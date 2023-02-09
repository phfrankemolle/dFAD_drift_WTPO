# -*- coding: utf-8 -*-
"""
@author: Philippe Frankemölle

Virtual surface drifter simulations 
"""

from parcels import (FieldSet, AdvectionRK4, ParticleFile, ParticleSet, JITParticle, Variable, ErrorCode)
from datetime import timedelta as delta
import numpy as np

class FADParticle(JITParticle):
    drift_time = Variable('drift_time', dtype=np.float32, to_write='True')
def OutOfBounds(particle, fieldset, time):
    particle.delete()
def DriftTime(particle, fieldset, time):
    particle.drift_time = particle.drift_time + 6*3600
    if particle.drift_time>= 180*24*3600:
        particle.delete()
def PeriodicBoundary1(particle, fieldset, time):
    if particle.lon > 180.:
        particle.lon += -360
def PeriodicBoundary2(particle, fieldset, time):
    if particle.lon < -180.:
        particle.lon += 360
def OffTheGrid1(particle, fieldset, time):
    if particle.lat > 31:
        particle.delete()

#load MOi data for horizontal velocities over 50 meters depth
def set_MOi_50m_fieldset():
    data_path = '/storage/shared/oceanparcels/input_data/'
    
    listdates= np.arange('2006-11-01', '2022-01-01', dtype='datetime64[D]') #time range for loaded data
    file_path_U = [] #list for data path to zonal velocity data
    file_path_V = [] #list for data path to meridional velocity data
    for i in range(len(listdates)):
        file_path_U.append(data_path + 'MOi/psy4v3r1/psy4v3r1-daily_U_' + str(listdates[i]) +'.nc')
        file_path_V.append(data_path + 'MOi/psy4v3r1/psy4v3r1-daily_V_' + str(listdates[i]) +'.nc')
    filenames = {'U': {'lon': data_path + 'NEMO-MEDUSA/ORCA0083-N006/domain/coordinates.nc',
                       'lat': data_path + 'NEMO-MEDUSA/ORCA0083-N006/domain/coordinates.nc',
                       'depth': data_path + 'MOi/psy4v3r1/psy4v3r1-daily_W_' + str(listdates[0]) +'.nc',
                       'data': file_path_U},
                 'V': {'lon': data_path + 'NEMO-MEDUSA/ORCA0083-N006/domain/coordinates.nc',
                       'lat': data_path + 'NEMO-MEDUSA/ORCA0083-N006/domain/coordinates.nc',
                       'depth': data_path + 'MOi/psy4v3r1/psy4v3r1-daily_W_' + str(listdates[0]) +'.nc',
                       'data': file_path_V}}
    variables = {'U': 'vozocrtx', 'V': 'vomecrty'}
    dimensions = {'lon': 'glamf', 'lat': 'gphif', 'depth': 'depthw','time': 'time_counter'}
    #indices that demarcate the tropical western pacific
    #indices = {'lon':range(500,1300),'lat': range(1100,1900), 'depth': range(0, 18)}
    indices = {'lat': range(1100,1900), 'depth': range(0, 18)}
    #indices = {'depth': range(0, 18)}
    fset = FieldSet.from_nemo(filenames, variables, dimensions, indices,allow_time_extrapolation=True)
    fset.add_constant('max_drift_time', delta(days=3).total_seconds())
    return fset

def run_particles():
    fieldset = set_MOi_50m_fieldset()
    dz = np.gradient(fieldset.U.depth)
    DZ = np.moveaxis(np.tile(dz, (fieldset.U.grid.ydim, fieldset.U.grid.xdim, 1)), [0, 1, 2], [1, 2, 0])

    def compute(fieldset):
        # Calculating vertical weighted average
        for f in [fieldset.U, fieldset.V]:
            for tind in f.loaded_time_indices:
                f.data[tind, :] = np.sum(f.data[tind, :] * DZ, axis=0) / sum(dz)

    fieldset.compute_on_defer = compute
    
    latrange = np.linspace(-15,15,num = 61) #0.5 graden
    lonrange = np.linspace(140,285,num = 291) #0.5 graden
    latmesh,lonmesh = np.meshgrid(latrange,lonrange)
    latlist,lonlist = np.reshape(latmesh,-1),np.reshape(lonmesh,-1)
    lonlist[lonlist >180] = lonlist[lonlist >180] - 360 
    repeatdt = delta(days=10)
    
    pset = ParticleSet.from_list(fieldset=fieldset, pclass=FADParticle, #Practice particle dataset.
                             lon = lonlist, lat = latlist, depth=0.5*np.ones(len(lonlist)), repeatdt=repeatdt)

    ofile = ParticleFile("20062021depth_wholepacific.nc", pset, outputdt=delta(days=1))

    kernels = pset.Kernel(AdvectionRK4) + PeriodicBoundary1 + PeriodicBoundary2 + OffTheGrid1 + DriftTime
    pset.execute(kernels, runtime=delta(days=5539), dt=delta(hours=6), output_file=ofile,
                 recovery={ErrorCode.ErrorOutOfBounds: OutOfBounds})

run_particles()