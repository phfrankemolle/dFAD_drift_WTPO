This repository features code and figures for the manuscript on "Assessing the drift of Fish Aggregating Devices in the tropical Pacific Ocean".

The 'Code' folder features several python codes ran with the module parcels (for more information, see https://oceanparcels.org/).

This consists of 4 codes:
- WTPO_d50_FAD.py - A virtual dFAD simulation with depth averaged currents (from MOi dataset)
- WTPO_d0_FAD.py  - A virtual particle simulation with only surface currents (from MOi dataset)
- Processing_metrics.ipynb - A notebook file that computes four different metrics (displacement, travel distance, distance ratio, loopiness)
  based on the simulated trajectories and exports them as .nc files
- Analysis and Plot Metrics.ipynb - A notebook file that performs data analysis on the previously calculated metrics and plots the outcome

The 'Figures' folder features the figures from the manuscript of the simulations in pdf format.
