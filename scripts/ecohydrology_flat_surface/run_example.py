"""
Created on Wed Jul 20 2016

This tutorial is on:
landlab/tutorials/ecohydrology/cellular_automaton_vegetation_flat_surface.ipynb

Creating a (.py) version of the same.

@author: Sai Nudurupati & Erkan Istanbulluoglu
"""

import os
import time

import numpy as np

from landlab import RasterModelGrid, load_params
from ecohyd_functions_flat import (initialize, empty_arrays,
                                   create_pet_lookup, save, plot)


grid1 = RasterModelGrid((100, 100), spacing=(5., 5.))
grid = RasterModelGrid((5, 4), spacing=(5., 5.))

# Create dictionary that holds the inputs
data = load_params('inputs_vegetation_ca.yaml')

(precip_dry, precip_wet, radiation, pet_tree, pet_shrub,
 pet_grass, soil_moisture, vegetation, vegca) = initialize(data, grid, grid1)

n_years = 2000 # Approx number of years for model to run

# Calculate approximate number of storms per year
fraction_wet = (data['doy__end_of_monsoon'] -
                data['doy__start_of_monsoon']) / 365.
fraction_dry = 1 - fraction_wet
no_of_storms_wet = 8760 * fraction_wet / (data['mean_interstorm_wet'] +
                                          data['mean_storm_wet'])
no_of_storms_dry = 8760 * fraction_dry / (data['mean_interstorm_dry'] +
                                          data['mean_storm_dry'])
n = int(n_years * (no_of_storms_wet + no_of_storms_dry))

(precip, inter_storm_dt, storm_dt, time_elapsed, veg_type, daily_pet,
 rad_factor, EP30, pet_threshold) = empty_arrays(n, grid, grid1)

create_pet_lookup(radiation, pet_tree, pet_shrub, pet_grass,  daily_pet,
                  rad_factor, EP30, grid)

# Represent current time in years
current_time = 0 # Start from first day of Jan

# Keep track of run time for simulation - optional
wallclock_start = time.clock() # Recording time taken for simulation

# declaring few variables that will be used in the storm loop
time_check = 0. # Buffer to store current_time at previous storm
yrs = 0 # Keep track of number of years passed
water_stress = 0. # Buffer for Water Stress
Tg = 270 # Growing season in days

# Run storm Loop
for i in range(n):
    # Update objects

    # Calculate Day of Year (DOY)
    julian = np.int(np.floor((current_time - np.floor(current_time)) * 365.))

    # Generate seasonal storms
    # Wet Season - Jul to Sep - NA Monsoon
    if data['doy__start_of_monsoon'] <= julian <= data['doy__end_of_monsoon']:
        precip_wet.update()
        precip[i] = precip_wet.storm_depth
        storm_dt[i] = precip_wet.storm_duration
        inter_storm_dt[i] = precip_wet.interstorm_duration
    else: # for Dry season
        precip_dry.update()
        precip[i] = precip_dry.storm_depth
        storm_dt[i] = precip_dry.storm_duration
        inter_storm_dt[i] = precip_dry.interstorm_duration

    # Spatially distribute PET and its 30-day-mean (analogous to degree day)
    grid.at_cell['surface__potential_evapotranspiration_rate'] = daily_pet[julian]
    grid.at_cell['surface__potential_evapotranspiration_30day_mean'] = EP30[julian]

    # Assign spatial rainfall data
    grid.at_cell['rainfall__daily_depth'] = np.full(grid.number_of_cells, precip[i])

    # Update soil moisture component
    current_time = soil_moisture.update(current_time, Tr=storm_dt[i],
                                        Tb=inter_storm_dt[i])

    # Decide whether its growing season or not
    if julian != 364:
        if EP30[julian + 1, 0] > EP30[julian, 0]:
            pet_threshold = 1
            # 1 corresponds to ETThresholdup (begin growing season)
        else:
            pet_threshold = 0
            # 0 corresponds to ETThresholddown (end growing season)

    # Update vegetation component
    vegetation.update(PETThreshold_switch=pet_threshold, Tb=inter_storm_dt[i],
                      Tr=storm_dt[i])

    # Update yearly cumulative water stress data
    water_stress += (grid.at_cell['vegetation__water_stress'] *
                     inter_storm_dt[i] / 24.)

    # Record time (optional)
    time_elapsed[i] = current_time

    # Update spatial PFTs with Cellular Automata rules
    if (current_time - time_check) >= 1.:
        if yrs % 100 == 0:
            print 'Elapsed time = ', yrs, ' years'
        veg_type[yrs] = grid1.at_cell['vegetation__plant_functional_type']
        WS_ = np.choose(veg_type[yrs], water_stress)
        grid1.at_cell['vegetation__cumulative_water_stress'] = WS_ / Tg
        vegca.update()
        time_check = current_time
        water_stress = 0
        yrs += 1

veg_type[yrs] = grid1.at_cell['vegetation__plant_functional_type']

wallclock_stop = time.clock()
walltime = (wallclock_stop - wallclock_start) / 60. # in minutes
print 'Time_consumed = ', walltime, ' minutes'

# Saving
try:
    os.mkdir('output')
except OSError:
    pass
finally:
    os.chdir('output')

save('veg', inter_storm_dt, storm_dt, precip, veg_type, yrs,
     walltime, time_elapsed)

plot('veg', grid1, veg_type, yrs, yr_step=100)
