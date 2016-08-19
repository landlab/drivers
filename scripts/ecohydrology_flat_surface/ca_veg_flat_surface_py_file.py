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
no_of_storms_dry = 8760 * (fraction_dry) / (data['mean_interstorm_dry'] +
                                            data['mean_storm_dry'])
n = int(n_years * (no_of_storms_wet + no_of_storms_dry))

(P, Tb, Tr, Time, veg_type, daily_pet, rad_factor,
 EP30, PET_threshold) = empty_arrays(n, grid, grid1)

create_pet_lookup(radiation, pet_tree, pet_shrub, pet_grass,  daily_pet,
                  rad_factor, EP30, grid)

# Represent current time in years
current_time = 0 # Start from first day of Jan

# Keep track of run time for simulation - optional
Start_time = time.clock() # Recording time taken for simulation

# declaring few variables that will be used in the storm loop
time_check = 0. # Buffer to store current_time at previous storm
yrs = 0 # Keep track of number of years passed
WS = 0. # Buffer for Water Stress
Tg = 270 # Growing season in days

# Run storm Loop
for i in range(n):
    # Update objects

    # Calculate Day of Year (DOY)
    Julian = np.int(np.floor((current_time - np.floor(current_time)) * 365.))

    # Generate seasonal storms
    # Wet Season - Jul to Sep - NA Monsoon
    if data['doy__start_of_monsoon'] <= Julian <= data['doy__end_of_monsoon']:
        precip_wet.update()
        P[i] = precip_wet.storm_depth
        Tr[i] = precip_wet.storm_duration
        Tb[i] = precip_wet.interstorm_duration
    else: # for Dry season
        precip_dry.update()
        P[i] = precip_dry.storm_depth
        Tr[i] = precip_dry.storm_duration
        Tb[i] = precip_dry.interstorm_duration

    # Spatially distribute PET and its 30-day-mean (analogous to degree day)
    grid.at_cell['surface__potential_evapotranspiration_rate'] = daily_pet[Julian]
    grid.at_cell['surface__potential_evapotranspiration_30day_mean'] = EP30[Julian]

    # Assign spatial rainfall data
    grid['cell']['rainfall__daily'] = P[i] * np.ones(grid.number_of_cells)

    # Update soil moisture component
    current_time = soil_moisture.update(current_time, Tr=Tr[i], Tb=Tb[i])

    # Decide whether its growing season or not
    if Julian != 364:
        if EP30[Julian + 1, 0] > EP30[Julian, 0]:
            PET_threshold = 1
            # 1 corresponds to ETThresholdup (begin growing season)
        else:
            PET_threshold = 0
            # 0 corresponds to ETThresholddown (end growing season)

    # Update vegetation component
    vegetation.update(PETThreshold_switch=PET_threshold, Tb=Tb[i], Tr=Tr[i])

    # Update yearly cumulative water stress data
    WS += (grid['cell']['vegetation__water_stress']) * Tb[i] / 24.

    # Record time (optional)
    Time[i] = current_time

    # Update spatial PFTs with Cellular Automata rules
    if (current_time - time_check) >= 1.:
        if yrs % 100 == 0:
            print 'Elapsed time = ', yrs, ' years'
        veg_type[yrs] = grid1.at_cell['vegetation__plant_functional_type']
        WS_ = np.choose(veg_type[yrs], WS)
        grid1.at_cell['vegetation__cumulative_water_stress'] = WS_ / Tg
        vegca.update()
        time_check = current_time
        WS = 0
        yrs += 1

veg_type[yrs] = grid1.at_cell['vegetation__plant_functional_type']

Final_time = time.clock()
Time_Consumed = (Final_time - Start_time) / 60. # in minutes
print 'Time_consumed = ', Time_Consumed, ' minutes'

# Saving
try:
    os.mkdir('output')
except OSError:
    pass
finally:
    os.chdir('output')

save('veg', Tb, Tr, P, veg_type, yrs, Time_Consumed, Time)

plot('veg', grid1, veg_type, yrs, yr_step=100)
