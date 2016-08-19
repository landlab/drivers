# Authors: Sai Nudurupati & Erkan Istanbulluoglu, 21May15
# Edited: 15Jul16 - to conform to Landlab version 1.
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

from landlab import load_params
from landlab.plot import imshow_grid
from landlab.components import (PrecipitationDistribution, Radiation,
                                PotentialEvapotranspiration, SoilMoisture,
                                Vegetation, VegCA)

GRASS = 0
SHRUB = 1
TREE = 2
BARE = 3
SHRUBSEEDLING = 4
TREESEEDLING = 5


def compose_veg_grid(grid, percent_bare=0.4, percent_grass=0.2,
                     percent_shrub=0.2, percent_tree=0.2):
    """Compose spatially distribute PFT."""
    no_cells = grid.number_of_cells
    shrub_point = int(percent_bare * no_cells)
    tree_point = int((percent_bare + percent_shrub) * no_cells)
    grass_point = int((1 - percent_grass) * no_cells)

    veg_grid = np.full(grid.number_of_cells, BARE, dtype=int)
    veg_grid[shrub_point:tree_point] = SHRUB
    veg_grid[tree_point:grass_point] = TREE
    veg_grid[grass_point:] = GRASS

    np.random.shuffle(veg_grid)
    return veg_grid


def initialize(data, grid, grid1):
    """Initialize random plant type field.

    Plant types are defined as the following:

    *  GRASS = 0
    *  SHRUB = 1
    *  TREE = 2
    *  BARE = 3
    *  SHRUBSEEDLING = 4
    *  TREESEEDLING = 5
    """
    grid1.at_cell['vegetation__plant_functional_type'] = compose_veg_grid(
        grid1, percent_bare=data['percent_bare_initial'],
        percent_grass=data['percent_grass_initial'],
        percent_shrub=data['percent_shrub_initial'],
        percent_tree=data['percent_tree_initial'])

    # Assign plant type for representative ecohydrologic simulations
    grid.at_cell['vegetation__plant_functional_type'] = np.arange(6)
    grid1.at_node['topographic__elevation'] = np.full(grid1.number_of_nodes,
                                                      1700.)
    grid.at_node['topographic__elevation'] = np.full(grid.number_of_nodes,
                                                     1700.)
    precip_dry = PrecipitationDistribution(
        mean_storm_duration=data['mean_storm_dry'],
        mean_interstorm_duration=data['mean_interstorm_dry'],
        mean_storm_depth=data['mean_storm_depth_dry'])
    precip_wet = PrecipitationDistribution(
        mean_storm_duration=data['mean_storm_wet'],
        mean_interstorm_duration=data['mean_interstorm_wet'],
        mean_storm_depth=data['mean_storm_depth_wet'])

    radiation = Radiation(grid)
    pet_tree = PotentialEvapotranspiration(grid, method=data['PET_method'],
                                           MeanTmaxF=data['MeanTmaxF_tree'],
                                           delta_d=data['DeltaD'])
    pet_shrub = PotentialEvapotranspiration(grid, method=data['PET_method'],
                                            MeanTmaxF=data['MeanTmaxF_shrub'],
                                            delta_d=data['DeltaD'])
    pet_grass = PotentialEvapotranspiration(grid, method=data['PET_method'],
                                            MeanTmaxF=data['MeanTmaxF_grass'],
                                            delta_d=data['DeltaD'])
    soil_moisture = SoilMoisture(grid, **data) # Soil Moisture object
    vegetation = Vegetation(grid, **data) # Vegetation object
    vegca = VegCA(grid1, **data) # Cellular automaton object

    # Initializing inputs for Soil Moisture object
    grid.at_cell['vegetation__live_leaf_area_index'] = (
        1.6 * np.ones(grid.number_of_cells))
    grid.at_cell['soil_moisture__initial_saturation_fraction'] = (
        0.59 * np.ones(grid.number_of_cells))

    return (precip_dry, precip_wet, radiation, pet_tree, pet_shrub,
            pet_grass, soil_moisture, vegetation, vegca)


def empty_arrays(n, grid, grid1):
    precip = np.empty(n) # Record precipitation
    inter_storm_dt = np.empty(n) # Record inter storm duration
    storm_dt = np.empty(n) # Record storm duration
    time_elapsed = np.empty(n) # To record time elapsed from the start of simulation

    # Cumulative Water Stress
    veg_type = np.empty([n / 55, grid1.number_of_cells], dtype=int)
    daily_pet = np.zeros([365, grid.number_of_cells])
    rad_factor = np.empty([365, grid.number_of_cells])
    EP30 = np.empty([365, grid.number_of_cells])

    # 30 day average PET to determine season
    pet_threshold = 0  # Initializing pet_threshold to ETThresholddown
    return (precip, inter_storm_dt, storm_dt, time_elapsed, veg_type,
            daily_pet, rad_factor, EP30, pet_threshold)


def create_pet_lookup(radiation, pet_tree, pet_shrub, pet_grass, daily_pet,
                      rad_factor, EP30, grid):
    for i in range(0, 365):
        pet_tree.update(float(i) / 365.25)
        pet_shrub.update(float(i) / 365.25)
        pet_grass.update(float(i) / 365.25)
        daily_pet[i] = [pet_grass._PET_value, pet_shrub._PET_value,
                   pet_tree._PET_value, 0., pet_shrub._PET_value,
                   pet_tree._PET_value]
        radiation.update(float(i) / 365.25)
        rad_factor[i] = grid.at_cell['radiation__ratio_to_flat_surface']

        if i < 30:
            if i == 0:
                EP30[0] = daily_pet[0]
            else:
                EP30[i] = np.mean(daily_pet[:i], axis=0)
        else:
            EP30[i] = np.mean(daily_pet[i - 30:i], axis=0)


def save(sim, inter_storm_dt, storm_dt, precip, veg_type, yrs,
         walltime, time_elapsed):
    np.save(sim + '_Tb', inter_storm_dt)
    np.save(sim + '_Tr', storm_dt)
    np.save(sim + '_P', precip)
    np.save(sim + '_VegType', veg_type)
    np.save(sim + '_Years', yrs)
    np.save(sim + '_Time_Consumed_minutes', walltime)
    np.save(sim + '_CurrentTime', time_elapsed)


def plot(sim, grid, veg_type, yrs, yr_step=10):
    pic = 0
    years = range(0, yrs)
    cmap = mpl.colors.ListedColormap(
        ['green', 'red', 'black', 'white', 'red', 'black'])
    bounds = [-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5]
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
    print 'Plotting cellular field of Plant Functional Type'
    print 'Green - Grass; Red - Shrubs; Black - Trees; White - Bare'

    # Plot images to make gif.
    for year in range(0, yrs, yr_step):
        filename = 'year_' + "%05d" % year
        pic += 1
        plt.figure(pic, figsize=(10, 8))
        imshow_grid(grid, veg_type[year], values_at='cell', cmap=cmap,
                    grid_units=('m', 'm'), norm=norm, limits=[0, 5],
                    allow_colorbar=False)
        plt.title(filename, weight='bold', fontsize=22)
        plt.xlabel('X (m)', weight='bold', fontsize=18)
        plt.ylabel('Y (m)', weight='bold', fontsize=18)
        plt.xticks(fontsize=14, weight='bold')
        plt.yticks(fontsize=14, weight='bold')
        plt.savefig(sim + '_' + filename)

    grass_cov = np.empty(yrs)
    shrub_cov = np.empty(yrs)
    tree_cov = np.empty(yrs)
    grid_size = float(veg_type.shape[1])

    for x in range(0, yrs):
        grass_cov[x] = (veg_type[x][veg_type[x] == GRASS].size / grid_size) * 100
        shrub_cov[x] = ((veg_type[x][veg_type[x] == SHRUB].size / grid_size) *
                        100 + (veg_type[x][veg_type[x] == SHRUBSEEDLING].size /
                        grid_size) * 100)
        tree_cov[x] = ((veg_type[x][veg_type[x] == TREE].size / grid_size) *
                       100 + (veg_type[x][veg_type[x] == TREESEEDLING].size /
                       grid_size) * 100)

    pic += 1
    plt.figure(pic, figsize=(10, 8))
    plt.plot(years, grass_cov, '-g', label='Grass', linewidth=4)
    plt.hold(True)
    plt.plot(years, shrub_cov, '-r', label='Shrub', linewidth=4)
    plt.hold(True)
    plt.plot(years, tree_cov, '-k', label='Tree', linewidth=4)
    plt.ylabel('% Area Covered by Plant Type', weight='bold', fontsize=18)
    plt.xlabel('Time in years', weight='bold', fontsize=18)
    plt.xticks(fontsize=12, weight='bold')
    plt.yticks(fontsize=12, weight='bold')
    plt.legend(loc=0, prop={'size': 16, 'weight': 'bold'})
    plt.savefig(sim + '_percent_cover')
