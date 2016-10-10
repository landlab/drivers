
# Authors: Sai Nudurupati & Erkan Istanbulluoglu, 21May15
# Edited: 15Jul16 - to conform to Landlab version 1.
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from landlab.plot import imshow_grid
from landlab.components import PrecipitationDistribution
from landlab.components import Radiation
from landlab.components import PotentialEvapotranspiration
from landlab.components import SoilMoisture
from landlab.components import Vegetation
from landlab.components import VegCA

GRASS = 0
SHRUB = 1
TREE = 2
BARE = 3
SHRUBSEEDLING = 4
TREESEEDLING = 5


# Function that converts text file to a dictionary
def txt_data_dict(InputFile):
    f = open(InputFile)
    data1 = {}
    for line in f:
        if line.strip() != '' and line[0] != '#':
            m, n = line.split(':')
            line = f.next()
            e = line[:].strip()
            if e[0].isdigit():
                if e.find('.') != -1:
                    data1[m.strip()] = float(line[:].strip())
                else:
                    data1[m.strip()] = int(line[:].strip())
            else:
                data1[m.strip()] = line[:].strip()
    f.close()
    return data1.copy()


# Function to compose spatially distribute PFT
def compose_veg_grid(grid, percent_bare=0.4, percent_grass=0.2,
                     percent_shrub=0.2, percent_tree=0.2):
    no_cells = grid.number_of_cells
    V = 3 * np.ones(grid.number_of_cells, dtype=int)
    shrub_point = int(percent_bare * no_cells)
    tree_point = int((percent_bare + percent_shrub) * no_cells)
    grass_point = int((1 - percent_grass) * no_cells)
    V[shrub_point:tree_point] = 1
    V[tree_point:grass_point] = 2
    V[grass_point:] = 0
    np.random.shuffle(V)
    return V


def Initialize_(data, grid, grid1, elevation):
    # Plant types are defined as following:
    # GRASS = 0; SHRUB = 1; TREE = 2; BARE = 3;
    # SHRUBSEEDLING = 4; TREESEEDLING = 5
    # Initialize random plant type field
    grid['cell']['vegetation__plant_functional_type'] = compose_veg_grid(
                grid, percent_bare=data['percent_bare_initial'],
                percent_grass=data['percent_grass_initial'],
                percent_shrub=data['percent_shrub_initial'],
                percent_tree=data['percent_tree_initial'])
    # Assign plant type for representative ecohydrologic simulations
    grid1['cell']['vegetation__plant_functional_type'] = np.arange(0, 6)
    grid1['node']['topographic__elevation'] = (1700. *
                                               np.ones(grid1.number_of_nodes))
    grid['node']['topographic__elevation'] = elevation
    PD_D = PrecipitationDistribution(
                        mean_storm_duration=data['mean_storm_dry'],
                        mean_interstorm_duration=data['mean_interstorm_dry'],
                        mean_storm_depth=data['mean_storm_depth_dry'])
    PD_W = PrecipitationDistribution(
                        mean_storm_duration=data['mean_storm_wet'],
                        mean_interstorm_duration=data['mean_interstorm_wet'],
                        mean_storm_depth=data['mean_storm_depth_wet'])
    Rad = Radiation(grid)
    Rad_PET = Radiation(grid1)
    PET_Tree = PotentialEvapotranspiration(grid1, method=data['PET_method'],
                                           MeanTmaxF=data['MeanTmaxF_tree'],
                                           delta_d=data['DeltaD'])
    PET_Shrub = PotentialEvapotranspiration(grid1, method=data['PET_method'],
                                            MeanTmaxF=data['MeanTmaxF_shrub'],
                                            delta_d=data['DeltaD'])
    PET_Grass = PotentialEvapotranspiration(grid1, method=data['PET_method'],
                                            MeanTmaxF=data['MeanTmaxF_grass'],
                                            delta_d=data['DeltaD'])
    SM = SoilMoisture(grid, **data)   # Soil Moisture object
    VEG = Vegetation(grid, **data)    # Vegetation object
    vegca = VegCA(grid, **data)      # Cellular automaton object

    # # Initializing inputs for Soil Moisture object
    grid['cell']['vegetation__live_leaf_area_index'] = (
                                    1.6 * np.ones(grid.number_of_cells))
    grid['cell']['soil_moisture__initial_saturation_fraction'] = (
                                    0.59 * np.ones(grid.number_of_cells))
    # Initializing Soil Moisture
    return PD_D, PD_W, Rad, Rad_PET, PET_Tree, PET_Shrub, PET_Grass, SM, \
        VEG, vegca


def Empty_arrays(n, n_years, grid, grid1):
    P = np.empty(n)    # Record precipitation
    Tb = np.empty(n)    # Record inter storm duration
    Tr = np.empty(n)    # Record storm duration
    Time = np.empty(n)  # To record time elapsed from the start of simulation
#    CumWaterStress = np.empty(grid.number_of_cells)  # Cum Water Stress
    VegType = np.empty([n_years+5, grid.number_of_cells], dtype=int)
    PET_ = np.zeros([365, grid1.number_of_cells])
    Rad_Factor = np.empty([365, grid.number_of_cells])
    EP30 = np.empty([365, grid1.number_of_cells])
    # 30 day average PET to determine season
    PET_threshold = 0  # Initializing PET_threshold to ETThresholddown
    return (P, Tb, Tr, Time, VegType,
            PET_, Rad_Factor, EP30, PET_threshold)


def Create_PET_lookup(Rad, PET_Tree, PET_Shrub, PET_Grass, PET_,
                      Rad_Factor, EP30, Rad_PET, grid):
    for i in range(0, 365):
        Rad_PET.update(float(i)/365.25)
        PET_Tree.update(float(i)/365.25)
        PET_Shrub.update(float(i)/365.25)
        PET_Grass.update(float(i)/365.25)
        PET_[i] = [PET_Grass._PET_value, PET_Shrub._PET_value,
                   PET_Tree._PET_value, 0., PET_Shrub._PET_value,
                   PET_Tree._PET_value]
        Rad.update(float(i)/365.25)
        Rad_Factor[i] = grid['cell']['radiation__ratio_to_flat_surface']
        if i < 30:
            if i == 0:
                EP30[0] = PET_[0]
            else:
                EP30[i] = np.mean(PET_[:i], axis=0)
        else:
            EP30[i] = np.mean(PET_[i-30:i], axis=0)


def Save_(sim, Tb, Tr, P, VegType, yrs, Time_Consumed, Time):
    np.save(sim+'Tb', Tb)
    np.save(sim+'Tr', Tr)
    np.save(sim+'P', P)
    np.save(sim+'VegType', VegType)
#    np.save(sim+'CumWaterStress', CumWaterStress)
    np.save(sim+'Years', yrs)
    np.save(sim+'Time_Consumed_minutes', Time_Consumed)
    np.save(sim+'CurrentTime', Time)


def Plot_(sim, grid, VegType, yrs, yr_step=10):
    # # Plotting
    elevation = grid['node']['topographic__elevation']
    ELEV = grid.node_vector_to_raster(elevation)
    Elev = np.zeros([ELEV.shape[0]-2, ELEV.shape[1]-2])
    xx = (grid.dx * (np.arange(Elev.shape[1])))
    yy = (grid.dy * (np.arange(Elev.shape[0])))
    for ii in range(1, int(ELEV.shape[0]-1)):
        for jj in range(1, int(ELEV.shape[1]-1)):
            Elev[ii-1][jj-1] = ELEV[ii][jj]

    pic = 0
    years = range(0, yrs)
    cmap = mpl.colors.ListedColormap(
                        ['green', 'red', 'black', 'white', 'red', 'black'])
    bounds = [-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5]
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
    print 'Plotting cellular field of Plant Functional Type'
    print 'Green - Grass; Red - Shrubs; Black - Trees; White - Bare'
    # # Plot images to make gif.
    for year in range(0, yrs, yr_step):
        filename = 'Year = ' + "%05d" % year
        pic += 1
        plt.figure(pic, figsize=(10, 8))
        imshow_grid(grid, VegType[year], values_at='cell', cmap=cmap,
                    grid_units=('m', 'm'), norm=norm, limits=[0, 5],
                    allow_colorbar=False)
        plt.hold(True)
        CLT = plt.contour(xx, yy, Elev, colors='b',
                          figsize=(10, 8), linewidths=3)
        plt.clabel(CLT, inline=1, fmt='%4.0f', fontsize=12)
        plt.title(filename, weight='bold', fontsize=14)
        plt.savefig(sim+filename)

    grass_cov = np.empty(yrs)
    shrub_cov = np.empty(yrs)
    tree_cov = np.empty(yrs)
    grid_size = float(VegType.shape[1])
    for x in range(0, yrs):
        grass_cov[x] = (VegType[x][VegType[x] == GRASS].size/grid_size) * 100
        shrub_cov[x] = ((VegType[x][VegType[x] == SHRUB].size/grid_size) *
                        100 + (VegType[x][VegType[x] == SHRUBSEEDLING].size /
                        grid_size) * 100)
        tree_cov[x] = ((VegType[x][VegType[x] == TREE].size/grid_size) *
                       100 + (VegType[x][VegType[x] == TREESEEDLING].size /
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
    plt.savefig(sim+'percentCover')
    # plt.show()
