# -*- coding: utf-8 -*-
"""
A driver for our version of AW's gFlex component.

Created on Fri Feb 20 11:17:52 2015

@author: danhobley
"""

from landlab.components.diffusion.diffusion import LinearDiffuser
import numpy as np
import pylab
from pylab import show
from landlab import RasterModelGrid
from landlab import ModelParameterDictionary
from landlab.plot.imshow import imshow_node_grid

inputs = ModelParameterDictionary('./diffusion_params.txt')
nrows = inputs.read_int('nrows')
ncols = inputs.read_int('ncols')
dx = inputs.read_float('dx')
dt = inputs.read_float('dt')
time_to_run = inputs.read_float('run_time')
init_elev = inputs.read_float('init_elev')
uplift_rate = inputs.read_float('uplift_rate')

mg = RasterModelGrid(nrows, ncols, dx)

#create the fields in the grid
mg.create_node_array_zeros('topographic__elevation')
z = mg.create_node_array_zeros() + init_elev
mg['node'][ 'topographic__elevation'] = z + np.random.rand(len(z))/1000.

mg.set_fixed_value_boundaries_at_grid_edges(True, True, True, True)

#instantiate:
dfn = LinearDiffuser(mg, './diffusion_params.txt')

#perform the loop:
elapsed_time = 0. #total time in simulation
while elapsed_time < time_to_run:
    print elapsed_time
    if elapsed_time+dt>time_to_run:
        print "Short step!"
        dt = time_to_run - elapsed_time
    dfn.diffuse(dt)
    mg.at_node['topographic__elevation'][mg.core_nodes] += uplift_rate*dt
    elapsed_time += dt

pylab.figure(1)
im = imshow_node_grid(mg, 'topographic__elevation')  # display a colored image

pylab.figure(2)
im2 = pylab.plot(mg.node_vector_to_raster(mg.at_node['topographic__elevation'])[:,ncols//2])
pylab.xlabel('Horizontal distance')
pylab.ylabel('Elevation')
pylab.title('Cross section')
