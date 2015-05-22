#! /usr/env/python

"""
test_script_for_route_flow_dn.py: 
    
Tests and illustrates use of route_flow_dn component.
"""

from landlab.components.flow_routing.route_flow_dn import FlowRouter
from landlab.io import read_esri_ascii
from landlab.plot.imshow import imshow_node_grid
import os
import pylab

dem_name = './west_bijou_gully.asc'
outlet_row = 82
outlet_column = 38

# Read in a DEM and set its boundaries
DATA_FILE = os.path.join(os.path.dirname(__file__), dem_name)
(grid, z) = read_esri_ascii(DATA_FILE, name='topographic__elevation')
grid.set_nodata_nodes_to_closed(z, 0.) # set nodata nodes to inactive bounds
outlet_node = grid.grid_coords_to_node_id(outlet_row, outlet_column)

# Route flow
flow_router = FlowRouter(grid)
flow_router.route_flow()

# Create a shaded image
pylab.close()  # clear any pre-existing plot
pylab.figure(1)
im = imshow_node_grid(grid, 'water__volume_flux', cmap = pylab.cm.RdBu)

# add a title and axis labels
pylab.title('Discharge')
pylab.xlabel('Distance (m)')
pylab.ylabel('Distance (m)')

pylab.figure(2)
im = imshow_node_grid(grid, 'topographic__elevation')
pylab.title('DEM')
pylab.xlabel('Distance (m)')
pylab.ylabel('Distance (m)')

# Display the plot
pylab.show()
    
import numpy as np
print np.sum(grid.node_status!=4)
