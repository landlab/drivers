import numpy
from landlab import RasterModelGrid
from landlab import ModelParameterDictionary
from landlab.plot.imshow import imshow_node_grid
#from landlab.components.nonlinear_diffusion.Perron_nl_diffuse import PerronNLDiffuse
from landlab.components.nonlinear_diffusion.explicit_nl_diffuse import NonlinearDiffuser
import pylab
import time

inputs = ModelParameterDictionary('./drive_perron_params.txt')
nrows = inputs.read_int('nrows')
ncols = inputs.read_int('ncols')
dx = inputs.read_float('dx')
dt = inputs.read_float('dt')
time_to_run = inputs.read_float('run_time')
# nt needs defining
uplift = inputs.read_float('uplift_rate')
init_elev = inputs.read_float('init_elev')

mg = RasterModelGrid(nrows, ncols, dx)
# mg.set_looped_boundaries(True, True)
mg.set_closed_boundaries_at_grid_edges(True,True,True,True)

#create the fields in the grid
mg.create_node_array_zeros('topographic__elevation')
z = mg.create_node_array_zeros() + init_elev
mg.at_node[ 'topographic__elevation'] = z + numpy.random.rand(len(z))/1000.

# Display a message
print( 'Running ...' )
start_time = time.time()

#instantiate the component:
diffusion_component = NonlinearDiffuser(mg, './drive_perron_params.txt')

#perform the loop:
elapsed_time = 0. #total time in simulation
while elapsed_time < time_to_run:
    print elapsed_time
    mg.at_node['topographic__elevation'][mg.active_nodes[:(mg.active_nodes.shape[0]//2.)]] += uplift*dt #half block uplift
    if elapsed_time>time_to_run:
        dt = elapsed_time - time_to_run
    mg = diffusion_component.diffuse(dt)
    elapsed_time += dt

print('Total run time = '+str(time.time()-start_time)+' seconds.')

# Clear previous plots
pylab.figure(1)
pylab.close()

# Plot topography
pylab.figure(1)
im = imshow_node_grid(mg, 'topographic__elevation')
pylab.title('Topography')

elev_r = mg.node_vector_to_raster(mg.at_node['topographic__elevation'])
pylab.figure(2)
im = pylab.plot(dx*numpy.arange(nrows), elev_r[:,int(ncols//2)])  # display a colored image
pylab.title('Vertical cross section')

pylab.show()

print('Done.')
