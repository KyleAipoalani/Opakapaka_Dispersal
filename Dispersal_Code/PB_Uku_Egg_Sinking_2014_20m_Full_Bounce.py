#!/usr/bin/env python
# coding: utf-8

# In[ ]:


#!/usr/bin/env python
# coding: utf-8


# In[26]:

# In[8]:


from parcels import FieldSet, Field, ParticleSet, JITParticle, AdvectionRK4, Variable, plotTrajectoriesFile, ErrorCode, DiffusionUniformKh, GeographicPolar,Geographic
from datetime import timedelta, datetime
from netCDF4 import Dataset
import numpy as np
import pandas as pd
import netCDF4
import matplotlib.pyplot as plt
import xarray
from scipy import interpolate


# In[27]:<br>
# Define the new Kernel that mimics egg hatching for snapper

# In[9]:


def EggHatchingMovement(particle, fieldset, time):
    eggdepth = 0.25 # egg depth in m
    larvadepth = 20  # larval drift depth in m
    drifttime =  45* 86400  # time of deep drift in seconds
    eggtime=1*86400
    vertical_speed = 0.02  # sink and rise speed in m/s
    driftdepth=20

    if particle.cycle_phase == 0:
        # Phase 0:particle is an egg and has yet to hatch
        particle.age += particle.dt
        particle.depth=eggdepth 
        if particle.age >= eggtime:
            particle.cycle_phase = 1
            
    elif particle.cycle_phase == 1:
        # Phase 0: Sinking with vertical_speed until depth is driftdepth
        particle.depth += vertical_speed * particle.dt
        if particle.depth >= driftdepth:
            particle.cycle_phase = 2

    elif particle.cycle_phase == 2:
        # Phase 1: Drifting at larval depth
        particle.age += particle.dt
        particle.depth=larvadepth 
        if particle.age > eggtime+drifttime:
            particle.delete()
        


# In[10]:


#break


# In[11]:


filenames = {'U': "/home/jsuca/MHI_ROMS_For_Uku_Spawning/2014/0.25_1_5_10_20_30_50m/*.nc",
             'V': "/home/jsuca/MHI_ROMS_For_Uku_Spawning/2014/0.25_1_5_10_20_30_50m/*.nc"}


variables = {'U': 'u',
             'V': 'v'}
dimensions = {'U': {'lat': 'latitude', 'lon': 'longitude', 'depth': 'depth', 'time': 'time'},
               'V': {'lat': 'latitude', 'lon': 'longitude', 'depth': 'depth', 'time': 'time'}}


# In[29]:

# In[14]:


fieldset = FieldSet.from_netcdf(filenames, variables, dimensions,interp_method={'U': 'freeslip', 'V': 'freeslip'})

# In[12]:

file_path_fine = '/home/jsuca/MHI_ROMS_For_Uku_Spawning/2014/0.25_1_5_10_20_30_50m/Yr_2014_04_05.nc'

def make_landmask(fielddata):
    """Returns landmask where land = 1 and ocean = 0
    fielddata is a netcdf file.
    """
    datafile = Dataset(fielddata)

    landmask = datafile.variables['u'][0, 0]
    landmask = np.ma.masked_invalid(landmask) #remove Nas? 
    landmask = landmask.mask.astype('int')

    return landmask

#


landmask_fine = make_landmask(file_path_fine)


# In[28]:
def get_coastal_nodes(landmask):
    """Function that detects the coastal nodes, i.e. the ocean nodes directly
    next to land. Computes the Laplacian of landmask.

    - landmask: the land mask built using `make_landmask`, where land cell = 1
                and ocean cell = 0.

    Output: 2D array array containing the coastal nodes, the coastal nodes are
            equal to one, and the rest is zero.
    """
    mask_lap = np.roll(landmask, -1, axis=0) + np.roll(landmask, 1, axis=0)
    mask_lap += np.roll(landmask, -1, axis=1) + np.roll(landmask, 1, axis=1)
    mask_lap -= 4*landmask
    coastal = np.ma.masked_array(landmask, mask_lap > 0)
    coastal = coastal.mask.astype('int')

    return coastal

def get_shore_nodes(landmask):
    """Function that detects the shore nodes, i.e. the land nodes directly
    next to the ocean. Computes the Laplacian of landmask.

    - landmask: the land mask built using `make_landmask`, where land cell = 1
                and ocean cell = 0.

    Output: 2D array array containing the shore nodes, the shore nodes are
            equal to one, and the rest is zero.
    """
    mask_lap = np.roll(landmask, -1, axis=0) + np.roll(landmask, 1, axis=0)
    mask_lap += np.roll(landmask, -1, axis=1) + np.roll(landmask, 1, axis=1)
    mask_lap -= 4*landmask
    shore = np.ma.masked_array(landmask, mask_lap < 0)
    shore = shore.mask.astype('int')

    return shore

# In[13]:

def get_coastal_nodes_diagonal(landmask):
    """Function that detects the coastal nodes, i.e. the ocean nodes where 
    one of the 8 nearest nodes is land. Computes the Laplacian of landmask
    and the Laplacian of the 45 degree rotated landmask.

    - landmask: the land mask built using `make_landmask`, where land cell = 1
                and ocean cell = 0.

    Output: 2D array array containing the coastal nodes, the coastal nodes are
            equal to one, and the rest is zero.
    """
    mask_lap = np.roll(landmask, -1, axis=0) + np.roll(landmask, 1, axis=0)
    mask_lap += np.roll(landmask, -1, axis=1) + np.roll(landmask, 1, axis=1)
    mask_lap += np.roll(landmask, (-1,1), axis=(0,1)) + np.roll(landmask, (1, 1), axis=(0,1))
    mask_lap += np.roll(landmask, (-1,-1), axis=(0,1)) + np.roll(landmask, (1, -1), axis=(0,1))
    mask_lap -= 8*landmask
    coastal = np.ma.masked_array(landmask, mask_lap > 0)
    coastal = coastal.mask.astype('int')
    
    return coastal
    
def get_shore_nodes_diagonal(landmask):
    """Function that detects the shore nodes, i.e. the land nodes where 
    one of the 8 nearest nodes is ocean. Computes the Laplacian of landmask 
    and the Laplacian of the 45 degree rotated landmask.

    - landmask: the land mask built using `make_landmask`, where land cell = 1
                and ocean cell = 0.

    Output: 2D array array containing the shore nodes, the shore nodes are
            equal to one, and the rest is zero.
    """
    mask_lap = np.roll(landmask, -1, axis=0) + np.roll(landmask, 1, axis=0)
    mask_lap += np.roll(landmask, -1, axis=1) + np.roll(landmask, 1, axis=1)
    mask_lap += np.roll(landmask, (-1,1), axis=(0,1)) + np.roll(landmask, (1, 1), axis=(0,1))
    mask_lap += np.roll(landmask, (-1,-1), axis=(0,1)) + np.roll(landmask, (1, -1), axis=(0,1))
    mask_lap -= 8*landmask
    shore = np.ma.masked_array(landmask, mask_lap < 0)
    shore = shore.mask.astype('int')

    return shore
#
coastal_fine = get_coastal_nodes_diagonal(landmask_fine)
shore_fine = get_shore_nodes_diagonal(landmask_fine)

#
def create_displacement_field(landmask, double_cell=False):
    """Function that creates a displacement field 1 m/s away from the shore.

    - landmask: the land mask dUilt using `make_landmask`.
    - double_cell: Boolean for determining if you want a double cell.
      Default set to False.

    Output: two 2D arrays, one for each camponent of the velocity.
    """
    shore = get_shore_nodes(landmask)
    shore_d = get_shore_nodes_diagonal(landmask) # bordering ocean directly and diagonally
    shore_c = shore_d - shore                    # corner nodes that only border ocean diagonally
    
    Ly = np.roll(landmask, -1, axis=0) - np.roll(landmask, 1, axis=0) # Simple derivative
    Lx = np.roll(landmask, -1, axis=1) - np.roll(landmask, 1, axis=1)
    
    Ly_c = np.roll(landmask, -1, axis=0) - np.roll(landmask, 1, axis=0)
    Ly_c += np.roll(landmask, (-1,-1), axis=(0,1)) + np.roll(landmask, (-1,1), axis=(0,1)) # Include y-component of diagonal neighbours
    Ly_c += - np.roll(landmask, (1,-1), axis=(0,1)) - np.roll(landmask, (1,1), axis=(0,1))
    
    Lx_c = np.roll(landmask, -1, axis=1) - np.roll(landmask, 1, axis=1)
    Lx_c += np.roll(landmask, (-1,-1), axis=(1,0)) + np.roll(landmask, (-1,1), axis=(1,0)) # Include x-component of diagonal neighbours
    Lx_c += - np.roll(landmask, (1,-1), axis=(1,0)) - np.roll(landmask, (1,1), axis=(1,0))
    
    v_x = -Lx*(shore)
    v_y = -Ly*(shore)
    
    v_x_c = -Lx_c*(shore_c)
    v_y_c = -Ly_c*(shore_c)
    
    v_x = v_x + v_x_c
    v_y = v_y + v_y_c

    magnitude = np.sqrt(v_y**2 + v_x**2)
    # the coastal nodes between land create a problem. Magnitude there is zero
    # I force it to be 1 to avoid problems when normalizing.
    ny, nx = np.where(magnitude == 0)
    magnitude[ny, nx] = 1

    v_x = v_x/magnitude
    v_y = v_y/magnitude

    return v_x, v_y



##
v_x_f, v_y_f = create_displacement_field(landmask_fine)
#

def distance_to_shore(landmask, dx=1):
    """Function that computes the distance to the shore. It is based in the
    the `get_coastal_nodes` algorithm.

    - landmask: the land mask dUilt using `make_landmask` function.
    - dx: the grid cell dimension. This is a crude approxsimation of the real
    distance (be careful).

    Output: 2D array containing the distances from shore.
    """
    ci = get_coastal_nodes(landmask) # direct neighbours
    dist = ci*dx                     # 1 dx away
    
    ci_d = get_coastal_nodes_diagonal(landmask) # diagonal neighbours
    dist_d = (ci_d - ci)*np.sqrt(2*dx**2)       # sqrt(2) dx away
        
    return dist+dist_d

#
d_2_s_f = distance_to_shore(landmask_fine)

#
def set_displacement(particle, fieldset, time):  
        particle.d2s = fieldset.distance2shore_fine[time, particle.depth,
                               particle.lat, particle.lon]
        if  particle.d2s < 0.5:
            dispUab = fieldset.dispUF[time, particle.depth, particle.lat,
                               particle.lon]
            dispVab = fieldset.dispVF[time, particle.depth, particle.lat,
                               particle.lon]
            
            particle.dU = dispUab
            particle.dV = dispVab
        else:
            particle.dU = 0.
            particle.dV = 0.


##
def displace(particle, fieldset, time):    
    if  particle.d2s < 0.5:
        particle.lon += particle.dU*particle.dt
        particle.lat += particle.dV*particle.dt
##
u_displacement_f = v_x_f
v_displacement_f = v_y_f
#
fieldset.add_field(Field('dispUF', data=u_displacement_f,
                         lon=fieldset.U.grid.lon, lat=fieldset.U.grid.lat,
                         mesh='spherical')) #have to index to choose which field we want to base it off of; 1 is choosing coarser

fieldset.add_field(Field('dispVF', data=v_displacement_f,
                         lon=fieldset.U.grid.lon, lat=fieldset.U.grid.lat,
                         mesh='spherical'))
fieldset.dispUF.units = GeographicPolar()
fieldset.dispVF.units = Geographic()
fieldset.add_field(Field('landmask_fine', landmask_fine,
                         lon=fieldset.U.grid.lon, lat=fieldset.U.grid.lat,
                         mesh='spherical'))
fieldset.add_field(Field('distance2shore_fine', d_2_s_f,
                         lon=fieldset.U.grid.lon, lat=fieldset.U.grid.lat,
                         mesh='spherical'))
class DisplacementParticle(JITParticle):
    dU = Variable('dU')
    dV = Variable('dV')
    d2s = Variable('d2s', initial=1e3)
    age = Variable('age', dtype=np.float32, initial=0.)
    cycle_phase=Variable('cycle_phase', dtype=np.float32, initial=0.)
    releaseSite = Variable('releaseSite', dtype=np.int32)
    distance = Variable('distance', dtype=np.int32, initial=0.) # not calculating distance for now but left this in
    prev_lat = Variable('prev_lat', initial=0., to_write=False)  
    prev_lon = Variable('prev_lon', initial=0., to_write=False)
    f = Variable('f', dtype=np.int32)

#


# In[30]:

# In[87]:


kh = 10.0   # This is the eddy diffusivity in m2/s


# Add even diffusivity to the fieldset; value of constant field is kh 

# In[88]:


fieldset.add_constant_field('Kh_zonal', kh, mesh='spherical')
#zonal follows lat
fieldset.add_constant_field('Kh_meridional', kh, mesh='spherical') 


# In[ ]:

# In[89]:


class AgeParticle(JITParticle):
    age = Variable('age', dtype=np.float32, initial=0.)
    cycle_phase=Variable('cycle_phase', dtype=np.float32, initial=0.)


# In[32]:

# In[90]:


source_loc = pd.read_csv('/home/jsuca/Uku_Dispersal_Code/Uku_General_Habitat_PB.csv', header=None)
# Number of particle released per location
npart = 3
# Release location from the file read in above
lon = np.repeat(source_loc[1], npart)
lat = np.repeat(source_loc[2],npart)
# Start date for release. Since we are releasing every set number of days the repeatdt version was simplest
#start_date = 0
dlayer = [0.25]*(len(source_loc)*npart)
repeatdt = timedelta(days=2)


# Start depth

# In[91]:


pset = ParticleSet.from_list(fieldset=fieldset, pclass=DisplacementParticle,lon=lon,lat=lat,
                             depth=dlayer, time=1245600, repeatdt=repeatdt)#. ,#this time adjusts it to start on June 1


# In[23]:

# In[92]:





# In[93]:


def DeleteParticle(particle, fieldset, time):
    particle.delete()


# In[24]:

# In[94]:


output_file = pset.ParticleFile(name="PB_Uku_Roms_2014_20m_Full_Diff_bounce.nc", outputdt=timedelta(hours=12))


# In[21]:

# In[95]:


kernels=EggHatchingMovement+pset.Kernel(displace)+pset.Kernel(AdvectionRK4)+ pset.Kernel(DiffusionUniformKh)+ pset.Kernel(set_displacement)# + pset.Kernel(smagdiff)
run_days=168
model_dt=timedelta(minutes=10)


# set.set_variable_write_status('depth', False)

# Execute and release daily during this timeframe

# In[96]:


pset.execute(kernels,
            runtime=timedelta(days=run_days),
            dt=model_dt, 
            recovery={ErrorCode.ErrorOutOfBounds: DeleteParticle},
            output_file=output_file)


# now stop the repeated release

# In[97]:


pset.repeatdt = None


# now continue running for the remaining length of the PLD

# In[98]:


pset.execute(kernels,
            runtime=timedelta(days=45+1),
            dt=model_dt, 
            recovery={ErrorCode.ErrorOutOfBounds: DeleteParticle},
            output_file=output_file)


# In[ ]:

# In[99]:


output_file.export()


# In[ ]:
