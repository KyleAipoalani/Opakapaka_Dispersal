#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os
os.environ["PARCELS_COMPILER"]="gcc"
os.environ["CC"]="gcc"



# In[ ]:


# ============================================================
# MASTER IMPORT CELL (Run this first)
# ============================================================

# -------------------------
# Core Python
# -------------------------
import os
import sys

# -------------------------
# Data handling
# -------------------------
import numpy as np
import pandas as pd
import xarray as xr
import zarr

# -------------------------
# Plotting
# -------------------------
import matplotlib.pyplot as plt

# -------------------------
# Geometry / spatial
# -------------------------
from shapely.geometry import Point, Polygon

# -------------------------
# Bathymetry / raster
# -------------------------
import rasterio

# -------------------------
# Projection (if used)
# -------------------------
from pyproj import Transformer

# -------------------------
# OceanParcels
# -------------------------
from parcels import (
    FieldSet,
    ParticleSet,
    JITParticle,
    AdvectionRK4,
    Field,
)

# -------------------------
# Time handling
# -------------------------
from datetime import timedelta

# ============================================================
# QUICK CHECK (optional but helpful)
# ============================================================

print("All imports loaded successfully ✅")

# Check Parcels version
import parcels
print("Parcels version:", parcels.__version__)


# In[2]:


import sys
import parcels
import numpy

print("Python path:", sys.executable)
print("Parcels version:", parcels.__version__)
print("NumPy version:", numpy.__version__)


# In[3]:


import parcels
print(parcels.__version__)


# In[4]:


import glob, numpy as np, xarray as xr, matplotlib.pyplot as plt
from datetime import timedelta

from parcels import (
    FieldSet, ParticleSet, ScipyParticle, JITParticle, AdvectionRK4,
    ParticleFile, Field, Variable, DiffusionUniformKh, GeographicPolar, Geographic
)

from netCDF4 import Dataset
import pandas as pd


# In[5]:


import glob
import os

# ============================================================
# YOUR FOLDER PATHS (EDIT THESE IF NEEDED)  this is what it should be named( yr_2009_06_depth_30m.nc)

# ============================================================

folders = [
    r"/Users/justin.suca/Documents/TPruitt/ROMS(0-50)/2009/0.25_10_20_30_40_50_60m",
    r"/Users/justin.suca/Documents/TPruitt/ROMS(0-50)/2010/0.25_10_20_30_40_50_60m"

]

# choose multiple depth here
depths = ["30m"]   # change as needed

# ============================================================
# AUTO FIND FILES FROM ALL FOLDERS
# ============================================================

all_files = []

for folder in folders:
    files_in_folder = glob.glob(os.path.join(folder, "*.nc"))
    print(f"\nScanning folder: {folder}")
    print("Files found in folder:", len(files_in_folder))

    all_files.extend(files_in_folder)

# ============================================================
# FILTER BY DEPTH
# ============================================================

roms_files_sorted = [
    f for f in all_files
    if any(f"_depth_{d}" in os.path.basename(f) for d in depths)
]

# sort files
roms_files_sorted = sorted(roms_files_sorted)

# ============================================================
# CHECK RESULTS
# ============================================================

print("\n====================================================")
print("FINAL RESULTS")
print("====================================================")
print("Folders searched:", len(folders))
print("Total .nc files found:", len(all_files))
print("Files matching depths:", len(roms_files_sorted))
print("Depths selected:", depths)

# count how many per depth (helpful debug)
depth_counts = {d: 0 for d in depths}
for f in roms_files_sorted:
    for d in depths:
        if f"_depth_{d}" in f:
            depth_counts[d] += 1

print("\nFiles per depth:")
for d, count in depth_counts.items():
    print(f"{d}: {count}")

print("\nFiles being used:")
for f in roms_files_sorted:
    print(f)


# In[6]:


import glob
import os
import xarray as xr

# ============================================================
# CHANGE THIS (your folders)
# ============================================================

folders = [
    r"/Users/justin.suca/Documents/TPruitt/ROMS(0-50)/2009/0.25_10_20_30_40_50_60m",
    r"/Users/justin.suca/Documents/TPruitt/ROMS(0-50)/2010/0.25_10_20_30_40_50_60m"
]

# ============================================================
# SCAN FILES
# ============================================================

all_files = []

for folder in folders:
    all_files.extend(glob.glob(os.path.join(folder, "*.nc")))

print("Total files found:", len(all_files))

# ============================================================
# CHECK DEPTH INSIDE FILE
# ============================================================

print("\nChecking depth info...\n")

for f in all_files:
    try:
        ds = xr.open_dataset(f)

        print("File:", os.path.basename(f))

        # check common depth variable names
        found = False

        for var in ["depth", "Depth", "z", "Z", "s_rho"]:
            if var in ds:
                print(f"  Found variable '{var}':", ds[var].values[:5])
                found = True

        # check attributes too
        for attr in ds.attrs:
            if "depth" in attr.lower():
                print(f"  Found attribute '{attr}':", ds.attrs[attr])
                found = True

        if not found:
            print("  ⚠️ No obvious depth info found")

        print("--------------------------------------------------")

        ds.close()

    except Exception as e:
        print("Error reading:", f)
        print(e)
        print("--------------------------------------------------")


# In[7]:


from parcels import FieldSet

# Make sure this exists (fails fast with a clear message)
assert "roms_files_sorted" in globals(), "roms_files_sorted is not defined. Run the file-list cell first."
assert len(roms_files_sorted) > 0, "roms_files_sorted is empty. Check your file paths."

# Use the first file as the grid/coord source (works for your monthly files)
roms_file = roms_files_sorted[0]

variables = {"U": "u", "V": "v"}
dimensions = {
    "U": {"lon": "longitude", "lat": "latitude", "time": "time", "depth": "depth"},
    "V": {"lon": "longitude", "lat": "latitude", "time": "time", "depth": "depth"},
}

filenames = {
    "U": {"lon": roms_file, "lat": roms_file, "depth": roms_file, "data": roms_files_sorted, "time": roms_files_sorted},
    "V": {"lon": roms_file, "lat": roms_file, "depth": roms_file, "data": roms_files_sorted, "time": roms_files_sorted},
}

fieldset = FieldSet.from_netcdf(
    filenames, variables, dimensions,
    mesh="spherical",
    allow_time_extrapolation=False
)

print("Depth levels:", fieldset.U.depth)
print("Time start:", fieldset.U.grid.time[0])
print("Time end:", fieldset.U.grid.time[-1])
print("Using ROMS files:")
for f in roms_files_sorted:
    print("  ", f)


# In[8]:


TOP_LEVEL = 30.0
fieldset.add_constant("TOP_LEVEL", TOP_LEVEL)
print("Running simulation at depth:", fieldset.TOP_LEVEL, "meters")


# In[10]:


import xarray as xr
import numpy as np

def all_times_from_files(files, time_name="time"):
    times = []
    for fp in files:
        with xr.open_dataset(fp, decode_times=True) as ds:
            t = ds[time_name].values
            t = np.atleast_1d(t)
            times.append(t)
    return np.concatenate(times)

all_times = all_times_from_files(roms_files_sorted)


print("First timestamp:", all_times[0])
print("Last timestamp:", all_times[-1])


# In[11]:


import numpy as np

# Parcels time is in SECONDS since a time origin
tsec = np.asarray(fieldset.U.grid.time, dtype=float)

# --- robustly extract a usable datetime64 origin ---
origin_obj = fieldset.U.grid.time_origin
print("time_origin raw type:", type(origin_obj))
print("time_origin raw value:", origin_obj)

# Try common Parcels shapes
if isinstance(origin_obj, np.datetime64):
    origin64 = origin_obj
elif isinstance(origin_obj, str):
    origin64 = np.datetime64(origin_obj)
elif hasattr(origin_obj, "origin"):          # e.g., TimeConverter(origin=...)
    origin64 = np.datetime64(origin_obj.origin)
elif hasattr(origin_obj, "time_origin"):     # alternative attribute name
    origin64 = np.datetime64(origin_obj.time_origin)
else:
    # last resort (works if it's a python datetime)
    origin64 = np.datetime64(origin_obj)

# Convert seconds -> datetime64 so we can mask by dates
tsec_int = np.round(tsec).astype("int64")  # safer than float -> timedelta directly
tdt = origin64 + tsec_int.astype("timedelta64[s]")

YEAR = 2009
start_target = np.datetime64(f"{YEAR}-06-01T00:00:00")
end_target   = np.datetime64(f"{YEAR}-10-01T00:00:00")  # June..Sep (end is Oct 1)

mask_win = (tdt >= start_target) & (tdt < end_target)

print("Times in window:", int(mask_win.sum()))
print("Window start dt:", tdt[mask_win][0] if mask_win.any() else None)
print("Window end dt  :", tdt[mask_win][-1] if mask_win.any() else None)

if not mask_win.any():
    raise ValueError("No timestamps found in June–Sep window. (Likely your FieldSet only includes 1 month.)")

# IMPORTANT: start_time/end_time must be SECONDS (for Parcels)
start_time = float(tsec[mask_win][0])
end_time   = float(tsec[mask_win][-1])

print("Chosen start_time (sec):", start_time)
print("Chosen end_time   (sec):", end_time)
print(" number of days/runtime days:", (end_time - start_time) / 86400.0)


# In[12]:


import numpy as np

t = fieldset.U.grid.time
print("t0 sec:", t[0], "tN sec:", t[-1], "days:", (t[-1]-t[0])/86400)

# If you have the datetime array you used earlier:
print("min dt:", all_times.min())
print("max dt:", all_times.max())
print("unique months:", np.unique(all_times.astype("datetime64[M]")))


# In[ ]:





# In[13]:


AGG_FACTOR = 5          # "convert resolution by factor of 5"
PR_THRESH  = 0.25       # below this -> 0 particles
WEIGHT_MULT = 10 
Kh = 10.0  # horizontal diffusion [m^2/s]
# particles = round(PR * 10)


# In[14]:


# ============================================================
# LOAD BATHYMETRY (ETOPO)
# ============================================================

import rioxarray as rxr

bathy_path = r"/Users/justin.suca/Documents/TPruitt/ROMS(0-50)/ETOPO_2022 (Bedrock; 15 arcseconds).tiff"

bathy = rxr.open_rasterio(bathy_path)

depth = bathy.values[0]
lon_bathy = bathy.x.values
lat_bathy = bathy.y.values

print("Bathymetry loaded")
print("Depth grid shape:", depth.shape)


# In[15]:


from scipy.ndimage import binary_dilation
import numpy as np

# select 40–60 m depth band
habitat_band = (depth <= -40) & (depth >= -60)

# 15 arc-second bathymetry is about 0.463 km per cell
km_per_cell = 0.463

# use ceil so the buffer is at least 4 km
buffer_cells = int(np.ceil(6/ km_per_cell))

print("km_per_cell =", km_per_cell)
print("buffer_cells =", buffer_cells)
print("actual buffer ≈", buffer_cells * km_per_cell, "km")

# expand habitat mask
habitat_buffered = binary_dilation(habitat_band, iterations=buffer_cells)

habitat_mask = habitat_buffered.astype(np.int32)

print("Settlement habitat mask created")
print("Habitat cells:", habitat_mask.sum())


# In[16]:


fieldset = FieldSet.from_netcdf(filenames, variables, dimensions,interp_method={'U': 'freeslip', 'V': 'freeslip'})

# In[12]:

file_path_fine = r"/Users/justin.suca/Documents/TPruitt/ROMS/2009/0.25_1_5_10_20_30_50m"
lon_grid = np.asarray(fieldset.U.grid.lon)
lat_grid = np.asarray(fieldset.U.grid.lat)

if lon_grid.ndim == 1 and lat_grid.ndim == 1:
    # lon[x], lat[y]  -> data must be [t, y, x]
    Kh_z = np.full((1, lat_grid.size, lon_grid.size), Kh, dtype=np.float32)
    Kh_m = np.full((1, lat_grid.size, lon_grid.size), Kh, dtype=np.float32)

    fieldset.add_field(Field("Kh_zonal", Kh_z, lon=lon_grid, lat=lat_grid, mesh="spherical"))
    fieldset.add_field(Field("Kh_meridional", Kh_m, lon=lon_grid, lat=lat_grid, mesh="spherical"))

else:
    # lon[y,x], lat[y,x] -> data must be [t, y, x] matching that grid
    Kh_z = np.full((1,) + lon_grid.shape, Kh, dtype=np.float32)
    Kh_m = np.full((1,) + lon_grid.shape, Kh, dtype=np.float32)

    fieldset.add_field(Field("Kh_zonal", Kh_z, lon=lon_grid, lat=lat_grid, mesh="spherical", transpose=True))
    fieldset.add_field(Field("Kh_meridional", Kh_m, lon=lon_grid, lat=lat_grid, mesh="spherical", transpose=True))


file_path_fine = roms_files_sorted[0]   # uses the same year/month/depth file
print("USING file_path_fine:", file_path_fine)


from parcels import Field

habitat_field = Field(
    name="habitat",
    data=habitat_mask,
    lon=lon_bathy,
    lat=lat_bathy
)

fieldset.add_field(habitat_field)

print("Habitat field added to fieldset")

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
    """Clamp to safe inner bounds before sampling static fields."""
    # use SAFE_* + EPS to stay one full cell in
    if particle.lon <= fieldset.SAFE_LON_MIN:
        particle.lon = fieldset.SAFE_LON_MIN + fieldset.LON_EPS
    if particle.lon >= fieldset.SAFE_LON_MAX:
        particle.lon = fieldset.SAFE_LON_MAX - fieldset.LON_EPS
    if particle.lat <= fieldset.SAFE_LAT_MIN:
        particle.lat = fieldset.SAFE_LAT_MIN + fieldset.LAT_EPS
    if particle.lat >= fieldset.SAFE_LAT_MAX:
        particle.lat = fieldset.SAFE_LAT_MAX - fieldset.LAT_EPS

    particle.d2s = fieldset.distance2shore_fine[time, particle.depth, particle.lat, particle.lon]

    if particle.d2s < fieldset.shore_threshold:
        particle.dU = fieldset.dispUF[time, particle.depth, particle.lat, particle.lon]
        particle.dV = fieldset.dispVF[time, particle.depth, particle.lat, particle.lon]
    else:
        particle.dU = 0.0
        particle.dV = 0.0



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
from parcels import JITParticle, Variable



#look at this warning


# In[17]:


print("depth exists:", "depth" in globals())
print("lon_bathy exists:", "lon_bathy" in globals())
print("lat_bathy exists:", "lat_bathy" in globals())
print("habitat_mask exists:", "habitat_mask" in globals())


# In[18]:


# ============================================================
# DEFINE HAWAII ISLAND POLYGONS
# ============================================================

from shapely.geometry import Polygon

island_polygons = {

    "Kauai": Polygon([
        (-159.70, 22.30),
        (-159.20, 22.30),
        (-159.20, 21.85),
        (-159.70, 21.85)
    ]),

    "Niihau": Polygon([
        (-160.30, 22.05),
        (-159.95, 22.05),
        (-159.95, 21.75),
        (-160.30, 21.75)
    ]),

    "Kaula": Polygon([
        (-160.75, 21.72),
        (-160.45, 21.72),
        (-160.45, 21.50),
        (-160.75, 21.50)
    ]),

    "Oahu": Polygon([
        (-158.40, 21.80),
        (-157.60, 21.80),
        (-157.60, 21.10),
        (-158.40, 21.10)
    ]),

    "Penguin_Bank": Polygon([
        (-157.85, 21.15),
        (-157.20, 21.15),
        (-157.20, 20.80),
        (-157.85, 20.80)
    ]),

    "Maui_Nui": Polygon([
        (-157.4, 21.2),
        (-156.0, 21.2),
        (-156.0, 20.4),
        (-157.4, 20.4)
    ]),

    "Hawaii": Polygon([
        (-156.0, 20.3),
        (-154.7, 20.3),
        (-154.7, 18.9),
        (-156.0, 18.9)
    ]),
}

print("Island polygons created:", list(island_polygons.keys()))
for a_name, a_poly in island_polygons.items():
    for b_name, b_poly in island_polygons.items():
        if a_name < b_name:
            if a_poly.intersects(b_poly):
                print(f"Overlap: {a_name} with {b_name}")


# In[19]:


from shapely.geometry import Point

region_masks = {}

LON2, LAT2 = np.meshgrid(lon_bathy, lat_bathy)

for region_name, poly in island_polygons.items():
    inside_poly = np.zeros_like(habitat_mask, dtype=np.int32)

    for j in range(LAT2.shape[0]):
        for i in range(LON2.shape[1]):
            p = Point(LON2[j, i], LAT2[j, i])
            if p.within(poly):
                inside_poly[j, i] = 1

    # combine polygon + habitat
    region_mask = ((habitat_mask == 1) & (inside_poly == 1)).astype(np.int32)

    region_masks[region_name] = region_mask
    print(region_name, "cells:", region_mask.sum())


# In[20]:


# ============================================================
# REMOVE PENGUIN BANK FROM OAHU AND MAUI_NUI MASKS
# ============================================================

pb = region_masks["Penguin_Bank"] == 1

region_masks["Oahu"][pb] = 0
region_masks["Maui_Nui"][pb] = 0

print("Removed Penguin_Bank cells from Oahu and Maui_Nui masks.")

# quick overlap check at MASK level
for a_name, a_mask in region_masks.items():
    for b_name, b_mask in region_masks.items():
        if a_name < b_name:
            overlap_cells = np.sum((a_mask == 1) & (b_mask == 1))
            if overlap_cells > 0:
                print(f"Mask overlap: {a_name} with {b_name} = {overlap_cells} cells")


# In[21]:


# ============================================================
# CELL 20A
# CONVERT STRICT REGION MASKS -> STRICT HABITAT POLYGONS
# ============================================================

from rasterio.features import shapes
from shapely.geometry import shape, MultiPolygon, Polygon
from shapely.ops import unary_union
from shapely.prepared import prep
import rasterio

strict_region_polygons = {}
strict_region_polygons_prepared = {}

# use the bathy raster transform so polygons line up with the grid
transform = bathy.rio.transform()

for region_name, region_mask in region_masks.items():
    geoms = []

    # extract polygons where mask == 1
    for geom, value in shapes(region_mask.astype(np.int16), mask=(region_mask == 1), transform=transform):
        if value == 1:
            geoms.append(shape(geom))

    if len(geoms) == 0:
        strict_region_polygons[region_name] = None
        strict_region_polygons_prepared[region_name] = None
        print(region_name, "-> no polygon created")
        continue

    merged = unary_union(geoms)

    # clean small geometry issues
    merged = merged.buffer(0)

    strict_region_polygons[region_name] = merged
    strict_region_polygons_prepared[region_name] = prep(merged)

    print(region_name, "polygon area created")

print("\nStrict habitat polygons created from region_masks.")


# In[22]:


from parcels import Field

for region_name, region_mask in region_masks.items():
    field_name = f"mask_{region_name}"

    mask_3d = region_mask[np.newaxis, :, :]

    fieldset.add_field(Field(
        name=field_name,
        data=mask_3d,
        lon=lon_bathy,
        lat=lat_bathy,
        mesh="spherical"
    ))

    print("Added:", field_name)


# In[23]:


print(hasattr(fieldset, "mask_Kaula"))
print(hasattr(fieldset, "mask_Kauai"))
print(hasattr(fieldset, "mask_Niihau"))
print([name for name in dir(fieldset) if name.startswith("mask_")])


# In[24]:


# ============================================================
# INITIALIZE SETTLEMENT COUNTER
# ============================================================

fieldset.add_constant("settlement_counter", 0)

print("Settlement counter initialized")


# In[25]:


# Quick domain bounds 
import numpy as np
LON_MIN, LON_MAX = float(np.min(fieldset.U.grid.lon)), float(np.max(fieldset.U.grid.lon))
LAT_MIN, LAT_MAX = float(np.min(fieldset.U.grid.lat)), float(np.max(fieldset.U.grid.lat))
print("Domain lon:", LON_MIN, "→", LON_MAX, " | lat:", LAT_MIN, "→", LAT_MAX)

# Detect available depth levels 
zlevels = np.asarray(fieldset.U.depth if hasattr(fieldset.U, "depth") else fieldset.U.grid.depth, dtype=float)
zlevels_sorted = np.sort(zlevels)
print("Depth levels:", zlevels_sorted)
TOP_LEVEL = float(zlevels_sorted[0])  # e.g., 50.0 for your 60–100 m data


# In[26]:


print("Kh_zonal:", fieldset.Kh_zonal)
print("Kh_meridional:", fieldset.Kh_meridional)


# In[ ]:





# In[27]:


# --- depth info from ROMS and safe constants ---

# Get all the depth levels that the ocean model uses.
# Think of this as a list of water depths where the model can place particles.
zlevels = np.asarray(getattr(fieldset.U, "depth", fieldset.U.grid.depth), dtype=float)

# Sort those depth values from shallowest to deepest.
# This makes it easy to grab the top (shallowest) and bottom (deepest) depths.
zlevels_sorted = np.sort(zlevels)

# Save the shallowest depth from the model.
# We use this so we know how close to the "surface" we’re allowed to start.
DEPTH_MIN = float(zlevels_sorted[0])  # e.g., 60.0

# Save the deepest depth from the model.
# This tells us how far down the model goes, so we don’t try to go past it.
DEPTH_MAX = float(zlevels_sorted[-1])  # e.g., 100.0

# Pick the shallowest depth as the depth where things will drift.
# This is used to put eggs/particles as close to the surface as the model allows.
DRIFT_DEPTH = DEPTH_MIN

# Store these values inside the fieldset so all parts of the code can use them.
# This keeps everything consistent and avoids hard-coding numbers in many places.
fieldset.add_constant("DEPTH_MIN",   DEPTH_MIN)
fieldset.add_constant("DEPTH_MAX",   DEPTH_MAX)
fieldset.add_constant("DRIFT_DEPTH", DRIFT_DEPTH)

# Use the drift depth as the "top level" for particles.
# This means eggs start at a depth that is valid for the model instead of an
# unrealistic depth that could crash the run.
TOP_LEVEL = DRIFT_DEPTH
fieldset.add_constant("TOP_LEVEL", TOP_LEVEL)


# In[28]:


# --- Edge nudges so particles don't sit exactly on the boundary ---
import numpy as np

# Grab all the longitude (left–right) positions from the model grid
LON_arr = np.asarray(fieldset.U.grid.lon, dtype=float)

# Grab all the latitude (up–down) positions from the model grid
LAT_arr = np.asarray(fieldset.U.grid.lat, dtype=float)

# Estimate a "typical" spacing between longitudes.
# This tells us roughly how far apart the grid points are left–to–right.
LON_DX = float(np.nanmedian(np.abs(np.diff(LON_arr))))

# Estimate a "typical" spacing between latitudes.
# This tells us roughly how far apart the grid points are up–down.
LAT_DY = float(np.nanmedian(np.abs(np.diff(LAT_arr))))

# Store a small nudge distance in longitude (1/4 of a grid step).
# We use this to push particles slightly away from the very edge,
# so they don't sit exactly on the border and cause "out of bounds" problems.
fieldset.add_constant("LON_EPS", 0.25 * LON_DX)

# Store a small nudge distance in latitude (1/4 of a grid step).
# Same idea: gently move particles away from the top/bottom borders
# so the model runs more safely.
fieldset.add_constant("LAT_EPS", 0.25 * LAT_DY)

# Find the shallowest depth level in the model.
# We save this so we know the highest (closest to surface) depth we can use.
fieldset.add_constant(
    "DEPTH_MIN",
    float(np.min(getattr(fieldset.U, "depth", fieldset.U.grid.depth)))
)

# Find the deepest depth level in the model.
# We save this so we know the lowest depth we can go to without breaking the model.
fieldset.add_constant(
    "DEPTH_MAX",
    float(np.max(getattr(fieldset.U, "depth", fieldset.U.grid.depth)))
)


# In[29]:


# --- Make a clean list of Opakapaka release locations that are only in the ocean ---

import numpy as np, pandas as pd, xarray as xr

# A) Load the habitat CSV file
# This file lists possible release spots for Opakapaka, with their longitude, latitude, and weights.
release_csv = r"/Users/justin.suca/Documents/TPruitt/ROMS(0-50)/Opakapaka_General_Habitat_PB.csv"
df = pd.read_csv(release_csv, header=None, names=["id", "lon", "lat", "weight"])

# Turn the longitude and latitude columns into regular arrays so we can work with them easily.
release_lons_raw = df["lon"].to_numpy(dtype=float)
release_lats_raw = df["lat"].to_numpy(dtype=float)

# Print how many total points we started with.
print(f"Loaded {release_lons_raw.size:,} candidate release points from CSV")

# B) Remove any points that fall outside the ocean model area (the ROMS grid)
# This makes sure we only keep points inside the part of the ocean your model covers.
LON_MIN, LON_MAX = float(np.min(fieldset.U.grid.lon)), float(np.max(fieldset.U.grid.lon))
LAT_MIN, LAT_MAX = float(np.min(fieldset.U.grid.lat)), float(np.max(fieldset.U.grid.lat))

# Check which points fall inside the model’s box.
in_box = (
    (release_lons_raw >= LON_MIN) & (release_lons_raw <= LON_MAX) &
    (release_lats_raw >= LAT_MIN) & (release_lats_raw <= LAT_MAX)
)

# Keep only those that are inside.
release_lons_box = release_lons_raw[in_box]
release_lats_box = release_lats_raw[in_box]
print(f"After domain box filter: {release_lons_box.size:,} points")

# C) Use the model data to tell which points are on land and which are in the ocean.
# The ROMS file stores ocean current data. If a point’s value is NaN (not a number),
# that means it’s on land, so we skip those.
roms_first = roms_files_sorted[0]
with xr.open_dataset(roms_first, decode_times=True) as ds:
    # Take a single time and depth layer from the U (east-west current) data.
    u0 = ds["u"].isel(time=0, depth=0).load()

    # Find the value closest to each release point.
    sampled = u0.interp(
        longitude=("points", release_lons_box),
        latitude=("points", release_lats_box),
        method="nearest"
    ).values

# Keep only the points where the data exists (those are the ocean points).
keep_ocean = np.isfinite(sampled)
release_lons = release_lons_box[keep_ocean]
release_lats = release_lats_box[keep_ocean]

# Print how many valid ocean release points are left.
print(f"After ocean mask: {release_lons.size:,} points kept (ocean only)")
# make a note where the land particles went to


# In[30]:


zlevels = np.asarray(getattr(fieldset.U, "depth", fieldset.U.grid.depth), dtype=float)
zlevels_sorted = np.sort(zlevels)
fieldset.add_constant('DEPTH_MIN', float(zlevels_sorted[0]))
fieldset.add_constant('DEPTH_MAX', float(zlevels_sorted[-1]))
fieldset.add_constant('TOP_LEVEL', float(zlevels_sorted[0]))


# In[31]:


fieldset.add_constant('TOP_LEVEL', float(TOP_LEVEL))


# In[32]:


release_csv = r"/Users/justin.suca/Documents/TPruitt/ROMS(0-50)/Opakapaka_General_Habitat_PB.csv"
release_df = pd.read_csv(release_csv, header=None, names=["id","lon","lat","weight"])
print("Total sites:", len(release_df))


# In[33]:


import numpy as np
import pandas as pd
from pyproj import Transformer

# --- Load habitat CSV ---
release_csv = r"/Users/justin.suca/Documents/TPruitt/ROMS(0-50)/Opakapaka_General_Habitat_PB.csv"

df = pd.read_csv(release_csv, header=None,
                 names=["id", "lon", "lat", "weight"])

# --- Keep only good habitat points ---
df = df[df["weight"] > PR_THRESH].copy()
print("After threshold:", len(df), "points")

# ============================================================
# ✅ BIN IN METERS (UTM) INSTEAD OF DEGREES
# ============================================================

BIN_M = 250  # <-- 250 meter bins 

# 1) Convert lon/lat (degrees) -> UTM meters (Zone 4N for Hawaii)
to_utm = Transformer.from_crs("EPSG:4326", "EPSG:32604", always_xy=True)

df["x_m"], df["y_m"] = to_utm.transform(
    df["lon"].to_numpy(),
    df["lat"].to_numpy()
)

# 2) Bin into 250m x 250m grid squares
df["x_bin"] = np.floor(df["x_m"] / BIN_M).astype(int)
df["y_bin"] = np.floor(df["y_m"] / BIN_M).astype(int)

agg = df.groupby(["x_bin", "y_bin"], as_index=False).agg(
    x_m=("x_m", "mean"),
    y_m=("y_m", "mean"),
    weight=("weight", "mean")  
)


print("After aggregation:", len(agg), "release cells")

# 4) Convert bin centers back -> lon/lat for Parcels
to_ll = Transformer.from_crs("EPSG:32604", "EPSG:4326", always_xy=True)

agg["lon"], agg["lat"] = to_ll.transform(
    agg["x_m"].to_numpy(),
    agg["y_m"].to_numpy()
)

# 5) Convert habitat weight -> particle counts
agg["n_particles"] = (agg["weight"] * WEIGHT_MULT).round().astype(int)
agg = agg[agg["n_particles"] > 0].copy()

print("Final release sites:", len(agg))
print("Total particles:", int(agg["n_particles"].sum()))

# 6) Final release arrays
release_lons = agg["lon"].to_numpy(float)
release_lats = agg["lat"].to_numpy(float)
release_counts = agg["n_particles"].to_numpy(int)

print("Spawn points created:", int(release_counts.sum()))


# In[34]:


lon = np.repeat(release_lons, release_counts)
lat = np.repeat(release_lats, release_counts)

print("Spawn points created:", lon.size)


# In[35]:


print("release sites:", len(release_lons))
print("total particles:", int(release_counts.sum()))
print("spawn points created:", lon.size)
print("min/max weight:", float(agg["weight"].min()), float(agg["weight"].max()))
print("min/max n_particles:", int(agg["n_particles"].min()), int(agg["n_particles"].max()))


# In[36]:


# --- Set up constants and settings for the model  ---
import numpy as np

# 0) Make sure we have release points ready
# This checks if the variables holding your release locations exist.
# If they don’t, it gives an error reminding you to run the earlier setup cell.
if 'release_lons' in globals() and 'release_lats' in globals():
    lon = np.asarray(release_lons, dtype=float)
    lat = np.asarray(release_lats, dtype=float)
elif 'lon' in globals() and 'lat' in globals():
    lon = np.asarray(lon, dtype=float)
    lat = np.asarray(lat, dtype=float)
else:
    raise NameError("No release positions found. Run the habitat/ocean-mask selection cell first.")

# 1) Get the domain limits (edges) of the model
# These values mark the minimum and maximum longitudes and latitudes
# in the ROMS grid — basically, the edges of your ocean model area.
LON_MIN, LON_MAX = float(np.min(fieldset.U.grid.lon)), float(np.max(fieldset.U.grid.lon))
LAT_MIN, LAT_MAX = float(np.min(fieldset.U.grid.lat)), float(np.max(fieldset.U.grid.lat))

# 2) Make these limits available to the particles
# This allows your particles to “see” the boundaries
# so kernels can keep them inside the domain if needed.
fieldset.add_constant('LON_MIN', LON_MIN)
fieldset.add_constant('LON_MAX', LON_MAX)
fieldset.add_constant('LAT_MIN', LAT_MIN)
fieldset.add_constant('LAT_MAX', LAT_MAX)

# 3) Define a shoreline buffer distance
# This sets how close a particle can get to the coast before being considered “too close.”
# The number 5.0 here means 5 grid cells from shore 
fieldset.add_constant('shore_threshold', 5.0)

# 4) Pick the top layer of the model for particle release depth
# The ROMS data includes multiple depths (e.g., 0m, 10m, 20m, etc.).
# This finds the shallowest one — where your Opakapaka eggs or larvae start.
try:
    TOP_LEVEL = float(zlevels_sorted[0])  # use the first (top) depth level if already sorted
except NameError:
    # if not yet defined, get the depth levels directly from the model and sort them
    zlevels = np.asarray(getattr(fieldset.U, "depth", fieldset.U.grid.depth), dtype=float)
    zlevels_sorted = np.sort(zlevels)
    TOP_LEVEL = float(zlevels_sorted[0])

# Print a short summary so you can confirm everything looks right
print(f"Release count: {lon.size} | TOP_LEVEL={TOP_LEVEL} m")
print(f"Domain lon:[{LON_MIN:.2f},{LON_MAX:.2f}] lat:[{LAT_MIN:.2f},{LAT_MAX:.2f}]")


# In[37]:


LON_MIN = float(np.min(fieldset.U.grid.lon))
LON_MAX = float(np.max(fieldset.U.grid.lon))
LAT_MIN = float(np.min(fieldset.U.grid.lat))
LAT_MAX = float(np.max(fieldset.U.grid.lat))

EPS = 0.02  # about ~2 km-ish buffer

SAFE_LON_MIN = LON_MIN + EPS
SAFE_LON_MAX = LON_MAX - EPS
SAFE_LAT_MIN = LAT_MIN + EPS
SAFE_LAT_MAX = LAT_MAX - EPS

fieldset.add_constant("LON_MIN", LON_MIN)
fieldset.add_constant("LON_MAX", LON_MAX)
fieldset.add_constant("LAT_MIN", LAT_MIN)
fieldset.add_constant("LAT_MAX", LAT_MAX)

fieldset.add_constant("SAFE_LON_MIN", SAFE_LON_MIN)
fieldset.add_constant("SAFE_LON_MAX", SAFE_LON_MAX)
fieldset.add_constant("SAFE_LAT_MIN", SAFE_LAT_MIN)
fieldset.add_constant("SAFE_LAT_MAX", SAFE_LAT_MAX)

fieldset.add_constant("LON_EPS", 1e-6)
fieldset.add_constant("LAT_EPS", 1e-6)

fieldset.add_constant("shore_threshold", 1.5)


# In[38]:


def DeleteOnError(particle, fieldset, time):
    particle.delete()


# In[39]:


# ===============================
# DEFINE PLD (settlement window)
# ===============================

PLD_MIN_DAYS = 60   # larvae must drift at least 60 days
PLD_MAX_DAYS = 180  # larvae stop settling after 180 days

fieldset.add_constant("PLD_MIN_SEC", PLD_MIN_DAYS * 86400)
fieldset.add_constant("PLD_MAX_SEC", PLD_MAX_DAYS * 86400)


# In[40]:


# ===============================
# MORTALITY SETTINGS
# ===============================

TARGET_MORTALITY = 0.10   # want about 10% dead overall by end of horizon
MORTALITY_HORIZON_DAYS = 180  # use 300 because your runtime is 300 days

fieldset.add_constant("TARGET_MORTALITY", TARGET_MORTALITY)
fieldset.add_constant("MORTALITY_HORIZON_SEC", MORTALITY_HORIZON_DAYS * 86400.0)

print("Mortality target:", TARGET_MORTALITY)
print("Mortality horizon (days):", MORTALITY_HORIZON_DAYS)


# In[41]:


from parcels import JITParticle, Variable
import numpy as np

class DisplacementParticle(JITParticle):
    age = Variable("age", dtype=np.float32, initial=0.0)
    d2s = Variable("d2s", dtype=np.float32, initial=9999.0)
    dU = Variable("dU", dtype=np.float32, initial=0.0)
    dV = Variable("dV", dtype=np.float32, initial=0.0)

    kill_reason = Variable("kill_reason", dtype=np.int32, initial=0)
    alive = Variable("alive", dtype=np.int32, initial=1)

    settled = Variable("settled", dtype=np.int32, initial=0)

    settle_lon = Variable("settle_lon", dtype=np.float32, initial=0.0)
    settle_lat = Variable("settle_lat", dtype=np.float32, initial=0.0)
    settle_time = Variable("settle_time", dtype=np.float32, initial=0.0)

    release_region = Variable("release_region", dtype=np.int32, initial=-1)
    settle_region = Variable("settle_region", dtype=np.int32, initial=-1)

    death_age = Variable("death_age", dtype=np.float32, initial=-1.0)

    release_group = Variable("release_group", dtype=np.int32, initial=-1)
    will_die = Variable("will_die", dtype=np.int32, initial=0)
    death_time_target = Variable("death_time_target", dtype=np.float32, initial=-1.0)
# ============================================================
# ADD SETTLEMENT VARIABLES TO PARTICLES
# ============================================================

from parcels import Variable
import numpy as np







def set_displacement(particle, fieldset, time):
    if particle.lon <= fieldset.SAFE_LON_MIN:
        particle.lon = fieldset.SAFE_LON_MIN + fieldset.LON_EPS
    if particle.lon >= fieldset.SAFE_LON_MAX:
        particle.lon = fieldset.SAFE_LON_MAX - fieldset.LON_EPS
    if particle.lat <= fieldset.SAFE_LAT_MIN:
        particle.lat = fieldset.SAFE_LAT_MIN + fieldset.LAT_EPS
    if particle.lat >= fieldset.SAFE_LAT_MAX:
        particle.lat = fieldset.SAFE_LAT_MAX - fieldset.LAT_EPS

    particle.d2s = fieldset.distance2shore_fine[time, particle.depth, particle.lat, particle.lon]

    if particle.d2s < fieldset.shore_threshold:
        particle.dU = fieldset.dispUF[time, particle.depth, particle.lat, particle.lon]
        particle.dV = fieldset.dispVF[time, particle.depth, particle.lat, particle.lon]
    else:
        particle.dU = 0.0
        particle.dV = 0.0


def displace(particle, fieldset, time):
    if particle.d2s < fieldset.shore_threshold:
        particle_dlon += particle.dU * particle.dt
        particle_dlat += particle.dV * particle.dt




def KillIfOutOfBounds(particle, fieldset, time):

    # Define safe domain
    if (
        particle.lon < fieldset.SAFE_LON_MIN or
        particle.lon > fieldset.SAFE_LON_MAX or
        particle.lat < fieldset.SAFE_LAT_MIN or
        particle.lat > fieldset.SAFE_LAT_MAX
    ):
        particle.delete()


def DeleteParticle(particle, fieldset, time):
    particle.delete()



def ClampDepth(particle, fieldset, time):

    if particle.depth < 0:
        particle.depth = 0.0

    if particle.depth > 200:
        particle.depth = 200.0

def KeepInDomain(particle, fieldset, time):

    eps = 1e-5   # small buffer inside the grid

    if particle.lon <= fieldset.SAFE_LON_MIN:
        particle_dlon += (fieldset.SAFE_LON_MIN + eps) - particle.lon

    elif particle.lon >= fieldset.SAFE_LON_MAX:
        particle_dlon += (fieldset.SAFE_LON_MAX - eps) - particle.lon

    if particle.lat <= fieldset.SAFE_LAT_MIN:
        particle_dlat += (fieldset.SAFE_LAT_MIN + eps) - particle.lat

    elif particle.lat >= fieldset.SAFE_LAT_MAX:
        particle_dlat += (fieldset.SAFE_LAT_MAX - eps) - particle.lat   
# ============================================================
# SETTLEMENT KERNEL
# ============================================================

# ============================================================

from parcels import ParcelsRandom
import math


def EggSurface_NoCap(particle, fieldset, time):
    particle.age += math.fabs(particle.dt)
    particle.depth = fieldset.TOP_LEVEL


from parcels import ParcelsRandom
import math

import math

def ApplyMortality(particle, fieldset, time):
    if particle.alive == 1:
        if particle.settled == 0:
            if particle.will_die == 1:
                if particle.age >= fieldset.PLD_MIN_SEC:
                    if particle.age <= fieldset.PLD_MAX_SEC:
                        if particle.death_time_target >= 0.0:
                            if particle.age >= particle.death_time_target:
                                particle.alive = 0
                                particle.kill_reason = 2
                                particle.death_age = particle.age


def CheckSettlement(particle, fieldset, time):
    if particle.alive == 1:
        if particle.settled == 0:
            if particle.age >= fieldset.PLD_MIN_SEC:
                if particle.age <= fieldset.PLD_MAX_SEC:

                    if fieldset.mask_Kaula[time, particle.depth, particle.lat, particle.lon] >= 0.5:
                        particle.settled = 1
                        particle.settle_region = 2

                    elif fieldset.mask_Niihau[time, particle.depth, particle.lat, particle.lon] >= 0.5:
                        particle.settled = 1
                        particle.settle_region = 1

                    elif fieldset.mask_Kauai[time, particle.depth, particle.lat, particle.lon] >= 0.5:
                        particle.settled = 1
                        particle.settle_region = 0

                    elif fieldset.mask_Oahu[time, particle.depth, particle.lat, particle.lon] >= 0.5:
                        particle.settled = 1
                        particle.settle_region = 3

                    elif fieldset.mask_Penguin_Bank[time, particle.depth, particle.lat, particle.lon] >= 0.5:
                        particle.settled = 1
                        particle.settle_region = 4

                    elif fieldset.mask_Maui_Nui[time, particle.depth, particle.lat, particle.lon] >= 0.5:
                        particle.settled = 1
                        particle.settle_region = 5

                    elif fieldset.mask_Hawaii[time, particle.depth, particle.lat, particle.lon] >= 0.5:
                        particle.settled = 1
                        particle.settle_region = 6

                    if particle.settled == 1:
                        particle.settle_time = time
                        particle.settle_lon = particle.lon
                        particle.settle_lat = particle.lat




def FreezeIfSettled(particle, fieldset, time):
    if particle.settled == 1:
        particle.lon = particle.settle_lon
        particle.lat = particle.settle_lat
        particle.depth = fieldset.TOP_LEVEL
        particle.dU = 0.0
        particle.dV = 0.0





def classify_region_from_lonlat(lon, lat):
    if island_polygons["Kaula"].contains(Point(lon, lat)):
        return 2
    elif island_polygons["Niihau"].contains(Point(lon, lat)):
        return 1
    elif island_polygons["Kauai"].contains(Point(lon, lat)):
        return 0
    elif island_polygons["Oahu"].contains(Point(lon, lat)):
        return 3
    elif island_polygons["Penguin_Bank"].contains(Point(lon, lat)):
        return 4
    elif island_polygons["Maui_Nui"].contains(Point(lon, lat)):
        return 5
    elif island_polygons["Hawaii"].contains(Point(lon, lat)):
        return 6
    else:
        return 7






# ======================================
# DEFINE ISLAND REGIONS (for connectivity)
# ======================================


def SetReleaseRegion(particle, fieldset, time):
    if particle.release_region == -1:
        particle.release_region = classify_region_from_lonlat(particle.lon, particle.lat)



# In[42]:


print([name for name in dir(fieldset) if name.startswith("mask_")])


# In[43]:


def build_kernels(pset):
    return (
        pset.Kernel(EggSurface_NoCap)
        + pset.Kernel(KeepInDomain)
        + pset.Kernel(AdvectionRK4)
        + pset.Kernel(set_displacement)
        + pset.Kernel(displace)
        + pset.Kernel(CheckSettlement)
        + pset.Kernel(FreezeIfSettled)
    )


# In[44]:


import numpy as np

rng = np.random.default_rng(42)

def apply_chunked_mortality(pset, target_fraction=0.10, pld_min_days=60, pld_max_days=180):
    particles = [p for p in pset]

    if len(particles) == 0:
        return

    alive = np.array([p.alive for p in particles], dtype=np.int32)
    settled = np.array([p.settled for p in particles], dtype=np.int32)
    age_days = np.array([p.age for p in particles], dtype=np.float64) / 86400.0
    kill_reason = np.array([p.kill_reason for p in particles], dtype=np.int32)
    release_group = np.array([p.release_group for p in particles], dtype=np.int32)

    unique_groups = np.unique(release_group)

    for g in unique_groups:
        group_idx = np.where(release_group == g)[0]
        total_g = len(group_idx)

        if total_g == 0:
            continue

        target_dead_final = int(round(target_fraction * total_g))

        dead_mask_g = (kill_reason[group_idx] == 2)
        current_dead = int(np.sum(dead_mask_g))

        # age should be basically the same within a release group
        current_day = float(np.max(age_days[group_idx]))

        if current_day < pld_min_days:
            target_dead_now = 0
        elif current_day >= pld_max_days:
            target_dead_now = target_dead_final
        else:
            frac = (current_day - pld_min_days) / (pld_max_days - pld_min_days)
            target_dead_now = int(round(target_dead_final * frac))

        need_to_kill = target_dead_now - current_dead

        if need_to_kill <= 0:
            continue

        eligible_mask = (
            (alive[group_idx] == 1) &
            (settled[group_idx] == 0) &
            (age_days[group_idx] >= pld_min_days) &
            (age_days[group_idx] <= pld_max_days)
        )

        eligible_local = np.where(eligible_mask)[0]

        if len(eligible_local) == 0:
            continue

        n_kill = min(need_to_kill, len(eligible_local))
        chosen_local = rng.choice(eligible_local, size=n_kill, replace=False)
        chosen_global = group_idx[chosen_local]

        for idx in chosen_global:
            particles[idx].alive = 0
            particles[idx].kill_reason = 2
            particles[idx].death_age = particles[idx].age

        print(
            f"Group {g}: day={current_day:.2f}, "
            f"dead_now={current_dead}, target_now={target_dead_now}, "
            f"killed_this_chunk={n_kill}"
        )


# In[45]:


from parcels import ParticleSet
import numpy as np

release_interval_days = 1   # 👈 changed to daily releases
n_release_periods = 120

release_interval_sec = release_interval_days * 86400

release_times = [
    start_time + release_interval_sec * i
    for i in range(n_release_periods)
]

print("Release interval (days):", release_interval_days)
print("Number of release periods:", n_release_periods)
print("Release times:")
for i, t in enumerate(release_times):
    print(i, t)

# ============================================================
# BUILD PARTICLES BY RELEASE COHORT
# ============================================================
all_lon = []
all_lat = []
all_time = []
all_release_group = []

group_id = 0

for rel_time in release_times:
    lon_group = np.repeat(release_lons, release_counts)
    lat_group = np.repeat(release_lats, release_counts)

    n_group = len(lon_group)

    all_lon.append(lon_group)
    all_lat.append(lat_group)
    all_time.append(np.full(n_group, rel_time))
    all_release_group.append(np.full(n_group, group_id, dtype=np.int32))

    print(f"Release group {group_id}: total={n_group}, release_time={rel_time}")

    group_id += 1

lon = np.concatenate(all_lon)
lat = np.concatenate(all_lat)
time_arr = np.concatenate(all_time)
release_group_arr = np.concatenate(all_release_group)

n_particles = len(lon)
print("\nTotal particles:", n_particles)
print("Number of release groups:", group_id)

# debug check
unique_groups, counts = np.unique(release_group_arr, return_counts=True)
print("\nParticles per release group:")
for g, c in zip(unique_groups, counts):
    print(f"Group {g}: {c}")

# ============================================================
# CREATE PARTICLE SET
# ============================================================
pset = ParticleSet.from_list(
    fieldset=fieldset,
    pclass=DisplacementParticle,
    lon=lon,
    lat=lat,
    depth=np.full(n_particles, float(fieldset.TOP_LEVEL), dtype=np.float32),
    time=time_arr,
    release_group=release_group_arr,
)


# In[46]:


fieldset.add_constant("shore_threshold", 5.0)
print("Shore push increased to 5 grid cells.")


# In[47]:


# Check if the dataset actually has multiple depth levels
z = np.asarray(getattr(fieldset.U, "depth", getattr(fieldset.U.grid, "depth", [])), dtype=float)
print("Depth array:", z)
print("Num depth levels:", len(np.atleast_1d(z)))


# In[48]:


print("===================================")
print("Running simulation at depth:", fieldset.TOP_LEVEL, "meters")
print("===================================")


# In[49]:


print("start_time (sec):", start_time)
print("end_time   (sec):", end_time)
print("runtime days:", (end_time - start_time)/86400)

print("start_dt:", all_times[all_times >= np.datetime64("2009-06-01")][0])
print("all_times max:", all_times.max())


# In[50]:


def DeleteParticle(particle, fieldset, time):
    particle.delete()


# In[51]:


# ============================================================
# FILTER RELEASE POINTS INSIDE SAFE ROMS DOMAIN
# ============================================================

lon_min = float(np.min(fieldset.U.grid.lon))
lon_max = float(np.max(fieldset.U.grid.lon))
lat_min = float(np.min(fieldset.U.grid.lat))
lat_max = float(np.max(fieldset.U.grid.lat))

buffer = 0.25   # ~5 km safety margin

safe_mask = (
    (release_lons > lon_min + buffer) &
    (release_lons < lon_max - buffer) &
    (release_lats > lat_min + buffer) &
    (release_lats < lat_max - buffer)
)

release_lons = release_lons[safe_mask]
release_lats = release_lats[safe_mask]
release_counts = release_counts[safe_mask]

print("Filtered release sites:", len(release_lons))


# In[52]:


# ============================================================
# IMPORT POLYGON TOOLS
# ============================================================

from shapely.geometry import Polygon, Point
import pandas as pd


# In[53]:


def DeleteParticle(particle, fieldset, time):
    particle.delete()


# In[54]:


from datetime import timedelta
from parcels import ErrorCode

# =========================
# CHANGE THESE
# =========================
total_runtime_days = 300      # test run first
chunk_days = 1                # mortality checked once per day
output_name = "OpakapakaOutput.zarr"
# =========================

kernels = build_kernels(pset)

ofile = pset.ParticleFile(
    name=output_name,
    outputdt=timedelta(hours=12),
)

def DeleteParticle(particle, fieldset, time):
    particle.delete()

elapsed_days = 0

while elapsed_days < total_runtime_days:
    step_days = min(chunk_days, total_runtime_days - elapsed_days)

    pset.execute(
        kernels,
        runtime=timedelta(days=step_days),
        dt=timedelta(minutes=2),
        output_file=ofile,
        verbose_progress=False,
        recovery={ErrorCode.ErrorOutOfBounds: DeleteParticle}
    )

    elapsed_days += step_days
    print(f"\nFinished model day {elapsed_days}")

    apply_chunked_mortality(
        pset,
        target_fraction=0.10,
        pld_min_days=0,
        pld_max_days=180
    )

# tiny final flush so last mortality changes get written
pset.execute(
    kernels,
    runtime=timedelta(minutes=2),
    dt=timedelta(minutes=2),
    output_file=ofile,
    verbose_progress=False,
    recovery={ErrorCode.ErrorOutOfBounds: DeleteParticle}
)

print("\nChunked run complete.")


# In[ ]:





# In[ ]:




