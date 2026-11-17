# OPAKAPAKA 2011: FIGURES FOR PROFESSOR
# Updated from plot_output (1).R. Run the whole script in RStudio.
# The Zarr and bathymetry TIFF must be available on this computer.

# 1. SETTINGS ---------------------------------------------------------------
packages <- c('reticulate', 'ggplot2', 'terra')
missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly=TRUE)]
if (length(missing)) stop('Install R packages first: install.packages(c(', paste(sprintf('"%s"',missing), collapse=', '), '))')
library(reticulate)
library(ggplot2)
library(terra)

# IMPORTANT: update these two paths to match your Mac/external drive.
zarr_path <- '/Volumes/One Touch/TPruitt/Output/OpakapakaOutput_2011_30m3.zarr'
bathy_path <- '/Users/justin.suca/Documents/TPruitt/ROMS(0-50)/ETOPO_2022 (Bedrock; 15 arcseconds).tiff'
output_dir <- '/Volumes/One Touch/TPruitt/Output/Opakapaka_2011_Figures'
dir.create(output_dir, recursive=TRUE, showWarnings=FALSE)

# This is the PLD window used in your old poster; change only if your 2011 model differs.
pld_min_days <- 60
pld_max_days <- 180
# Set to TRUE only if the 2011 simulation actually used the 40-60m + 6km mask.
show_habitat_outline <- TRUE
max_map_points <- 150000
batch_size <- 1000L # small batches avoid loading entire 300-day Zarr into R
set.seed(11)

if (!dir.exists(zarr_path)) stop('Zarr folder not found: ', zarr_path, '\nCheck zarr_path.')
if (!file.exists(bathy_path)) stop('Bathymetry TIFF not found: ', bathy_path, '\nCheck bathy_path.')

# 2. LOAD THE ZARR IN BATCHES -----------------------------------------------
# reticulate uses the Python environment where xarray, zarr, numpy are installed.
py_run_string("import numpy as np\nimport xarray as xr")
py$zarr_file <- zarr_path
py_run_string("
ds = xr.open_zarr(zarr_file, consolidated=False)
print('Zarr dimensions:', dict(ds.sizes))
print('Zarr variables:', list(ds.data_vars))
if 'lon' not in ds or 'lat' not in ds:
    raise ValueError('Expected lon and lat arrays. Check the Zarr structure.')
trajectory_dim = next((d for d in ds['lon'].dims if d != 'obs'), None)
if trajectory_dim is None or 'obs' not in ds['lon'].dims:
    raise ValueError('Expected lon and lat with trajectory and obs dimensions.')
n_particles = ds.sizes[trajectory_dim]

def load_batch(start, end):
    sub = ds.isel({trajectory_dim: slice(start, end)})
    def get(name):
        if name not in sub: return None
        arr = np.asarray(sub[name].transpose(trajectory_dim, 'obs').values)
        return arr
    lon, lat = get('lon'), get('lat')
    age, settled, region = get('age'), get('settled'), get('settle_region')
    if settled is None and region is None:
        raise ValueError('Neither settled nor settle_region exists. Cannot calculate connectivity.')
    if region is None:
        raise ValueError('settle_region missing. Cannot identify destination islands.')
    # Detect first saved position for each trajectory, including staggered releases.
    good = np.isfinite(lon) & np.isfinite(lat)
    has_start = good.any(axis=1)
    start_idx = good.argmax(axis=1)
    rows = np.arange(lon.shape[0])
    release_lon = np.where(has_start, lon[rows, start_idx], np.nan)
    release_lat = np.where(has_start, lat[rows, start_idx], np.nan)
    # Prefer explicit settled flag; do NOT treat region alone as proof of settlement
    # if settled is available. Only the first settlement event counts.
    if settled is not None:
        settled_now = (settled == 1)
    else:
        settled_now = np.isfinite(region) & (region >= 0)
    events = settled_now & np.isfinite(region) & (region >= 0) & good
    has_event = events.any(axis=1)
    event_idx = events.argmax(axis=1)
    dest = np.where(has_event, region[rows, event_idx], -1)
    first_lon = np.where(has_event, lon[rows, event_idx], np.nan)
    first_lat = np.where(has_event, lat[rows, event_idx], np.nan)
    first_age = np.where(has_event, age[rows, event_idx] / 86400.0, np.nan) if age is not None else np.full(len(rows), np.nan)
    # Take a modest, reproducible subset of positions from each batch for the dispersal map.
    step = max(1, lon.shape[1] // 35)
    tracks_lon = lon[:, ::step].ravel()
    tracks_lat = lat[:, ::step].ravel()
    valid_tracks = np.isfinite(tracks_lon) & np.isfinite(tracks_lat)
    return dict(release_lon=release_lon, release_lat=release_lat,
                first_lon=first_lon, first_lat=first_lat,
                first_age=first_age, destination=dest, ever_settled=has_event,
                track_lon=tracks_lon[valid_tracks], track_lat=tracks_lat[valid_tracks])
")

n_particles <- as.integer(py$n_particles)
if (n_particles < 1) stop('Zarr contains no trajectories.')
parts <- list(); track_parts <- list()
for (start in seq.int(0L, n_particles - 1L, by=batch_size)) {
  end <- min(n_particles, start + batch_size)
  b <- py$load_batch(as.integer(start), as.integer(end))
  parts[[length(parts)+1L]] <- data.frame(
    release_lon=as.numeric(b$release_lon), release_lat=as.numeric(b$release_lat),
    settle_lon=as.numeric(b$first_lon), settle_lat=as.numeric(b$first_lat),
    age_days=as.numeric(b$first_age), destination=as.integer(b$destination),
    ever_settled=as.logical(b$ever_settled))
  tr <- data.frame(lon=as.numeric(b$track_lon), lat=as.numeric(b$track_lat))
  if (nrow(tr) > 1000) tr <- tr[sample.int(nrow(tr), 1000), ,drop=FALSE]
  track_parts[[length(track_parts)+1L]] <- tr
  if (start %% (batch_size*10L) == 0L) cat('Processed', end, 'of', n_particles, 'trajectories\n')
}
particles <- do.call(rbind, parts)
tracks <- do.call(rbind, track_parts)
rm(parts, track_parts)
if (nrow(tracks) > max_map_points) tracks <- tracks[sample.int(nrow(tracks), max_map_points), ,drop=FALSE]

# 3. REGIONS AND CONNECTIVITY ----------------------------------------------
# These approximate release boxes come from your OLD R script. If your 2011
# model has actual release-region IDs or polygon boundaries, use those instead.
classify_release <- function(lon, lat) {
  out <- rep('Ocean / Unclassified', length(lon))
  set <- function(name, x1, x2, y1, y2) {
    hit <- is.finite(lon) & is.finite(lat) & lon >= x1 & lon <= x2 & lat >= y1 & lat <= y2 & out == 'Ocean / Unclassified'
    out[hit] <<- name
  }
  set('Kaula', -160.75, -160.45, 21.50, 21.72)
  set('Niihau', -160.30, -159.95, 21.75, 22.05)
  set('Kauai', -159.70, -159.20, 21.85, 22.30)
  set('Oahu', -158.40, -157.60, 21.10, 21.80)
  set('Penguin Bank', -157.85, -157.20, 20.80, 21.15)
  set('Maui Nui', -157.40, -156.00, 20.40, 21.20)
  set('Hawaii', -156.00, -154.70, 18.90, 20.30)
  out
}
regions <- c('Kauai','Niihau','Kaula','Oahu','Penguin Bank','Maui Nui','Hawaii')
