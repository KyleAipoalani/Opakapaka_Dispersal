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
# Destination ID mapping follows your previous R script, NOT independently verified
# against the 2011 model. Verify this mapping before reporting numerical results.
region_id <- c('0'='Kauai','1'='Niihau','2'='Kaula','3'='Oahu',
               '4'='Penguin Bank','5'='Maui Nui','6'='Hawaii')
particles$release <- classify_release(particles$release_lon, particles$release_lat)
particles$settlement <- unname(region_id[as.character(particles$destination)])
particles$valid_pld <- particles$ever_settled & is.finite(particles$age_days) &
  particles$age_days >= pld_min_days & particles$age_days <= pld_max_days &
  !is.na(particles$settlement)
valid <- particles[particles$valid_pld & particles$release %in% regions, ,drop=FALSE]
counts <- table(factor(valid$release, levels=regions), factor(valid$settlement, levels=regions))
row_totals <- rowSums(counts)
percent <- sweep(counts, 1, pmax(row_totals,1), '/') * 100
percent[row_totals == 0, ] <- NA_real_ # no settled particles: undefined, NOT 0%
write.csv(as.data.frame.matrix(counts), file.path(output_dir,'2011_connectivity_counts.csv'))
write.csv(round(as.data.frame.matrix(percent),2), file.path(output_dir,'2011_connectivity_percent.csv'))
write.csv(particles[,c('release','settlement','age_days','ever_settled','valid_pld')],
          file.path(output_dir,'2011_particle_settlement_summary.csv'), row.names=FALSE)

# 4. BATHYMETRY -------------------------------------------------------------
bathy <- terra::rast(bathy_path)
if (is.na(terra::crs(bathy)) || terra::crs(bathy) == '') stop('Bathymetry TIFF needs a valid CRS.')
# Crop in the raster's native CRS; transform to longitude/latitude for plots.
area <- terra::vect(matrix(c(-161.5,18.0,-154.5,18.0,-154.5,22.5,-161.5,22.5,-161.5,18.0), ncol=2, byrow=TRUE), type='polygons', crs='EPSG:4326')
area_native <- terra::project(area, terra::crs(bathy))
bathy <- terra::crop(bathy, area_native)
if (!terra::same.crs(bathy, 'EPSG:4326')) bathy <- terra::project(bathy, 'EPSG:4326', method='bilinear')
# Reduce plotting cost only; keep full raster for the habitat outline.
bathy_plot <- if (terra::ncell(bathy) > 350000) terra::aggregate(bathy, fact=ceiling(sqrt(terra::ncell(bathy)/350000)), fun=mean, na.rm=TRUE) else bathy
bathy_df <- as.data.frame(bathy_plot, xy=TRUE, na.rm=TRUE)
names(bathy_df) <- c('lon','lat','depth')

habitat <- NULL
if (show_habitat_outline) {
  # Approximate 40-60m bathymetric band plus 6 km, shown for visual reference.
  # This does NOT reconstruct the exact model settlement mask.
  band <- terra::ifel(bathy >= -60 & bathy <= -40, 1, 0)
  # Project to meters before buffering; EPSG:32604 covers the main islands.
  band_utm <- terra::project(band, 'EPSG:32604', method='near')
  # Fast focal dilation using a circular 6 km kernel (can take time).
  w <- terra::focalMat(band_utm, 6000, type='circle')
  w[!is.na(w)] <- 1
  buffered <- terra::focal(band_utm, w=w, fun='max', na.policy='omit', fillvalue=0)
  habitat <- terra::project(buffered, bathy_plot, method='near')
  habitat <- terra::ifel(bathy_plot < 0, habitat, 0)
}

# 5. SAVE FIGURES -----------------------------------------------------------
base_map <- function(title, subtitle=NULL) {
  ggplot(bathy_df, aes(lon,lat)) +
    geom_raster(aes(fill=depth)) +
    scale_fill_gradientn(colours=c('#24134f','#276b99','#40a9a0','#d6d9a7'), name='Elevation (m)') +
    coord_quickmap(xlim=c(-161.5,-154.5), ylim=c(18,22.5), expand=FALSE) +
    labs(title=title, subtitle=subtitle, x='Longitude', y='Latitude') +
    theme_minimal(base_size=13) + theme(panel.grid=element_blank(), plot.title=element_text(face='bold'))
}
add_habitat <- function(p) {
  if (is.null(habitat)) return(p)
  # Use a contour line to show approximate habitat extent.
  h <- as.data.frame(habitat, xy=TRUE, na.rm=TRUE)
  names(h) <- c('lon','lat','mask')
  p + geom_contour(data=h, aes(x=lon,y=lat,z=mask), breaks=0.5,
                   colour='#ff9b21', linewidth=0.5, inherit.aes=FALSE)
}
save_map <- function(p, name) ggsave(file.path(output_dir,name), p, width=12, height=7.5, dpi=220, bg='white')

p1 <- base_map('2011 Hawaiian Islands Bathymetry', 'ETOPO 2022 bedrock elevation')
save_map(p1, '2011_01_bathymetry_heatmap.png')

release_points <- particles[is.finite(particles$release_lon) & is.finite(particles$release_lat),]
if (nrow(release_points)>max_map_points) release_points <- release_points[sample.int(nrow(release_points),max_map_points),]
p2 <- add_habitat(base_map('2011 Opakapaka Release Locations', 'Orange outline: approximate 40-60 m habitat with 6 km buffer')) +
  geom_point(data=release_points, aes(x=release_lon,y=release_lat), inherit.aes=FALSE,
             colour='#fff07a', alpha=0.4, size=0.22)
save_map(p2, '2011_02_release_locations.png')

settle_points <- particles[particles$valid_pld & is.finite(particles$settle_lon) & is.finite(particles$settle_lat),]
if (nrow(settle_points)>max_map_points) settle_points <- settle_points[sample.int(nrow(settle_points),max_map_points),]
p3 <- add_habitat(base_map('2011 First-Entry Settlement Locations',
                            sprintf('First valid settlement, %d-%d days PLD',pld_min_days,pld_max_days))) +
  geom_point(data=settle_points, aes(x=settle_lon,y=settle_lat), inherit.aes=FALSE,
             colour='yellow', alpha=0.65, size=0.32)
save_map(p3, '2011_03_first_entry_settlement.png')

heat_df <- as.data.frame(as.table(percent), stringsAsFactors=FALSE)
names(heat_df) <- c('Release','Settlement','Percent')
heat_df$Release <- factor(heat_df$Release, levels=rev(regions))
heat_df$Settlement <- factor(heat_df$Settlement, levels=regions)
p4 <- ggplot(heat_df, aes(Settlement,Release,fill=Percent)) +
  geom_tile(colour='white',linewidth=0.8) +
  geom_text(aes(label=ifelse(is.na(Percent),'N/A',sprintf('%.1f%%',Percent))),size=3) +
  scale_fill_gradient(low='white',high='#6A0DAD',na.value='#e4e4e4',name='Settlement (%)') +
  labs(title='2011 Opakapaka Larval Connectivity',
       subtitle='Rows = release region; columns = first-entry settlement region; percentages among valid settled particles per release region',
       x='Settlement Region',y='Release Region') +
  theme_minimal(base_size=12) + theme(panel.grid=element_blank(),axis.text.x=element_text(angle=35,hjust=1),plot.title=element_text(face='bold'))
ggsave(file.path(output_dir,'2011_04_connectivity_heatmap.png'),p4,width=12,height=7.5,dpi=250,bg='white')

if (nrow(tracks)>0) {
  p5 <- base_map('2011 Opakapaka Particle Dispersal', 'Sampled particle positions from saved model observations') +
    geom_point(data=tracks,aes(x=lon,y=lat),inherit.aes=FALSE,colour='#ef8b2c',size=0.1,alpha=0.25)
  save_map(p5,'2011_05_particle_dispersal.png')
}

# Counts and percentage matrix also saved as CSV for professor/poster.
summary_df <- data.frame(
  Metric=c('Total trajectories','Ever settled (first-entry recorded)',
           sprintf('First settlement at %d-%d days',pld_min_days,pld_max_days),
           'Valid settlements with classified release region',
           'No recorded first-entry settlement'),
  Count=c(nrow(particles),sum(particles$ever_settled,na.rm=TRUE),
          sum(particles$valid_pld,na.rm=TRUE),nrow(valid),
          sum(!particles$ever_settled,na.rm=TRUE)))
# No settlement does NOT mean death. Mortality requires separate kill_reason analysis.
write.csv(summary_df,file.path(output_dir,'2011_simulation_summary.csv'),row.names=FALSE)
print(summary_df)
cat('\nAll plots and tables saved to:', output_dir, '\n')
cat('IMPORTANT: Confirm destination region IDs and release polygons against the 2011 Python run.\n')
cat('IMPORTANT: This script cannot verify that the new Zarr data or TIFF exists until you run it on the Mac.\n')
