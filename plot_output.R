# note: run separately in sections by comments
#      ie.  highlight library first then run, 
#           highlight virtual env next then run, dataset next, etc.


library(ggplot2)
library(dplyr)
# install.packages("reticulate")
library(reticulate) 
library(terra)
library(sf)
library(viridis)


xr <- import("xarray")
np <- import("numpy")
pd <- import("pandas")
library(reticulate)

zarr <- import("zarr")

zarr_path <- "/Users/justin.suca/Documents/TPruitt/Python Script/OpakapakaOutput_2010_30m.zarr"

zg <- zarr$open_group(zarr_path, mode = "r")

lon_arr <- py_to_r(zg[["lon"]]$`__getitem__`(tuple()))
lat_arr <- py_to_r(zg[["lat"]]$`__getitem__`(tuple()))
age_arr <- py_to_r(zg[["age"]]$`__getitem__`(tuple()))
settled_arr <- py_to_r(zg[["settled"]]$`__getitem__`(tuple()))
settle_region_arr <- py_to_r(zg[["settle_region"]]$`__getitem__`(tuple()))

min_obs <- min(
  ncol(lon_arr),
  ncol(lat_arr),
  ncol(age_arr),
  ncol(settled_arr),
  ncol(settle_region_arr)
)

lon_arr <- lon_arr[, 1:min_obs]
lat_arr <- lat_arr[, 1:min_obs]
age_arr <- age_arr[, 1:min_obs]
settled_arr <- settled_arr[, 1:min_obs]
settle_region_arr <- settle_region_arr[, 1:min_obs]

cat("Loaded fixed arrays with", min_obs, "time steps\n")

# ============================================================
# RELEASE LOCATIONS
# ============================================================

release_lon <- lon_arr[, 1]
release_lat <- lat_arr[, 1]

# ============================================================
# FIRST-ENTRY SETTLEMENT POINTS
# ============================================================

first_settle_obs <- apply(settled_arr, 1, function(x) {
  idx <- which(x == 1)
  if (length(idx) == 0) NA_integer_ else idx[1]
})

settled_particles <- which(!is.na(first_settle_obs))
first_idx <- first_settle_obs[settled_particles]

first_lon <- lon_arr[cbind(settled_particles, first_idx)]
first_lat <- lat_arr[cbind(settled_particles, first_idx)]
first_age_days <- age_arr[cbind(settled_particles, first_idx)] / 86400
first_region <- settle_region_arr[cbind(settled_particles, first_idx)]

valid_pld <- first_age_days >= 60 &
  first_age_days <= 180 &
  is.finite(first_lon) &
  is.finite(first_lat) &
  first_region >= 0

first_lon_plot <- first_lon[valid_pld]
first_lat_plot <- first_lat[valid_pld]

cat("Total particles:", nrow(lon_arr), "\n")
cat("Ever settled:", length(settled_particles), "\n")
cat("Settled inside 60–180 days:", length(first_lon_plot), "\n")

# ============================================================
# BATHYMETRY + HABITAT MASK
# ============================================================

bathy_path <- "/Users/justin.suca/Documents/TPruitt/ROMS(0-50)/ETOPO_2022 (Bedrock; 15 arcseconds).tiff"

bathy <- rast(bathy_path)

hawaii_extent <- ext(-161.5, -154.5, 18.0, 22.5)
bathy_crop <- crop(bathy, hawaii_extent)

habitat_band <- bathy_crop <= -40 & bathy_crop >= -60

km_per_cell <- 0.463
buffer_km <- 6
buffer_cells <- ceiling(buffer_km / km_per_cell)

r <- buffer_cells
x <- -r:r
y <- -r:r

w <- outer(x, y, function(a, b) sqrt(a^2 + b^2) <= r)
w <- matrix(as.numeric(w), nrow = length(x), ncol = length(y))

habitat_mask <- focal(
  habitat_band,
  w = w,
  fun = max,
  na.policy = "omit",
  fillvalue = 0
)

habitat_mask <- habitat_mask & (bathy_crop < 0)

# ============================================================
# PLOT 1: FIRST-ENTRY SETTLEMENT MAP
# ============================================================

png("~/Documents/2010_first_entry_settlement_locations.png", width = 2200, height = 1500, res = 250)

plot(
  bathy_crop,
  col = hcl.colors(120, "viridis"),
  main = "2010 First-Entry Settlement Locations",
  xlab = "Longitude",
  ylab = "Latitude",
  axes = TRUE,
  legend = FALSE
)

contour(
  habitat_mask,
  levels = 0.5,
  add = TRUE,
  drawlabels = FALSE,
  col = "orange",
  lwd = 2
)

points(
  first_lon_plot,
  first_lat_plot,
  pch = 16,
  cex = 0.19,
  col = rgb(1, 1, 0, 0.75)
)

text(-158.95, 22.12, "Kauai", col = "white", cex = 0.9)
text(-159.90, 21.65, "Niihau", col = "white", cex = 0.9)
text(-160.60, 21.20, "Kaula", col = "white", cex = 0.9)
text(-158.05, 21.95, "Oahu", col = "white", cex = 0.9)
text(-157.85, 20.60, "Penguin Bank", col = "white", cex = 0.8)
text(-156.95, 21.45, "Maui Nui", col = "white", cex = 0.9)
text(-155.70, 19.80, "Hawaii", col = "white", cex = 0.9)

dev.off()

# ============================================================
# CONNECTIVITY MATRIX
# ============================================================

classify_release <- function(lon, lat) {
  if (lon >= -160.75 && lon <= -160.45 && lat >= 21.50 && lat <= 21.72) return("Kaula")
  if (lon >= -160.30 && lon <= -159.95 && lat >= 21.75 && lat <= 22.05) return("Niihau")
  if (lon >= -159.70 && lon <= -159.20 && lat >= 21.85 && lat <= 22.30) return("Kauai")
  if (lon >= -158.40 && lon <= -157.60 && lat >= 21.10 && lat <= 21.80) return("Oahu")
  if (lon >= -157.85 && lon <= -157.20 && lat >= 20.80 && lat <= 21.15) return("Penguin_Bank")
  if (lon >= -157.40 && lon <= -156.00 && lat >= 20.40 && lat <= 21.20) return("Maui_Nui")
  if (lon >= -156.00 && lon <= -154.70 && lat >= 18.90 && lat <= 20.30) return("Hawaii")
  return("Ocean")
}

release_name <- mapply(classify_release, release_lon, release_lat)

region_names <- c(
  "0" = "Kauai",
  "1" = "Niihau",
  "2" = "Kaula",
  "3" = "Oahu",
  "4" = "Penguin_Bank",
  "5" = "Maui_Nui",
  "6" = "Hawaii"
)

settle_region_first <- first_region[valid_pld]
settle_name <- region_names[as.character(settle_region_first)]

release_valid <- release_name[settled_particles][valid_pld]

connectivity_counts <- table(release_valid, settle_name)

connectivity_percent <- round(
  prop.table(connectivity_counts, 1) * 100,
  2
)

print(connectivity_percent)

# ============================================================
# PLOT 2: HEATMAP
# ============================================================

heat_df <- as.data.frame(as.table(connectivity_percent))
names(heat_df) <- c("Release", "Settlement", "Percent")

heat_df <- heat_df %>%
  filter(Release != "Ocean", Settlement != "Ocean")

heat_df$Release <- gsub("_", " ", heat_df$Release)
heat_df$Settlement <- gsub("_", " ", heat_df$Settlement)

region_order <- c("HawaiiKauai", "Maui Nui", "Penguin Bank", "Oahu", "Kauai", "Niihau", "Kaula")

heat_df$Release <- factor(heat_df$Release, levels = region_order)
heat_df$Settlement <- factor(heat_df$Settlement, levels = region_order)

heatmap_plot <- ggplot(heat_df, aes(x = Settlement, y = Release, fill = Percent)) +
  geom_tile(color = "white", linewidth = 0.8) +
  scale_fill_gradient(
    name = "Settlement (%)",
    low = "white",
    high = "#6A0DAD"
  ) +
  labs(
    title = "Larval Connectivity Between Hawaiian Islands",
    subtitle = "2010 Opakapaka (June–September, 30 m depth, 60–180 day PLD)",
    x = "Settlement Region",
    y = "Release Region"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    plot.title = element_text(size = 26, face = "bold"),
    plot.subtitle = element_text(size = 16),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 13),
    panel.grid = element_blank()
  )

ggsave(
  "2010_connectivity_heatmap_clean.png",
  heatmap_plot,
  width = 10,
  height = 7,
  dpi = 300,
  bg = "white"
)

# ============================================================
# PLOT 3: PARTICLE TRAJECTORY / DISPERSAL SNAPSHOT
# ============================================================


particle_lon <- as.vector(lon_arr)
particle_lat <- as.vector(lat_arr)

valid_particles <- is.finite(particle_lon) & is.finite(particle_lat)

particle_df <- data.frame(
  lon = particle_lon[valid_particles],
  lat = particle_lat[valid_particles]
)
# ============================================================
# PLOT 3: DENSE PARTICLE DISPERSAL MAP - OCEAN ONLY
# ============================================================

# flatten all particle positions across all saved timesteps
particle_lon <- as.vector(lon_arr)
particle_lat <- as.vector(lat_arr)

# remove bad values
valid_particles <- is.finite(particle_lon) & is.finite(particle_lat)

particle_df <- data.frame(
  lon = particle_lon[valid_particles],
  lat = particle_lat[valid_particles]
)

# remove particles on land using bathymetry
bathy_vals <- terra::extract(
  bathy_crop,
  particle_df[, c("lon", "lat")]
)

ocean_keep <- is.finite(bathy_vals[, 2]) & bathy_vals[, 2] < 0

particle_df_ocean <- particle_df[ocean_keep, ]

cat("Ocean particle points plotted:", nrow(particle_df_ocean), "\n")

# optional: thin only if too many points
set.seed(1)
max_points <- 1000000

if (nrow(particle_df_ocean) > max_points) {
  particle_df_ocean <- particle_df_ocean[
    sample(seq_len(nrow(particle_df_ocean)), max_points),
  ]
}

png("~/Documents/2010_particle_dispersal_dense_ocean_only.png",
    width = 2200, height = 1500, res = 250)

plot(
  particle_df_ocean$lon,
  particle_df_ocean$lat,
  pch = 16,
  cex = 0.05,
  col = rgb(0.75, 0.35, 0.05, 0.25),
  xlim = c(-161.3, -154.6),
  ylim = c(18.4, 22.5),
  xlab = "Longitude",
  ylab = "Latitude",
  main = "2010 Opakapaka Particle Dispersal"
)

points(
  release_lon,
  release_lat,
  pch = 16,
  cex = 0.12,
  col = rgb(0.1, 0.35, 0.9, 0.55)
)

legend(
  "topright",
  legend = c("Release", "Particles"),
  col = c(rgb(0.1, 0.35, 0.9, 0.6), rgb(0.75, 0.35, 0.05, 0.6)),
  pch = 16,
  pt.cex = 1.2,
  bty = "n"
)

dev.off()

























##############################


target_obs <- min_obs

lon_target <- lon_arr[, target_obs]
lat_target <- lat_arr[, target_obs]

valid_target <- is.finite(lon_target) & is.finite(lat_target)

set.seed(1)
keep_release <- sample(seq_along(release_lon), min(500000, length(release_lon)))
keep_target <- sample(which(valid_target), min(500000, sum(valid_target)))

png("2010_particle_dispersal_snapshot1.png", width = 2200, height = 1500, res = 250)

plot(
  release_lon[keep_release],
  release_lat[keep_release],
  pch = 16,
  cex = 0.18,
  col = rgb(0.1, 0.35, 0.9, 0.35),
  xlim = c(-161.3, -154.6),
  ylim = c(18.4, 22.5),
  xlab = "Longitude",
  ylab = "Latitude",
  main = "2010 Opakapaka Particle Dispersal Snapshot"
)

points(
  lon_target[keep_target],
  lat_target[keep_target],
  pch = 16,
  cex = 0.25,
  col = rgb(1, 0.45, 0, 0.65)
)

legend(
  "topright",
  legend = c("Release", "Particles"),
  col = c(rgb(0.1, 0.35, 0.9, 0.6), rgb(1, 0.45, 0, 0.7)),
  pch = 16,
  pt.cex = 1.2,
  bty = "n"
)

dev.off()








library(gridExtra)
library(grid)

# ============================================================
# CLEAN POSTER CONNECTIVITY TABLE
# ============================================================

poster_table <- as.data.frame.matrix(connectivity_percent_60_180)

# Add release region as first column
poster_table <- cbind(
  "Release Region" = rownames(poster_table),
  poster_table
)

# Clean underscores
names(poster_table) <- gsub("_", " ", names(poster_table))
poster_table$`Release Region` <- gsub("_", " ", poster_table$`Release Region`)

# Add percent signs
for (col in 2:ncol(poster_table)) {
  poster_table[[col]] <- paste0(poster_table[[col]], "%")
}

# ============================================================
# MAKE TABLE WITH LIGHT POSTER STYLING
# ============================================================

table_plot <- tableGrob(
  poster_table,
  rows = NULL,
  theme = ttheme_default(
    core = list(
      fg_params = list(
        fontsize = 13
      ),
      bg_params = list(
        fill = rep(c("white", "gray95"), length.out = nrow(poster_table)),
        col = "gray70"
      ),
      padding = unit(c(5, 5), "mm")
    ),
    colhead = list(
      fg_params = list(
        fontsize = 14,
        fontface = "bold"
      ),
      bg_params = list(
        fill = "gray85",
        col = "gray60"
      ),
      padding = unit(c(5, 5), "mm")
    )
  )
)

# ============================================================
# SAVE CLEAN POSTER TABLE
# ============================================================

png(
  "2010_connectivity_matrix_percent_poster_clean.png",
  width = 2400,
  height = 950,
  res = 250
)

grid.newpage()
grid.draw(table_plot)

dev.off()




# Trevellʻs code is under here


# ============================================================
# LOAD ZARR
# ============================================================

xr <- import("xarray", convert = TRUE)
ds <- xr$open_zarr(zarr_path)

# ============================================================
# GET FIRST-ENTRY SETTLEMENT POINTS
# ============================================================

settled_arr <- py_to_r(ds$settled$values)
lon_arr <- py_to_r(ds$lon$values)
lat_arr <- py_to_r(ds$lat$values)
age_arr <- py_to_r(ds$age$values)
settle_region_arr <- py_to_r(ds$settle_region$values)

first_settle_obs <- apply(settled_arr, 1, function(x) {
  idx <- which(x == 1)
  if (length(idx) == 0) NA_integer_ else idx[1]
})

settled_particles <- which(!is.na(first_settle_obs))

first_lon <- mapply(function(i, j) lon_arr[i, j], settled_particles, first_settle_obs[settled_particles])
first_lat <- mapply(function(i, j) lat_arr[i, j], settled_particles, first_settle_obs[settled_particles])
first_age_days <- mapply(function(i, j) age_arr[i, j], settled_particles, first_settle_obs[settled_particles]) / 86400
first_region <- mapply(function(i, j) settle_region_arr[i, j], settled_particles, first_settle_obs[settled_particles])

valid_pld <- first_age_days >= 60 &
  first_age_days <= 180 &
  is.finite(first_lon) &
  is.finite(first_lat) &
  first_region != -1

first_lon_plot <- first_lon[valid_pld]
first_lat_plot <- first_lat[valid_pld]

cat("\nTotal particles:", length(first_settle_obs), "\n")
cat("Ever settled:", length(settled_particles), "\n")
cat("Settled inside 60–180 days:", length(first_lon_plot), "\n")

# ============================================================
# ============================================================
# LOAD BATHYMETRY
# ============================================================

bathy_path <- "/Users/justin.suca/Documents/TPruitt/ROMS(0-50)/ETOPO_2022 (Bedrock; 15 arcseconds).tiff"

bathy <- rast(bathy_path)

hawaii_extent <- ext(-161.5, -154.5, 18.0, 22.5)
bathy_crop <- crop(bathy, hawaii_extent)

# ============================================================
# MAKE SETTLEMENT HABITAT MASK: 40-60 m + 6 km CIRCULAR BUFFER
# ============================================================

habitat_band <- bathy_crop <= -40 & bathy_crop >= -60

km_per_cell <- 0.463
buffer_km <- 6
buffer_cells <- ceiling(buffer_km / km_per_cell)

cat("Buffer cells:", buffer_cells, "\n")
cat("Actual buffer km:", buffer_cells * km_per_cell, "\n")

# circular buffer window instead of square window
r <- buffer_cells
x <- -r:r
y <- -r:r

w <- outer(
  x,
  y,
  function(a, b) sqrt(a^2 + b^2) <= r
)

w <- matrix(
  as.numeric(w),
  nrow = length(x),
  ncol = length(y)
)

habitat_mask <- focal(
  habitat_band,
  w = w,
  fun = max,
  na.policy = "omit",
  fillvalue = 0
)

# remove land from buffered habitat mask
habitat_mask <- habitat_mask & (bathy_crop < 0)

# testing habitat mask only
# plot (
#   habitat_mask,
#   main = "Habitat mask only"
# )

# ============================================================
# PLOT NICE MAP
# ============================================================

plot(
  bathy_crop,
  col = hcl.colors(120, "viridis"),
  main = "Opakapaka Settlement Locations (2010 / 30m)",
  xlab = "Longitude",
  ylab = "Latitude",
  axes = TRUE,
  legend = FALSE
)

# remove the white coastline contour because it creates weird inner squiggles
# contour(bathy_crop, levels = 0, add = TRUE)

# buffered settlement habitat outline
contour(
  habitat_mask,
  levels = 0.5,
  add = TRUE,
  drawlabels = FALSE,
  col = "orange",
  lwd = 2
)


keep <- seq_along(first_lon_plot)


cat("Points being plotted:", length(keep), "\n")

# settlement particles
points(
  first_lon_plot[keep],
  first_lat_plot[keep],
  pch = 16,
  cex = 0.19,
  col = rgb(1, 1, 0, 0.75)
)

# ============================================
# TEST: drawing connectivity regions
# output looks rectangular so maybe
# ============================================
region_boxes <- data.frame(
  region = c("Kauai", "Niihau", "Kaula", "Oahu", "Penguin Bank", "Maui Nui", "Hawaii"),
  xmin = c(-159.70, -160.30, -160.75, -158.40, -157.85, -157.40, -156.00),
  xmax = c(-159.20, -159.95, -160.45, -157.60, -157.20, -156.00, -154.70),
  ymin = c(21.85, 21.75, 21.50, 21.10, 20.80, 20.40, 18.90),
  ymax = c(22.30, 22.05, 21.72, 21.80, 21.15, 21.20, 20.30)
)
# coords taken directly from run code

for (i in seq_len(nrow(region_boxes))) {
  rect(
    xleft = region_boxes$xmin[i],
    ybottom = region_boxes$ymin[i],
    xright = region_boxes$xmax[i],
    ytop = region_boxes$ymax[i],
    border = "red",
    lwd = 2
  )
}


# labels
# labels (adjusted positions)

# Kauai → move right
text(-158.95, 22.12, "Kauai", col = "white", cex = 0.9)

# Niihau → move down + right
text(-159.90, 21.65, "Niihau", col = "white", cex = 0.9)

# Kaula → move down
text(-160.60, 21.20, "Kaula", col = "white", cex = 0.9)

# Oahu → move slightly up
text(-158.05, 21.95, "Oahu", col = "white", cex = 0.9)

# Penguin Bank → move slightly down
text(-157.85, 20.60, "Penguin Bank", col = "white", cex = 0.8)

# Maui Nui → move slightly up
text(-156.95, 21.45, "Maui Nui", col = "white", cex = 0.9)

# Hawaii (this one was fine)
text(-155.70, 19.80, "Hawaii", col = "white", cex = 0.9)




















#STATS===============

library(gridExtra)
library(grid)

# ============================================================
# SUMMARY NUMBERS
# ============================================================

total_particles <- length(first_settle_obs)
ever_settled_n <- length(settled_particles)
settled_60_180_n <- length(first_lon_plot)

# If your Zarr has kill_reason, use it for died count
if (exists("kill_reason_arr")) {
  died_n <- sum(apply(kill_reason_arr, 1, max, na.rm = TRUE) == 2, na.rm = TRUE)
} else {
  died_n <- total_particles - ever_settled_n
}

summary_df <- data.frame(
  Metric = c(
    "Total Particles",
    "Ever Settled",
    "Settled Within 60–180 Days",
    "Did Not Settle / Died"
  ),
  Count = c(
    total_particles,
    ever_settled_n,
    settled_60_180_n,
    died_n
  )
)

summary_df$Count <- format(summary_df$Count, big.mark = ",")

# ============================================================
# CONNECTIVITY COUNTS TABLE
# ============================================================

counts_df <- as.data.frame.matrix(connectivity_counts_60_180)

counts_df <- cbind(
  "Release Region" = rownames(counts_df),
  counts_df
)

names(counts_df) <- gsub("_", " ", names(counts_df))
counts_df$`Release Region` <- gsub("_", " ", counts_df$`Release Region`)

# Add commas to counts
for (col in 2:ncol(counts_df)) {
  counts_df[[col]] <- format(counts_df[[col]], big.mark = ",")
}

# ============================================================
# MAKE SUMMARY TABLE
# ============================================================

summary_table <- tableGrob(
  summary_df,
  rows = NULL,
  theme = ttheme_default(
    core = list(
      fg_params = list(fontsize = 13),
      bg_params = list(fill = "white", col = "gray70"),
      padding = unit(c(5, 5), "mm")
    ),
    colhead = list(
      fg_params = list(fontsize = 14, fontface = "bold", col = "white"),
      bg_params = list(fill = "#6A0DAD", col = "#6A0DAD"),
      padding = unit(c(5, 5), "mm")
    )
  )
)

# ============================================================
# MAKE CONNECTIVITY TABLE
# ============================================================

counts_table <- tableGrob(
  counts_df,
  rows = NULL,
  theme = ttheme_default(
    core = list(
      fg_params = list(fontsize = 11),
      bg_params = list(
        fill = rep(c("white", "gray95"), length.out = nrow(counts_df)),
        col = "gray70"
      ),
      padding = unit(c(4, 4), "mm")
    ),
    colhead = list(
      fg_params = list(fontsize = 12, fontface = "bold", col = "white"),
      bg_params = list(fill = "#6A0DAD", col = "#6A0DAD"),
      padding = unit(c(4, 4), "mm")
    )
  )
)

# ============================================================
# SAVE SEPARATE PNG FILES
# ============================================================

png("2010_simulation_summary.png", width = 1400, height = 500, res = 250)
grid.newpage()
grid.draw(summary_table)
dev.off()

png("2010_connectivity_matrix_counts_60_180.png", width = 2600, height = 1100, res = 250)
grid.newpage()
grid.draw(counts_table)
dev.off()

# ============================================================
# SAVE BOTH TOGETHER AS ONE PNG
# ============================================================

combined_plot <- arrangeGrob(
  summary_table,
  counts_table,
  ncol = 1,
  heights = c(1, 2.4)
)

png("2010_combined_connectivity_summary_counts.png", width = 2800, height = 1700, res = 250)
grid.newpage()
grid.draw(combined_plot)
dev.off()


# HEAT MAP=========================================

library(ggplot2)

# ============================================================

heat_df <- as.data.frame(as.table(connectivity_percent_60_180))
names(heat_df) <- c("Release", "Settlement", "Percent")

# remove Ocean
heat_df <- heat_df %>%
  dplyr::filter(Release != "Ocean", Settlement != "Ocean")

# clean labels
heat_df$Release <- gsub("_", " ", heat_df$Release)
heat_df$Settlement <- gsub("_", " ", heat_df$Settlement)

# ============================================================
# ORDER REGIONS (geographic)
# ============================================================

region_order <- c("Kauai", "Niihau", "Kaula", "Oahu", "Penguin Bank", "Maui Nui", "Hawaii")

heat_df$Release <- factor(heat_df$Release, levels = region_order)
heat_df$Settlement <- factor(heat_df$Settlement, levels = region_order)

# ============================================================
# HEATMAP
# ============================================================

heatmap_plot <- ggplot(heat_df, aes(x = Settlement, y = Release, fill = Percent)) +
  geom_tile(color = "white", linewidth = 0.8) +
  
  scale_fill_gradient(
    name = "Settlement (%)",
    low = "white",
    high = "#6A0DAD"
  ) +
  
  labs(
    title = "Larval Connectivity Between Hawaiian Islands",
    subtitle = "2010 Opakapaka (June–September, 30 m depth, 60–180 day PLD)",
    x = "Settlement Region",
    y = "Release Region"
  ) +
  
  theme_minimal(base_size = 16) +
  
  theme(
    plot.title = element_text(size = 26, face = "bold"),
    plot.subtitle = element_text(size = 16),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 14),
    
    axis.text.x = element_text(angle = 45, hjust = 1),
    
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 13),
    
    panel.grid = element_blank()
  )

print(heatmap_plot)

ggsave(
  "2010_connectivity_heatmap_clean.png",
  heatmap_plot,
  width = 10,
  height = 7,
  dpi = 300,
  bg = "white"
)

#SPECIES TABLE =======================
#library(gridExtra)
#library(grid)

#species_table <- data.frame(
#  `Spawning Period` = "Jun–Sep",
#  `OVM / Depth` = "Surface to mid-water; 30 m (ROMS)",
#  PLD = "60–180 days",
#  `Settlement Habitat` = "Soft substrate, 60–100 m",
#  `Release Strategy` = "Daily, habitat-weighted",
#  `Settlement Rule` = "First-entry",
#  check.names = FALSE
#)


#table_plot <- tableGrob(
#  species_table,
#  rows = NULL,
#  theme = ttheme_default(
#    core = list(
#      fg_params = list(fontsize = 14),
#      bg_params = list(fill = "white", col = "gray70"),
#      padding = unit(c(5, 5), "mm")
#    ),
#    colhead = list(
#      fg_params = list(fontsize = 15, fontface = "bold", col = "white"),
#      bg_params = list(fill = "#6A0DAD", col = "#6A0DAD"),
#      padding = unit(c(5, 5), "mm")
#    )
#  )
#)
#png("species_table_compact2.png", width = 3500, height = 350, res = 250)
#grid.newpage()
#grid.draw(table_plot)
#dev.off()











#plot OG================

# ============================================================
# DISPERSAL PLOT USING REAL MODEL TIMESTAMP
# Release locations vs October 1, 2009 positions
# ============================================================

# library(reticulate)

# Make sure ds is already loaded
# xr <- import("xarray", convert = TRUE)
# ds <- xr$open_zarr(zarr_path)

# lon_arr <- py_to_r(ds$lon$values)
# lat_arr <- py_to_r(ds$lat$values)
# time_arr <- py_to_r(ds$time$values)

# release locations = first saved position
# release_lon <- lon_arr[, 1]
# release_lat <- lat_arr[, 1]

# target timestamp based on your release window
# target_time <- as.POSIXct("2009-10-01 00:00:00", tz = "UTC")

# convert model time to POSIXct if needed
# time_posix <- time_arr
# find each particle's location closest to Oct 1, 2009
# idx_target <- apply(time_posix, 1, function(x) {
#   good <- which(!is.na(x))
#   if (length(good) == 0) return(NA_integer_)
#   good[which.min(abs(as.numeric(difftime(x[good], target_time, units = "secs"))))]
# })

# valid_target <- !is.na(idx_target)

# lon_target <- mapply(function(i, j) lon_arr[i, j], which(valid_target), idx_target[valid_target])
# lat_target <- mapply(function(i, j) lat_arr[i, j], which(valid_target), idx_target[valid_target])

# sample points so plot is not too heavy
# set.seed(1)
# keep_release <- sample(seq_along(release_lon), min(50000, length(release_lon)))
# keep_target <- sample(seq_along(lon_target), min(50000, length(lon_target)))
# ============================================================
# REMOVE PARTICLES THAT ARE ON LAND USING BATHYMETRY
# ============================================================

# target_points <- data.frame(
#   lon = lon_target,
#   lat = lat_target
# )

# Extract bathymetry value at each particle location
# bathy_values <- terra::extract(
#   bathy_crop,
#   target_points[, c("lon", "lat")]
# )

# Keep only ocean points where depth is below 0
# ocean_keep <- bathy_values[, 2] < 0 & is.finite(bathy_values[, 2])

# lon_target_ocean <- lon_target[ocean_keep]
# lat_target_ocean <- lat_target[ocean_keep]

# cat("Original target particles:", length(lon_target), "\n")
# cat("Ocean-only target particles:", length(lon_target_ocean), "\n")
# ============================================================
# PLOT
# ============================================================
# ============================================================
# SET ZOOM (DO THIS FIRST)
# ============================================================
# xlim <- c(-161.3, -154.6)
# ylim <- c(18.4, 22.5)

# ============================================================
# BASE PLOT (RELEASE POINTS)
# ============================================================
# make space on the right side
# par(mar = c(5, 5, 4, 8))  # extra space on right

# your plot
# plot(
#   release_lon[keep_release],
#   release_lat[keep_release],
#   pch = 16,
#   cex = 0.18,
#   col = rgb(0.1, 0.35, 0.9, 0.35),
#   xlim = xlim,
#   ylim = ylim,
#   xlab = "Longitude",
#   ylab = "Latitude",
#   main = "Opakapaka Dispersal — June 1 to October 1, 2009"
# )

# particles
# points(
#   lon_target_ocean[keep_target],
#   lat_target_ocean[keep_target],
#   pch = 16,
#   cex = 0.28,
#   col = rgb(1, 0.45, 0, 0.65)
# )

# 👇 legend OUTSIDE plot (this fixes everything)
# legend(
#   "topright",
#   inset = c(-0.25, 0),   # pushes it outside right
#   legend = c("Release", "Particles"),
#   col = c(rgb(0.1, 0.35, 0.9, 0.6), rgb(1, 0.45, 0, 0.7)),
#   pch = 16,
#   pt.cex = 1.2,
#   bty = "n",
#   xpd = TRUE             # allows drawing outside plot
# )

#===============================

# library(DiagrammeR)
# library(DiagrammeRsvg)
# library(rsvg)

# pipeline <- grViz("
# digraph pipeline {

#   graph [layout = dot, rankdir = LR, bgcolor = white]

#   node [
#     shape = box,
#     style = 'rounded,filled',
#     fillcolor = '#6A0DAD',
#     color = '#4B0082',
#     fontcolor = white,
#     fontsize = 18,
#     fontname = Helvetica,
#     margin = 0.18
#   ]

#   edge [
#     color = '#555555',
#     penwidth = 2,
#     arrowsize = 0.8
#   ]

#   A [label = 'Release Sites\\nHabitat-based seeding']
#   B [label = 'ROMS Currents\\n2009 velocity fields']
#   C [label = 'Particle Drift\\nOceanParcels at 30 m']
#   D [label = 'Settlement Habitat\\n40–60 m + 6 km buffer']
#   E [label = 'Connectivity Matrix\\nFirst-entry, 60–180 days']

#   A -> B -> C -> D -> E
# }
# ")

# # Save as PNG
# svg <- export_svg(pipeline)
# rsvg_png(charToRaw(svg), "pipeline_overview.png", width = 2200, height = 500)


# ============================================================
# 2009 DATASET
# ============================================================

zarr_path <- "/Users/justin.suca/Documents/TPruitt/Python Script/OpakapakaOutput.zarr"
ds <- xr$open_zarr(zarr_path)

# ============================================================
# 60–180 DAY CONNECTIVITY MATRIX
# ============================================================

release_lon <- as.numeric(ds$lon$isel(list(obs = 0L))$values)
release_lat <- as.numeric(ds$lat$isel(list(obs = 0L))$values)

py$ds_r <- ds

py_run_string("
import numpy as np

ever_settled = ds_r['settled'].max(dim='obs').values
settle_region_any = ds_r['settle_region'].max(dim='obs').values
age_at_settle = ds_r['age'].where(ds_r['settled'] == 1).min(dim='obs').values
")

settled <- as.integer(py$ever_settled)
settle_region <- as.integer(py$settle_region_any)
age_at_settle_days <- as.numeric(py$age_at_settle) / 86400

valid <- (settled == 1) &
  !is.na(age_at_settle_days) &
  age_at_settle_days >= 60 &
  age_at_settle_days <= 180

classify_release <- function(lon, lat) {
  if (lon >= -160.75 && lon <= -160.45 && lat >= 21.50 && lat <= 21.72) return("Kaula")
  if (lon >= -160.30 && lon <= -159.95 && lat >= 21.75 && lat <= 22.05) return("Niihau")
  if (lon >= -159.70 && lon <= -159.20 && lat >= 21.85 && lat <= 22.30) return("Kauai")
  if (lon >= -158.40 && lon <= -157.60 && lat >= 21.10 && lat <= 21.80) return("Oahu")
  if (lon >= -157.85 && lon <= -157.20 && lat >= 20.80 && lat <= 21.15) return("Penguin_Bank")
  if (lon >= -157.40 && lon <= -156.00 && lat >= 20.40 && lat <= 21.20) return("Maui_Nui")
  if (lon >= -156.00 && lon <= -154.70 && lat >= 18.90 && lat <= 20.30) return("Hawaii")
  return("Ocean")
}

release_name <- mapply(classify_release, release_lon, release_lat)

region_names <- c(
  "0" = "Kauai",
  "1" = "Niihau",
  "2" = "Kaula",
  "3" = "Oahu",
  "4" = "Penguin_Bank",
  "5" = "Maui_Nui",
  "6" = "Hawaii"
)

settle_name <- region_names[as.character(settle_region)]

release_valid <- release_name[valid]
settle_valid <- settle_name[valid]

connectivity_counts_60_180 <- table(release_valid, settle_valid)

connectivity_percent_60_180 <- round(
  prop.table(connectivity_counts_60_180, 1) * 100,
  2
)

print(connectivity_percent_60_180)

# ============================================================
# HEATMAP
# ============================================================

heat_df <- as.data.frame(as.table(connectivity_percent_60_180))
names(heat_df) <- c("Release", "Settlement", "Percent")

heat_df <- heat_df %>%
  dplyr::filter(Release != "Ocean", Settlement != "Ocean")

heat_df$Release <- gsub("_", " ", heat_df$Release)
heat_df$Settlement <- gsub("_", " ", heat_df$Settlement)

region_order <- c("Hawaii", "Maui Nui", "Penguin Bank", "Oahu", "Kauai", "Niihau", "Kaula")

heat_df$Release <- factor(heat_df$Release, levels = region_order)
heat_df$Settlement <- factor(heat_df$Settlement, levels = region_order)

heatmap_plot <- ggplot(heat_df, aes(x = Release, y = Settlement, fill = Percent)) +
  geom_tile(color = "white", linewidth = 0.8) +
  scale_fill_gradient(
    name = "Settlement (%)",
    low = "white",
    high = "#6A0DAD"
  ) +
  labs(
    title = "Larval Connectivity Between Hawaiian Islands",
    subtitle = "2009 Opakapaka (June–September, 30 m depth, 60–180 day PLD)",
    x = "Release Region",
    y = "Settlement Region"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    plot.title = element_text(size = 26, face = "bold"),
    plot.subtitle = element_text(size = 16),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 13),
    panel.grid = element_blank()
  )

print(heatmap_plot)

ggsave(
  "2009_connectivity_heatmap_clean.png",
  heatmap_plot,
  width = 10,
  height = 7,
  dpi = 300,
  bg = "white"
)




