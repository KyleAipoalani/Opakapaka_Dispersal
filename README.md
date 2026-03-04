# Opakapaka Larval Dispersal Model

## Overview
This project models the larval dispersal of *Opakapaka* (Pristipomoides filamentosus), a deep-water snapper species found throughout the Hawaiian Islands. Using ocean current data from the Regional Ocean Modeling System (ROMS), simulated larvae are released from known habitat locations and tracked through time to simulate how ocean circulation could transport them between islands and potential settlement areas.

## Data Sources
This project integrates several datasets to simulate Opakapaka larval dispersal and potential settlement:
- **ROMS Ocean Current Data**: Output from the Regional Ocean Modeling System (ROMS) used to drive particle advection and simulate larval transport through ocean circulation.
- **ETOPO 2022 (15 arc-second Bedrock Bathymetry)**: Global bathymetry dataset from NOAA used to determine seafloor depth and define potential settlement habitat ranges used in the model.
- **Parcels Particle Tracking Framework**: The OceanParcels framework is used to simulate particle trajectories and model larval movement through time.

## Dependancies
- **parcels**: Used to simulate larval particles moving through ocean currents.
- **numpy**: Used for numerical calculations and working with arrays of model data.
- **matplotlib**: Used to create figures and visualize bathymetry, particle tracks, and settlement areas.
- **rasterio**: Used to read bathymetry data from GeoTIFF files.
- **geopandas**: Used to work with spatial data and generate settlement regions.
- **shapely**: Used for geometric operations such as buffering habitat areas and creating spatial boundaries.
