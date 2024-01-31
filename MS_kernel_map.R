### Plot maps for Figure 1 

# Load functions & packages ----------------------------------------------------

# getLevel function for kernel density estimation
getLevel <- function(x,y,prob) { 
  require(MASS)
  kk <- MASS::kde2d(x,y) # Compute the 2D kernel density estimate
  dx <- diff(kk$x[1:2]) # Calculate the differences in x and y values to determine the grid cell sizes
  dy <- diff(kk$y[1:2])  
  sz <- sort(kk$z) # Sort the density values in ascending order
  c1 <- cumsum(sz) * dx * dy # Calculate the cumulative sum of the sorted density values multiplied by grid cell areas
  approx(c1, sz, xout = 1 - prob)$y  # Use linear interpolation to estimate the contour level corresponding to the specified probability
}

# Define packages
packages <- c("ggplot2", "raster", "dplyr", "data.table", "rnaturalearth", 
              "rnaturalearthdata", "ggpubr", "sf", "readxl", "magrittr", "MASS")

# Install packages not yet installed - change lib to library path
# installed_packages <- packages %in% rownames(installed.packages())
# 
# if (any(installed_packages == FALSE)) {
#   install.packages(packages[!installed_packages])
# }

# Load packages
invisible(lapply(packages, library, character.only = TRUE))

# Load the data ----------------------------------------------------------------

gps_waal <- read.csv("Data_inputs/WAAL_gpsLocations_1989-2020.csv")
gps_waal <- subset(gps_waal, !is.na(DateTime))
gps_waal$DateTime <- as.POSIXct(gps_waal$DateTime, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")

# Isolate period 2010 - 2020
gps_waal <- subset(gps_waal, Year >= 2010)

# Get months
gps_waal$month <- as.numeric(format(gps_waal$DateTime, "%m"))

# Make year a character
gps_waal$Year <- as.character(gps_waal$Year)
gps_waal$RingYr <- paste0(gps_waal$Ring, gps_waal$Year, sep = "_")


## Set projections, project world and colony -----------------------------------

# Set projections Lambert Azimuthal Equal Area projection (proj = laea)
proj.laea <- "+proj=laea +lat_0=-46.358639 +lon_0=51.706972 +x_0=0 +y_0=0 +datum=WGS84 
              +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0" 
proj.dec <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

# Get a base world plot
world <- ne_countries(scale = "medium", returnclass = "sf")
world2 <- sf::st_transform(world, crs = proj.laea)

# Set colony location for plotting 
colony_sf <- st_as_sf(data.frame(Lon = 51.706972, Lat = -46.358639), coords = c("Lon", "Lat"))

# Set the CRS for the colony sf object
st_crs(colony_sf) <- st_crs(proj.dec) 
colony_transformed <- st_transform(colony_sf, crs = proj.laea) 
colony_transformed$Lon <- st_coordinates(colony_transformed)[, 1]
colony_transformed$Lat <- st_coordinates(colony_transformed)[, 2]


## Project the GPS data --------------------------------------------------------

gps_waal %<>% dplyr::select(c(Ring, Sex, RingYr, Longitude, Latitude))
gps_waal.sf <- st_as_sf(gps_waal, coords = c("Longitude", "Latitude"))
st_crs(gps_waal.sf) <- st_crs(proj.dec) 
gps_waal.sf.transformed <- st_transform(gps_waal.sf, crs = proj.laea) 

# Get coordinates
gps_waal.sf.transformed$Longitude <- st_coordinates(gps_waal.sf.transformed)[, 1]
gps_waal.sf.transformed$Latitude <- st_coordinates(gps_waal.sf.transformed)[, 2]
gps_waal.sf.df <- as.data.frame(gps_waal.sf.transformed)

## Separate males and females
gps_waal_f <- subset(gps_waal.sf.transformed, Sex == "F")
gps_waal_m <- subset(gps_waal.sf.transformed, Sex == "M")


# Calculate kernel density estimates -------------------------------------------

## Females ##
L90_F <- getLevel(gps_waal_f$Longitude, gps_waal_f$Latitude, 0.9)
L75_F <- getLevel(gps_waal_f$Longitude, gps_waal_f$Latitude, 0.75)
L50_F <- getLevel(gps_waal_f$Longitude, gps_waal_f$Latitude, 0.5)

kk <- MASS::kde2d(gps_waal_f$Longitude, gps_waal_f$Latitude)
dimnames(kk$z) <- list(kk$x, kk$y)
# Transform data for plotting
dc_F <- melt(kk$z)


## Males ##
L90_M <- getLevel(gps_waal_m$Longitude, gps_waal_m$Latitude, 0.9)
L75_M <- getLevel(gps_waal_m$Longitude, gps_waal_m$Latitude, 0.75)
L50_M <- getLevel(gps_waal_m$Longitude, gps_waal_m$Latitude, 0.5)

kk <- MASS::kde2d(gps_waal_m$Longitude, gps_waal_m$Latitude)
dimnames(kk$z) <- list(kk$x, kk$y)
# Transform data for plotting
dc_M <- melt(kk$z)


# Create the plots -------------------------------------------------------------

## Females ##
plot_female <-
  ggplot(data = world2) + 
  # Build the base world plot
  geom_sf(fill = "cadetblue", colour = "grey") +
  ggspatial::annotation_scale(location = "bl", width_hint = 0.25, style = "bar") +
  coord_sf(crs = proj.laea, xlim = c(-3500000, 2500000), ylim = c(-3500000, 2500000),
           label_axes = list(top = "E", left = "N", bottom = "E", right = "N")) +
  # Set title and subtitle
  ggtitle("Females") +
  # Add female contour lines
  geom_path(data = gps_waal_f, aes(x = Longitude, y = Latitude, group = RingYr), alpha = 0.2, col = "grey") + 
  geom_contour(data = dc_F, aes(x = Var1, y = Var2, z = value), breaks = L50_F, colour = "#FDD835", linewidth = 1) + 
  geom_contour(data = dc_F, aes(x = Var1, y = Var2, z = value), breaks = L75_F, colour = "#FBC02D", linewidth = 1) + 
  geom_contour(data = dc_F, aes(x = Var1, y = Var2, z = value), breaks = L90_F, colour = "#FF8F00", linewidth = 1) +
  # Location of Crozet
  annotate("point", shape = 17, (colony_transformed$Lon + 100), colony_transformed$Lat) +
  annotate("text", label = "Crozet", colony_transformed$Lon, (colony_transformed$Lat - 100000)) +
  # Plot theme
  scale_x_continuous(breaks = seq(-180, 180, by = 10)) +
  scale_y_continuous(breaks = seq(180, -180, by = -10)) + 
  theme(panel.background = element_rect(fill = "white"), 
        panel.grid.major = element_line(colour = "grey80"),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text.x.bottom = element_blank(), 
        axis.title.x.bottom = element_blank(),
        axis.text.y.right = element_blank(), 
        axis.title.y.right = element_blank(),
        axis.title.y.left = element_blank(),
        plot.title = element_text(size = 18, face = "bold"),
        plot.subtitle = element_text(size = 14),
        legend.position = "none") 


## Males ##

plot_male <-
  ggplot(data = world2) + 
  # Build the base world plot
  geom_sf(fill = "cadetblue", colour = "grey") +
  ggspatial::annotation_scale(location = "bl", width_hint = 0.25, style = "bar") +
  coord_sf(crs = proj.laea, xlim = c(-3500000, 2500000), ylim = c(-3500000, 2500000),
           label_axes = list(top = "E", left = "N", bottom = "E", right = "N")) +
  # Set title and subtitle
  ggtitle("Males") +
  # Add male contour lines
  geom_path(data = gps_waal_m, aes(x = Longitude, y = Latitude, group = RingYr), alpha = 0.2, col = "grey") + 
  geom_contour(data = dc_M, aes(x = Var1, y = Var2, z = value), breaks = L50_M, 
               colour = "#81D4FA", linewidth = 1) + 
  geom_contour(data = dc_M, aes(x = Var1, y = Var2, z = value), breaks = L75_M, 
               colour = "#1E88E5", linewidth = 1) + 
  geom_contour(data = dc_M, aes(x = Var1, y = Var2, z = value), breaks = L90_M, 
               colour = "#0D47A1", linewidth = 1) +
  # Location of Crozet
  annotate("point", shape = 17, (colony_transformed$Lon + 100), colony_transformed$Lat) +
  annotate("text", label = "Crozet", colony_transformed$Lon, (colony_transformed$Lat - 100000)) +
  # Plot theme
  scale_x_continuous(breaks = seq(-180, 180, by = 10)) +
  scale_y_continuous(breaks = seq(180, -180, by = -10)) + 
  theme(panel.background = element_rect(fill = "white"), 
        panel.grid.major = element_line(colour = "grey80"),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text.x.bottom = element_blank(), 
        axis.title.x.bottom = element_blank(),
        axis.text.y.right = element_blank(), 
        axis.title.y.right = element_blank(),
        axis.title.y.left = element_blank(),
        plot.title = element_text(size = 18, face = "bold"),
        plot.subtitle = element_text(size = 14),
        legend.position = "none") 



# ~ # ANNUAL UD FIGURE # ~ ------------------------------------------------

png("Figures/male_female_UD.png", width = 12, height = 8, units = "in", res = 300)
ggpubr::ggarrange(plot_female, plot_male,
                  ncol = 2,
                  align = "hv")
dev.off()





# ADEHABITAT - LEGACY CODE =====================================================

# Load packages ----------------------------------------------------------------

# Define packages
packages <- c("ggplot2", "raster", "dplyr", "data.table", "rnaturalearth", 
              "rnaturalearthdata", "ggpubr", "sp", "adehabitatHR", "sf", "readxl",
              "polyggon")

# Install packages not yet installed - change lib to library path
# installed_packages <- packages %in% rownames(installed.packages())
# 
# if (any(installed_packages == FALSE)) {
#   install.packages(packages[!installed_packages])
# }

# Load packages
invisible(lapply(packages, library, character.only = TRUE))

# Set colours
female_col <- "#FFC20A"
male_col <- "#0C7BDC"

#  Set UD sizes
UDs <- c(50, 75, 90)


# Load the data ----------------------------------------------------------------

gps_waal <- read.csv("Data_original/WAAL_GPS_1990-2020.csv")
gps_waal <- subset(gps_waal, !is.na(DateTime))
gps_waal$DateTime <- as.POSIXct(gps_waal$DateTime, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")

# Isolate period after 2010
gps_waal <- subset(gps_waal, Year >= 2010)

# Subset to relevant columns
gps_waal <- gps_waal[,c("Ring", "DateTime", "Sex", "Year", "Longitude", "Latitude"),]

# Get months
gps_waal$month <- as.numeric(format(gps_waal$DateTime, "%m"))

# Make year a character
gps_waal$Year <- as.character(gps_waal$Year)


# Set projections, project world and colony ------------------------------------

# Set projections 
proj.laea <- "+proj=laea +lat_0=-46.358639 +lon_0=51.706972 +x_0=0 +y_0=0 +datum=WGS84 
              +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0" 
proj.dec <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

# Set and project colony location for plotting 
colony <- data.frame(Lon = 51.706972, Lat = -46.358639)
coordinates(colony) <- ~Lon+Lat
proj4string(colony) <- proj.dec
colony.sp <- spTransform(colony, CRS(proj.laea))
colony.df <- as.data.frame(colony.sp)
names(colony.df) <- c("Longitude", "Latitude")

# Get a base world plot
world <- ne_countries(scale = "medium", returnclass = "sf")
world2 <- sf::st_transform(world, crs = proj.laea)

## Project the data ------------------------------------------------------------

# Data for kernel estimation

# Make a spatial points dataframe of dataset for kernels
mydat.sp <- gps_waal[,c("Sex", "Longitude", "Latitude")]

# Set coordinates and project
coordinates(mydat.sp) <- ~Longitude + Latitude
proj4string(mydat.sp) <- proj.dec
mydat.sp <- spTransform(mydat.sp, CRS(proj.laea))
mydat.df <- as.data.frame(mydat.sp)


# Data for plotting (incudes ring)
gps_waal.proj <- gps_waal
coordinates(gps_waal.proj) <- ~ Longitude+Latitude
proj4string(gps_waal.proj) <- proj.dec

gps_waal.sp <- spTransform(gps_waal.proj, CRS(proj.laea))
gps_waal.df <- as.data.frame(gps_waal.sp)
names(gps_waal.df)[c(6,7)] <- c("Longitude", "Latitude")

# Estimate each kernel ---------------------------------------------------------

# Estimate utilisation distributions for each sex
kud_annual <- kernelUD(mydat.sp[,1], grid = 400, same4all = TRUE)

# Extract vertices for each UD
kdareas_list <- list()
sexes <- c("F", "M")

for (i in 1:length(UDs)) {
  
  kareas_ud <- getverticeshr(kud_annual, UDs[i])
  
  # Fortify the raster for plotting
  kdareas_ud <- fortify(kareas_ud)
  kdareas_ud$UD <- UDs[i]
  
  kdareas_list[[i]] <- kdareas_ud
  
}

# Combine the UD datasets
kdareas_annual <- do.call("rbind", kdareas_list)

#save(kdareas_annual, file = "annual_kernels.RData")
load("Data_outputs/annual_kernels.RData")


# Create the plots -------------------------------------------------------------

## FEMALE ##

annual_F <- subset(kdareas_annual, id == "F")

annual_kernelPlot_F <-
  # Set up the base map
  ggplot(data = world2) + 
  geom_sf(fill = "cadetblue", colour = "grey") +
  ggspatial::annotation_scale(location = "bl", width_hint = 0.25, style = "bar") +
  coord_sf(crs = proj.laea, xlim = c(-3500000, 2500000), ylim = c(-3500000, 2500000),
            label_axes = list(top = "E", left = "N", bottom = "E", right = "N")) +
  # Add the GPS tracks
  geom_path(aes(x = Longitude, y = Latitude, group = Ring), 
            alpha = 0.35, size = 0.25, col = "black",
            dat = subset(gps_waal.df, Sex == "F")) +
  # Plot the kernels
  geom_polygon(aes(x = long, y = lat, group = group), 
                         size = 0.8, fill = female_col, col = "#ff800a",
                         dat = subset(annual_F, UD == 50),
                         alpha = .3) +
  geom_polygon(aes(x = long, y = lat, group = group), 
                         size = 0.8, fill = female_col, col = "#ff800a",
                         dat = subset(annual_F, UD == 75),
                         alpha = .3) +
  geom_polygon(aes(x = long, y = lat, group = group), 
                         size = 0.8, fill = female_col, col = "#ff800a",
                         dat = subset(annual_F, UD == 90),
                         alpha = .3) +
  # Add point for Crozet
  annotate("point", shape = 17, (colony.df$Lon + 100), colony.df$Lat) +
  annotate("text", label = "Crozet", colony.df$Lon, (colony.df$Lat - 100000)) +
  # Set plot parameters
  scale_x_continuous(breaks = seq(-180, 180, by = 10)) +
  scale_y_continuous(breaks = seq(180, -180, by = -10)) + 
  scale_fill_manual(values = female_col) +
  ggtitle("Females") +
  theme(panel.background = element_rect(fill = "white"), 
        panel.grid.major = element_line(colour = "grey80"),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text.x.bottom = element_blank(), 
        axis.title.x.bottom = element_blank(),
        axis.text.y.right = element_blank(), 
        axis.title.y.right = element_blank(),
        axis.title.y.left = element_blank(),
        plot.title = element_text(face = "bold"),
        legend.position = "none") 

  
## Males ##

annual_M <- subset(kdareas_annual, id == "M")

annual_kernelPlot_M <-
  # Set up the base map
  ggplot(data = world2) + 
  geom_sf(fill = "cadetblue", colour = "grey") +
  ggspatial::annotation_scale(location = "bl", width_hint = 0.25, style = "bar") +
  coord_sf(crs = proj.laea, xlim = c(-3500000, 2500000), ylim = c(-3500000, 2500000),
           label_axes = list(top = "E", left = "N", bottom = "E", right = "N")) +
  # Add the GPS tracks
  geom_path(aes(x = Longitude, y = Latitude, group = Ring), 
            alpha = 0.35, size = 0.25, col = "black",
            dat = subset(gps_waal.df, Sex == "M")) +
  # Plot the kernels
  geom_polygon(aes(x = long, y = lat, group = group), 
               size = 0.8, fill = male_col, col = "#2b0cdc",
               dat = subset(annual_M, UD == 50),
               alpha = .3) +
  geom_polygon(aes(x = long, y = lat, group = group), 
               size = 0.8, fill = male_col, col = "#2b0cdc",
               dat = subset(annual_M, UD == 75),
               alpha = .3) +
  geom_polygon(aes(x = long, y = lat, group = group), 
               size = 0.8, fill = male_col, col = "#2b0cdc",
               dat = subset(annual_M, UD == 90),
               alpha = .3) +
  # Add point for Crozet
  annotate("point", shape = 17, (colony.df$Lon + 100), colony.df$Lat) +
  annotate("text", label = "Crozet", colony.df$Lon, (colony.df$Lat - 100000)) +
  # Set plot parameters
  scale_x_continuous(breaks = seq(-180, 180, by = 10)) +
  scale_y_continuous(breaks = seq(180, -180, by = -10)) + 
  scale_fill_manual(values = female_col) +
  ggtitle("Males") +
  theme(panel.background = element_rect(fill = "white"), 
        panel.grid.major = element_line(colour = "grey80"),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text.x.bottom = element_blank(), 
        axis.title.x.bottom = element_blank(),
        axis.text.y.right = element_blank(), 
        axis.title.y.right = element_blank(),
        axis.title.y.left = element_blank(),
        plot.title = element_text(face = "bold"),
        legend.position = "none") 


# FIGURE 1: Breeding UD ========================================================
png("Figures/FIG1_male_female_UD_legacy.png", 
    width = 12, height = 6, units = "in", res = 300)
ggpubr::ggarrange(annual_kernelPlot_F, annual_kernelPlot_M,
                  ncol = 2, nrow = 1,
                  widths = c(1, 0.92))
dev.off()




