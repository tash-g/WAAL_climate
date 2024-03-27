### Plot maps for Figure 1 

# NB. LEGACY CODE USING SP PACKAGES ============================================

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

gps_waal <- read.csv("Data_inputs/WAAL_gpsLocations_1989-2020.csv")
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



# Get polar front ---------------------------------------------------------

apf_sf <- st_read(file.path("Data_inputs/GIS/shapefile/antarctic_circumpolar_current_fronts.shp"))
apf_sf <- apf_sf[apf_sf$NAME == "Polar Front (PF)",]
plot(st_geometry(apf_sf), axes = T)

## Change the CRS
apf_sf2 <- sf::st_transform(apf_sf, crs = proj.laea)


# Create the plots -------------------------------------------------------------

## FEMALE ##

annual_F <- subset(kdareas_annual, id == "F")

annual_kernelPlot_F <-
  # Set up the base map
  ggplot() + 
  geom_sf(data = world2, fill = "cadetblue", colour = "grey") +
  # Add the polar front
  geom_sf(data = apf_sf2, col = "#CF00BE", linewidth = 1, linetype = "dashed") +
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
  ggplot() + 
  geom_sf(data = world2, fill = "cadetblue", colour = "grey") +
  # Add the polar front
  geom_sf(data = apf_sf2, col = "#CF00BE", linewidth = 1, linetype = "dashed") +
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




