
# Load data ---------------------------------------------------------------

femalegps <- read.csv("Data_inputs/WAAL_foraging_2010-2020_F.csv")
malegps <- read.csv("Data_inputs/WAAL_foraging_2010-2020_M.csv")

allgps <- rbind(femalegps, malegps)

# Individuals
allgps %>% group_by(sex) %>% summarise(n_distinct(id))

# sex    `n_distinct(id)`
# <chr>             <int>
#   1 Female              171
# 2 Male                175

# Number of trips
length(unique(allgps$trackid))

trip_summ <- allgps %>% group_by(id) %>% summarise(n_trips = n_distinct(trackid))
mean(trip_summ$n_trips); sd(trip_summ$n_trips)

# Duration of trips
mean(allgps$tripduration.days); sd(allgps$tripduration.days)

# Mean boldness -----------------------------------------------------------

mean(allgps$boldness, na.rm = T); sd(allgps$boldness, na.rm = T)

# Plot histogram ----------------------------------------------------------

## Set personality quantiles
upperQ <- 0.90 # 90
lowerQ <- 0.10 # 10

# Isolate individuals
boldness <- allgps %>% 
  dplyr::group_by(id) %>%
  dplyr::summarise(boldness = boldness[1], Sex = sex[1]) %>%
  filter(!is.na(boldness))

# Extract quantiles
quantiles <- boldness %>% 
  dplyr::group_by(Sex) %>% 
  dplyr::summarise(upperQuan = quantile(boldness, prob = upperQ),
                   lowerQuan = quantile(boldness, prob = lowerQ),
                   midQuan = quantile(boldness, prob = 0.5))



# Plot personality distribution -------------------------------------------

# Set colours
shy_col <- "#00DD2F"
bold_col <- "purple"

boldness$Sex <- factor(boldness$Sex, labels = c("Females", "Males"))

# Colour based on quantiles
lowQ <- data.frame(Sex = c("Females", "Males"), 
                   xmin = c(min(subset(boldness, Sex == "Females")$boldness), 
                            min(subset(boldness, Sex == "Males")$boldness)),
                   xmax = c(quantiles$lowerQuan[quantiles$Sex == "Female"], 
                            quantiles$lowerQuan[quantiles$Sex == "Male"]),
                   ymin = c(0, 0),
                   ymax = c(Inf, Inf))

uppQ <- data.frame(Sex = c("Females", "Males"), 
                   xmin = c(quantiles$upperQuan[quantiles$Sex == "Female"], 
                            quantiles$upperQuan[quantiles$Sex == "Male"]),
                   xmax = c(max(subset(boldness, Sex == "Females")$boldness),
                            max(subset(boldness, Sex == "Males")$boldness)),
                   ymin = c(0, 0),
                   ymax = c(Inf, Inf))

## Binning means shading might not match lower/upper limits; fix that here:
lowQ$xmin[1] <- lowQ$xmin[1]-0.05
uppQ$xmax[1] <- uppQ$xmax[1]+0.08

lowQ$xmin[2] <- lowQ$xmin[2]-0.04
uppQ$xmax[2] <- uppQ$xmax[2]+0.09

png("Figures/SX_boldness_histogram.png", width = 12, height = 6, units = "in", res = 600)
ggplot(aes(x = boldness), data = boldness) + 
  geom_histogram(bins = 30) + 
  facet_grid(.~Sex, labeller = label_value) + 
  theme_bw(base_size=15) + theme(strip.placement = "outside",
                                 strip.background = element_blank(),
                                 strip.text.x = element_text(size = 16)) +
  geom_rect(data = lowQ, 
            aes(x = NULL, y = NULL,
                xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), 
            alpha = 0.2, fill = shy_col) +
  geom_rect(data = uppQ, 
            aes(x = NULL, y = NULL,
                xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), 
            alpha = 0.2, fill = bold_col) +
  xlab("Boldness Estimate") + ylab("Frequency")
dev.off()
