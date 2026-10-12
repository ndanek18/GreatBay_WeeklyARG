library(tidyverse)
library(here)
library(readxl)
library(vegan)

#### ====================================================================== ####
#LOADING IN AND CONCATENTING .MAPSTAT.FILTERED AND .MAPSTAT FILES FOR EACH SAMPLE

#load in and concat all of the mapstat.filtered kma_panres (ARG hits) files that I downloaded from premise

input_directory <- "C:/Users/ndane/OneDrive - USNH/Documents/R/Harvey_Lab/ARGprofiler/raw_data/panres_files" # file path


list.files(input_directory, recursive = TRUE)

file_pattern <- ".mapstat.filtered"


panres_files <- list.files(input_directory,
                           pattern = file_pattern,
                           full.names = TRUE,
                           recursive = TRUE)


# Fix the name of the files to match the sample id
panres_file_names <- gsub(input_directory, "", panres_files) # replace input directory with blank space
panres_file_names <- gsub(file_pattern, "", panres_file_names) # replace file pattern with blank space
panres_file_names <- gsub("\\/", "", panres_file_names) # remove forward slashes, need to escape the forward 

names(panres_files) <- panres_file_names


panres_files # check they are all there

# Read in the panres files and concatenate them by sample id
panres_file_concat <- map_dfr(panres_files, 
                              ~read_tsv(.x,
                                        skip = 6,  col_types = cols(.default = "c")),  # read everything as character first # skip = number of hashtag lines (excluding column titles)
                              .id = "sampleid")

#infer the correct type per column after everything is combined as character
panres_file_concat <- panres_file_concat %>%
  type_convert()

#rename column 2 
names(panres_file_concat)[2] <- "RefSequence"


#now do the same for the .mapstat kma_mOTU files (bacteria hits)

input_directory_mOTUs <- "C:/Users/ndane/OneDrive - USNH/Documents/R/Harvey_Lab/ARGprofiler/raw_data/mOTU_files" # file path


list.files(input_directory_mOTUs, recursive = TRUE)

file_pattern2 <- "\\.mapstat$"


mOTU_files <- list.files(input_directory_mOTUs,
                           pattern = file_pattern2,
                           full.names = TRUE,
                           recursive = TRUE)


# Fix the name of the files to match the sample id
mOTU_file_names <- gsub(input_directory_mOTUs, "", mOTU_files) # replace input directory with blank space
mOTU_file_names <- gsub(file_pattern2, "", mOTU_file_names) # replace file pattern with blank space
mOTU_file_names <- gsub("^/", "", mOTU_file_names) # remove forward slashes, need to escape the forward 


names(mOTU_files) <- mOTU_file_names


mOTU_files # check they are all there

# Read in the files and concatenate them by sample id
mOTU_file_concat <- map_dfr(mOTU_files, 
                              ~read_tsv(.x,
                                        skip = 6,  col_types = cols(.default = "c")),  # read everything as character first # skip = number of hashtag lines (excluding column titles)
                              .id = "sampleid")


#check it worked
unique(mOTU_file_concat$sampleid)

#check number of rows per sample
table(mOTU_file_concat$sampleid)


#infer the correct type per column after everything is combined as character
mOTU_file_concat <- mOTU_file_concat %>%
  type_convert()

#rename column 2 
names(mOTU_file_concat)[2] <- "RefSequence"

# add "mOTU_" to all column names except sampleid to be able to differentiate it
#from the panres_concat file
colnames(mOTU_file_concat)[-1] <- paste0("mOTU_", colnames(mOTU_file_concat)[-1])


#### ====================================================================== ####
#LOADING IN METADATA

#load in metadata file that contains the PanRes info
PanRes_metadata <- read_xlsx("C:/Users/ndane/OneDrive - USNH/Documents/R/Harvey_Lab/ARGprofiler/panres_annotations.xlsx", 
                             sheet = 1, 
                             col_types = c("text", "text", "text", "text", "text", "numeric", "text", "text"))


#load in weekly metadata csv file (contains dates, sites, chlorophyll, etc for each sample id)
weekly_metadata <- read_csv("C:/Users/ndane/OneDrive - USNH/Documents/R/Harvey_Lab/ARGprofiler/raw_data/Metadata_ARG_metagenome_samples.csv")

#remove NAs
weekly_metadata <- weekly_metadata %>% 
  filter(if_any(everything(), ~ !is.na(.)))

#join all metadata files to the panres concatenated files 
all_joined <- panres_file_concat %>%
  left_join(weekly_metadata, by = "sampleid") %>% #left join in my weekly metadata
  left_join(PanRes_metadata, by = c("RefSequence" = "cluster_representative")) %>% #left join in the PanRes metadata
  left_join(panres_fragment_sums, by = "sampleid") %>%  #left join in calculated panres fragment count sums for each sample id
  left_join(mOTU_fragment_sums, by = "sampleid")#left join in mOTU sums to have bacteria data to standardize to 

#remove NAs
all_joined_NoNA <- all_joined[!is.na(all_joined$class), ]

#set the order of the locations
all_joined$Location <- factor(
  all_joined$Location,
  levels = c("CML", "JEL", "HILT", "NBL", "WIS", "LAMP", "MAIN"))

#name the months according to their number
all_joined <- all_joined %>%
  mutate(
    Month = factor(
      month_pull,
      levels = c(9, 10, 11, 12),
      labels = c("September", "October", "November", "December")))


#### ====================================================================== ####
#CALCULATIONS FOR STANDARDIZED COUNTS, MEANS, STANDARD ERROR 

#calculate panres (ARG) total fragment sum for each sampleid (produces one sum per sampleid)
panres_fragment_sums <- panres_file_concat %>%
  group_by(sampleid) %>%
  summarise(sum_panres_fragmentCount = sum(as.numeric(fragmentCount), na.rm = TRUE))


#calculate mOTU (bacterial) fragment sum for each sampleid
mOTU_fragment_sums <- mOTU_file_concat %>%
  group_by(sampleid) %>%
  summarise(sum_mOTU_fragmentCount = sum(as.numeric(mOTU_fragmentCount), na.rm = TRUE))


#calculate each ARG class fragment count / mOTU fragment sum for that sample id
#(it will become a new column in all_joined)
all_joined <- all_joined %>%
  ungroup() %>% 
  group_by(sampleid) %>% 
  mutate(fragmentCount_prop_eachARG = 
           as.numeric(fragmentCount) / as.numeric(sum_mOTU_fragmentCount)*10000) 
                                    #multiply by 10,000 so the ouput isn't so tiny


#calculate ARG fragment count sum / mOTU fragment count sum for each sampleid
all_joined <- all_joined %>%
  mutate(fragmentCount_prop_all = 
           as.numeric(sum_panres_fragmentCount) / as.numeric(sum_mOTU_fragmentCount)*10000)
                                    #multiply by 10,000 so the ouput isn't so tiny


#means and standard error for fragment counts for each class per month per site
ARG_summary <- all_joined %>%
  group_by(Location, Month, class) %>%
  summarise(
    mean_abundance = mean(fragmentCount_prop_eachARG, na.rm = TRUE),
    sd = sd(fragmentCount_prop_eachARG, na.rm = TRUE),
    n = sum(!is.na(fragmentCount_prop_eachARG)),
    SE = sd / sqrt(n),
    .groups = "drop")

#means and standard error per month and site
ARG_summary2 <- all_joined %>%
  group_by(Location, Month) %>%
  summarise(
    mean_frags = mean(fragmentCount_prop_eachARG, na.rm = TRUE),
    sd_frags = sd(fragmentCount_prop_eachARG, na.rm = TRUE),
    .groups = "drop")

#chlorophyll means per site and month 
chlorophyll_summary <- weekly_metadata %>%
  group_by(Location, month_pull) %>%
  summarise(
    mean_chlor = mean(chlorophyll_ug_per_L, na.rm = TRUE),
    sd_chlor = sd(chlorophyll_ug_per_L, na.rm = TRUE),
    .groups = "drop") %>%
  mutate(Month = factor(
      month_pull,
      levels = c(9, 10, 11, 12),
      labels = c("September", "October",
                 "November", "December"))) #name the months according to their number


#### ====================================================================== ####
#GRAPHS

#set color palletes for ARGs, sites, and months

ARG_palette <- c(
  "#B8B8E8", 
  "#A8D5BA", 
  "#F6C98B",  
  "#D9B8D8", 
  "#A6CEE3",  
  "#6FB5AA", 
  "#F4A582",  
  "#C6B4CE")


#bolder option if needed 
ARG_palette2 <- c(
  "#7A7AC7",  # purple-blue
  "#4FAF7B",  # green
  "#E6A23C",  # golden orange
  "#B66BB3",  # purple
  "#4C9BC6",  # blue
  "#35A89A",  # teal
  "#E76F51",  # coral
  "#8E6FAE"   # violet
)

#test the ARG colors 
ggplot(all_joined, aes(x = Location, y = fragmentCount_prop_eachARG, fill = class)) +
  geom_col(position = "stack") + #scale_y_continuous(transform = scales:: +
  labs( x = "Site", y = "Standardized fragment count", title = "ARG abundance by month", fill = "ARG Class") +
  all_bars + scale_fill_manual(values = ARG_palette)


Site_palette <- c(
  "LAMP" = "#332288",  
  "MAIN" = "#882255", 
  "WIS"  = "#D55E00",  
  "JEL"  = "#117733",  
  "NBL"  = "#E69F00", 
  "HILT" = "#CC79A7",  
  "CML"  = "#994F00")

#test the site colors
ggplot(pcoa_points, aes(x = PCoA1, y = PCoA2, color = Location)) +
  geom_point(size = 4) +
  theme_classic() +
  labs(
    x = paste0("PCoA Axis 1 (", PCoA1_percent, "%)"),
    y = paste0("PCoA Axis 2 (", PCoA2_percent, "%)"),
    color = "Location") + stat_ellipse() + theme_classic() +scale_color_manual(values = Site_palette) + xlim(-1.5, 1.5) + ylim(-2, 1.5)


Month_palette <- c(
  "October" = "#E69F00",  
  "November"   = "#009E73",  
  "September"  = "#56B4E9",  
  "December"  = "#332288")  

#test the month colors 
ggplot(pcoa_points_location, aes(x = PCoA1, y = PCoA2, color = Month)) +
  geom_point(size = 4) +
  theme_classic() +
  labs(
    x = paste0("PCoA Axis 1 (", PCoA1_percent, "%)"),
    y = paste0("PCoA Axis 2 (", PCoA2_percent, "%)"),
    color = "Month") + stat_ellipse() + theme_classic() + scale_color_manual(values = Month_palette) + xlim(-1.5, 1.5) + ylim(-2, 1.5)


#### ====================================================================== ####
#STACKED BAR GRAPHS 

#all_bars is the same for the first two stacked bar graphs, so can just add <+ all_bars> to the ggplot script 
all_bars <- list(geom_col(position = "stack"),
                       theme_classic(),
                       theme(legend.position = "right",
                             axis.text.x = element_text(angle = 45, hjust = 1)))


#stacked bar plot of site vs raw fragment count (unstandardized counts)
ggplot(all_joined, aes(x = Location, y = fragmentCount, fill = class)) +
  labs( x = "Site", y = "Fragment Count", title = "ARG abundance by month",
  fill = "ARG Class") + all_bars + scale_fill_manual(values = ARG_palette) + 
  facet_wrap(~Month) + 
  theme_classic() + theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.text = element_text(size = 8))

  
#stacked bar plot of site vs standardized fragment count(to mOTU sums)
ggplot(all_joined, aes(x = Location, y = fragmentCount_prop_eachARG, fill = class)) +
    geom_col(position = "stack") + 
    labs( x = "Site", y = "Standardized ARG fragment count", fill = "ARG Class") +
    all_bars + scale_fill_manual(values = ARG_palette) + facet_wrap(~Month) + 
    theme_classic() + theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.text = element_text(size = 8))
  
 
#stacked bar plot of Date vs fragment count
ggplot(all_joined, aes(x = Month, y = fragmentCount, fill = class)) +
    geom_col(position = "stack") +
    theme_classic() +
    theme(legend.position = "right",
          axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs( x = "Month", y = "Fragment Count",
          fill = "ARG Class") +
    facet_wrap("Location", nrow = 1, ncol = 7) + scale_fill_manual(values = ARG_palette)
  
 
#stacked bar plot of Date vs standardized fragment count
ggplot(all_joined, aes(x = Month, y = fragmentCount_prop_eachARG, fill = class)) +
    geom_col(position = "stack") +
    theme_classic() +
    theme(legend.position = "right",
          axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs( x = "Month", y = "Standardized Fragment Count", title = "ARG abundance over time by site",
          fill = "ARG Class") +
    facet_wrap("Location", nrow = 1, ncol = 7) + scale_fill_manual(values = ARG_palette)

#stacked bar plot of Date vs standardized fragment count
#work in progress
ggplot(all_joined, aes(x = Month, y = fragmentCount_prop_eachARG)) +
  geom_point() + geom_line() + 
  theme_classic() +
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs( x = "Month", y = "Standardized Fragment Count", title = "ARG abundance over time by site",
        fill = "ARG Class") +
  facet_wrap("Location", nrow = 7, ncol = 1) + scale_fill_manual(values = ARG_palette)

#### ====================================================================== ####
#scatterplots

#point graph of mean frag count faceted by site
ggplot(ARG_summary, aes(x = Month, y = mean_abundance, color = class, group = class)) +
  geom_point(size = 3) + 
  geom_errorbar(aes(ymin = mean_abundance - SE,ymax = mean_abundance + SE),width = 0.2) +
  facet_wrap(~Location, ncol = 1, nrow = 7) +
  theme_classic() + labs(y = "Mean standardized ARG fragment count", color= "ARG class") +
  scale_color_manual(values = ARG_palette) +
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.text = element_text(size = 7))


#mean frag count over time by site
ggplot(ARG_summary2, aes(x = Month, y = mean_frags, color = Location, group = Location)) +
  geom_point(size = 2) + labs(y = "Mean fragment count") + 
  geom_errorbar(aes(ymin = mean_frags - sd_frags,ymax = mean_frags + sd_frags),width = 0.2) +
  theme_classic() + scale_color_manual(values = Site_palette)
  
#### ====================================================================== ####
#CHLOROPHYLL GRAPHS 

#scatterplot of chlorophyll over time, colored by site 
ggplot(chlorophyll_summary, aes(x = Month, y = mean_chlor, color = Location, group = Location)) +
  geom_point(size = 3) + labs(y = "Mean chlorophyll(ug/L)") + 
  geom_errorbar(aes(ymin = mean_chlor - sd_chlor,ymax = mean_chlor + sd_chlor),width = 0.2) +
  theme_classic() + scale_color_manual(values = Site_palette)


#scatterplot of standardized fragment count vs log chlorophyll by location
ggplot( all_joined, aes(x = log(as.numeric(chlorophyll_ug_per_L)), y = as.numeric(fragmentCount_prop_all), color = Location)) + 
  geom_point(size = 3)  + 
  theme_classic() + labs(x = "log Chlorophyll (ug/L))", y = "Standardized Total ARG Fragment Count") + scale_color_manual(values = Site_palette)

#scatterplot of standardized fragment count vs log chlorophyll by class 
ggplot( all_joined, aes(x = log(as.numeric(chlorophyll_ug_per_L)), y = as.numeric(fragmentCount_prop_all), color = class)) + 
  geom_point(size = 2)  + 
  theme_classic() + labs(x = "log Chlorophyll (ug/L))", y = "Standardized Total ARG Fragment Count") + scale_color_manual(values = ARG_palette)
#+ geom_smooth(method = "lm", se = TRUE)


#### ====================================================================== ####
#Just CML, JEL data

#create a data table from all_joined for just CML and JEL
CML_JEL <- all_joined %>% 
  filter(Location %in% c("CML", "JEL")) %>% 
  mutate(
    Month = factor(
      month_pull,
      levels = c(9, 10, 11, 12),
      labels = c("September", "October", "November", "December"))) #name the months according to their number

#pull out CML and JEL from monthyl means data
CML_JEL_means <- ARG_summary2 %>% 
  filter(Location %in% c("CML", "JEL"))

# Calculate mean temperature for each location and month
temp_summary <- CML_JEL %>%
  group_by(Location, Month) %>%
  summarise(
    mean_temp = mean(as.numeric(temp), na.rm = TRUE),
    .groups = "drop")

# Scaling factor for the secondary axis
count_range <- range(as.numeric(CML_JEL$fragmentCount_prop_all), na.rm = TRUE)
temp_range <- range(temp_summary$mean_temp, na.rm = TRUE)

scale_factor <- diff(count_range) / diff(temp_range)

#plot fragment counts of all sites and temperature as 2nd y axis
ggplot(ARG_summary2, aes(x = Month)) +
  geom_point(
    aes(y = mean_frags, color = Location),
    size = 3) +
  geom_line(
    data = temp_summary,
    aes(
      y = (mean_temp - temp_range[1]) * scale_factor + count_range[1],
      color = Location,
      group = Location),linewidth = 1,linetype = "dashed") +
  scale_y_continuous(
    name = "Mean Standardized Fragment Count",
    sec.axis = sec_axis(
      ~ (. - count_range[1]) / scale_factor + temp_range[1],
      name = "Mean Temperature (C)")) +
  theme_classic() + labs(x = "Month",
    color = "Location") + scale_color_manual(values = Site_palette)


#plot fragment counts of just CML, JEL and temperature as 2nd y axis
ggplot(CML_JEL_means, aes(x = Month)) +
  geom_point(
    aes(y = mean_frags, color = Location),
    size = 3) +
  geom_line(
    data = temp_summary,
    aes(
      y = (mean_temp - temp_range[1]) * scale_factor + count_range[1],
      color = Location,
      group = Location),linewidth = 1,linetype = "dashed") +
  scale_y_continuous(
    name = "Mean Standardized Fragment Count",
    sec.axis = sec_axis(
      ~ (. - count_range[1]) / scale_factor + temp_range[1],
      name = "Mean Temperature (C)")) +
  theme_classic() + labs(x = "Month",
                         color = "Location") + scale_color_manual(values = Site_palette)


#plot of temperature over time at CML and JEL
ggplot(CML_JEL,
       aes(x = Month,
           y = as.numeric(temp),
           color = Location)) +
  geom_point(size = 2) +
  theme_classic() +
  labs(x = "Month", y = "Temperature") + scale_color_manual(values = Site_palette)

#scatterplot of standardized fragment count vs salinity with linear regression line 
ggplot(CML_JEL, aes(x = as.numeric(sal), y = as.numeric(fragmentCount_prop_all), color = Location)) +
  geom_point(size = 2) + geom_smooth(method = "lm", color = "black", se = TRUE) +
  theme_classic() + labs(x = "Salinity (ppt)", y = "Standardized Total Fragment Count") + 
  scale_color_manual(values = Site_palette)

#scatterplot of standardized fragment count vs temperature by site  with linear regression line 
ggplot(CML_JEL, aes(x = as.numeric(temp), y = as.numeric(fragmentCount_prop_all), color = Location)) + 
  geom_point(size = 2) + geom_smooth(method = "lm", color = "black", se = TRUE) +
  theme_classic() + labs(x = "Temperature (°C)", y = "Standardized Total Fragment Count") + 
  scale_color_manual(values = Site_palette)

#scatterplot of standardized fragment count vs temperature by class with linear regression line 
ggplot(CML_JEL, aes(x = as.numeric(temp), y = as.numeric(fragmentCount_prop_all), color = class)) + 
  geom_point(size = 2) + geom_smooth(method = "lm", se = TRUE) +
  theme_classic() + labs(x = "Temperature (°C)", y = "Standardized Total Fragment Count")  + 
  scale_color_manual(values = ARG_palette)


#### ====================================================================== ####
#PCoA plots

#by location 
#reshape data

arg_matrix <- all_joined %>%
  select(sampleid, RefSequence, fragmentCount_prop_eachARG) %>%
  pivot_wider(
    names_from = RefSequence,
    values_from = fragmentCount_prop_eachARG,
    values_fill = 0)


#values_fn = sum, #sum values within each ARG class  

# Remove sample ID from the abundance matrix
arg_abundance <- arg_matrix %>%
  column_to_rownames("sampleid")

# Calculate Bray-Curtis distance
arg_dist <- vegdist(arg_abundance, method = "bray")

# PCoA 
pcoa <- cmdscale(arg_dist, eig = TRUE, k = 2)

#get the PCoA coordinates
pcoa_points <- as.data.frame(pcoa$points)

colnames(pcoa_points) <- c("PCoA1", "PCoA2")

pcoa_points$sampleid <- rownames(pcoa_points)

#add in location and month data
pcoa_points <- pcoa_points %>%
  left_join(
    all_joined %>%
      mutate(sampleid = as.character(sampleid)) %>%
      select(sampleid, Location) %>%
      distinct(),
    by = "sampleid") %>%
  left_join(
    all_joined %>%
      mutate(sampleid = as.character(sampleid)) %>%
      select(sampleid, Month) %>%
      distinct(),
    by = "sampleid")


#calcualte % variation explained by each axis
eig_percent <- pcoa$eig / sum(pcoa$eig) * 100


PCoA1_percent <- round(eig_percent[1], 1)
PCoA2_percent <- round(eig_percent[2], 1)

#plot PCoA by Location
ggplot(pcoa_points, aes(x = PCoA1, y = PCoA2, color = Location)) +
  geom_point(size = 4) +
  theme_classic() +
  labs(
    x = paste0("PCoA Axis 1 (", PCoA1_percent, "%)"),
    y = paste0("PCoA Axis 2 (", PCoA2_percent, "%)"),
    color = "Location") + stat_ellipse() + theme_classic() + 
  scale_color_manual(values = Site_palette) + xlim(-1.5, 1.5) + ylim(-2, 1.5)

#PERMNOVA test for clustering
adonis2(arg_dist ~ Location, data = pcoa_points, permutations = 999)
#Location explained 64% of the variation in ARG composition 
#(PERMANOVA, pseudo-F = 9.185, R² = 0.64, p = 0.001, 999 permutations)

location_month <- adonis2(arg_dist ~ Location + Month, data = pcoa_points, permutations = 999, by = "terms")
location_month
#both are significant

#plot PCoA by month
ggplot(pcoa_points, aes(x = PCoA1, y = PCoA2, color = Month)) +
  geom_point(size = 4) +
  theme_classic() +
    labs(
    x = paste0("PCoA Axis 1 (", PCoA1_percent, "%)"),
    y = paste0("PCoA Axis 2 (", PCoA2_percent, "%)"),
    color = "Month") + stat_ellipse() + theme_classic() +
  scale_color_manual(values = Month_palette) + xlim(-1.5, 1.5) + ylim(-2, 1.5)

#PERMNOVA test for clustering
adonis2(arg_dist ~ Month, data = pcoa_points, permutations = 999)
#month explained 25.4% of the variation in ARG composition 
#(PERMANOVA, pseudo-F = 3.8639, R² = 0.254, p = 0.002, 999 permutations)

#PCoA for salty sites

#create a data table from all_joined for just CML, JEL, HILT, NBL
salty <- all_joined %>% 
  filter(Location %in% c("CML", "JEL", "HILT", "NBL"))

#reshape data
arg_matrix_salty <- salty %>%
  select(sampleid, RefSequence, fragmentCount_prop_eachARG) %>%
  pivot_wider(
    names_from = RefSequence,
    values_from = fragmentCount_prop_eachARG,
    values_fill = 0)

# Remove sample ID from the abundance matrix
arg_abundance_salty <- arg_matrix_salty %>%
  column_to_rownames("sampleid")

# Calculate Bray-Curtis distance
arg_dist_salty <- vegdist(arg_abundance_salty, method = "bray")

# PCoA 
pcoa_salty <- cmdscale(arg_dist_salty, eig = TRUE, k = 2)

#get the PCoA coordinates
pcoa_points_salty <- as.data.frame(pcoa_salty$points)

colnames(pcoa_points_salty) <- c("PCoA1", "PCoA2")


pcoa_points_salty$sampleid <- rownames(pcoa_points_salty)

#add in location and month data
pcoa_points_salty <- pcoa_points_salty %>%
  left_join(
    salty %>%
      mutate(sampleid = as.character(sampleid)) %>%
      select(sampleid, Location) %>%
      distinct(),
    by = "sampleid") %>%
  left_join(
    salty %>%
      mutate(sampleid = as.character(sampleid)) %>%
      select(sampleid, Month) %>%
      distinct(),
    by = "sampleid")


#calculate % variation explained by each axis
eig_percent_salty <- pcoa_salty$eig / sum(pcoa_salty$eig) * 100

PCoA1_percent_salty <- round(eig_percent_salty[1], 1)
PCoA2_percent_salty <- round(eig_percent_salty[2], 1)

#plot salty PCoA by Location
ggplot(pcoa_points_salty, aes(x = PCoA1, y = PCoA2, color = Location)) +
  geom_point(size = 4) +
  theme_classic() +
  labs(title = "Saltwater sites", 
    x = paste0("PCoA Axis 1 (", PCoA1_percent_salty, "%)"),
    y = paste0("PCoA Axis 2 (", PCoA2_percent_salty, "%)"),
    color = "Location") + stat_ellipse() + theme_classic() +scale_color_manual(values = Site_palette) + xlim(-2, 1.5) + ylim(-2, 1.5)

#plot salty pcoa by month
ggplot(pcoa_points_salty, aes(x = PCoA1, y = PCoA2, color = Month)) +
  geom_point(size = 4) +
  theme_classic() +
  labs(title= "Saltwater sites",
    x = paste0("PCoA Axis 1 (", PCoA1_percent_salty, "%)"),
    y = paste0("PCoA Axis 2 (", PCoA2_percent_salty, "%)"),
    color = "Month") + stat_ellipse() + theme_classic() +scale_color_manual(values = Month_palette) + xlim(-2, 1.5) + ylim(-2, 1.5)


#PCoA for fresh sites 

#create a data table from all_joined for just LAMP, WIS, MAIN
fresh <- all_joined %>% 
  filter(Location %in% c("LAMP", "WIS", "MAIN"))

#remove NAs
fresh <- fresh[!is.na(fresh$class), ]


#reshape data
arg_matrix_fresh <- fresh %>%
  select(sampleid, RefSequence, fragmentCount_prop_eachARG) %>%
  pivot_wider(
    names_from = RefSequence,
    values_from = fragmentCount_prop_eachARG,
    values_fill = 0)

# Remove sample ID from the abundance matrix
arg_abundance_fresh <- arg_matrix_fresh %>%
  column_to_rownames("sampleid")

# Calculate Bray-Curtis distance
arg_dist_fresh <- vegdist(arg_abundance_fresh, method = "bray")

# PCoA 
pcoa_fresh <- cmdscale(arg_dist_fresh, eig = TRUE, k = 2)

#get the PCoA coordinates
pcoa_points_fresh <- as.data.frame(pcoa_fresh$points)

colnames(pcoa_points_fresh) <- c("PCoA1", "PCoA2")


pcoa_points_fresh$sampleid <- rownames(pcoa_points_fresh)


#add in location and month data
pcoa_points_fresh <- pcoa_points_fresh %>%
  left_join(
    fresh %>%
      mutate(sampleid = as.character(sampleid)) %>%
      select(sampleid, Location) %>%
      distinct(),
    by = "sampleid") %>%
  left_join(
    fresh %>%
      mutate(sampleid = as.character(sampleid)) %>%
      select(sampleid, Month) %>%
      distinct(),
    by = "sampleid")


#calcualte % variation explained by each axis
eig_percent_fresh <- pcoa_fresh$eig / sum(pcoa_fresh$eig) * 100

PCoA1_percent_fresh <- round(eig_percent_fresh[1], 1)
PCoA2_percent_fresh <- round(eig_percent_fresh[2], 1)

#plot fresh PCoA by Location
ggplot(pcoa_points_fresh, aes(x = PCoA1, y = PCoA2, color = Location)) +
  geom_point(size = 4) +
  theme_classic() +
  labs(title= "Freshwater sites",
    x = paste0("PCoA Axis 1 (", PCoA1_percent_fresh, "%)"),
    y = paste0("PCoA Axis 2 (", PCoA2_percent_fresh, "%)"),
    color = "Location") + stat_ellipse() + theme_classic() +
  scale_color_manual(values = Site_palette) + xlim(-2, 1.5) + ylim(-2, 1.5)

#plot fresh pcoa by month
ggplot(pcoa_points_fresh, aes(x = PCoA1, y = PCoA2, color = Month)) +
  geom_point(size = 4) +
  theme_classic() +
  labs(title= "Freshwater sites",
    x = paste0("PCoA Axis 1 (", PCoA1_percent_fresh, "%)"),
    y = paste0("PCoA Axis 2 (", PCoA2_percent_fresh, "%)"),
    color = "Month") + stat_ellipse() + theme_classic() +
  scale_color_manual(values = Month_palette) + xlim(-2, 1.5) + ylim(-2, 1.5)

