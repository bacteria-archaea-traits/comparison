
###################
# Load data files #
###################

# Main data frame (see Madin et al. 2020)
df <- read.csv("data/Madin et al.2020/condensed_species_GTDB[NCBI_fill]_16102020.csv", as.is=TRUE)
df <- df[!is.na(df$species),]

#MIST transporter db (MIST3.0 database [mistdb.com])
ms <- read.csv("data/MIST/clean_mist01022019.csv", as.is=TRUE)

# Data from Lever et al. 2015
lv <- read_xlsx("data/Lever et al.2015/Lever_MicroEnergySupp_TableA2.xlsx", skip = 2)


#####################
# Layout parameters #
#####################

save_path <- "figs"

# Define fixed formats for all figures
text_size <- 9
text_color <- "black"
font_family <- "sans"
plot_line_color <- "black"
plot_line_width <- 0.5 #mm

# Define theme standards 
basic_layout <- 
  theme_bw() + 
  theme(
    panel.border = element_rect(color = plot_line_color, size = plot_line_width),
    panel.grid = element_blank(),
    axis.title = element_text(size = text_size, family=font_family, color=text_color),
    axis.text = element_text(size = text_size, family=font_family, color=text_color),
    axis.line = element_blank(),
    axis.ticks = element_line(color = plot_line_color, size = plot_line_width)
  )

#Define colour sets
colours_raw <- c("#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FB9A99","#E31A1C","#FDBF6F","#FF7F00")
#Rearrange colours
colours <- c(colours_raw[1],colours_raw[2],colours_raw[4],colours_raw[6],colours_raw[5],colours_raw[8],"#BDBDBD")

# Create special axes
growth_rate_axis_text <- expression(Max~growth~rate~(hr^{-1}))
growth_rate_norm_axis_text <- expression(Normalised~max~growth~rate~(hr^{-1}))
volume_axis_text <- expression(paste("Cell volume (",mu,m^3,")", sep = ""))


#######################
# Fix any data errors #
#######################

#Fix habitat type for "Halonatronum saccharophilum". Listed as 'soil' but should be 'sediment_hypersaline'
df$isolation_source[df$species == "Halonatronum saccharophilum"] <- "sediment_hypersaline"


####################################
# Merge with Mist transporter data #
####################################

# Merge and create new groups

df <- ms %>% group_by(species_tax_id) %>% 
  mutate(tcp_tot = sum(tcp.hk,tcp.hhk,tcp.rr,tcp.hrr,tcp.other, na.rm = TRUE), 
         tcp_hk_hkk_other_chemotaxis = sum(tcp.hk,tcp.hhk,tcp.other,tcp.chemotaxis, na.rm = TRUE),
         tcp_hk_hkk = sum(tcp.hk,tcp.hhk, na.rm = TRUE),
         hk_tot = sum(tcp.hk,tcp.hhk, na.rm = TRUE), 
         stp = sum(tcp_tot, ocp, na.rm = TRUE)) %>%
  select(species_tax_id,ocp, tcp_tot, tcp_hk_hkk_other_chemotaxis, tcp_hk_hkk,  stp, hk_tot, majormodes_total) %>% 
  right_join(df, by = "species_tax_id")
rm(ms)


######################################
# Calculate mid diameter and volumes #
######################################

df$d1_mid <- ifelse(!is.na(df$d1_up), (df$d1_lo + df$d1_up)/2, df$d1_lo)
df$d2_mid <- ifelse(!is.na(df$d2_up), (df$d2_lo + df$d2_up)/2, df$d2_lo)


###################
# Add growth rate #
###################

df$growth_rate <- log(2)/df$doubling_h


###########################
# Add Cobo Simon habitats #
###########################

en <- read.csv("data/Madin et al.2020/environments.csv", as.is=TRUE)

df <- df %>% left_join(en[,c("Type","Cobo.Simon.habitat")], by = c("isolation_source"="Type"))
rm(en)


###############################
# Create new habitat category #
###############################

#First letter upper case for environment
df$Cobo.Simon.habitat <- str_to_title(df$Cobo.Simon.habitat, locale = "en")

#Use Cobo Simon habitat scheme
df$habitat <- as.character(df$Cobo.Simon.habitat)

#Set any isolation_source that does not fit into the Cobo Simon scheme to "Other"
df$habitat[!is.na(df$isolation_source) & is.na(df$habitat)] <- "Other"

#Fix particular types of Cobo Simon habitats
df$habitat[df$habitat == "Rock"] <- "Other"
df$habitat[df$habitat == "Therm"] <- "Thermal"
df$habitat[grepl("plant|fungus|algae",df$isolation_source)] <- "Other"
df$habitat[grepl("feces|endotherm_surface|ectotherm_surface",df$isolation_source)] <- "Other"
df$habitat[df$isolation_source %in% c("host","host_animal")] <- "Other"

#Where habitat is not "Other" split host into 'endo' and 'ecto'
df$habitat[grepl("endotherm",df$isolation_source) & !(df$habitat == "Other")] <- "Endotherm"
df$habitat[grepl("ectotherm",df$isolation_source) & !(df$habitat == "Other")] <- "Ectotherm"

#Set all species with no isolation_source/habitat information to habitat = "Other"
df$habitat[is.na(df$habitat)] <- "Other"


####################################################
# Temperature adjusted growth rate using residuals #
####################################################

# With arrhenius temperature

df <- df %>% mutate(arrhenius_tmp = ifelse(!is.na(growth_tmp), 1/(growth_tmp + 273), NA))

#Get data set for creating regression model on growth rate against arrhenius growt temperature
sub <- df %>% filter(!is.na(growth_rate) & !is.na(growth_tmp))

#Get random 75% subset for training model (never use all data for model)
set.seed(125) 
sample = sample.split(sub$growth_rate, SplitRatio = .75)
train = subset(sub, sample == TRUE)

#Ordinary least squares using lm() 
model <- lm(log10(growth_rate)~arrhenius_tmp, data = train)

#Use model to calculate residuals in original data frame
#measured growth rate minus estimated growth rate = residual
df <- df %>% mutate(
  temp_adjusted_maxgrowth_arrhenius = ifelse(!is.na(growth_rate), log10(growth_rate) - (model$coefficients['arrhenius_tmp']*arrhenius_tmp+model$coefficients['(Intercept)']), NA)) 

# With growth temperature:

#Ordinary least squares using lm() 
model <- lm(log10(growth_rate)~growth_tmp, data = train)

#Use model to calculate residuals in original data frame
#measured growth rate minus estimated growth rate = residual
df <- df %>% mutate(
  temp_adjusted_maxgrowth = ifelse(!is.na(growth_rate), log10(growth_rate) - (model$coefficients['growth_tmp']*growth_tmp+model$coefficients['(Intercept)']), NA)) 

# Clean up
rm(train,model,sub)


###########################################
# Save full df for plots on original data #
###########################################

full_df <- df %>% ungroup()


####################################
# Create redundant data frame with #
# species in selected pathways     #
####################################

# This creates a redundant data frame 
# with each microbe potentially represented by multiple rows

#Extract all microbes in our data frame that fits into either of these pathways
selected_pathways <- c("hydrogen_oxidation_dark",	"cellulose_degradation",	"chitin_degradation",	"ligninolysis",	"xylan_degradation",	"sulfate_reduction",	"thiosulfate_reduction")

all <- df[0,]

for(i in 1:length(selected_pathways)) {
  
  pathway <- selected_pathways[i]

  sub <- df[grepl(sprintf("\\b%s\\b",pathway),df$pathways),]
  
  sub$pathway <- pathway

  all <- all %>% bind_rows(sub)
  print(pathway);
}

# Overwrite df with this
df <- all


##################################
# Create combined pathway groups #
##################################

df$pathway_final <- df$pathway

# Merger degradation groups
degraders <- c("cellulose_degradation","ligninolysis","chitin_degradation","xylan_degradation")
df$pathway_final[df$pathway %in% degraders] <- "cclx_degraders"

# Merge sulfate and thiosulfate reducers
df$pathway_final[df$pathway %in% c("sulfate_reduction","thiosulfate_reduction")] <- "Sulfate, Thiosulfate red."

#Ensure that each species is only included once in each pathway category 
df <- df %>% distinct(species, pathway_final, .keep_all = TRUE)
#(a species with methylotrophy and methanol oxidation would end up in the same category twice)

# Remove any pathway groups where n < 20
df <- df %>% group_by(pathway_final) %>% mutate(n = n()) %>% 
  filter(n >= 20) %>% 
  ungroup()

# Create a pathway name column (more readable names)
df$pathway_name <- NA
df$pathway_name <-str_to_sentence(gsub("_"," ",df$pathway_final))


############
# Finalise #
############

# Change factor levels of habitat for plots
df$habitat <- factor(df$habitat, levels = c("Fresh water","Fresh sediment","Marine water","Marine sediment","Soil","Thermal","Endotherm","Ectotherm","Intracellular","Other"))

# Convert genome size to Mbp
df <- df %>% mutate(genome_size = genome_size/1000000) 
full_df <- full_df %>% mutate(genome_size = genome_size/1000000) 
