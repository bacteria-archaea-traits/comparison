
############
# Figure 1 #
############

# Transform data according to data type
# Note: temp_adjusted_maxgrowth is the residuals of the log10 transformed growth rate against growth temperature (see prep.R),
# and so should not be transformed again here.

org <- full_df %>%
  filter(!is.na(genome_size) & 
           !is.na(d1_mid) & 
           !is.na(rRNA16S_genes) & 
           !is.na(growth_tmp) & 
           !is.na(stp),
         !is.na(ocp)) %>%
  mutate("Genome size" = log10(genome_size), 
         D1 = log10(d1_mid),
         RRN = sqrt(rRNA16S_genes),
         STP = sqrt(stp),
         "HK tot" = sqrt(hk_tot), 
         "Growth tmp" = growth_tmp,
         "Max growth rate" = temp_adjusted_maxgrowth)

org <- org %>% select(species, `Genome size`, `HK tot`, `Growth tmp`, `Max growth rate`, D1, RRN, STP)

# Restrict data to rows with max growth rate
dat <- org %>% filter(!is.na(`Max growth rate`)) %>% 
  select(-species)

# Chose PC axes to output
plot_pcs <- c("PC1","PC2")

# Set rownames as number values
rownames(dat) <- seq(1,nrow(dat),1)

# Run PCA function to generate data for plots
pca <- pca_format(dat)

# Get % and vectors for specified PC axes (above)
props <- pca_prop_explained(dat,plot_pcs)
vectors<- pca_vectors(dat,plot_pcs)

# Plot data
p <- ggplot(pca, aes_string(x = plot_pcs[1], y = plot_pcs[2])) + 
  geom_point(colour = colours[2], alpha = 0.3, size = 1) +
  scale_x_continuous(limits = c(-5,5)) +
  scale_y_continuous(limits = c(-3.5,3.5)) +
  basic_layout +
  theme(
    legend.justification = c(0,1),
    legend.text = element_text(size = 8),
    legend.key.height = unit(0.3, "cm"),
    legend.spacing.y = unit(0.2, "cm")
  ) +
  labs(x = paste0(sprintf("%s (%s",plot_pcs[1], props[1]),"%)"), y = paste0(sprintf("%s (%s",plot_pcs[2], props[2]),"%)"),
       title = "")

# Add PCA vectors
p <- p + 
  coord_equal() + 
  geom_segment(data = vectors, aes(x = 0, y = 0, xend = v1, yend = v2), arrow = arrow(length = unit(0.1,"cm")), size = 0.3, alpha = 1, color = "black")
p

# Save version without vector labels
#ggsave(filename = "pca_all_with_tmp.adj.max.growth[PC1,PC2]_noVNames.png", plot = p, device = "png", path = save_path, units = "cm", width = 8, height = 6, dpi = 600, limitsize = TRUE)

# Add vector names
p <- p + 
  geom_text(data = vectors, aes(x = v1*1.1, y = v2*1.1, label = varnames), size = 2, vjust = 0, hjust = 0, alpha = 1, colour = "black")
p

# Save complete version
ggsave(filename = "Fig1.png", plot = p, device = "png", path = save_path, units = "cm", width = 8, height = 6, dpi = 600, limitsize = TRUE)

#############
# Figure S1 #
#############

# Use PCA data from above (Fig 1)

# Chose PC axes to output
plot_pcs <- c("PC1","PC3")

# Get % and vectors for specified PC axes (above)
props <- pca_prop_explained(dat,plot_pcs)
vectors<- pca_vectors(dat,plot_pcs)

#PCA plot layout
p <- ggplot(pca, aes_string(x = plot_pcs[1], y = plot_pcs[2])) + 
  geom_point(colour = colours[2], alpha = 0.3, size = 1) +
  scale_x_continuous(limits = c(-5,5)) +
  #scale_y_continuous(limits = c(-5.8,4)) +
  basic_layout +
  theme(
    legend.justification = c(0,1),
    legend.text = element_text(size = 8),
    legend.key.height = unit(0.3, "cm"),
    legend.spacing.y = unit(0.2, "cm")
  ) +
  labs(x = paste0(sprintf("%s (%s",plot_pcs[1], props[1]),"%)"), y = paste0(sprintf("%s (%s",plot_pcs[2], props[2]),"%)"),
       title = "")
#PCA vector layout
p <- p + 
  coord_equal() + 
  geom_segment(data = vectors, aes(x = 0, y = 0, xend = v1, yend = v2), arrow = arrow(length = unit(0.1,"cm")), size = 0.3, alpha = 1, color = "black")
p

#Save version without vector labels
#ggsave(filename = "pca_all_with_tmp.adj.max.growth[PC1,PC3]_noVNames.png", plot = p, device = "png", path = save_path, units = "cm", width = 8, height = 6, dpi = 600, limitsize = TRUE)

p <- p + 
  geom_text(data = vectors, aes(x = v1*1.1, y = v2*1.1, label = varnames), size = 2, vjust = 0, hjust = 0.3, alpha = 1, colour = "black")
p

#Save complete version
ggsave(filename = "FigS1.png", plot = p, device = "png", path = save_path, units = "cm", width = 8, height = 6, dpi = 600, limitsize = TRUE)


############
# Figure 2 #
############

sub <- full_df %>% filter(!is.na(genome_size) & !is.na(majormodes_total))

colours_tcp <- c(colours_raw[8],colours_raw[5],colours_raw[1],colours_raw[2],"#BDBDBD",colours_raw[4],colours_raw[6])

sub2 <- sub %>% filter(!is.na(tcp_hk_hkk))
sub2_no_other <- sub2 %>% filter(habitat == "Other")

p <- ggplot(sub2_no_other, aes(x = genome_size, y = tcp_hk_hkk)) + 
  geom_jitter(aes(fill = habitat), alpha = 0.7, colour = "#656565", size = 2, shape = 21, height = 0.03) + 
  geom_jitter(data = sub2 %>% filter(!(habitat == "Other")),  aes(x = genome_size, y = tcp_hk_hkk, fill = habitat), alpha = 0.7, colour = "#656565", size = 2, shape = 21, height = 0.03) + 
  geom_abline(intercept = log10(10), slope = 1) +
  geom_abline(intercept = log10(2), slope = 2, col = "red") +
  scale_fill_manual(values = colours_tcp) + 
  scale_x_log10() + 
  scale_y_log10() + 
  annotation_logticks(sides="lb", short = unit(1,"mm"), mid = unit(1,"mm"), long = unit(1,"mm")) +
  basic_layout + 
  theme(
    legend.justification = c(0,1),
    legend.text = element_text(size = 8),
    legend.key.height = unit(0.3, "cm"),
    legend.spacing.y = unit(0.2, "cm")
  )+
  #guides(fill = guide_legend(override.aes = list(shape = 21))) +
  labs(x = "Genome size (Mb)", y = "hk/hhk sensor genes", fill = "Habitat")
p

ggsave(filename = "Fig2.png", plot = p, device = "png", path = save_path, units = "cm", width = 13.6, height = 10, dpi = 600, limitsize = TRUE)


###############
# Figure B1.1 #
###############

# Data from Lever et al. 2015 in 'lv'

#Rename columns
names(lv) <- c("organism","width","length","shape_factor")
#Make starved/non-starved column
lv$starved <- ifelse(grepl("non-starved",lv$organism), "no", "yes")
#Remove starvation information from names
lv$organism <- word(lv$organism, 1, 2)

#Split into starved and non-starved
yes <- lv[lv$starved == "yes",]
yes <- subset(yes, select = -starved)
no <- lv[lv$starved == "no",]
no <- subset(no, select = -starved)
#Change column names of starved
names(yes) <- c("organism","width_st","length_st","shape_factor_st")
#Combine

lv <- no %>% left_join(yes, by = "organism") %>% 
 # select(-organism1) %>% 
  mutate(diff_width = width_st - width) %>% 
  mutate(perc_diff_width = diff_width/width*100) %>%
  mutate(diff_length = length_st - length) %>% 
  mutate(perc_diff_length = diff_length/length*100) 

sub <- lv

starved_width_axis_text <- expression(paste("Starved width (",mu,"m)", sep = ""))
non_starved_width_axis_text <- expression(paste("Non-starved width (",mu,"m)", sep = ""))

p1 <- ggplot(sub, aes(x = width, y = width_st)) + 
  geom_point(aes(fill = organism), shape = 21, alpha = 0.8, colour = "#656565", show.legend = FALSE) + 
  scale_x_continuous(limits = c(0,1)) + 
  scale_y_continuous(limits = c(0,1)) + 
  geom_abline(intercept = 0, slope = 1) +
  basic_layout + 
  labs(title = "Cell width", x = non_starved_width_axis_text, y = starved_width_axis_text)
p1

#Starved length by non-starved length
starved_length_axis_text <- expression(paste("Starved length (",mu,"m)", sep = ""))
non_starved_length_axis_text <- expression(paste("Non-starved length (",mu,"m)", sep = ""))

p2 <- ggplot(sub, aes(x = length, y = length_st), group = 1) + 
  geom_point(aes(fill = organism), shape = 21, alpha = 0.8, colour = "#656565", show.legend = TRUE) + 
  scale_x_continuous(limits = c(0,2)) + 
  scale_y_continuous(limits = c(0,2)) + 
  geom_abline(intercept = 0, slope = 1) +
  basic_layout + 
  theme(
    legend.justification = c(0,1),
    legend.text = element_text(size = 8),
    legend.key.height = unit(0.3, "cm"),
    legend.spacing.y = unit(0.2, "cm"),
    axis.title.y = element_blank()
  )+
  labs(title = "Cell length", x = non_starved_length_axis_text, y = starved_length_axis_text, fill = "Species")
p2

fig <- ggarrange(p1, p2, widths = c(1,1.57))
fig

ggsave(filename = "FigureB1.1.png", plot = fig, device = "png", path = save_path, units = "cm", width = 18, height = 8, dpi = 600, limitsize = TRUE)


############
# Figure 4 #
############

# Histograms of specific resource use groups

# 3x3 panel of traits: genome_size, d1_mid, growth_rate_residuals
# resource use groups: Sulfate and thiosulfate reducers, Hydrogen oxidizers, cclx degraders

# NOTE: For clarity, some axes are scaled to leave out taling data. 
# This does not affect the conclusions of the figures in the context of this paper.

df2 <- df %>% filter(pathway_final %in% c("Sulfate, Thiosulfate red.","hydrogen_oxidation_dark","cclx_degraders")) %>% 
  select(pathway_final, pathway_name, genome_size, d1_mid, temp_adjusted_maxgrowth) %>%
  mutate(pathway_final = as.factor(pathway_final))

#Order pathways as factors
df2$pathway_final <- factor(df2$pathway_final, levels = c("Sulfate, Thiosulfate red.","hydrogen_oxidation_dark","cclx_degraders"))


# genome_size
#########

trait <- "genome_size"
title <- "Genome size (Mbp)"

sub <- df2 %>% filter(!is.na(.data[[trait]]))
sub_non <- full_df %>% filter(!is.na(.data[[trait]]))

# Create data frame with total count within each group
dat_text <- sub %>% group_by(pathway_final) %>% summarise(n = n())

# Get names of pathways in the order of appearance for facet labels
pathway_names <- unique(sub$pathway_name)
names(pathway_names) <- unique(sub$pathway_final)

p1 <- sub %>% 
  ggplot(aes(x = .data[[trait]])) + 
  geom_histogram(data = sub_non, bins = 30, aes(x = .data[[trait]], y= ..ncount..), fill = "gray", alpha = 1) +
  geom_histogram(bins = 30, aes(x = .data[[trait]], y= ..ncount..), fill = colours[2], alpha = 1) +
  scale_y_continuous() +
  scale_x_log10(limits = c(0.5,18)) +
  basic_layout +
  theme(
    strip.background = element_blank(),
    strip.text.y = element_blank()
  ) +
  geom_text(data = dat_text, mapping = aes(label = sprintf("n = %s",n), x = 0.2, y = 0.9), size = 3.5) +
  facet_grid(pathway_final~., labeller = labeller(pathway_final = pathway_names)) + 
  labs(x = title, y = "Scaled proportion")
p1

print("Note: For clarity, X-axis scaled to min 0.5Mb (cuts off tail).")

# d1_mid
#########

trait <- "d1_mid"
title <- expression(paste("Mean diameter (", mu, "m)"))

sub <- df2 %>% filter(!is.na(.data[[trait]]))
sub_non <- full_df %>% filter(!is.na(.data[[trait]]))

# Create data frame with total count within each group
dat_text <- sub %>% group_by(pathway_final) %>% summarise(n = n())

# Get names of pathways in the order of appearance for facet labels
pathway_names <- unique(sub$pathway_name)
names(pathway_names) <- unique(sub$pathway_final)

p2 <- sub %>% 
  ggplot(aes(x = .data[[trait]])) + 
  geom_histogram(data = sub_non, bins = 30, aes(x = .data[[trait]], y= ..ncount..), fill = "gray", alpha = 1) +
  geom_histogram(bins = 30, aes(x = .data[[trait]], y= ..ncount..), fill = colours[2], alpha = 1) +
  scale_y_continuous() +
  scale_x_log10(limits = c(0.1,10)) +
  basic_layout +
  theme(
    strip.background = element_blank(),
    strip.text.y = element_blank()
  ) +
  geom_text(data = dat_text, mapping = aes(label = sprintf("n = %s",n), x = 0.2, y = 0.9), size = 3.5) +
  facet_grid(pathway_final~., labeller = labeller(pathway_final = pathway_names)) + 
  labs(x = title, y = "Scaled proportion")
p2

print("Note: For clarity, X-axis scaled to max 10um (cuts off tail).")

# temp_adjusted_maxgrowth
#########

trait <- "temp_adjusted_maxgrowth"
title <- "Temp-adj max growth (log10 resid)"

sub <- df2 %>% filter(!is.na(.data[[trait]]))
sub_non <- full_df %>% filter(!is.na(.data[[trait]]))

# Create data frame with total count within each group
dat_text <- sub %>% group_by(pathway_final) %>% summarise(n = n())

# Get names of pathways in the order of appearance for facet labels
pathway_names <- unique(sub$pathway_name)
names(pathway_names) <- unique(sub$pathway_final)

p3 <- sub %>% 
  ggplot(aes(x = .data[[trait]])) + 
  geom_histogram(data = sub_non, bins = 30, aes(x = .data[[trait]], y= ..ncount..), fill = "gray", alpha = 1) +
  geom_histogram(bins = 30, aes(x = .data[[trait]], y= ..ncount..), fill = colours[2], alpha = 1) +
  scale_y_continuous() +
  basic_layout +
  geom_text(data = dat_text, mapping = aes(label = sprintf("n = %s",n), x = -1.5, y = 0.9), size = 3.5) +
  facet_grid(pathway_final~., labeller = labeller(pathway_final = pathway_names)) + 
  labs(x = title, y = "Scaled proportion")
p3

#Save as combined figures
fig <- ggarrange(p1, p2, p3, ncol = 3, nrow = 1,  widths = c(1, 1, 1.1))
fig

ggsave(filename = "Fig4.png", plot = fig, device = "png", path = save_path, units = "cm", width = 18, height = 14, dpi = 600, limitsize = TRUE) 

