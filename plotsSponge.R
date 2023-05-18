#
library(vegan)

qualitative_colors <- 
  c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


#
spongeMicro <- sponge_clean %>% 
  group_by(Sample) %>% 
  mutate(RelativeAbundance = Abundance*100/sum(Abundance)) %>% 
  ungroup()

#
metadata <- read.table("metadata", sep = "\t", header = TRUE)
metadata <- metadata %>% 
  rename(Sample = Run...Assembly.accession, Biome = Sample.description) %>% 
  select(Sample, Biome)

# Merge microbiome data with metadata
spongeMicroFull <- spongeMicro %>% 
  full_join(metadata, by = "Sample") %>% 
  mutate(Biome = case_when(Biome == "sponge-derived microbial pellet" ~ "Sponge",
                           Biome == "Sediment microbiome" ~ "Sediment",
                           Biome == "Seawater microbiome" ~ "Seawater"))

##
spongeMicroFull <- spongeMicroFull %>% 
  group_by(Sample) %>% 
  mutate(Species_richness = specnumber(Abundance),
         Shannon_index = diversity(Abundance, index = "shannon"),
         Simpson_index = diversity(Abundance, index = "inv"))
  

## Plots ##

## basic level
spongeMicroFull %>% 
  ggplot(aes(x = Biome, y = Species_richness)) + 
  geom_point() + 
  theme_bw() + 
  theme(panel.grid.minor = element_blank()) + 
  labs(title = "Boring, simplistic species richness \U1F634")


## but what about other diversity metrics?
spongeMicroFull %>% 
  pivot_longer(cols = c("Species_richness", "Shannon_index", "Simpson_index"), values_to = "Diversity", names_to = "Index") %>% 
  ggplot(aes(x = Biome, y = Diversity)) + 
  geom_point() + 
  theme_bw() + 
  theme(panel.grid.minor = element_blank()) + 
  facet_grid(~Index) + 
  labs(title = "PROBLEMATIC AXIS! \U2620")

# ok, but boring
spongeMicroFull %>% 
  pivot_longer(cols = c("Species_richness", "Shannon_index", "Simpson_index"), values_to = "Diversity", names_to = "Index") %>% 
  ggplot(aes(x = Biome, y = Diversity)) + 
  geom_point() + 
  theme_bw() + 
  theme(panel.grid.minor = element_blank()) + 
  facet_wrap(~Index,scales = "free") + 
  labs(title = "An ok, basic plot \U1F60C")

###

spongeMicroFull %>% 
  ggplot(aes(Biome, RelativeAbundance, fill = Domain)) + 
  geom_col() + 
  theme_bw() + 
  labs(title = "Why did this happen? \U1F622", 
       y = "Relative Abundance (%)",
       subtitle = "\U26A0 Don't do this \U26A0") + 
  theme(legend.position = "top")


# By each sample
spongeMicroFull %>% 
  ggplot(aes(Sample, RelativeAbundance, fill = Domain)) + 
  geom_col() + 
  theme_bw() + 
  labs(title = "Why did this happen?", 
       y = "Relative Abundance (%)",
       subtitle = "Relative abundance was calculated by sample, not by biome") + 
  theme(legend.position = "top")


# By each sample (Domain)
spongeMicroFull %>% 
  ggplot(aes(Sample, RelativeAbundance, fill = Domain)) + 
  geom_col() + 
  theme_bw() + 
  labs(title = "Some basic options", 
       y = "Relative Abundance (%)",
       subtitle = "Divide by pannels\n Works, because we have two categorical/qualitative variables ") + 
  theme(legend.position = "top") + 
  facet_grid(~Biome, scales = "free")


# By each sample (Phylum -- it no longer works)
spongeMicroFull %>% 
  ggplot(aes(Sample, RelativeAbundance, fill = Phylum)) + 
  geom_col() + 
  theme_bw() + 
  labs(title = "What if we extend to Phylum level??", 
       y = "Relative Abundance (%)",
       subtitle = "\U1F616 Please stop doing this \U1F616") + 
  facet_grid(~Biome, scales = "free") 

### Common workarounds

# main problem is colloring, but number of samples is also a problem

# the coloring is misleading, because it gives the sense of connection, where there is none

## Try to give some biological meaning
## e.g., order them by evolutionary relationship
## took a phylogenetic tree
library(RColorBrewer)
spongeMicroFull %>% 
  mutate(Phylum = factor(Phylum, levels = c("Proteobacteria",
                                            "Nitrospinae",
                                            "Nitrospirae",
                                            "Acidobacteria", # group1
                                            "Fusobacteria",
                                            "Spirochaetes",
                                            "Lentisphaerae",
                                            "Verrucomicrobia",
                                            "Planctomycetes",
                                            "Chlamydiae",
                                            "Ignavibacteriae",
                                            "Bacteroidetes", # group 2
                                            "Synergistetes",
                                            "Cyanobacteria",
                                            "Deinococcus-Thermus",
                                            "Actinobacteria",
                                            "Tenericutes",
                                            "Firmicutes", 
                                            "Armatimonadetes",             
                                            "Balneolaeota", 
                                            "Calditrichaeota",  
                                            "Chloroflexi",                 
                                            "Elusimicrobia",
                                            "Fibrobacteres", 
                                            "Gemmatimonadetes",# group 3
                                            "Thaumarchaeota", 
                                            "Euryarchaeota",
                                            "Kiritimatiellaeota",
                                            "Rhodothermaeota", 
                                            "Candidatus_Falkowbacteria", # group 4   
                                            "Candidatus_Kaiserbacteria",
                                            "Candidatus_Kerfeldbacteria",
                                            "Candidatus_Latescibacteria",  
                                            "Candidatus_Levybacteria",
                                            "Candidatus_Magasanikbacteria",
                                            "Candidatus_Nomurabacteria",   
                                            "Candidatus_Omnitrophica",
                                            "Candidatus_Pacebacteria",
                                            "Candidatus_Peregrinibacteria",
                                            "Candidatus_Poribacteria", 
                                            "Candidatus_Raymondbacteria",
                                            "Candidatus_Tectomicrobia",    
                                            "Candidatus_Uhrbacteria", 
                                            "Candidatus_Yonathbacteria",
                                            "NA"
                                            ))) %>% # groups 4 and 5 
  ggplot(aes(Sample, RelativeAbundance, fill = Phylum)) + 
  geom_col() + 
  theme_bw() + 
  labs(title = "Organizes color and ordering by biological groups", 
       y = "Relative Abundance (%)",
       subtitle = "To much work for little reward \U1F62B, but is an improvement \U1F60C") + 
  facet_grid(~Biome, scales = "free") + 
  scale_fill_manual(values = c(brewer.pal(4, "YlOrBr"),
                               brewer.pal(8, "Greens"),
                               brewer.pal(9, "Purples"), brewer.pal(3, "Purples"),
                               brewer.pal(5, "Blues"),
                               brewer.pal(9, "Reds"),brewer.pal(6, "Reds"), "grey"))



## filter less important taxa -- this is tricky, because data is compositional
spongeMicroFull %>% 
  filter(RelativeAbundance > 1) %>% 
  ggplot(aes(Sample, RelativeAbundance, fill = Phylum)) + 
  geom_col() + 
  theme_bw() + 
  labs(title = "Remove taxonomic lineages below 1% relative abundance (per sample)", 
       y = "Relative Abundance (%)",
       subtitle = "Also an improvement, still not aceptable \U2620, but why not? \U1F611") +
  facet_grid(~Biome, scales = "free") + 
  scale_fill_manual(values = brewer.pal(n = 11, "Set3"))


## Workaround is to explicitly indicate whats missing (but still misleading)
spongeMicroFull %>% 
  mutate(PhylumModified = ifelse(RelativeAbundance > 1, Phylum, "Taxa below 1%")) %>% 
  ggplot(aes(Sample, RelativeAbundance, fill = PhylumModified)) + 
  geom_col() + 
  theme_bw() + 
  labs(title = "Remove taxonomic lineages below 1% relative abundance (per sample)", 
       y = "Relative Abundance (%)",
       subtitle = "Getting better! But not good enough, why? \U1F624") +
  facet_grid(~Biome, scales = "free") + 
  scale_fill_manual(values = brewer.pal(n = 11, "Set3"))



### We love heatmaps, because they are cute
### but more often than not, they are bad visualization tools

# My chalenge is for you to realize that you can change the place of things
## this is wrong
spongeMicroFull %>%
  mutate(PhylumModified = ifelse(RelativeAbundance > 1, Phylum, "Taxa below 1%")) %>%
  ggplot(aes(x = Sample, y = RelativeAbundance, fill = PhylumModified)) + 
  geom_col(position = "dodge") + 
  facet_grid(~Biome, scales = "free") + 
  scale_fill_manual(values = brewer.pal(n = 11, "Set3")) +
  theme_bw() + 
  labs(title = "On the right track, but wrong \U1F635 (why??)", 
       y = "Relative Abundance (%)") + 
  theme(legend.position = "top")


## this was corrected
spongeMicroFull %>%
  mutate(PhylumModified = ifelse(RelativeAbundance > 1, Phylum, "Taxa below 1%")) %>%
  group_by(Biome, PhylumModified, Sample) %>% 
  summarise(RelativeAbundance = sum(RelativeAbundance)) %>% 
  ggplot(aes(x = Sample, y = RelativeAbundance, fill = PhylumModified)) + 
  geom_col(position = "dodge") + 
  facet_grid(~Biome, scales = "free") + 
  scale_fill_manual(values = brewer.pal(n = 11, "Set3")) +
  theme_bw() + 
  labs(title = "Maybe ugly, but best so far \U1F60A", 
       y = "Relative Abundance (%)") + 
  theme(legend.position = "top")



### You can try and forget relative abundance,
### but then be carefull and normalize data (so redundant with relative abundance)
library(purrr)
spongeMicroFullRarefied <- spongeMicroFull %>% 
  group_by(Sample) %>% 
  nest() %>% 
  mutate(RarefiedAbundance = map(.x = data, ~as.data.frame(t(rrarefy(.x$Abundance,sample = 1000))))) %>% 
  unnest(c(data,RarefiedAbundance)) %>% 
  rename(RarefiedAbundance = "V1")


spongeMicroFullRarefied %>% 
  group_by(Sample) %>% 
  mutate(RarefiedRelativeAbundance = RarefiedAbundance*100/sum(RarefiedAbundance)) %>% 
  mutate(PhylumModified = ifelse(RarefiedRelativeAbundance > 1, Phylum, "Taxa below 1%")) %>%
  ggplot(aes(x = Sample, y = RarefiedAbundance, fill = PhylumModified)) + 
  #geom_col(position = "dodge") +
  geom_col()+
  facet_grid(~Biome, scales = "free") + 
  scale_fill_manual(values = c(brewer.pal(n = 12, "Set3"), "grey")) +
  theme_bw() + 
  labs(title = "Absolute abundance (1000 reads per sample)", 
       y = "Relative Abundance (%)",
       subtitle = "We are walking in circles, the main problem remains, what is it? \U1F62C") + 
  theme(legend.position = "top")


## the main problem is that we are still using color to distinguish phylum
spongeMicroFullRarefied %>% 
  group_by(Sample) %>% 
  mutate(RarefiedRelativeAbundance = RarefiedAbundance*100/sum(RarefiedAbundance)) %>% 
  mutate(PhylumModified = ifelse(RarefiedRelativeAbundance > 1, Phylum, "Taxa below 1%")) %>%
  ungroup() %>% 
  group_by(Biome, PhylumModified, Sample) %>% 
  summarise(RarefiedAbundance = sum(RarefiedAbundance)) %>% 
  ggplot(aes(x = Sample, y = RarefiedAbundance, fill = PhylumModified)) + 
  geom_col(position = "dodge") +
  facet_grid(~Biome, scales = "free") + 
  scale_fill_manual(values = c(brewer.pal(n = 12, "Set3"), "grey")) +
  theme_bw() + 
  labs(title = "Absolute abundance (1000 reads per sample)", 
       y = "Relative Abundance (%)",
       subtitle = "We are walking in circles, the main problem remains, what is it? \U1F62C") + 
  theme(legend.position = "top")


## boxplot approach

spongeMicroFull %>% 
  ggplot(aes(x = reorder(Phylum, RelativeAbundance), 
             y = RelativeAbundance, fill = Biome)) + 
  geom_col(position = "dodge") + 
  coord_flip() + 
  theme_bw() + 
  labs(x = "Phylum",
       y = "Relative Abundance (%)",
       title = "Not perfect, but best so far \U1F60A",
       subtitle = "All taxa. One small problem remains \U1F613") + 
  theme(legend.position = "top",
        panel.grid.minor = element_blank()) + 
  scale_fill_manual(values = qualitative_colors[c(2,4,6)])

## we can still improve it by removing rare taxa
spongeMicroFull %>% 
  mutate(PhylumModified = ifelse(RelativeAbundance > 1, Phylum, "Taxa below 1%")) %>%
  ggplot(aes(x = reorder(PhylumModified, RelativeAbundance), 
             y = RelativeAbundance, fill = Biome)) + 
  geom_col(position = "dodge") + 
  coord_flip() + 
  theme_bw() + 
  labs(x = "Phylum",
       y = "Relative Abundance (%)",
       title = "Not perfect, but best so far \U1F60A",
       subtitle = "We can remove low abundance taxa again. \nOne small problem remains \U1F613") + 
  theme(legend.position = "top",
        panel.grid.minor = element_blank()) + 
  scale_fill_manual(values = qualitative_colors[c(2,4,6)])


## The sample information was lost.


#### What if we try to show everything?

## group taxa within samples
spongeMicroFull %>% 
  group_by(Phylum, Biome, Sample) %>% 
  summarise(PhylumAbundance = sum(Abundance)) %>%
  ungroup() %>% 
  group_by(Sample) %>% 
  mutate(RelativeAbundance = PhylumAbundance*100/sum(PhylumAbundance)) %>% 
  ggplot(aes(x = reorder(Phylum, PhylumAbundance), 
             y = PhylumAbundance)) + 
  stat_summary(size = 0.5) + 
  facet_grid(~Biome) + 
  coord_flip() + 
  theme_bw() + 
  labs(x = "Phylum",
       y = "Relative Abundance (%)",
       title = "Centrality metric --- mean with SD",
       subtitle = "Too few samples for mean to be informative \ndifficult to see missing phyla across biomes") + 
  theme(legend.position = "top",
        panel.grid.minor = element_blank()) + 
  scale_fill_manual(values = qualitative_colors[c(2,4,6)])

###
## group taxa within samples
spongeMicroFull %>% 
  group_by(Phylum, Biome, Sample) %>% 
  summarise(PhylumAbundance = sum(Abundance)) %>%
  ungroup() %>% 
  group_by(Sample) %>% 
  mutate(RelativeAbundance = PhylumAbundance*100/sum(PhylumAbundance)) %>%
  ggplot(aes(x = reorder(Phylum, PhylumAbundance), 
             y = RelativeAbundance)) + 
  geom_point(size = 0.5) + 
  facet_grid(~Biome) + 
  coord_flip() + 
  theme_bw() + 
  labs(x = "Phylum",
       y = "Relative Abundance (%)",
       title = "Each sample is a point",
       subtitle = "Small values are hard to interpret, are they missing or are they very few? ") + 
  theme(legend.position = "top",
        panel.grid.minor = element_blank()) + 
  scale_fill_manual(values = qualitative_colors[c(2,4,6)])

### we can try log10
spongeMicroFull %>% 
  group_by(Phylum, Biome, Sample) %>% 
  summarise(PhylumAbundance = sum(Abundance)) %>%
  ungroup() %>% 
  group_by(Sample) %>% 
  mutate(RelativeAbundance = PhylumAbundance*100/sum(PhylumAbundance)) %>%
  ggplot(aes(x = reorder(Phylum, PhylumAbundance), 
             y = RelativeAbundance)) + 
  geom_point(size = 0.5) + 
  facet_grid(~Biome) + 
  coord_flip() + 
  theme_bw() + 
  labs(x = "Phylum",
       y = "Relative Abundance (%)",
       title = "Each sample is a point",
       subtitle = "Small values are hard to interpret, are they missing or are they very few? ") + 
  theme(legend.position = "top",
        panel.grid.minor = element_blank()) + 
  scale_fill_manual(values = qualitative_colors[c(2,4,6)]) + 
  scale_y_log10()

## even better, remove zeros
spongeMicroFull %>% 
  group_by(Phylum, Biome, Sample) %>% 
  summarise(PhylumAbundance = sum(Abundance)) %>%
  ungroup() %>% 
  group_by(Sample) %>% 
  mutate(RelativeAbundance = PhylumAbundance*100/sum(PhylumAbundance)) %>%
  filter(RelativeAbundance > 0) %>% 
  ggplot(aes(x = reorder(Phylum, PhylumAbundance), 
             y = RelativeAbundance)) + 
  geom_point(size = 0.5) + 
  facet_grid(~Biome) + 
  coord_flip() + 
  theme_bw() + 
  labs(x = "Phylum",
       y = "Relative Abundance (%)",
       title = "Each sample is a point",
       subtitle = "Log10 and zeros removed") + 
  theme(legend.position = "top",
        panel.grid.minor = element_blank()) + 
  scale_fill_manual(values = qualitative_colors[c(2,4,6)]) + 
  scale_y_log10()


### again, we can remove low abundance taxa
spongeMicroFull %>% 
  group_by(Phylum, Biome, Sample) %>% 
  summarise(PhylumAbundance = sum(Abundance)) %>%
  ungroup() %>% 
  group_by(Sample) %>% 
  mutate(RelativeAbundance = PhylumAbundance*100/sum(PhylumAbundance)) %>%
  filter(RelativeAbundance > 1) %>% 
  ggplot(aes(x = reorder(Phylum, PhylumAbundance), 
             y = RelativeAbundance)) + 
  geom_point() + 
  facet_grid(~Biome) + 
  coord_flip() + 
  theme_bw() + 
  labs(x = "Phylum",
       y = "Relative Abundance (%)",
       title = "Each sample is a point",
       subtitle = "Log10 and zeros removed and relative abundance > 1%") + 
  theme(legend.position = "top",
        panel.grid.minor = element_blank())+
  scale_fill_manual(values = qualitative_colors[c(2,4,6)]) + 
  scale_y_log10()

### We can think on yet other approaches, be creative and don't be affraid to test
gridExtra::grid.arrange(
spongeMicroFull %>% 
  mutate(PhylumModified = ifelse(RelativeAbundance > 1, Phylum, "Taxa below 1%")) %>%
  filter(RelativeAbundance > 0) %>% 
  ggplot(aes(x = reorder(PhylumModified, Abundance), 
             y = Abundance)) + 
  geom_boxplot(aes(fill = Biome)) +
  coord_flip() + 
  theme_bw() + 
  labs(x = "Phylum",
       y = "Abundance",
       title = "Boxplots (median, etc) --- all taxonomic lineages",
       subtitle = "") + 
  theme(legend.position = "top",
        panel.grid.minor = element_blank())+
  scale_fill_manual(values = qualitative_colors[c(2,4,6)]) + 
  scale_y_log10()
,
spongeMicroFull %>% 
  mutate(PhylumModified = ifelse(RelativeAbundance > 1, Phylum, "Taxa below 1%")) %>%
  filter(RelativeAbundance > 0) %>% 
  ggplot(aes(x = reorder(PhylumModified, Abundance), 
             y = Abundance)) + 
  stat_summary(aes(col = Biome),
               position = position_jitter(width = 0.2)) +
  coord_flip() + 
  theme_bw() + 
  labs(x = "Phylum",
       y = "Abundance",
       title = "Mean \U2213 SD --- all taxonomic lineages",
       subtitle = "Zeros removed and relative abundance > 1%") + 
  theme(legend.position = "top",
        panel.grid.minor = element_blank())+
  scale_fill_manual(values = qualitative_colors[c(2,4,6)]) + 
  scale_color_manual(values = qualitative_colors[c(2,4,6)]) + 
  scale_y_log10(), ncol = 2
)

### we can go even further ####
# don't forget that facets are also aesthetics!
spongeMicroFull %>% 
  mutate(PhylumModified = ifelse(RelativeAbundance > 1, Phylum, "Taxa below 1%")) %>%
  filter(RelativeAbundance > 0) %>% 
  ggplot(aes(y = Species_richness, x = Biome)) + 
  geom_point() + 
  facet_wrap(~PhylumModified, scales = "free_y") + 
  theme_bw() + 
  theme(panel.grid.minor = element_blank()) + 
  labs(title = "Species richness per Phylum",
       y = "Species richness",
       x = "Biome",
       subtitle = "Most economic approach so far -- carefull with axis\nNote that we still have: \n-Shape;\n-color;\n-and size.")






