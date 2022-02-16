##################################################################
## Calculating karyotype proportion across labrador
##
## Matt Brachmann (PhDMattyB)
##
## 2021-06-02
##
##################################################################

## Load the data manipulation work horse
library(tidyverse)
library(maps)
library(scatterpie)


# Functions ---------------------------------------------------------------

`%!in%` = Negate(`%in%`)
# Charr project -----------------------------------------------------------

setwd('~/Charr_Adaptive_Introgression/Charr_Project_1/GeneticData/Fst_sliding_window')


inversion = read_csv('Lab_AC08_Inversion_grouping.csv')
south_ancestral = read_csv('Lab_AC08_South_Ancestral_Karyotype.csv')
North_ancestral = read_csv('Lab_AC08_North_Ancestral_Karyotype.csv')
hetero = read_csv('Lab_AC08_heterozygous_grouping.csv')

ancestral = bind_rows(south_ancestral, 
                      North_ancestral)


# Begin data processing ---------------------------------------------------

## Need to make labels so we know what everything is when we combine into
## large data frame (big pappi)

inver_lab = rep('rearranged homozygous', 
                length(inversion$Population)) %>% 
  as_tibble()
inversion = bind_cols(inversion, 
                      inver_lab)

View(inversion)
anc_lab = rep('non-rearranged homozygous', 
              length(ancestral$Population)) %>% 
  as_tibble()
ancestral = bind_cols(ancestral, 
                      anc_lab)

het_lab = rep('rearranged heterozygous', 
              length(hetero$Population)) %>% 
  as_tibble()
hetero = bind_cols(hetero, 
                   het_lab)


big_pappi = bind_rows(inversion, 
                      ancestral, 
                      hetero)

View(big_pappi)
small_pappi = big_pappi %>% 
  dplyr::select(Population, 
         Latitude, 
         Longitude) %>% 
  distinct(Population, 
           .keep_all = T)

big_pappi_freq = big_pappi %>% 
  group_by(Population, 
           value) %>% 
  summarise(n = n())%>%
  mutate(freq = n / sum(n)) 


big_pappi %>% 
  # filter(Population %in% c('FRD', 
  #                         'FRN')) %>% 
  # filter(value == 'rearranged heterozygous') %>%
  dplyr::select(Population, 
                IndividualID, 
                Name, 
                Latitude, 
                Longitude, 
                Loc2, 
                value) %>% 
  filter(value %in% c('rearranged homozygous',
                      'non-rearranged homozygous',
                      'rearranged heterozygous')) %>%
  # filter(value == 'rearranged homozygous') %>%
  filter(Population == 'NOR') %>% 
  # filter(Population %!in% c('BLD', 
  #                           'BRG', 
  #                           'ENG', 
  #                           'GDL')) %>%
  arrange(Population, 
          value) %>% 
  group_by(Population, 
           Name,
           value) %>% 
  # summarise(num = n()) %>% 
  write_csv('NOR_AC08_SV.csv')

clean_data = inner_join(big_pappi_freq, 
           small_pappi, 
           by = 'Population')

clean_data %>% 
  filter(value == 'rearranged homozygous') %>% 
  mutate(percent = freq*100) %>% 
  View()

spread_data = clean_data %>% 
  group_by(value) %>% 
  mutate(i1 = row_number()) %>% 
  spread(value, 
         freq) %>% 
  replace(is.na(.), 0) %>% 
  group_by(Population) %>% 
  mutate(Lat_jitter = Latitude + (Latitude[1] - jitter(Latitude[1], 
                                           amount = 1,
                                           factor = 1)), 
       Long_jitter = Longitude + (Longitude[1] - jitter(Longitude[1], 
                                              amount = 1,
                                              factor = 0.9)))

View(spread_data)
# Map data ----------------------------------------------------------------

East_coastish = map_data('world') %>% 
  filter(region == 'Canada') %>% 
  as_tibble() %>% 
  rename(Latitude = lat, 
         Longitude = long) %>% 
  filter(Longitude < -47.0, 
         Longitude > -65.0, 
         Latitude < 60.0, 
         Latitude > 47.5) %>% 
  rename(lat = Latitude, 
         long = Longitude)


# Map plot ----------------------------------------------------------------

# map_palette = c('#0583F2', 
#                   '#F28705',
#                   '#F20530')

new_palette = c('#03588C',
                '#F25C5C',
                '#04BFAD')

theme_set(theme_bw())
# theme_set(theme_void())

karyotype_map = ggplot(East_coastish) +
    geom_map(data = East_coastish, 
             map = East_coastish, 
             aes(x = long, 
                 y = lat, 
                 map_id = region), 
             col = 'white', 
             fill = 'black')+
  labs(x = 'Longitude', 
       y = 'Latitude', 
       color = 'Karyotype', 
       fill = 'Karyotype')+
  # scale_fill_manual(values = map_palette)+
  scale_fill_manual(values = new_palette)+
  theme(axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12), 
        legend.title = element_text(size = 14), 
        legend.text = element_text(size = 12))+
    geom_scatterpie(data = spread_data, 
                    aes(x = Long_jitter, 
                        y = Lat_jitter, 
                        group = Population), 
                    pie_scale = 2, 
                    cols = colnames(spread_data[,c(6:8)]))

karyotype_map
  
ggsave('~/Charr_Adaptive_Introgression/Charr_Project_1/Figures/Karyotypes_map_Newcols.tiff', 
       plot = karyotype_map, 
       dpi = 'retina', 
       units = 'cm', 
       width = 20, 
       height = 15)


##
# All populations on AC08 -------------------------------------------------
ancestral_all = read_csv('All_Pops_AC08_Ancestral_karyotype.csv')
Hetero_all = read_csv('All_Pops_AC08_Heterozygous_karyotype.csv')
Inversion_all = read_csv('All_Pops_AC08_Derived_Karyotype.csv')



# All pops AC08 data cleaning ---------------------------------------------


inver_lab = rep('rearranged homozygous', 
                length(Inversion_all$Population)) %>% 
  as_tibble()
inversion_all = bind_cols(Inversion_all, 
                      inver_lab)

anc_lab = rep('non-rearranged homozygous', 
              length(ancestral_all$Population)) %>% 
  as_tibble()
ancestral_all = bind_cols(ancestral_all, 
                      anc_lab)

het_lab = rep('rearranged heterozygous', 
              length(Hetero_all$Population)) %>% 
  as_tibble()
hetero_all = bind_cols(Hetero_all, 
                   het_lab)


bigger_pappi = bind_rows(inversion_all, 
                      ancestral_all, 
                      hetero_all)

smallish_pappi = bigger_pappi %>% 
  dplyr::select(Population, 
                Latitude, 
                Longitude) %>% 
  distinct(Population, 
           .keep_all = T)

bigger_pappi_freq = bigger_pappi %>% 
  group_by(Population, 
           value) %>% 
  summarise(n = n())%>%
  mutate(freq = n / sum(n)) 

clean_data = inner_join(bigger_pappi_freq, 
                        smallish_pappi, 
                        by = 'Population')

spread_data = clean_data %>% 
  group_by(value) %>% 
  mutate(i1 = row_number()) %>% 
  spread(value, 
         freq) %>% 
  replace(is.na(.), 0)


# All pops AC08 map data --------------------------------------------------

spread_data %>% 
  arrange(-Longitude)
  
maine = map_data('state') %>%
  filter(region == 'maine') %>% 
  as_tibble()

other_pops = map_data('world') %>% 
  filter(region %in% c('Canada', 
                       'Greenland', 
                       'Norway', 
                       'Iceland', 
                       'UK', 
                       'state'))%>% 
  as_tibble()

all_map_data = bind_rows(other_pops, 
          maine)

please_work =  all_map_data %>% 
  rename(Latitude = lat, 
         Longitude = long) %>% 
  filter(Longitude < 20.00, 
         Longitude > -80.00, 
         Latitude < 70.00, 
         Latitude > 43.0) %>% 
  rename(lat = Latitude, 
         long = Longitude)


# All pops AC08 map -------------------------------------------------------
map_palette = c('#D93E30', 
                '#D91895',
                '#0A21A6')

karyotype_map_all = ggplot(please_work) +
  geom_map(data = please_work, 
           map = please_work, 
           aes(x = long, 
               y = lat, 
               map_id = region), 
           col = 'white', 
           fill = 'black')+
  labs(x = 'Longitude', 
       y = 'Latitude', 
       color = 'Karyotype', 
       fill = 'Karyotype')+
  scale_fill_manual(values = map_palette)+
  theme(axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12), 
        legend.title = element_text(size = 14), 
        legend.text = element_text(size = 12))+
  geom_scatterpie(data = spread_data, 
                  aes(x = Longitude, 
                      y = Latitude, 
                      group = Population), 
                  pie_scale = 0.75, 
                  cols = colnames(spread_data[,c(6:8)]))+
  coord_fixed()

karyotype_map_all

ggsave('~/Charr_Adaptive_Introgression/Charr_Project_1/Figures/Karyotypes_map_All_populations.tiff', 
       plot = karyotype_map_all, 
       dpi = 'retina', 
       units = 'cm')

