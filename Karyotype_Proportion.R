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

clean_data = inner_join(big_pappi_freq, 
           small_pappi, 
           by = 'Population')

spread_data = clean_data %>% 
  group_by(value) %>% 
  mutate(i1 = row_number()) %>% 
  spread(value, 
         freq) %>% 
  replace(is.na(.), 0)


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

map_palette = c('#D93E30', 
                  '#D91895',
                  '#0A21A6')

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
  scale_fill_manual(values = map_palette)+
  theme(axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12), 
        legend.title = element_text(size = 14), 
        legend.text = element_text(size = 12))+
    geom_scatterpie(data = spread_data, 
                    aes(x = Longitude, 
                        y = Latitude, 
                        group = Population), 
                    pie_scale = 3, 
                    cols = colnames(spread_data[,c(6:8)]))
  
ggsave('~/Charr_Adaptive_Introgression/Charr_Project_1/Figures/Karyotypes_map.tiff', 
       plot = karyotype_map, 
       dpi = 'retina', 
       units = 'cm')
