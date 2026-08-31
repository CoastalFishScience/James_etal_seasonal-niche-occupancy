#' """ Data analysis, tables, and figures for 
#'     "Seasonality in species trophic niche width and 
#'     position alters community niche occupancy"
#'     @authors: W. Ryan James, Mack White, Herbert Leavitt, 
#'     Jonathan R. Rodemann, Gina A. Badlowski, James W. Sturges, 
#'     Justin S. Lesser, Valentina Bautista, Nicholas A. Castillo, 
#'     Sophia V. Costa, Andy Distrubell, Cody W. Eggenberger, 
#'     Lauren J. Kabat, Joshua O. Linenfelser, Nicolas Rivas, 
#'     Shakira Trabelsi, Jennifer S. Rehage, and Rolando O. Santos 
#'     Date : 8/28/16

#load necesary packages
librarian::shelf(tidyverse, hypervolume, vegan, ggpubr, 
                 grid, gtable, factoextra,viridis, ggforce)

# Calculate species level trophic niches ----
# Function to make random points of n length
# data from random sample with mean and sd but 
# can be generated between a high and low value if end_points = T
# *** Note chose either column names or column numbers for ID_rows and names
# either work but must be the same
# df = dataframe or tibble with each row containing 
#        unique entry for making random points 
# ID_rows = vector of column names or numbers with id information
# names_c = column name or number of name of measure variables
# mean_c = column name or column number of df with mean 
# sd_c = column name or column number of df with sd 
# n = number of points to randomly generate
# z_score = T or F, if T z-scores values
# end_points = T or F for if random points need to be generated between
#        a high and low end point (e.g. 5% and 95% interval)
#        low and high required if end_points = T
# low_c = column name or column number of df with lower bound to sample in
# high_c = column name or column number of df with upper bound to sample in
HVvalues = function(df, ID_rows, names_c, mean_c, sd_c, n, z_score = F,
                    end_points = F, low_c = NULL, high_c = NULL){
  require(tidyverse)
  require(truncnorm)
  
  # check to see if information is needed to restrict where points are
  if (end_points){
    if (is_empty(df[,low_c]) | is_empty(df[,high_c])){
      return(cat('Warning: low and/or high columns not specified \n
                  Specific and run again or end_points = F \n'))
    }
  }
  
  # check to see if there are more 
  if(T %in% duplicated(df[,c(ID_rows,names_c)])){
    return(cat('Warning: some of the rows contain duplicated information \n
                make sure data is correct \n'))
  }
  
  # rename variables to make code work
  if (is.numeric(mean_c)){
    names(df)[mean_c] = 'Mean'
  }else {
    df = df |>  rename(Mean = all_of(mean_c))
  }
  
  if (is.numeric(sd_c)){
    names(df)[sd_c] = 'SD'
  }else {
    df = df |>  rename(SD = all_of(sd_c))
  }
  
  if (end_points){
    if (is.numeric(low_c)){
      names(df)[low_c] = 'lower'
    }else {
      df = df |>  rename(lower = all_of(low_c))
    }
    
    if (is.numeric(high_c)){
      names(df)[high_c] = 'upper'
    }else {
      df = df |>  rename(upper = all_of(high_c))
    }
  }
  
  # make sure the names is not numeric 
  if (is.numeric(names_c)){
    names_c = names(df)[names_c]
  }
  
  labs = unique(as.vector(df[,names_c])[[1]])
  # generate random points within bounds
  if(end_points){
    
    df_tot = df |> slice(rep(1:n(), each=n))|> 
      mutate(point = 
               truncnorm::rtruncnorm(n(), a = lower, b = upper,
                                     mean = Mean, sd = SD)) |> 
      ungroup() |> 
      mutate(num = rep(1:n, times=nrow(df))) |>
      dplyr::select(-Mean, -SD, -lower, -upper)|>
      pivot_wider(names_from = all_of(names_c), values_from = point)|> 
      dplyr::select(-num)
  }else {
    # generate random points outside of bounds
    df_tot = df |> slice(rep(1:n(), each=n))|>
      mutate(point = 
               truncnorm::rtruncnorm(n(), mean = Mean, sd = SD)) |> 
      ungroup() |> 
      mutate(num = rep(1:n, times=nrow(df))) |>
      dplyr::select(-Mean, -SD)|>
      pivot_wider(names_from = all_of(names_c), values_from = point)|> 
      dplyr::select(-num)
  }
  if (z_score){
    df_tot = df_tot  |>  
      mutate(across(all_of(labs), scale))
  }
  
  return(df_tot)
  
}

# load mixing model data
d = read_csv('CESImixResults.csv')

# number or iterations
reps = 100

# generate points and z-score across iterations
set.seed(14)
df = d |> 
  slice(rep(1:n(), each=reps))|> 
  mutate(i = rep(1:reps, times=nrow(d))) |> 
  group_by(i) |> 
  nest() |> 
  mutate(points = map(data, \(data) HVvalues(df = data, ID_rows = c('species', 'season'),
                                             names_c = c('source'), 
                                             mean_c = 'mean', sd_c = 'sd', n = 30,
                                             end_points = T, 
                                             low_c = 'lowend', high_c = 'highend',
                                             z_score = T))) |> 
  select(i, points) |> 
  unnest(points)


# generate trophic niche hypervolumes
df = df |> 
  group_by(species, season, i) |> 
  nest() |> 
  mutate(hv = map(data, \(data) hypervolume_gaussian(data, name = paste(species, season, i,sep = '_'),
                                                     samples.per.point = 1000,
                                                     kde.bandwidth = estimate_bandwidth(data), 
                                                     sd.count = 3, 
                                                     quantile.requested = 0.95, 
                                                     quantile.requested.type = "probability", 
                                                     chunk.size = 1000, 
                                                     verbose = F)),
         hv_size = map_dbl(hv, \(hv) get_volume(hv)),
         centroid = map(hv, \(hv) get_centroid(hv)))

# save output so don't need to run
#write_rds(df, 'data/CESIhvAll.rds', compress = 'gz')

# Trophic niche comparisons ----
# read in saved hvs 
# df = read_rds('data/CESIhvAll.rds')

## Intraspecies across seasons ----
ov_sn = df |> 
  select(species, season, hv, hv_size) |> 
  pivot_wider(names_from = season, values_from = c(hv,hv_size)) |> 
  mutate(size_rat = hv_size_Dry/hv_size_Wet,
         set = map2(hv_Wet,hv_Dry, \(hv1, hv2) hypervolume_set(hv1, hv2, check.memory = F, verbose = F)),
         ov = map(set, \(set) hypervolume_overlap_statistics(set)),
         dist_cent = map2_dbl(hv_Wet,hv_Dry, \(hv1,hv2) hypervolume_distance(hv1, hv2, type = 'centroid', check.memory=F))) |> 
  unnest_wider(ov) |> 
  select(species, i, hv_size_Wet, hv_size_Dry, 
         size_rat, jaccard, sorensen,
         uniq_Wet = frac_unique_1, uniq_Dry = frac_unique_2, 
         dist_cent)

#write_csv(ov_sn, 'data/hvOv_season.csv')

### Seasonal niche strategies ----

# calculate brown and green seasonal differences in resource use
df_dif = df |> 
  select(species, season, i, centroid) |> 
  unnest_longer(centroid) |> 
  pivot_wider(names_from = season, values_from = centroid) |> 
  ungroup() |> 
  rename(basal_resource = centroid_id) |> 
  mutate(seasonal_centroid_diff = Dry - Wet) |> 
  group_by(species, i) |> 
  mutate(mean_seasonal_centroid_diff = mean(seasonal_centroid_diff)) |>
  ungroup() |> 
  ### pivot wider to make that seasonal change for each basal resource a new column
  pivot_wider(names_from = basal_resource, values_from = seasonal_centroid_diff, 
              id_cols = c(species, i, mean_seasonal_centroid_diff),
              ### add "_centroid_diff" to each of the new columns
              names_glue = "{basal_resource}_centroid_diff") |> 
  group_by(species, i) |> 
  mutate(brown = sum(Mangrove_centroid_diff, Seagrass_centroid_diff),
         green = sum(Epiphytes_centroid_diff, Algae_centroid_diff))

# join with intraspecies seasonal comparison metrics
dt_all = left_join(ov_sn, df_dif, by = c("species", "i")) |> 
  janitor::clean_names() |> 
  mutate(species = as.factor(species),
         i = as.factor(i)) |> 
  ### removing for now - think it may be better to focus on mean change if trying to generalize findings
  dplyr::select(species, i, brown, green,
                size_rat, dist_cent, uniq_wet, uniq_dry, sorensen)
glimpse(dt_all)

#### pca ----
# create tibble of numeric only values
dt_numeric = dt_all |> 
  select(brown:sorensen)

# perform PCA using prcomp()  
pca_result = prcomp(dt_numeric, center = TRUE, scale. = TRUE)
summary(pca_result)

# pull pca principal component calculations 
pca_scores = as.data.frame(pca_result$x) 


#### k-means clustering ----
### calculate the total within-cluster sum of squares for differnt values of "k"
set.seed(123)
wss = sapply(1:10, function(k) {
  sum(kmeans(pca_scores, centers = k, nstart = 20)$withinss)
})

# examine the wss plot to determine appropriate number of clusters
plot(1:10, wss, type = "b", xlab = "Number of Clusters", ylab = "Within-Cluster Sum of Squares")

# calculate k-means of 3 clusters based on wss
set.seed(123)
kmeans_result = kmeans(pca_scores, centers = 3, nstart = 25)

# add cluster information for plotting
pca_scores$cluster = as.factor(kmeans_result$cluster)
pca_scores$species = dt_all$species  # Adding species information for coloring


## Occupancy ----
# generate community trophic niche for each season
df_occ = df |> 
  select(species, season, i, hv) |> 
  pivot_wider(names_from = c(season, species), values_from = hv) |> 
  mutate(hv_join = pmap(list(`Wet_Bay anchovy`,Wet_Mojarra,
                             Wet_Pigfish, Wet_Pinfish, `Wet_Pink shrimp`,
                             `Wet_Rainwater killifish`,`Wet_Silver perch`,
                             `Dry_Bay anchovy`,Dry_Mojarra,
                             Dry_Pigfish, Dry_Pinfish, `Dry_Pink shrimp`,
                             `Dry_Rainwater killifish`,`Dry_Silver perch`), 
                        hypervolume_join),
         hv_occ = map(hv_join, \(data) hypervolume_n_occupancy(data, 
                                                               method = "box", 
                                                               classification = rep(c('Wet','Dry'), each = 7), 
                                                               box_density = 5000, 
                                                               FUN = mean,
                                                               verbose = F)),
         size = map(hv_occ, \(hv) get_volume(hv)),
         cent = map(hv_occ, \(hv) get_centroid_weighted(hv)),
         hv_df = map(hv_occ, \(data) hypervolume_to_data_frame(data))) 

#write_rds(df_occ, 'data/hvOCC.rds')
#df_occ = read_rds('data/hvOCC.rds')

# calculate the proportion of total niche occupied by single vs multiple spp
oc_p = df_occ  |> 
  select(i,hv_df) |> 
  unnest(hv_df) |> 
  rename(season = Name, occupancy = ValueAtRandomPoints) |> 
  mutate(cat = if_else(occupancy > 1/7, 'Multi species', 'Single species')) |> 
  group_by(i, season, cat) |> 
  summarize(n = n(), .groups = 'drop') |> 
  group_by(i, season) |> 
  mutate(prop = n/sum(n))

# calculate size of multi- and single- spp community niche breadth
oc_s = df_occ |> 
  select(i, size) |> 
  unnest_longer(size, indices_to = 'season') |> 
  left_join(oc_p) |> 
  mutate(size_occ = prop*size)

#write_csv(oc_s, 'data/occ_size.csv')

## Interspecies within season ----
spc = tibble(sp1 = unique(df$species),
             sp2 = unique(df$species))

df_sp = spc |>  expand(sp1,sp2)

df_sp = df_sp[!duplicated(t(apply(df_sp,1,sort))),] %>%
  filter(!(sp1 == sp2))

df_spw = df_sp %>%
  mutate(season = 'Wet')

df_spd = df_sp %>%
  mutate(season = 'Dry')

df_spp = bind_rows(df_spw, df_spd)

# make hv to bind
df1 = df |> 
  select(sp1 = species, season, i, hv1 = hv)

df2 = df |> 
  select(sp2 = species, season, i, hv2 = hv)

# create large df to store all data
ov_sp = df_spp |> 
  slice(rep(1:n(), each=reps))|> 
  mutate(i = rep(1:reps, times=nrow(df_spp))) |> 
  left_join(df1, by = c('sp1', 'season', 'i')) |> 
  left_join(df2, by = c('sp2', 'season', 'i')) |> 
  mutate(set = map2(hv1, hv2, \(hv1, hv2) hypervolume_set(hv1, hv2, check.memory = F, verbose = F)),
         ov = map(set, \(set) hypervolume_overlap_statistics(set)),
         dist_cent = map2_dbl(hv1, hv2, \(hv1,hv2) hypervolume_distance(hv1, hv2, type = 'centroid', check.memory=F))) |> 
  unnest_wider(ov) |> 
  select(sp1, sp2, season, i,
         jaccard, sorensen,
         uniq_1 = frac_unique_1, uniq_2 = frac_unique_2, 
         dist_cent)

#write_csv(ov_sp, 'data/hvOv_spp.csv')


# Figures and tables ----

## Table 1----
df_t1 = ov_sn |>
  select(species, i, vol_Wet = hv_size_Wet, vol_Dry = hv_size_Dry,
         uniq_Wet, uniq_Dry) |>
  pivot_longer(
    cols = vol_Wet:uniq_Dry,
    names_to = c('.value', 'Season'),
    names_sep = '_'
  ) |>
  mutate(Season = factor(Season, levels = c('Wet', 'Dry'))) |>
  group_by(species, Season) |>
  summarize(
    vmean = mean(vol),
    vlow = quantile(vol, 0.025),
    vup = quantile(vol, 0.975),
    umean = mean(uniq),
    ulow = quantile(uniq, 0.025),
    uup = quantile(uniq, 0.975),
    .groups = 'drop'
  ) |>
  mutate(
    `Niche volume` = paste0(
      format(round(vmean, digits = 2), nsmall = 2,trim = T),
      ' (',
      format(round(vlow, digits = 2), nsmall = 2,trim = T),
      '-',
      format(round(vup, digits = 2), nsmall = 2,trim = T),
      ')'
    ),
    `Percent unique` = paste0(
      format(round(umean, digits = 2), nsmall = 2),
      ' (',
      format(round(ulow, digits = 2), nsmall = 2),
      '-',
      format(round(uup, digits = 2), nsmall = 2),
      ')'
    )
  ) |>
  select(Species = species,
         Season,
         `Niche volume`,
         `Percent unique`)

d_t1 = ov_sn |>
  select(i, species, sorensen, dist_cent) |>
  group_by(species) |>
  summarize(across(c(sorensen, dist_cent),
                   list(
                     mean = ~ mean(.x),
                     low = ~ quantile(.x, 0.025),
                     up = ~ quantile(.x, 0.975)
                   ))) |> 
  mutate(Overlap = paste0(
    format(round(sorensen_mean, digits = 2), nsmall = 2),
    ' (',
    format(round(sorensen_low, digits = 2), nsmall = 2),
    '-',
    format(round(sorensen_up, digits = 2), nsmall = 2),
    ')'),
    `Centroid distance` = paste0(
      format(round(dist_cent_mean, digits = 2), nsmall = 2),
      ' (',
      format(round(dist_cent_low, digits = 2), nsmall = 2),
      '-',
      format(round(dist_cent_up, digits = 2), nsmall = 2),
      ')'),
    Season = 'Wet'
  ) |> 
  select(Species = species, Season, Overlap, `Centroid distance`)

df_t1 = left_join(df_t1,d_t1)

write_excel_csv(df_t1, 'tables/table1.csv', na = '')

## Fig 1----
# a) overlap
df_1a = ov_sn |> 
  group_by(species) |> 
  summarize(mean = mean(sorensen),
            median = median(sorensen),
            low = quantile(sorensen, 0.025),
            up = quantile(sorensen, 0.975))

cols = c("Pinfish" = 'yellow2',
         "Mojarra" = 'slategray4',
         "Silver perch" = 'snow3',
         "Bay anchovy" = 'deepskyblue1',
         "Pigfish" = 'orange', 
         "Pink shrimp" = 'pink',
         "Rainwater killifish" = 'firebrick',
         'All' = 'black')

a = ggplot(df_1a, aes(x = species, y = mean, color = species))+
  geom_pointrange(aes(ymin = low, ymax = up),
                  size = 1.5, linewidth = 1.5, fatten = 2,  
                  position=position_dodge(width = 0.5))+
  labs(x = NULL, y = 'Niche overlap') +
  scale_x_discrete(labels = c("Bay \nanchovy",
                              "Mojarra",
                              "Pigfish",
                              "Pinfish",
                              "Pink \nshrimp",
                              "Rainwater \nkillifish",
                              "Silver \nperch" ))+
  scale_color_manual(values = cols)+
  theme_bw()+
  theme(axis.title = element_text(size = 14), 
        axis.text.y = element_text(size = 14, colour = "black"), 
        axis.text.x = element_text(size = 12, colour = "black"),
        plot.title = element_text(size = 14, hjust=0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = 'right',
        legend.title = element_text(size = 14),
        strip.text.x = element_text(size = 14),
        legend.text = element_text(size = 12))

# b) distance
df_1b = ov_sn |> 
  group_by(species) |> 
  summarize(mean = mean(dist_cent),
            median = median(dist_cent),
            low = quantile(dist_cent, 0.025),
            up = quantile(dist_cent, 0.975))

b = ggplot(df_1b, aes(x = species, y = mean, color = species))+
  geom_pointrange(aes(ymin = low, ymax = up),
                  size = 1.5, linewidth = 1.5, fatten = 2, 
                  position=position_dodge(width = 0.5))+
  labs(x = 'Species', y = 'Centroid distance') +
  scale_x_discrete(labels = c("Bay \nanchovy",
                              "Mojarra",
                              "Pigfish",
                              "Pinfish",
                              "Pink \nshrimp",
                              "Rainwater \nkillifish",
                              "Silver \nperch" ))+
  scale_color_manual(values = cols)+
  theme_bw()+
  theme(axis.title = element_text(size = 14), 
        axis.text.y = element_text(size = 14, colour = "black"), 
        axis.text.x = element_text(size = 12, colour = "black"),
        plot.title = element_text(size = 14, hjust=0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = 'right',
        legend.title = element_text(size = 14),
        strip.text.x = element_text(size = 14),
        legend.text = element_text(size = 12))

ggarrange(a,b,
          labels = c('a)','b)'),
          nrow = 2, 
          align = 'hv',
          legend = 'none')

ggsave("figs/fig1.png", units="in", width=7, height=8, dpi=600)
## Fig 2 ----
df_2 = df |> 
  select(species, season, i, centroid) |> 
  unnest_longer(centroid) |> 
  pivot_wider(names_from = season, values_from = centroid)|> 
  ungroup() |> 
  mutate(dif = Dry - Wet) |> 
  group_by(species, centroid_id) |> 
  summarize(mean = mean(dif),
            low = quantile(dif, 0.025),
            up = quantile(dif, 0.975),
            .groups = 'drop')

cols = c("Pinfish" = 'yellow2',
         "Mojarra" = 'slategray4',
         "Silver perch" = 'snow3',
         "Bay anchovy" = 'deepskyblue1',
         "Pigfish" = 'orange', 
         "Pink shrimp" = 'Pink',
         "Rainwater killifish" = 'firebrick')

ggplot(df_2, aes(x = centroid_id, y = mean, color = species))+
  geom_hline(aes(yintercept = 0), 
             linetype = 'dashed', linewidth = 1.5)+
  geom_pointrange(aes(ymin = low, ymax = up),
                  size = 1.5, linewidth = 1.5, fatten = 2, 
                  position=position_dodge(width = 1))+
  labs(x = NULL, y = 'Centroid difference',
       color = 'Species') +
  scale_color_manual(values = cols)+
  theme_bw()+
  facet_grid(cols = vars(centroid_id), scales = 'free_x')+
  #facet_wrap(~centroid_id, nrow = 2)+
  theme(axis.title = element_text(size = 14), 
        axis.text.y = element_text(size = 14, colour = "black"), 
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        # axis.text.x = element_text(size = 12, colour = "black"),
        plot.title = element_text(size = 14, hjust=0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = 'right',
        legend.title = element_text(size = 14),
        strip.text = element_text(size = 14),
        legend.text = element_text(size = 12))

ggsave("figs/fig2.png", units="in", width=8, height=4, dpi=600)

## Fig 3----
### species niches ----
df_3ab = df |> 
  mutate(hv_df = map(hv, \(data) hypervolume_to_data_frame(data))) |> 
  select(i,species, season, hv_df) |> 
  unnest(hv_df) |> 
  mutate(across(Algae:Seagrass, ~round(.x, digits = 2))) |> 
  group_by(species, season, Algae, Epiphytes, Mangrove, Seagrass) |> 
  summarize(ValueAtRandomPoints = mean(ValueAtRandomPoints), .groups = 'drop')

# get points for plotting
do = df_3ab |> 
  group_by(species, season) |> 
  sample_n(5000) |>
  mutate(num = row_number()) |> 
  pivot_longer(Algae:Seagrass, names_to = 'axis', values_to = 'value')

do1 = do |> 
  rename(axis1 = axis, val1 = value)

do2 = do |> 
  rename(axis2 = axis, val2 = value)


spc = tibble(axis1 = unique(do$axis),
             axis2 = unique(do$axis))

df_ax = spc |>  expand(axis1,axis2)

df_ax = df_ax[!duplicated(t(apply(df_ax,1,sort))),] |>
  filter(!(axis1 == axis2))

df_axw = df_ax |>
  mutate(season = 'Wet') |> 
  slice(rep(1:n(), each=length(unique(do$species)))) |> 
  mutate(species = rep(unique(do$species), times = nrow(df_ax)))

df_axd = df_ax |>
  mutate(season = 'Dry')|> 
  slice(rep(1:n(), each=length(unique(do$species)))) |> 
  mutate(species = rep(unique(do$species), times = nrow(df_ax)))

set.seed(14)
df_axx = bind_rows(df_axw, df_axd) |> 
  slice(rep(1:n(), each=60000))|> 
  group_by(species, season, axis1, axis2) |> 
  mutate(num = row_number()) |> 
  left_join(do1) |> 
  left_join(do2)

# centroid
ax = bind_rows(df_axw, df_axd) 

cent = df |> 
  select(species, season, i, centroid) |> 
  unnest_wider(centroid) |> 
  group_by(species, season) |> 
  summarize(across(Algae:Seagrass, mean),
            .groups = 'drop') |> 
  pivot_longer(Algae:Seagrass, names_to = 'axis', values_to = 'value')

c1 = cent |> 
  rename(axis1 = axis, val1 = value)

c2 = cent |> 
  rename(axis2 = axis, val2 = value)

df_cent = ax |>  
  left_join(c1) |> 
  left_join(c2)

cols = c("Pinfish" = 'yellow2',
         "Mojarra" = 'slategray4',
         "Silver perch" = 'snow3',
         "Bay anchovy" = 'deepskyblue1',
         "Pigfish" = 'orange', 
         "Pink shrimp" = 'Pink',
         "Rainwater killifish" = 'firebrick')

#### a) wet season ---- 
hw = ggplot()+
  geom_point(data = df_axx[sample(nrow(df_axx)), ] |>  filter(season == 'Wet'),
             aes(val2, val1, color = species), #alpha = 0.5, 
             size = 1.5) +
  geom_point(data = df_cent |> filter(season == 'Wet'), 
             aes(val2, val1, fill = species), shape = 21, 
             color ='white', size = 5, stroke = 1)+
  scale_x_continuous(limits = c(-5,5), breaks = c(-4,-2,0,2,4))+
  scale_y_continuous(limits = c(-5,5), breaks = c(-4,-2,0,2,4))+
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols)+
  labs(x = NULL, y = NULL, fill = 'Species', color = 'Occupancy')+
  theme_bw()+
  facet_grid(cols = vars(axis2), rows = vars(axis1), switch = 'both')+
  guides(color = 'none')+
  theme(axis.title = element_text(size = 14), 
        axis.text.y = element_text(size = 13, colour = "black"), 
        axis.text.x = element_text(size = 13, colour = "black"),
        plot.title = element_text(size = 14, hjust=0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = 'right',
        legend.title = element_text(size = 17),
        strip.text = element_text(size = 17),
        legend.text = element_text(size = 16),
        strip.text.y.left = element_text(angle = 0),
        strip.background = element_blank(),
        strip.placement = "outside")


# fix all labels
gt = ggplotGrob(hw)


del = which(gt$layout$name %in% c('panel-1-2', 'panel-1-3', 'panel-2-3',
                                  'strip-b-1', 'strip-b-2'))
gt$grobs[del] = NULL
gt$layout = gt$layout[-del, ]

# y axis row 2
gt$layout$l[gt$layout$name == 'axis-l-2'] = gt$layout$l[gt$layout$name == 'axis-l-2'] + 2
gt$layout$r[gt$layout$name == 'axis-l-2'] = gt$layout$r[gt$layout$name == 'axis-l-2'] + 2

# y axis row 3
gt$layout$l[gt$layout$name == 'axis-l-3'] = gt$layout$l[gt$layout$name == 'axis-l-3'] + 4
gt$layout$r[gt$layout$name == 'axis-l-3'] = gt$layout$r[gt$layout$name == 'axis-l-3'] + 4

# x axis column 1
gt$layout$t[gt$layout$name == 'axis-b-1'] = gt$layout$t[gt$layout$name == 'axis-b-1'] - 4
gt$layout$b[gt$layout$name == 'axis-b-1'] = gt$layout$b[gt$layout$name == 'axis-b-1'] - 4

# x axis column 2
gt$layout$t[gt$layout$name == 'axis-b-2'] = gt$layout$t[gt$layout$name == 'axis-b-2'] - 2
gt$layout$b[gt$layout$name == 'axis-b-2'] = gt$layout$b[gt$layout$name == 'axis-b-2'] - 2


# y strip row 2
gt$layout$l[gt$layout$name == 'strip-l-2'] = gt$layout$l[gt$layout$name == 'strip-l-2'] + 3
gt$layout$r[gt$layout$name == 'strip-l-2'] = gt$layout$r[gt$layout$name == 'strip-l-2'] + 3

# y strip row 3
gt$layout$l[gt$layout$name == 'strip-l-3'] = gt$layout$l[gt$layout$name == 'strip-l-3'] + 5
gt$layout$r[gt$layout$name == 'strip-l-3'] = gt$layout$r[gt$layout$name == 'strip-l-3'] + 5


# save 
png("figs/fig3a.png", width = 8.5, height = 5.5, 
     units = 'in', res = 600)

grid.newpage()
grid.draw(gt)

dev.off()



#### b) dry season ----
hd = ggplot()+
  geom_point(data = df_axx[sample(nrow(df_axx)), ] |>  
               filter(season == 'Dry'),
             aes(val2, val1, color = species), #alpha = 0.5, 
             size = 1.5) +
  geom_point(data = df_cent |> filter(season == 'Dry'), 
             aes(val2, val1, fill = species), shape = 21, 
             color ='white', size = 5, stroke = 1)+
  scale_x_continuous(limits = c(-5,5), breaks = c(-4,-2,0,2,4))+
  scale_y_continuous(limits = c(-5,5), breaks = c(-4,-2,0,2,4))+
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols)+
  labs(x = NULL, y = NULL, fill = 'Species', color = 'Occupancy')+
  theme_bw()+
  facet_grid(cols = vars(axis2), rows = vars(axis1), switch = 'both')+
  guides(color = 'none')+
  theme(axis.title = element_text(size = 14), 
        axis.text.y = element_text(size = 13, colour = "black"), 
        axis.text.x = element_text(size = 13, colour = "black"),
        plot.title = element_text(size = 14, hjust=0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = 'right',
        legend.title = element_text(size = 17),
        strip.text = element_text(size = 17),
        legend.text = element_text(size = 16),
        strip.text.y.left = element_text(angle = 0),
        strip.background = element_blank(),
        strip.placement = "outside")

# fix all labels
gt = ggplotGrob(hd)


del = which(gt$layout$name %in% c('panel-1-2', 'panel-1-3', 'panel-2-3',
                                  'strip-b-1', 'strip-b-2'))
gt$grobs[del] = NULL
gt$layout = gt$layout[-del, ]

# y axis row 2
gt$layout$l[gt$layout$name == 'axis-l-2'] = gt$layout$l[gt$layout$name == 'axis-l-2'] + 2
gt$layout$r[gt$layout$name == 'axis-l-2'] = gt$layout$r[gt$layout$name == 'axis-l-2'] + 2

# y axis row 3
gt$layout$l[gt$layout$name == 'axis-l-3'] = gt$layout$l[gt$layout$name == 'axis-l-3'] + 4
gt$layout$r[gt$layout$name == 'axis-l-3'] = gt$layout$r[gt$layout$name == 'axis-l-3'] + 4

# x axis column 1
gt$layout$t[gt$layout$name == 'axis-b-1'] = gt$layout$t[gt$layout$name == 'axis-b-1'] - 4
gt$layout$b[gt$layout$name == 'axis-b-1'] = gt$layout$b[gt$layout$name == 'axis-b-1'] - 4

# x axis column 2
gt$layout$t[gt$layout$name == 'axis-b-2'] = gt$layout$t[gt$layout$name == 'axis-b-2'] - 2
gt$layout$b[gt$layout$name == 'axis-b-2'] = gt$layout$b[gt$layout$name == 'axis-b-2'] - 2


# y strip row 2
gt$layout$l[gt$layout$name == 'strip-l-2'] = gt$layout$l[gt$layout$name == 'strip-l-2'] + 3
gt$layout$r[gt$layout$name == 'strip-l-2'] = gt$layout$r[gt$layout$name == 'strip-l-2'] + 3

# y strip row 3
gt$layout$l[gt$layout$name == 'strip-l-3'] = gt$layout$l[gt$layout$name == 'strip-l-3'] + 5
gt$layout$r[gt$layout$name == 'strip-l-3'] = gt$layout$r[gt$layout$name == 'strip-l-3'] + 5


# save 
png("figs/fig3b.png", width = 8.5, height = 5.5, 
     units = 'in', res = 600)

grid.newpage()
grid.draw(gt)

dev.off()

### occupancy ----
df_oc = df_occ  |> 
  select(i,hv_df) |> 
  unnest(hv_df) |> 
  rename(season = Name, occupancy = ValueAtRandomPoints) |> 
  mutate(across(Algae:Seagrass, ~round(.x, digits = 2))) |> 
  group_by(season, Algae, Epiphytes, Mangrove, Seagrass) |> 
  summarize(occupancy = mean(occupancy), .groups = 'drop')

# get points for plotting
set.seed(14)
do = df_oc |> 
  mutate(cat = if_else(occupancy > 1/7, 'Multi species', 'Single species')) |> 
  #filter(ValueAtRandomPoints > 1/7) |> 
  group_by(season) |> 
  sample_n(5000) |>
  mutate(num = row_number()) |> 
  pivot_longer(Algae:Seagrass, names_to = 'axis', values_to = 'value')

do1 = do |> 
  rename(axis1 = axis, val1 = value)

do2 = do |> 
  rename(axis2 = axis, val2 = value)


spc = tibble(axis1 = unique(do$axis),
             axis2 = unique(do$axis))

df_ax = spc |>  expand(axis1,axis2)

df_ax = df_ax[!duplicated(t(apply(df_ax,1,sort))),] |>
  filter(!(axis1 == axis2))

df_axw = df_ax |>
  mutate(season = 'Wet')

df_axd = df_ax |>
  mutate(season = 'Dry')

set.seed(14)
df_axx = bind_rows(df_axw, df_axd) |> 
  slice(rep(1:n(), each=60000))|> 
  group_by(season, axis1, axis2) |> 
  mutate(num = row_number()) |> 
  left_join(do1) |> 
  left_join(do2)

#### c) wet season---- 
pw = ggplot()+
  geom_point(data = df_axx |>  filter(cat == 'Single species',
                                      season == 'Wet'),
             aes(val2, val1, color = cat), size = 1.5) +
  geom_point(data = df_axx |>  filter(cat == 'Multi species',
                                      season == 'Wet'),
             aes(val2, val1, color = cat), size = 1.5) +
  # geom_point(data = df_cent |> filter(season == 'Wet'), 
  #            aes(val2, val1, fill = species), shape = 21, 
  #            color ='white', size = 5, stroke = 1)+
  scale_x_continuous(limits = c(-5,5), breaks = c(-4,-2,0,2,4))+
  scale_y_continuous(limits = c(-5,5), breaks = c(-4,-2,0,2,4))+
  scale_color_manual(values = c('Multi species' = "#132B43",
                                'Single species' ="#add8e6")) +
  #scale_fill_manual(values = cols)+
  labs(x = NULL, y = NULL, fill = 'Species', color = 'Occupancy')+
  theme_bw()+
  facet_grid(cols = vars(axis2), rows = vars(axis1), switch = 'both')+
  theme(axis.title = element_text(size = 14), 
        axis.text.y = element_text(size = 13, colour = "black"), 
        axis.text.x = element_text(size = 13, colour = "black"),
        plot.title = element_text(size = 14, hjust=0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = 'right',
        legend.title = element_text(size = 17),
        strip.text = element_text(size = 17),
        legend.text = element_text(size = 16),
        strip.text.y.left = element_text(angle = 0),
        strip.background = element_blank(),
        strip.placement = "outside")


# fix all labels
gt = ggplotGrob(pw)


del = which(gt$layout$name %in% c('panel-1-2', 'panel-1-3', 'panel-2-3',
                                  'strip-b-1', 'strip-b-2'))
gt$grobs[del] = NULL
gt$layout = gt$layout[-del, ]

# y axis row 2
gt$layout$l[gt$layout$name == 'axis-l-2'] = gt$layout$l[gt$layout$name == 'axis-l-2'] + 2
gt$layout$r[gt$layout$name == 'axis-l-2'] = gt$layout$r[gt$layout$name == 'axis-l-2'] + 2

# y axis row 3
gt$layout$l[gt$layout$name == 'axis-l-3'] = gt$layout$l[gt$layout$name == 'axis-l-3'] + 4
gt$layout$r[gt$layout$name == 'axis-l-3'] = gt$layout$r[gt$layout$name == 'axis-l-3'] + 4

# x axis column 1
gt$layout$t[gt$layout$name == 'axis-b-1'] = gt$layout$t[gt$layout$name == 'axis-b-1'] - 4
gt$layout$b[gt$layout$name == 'axis-b-1'] = gt$layout$b[gt$layout$name == 'axis-b-1'] - 4

# x axis column 2
gt$layout$t[gt$layout$name == 'axis-b-2'] = gt$layout$t[gt$layout$name == 'axis-b-2'] - 2
gt$layout$b[gt$layout$name == 'axis-b-2'] = gt$layout$b[gt$layout$name == 'axis-b-2'] - 2


# y strip row 2
gt$layout$l[gt$layout$name == 'strip-l-2'] = gt$layout$l[gt$layout$name == 'strip-l-2'] + 3
gt$layout$r[gt$layout$name == 'strip-l-2'] = gt$layout$r[gt$layout$name == 'strip-l-2'] + 3

# y strip row 3
gt$layout$l[gt$layout$name == 'strip-l-3'] = gt$layout$l[gt$layout$name == 'strip-l-3'] + 5
gt$layout$r[gt$layout$name == 'strip-l-3'] = gt$layout$r[gt$layout$name == 'strip-l-3'] + 5

# save
png("figs/fig3c.png", width = 8.5, height = 5.5, 
     units = 'in', res = 600)

grid.newpage()
grid.draw(gt)

dev.off()

#### d) dry season ----
pd = ggplot()+
  geom_point(data = df_axx |>  filter(cat == 'Single species',
                                      season == 'Dry'),
             aes(val2, val1, color = cat), size = 1.5) +
  geom_point(data = df_axx |>  filter(cat == 'Multi species',
                                      season == 'Dry'),
             aes(val2, val1, color = cat), size = 1.5) +
  # geom_point(data = df_cent |> filter(season == 'Dry'), 
  #            aes(val2, val1, fill = species), shape = 21, 
  #            color ='white', size = 5, stroke = 1)+
  scale_x_continuous(limits = c(-5,5), breaks = c(-4,-2,0,2,4))+
  scale_y_continuous(limits = c(-5,5), breaks = c(-4,-2,0,2,4))+
  scale_color_manual(values = c('Multi species' = "#132B43",
                                'Single species' ="#add8e6")) +
  #scale_fill_manual(values = cols)+
  labs(x = NULL, y = NULL, fill = 'Species', color = 'Occupancy')+
  theme_bw()+
  facet_grid(cols = vars(axis2), rows = vars(axis1), switch = 'both')+
  theme(axis.title = element_text(size = 14), 
        axis.text.y = element_text(size = 13, colour = "black"), 
        axis.text.x = element_text(size = 13, colour = "black"),
        plot.title = element_text(size = 14, hjust=0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = 'right',
        legend.title = element_text(size = 17),
        strip.text = element_text(size = 17),
        legend.text = element_text(size = 16),
        strip.text.y.left = element_text(angle = 0),
        strip.background = element_blank(),
        strip.placement = "outside")

# fix all labels
gt = ggplotGrob(pd)

del = which(gt$layout$name %in% c('panel-1-2', 'panel-1-3', 'panel-2-3',
                                  'strip-b-1', 'strip-b-2'))
gt$grobs[del] = NULL
gt$layout = gt$layout[-del, ]

# y axis row 2
gt$layout$l[gt$layout$name == 'axis-l-2'] = gt$layout$l[gt$layout$name == 'axis-l-2'] + 2
gt$layout$r[gt$layout$name == 'axis-l-2'] = gt$layout$r[gt$layout$name == 'axis-l-2'] + 2

# y axis row 3
gt$layout$l[gt$layout$name == 'axis-l-3'] = gt$layout$l[gt$layout$name == 'axis-l-3'] + 4
gt$layout$r[gt$layout$name == 'axis-l-3'] = gt$layout$r[gt$layout$name == 'axis-l-3'] + 4

# x axis column 1
gt$layout$t[gt$layout$name == 'axis-b-1'] = gt$layout$t[gt$layout$name == 'axis-b-1'] - 4
gt$layout$b[gt$layout$name == 'axis-b-1'] = gt$layout$b[gt$layout$name == 'axis-b-1'] - 4

# x axis column 2
gt$layout$t[gt$layout$name == 'axis-b-2'] = gt$layout$t[gt$layout$name == 'axis-b-2'] - 2
gt$layout$b[gt$layout$name == 'axis-b-2'] = gt$layout$b[gt$layout$name == 'axis-b-2'] - 2


# y strip row 2
gt$layout$l[gt$layout$name == 'strip-l-2'] = gt$layout$l[gt$layout$name == 'strip-l-2'] + 3
gt$layout$r[gt$layout$name == 'strip-l-2'] = gt$layout$r[gt$layout$name == 'strip-l-2'] + 3

# y strip row 3
gt$layout$l[gt$layout$name == 'strip-l-3'] = gt$layout$l[gt$layout$name == 'strip-l-3'] + 5
gt$layout$r[gt$layout$name == 'strip-l-3'] = gt$layout$r[gt$layout$name == 'strip-l-3'] + 5

# save
png("figs/fig3d.png", width = 8.5, height = 5.5, 
     units = 'in', res = 600)

grid.newpage()
grid.draw(gt)

dev.off()

## Fig 4----
### a) niche size ----
sp = ov_sn |> 
  rename(Wet = hv_size_Wet, Dry = hv_size_Dry) |> 
  pivot_longer(Wet:Dry,names_to = 'season',
               values_to = 'vol') |> 
  select(species, season, i, vol)

ms = oc_s |> 
  filter(cat == 'Multi species') |> 
  select(species = cat, season, i, vol = size_occ)

all = oc_s |> 
  filter(cat == 'Multi species') |> 
  mutate(species = 'Community') |> 
  select(species, season, i, vol = size)

ctn = all |> 
  select(season, i, cnb = vol)

df_4a = bind_rows(sp, ms, all) |> 
  group_by(species, season) |> 
  summarize(mean = mean(vol),
            median = median(vol),
            low = quantile(vol, 0.025),
            up = quantile(vol, 0.975),
            .groups = 'drop') |> 
  mutate(season = factor(season, levels = c('Wet', 
                                            'Dry')),
         species = factor(species, levels = c('Community',
                                              'Multi species',
                                              "Bay anchovy",
                                              "Mojarra",
                                              "Pigfish",
                                              "Pinfish",
                                              "Pink shrimp",
                                              "Rainwater killifish",
                                              "Silver perch")))

a = ggplot(df_4a, aes(x = species, y = mean, color = season))+
  geom_pointrange(aes(ymin = low, ymax = up),
                  size = 1.5, linewidth = 1.5, fatten = 2, 
                  position=position_dodge(width = 0.5))+
  labs(x = 'Species', y = 'Trophic niche width',
       color = 'Season') +
  # scale_fill_manual(values = c('Wet' = 'skyblue3', 
  #                              'Dry' = 'indianred3')) +
  scale_color_manual(values = c('Wet' = 'skyblue3', 
                                'Dry' = 'indianred3')) +
  scale_x_discrete(labels = c("Community",
                              "Multi \nspecies",
                              "Bay \nanchovy",
                              "Mojarra",
                              "Pigfish",
                              "Pinfish",
                              "Pink \nshrimp",
                              "Rainwater \nkillifish",
                              "Silver \nperch" ))+
  theme_bw()+
  theme(axis.title = element_text(size = 14), 
        axis.text.y = element_text(size = 14, colour = "black"), 
        axis.text.x = element_text(size = 12, colour = "black"),
        plot.title = element_text(size = 14, hjust=0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = 'right',
        legend.title = element_text(size = 14),
        strip.text.x = element_text(size = 14),
        legend.text = element_text(size = 12))

### b) proportion of community niche----
sp = ov_sn |> 
  rename(Wet = hv_size_Wet, Dry = hv_size_Dry) |> 
  pivot_longer(Wet:Dry,names_to = 'season',
               values_to = 'vol') |> 
  select(species, season, i, vol)

ms = oc_s |> 
  filter(cat == 'Multi species') |> 
  select(species = cat, season, i, vol = size_occ)

all = oc_s |> 
  filter(cat == 'Multi species') |> 
  mutate(species = 'Community') |> 
  select(species, season, i, vol = size)

ctn = all |> 
  select(season, i, cnb = vol)

df_4b = bind_rows(sp, ms, all) |> 
  left_join(ctn) |> 
  mutate(prop = vol/cnb) |> 
  group_by(species, season) |> 
  filter(species != 'Community') |> 
  summarize(mean = mean(prop),
            median = median(prop),
            low = quantile(prop, 0.025),
            up = quantile(prop, 0.975),
            .groups = 'drop') |> 
  mutate(season = factor(season, levels = c('Wet', 
                                            'Dry')),
         species = factor(species, levels = c('Multi species',
                                              "Bay anchovy",
                                              "Mojarra",
                                              "Pigfish",
                                              "Pinfish",
                                              "Pink shrimp",
                                              "Rainwater killifish",
                                              "Silver perch")))

b = ggplot(df_4b, aes(x = species, y = mean, color = season))+
  geom_pointrange(aes(ymin = low, ymax = up),
                  size = 1.5, linewidth = 1.5, fatten = 2, 
                  position=position_dodge(width = 0.5))+
  labs(x = 'Species', y = 'Proportion of community \ntrophic niche breadth',
       color = 'Season') +
  # scale_fill_manual(values = c('Wet' = 'skyblue3', 
  #                              'Dry' = 'indianred3')) +
  scale_color_manual(values = c('Wet' = 'skyblue3', 
                                'Dry' = 'indianred3')) +
  scale_x_discrete(labels = c("Multi \nspecies",
                              "Bay \nanchovy",
                              "Mojarra",
                              "Pigfish",
                              "Pinfish",
                              "Pink \nshrimp",
                              "Rainwater \nkillifish",
                              "Silver \nperch" ))+
  theme_bw()+
  theme(axis.title = element_text(size = 14), 
        axis.text.y = element_text(size = 14, colour = "black"), 
        axis.text.x = element_text(size = 12, colour = "black"),
        plot.title = element_text(size = 14, hjust=0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = 'right',
        legend.title = element_text(size = 14),
        strip.text.x = element_text(size = 14),
        legend.text = element_text(size = 12))


ggarrange(a,b,
          labels = c('a)','b)'),
          nrow = 2, 
          align = 'hv')

ggsave("figs/fig4.png", units="in", width=9, height=9, dpi=600)

## Fig 5 ----
### a) overlap ----
df_5a = ov_sp |>
  group_by(sp1, sp2, season) |> 
  summarize(mean = mean(sorensen),
            median = median(sorensen),
            low = quantile(sorensen, 0.025),
            up = quantile(sorensen, 0.975),
            .groups = 'drop')|> 
  mutate(season = factor(season, levels = c('Wet', 
                                            'Dry')))

d_5a = df_5a |> rename(Sp2 = sp1, Sp1 = sp2) 

df_5a = df_5a |> rename(Sp1 = sp1, Sp2 = sp2) |> 
  bind_rows(d_5a)

m_5a = ov_sp |>
  group_by(season) |> 
  summarize(mean = mean(sorensen),
            median = median(sorensen),
            low = quantile(sorensen, 0.025),
            up = quantile(sorensen, 0.975))|> 
  mutate(season = factor(season, levels = c('Wet', 
                                            'Dry')))

cols = c("Pinfish" = 'yellow2',
         "Mojarra" = 'slategray4',
         "Silver perch" = 'snow3',
         "Bay anchovy" = 'deepskyblue1',
         "Pigfish" = 'orange', 
         "Pink shrimp" = 'Pink',
         "Rainwater killifish" = 'firebrick')

a = ggplot(df_5a, aes(x = Sp1, y = mean, color = Sp2))+
  geom_hline(data = m_5a, aes(yintercept = mean), 
             linetype = 'dashed', linewidth = 1.5)+
  geom_pointrange(aes(ymin = low, ymax = up),
                  size = 1.5, linewidth = 1.5, fatten = 2, 
                  position=position_dodge(width = 0.5))+
  labs(x = NULL, y = 'Niche overlap', color = 'Species') +
  scale_x_discrete(labels = c("Bay \nanchovy",
                              "Mojarra",
                              "Pigfish",
                              "Pinfish",
                              "Pink \nshrimp",
                              "Rainwater \nkillifish",
                              "Silver \nperch" ))+
  scale_color_manual(values = cols)+
  theme_bw()+
  facet_grid(row = vars(season))+
  theme(axis.title = element_text(size = 14), 
        axis.text.y = element_text(size = 14, colour = "black"), 
        axis.text.x = element_text(size = 12, colour = "black"),
        plot.title = element_text(size = 14, hjust=0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = 'right',
        legend.title = element_text(size = 14),
        strip.text.y = element_text(size = 14),
        legend.text = element_text(size = 12))

### b) centroid distance ----
df_5b = ov_sp |>
  group_by(sp1, sp2, season) |> 
  summarize(mean = mean(dist_cent),
            median = median(dist_cent),
            low = quantile(dist_cent, 0.025),
            up = quantile(dist_cent, 0.975),
            .groups = 'drop')|> 
  mutate(season = factor(season, levels = c('Wet', 
                                            'Dry')))

d_5b = df_5b |> rename(Sp2 = sp1, Sp1 = sp2) 

df_5b = df_5b |> rename(Sp1 = sp1, Sp2 = sp2) |> 
  bind_rows(d_5b)

m_5b = ov_sp |>
  group_by(season) |> 
  summarize(mean = mean(dist_cent),
            median = median(dist_cent),
            low = quantile(dist_cent, 0.025),
            up = quantile(dist_cent, 0.975))|> 
  mutate(season = factor(season, levels = c('Wet', 
                                            'Dry')))

cols = c("Pinfish" = 'yellow2',
         "Mojarra" = 'slategray4',
         "Silver perch" = 'snow3',
         "Bay anchovy" = 'deepskyblue1',
         "Pigfish" = 'orange', 
         "Pink shrimp" = 'Pink',
         "Rainwater killifish" = 'firebrick')

b = ggplot(df_5b, aes(x = Sp1, y = mean, color = Sp2))+
  geom_hline(data = m_5b, aes(yintercept = mean), 
             linetype = 'dashed', linewidth = 1.5)+
  geom_pointrange(aes(ymin = low, ymax = up),
                  size = 1.5, linewidth = 1.5, fatten = 2,
                  position=position_dodge(width = 0.5))+
  labs(x = 'Species', y = 'Centroid distance', color = 'Species') +
  scale_x_discrete(labels = c("Bay \nanchovy",
                              "Mojarra",
                              "Pigfish",
                              "Pinfish",
                              "Pink \nshrimp",
                              "Rainwater \nkillifish",
                              "Silver \nperch" ))+
  scale_color_manual(values = cols)+
  theme_bw()+
  facet_grid(row = vars(season))+
  theme(axis.title = element_text(size = 14), 
        axis.text.y = element_text(size = 14, colour = "black"), 
        axis.text.x = element_text(size = 12, colour = "black"),
        plot.title = element_text(size = 14, hjust=0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = 'right',
        legend.title = element_text(size = 14),
        strip.text.y = element_text(size = 14),
        legend.text = element_text(size = 12))

ggarrange(a,b,
          labels = c('a)','b)'),
          nrow = 2, 
          align = 'hv',
          common.legend = T,
          legend = 'right')

ggsave("figs/fig5.png", units="in", width=10, height=9, dpi=600)

#Supplemental ----
## Table S1 ----
df_ts1 = read_csv('CESImixResults.csv') |> 
  mutate(val = paste0(format(round(mean, digits = 2),nsmall = 2),
                      ' \u00B1 ',
                      format(round(sd, digits = 2),nsmall = 2),
                      ' (',
                      format(round(lowend, digits = 2),nsmall = 2),
                      '-',
                      format(round(highend, digits = 2),nsmall = 2),
                      ')')) |> 
  select(Species = species, Season = season, source, val) |> 
  pivot_wider(names_from = 'source', values_from = 'val') |> 
  arrange(Species)

write_excel_csv(df_ts1, 'tables/tables1.csv')

## Table S2----
df_ts2 = ov_sp |> 
  mutate(season = factor(season, levels = c('Wet', 'Dry'))) |> 
  group_by(sp1,sp2,season) |> 
  summarize(across(c(sorensen, dist_cent),
                   list(
                     mean = ~ mean(.x),
                     low = ~ quantile(.x, 0.025),
                     up = ~ quantile(.x, 0.975)
                   )),
            .groups = 'drop')|> 
  mutate(Overlap = paste0(
    format(round(sorensen_mean, digits = 2), nsmall = 2),
    ' (',
    format(round(sorensen_low, digits = 2), nsmall = 2),
    '-',
    format(round(sorensen_up, digits = 2), nsmall = 2),
    ')'),
    `Centroid distance` = paste0(
      format(round(dist_cent_mean, digits = 2), nsmall = 2),
      ' (',
      format(round(dist_cent_low, digits = 2), nsmall = 2),
      '-',
      format(round(dist_cent_up, digits = 2), nsmall = 2),
      ')'),
    Season = 'Wet'
  ) |> 
  select(`Species 1` = sp1, `Species 2` = sp2,
         Season = season,Overlap, `Centroid distance`)

write_excel_csv(df_ts2, 'tables/tableS2.csv', na = '')

## Table S3----
df |> 
  select(species, season, i, centroid) |> 
  unnest_longer(centroid) |> 
  pivot_wider(names_from = season, values_from = centroid) |> 
  ungroup() |> 
  mutate(dist = Dry - Wet) |> 
  group_by(species, centroid_id) |> 
  summarize(mean = mean(dist),
            low = quantile(dist, 0.025),
            up = quantile(dist, 0.975),
            .groups = 'drop') |> 
  mutate(
    `diff` = paste0(
      format(round(mean, digits = 2), nsmall = 2,trim = T),
      ' (',
      format(round(low, digits = 2), nsmall = 2,trim = T),
      '-',
      format(round(up, digits = 2), nsmall = 2,trim = T),
      ')'
    )) |> 
  select(species,centroid_id, diff) |> 
  pivot_wider(names_from = centroid_id, values_from = diff) |> 
  write_excel_csv('tables/tableS3.csv', na = '')

## Table S3 ----
dt_all |> 
  mutate(cluster = as.factor(kmeans_result$cluster)) |> 
  group_by(cluster) |> 
  summarize(across(where(is.numeric), mean),
            .groups = 'drop') |> 
  mutate(across(where(is.numeric), \(x) format(round(x, digits = 2), nsmall = 2,trim = T)))|> 
  select(Cluster = cluster, `Brown difference` = brown, `Green difference` = green,
         `Size ratio` = size_rat, `Centroid distance` = dist_cent, 
         `Wet unique` = uniq_wet, `Dry unique` = uniq_dry, Overlap = sorensen) |> 
  write_excel_csv('tables/tableS3.csv', na = '')

## Fig S2 ----
pca_scores = pca_scores |> 
  mutate(cluster = as.factor(cluster))

# loadings 
loadings = pca_result$rotation
loadings_scaled =loadings * max(abs(pca_scores$PC1), abs(pca_scores$PC2)) / max(abs(loadings))


### set species colors
cols = c("Pinfish" = 'yellow2',
         "Mojarra" = 'slategray4',
         "Silver perch" = 'snow3',
         "Bay anchovy" = 'deepskyblue1',
         "Pigfish" = 'orange', 
         "Pink shrimp" = 'Pink',
         "Rainwater killifish" = 'firebrick')

labels = c("Brown difference", "Green difference", 
           "Size ratio", "Centroid distance",
           "Wet unique", "Dry unique", 'Overlap')

### plot it out
ggplot(pca_scores, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = species)) +
  geom_mark_ellipse(
    aes(fill = cluster, group = cluster),
    alpha = 0.1,
    color = NA,
    show.legend = FALSE,
    expand = unit(2, "mm")
  ) +
  labs(x = "PC1", y = "PC2") +
  scale_color_manual(values = cols) +
  scale_fill_viridis_d(option = 'viridis')+
  theme_bw() +
  labs(color = 'Species')+
  scale_x_continuous(limits = c(-4.5,6))+
  # guides(color = guide_legend("Species")) +
  geom_segment(aes(x = 0, y = 0, xend = PC1, yend = PC2), 
               data = as.data.frame(loadings_scaled), 
               arrow = arrow(type = "closed", length = unit(0.1, "inches")),
               color = "black") +
  geom_text(data = as.data.frame(loadings_scaled), 
            aes(x = PC1, y = PC2, label = labels), 
            nudge_x = -0.35, nudge_y = 0.4, check_overlap = TRUE, color = "black") +
  theme(axis.title = element_text(size = 14), 
        axis.text.y = element_text(size = 14, colour = "black"), 
        # axis.text.x = element_blank(),
        # axis.ticks.x = element_blank(),
        axis.text.x = element_text(size = 12, colour = "black"),
        plot.title = element_text(size = 14, hjust=0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = 'right',
        legend.title = element_text(size = 14),
        strip.text = element_text(size = 14),
        legend.text = element_text(size = 12))

# ggsave("figs/clustered_grouped.tiff", units="in", width=8, height=4, dpi=600,compression = 'lzw')
ggsave("figs/figS2.png", units="in", width=7, height=4, dpi=600)
