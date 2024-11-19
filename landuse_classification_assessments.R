#Code written by Cynthia L. Norton, Univserity of Arizona, 2021 
#The code was used to classify NEON scenes for species specific woody vegetation using LiDAR and Hyperspectral data.
#Script written by Cynthia Norton
# load the raster and rgdal libraries
library(raster)
library(tidyverse)
library(rgdal)
library(dplyr)
library(exactextractr)
library(sf)
library(tidyr)
library(terra)
library(tibble)
library(reshape2)
library(ggplot2)
library(ExtractTrainData)
library(rLiDAR)
library(lidR)
#loadlibraries
library(doParallel)
library(raster)
library(foreach)
library("ForestTools")
#install.packages('ggplot')
#library(ggplot)

setwd("//gaea/projects/RaBET/RaBET_landuse/landuse/")
options(scipen = 100, digits = 4)

# read data    
tif <- raster("//snow/projects/RaBET/RaBET_species/RESULTS/WGEW/speciesspecific_RF_082123_WGEW_2018.tif")
CHM <- raster("//snow/projects/RaBET/RaBET_species/RESULTS/WGEW/wger_chm_mask.tif")
pasture_wgew <- shapefile("//gaea/projects/RaBET/RaBET_species/landuse/WGEW_studyarea_08102023.shp")  %>% st_as_sf()
pasture_wgew_df <- pasture_wgew %>% as.data.frame() %>% rownames_to_column() %>% rename(FID = 1)%>%mutate_at('FID', as.integer)

plot(pasture_wgew)

# merge on common variable, here called 'key'
merge <- merge(shp, csv, by='pasture')

# perhaps save as shapefile again
#e <- extract(tif, shp)


# create classification matrix
reclass_df <- c(864, 950, 1,
                950, 1000, 2,
                1000, 1100, 3,
                1100,1200, 4,
                1200,1300,5,
                1300,1400,6,
                1400,Inf,7)


# reshape the object into a matrix with columns and rows
reclass_m <- matrix(reclass_df,
                    ncol = 3,
                    byrow = TRUE)
reclass_m


# reclassify the raster using the reclass object - reclass_m
dtm_classified <- reclassify(DTM,
                             reclass_m)


p<- as.polygons(dtm_classified == 4)
writeVector(p, "//gaea/projects/RaBET/RaBET_species/grazing/Temp/highelevation_srer.gpkg", filetype="ESRI Shapefile",
            overwrite=TRUE, options="ENCODING=UTF-8")






##Species frequencies###
#extract pixel frequencies
freq_species <- exact_extract(tif, pasture_wgew, function(value, coverage_fraction) {
  data.frame(value = value,
             frac = sum(coverage_fraction)) %>%
    group_by(value) %>%
    summarize(freq = sum(frac), .groups = 'drop') %>%
    pivot_wider(names_from = 'value',
                names_prefix = 'freq_',
                values_from = 'freq')
}) %>% 
  mutate(across(starts_with('freq'), replace_na, 0)) 
##rename rows and columns

freq_species_fin <- freq_species %>% rownames_to_column() %>% rename(FID = rowname, 
                                                                          bareground=freq_1, 
                                                                          grass = freq_2, 
                                                                          mesquite = freq_3,
                                                                          cactus = freq_4,
                                                                          creosote = freq_5,
                                                                          whitethorn = freq_6,
                                                                          urban = freq_7,
                                                                          noclass = freq_NaN) %>% mutate_at('FID', as.integer) %>% 
  mutate(total = rowSums(.[2:8]))%>%
  transform(bareground = (bareground/total)*100,
            grass = (grass/total)*100,
            mesquite = (mesquite/total)*100,
            cactus = (cactus/total)*100,
            creosote = (creosote/total)*100,
            whitethorn = (whitethorn/total)*100,
            urban = (urban/total)*100) %>%
  left_join(pasture_wgew_df, by = "FID") %>% mutate(area = st_area(geometry)/10000) %>%as.data.frame()%>% select(-geometry)


woody_freq_species_fin = freq_species_fin %>% select(-noclass,-total,-bareground,-grass)
baregrass_freq_species_fin = freq_species_fin %>% select(-noclass,-total,-mesquite,-whitethorn,-cactus,-creosote)

write.csv(freq_species_fin, file = 'wgew_species_freq_082223.csv')






getwd()





##CHM and species###
urban <- tif
urban[urban != 7] = NA
urban_project<- projectRaster(urban,CHM)
urban_chm <- mask(CHM, urban_project)
plot(urban_chm)
plot(urban_project)
writeRaster(urban_chm, file = "urban_chm_SRER.tif")

whitethorn <- tif
whitethorn[whitethorn != 6] = NA
whitethorn_project<- projectRaster(whitethorn,CHM)
whitethorn_chm <- mask(CHM, whitethorn_project)
writeRaster(whitethorn_chm, file = "whitethorn_chm_WGEW.tif")


creosote <- tif
creosote[creosote != 5] = NA
creosote_project<- projectRaster(creosote,CHM)
creosote_chm <- mask(CHM, creosote_project)
writeRaster(creosote_chm, file = "creosote_chm_WGEW.tif")

cactus <- tif
cactus[cactus != 4] = NA
cactus_project<- projectRaster(cactus,CHM)
cactus_chm <- mask(CHM, cactus_project)
writeRaster(cactus_chm, file = "cactus_chm_WGEW.tif")

mesquite <- tif
mesquite[mesquite != 3] = NA
mesquite_project<- projectRaster(mesquite,CHM)
mesquite_chm <- mask(CHM, mesquite_project)
writeRaster(mesquite_chm, file = "mesquite_chm_WGEW.tif")

# extract the data
chm_creosote_mean <- exact_extract(creosote_chm, pasture_wgew, 'mean') %>% 
  as.data.frame() %>%  
  rownames_to_column() %>% 
  rename(FID = 1, creosote = 2) %>% 
  transform(FID = as.numeric(FID))
chm_whitethorn_mean <- exact_extract(whitethorn_chm, pasture_wgew, 'mean') %>% 
  as.data.frame()  %>% 
  rownames_to_column() %>% 
  rename(FID = 1, whitethorn = 2)%>% 
  transform(FID = as.numeric(FID))
chm_whitethorn_mean <- exact_extract(whitethorn_chm, pasture_wgew, 'mean') %>% 
  as.data.frame()  %>% 
  rownames_to_column() %>% 
  rename(FID = 1, whitethorn = 2)%>% 
  transform(FID = as.numeric(FID))
chm_cactus_mean <- exact_extract(cactus_chm, pasture_wgew, 'mean') %>% 
  as.data.frame()  %>% 
  rownames_to_column() %>% 
  rename(FID = 1, cactus = 2)%>% 
  transform(FID = as.numeric(FID))
chm_mesquite_mean <- exact_extract(mesquite_chm, pasture_wgew, 'mean') %>%
  as.data.frame()  %>% 
  rownames_to_column() %>% 
  rename(FID = 1, mesquite = 2)%>% 
  transform(FID = as.numeric(FID))


pasture_wgew_df$FID <- as.numeric(pasture_wgew_df$FID)
chm_species_mean_final <- chm_mesquite_mean %>% 
  left_join(chm_cactus_mean, by = "FID") %>%
  left_join(chm_whitethorn_mean, by = "FID") %>%
  left_join(chm_creosote_mean, by = "FID") %>% 
  left_join(pasture_wgew_df, by = "FID")
  
write.csv(chm_species_mean_final, "chm_species_mean_final_wgew.csv")






e <- extent(pasture_wgew)



mesquite_chm = raster('//gaea/projects/RaBET/RaBET_species/landuse/mesquite_chm_wgew.tif')
mesquite_chm_crop <- crop(mesquite_chm, e) 
#mask based on study area extent
mesquite_chm_mask <- mask(mesquite_chm_crop, pasture_wgew)

cactus_chm = raster('//gaea/projects/RaBET/RaBET_species/landuse/cactus_chm_wgew.tif')
cactus_chm_crop <- crop(cactus_chm, e) 
#mask based on study area extent
cactus_chm_mask <- mask(cactus_chm_crop, pasture_wgew)

creosote_chm = raster('//gaea/projects/RaBET/RaBET_species/landuse/creosote_chm_wgew.tif')
creosote_chm_crop <- crop(creosote_chm, e) 
#mask based on study area extent
creosote_chm_mask <- mask(creosote_chm_crop, pasture_wgew)

whitethorn_chm = raster('//gaea/projects/RaBET/RaBET_species/landuse/whitethorn_chm_wgew.tif')
whitethorn_chm_crop <- crop(whitethorn_chm, e) 
#mask based on study area extent
whitethorn_chm_mask <- mask(whitethorn_chm_crop, pasture_wgew)


writeRaster(mesquite_chm_mask, file = "mesquite_chm_mask_wgew_082223.tif")
writeRaster(cactus_chm_mask, file = "cactus_chm_mask_wgew_082223.tif")
writeRaster(creosote_chm_mask, file = "creosote_chm_mask_wgew_082223.tif")
writeRaster(whitethorn_chm_mask, file = "whitethorn_chm_mask_wgew_082223.tif")


###ITD and Crown area#####
# Smoothing CHM

ft_win <- function(x){x * 0.06 + 0.5}
# Setting the fws:
#fws<-7 # dimention 5x5
# Setting the specified height above ground for detectionbreak
#minht<-1.0
# Getting the individual tree detection list
schm_mesquite<-CHMsmoothing(mesquite_chm_mask, "mean", 5)
schm_cactus<-CHMsmoothing(cactus_chm_mask, "mean", 5)
schm_creosote<-CHMsmoothing(creosote_chm_mask, "mean", 5)
schm_whitethorn<-CHMsmoothing(whitethorn_chm_mask, "mean", 5)

mesquitetops_ftools <- vwf(CHM = schm_mesquite, winFun = ft_win, minHeight = 1)
#mesquitetops_rlidr<- FindTreesCHM(schm_mesquite, fws, minht)
mesquitetops_lidr <- locate_trees(schm_mesquite,lmf(5, hmin = 0, shape = "circular"))
mesquitetops_lidr = sf::st_zm(mesquitetops_lidr)
st_write(mesquitetops_lidr, "mesquite_treetops_lidr_5lmf_wgew_082223.shp", overwrite=TRUE)

cactustops_lidr <- locate_trees(schm_cactus,lmf(5, hmin = 0, shape = "circular"))
cactustops_lidr = sf::st_zm(cactustops_lidr)
st_write(cactustops_lidr, "cactus_treetops_lidr_5lmf_wgew_082223.shp", overwrite=TRUE)

creosotetops_lidr <- locate_trees(schm_creosote,lmf(5, hmin = 0, shape = "circular"))
creosotetops_lidr = sf::st_zm(creosotetops_lidr)
st_write(creosotetops_lidr, "creosote_treetops_lidr_5lmf_wgew_082223.shp", overwrite=TRUE)

whitethorntops_lidr <- locate_trees(schm_whitethorn,lmf(5, hmin = 0, shape = "circular"))
whitethorntops_lidr = sf::st_zm(whitethorntops_lidr)
st_write(whitethorntops_lidr, "whitethorn_treetops_lidr_5lmf_wgew_082223.shp", overwrite=TRUE)


#mesquitetops_ftools = mesquitetops_ftools %>% st_set_crs(st_crs(pasture_wgew))
# Plotting the individual tree location on the CHM

#XY<-SpatialPoints(mesquitetops_rlidr[,1:2]) %>% st_as_sf() # Spatial points

mesquitetops_lidr$geometry[1:2,]
#st_write(XY, "mesquite_treetops_rlidar_7fws.shp", overwrite=TRUE,append=FALSE)
#st_write(mesquitetops_ftools, "mesquite_treetops_ftools_5vwf.shp", overwrite=TRUE)



mesquite_lidR_ttops <- shapefile("//gaea/projects/RaBET/RaBET_species/landuse/mesquite_treetops_lidr_5lmf_wgew_082223.shp") %>% st_as_sf() %>% st_set_crs(st_crs(pasture_wgew))
pasture_wgew_area = pasture_wgew %>% mutate(area = st_area(geometry)/10000) %>%as.data.frame()%>% select(-geometry) %>% rownames_to_column() %>% rename(FID = 1)%>%mutate_at('FID', as.integer)
mesquite_wgew_points <- st_join(pasture_wgew,mesquite_lidR_ttops) 
mesquite_wgew_points_count <- mesquite_wgew_points %>% as.data.frame() %>% rownames_to_column() %>% rename(FID = 1)%>% mutate_at('FID', as.integer)%>% group_by(CALCACRES,MUSYM) %>% count() 
mesquite_wgew_points_density <- mesquite_wgew_points_count %>% left_join(pasture_wgew_area, by = c("CALCACRES","MUSYM")) %>% as.data.frame() %>%
  mutate(density = as.numeric(n/area)) %>% mutate(species = "mesquite")


cactus_lidR_ttops <- shapefile("//gaea/projects/RaBET/RaBET_species/landuse/cactus_treetops_lidr_5lmf_wgew_082223.shp") %>% st_as_sf() %>% st_set_crs(st_crs(pasture_wgew))
pasture_wgew_area = pasture_wgew %>% mutate(area = st_area(geometry)/10000) %>%as.data.frame()%>% select(-geometry) %>% rownames_to_column() %>% rename(FID = 1)%>%mutate_at('FID', as.integer)
cactus_wgew_points <- st_join(pasture_wgew,cactus_lidR_ttops) 
cactus_wgew_points_count <- cactus_wgew_points %>% as.data.frame() %>% rownames_to_column() %>% rename(FID = 1)%>% mutate_at('FID', as.integer)%>% group_by(CALCACRES,MUSYM) %>% count() 
cactus_wgew_points_density <- cactus_wgew_points_count %>% left_join(pasture_wgew_area, by = c("CALCACRES","MUSYM")) %>% as.data.frame() %>%
  mutate(density = as.numeric(n/area)) %>% mutate(species = "cactus")


creosote_lidR_ttops <- shapefile("//gaea/projects/RaBET/RaBET_species/landuse/creosote_treetops_lidr_5lmf_wgew_082223.shp") %>% st_as_sf() %>% st_set_crs(st_crs(pasture_wgew))
pasture_wgew_area = pasture_wgew %>% mutate(area = st_area(geometry)/10000) %>%as.data.frame()%>% select(-geometry) %>% rownames_to_column() %>% rename(FID = 1)%>%mutate_at('FID', as.integer)
creosote_wgew_points <- st_join(pasture_wgew,creosote_lidR_ttops) 
creosote_wgew_points_count <- creosote_wgew_points %>% as.data.frame() %>% rownames_to_column() %>% rename(FID = 1)%>% mutate_at('FID', as.integer)%>% group_by(CALCACRES,MUSYM) %>% count() 
creosote_wgew_points_density <- creosote_wgew_points_count %>% left_join(pasture_wgew_area, by = c("CALCACRES","MUSYM")) %>% as.data.frame() %>%
  mutate(density = as.numeric(n/area)) %>% mutate(species = "creosote")

whitethorn_lidR_ttops <- shapefile("//gaea/projects/RaBET/RaBET_species/landuse/whitethorn_treetops_lidr_5lmf_wgew_082223.shp") %>% st_as_sf() %>% st_set_crs(st_crs(pasture_wgew))
pasture_wgew_area = pasture_wgew %>% mutate(area = st_area(geometry)/10000) %>%as.data.frame()%>% select(-geometry) %>% rownames_to_column() %>% rename(FID = 1)%>%mutate_at('FID', as.integer)
whitethorn_wgew_points <- st_join(pasture_wgew,whitethorn_lidR_ttops) 
whitethorn_wgew_points_count <- whitethorn_wgew_points %>% as.data.frame() %>% rownames_to_column() %>% rename(FID = 1)%>% mutate_at('FID', as.integer)%>% group_by(CALCACRES,MUSYM) %>% count() 
whitethorn_wgew_points_density <- whitethorn_wgew_points_count %>% left_join(pasture_wgew_area, by = c("CALCACRES","MUSYM")) %>% as.data.frame() %>%
  mutate(density = as.numeric(n/area)) %>% mutate(species = "whitethorn")


densities_final = mesquite_wgew_points_density %>% 
  rbind(cactus_wgew_points_density)%>%
  rbind(creosote_wgew_points_density)%>%
  rbind(whitethorn_wgew_points_density) 
  
write.csv(densities_final, "densities_final_wgew_082223.csv")

ggplot(densities_final, aes(fill=Grazing, y=density, x=species)) + 
  geom_bar(position="dodge", stat="identity")+
  scale_fill_manual(values=c('dark green', 'dark red'))+
  facet_wrap(~Grazing) 




















####crown area assessment###
mesquite_lidR_ttops <- vect("//gaea/projects/RaBET/RaBET_landuse/landuse/mesquite_treetops_lidr_5lmf_wgew_082223.shp")%>% sf::st_as_sf()
mesquite_chm_mask <- rast("//gaea/projects/RaBET/RaBET_landuse/landuse/mesquite_chm_mask_wgew.tif")
mesquite_lidR_crown <- mcws(treetops = mesquite_lidR_ttops, CHM = mesquite_chm_mask, format = "polygons", minHeight = 1) 
mesquite_lidR_crown = mesquite_lidR_crown%>% st_as_sf() %>% st_set_crs(st_crs(pasture_wgew))
plot(mesquite_lidR_crown)

creosote_lidR_ttops <- shapefile("//gaea/projects/RaBET/RaBET_species/landuse/creosote_treetops_lidr_5lmf_wgew_082223.shp")
creosote_chm_mask <- rast("//gaea/projects/RaBET/RaBET_species/landuse/creosote_chm_mask_wgew.tif")
creosote_lidR_crown <- mcws(treetops = creosote_lidR_ttops, CHM = creosote_chm_mask, format = "polygons", minHeight = 1) 
creosote_lidR_crown = creosote_lidR_crown%>% st_as_sf() %>% st_set_crs(st_crs(pasture_wgew))

cactus_lidR_ttops <- shapefile("//gaea/projects/RaBET/RaBET_species/landuse/cactus_treetops_lidr_5lmf_wgew_082223.shp")
cactus_chm_mask <- raster("//gaea/projects/RaBET/RaBET_species/landuse/cactus_chm_mask_wgew.tif")
cactus_lidR_crown <- mcws(treetops = cactus_lidR_ttops, CHM = cactus_chm_mask, format = "polygons", minHeight = 1, verbose = FALSE) 
cactus_lidR_crown = cactus_lidR_crown%>% st_as_sf() %>% st_set_crs(st_crs(pasture_wgew))

whitethorn_lidR_ttops <- shapefile("//gaea/projects/RaBET/RaBET_species/landuse/whitethorn_treetops_lidr_5lmf_wgew_082223.shp")
whitethorn_chm_mask <- raster("//gaea/projects/RaBET/RaBET_species/landuse/whitethorn_chm_mask_wgew.tif")
whitethorn_lidR_crown <- mcws(treetops = whitethorn_lidR_ttops, CHM = whitethorn_chm_mask, format = "polygons", minHeight = 1, verbose = FALSE) 
whitethorn_lidR_crown = whitethorn_lidR_crown%>% st_as_sf() %>% st_set_crs(st_crs(pasture_wgew))

st_write(mesquite_lidR_crown, "wgew_mesquite_082223_crown.shp")
st_write(creosote_lidR_crown, "wgew_creosote_082223_crown.shp")
st_write(cactus_lidR_crown, "wgew_cactus_082223_crown.shp")
st_write(whitethorn_lidR_crown, "wgew_whitethorn_082223_crown.shp")

wgew_mesquite_lidR_crown <- shapefile("//gaea/projects/RaBET/RaBET_species/landuse/wgew_mesquite_082223_crown.shp") %>% st_as_sf()
wgew_creosote_lidR_crown <- shapefile("//gaea/projects/RaBET/RaBET_species/landuse/wgew_creosote_082223_crown.shp")%>% st_as_sf()
wgew_cactus_lidR_crown <- shapefile("//gaea/projects/RaBET/RaBET_species/landuse/wgew_cactus_082223_crown.shp")%>% st_as_sf()
wgew_whitethorn_lidR_crown <- shapefile("//gaea/projects/RaBET/RaBET_species/landuse/wgew_whitethorn_082223_crown.shp")%>% st_as_sf()
windows()
plot(wgew_creosote_lidR_crown)
















###MASK PARALLEL
#mask training rasters with rasterized training Points
ptm <- proc.time()
# Start by defining the function that does the work 
join_poly <- function(crown, pasture) {
  joined_poly <- sf::st_join(crown, pasture)
  return(joined_poly)
}

# We need the RasterStack to be a list of RasterLayers
layer_list <- list(pasture_wgew)


# Set up parallel processing
num_cores = 4

# Start the cluster
cl <- makeCluster(num_cores)

# Load required libraries
clusterEvalQ(cl, {
  library(sf)
  library(raster)
})

# Run the function in parallel
masked_list <- parallel::parLapply(cl = cl,
                                   X = layer_list,
                                   fun = join_poly,
                                   crown = wgew_whitethorn_lidR_crown) 


gc()

#NEON_indices_mask <- raster::stack(masked_list)
#writeRaster(x = NEON_indices_mask, filename = "neon_stack.grd", overwrite=TRUE)
#neon_mask <- raster(x = "neon_stack.grd")

masked_list_df <-  masked_list %>% as.data.frame() %>% rename(FID =1)


stopCluster(cl)




whitethorn_crownArea_df <- masked_list_df %>% rename(FID =1)%>%
  group_by(FID,MUSYM, CALCACRES) %>%
  summarise_at(vars(crownArea,Z), list(mean = mean)) %>% 
  as.data.frame() %>%
  mutate(vol_mean = Z_mean*crownArea_mean)%>%
  left_join(select(pasture_wgew_df, c('FID','geometry')),by="FID")%>% 
  mutate(area = st_area(geometry)/10000) %>%as.data.frame()%>% select(-geometry)



write.csv(whitethorn_crownArea_df, file = "whitethorn_crownArea_df_wgew_082223.csv")










ggplot2(mesquite_crownArea_df, aes(fill=Grazing, y=crownArea_mean, x=Grazing)) + 
  geom_bar(position="dodge", stat="identity")+
  scale_fill_manual(values=c('dark green', 'dark red'))























































































































































































#########################SRER####################################
tif <- raster("//snow/projects/RaBET/RaBET_species/RESULTS/final_results/Classification/RF_allyr_SRER.tif")
CHM <- raster("//snow/projects/RaBET/RaBET_species/RESULTS/SRER/NEON_CHM_2018_mask_SRER.tif")
pasture_srer <- shapefile("//gaea/projects/RaBET/RaBET_species/landuse/srer_grid_landuse_100m_inter_soil.shp") %>% st_as_sf()
 #   mutate(area = as.numeric(st_area(geometry)/10000)) %>%filter(area >=1)
pasture_srer_df <- pasture_srer %>% as.data.frame() %>% rownames_to_column() %>% rename(FID = rowname)%>% mutate_at('FID', as.integer)

# create classification matrix
reclass_df <- c(864, 950, 1,
                950, 1000, 2,
                1000, 1100, 3,
                1100,1200, 4,
                1200,1300,5,
                1300,1400,6,
                1400,Inf,7)


# reshape the object into a matrix with columns and rows
reclass_m <- matrix(reclass_df,
                    ncol = 3,
                    byrow = TRUE)
reclass_m


# reclassify the raster using the reclass object - reclass_m
dtm_classified <- reclassify(DTM,
                             reclass_m)


p<- as.polygons(dtm_classified == 4)
writeVector(p, "//gaea/projects/RaBET/RaBET_species/grazing/Temp/highelevation_srer.gpkg", filetype="ESRI Shapefile",
            overwrite=TRUE, options="ENCODING=UTF-8")






##Species frequencies###
#extract pixel frequencies
freq_species <- exact_extract(tif, pasture_srer, function(value, coverage_fraction) {
  data.frame(value = value,
             frac = sum(coverage_fraction)) %>%
    group_by(value) %>%
    summarize(freq = sum(frac), .groups = 'drop') %>%
    pivot_wider(names_from = 'value',
                names_prefix = 'freq_',
                values_from = 'freq')
}) %>% 
  mutate(across(starts_with('freq'), replace_na, 0)) 
##rename rows and columns

freq_species_fin <- freq_species %>% rownames_to_column() %>% rename(FID = rowname, 
                                                                     bareground=freq_1, 
                                                                     grass = freq_2, 
                                                                     mesquite = freq_3,
                                                                     cactus = freq_4,
                                                                     lotebush = freq_5,
                                                                     paloverde = freq_6,
                                                                     creosote = freq_7,
                                                                     noclass = freq_NA) %>% mutate_at('FID', as.integer) %>% 
  mutate(total = rowSums(.[2:8]))%>%
  transform(bareground = (bareground/total)*100,
            grass = (grass/total)*100,
            mesquite = (mesquite/total)*100,
            cactus = (cactus/total)*100,
            creosote = (creosote/total)*100,
            lotebush = (lotebush/total)*100,
            paloverde = (paloverde/total)*100) %>%
  left_join(pasture_srer_df, by = "FID")%>%
  mutate(grazing = ifelse(AUY_per_ha > 300, '1','0'))%>% mutate(area = st_area(geometry)/10000) %>%as.data.frame()%>% select(-geometry)
  

#with or without grazing

woody_freq_species_fin = freq_species_fin %>% select(-noclass,-total,-bareground,-grass)
baregrass_freq_species_fin = freq_species_fin %>% select(-noclass,-total,-mesquite,-lotebush,-cactus,-creosote, -paloverde)
write.csv(freq_species_fin, file = 'srer_species_freq_082223.csv')














getwd()





##CHM and species###
creosote <- tif
creosote[creosote != 7] = NA
creosote_project<- projectRaster(creosote,CHM)
creosote_chm <- mask(CHM, creosote_project)
plot(creosote_chm)
plot(creosote_project)
writeRaster(creosote_chm, file = "creosote_chm_SRER.tif")

paloverde <- tif
paloverde[paloverde != 6] = NA
paloverde_project<- projectRaster(paloverde,CHM)
paloverde_chm <- mask(CHM, paloverde_project)
writeRaster(paloverde_chm, file = "paloverde_chm_srer.tif")


lotebush <- tif
lotebush[lotebush != 5] = NA
lotebush_project<- projectRaster(lotebush,CHM)
lotebush_chm <- mask(CHM, lotebush_project)
writeRaster(lotebush_chm, file = "lotebush_chm_srer.tif")

cactus <- tif
cactus[cactus != 4] = NA
cactus_project<- projectRaster(cactus,CHM)
cactus_chm <- mask(CHM, cactus_project)
writeRaster(cactus_chm, file = "cactus_chm_srer.tif")

mesquite <- tif
mesquite[mesquite != 3] = NA
mesquite_project<- projectRaster(mesquite,CHM)
mesquite_chm <- mask(CHM, mesquite_project)
writeRaster(mesquite_chm, file = "mesquite_chm_srer.tif")

# extract the data
chm_creosote_mean <- exact_extract(creosote_chm, pasture_srer, 'mean') %>% 
  as.data.frame() %>%  
  rownames_to_column() %>% 
  rename(FID = 1, creosote = 2) %>% 
  transform(FID = as.numeric(FID))
chm_paloverde_mean <- exact_extract(paloverde_chm, pasture_srer, 'mean') %>% 
  as.data.frame()  %>% 
  rownames_to_column() %>% 
  rename(FID = 1, paloverde = 2)%>% 
  transform(FID = as.numeric(FID))
chm_lotebush_mean <- exact_extract(lotebush_chm, pasture_srer, 'mean') %>% 
  as.data.frame()  %>% 
  rownames_to_column() %>% 
  rename(FID = 1, lotebush = 2)%>% 
  transform(FID = as.numeric(FID))
chm_cactus_mean <- exact_extract(cactus_chm, pasture_srer, 'mean') %>% 
  as.data.frame()  %>% 
  rownames_to_column() %>% 
  rename(FID = 1, cactus = 2)%>% 
  transform(FID = as.numeric(FID))
chm_mesquite_mean <- exact_extract(mesquite_chm, pasture_srer, 'mean') %>%
  as.data.frame()  %>% 
  rownames_to_column() %>% 
  rename(FID = 1, mesquite = 2)%>% 
  transform(FID = as.numeric(FID))


pasture_srer_df$FID <- as.numeric(pasture_srer_df$FID)
chm_species_mean_final <- chm_mesquite_mean %>% 
  left_join(chm_cactus_mean, by = "FID") %>%
  left_join(chm_paloverde_mean, by = "FID") %>%
  left_join(chm_creosote_mean, by = "FID") %>% 
  left_join(chm_lotebush_mean, by = "FID") %>% 
  left_join(pasture_srer_df, by = "FID")

write.csv(chm_species_mean_final, "chm_species_mean_final_srer.csv")









e <- extent(pasture_srer)

mesquite_chm = raster('//gaea/projects/RaBET/RaBET_species/landuse/datasets/mesquite_chm_srer.tif')
mesquite_chm_crop <- crop(mesquite_chm, e) 
#mask based on study area extent
mesquite_chm_mask <- mask(mesquite_chm_crop, pasture_srer)

cactus_chm = raster('//gaea/projects/RaBET/RaBET_species/landuse/datasets/cactus_chm_srer.tif')
cactus_chm_crop <- crop(cactus_chm, e) 
#mask based on study area extent
cactus_chm_mask <- mask(cactus_chm_crop, pasture_srer)

creosote_chm = raster('//gaea/projects/RaBET/RaBET_species/landuse/datasets/creosote_chm_srer.tif')
creosote_chm_crop <- crop(creosote_chm, e) 
#mask based on study area extent
creosote_chm_mask <- mask(creosote_chm_crop, pasture_srer)

paloverde_chm = raster('//gaea/projects/RaBET/RaBET_species/landuse/datasets/paloverde_chm_srer.tif')
paloverde_chm_crop <- crop(paloverde_chm, e) 
#mask based on study area extent
paloverde_chm_mask <- mask(paloverde_chm_crop, pasture_srer)

lotebush_chm = raster('//gaea/projects/RaBET/RaBET_species/landuse/datasets/lotebush_chm_srer.tif')
lotebush_chm_crop <- crop(lotebush_chm, e) 
#mask based on study area extent
lotebush_chm_mask <- mask(lotebush_chm_crop, pasture_srer)

writeRaster(mesquite_chm_mask, file = "mesquite_chm_mask_srer_082223.tif")
writeRaster(cactus_chm_mask, file = "cactus_chm_mask_srer_082223.tif")
writeRaster(creosote_chm_mask, file = "creosote_chm_mask_srer_082223.tif")
writeRaster(paloverde_chm_mask, file = "paloverde_chm_mask_srer_082223.tif")
writeRaster(lotebush_chm_mask, file = "lotebush_chm_mask_srer_082223.tif")


###ITD and Crown area#####
# Smoothing CHM

ft_win <- function(x){x * 0.06 + 0.5}
# Setting the fws:
#fws<-7 # dimention 5x5
# Setting the specified height above ground for detectionbreak
#minht<-1.0
# Getting the individual tree detection list
schm_mesquite<-CHMsmoothing(mesquite_chm_mask, "mean", 5)
schm_cactus<-CHMsmoothing(cactus_chm_mask, "mean", 5)
schm_creosote<-CHMsmoothing(creosote_chm_mask, "mean", 5)
schm_paloverde<-CHMsmoothing(paloverde_chm_mask, "mean", 5)
schm_lotebush<-CHMsmoothing(lotebush_chm_mask, "mean", 5)

mesquitetops_ftools <- vwf(CHM = schm_mesquite, winFun = ft_win, minHeight = 1)
#mesquitetops_rlidr<- FindTreesCHM(schm_mesquite, fws, minht)
mesquitetops_lidr <- locate_trees(schm_mesquite,lmf(5, hmin = 0, shape = "circular"))
mesquitetops_lidr = sf::st_zm(mesquitetops_lidr)
st_write(mesquitetops_lidr, "mesquite_treetops_lidr_5lmf_srer_082223.shp", overwrite=TRUE)

cactustops_lidr <- locate_trees(schm_cactus,lmf(5, hmin = 0, shape = "circular"))
cactustops_lidr = sf::st_zm(cactustops_lidr)
st_write(cactustops_lidr, "cactus_treetops_lidr_5lmf_srer_082223.shp", overwrite=TRUE)

creosotetops_lidr <- locate_trees(schm_creosote,lmf(5, hmin = 0, shape = "circular"))
creosotetops_lidr = sf::st_zm(creosotetops_lidr)
st_write(creosotetops_lidr, "creosote_treetops_lidr_5lmf_srer_082223.shp", overwrite=TRUE)

paloverdetops_lidr <- locate_trees(schm_paloverde,lmf(5, hmin = 0, shape = "circular"))
paloverdetops_lidr = sf::st_zm(paloverdetops_lidr)
st_write(paloverdetops_lidr, "paloverde_treetops_lidr_5lmf_srer_082223.shp", overwrite=TRUE)

lotebushtops_lidr <- locate_trees(schm_lotebush,lmf(5, hmin = 0, shape = "circular"))
lotebushtops_lidr = sf::st_zm(lotebushtops_lidr)
st_write(lotebushtops_lidr, "lotebush_treetops_lidr_5lmf_srer_082223.shp", overwrite=TRUE)
#mesquitetops_ftools=st_as_sf(mesquitetops_ftools)
plot(XY)
#mesquitetops_ftools = mesquitetops_ftools %>% st_set_crs(st_crs(pasture_srer))
# Plotting the individual tree location on the CHM

#XY<-SpatialPoints(mesquitetops_rlidr[,1:2]) %>% st_as_sf() # Spatial points

mesquitetops_lidr$geometry[1:2,]
#st_write(XY, "mesquite_treetops_rlidar_7fws.shp", overwrite=TRUE,append=FALSE)
#st_write(mesquitetops_ftools, "mesquite_treetops_ftools_5vwf.shp", overwrite=TRUE)

pasture_srer_area = pasture_srer %>% mutate(area = as.numeric(st_area(geometry)/10000)) %>%as.data.frame()%>% select(-geometry)%>%
  mutate(grazing = ifelse(AUY_per_ha > 300, '1','0'))%>%filter(area >=1)

mesquite_lidR_ttops <- shapefile("//gaea/projects/RaBET/RaBET_species/landuse/datasets/mesquite_treetops_lidr_5lmf.shp") %>% 
  st_as_sf() %>% 
  st_set_crs(st_crs(pasture_srer))
mesquite_srer_points <- st_join(pasture_srer,mesquite_lidR_ttops) %>%
  mutate(grazing = as.numeric(ifelse(AUY_per_ha > 300, '1','0')))
mesquite_srer_points_count <- mesquite_srer_points %>% as.data.frame() %>% 
  group_by(layer, MUSYM,grazing,FID_srer_g) %>% 
  count() %>%
  mutate(species = "mesquite")
mesquite_srer_points_density <- mesquite_srer_points_count %>% 
  group_by(layer, MUSYM,grazing) %>%
  summarise_at(vars(n), list(density = mean))%>%
  mutate(species = "mesquite")

#plot(mesquite_lidR_ttops)

cactus_lidR_ttops <- shapefile("//gaea/projects/RaBET/RaBET_species/landuse/datasets/cactus_treetops_lidr_5lmf.shp") %>% st_as_sf() %>% st_set_crs(st_crs(pasture_srer))
cactus_srer_points <- st_join(pasture_srer,cactus_lidR_ttops) %>%
  mutate(grazing = as.numeric(ifelse(AUY_per_ha > 300, '1','0')))
cactus_srer_points_count <- cactus_srer_points %>% as.data.frame() %>% 
  group_by(layer, MUSYM,grazing,FID_srer_g) %>% 
  count() %>%
 mutate(species = "cactus")
cactus_srer_points_density <- cactus_srer_points_count %>% 
  group_by(layer, MUSYM,grazing) %>%
  summarise_at(vars(n), list(density = mean))%>%
  mutate(species = "mesquite")

creosote_lidR_ttops <- shapefile("//gaea/projects/RaBET/RaBET_species/landuse/datasets/creosote_treetops_lidr_5lmf.shp") %>% st_as_sf() %>% st_set_crs(st_crs(pasture_srer))
creosote_srer_points <- st_join(pasture_srer,creosote_lidR_ttops) %>%
  mutate(grazing = as.numeric(ifelse(AUY_per_ha > 300, '1','0')))
creosote_srer_points_count <- creosote_srer_points %>% as.data.frame() %>% 
  group_by(layer, MUSYM,grazing,FID_srer_g) %>% 
  count() %>%
mutate(species = "creosote")
creosote_srer_points_density <- creosote_srer_points_count %>% 
  group_by(layer, MUSYM,grazing) %>%
  summarise_at(vars(n), list(density = mean))%>%
  mutate(species = "mesquite")


paloverde_lidR_ttops <- shapefile("//gaea/projects/RaBET/RaBET_species/landuse/datasets/paloverde_treetops_lidr_5lmf.shp") %>% st_as_sf() %>% st_set_crs(st_crs(pasture_srer))
paloverde_srer_points <- st_join(pasture_srer,paloverde_lidR_ttops) %>%
  mutate(grazing = as.numeric(ifelse(AUY_per_ha > 300, '1','0')))
paloverde_srer_points_count <- paloverde_srer_points %>% as.data.frame() %>% 
  group_by(layer, MUSYM,grazing,FID_srer_g) %>% 
  count() %>%
mutate(species = "paloverde")
paloverde_srer_points_density <- paloverde_srer_points_count %>% 
  group_by(layer, MUSYM,grazing) %>%
  summarise_at(vars(n), list(density = mean))%>%
  mutate(species = "mesquite")

lotebush_lidR_ttops <- shapefile("//gaea/projects/RaBET/RaBET_species/landuse/datasets/lotebush_treetops_lidr_5lmf.shp") %>% st_as_sf() %>% st_set_crs(st_crs(pasture_srer))
lotebush_srer_points <- st_join(pasture_srer,lotebush_lidR_ttops) %>%
  mutate(grazing = as.numeric(ifelse(AUY_per_ha > 300, '1','0')))
lotebush_srer_points_count <- lotebush_srer_points %>% as.data.frame() %>% 
  group_by(layer, MUSYM,grazing,FID_srer_g) %>% 
  count() %>%
  mutate(species = "lotebush")
lotebush_srer_points_density <- lotebush_srer_points_count %>% 
  group_by(layer, MUSYM,grazing) %>%
  summarise_at(vars(n), list(density = mean))%>%
  mutate(species = "mesquite")


densities_final = mesquite_srer_points_density %>% 
  rbind(cactus_srer_points_density)%>%
  rbind(creosote_srer_points_density)%>%
  rbind(paloverde_srer_points_density) %>%
  rbind(lotebush_srer_points_density) 
  
count_final = mesquite_srer_points_count %>% 
  rbind(cactus_srer_points_count)%>%
  rbind(creosote_srer_points_count)%>%
  rbind(paloverde_srer_points_count) %>%
  rbind(lotebush_srer_points_count) 

write.csv(densities_final, "densities_final_srer_083123.csv")
write.csv(count_final, "count_final_srer_083123.csv")


ggplot(densities_final, aes(fill=Grazing, y=density, x=species)) + 
  geom_bar(position="dodge", stat="identity")+
  scale_fill_manual(values=c('dark green', 'dark red'))+
  facet_wrap(~layer) 




















####crown area assessment###
srer_mesquite_lidR_ttops <- shapefile("//gaea/projects/RaBET/RaBET_species/landuse/datasets/mesquite_treetops_lidr_5lmf.shp")
srer_mesquite_chm_mask <- raster("//gaea/projects/RaBET/RaBET_species/landuse/datasets/mesquite_chm_mask_SRER.tif")
srer_mesquite_lidR_crown <- mcws(treetops = srer_mesquite_lidR_ttops, CHM = srer_mesquite_chm_mask, format = "polygons", minHeight = 1, verbose = FALSE) 
srer_mesquite_lidR_crown = srer_mesquite_lidR_crown%>% st_as_sf() %>% st_set_crs(st_crs(pasture_srer))


srer_creosote_lidR_ttops <- shapefile("//gaea/projects/RaBET/RaBET_species/landuse/datasets/creosote_treetops_lidr_5lmf.shp")
srer_creosote_chm_mask <- raster("//gaea/projects/RaBET/RaBET_species/landuse/datasets/creosote_chm_mask_SRER.tif")
srer_creosote_lidR_crown <- mcws(treetops = srer_creosote_lidR_ttops, CHM = srer_creosote_chm_mask, format = "polygons", minHeight = 1, verbose = FALSE) 
srer_creosote_lidR_crown = srer_creosote_lidR_crown%>% st_as_sf() %>% st_set_crs(st_crs(pasture_srer))

srer_cactus_lidR_ttops <- shapefile("//gaea/projects/RaBET/RaBET_species/landuse/datasets/cactus_treetops_lidr_5lmf.shp")
srer_cactus_chm_mask <- raster("//gaea/projects/RaBET/RaBET_species/landuse/datasets/cactus_chm_mask_SRER.tif")
srer_cactus_lidR_crown <- mcws(treetops = srer_cactus_lidR_ttops, CHM = srer_cactus_chm_mask, format = "polygons", minHeight = 1, verbose = FALSE) 
srer_cactus_lidR_crown = srer_cactus_lidR_crown%>% st_as_sf() %>% st_set_crs(st_crs(pasture_srer))

srer_paloverde_lidR_ttops <- shapefile("//gaea/projects/RaBET/RaBET_species/landuse/datasets/paloverde_treetops_lidr_5lmf.shp")
srer_paloverde_chm_mask <- raster("//gaea/projects/RaBET/RaBET_species/landuse/datasets/paloverde_chm_mask_SRER.tif")
srer_paloverde_lidR_crown <- mcws(treetops = srer_paloverde_lidR_ttops, CHM = srer_paloverde_chm_mask, format = "polygons", minHeight = 1, verbose = FALSE) 
srer_paloverde_lidR_crown = srer_paloverde_lidR_crown%>% st_as_sf() %>% st_set_crs(st_crs(pasture_srer))

srer_lotebush_lidR_ttops <- shapefile("//gaea/projects/RaBET/RaBET_species/landuse/datasets/lotebush_treetops_lidr_5lmf.shp")
srer_lotebush_chm_mask <- raster("//gaea/projects/RaBET/RaBET_species/landuse/datasets/lotebush_chm_mask_SRER.tif")
srer_lotebush_lidR_crown <- mcws(treetops = srer_lotebush_lidR_ttops, CHM = srer_lotebush_chm_mask, format = "polygons", minHeight = 1, verbose = FALSE) 
srer_lotebush_lidR_crown = srer_lotebush_lidR_crown%>% st_as_sf() %>% st_set_crs(st_crs(pasture_srer))

st_write(srer_mesquite_lidR_crown, "srer_mesquite_crown_082223.shp")
st_write(srer_creosote_lidR_crown, "srer_creosote_crown_082223.shp")
st_write(srer_cactus_lidR_crown, "srer_cactus_crown_082223.shp")
st_write(srer_paloverde_lidR_crown, "srer_paloverde_crown_082223.shp")
st_write(srer_lotebush_lidR_crown, "srer_lotebush_crown_082223.shp")

srer_mesquite_lidR_crown <- shapefile("//gaea/projects/RaBET/RaBET_species/landuse/srer_mesquite_crown_082223.shp") %>% st_as_sf()
srer_creosote_lidR_crown <- shapefile("//gaea/projects/RaBET/RaBET_species/landuse/srer_creosote_crown_082223.shp")%>% st_as_sf()
srer_cactus_lidR_crown <- shapefile("//gaea/projects/RaBET/RaBET_species/landuse/srer_cactus_crown_082223.shp")%>% st_as_sf()
srer_paloverde_lidR_crown <- shapefile("//gaea/projects/RaBET/RaBET_species/landuse/srer_paloverde_crown_082223.shp")%>% st_as_sf()
srer_lotebush_lidR_crown <- shapefile("//gaea/projects/RaBET/RaBET_species/landuse/srer_lotebush_crown_082223.shp")%>% st_as_sf()
windows()
plot(srer_mesquite_lidR_crown)

###MASK PARALLEL
#mask training rasters with rasterized training Points
ptm <- proc.time()
# Start by defining the function that does the work 
join_poly <- function(crown, pasture) {
  joined_poly <- sf::st_join(crown, pasture)
  return(joined_poly)
}

# We need the RasterStack to be a list of RasterLayers
layer_list <- list(pasture_srer)


# Set up parallel processing
num_cores = 4

# Start the cluster
cl <- makeCluster(num_cores)

# Load required libraries
clusterEvalQ(cl, {
  library(sf)
  library(raster)
})

# Run the function in parallel
masked_list <- parallel::parLapply(cl = cl,
                                   X = layer_list,
                                   fun = join_poly,
                                   crown = srer_paloverde_lidR_crown) 


gc()

#NEON_indices_mask <- raster::stack(masked_list)
#writeRaster(x = NEON_indices_mask, filename = "neon_stack.grd", overwrite=TRUE)
#neon_mask <- raster(x = "neon_stack.grd")

masked_list_df <-  masked_list %>% as.data.frame()


stopCluster(cl)


masked_list_df2 <- masked_list_df

paloverde_crownArea_df <- masked_list_df2 %>%
  as.data.frame() %>%
  mutate(grazing = ifelse(AUY_per_ha > 300, "2",'1'))%>%
  select(-geometry)
  
#left_join(masked_list_df, c('FID_srer_g','layer','MUSYM'),by="FID_srer_g")%>%
#  mutate(area = st_area(geometry)/10000)%>%as.data.frame()%>% select(-geometry)%>%
#  mutate(species = "paloverde")


write.csv(paloverde_crownArea_df, file = "paloverde_crownArea_srer_df_090123.csv")



