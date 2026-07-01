# bry gbif

library(ape)
library(dplyr)
library(phangorn)
library(raster)
library(sf)
library(sp)
library(phytools)
library(purrr)
#install.packages("sf")
#install.packages("sp")
#install.packages("phangorn")
#install.packages("purrr")

setwd("")

# Read tree and coordinates data
bry_tree = force.ultrametric(read.tree("A.tre"))
bry_data = read.csv("distribution_data.csv") # The curated GBIF coordinates

# Write list of species sampled to a file for easy reference later
bry_taxa = as.character(bry_tree$tip.label)
write(bry_taxa,"bry_tree_taxon_list.txt")

# Write species sampled in GBIF for later use
unique_bry_species = unique(bry_data$species)
write(unique_bry_species,"bry_gbif_unique_species.txt")

# Find spaces sampled by both the tree and GBIF
bry_tree_taxa_in_gbif = data.frame(species = intersect(bry_data$species, bry_tree$tip.label))
bry_gbif = semi_join(bry_data,bry_tree_taxa_in_gbif,by = "species")

# Prune tree to just taxa overlapping with GBIF
bry_gbif_tree = keep.tip(bry_tree,unique(bry_gbif$species))

# Now keep only one occurrence per species per 1km grid
bry_sp <- SpatialPointsDataFrame(coords = bry_gbif[, c("decimalLongitude", "decimalLatitude")], 
                                    data = bry_gbif)
grid <- raster(extent(bry_sp), resolution = 1) # 1km x 1km resolution
bry_gbif$grid_cell <- cellFromXY(grid, coordinates(bry_sp))
sampled_data = bry_gbif %>%
  group_by(species, grid_cell) %>%
  sample_n(1)
write.csv(sampled_data, "1_per_grid_bry_gbif.csv", row.names = FALSE)

# Prepare GBIF data for extracting climatic variables
sp = sampled_data[,c("species","decimalLatitude","decimalLongitude")]
xy <- sp[c("decimalLongitude","decimalLatitude")]
sp_xy <- sp[c("species", "decimalLongitude","decimalLatitude")]

###### Below here: need to replace these with extracting the new CHELSA stuff
bio1 = raster("tiff/wc2.1_30s_bio_1.tif")
xy <- sp[c("decimalLongitude","decimalLatitude")]
sp_xy <- sp[c("species", "decimalLongitude","decimalLatitude")]
xy_bio <- raster::extract(bio1, xy)
bio1 <- cbind(sp_xy, xy_bio)
write.csv(bio1,"var/bry_bio1.csv")

bio2 = raster("tiff/wc2.1_30s_bio_2.tif")
xy <- sp[c("decimalLongitude","decimalLatitude")]
sp_xy <- sp[c("species", "decimalLongitude","decimalLatitude")]
xy_bio <- raster::extract(bio2, xy)
bio2 <- cbind(sp_xy, xy_bio)
write.csv(bio2,"var/bry_bio2.csv")

bio3 = raster("tiff/wc2.1_30s_bio_3.tif")
xy <- sp[c("decimalLongitude","decimalLatitude")]
sp_xy <- sp[c("species", "decimalLongitude","decimalLatitude")]
xy_bio <- raster::extract(bio3, xy)
bio3 <- cbind(sp_xy, xy_bio)
write.csv(bio3,"var/bry_bio3.csv")

bio4 = raster("tiff/wc2.1_30s_bio_4.tif")
xy <- sp[c("decimalLongitude","decimalLatitude")]
sp_xy <- sp[c("species", "decimalLongitude","decimalLatitude")]
xy_bio <- raster::extract(bio4, xy)
bio4 <- cbind(sp_xy, xy_bio)
write.csv(bio4,"var/bry_bio4.csv")

bio5 = raster("tiff/wc2.1_30s_bio_5.tif")
xy <- sp[c("decimalLongitude","decimalLatitude")]
sp_xy <- sp[c("species", "decimalLongitude","decimalLatitude")]
xy_bio <- raster::extract(bio5, xy)
bio5 <- cbind(sp_xy, xy_bio)
write.csv(bio5,"var/bry_bio5.csv")

bio6 = raster("tiff/wc2.1_30s_bio_6.tif")
xy <- sp[c("decimalLongitude","decimalLatitude")]
sp_xy <- sp[c("species", "decimalLongitude","decimalLatitude")]
xy_bio <- raster::extract(bio6, xy)
bio6 <- cbind(sp_xy, xy_bio)
write.csv(bio6,"var/bry_bio6.csv")

bio7 = raster("tiff/wc2.1_30s_bio_7.tif")
xy <- sp[c("decimalLongitude","decimalLatitude")]
sp_xy <- sp[c("species", "decimalLongitude","decimalLatitude")]
xy_bio <- raster::extract(bio7, xy)
bio7 <- cbind(sp_xy, xy_bio)
write.csv(bio7,"var/bry_bio7.csv")

bio8 = raster("tiff/wc2.1_30s_bio_8.tif")
xy <- sp[c("decimalLongitude","decimalLatitude")]
sp_xy <- sp[c("species", "decimalLongitude","decimalLatitude")]
xy_bio <- raster::extract(bio8, xy)
bio8 <- cbind(sp_xy, xy_bio)
write.csv(bio8,"var/bry_bio8.csv")

bio9 = raster("tiff/wc2.1_30s_bio_9.tif")
xy <- sp[c("decimalLongitude","decimalLatitude")]
sp_xy <- sp[c("species", "decimalLongitude","decimalLatitude")]
xy_bio <- raster::extract(bio9, xy)
bio9 <- cbind(sp_xy, xy_bio)
write.csv(bio9,"var/bry_bio9.csv")

bio10 = raster("tiff/wc2.1_30s_bio_10.tif")
xy <- sp[c("decimalLongitude","decimalLatitude")]
sp_xy <- sp[c("species", "decimalLongitude","decimalLatitude")]
xy_bio <- raster::extract(bio10, xy)
bio10 <- cbind(sp_xy, xy_bio)
write.csv(bio10,"var/bry_bio10.csv")

bio11 = raster("tiff/wc2.1_30s_bio_11.tif")
xy <- sp[c("decimalLongitude","decimalLatitude")]
sp_xy <- sp[c("species", "decimalLongitude","decimalLatitude")]
xy_bio <- raster::extract(bio11, xy)
bio11 <- cbind(sp_xy, xy_bio)
write.csv(bio11,"var/bry_bio11.csv")

bio12 = raster("tiff/wc2.1_30s_bio_12.tif")
xy <- sp[c("decimalLongitude","decimalLatitude")]
sp_xy <- sp[c("species", "decimalLongitude","decimalLatitude")]
xy_bio <- raster::extract(bio12, xy)
bio12 <- cbind(sp_xy, xy_bio)
write.csv(bio12,"var/bry_bio12.csv")

bio13 = raster("tiff/wc2.1_30s_bio_13.tif")
xy <- sp[c("decimalLongitude","decimalLatitude")]
sp_xy <- sp[c("species", "decimalLongitude","decimalLatitude")]
xy_bio <- raster::extract(bio13, xy)
bio13 <- cbind(sp_xy, xy_bio)
write.csv(bio13,"var/bry_bio13.csv")

bio14 = raster("tiff/wc2.1_30s_bio_14.tif")
xy <- sp[c("decimalLongitude","decimalLatitude")]
sp_xy <- sp[c("species", "decimalLongitude","decimalLatitude")]
xy_bio <- raster::extract(bio14, xy)
bio14 <- cbind(sp_xy, xy_bio)
write.csv(bio14,"var/bry_bio14.csv")

bio15 = raster("tiff/wc2.1_30s_bio_15.tif")
xy <- sp[c("decimalLongitude","decimalLatitude")]
sp_xy <- sp[c("species", "decimalLongitude","decimalLatitude")]
xy_bio <- raster::extract(bio15, xy)
bio15 <- cbind(sp_xy, xy_bio)
write.csv(bio15,"var/bry_bio15.csv")

bio16 = raster("tiff/wc2.1_30s_bio_16.tif")
xy <- sp[c("decimalLongitude","decimalLatitude")]
sp_xy <- sp[c("species", "decimalLongitude","decimalLatitude")]
xy_bio <- raster::extract(bio16, xy)
bio16 <- cbind(sp_xy, xy_bio)
write.csv(bio16,"var/bry_bio16.csv")

bio17 = raster("tiff/wc2.1_30s_bio_17.tif")
xy <- sp[c("decimalLongitude","decimalLatitude")]
sp_xy <- sp[c("species", "decimalLongitude","decimalLatitude")]
xy_bio <- raster::extract(bio17, xy)
bio17 <- cbind(sp_xy, xy_bio)
write.csv(bio17,"var/bry_bio17.csv")

bio18 = raster("tiff/wc2.1_30s_bio_18.tif")
xy <- sp[c("decimalLongitude","decimalLatitude")]
sp_xy <- sp[c("species", "decimalLongitude","decimalLatitude")]
xy_bio <- raster::extract(bio18, xy)
bio18 <- cbind(sp_xy, xy_bio)
write.csv(bio18,"var/bry_bio18.csv")

bio19 = raster("tiff/wc2.1_30s_bio_19.tif")
xy <- sp[c("decimalLongitude","decimalLatitude")]
sp_xy <- sp[c("species", "decimalLongitude","decimalLatitude")]
xy_bio <- raster::extract(bio19, xy)
bio19 <- cbind(sp_xy, xy_bio)
write.csv(bio19,"var/bry_bio19.csv")

#Global Aridity Index and Potential Evapotranspiration (ET0) Database: Version 3
ai = raster("tiff/ai_v3_yr.tif") 
xy <- sp[c("decimalLongitude","decimalLatitude")]
sp_xy <- sp[c("species", "decimalLongitude","decimalLatitude")]
xy_bio <- raster::extract(ai, xy)
ai <- cbind(sp_xy, xy_bio)
write.csv(ai,"var/bry_ai.csv")

pet = raster("tiff/et0_v3_yr.tif") #https://figshare.com/articles/dataset/Global_Aridity_Index_and_Potential_Evapotranspiration_ET0_Climate_Database_v2/7504448/5?file=34377269
xy <- sp[c("decimalLongitude","decimalLatitude")]
sp_xy <- sp[c("species", "decimalLongitude","decimalLatitude")]
xy_bio <- raster::extract(pet, xy)
pet <- cbind(sp_xy, xy_bio)
write.csv(pet,"var/bry_pet.csv")


# Calculate median, rename column and remove lat and long columns
all_data = list(bio1=bio1, bio2=bio2, bio3=bio3, bio4=bio4, bio5=bio5,
                 bio6=bio6, bio7=bio7, bio8=bio8, bio9=bio9, bio10=bio10,
                 bio11=bio11, bio12=bio12, bio13=bio13, bio14=bio14, bio15=bio15,
                 bio16=bio16, bio17=bio17, bio18=bio18, bio19=bio19, ai=ai,
                 pet=pet)

median_results = map2(all_data, names(all_data), ~{
  .x %>% 
    group_by(species) %>%
    summarize(median_val = median(`...4`, na.rm = TRUE)) %>%
    rename(!!paste0(.y) := median_val)
})

combined_result = reduce(median_results, function(df1, df2) {
  left_join(df1, df2, by = "species")
})

# Calculate median and add
median_latitude = bio1 %>%
  group_by(species) %>%
  summarize(latitude = median(decimalLatitude, na.rm = TRUE))

median_chelsa = left_join(combined_result, median_latitude, by = "species")
write.csv(median_chelsa,"bry_median_CHELSA.csv")

# Biomes, with code from https://gis.stackexchange.com/questions/376280/converting-a-component-of-a-shapefile-into-a-geotiff-file-in-r
# Load packages
packs = list("tidyverse", "raster", "sf", "fasterize")
lapply(packs, require, character.only = T)

# Load shapefile # Data from https://globil-panda.opendata.arcgis.com/datasets/af92f20b1b43479581c819941e0f75ea_0/geoservice?geometry=11.953%2C-89.797%2C-11.953%2C78.032

wwf_shp = st_read("Terrestrial_Ecoregions_World.shp", stringsAsFactors = F) 
# Rasterize biome field and write to disk
raster(resolution = 1, crs = "+proj=longlat +datum=WGS84", vals = 0) %>% 
  fasterize(wwf_shp, ., field = "BIOME") %>% 
  writeRaster(., "biomes.tif")

biome = raster("tiff/biomes.tif") # http://worldmap.harvard.edu/data/geonode:wwf_terr_ecos_oRn

xy <- sp[c("decimalLongitude","decimalLatitude")]
sp_xy <- sp[c("species", "decimalLongitude","decimalLatitude")]
xy_bio <- raster::extract(biome, xy)
biome_bio <- na.omit(cbind(sp_xy, xy_bio))
write.csv(biome_bio, "bry_biome.csv")

getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}
all_biome = aggregate(biome_bio$'...4',list(biome_bio$species),getmode)
write.csv(all_biome,"biome_mode.csv",quote=FALSE)

# Elevation - from https://www.earthenv.org/topography
elev = raster("tiff/wc2.1_30s_elev.tif")
xy <- sp[c("decimalLongitude","decimalLatitude")]
sp_xy <- sp[c("species", "decimalLongitude","decimalLatitude")]
xy_bio <- raster::extract(elev, xy)
elev_bio <- cbind(sp_xy, xy_bio)
write.csv(elev_bio,"bry_elev.csv")

all_elev = aggregate(elev_bio$'...4',list(elev_bio$species),median)
write.csv(all_elev,"median_elev.csv")

# Sand, texture, water ph and organic carbon content 
#from https://github.com/ISRICWorldSoil/SoilGrids250m
sand = raster("tiff/SNDPPT_M_sl1_250m_ll.tif")
texture = raster("tiff/TEXMHT_M_sl1_250m_ll.tif")
water = raster("tiff/AWCh1_M_sl1_250m_ll.tif")
ph = raster("tiff/PHIHOX_M_sl1_250m_ll.tif")
organic_carbon = raster("tiff/ORCDRC_M_sl1_250m_ll.tif")

xy <- sp[c("decimalLongitude","decimalLatitude")]
sp_xy <- sp[c("species", "decimalLongitude","decimalLatitude")]

xy_bio <- raster::extract(sand, xy)
sand_bio <- cbind(sp_xy, xy_bio)

xy_bio <- raster::extract(texture, xy)
texture_bio <- cbind(sp_xy, xy_bio)

xy_bio <- raster::extract(water, xy)
water_bio <- cbind(sp_xy, xy_bio)

xy_bio <- raster::extract(ph, xy)
ph_bio <- cbind(sp_xy, xy_bio)

xy_bio <- raster::extract(organic_carbon, xy)
oc_bio <- cbind(sp_xy, xy_bio)

sand_bio = na.omit(sand_bio)
texture_bio = na.omit(texture_bio)
water_bio = na.omit(water_bio)
ph_bio = na.omit(ph_bio)
oc_bio = na.omit(oc_bio)

write.csv(sand_bio,"bry_sand_0cm.csv")
write.csv(water_bio,"bry_water_0cm.csv")
write.csv(texture_bio,"bry_texture.csv")
write.csv(ph_bio,"bry_ph_0cm.csv")
write.csv(oc_bio,"bry_oc.csv")

getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

all_sand = aggregate(sand_bio$'...4',list(sand_bio$species),median)
all_texture = aggregate(texture_bio$'...4',list(texture_bio$species),getmode)
all_water = aggregate(water_bio$'...4',list(water_bio$species),median)
all_ph = aggregate(ph_bio$'...4',list(ph_bio$species),median)
all_oc = aggregate(oc_bio$'...4',list(oc_bio$species),median)

write.csv(all_sand,"median_sand_0cm.csv")
write.csv(all_water,"median_water_0cm.csv")
write.csv(all_texture,"mode_texture.csv")
write.csv(all_ph,"median_ph.csv")
write.csv(all_oc,"median_oc.csv")

#合并海拔、土壤数据
water_bio = read.csv("median_water_0cm.csv") 
texture_bio = read.csv("mode_texture.csv")
ph_bio = read.csv("median_ph.csv") 
oc_bio = read.csv("median_oc.csv")
sand_bio = read.csv("median_sand_0cm.csv")
elev_bio = read.csv("median_elev.csv")

all_data = list(water_bio=water_bio, texture_bio=texture_bio, ph_bio=ph_bio,
                oc_bio=oc_bio, sand_bio=sand_bio, elev_bio=elev_bio)

median_results = map2(all_data, names(all_data), ~{
  .x %>% 
    group_by(Group.1) %>%
    summarize(median_val = median(`x`, na.rm = TRUE)) %>%
    rename(!!paste0(.y) := median_val)
})

combined_result = reduce(median_results, function(df1, df2) {
  left_join(df1, df2, by = "Group.1")
})

combined_result = combined_result %>%
  rename(species = Group.1)
         #NewName2 = OldName2,
         #NewName3 = OldName3,
write.csv(combined_result, "bry_median_CHELSA2.csv", row.names = FALSE)

#Climate data from WorldClim version 2.1
prec1 = raster("tiff/wc2.1_30s_prec_01.tif")
prec2 = raster("tiff/wc2.1_30s_prec_02.tif")
prec3 = raster("tiff/wc2.1_30s_prec_03.tif")
prec4 = raster("tiff/wc2.1_30s_prec_04.tif")
prec5 = raster("tiff/wc2.1_30s_prec_05.tif")
prec6 = raster("tiff/wc2.1_30s_prec_06.tif")
prec7 = raster("tiff/wc2.1_30s_prec_07.tif")
prec8 = raster("tiff/wc2.1_30s_prec_08.tif")
prec9 = raster("tiff/wc2.1_30s_prec_09.tif")
prec10 = raster("tiff/wc2.1_30s_prec_10.tif")
prec11 = raster("tiff/wc2.1_30s_prec_11.tif")
prec12 = raster("tiff/wc2.1_30s_prec_12.tif")
srad1 = raster("tiff/wc2.1_30s_srad_01.tif")
srad2 = raster("tiff/wc2.1_30s_srad_02.tif")
srad3 = raster("tiff/wc2.1_30s_srad_03.tif")
srad4 = raster("tiff/wc2.1_30s_srad_04.tif")
srad5 = raster("tiff/wc2.1_30s_srad_05.tif")
srad6 = raster("tiff/wc2.1_30s_srad_06.tif")
srad7 = raster("tiff/wc2.1_30s_srad_07.tif")
srad8 = raster("tiff/wc2.1_30s_srad_08.tif")
srad9 = raster("tiff/wc2.1_30s_srad_09.tif")
srad10 = raster("tiff/wc2.1_30s_srad_10.tif")
srad11 = raster("tiff/wc2.1_30s_srad_11.tif")
srad12 = raster("tiff/wc2.1_30s_srad_12.tif")
tavg1 = raster("tiff/wc2.1_30s_tavg_01.tif")
tavg2 = raster("tiff/wc2.1_30s_tavg_02.tif")
tavg3 = raster("tiff/wc2.1_30s_tavg_03.tif")
tavg4 = raster("tiff/wc2.1_30s_tavg_04.tif")
tavg5 = raster("tiff/wc2.1_30s_tavg_05.tif")
tavg6 = raster("tiff/wc2.1_30s_tavg_06.tif")
tavg7 = raster("tiff/wc2.1_30s_tavg_07.tif")
tavg8 = raster("tiff/wc2.1_30s_tavg_08.tif")
tavg9 = raster("tiff/wc2.1_30s_tavg_09.tif")
tavg10 = raster("tiff/wc2.1_30s_tavg_10.tif")
tavg11 = raster("tiff/wc2.1_30s_tavg_11.tif")
tavg12 = raster("tiff/wc2.1_30s_tavg_12.tif")
tmax1 = raster("tiff/wc2.1_30s_tmax_01.tif")
tmax2 = raster("tiff/wc2.1_30s_tmax_02.tif")
tmax3 = raster("tiff/wc2.1_30s_tmax_03.tif")
tmax4 = raster("tiff/wc2.1_30s_tmax_04.tif")
tmax5 = raster("tiff/wc2.1_30s_tmax_05.tif")
tmax6 = raster("tiff/wc2.1_30s_tmax_06.tif")
tmax7 = raster("tiff/wc2.1_30s_tmax_07.tif")
tmax8 = raster("tiff/wc2.1_30s_tmax_08.tif")
tmax9 = raster("tiff/wc2.1_30s_tmax_09.tif")
tmax10 = raster("tiff/wc2.1_30s_tmax_10.tif")
tmax11 = raster("tiff/wc2.1_30s_tmax_11.tif")
tmax12 = raster("tiff/wc2.1_30s_tmax_12.tif")
tmin1 = raster("tiff/wc2.1_30s_tmin_01.tif")
tmin2 = raster("tiff/wc2.1_30s_tmin_02.tif")
tmin3 = raster("tiff/wc2.1_30s_tmin_03.tif")
tmin4 = raster("tiff/wc2.1_30s_tmin_04.tif")
tmin5 = raster("tiff/wc2.1_30s_tmin_05.tif")
tmin6 = raster("tiff/wc2.1_30s_tmin_06.tif")
tmin7 = raster("tiff/wc2.1_30s_tmin_07.tif")
tmin8 = raster("tiff/wc2.1_30s_tmin_08.tif")
tmin9 = raster("tiff/wc2.1_30s_tmin_09.tif")
tmin10 = raster("tiff/wc2.1_30s_tmin_10.tif")
tmin11 = raster("tiff/wc2.1_30s_tmin_11.tif")
tmin12 = raster("tiff/wc2.1_30s_tmin_12.tif")
vapr1 = raster("tiff/wc2.1_30s_vapr_01.tif")
vapr2 = raster("tiff/wc2.1_30s_vapr_02.tif")
vapr3 = raster("tiff/wc2.1_30s_vapr_03.tif")
vapr4 = raster("tiff/wc2.1_30s_vapr_04.tif")
vapr5 = raster("tiff/wc2.1_30s_vapr_05.tif")
vapr6 = raster("tiff/wc2.1_30s_vapr_06.tif")
vapr7 = raster("tiff/wc2.1_30s_vapr_07.tif")
vapr8 = raster("tiff/wc2.1_30s_vapr_08.tif")
vapr9 = raster("tiff/wc2.1_30s_vapr_09.tif")
vapr10 = raster("tiff/wc2.1_30s_vapr_10.tif")
vapr11 = raster("tiff/wc2.1_30s_vapr_11.tif")
vapr12 = raster("tiff/wc2.1_30s_vapr_12.tif")
wind1 = raster("tiff/wc2.1_30s_wind_01.tif")
wind2 = raster("tiff/wc2.1_30s_wind_02.tif")
wind3 = raster("tiff/wc2.1_30s_wind_03.tif")
wind4 = raster("tiff/wc2.1_30s_wind_04.tif")
wind5 = raster("tiff/wc2.1_30s_wind_05.tif")
wind6 = raster("tiff/wc2.1_30s_wind_06.tif")
wind7 = raster("tiff/wc2.1_30s_wind_07.tif")
wind8 = raster("tiff/wc2.1_30s_wind_08.tif")
wind9 = raster("tiff/wc2.1_30s_wind_09.tif")
wind10 = raster("tiff/wc2.1_30s_wind_10.tif")
wind11 = raster("tiff/wc2.1_30s_wind_11.tif")
wind12 = raster("tiff/wc2.1_30s_wind_12.tif")

xy <- sp[c("decimalLongitude","decimalLatitude")]
sp_xy <- sp[c("species", "decimalLongitude","decimalLatitude")]

xy_prec1 <- raster::extract(prec1, xy)
xy_prec2 <- raster::extract(prec2, xy)
xy_prec3 <- raster::extract(prec3, xy)
xy_prec4 <- raster::extract(prec4, xy)
xy_prec5 <- raster::extract(prec5, xy)
xy_prec6 <- raster::extract(prec6, xy)
xy_prec7 <- raster::extract(prec7, xy)
xy_prec8 <- raster::extract(prec8, xy)
xy_prec9 <- raster::extract(prec9, xy)
xy_prec10 <- raster::extract(prec10, xy)
xy_prec11 <- raster::extract(prec11, xy)
xy_prec12 <- raster::extract(prec12, xy)
xy_srad1 <- raster::extract(srad1, xy)
xy_srad2 <- raster::extract(srad2, xy)
xy_srad3 <- raster::extract(srad3, xy)
xy_srad4 <- raster::extract(srad4, xy)
xy_srad5 <- raster::extract(srad5, xy)
xy_srad6 <- raster::extract(srad6, xy)
xy_srad7 <- raster::extract(srad7, xy)
xy_srad8 <- raster::extract(srad8, xy)
xy_srad9 <- raster::extract(srad9, xy)
xy_srad10 <- raster::extract(srad10, xy)
xy_srad11 <- raster::extract(srad11, xy)
xy_srad12 <- raster::extract(srad12, xy)
xy_tavg1 <- raster::extract(tavg1, xy)
xy_tavg2 <- raster::extract(tavg2, xy)
xy_tavg3 <- raster::extract(tavg3, xy)
xy_tavg4 <- raster::extract(tavg4, xy)
xy_tavg5 <- raster::extract(tavg5, xy)
xy_tavg6 <- raster::extract(tavg6, xy)
xy_tavg7 <- raster::extract(tavg7, xy)
xy_tavg8 <- raster::extract(tavg8, xy)
xy_tavg9 <- raster::extract(tavg9, xy)
xy_tavg10 <- raster::extract(tavg10, xy)
xy_tavg11 <- raster::extract(tavg11, xy)
xy_tavg12 <- raster::extract(tavg12, xy)
xy_tmax1 <- raster::extract(tmax1, xy)
xy_tmax2 <- raster::extract(tmax2, xy)
xy_tmax3 <- raster::extract(tmax3, xy)
xy_tmax4 <- raster::extract(tmax4, xy)
xy_tmax5 <- raster::extract(tmax5, xy)
xy_tmax6 <- raster::extract(tmax6, xy)
xy_tmax7 <- raster::extract(tmax7, xy)
xy_tmax8 <- raster::extract(tmax8, xy)
xy_tmax9 <- raster::extract(tmax9, xy)
xy_tmax10 <- raster::extract(tmax10, xy)
xy_tmax11 <- raster::extract(tmax11, xy)
xy_tmax12 <- raster::extract(tmax12, xy)
xy_tmin1 <- raster::extract(tmin1, xy)
xy_tmin2 <- raster::extract(tmin2, xy)
xy_tmin3 <- raster::extract(tmin3, xy)
xy_tmin4 <- raster::extract(tmin4, xy)
xy_tmin5 <- raster::extract(tmin5, xy)
xy_tmin6 <- raster::extract(tmin6, xy)
xy_tmin7 <- raster::extract(tmin7, xy)
xy_tmin8 <- raster::extract(tmin8, xy)
xy_tmin9 <- raster::extract(tmin9, xy)
xy_tmin10 <- raster::extract(tmin10, xy)
xy_tmin11 <- raster::extract(tmin11, xy)
xy_tmin12 <- raster::extract(tmin12, xy)
xy_vapr1 <- raster::extract(vapr1, xy)
xy_vapr2 <- raster::extract(vapr2, xy)
xy_vapr3 <- raster::extract(vapr3, xy)
xy_vapr4 <- raster::extract(vapr4, xy)
xy_vapr5 <- raster::extract(vapr5, xy)
xy_vapr6 <- raster::extract(vapr6, xy)
xy_vapr7 <- raster::extract(vapr7, xy)
xy_vapr8 <- raster::extract(vapr8, xy)
xy_vapr9 <- raster::extract(vapr9, xy)
xy_vapr10 <- raster::extract(vapr10, xy)
xy_vapr11 <- raster::extract(vapr11, xy)
xy_vapr12 <- raster::extract(vapr12, xy)
xy_wind1 <- raster::extract(wind1, xy)
xy_wind2 <- raster::extract(wind2, xy)
xy_wind3 <- raster::extract(wind3, xy)
xy_wind4 <- raster::extract(wind4, xy)
xy_wind5 <- raster::extract(wind5, xy)
xy_wind6 <- raster::extract(wind6, xy)
xy_wind7 <- raster::extract(wind7, xy)
xy_wind8 <- raster::extract(wind8, xy)
xy_wind9 <- raster::extract(wind9, xy)
xy_wind10 <- raster::extract(wind10, xy)
xy_wind11 <- raster::extract(wind11, xy)
xy_wind12 <- raster::extract(wind12, xy)

prec1 <- cbind(sp_xy, xy_prec1)
prec2 <- cbind(sp_xy, xy_prec2)
prec3 <- cbind(sp_xy, xy_prec3)
prec4 <- cbind(sp_xy, xy_prec4)
prec5 <- cbind(sp_xy, xy_prec5)
prec6 <- cbind(sp_xy, xy_prec6)
prec7 <- cbind(sp_xy, xy_prec7)
prec8 <- cbind(sp_xy, xy_prec8)
prec9 <- cbind(sp_xy, xy_prec9)
prec10 <- cbind(sp_xy, xy_prec10)
prec11 <- cbind(sp_xy, xy_prec11)
prec12 <- cbind(sp_xy, xy_prec12)
srad1 <- cbind(sp_xy, xy_srad1)
srad2 <- cbind(sp_xy, xy_srad2)
srad3 <- cbind(sp_xy, xy_srad3)
srad4 <- cbind(sp_xy, xy_srad4)
srad5 <- cbind(sp_xy, xy_srad5)
srad6 <- cbind(sp_xy, xy_srad6)
srad7 <- cbind(sp_xy, xy_srad7)
srad8 <- cbind(sp_xy, xy_srad8)
srad9 <- cbind(sp_xy, xy_srad9)
srad10 <- cbind(sp_xy, xy_srad10)
srad11 <- cbind(sp_xy, xy_srad11)
srad12 <- cbind(sp_xy, xy_srad12)
tavg1 <- cbind(sp_xy, xy_tavg1)
tavg2 <- cbind(sp_xy, xy_tavg2)
tavg3 <- cbind(sp_xy, xy_tavg3)
tavg4 <- cbind(sp_xy, xy_tavg4)
tavg5 <- cbind(sp_xy, xy_tavg5)
tavg6 <- cbind(sp_xy, xy_tavg6)
tavg7 <- cbind(sp_xy, xy_tavg7)
tavg8 <- cbind(sp_xy, xy_tavg8)
tavg9 <- cbind(sp_xy, xy_tavg9)
tavg10 <- cbind(sp_xy, xy_tavg10)
tavg11 <- cbind(sp_xy, xy_tavg11)
tavg12 <- cbind(sp_xy, xy_tavg12)
tmax1 <- cbind(sp_xy, xy_tmax1)
tmax2 <- cbind(sp_xy, xy_tmax2)
tmax3 <- cbind(sp_xy, xy_tmax3)
tmax4 <- cbind(sp_xy, xy_tmax4)
tmax5 <- cbind(sp_xy, xy_tmax5)
tmax6 <- cbind(sp_xy, xy_tmax6)
tmax7 <- cbind(sp_xy, xy_tmax7)
tmax8 <- cbind(sp_xy, xy_tmax8)
tmax9 <- cbind(sp_xy, xy_tmax9)
tmax10 <- cbind(sp_xy, xy_tmax10)
tmax11 <- cbind(sp_xy, xy_tmax11)
tmax12 <- cbind(sp_xy, xy_tmax12)
tmin1 <- cbind(sp_xy, xy_tmin1)
tmin2 <- cbind(sp_xy, xy_tmin2)
tmin3 <- cbind(sp_xy, xy_tmin3)
tmin4 <- cbind(sp_xy, xy_tmin4)
tmin5 <- cbind(sp_xy, xy_tmin5)
tmin6 <- cbind(sp_xy, xy_tmin6)
tmin7 <- cbind(sp_xy, xy_tmin7)
tmin8 <- cbind(sp_xy, xy_tmin8)
tmin9 <- cbind(sp_xy, xy_tmin9)
tmin10 <- cbind(sp_xy, xy_tmin10)
tmin11 <- cbind(sp_xy, xy_tmin11)
tmin12 <- cbind(sp_xy, xy_tmin12)
vapr1 <- cbind(sp_xy, xy_vapr1)
vapr2 <- cbind(sp_xy, xy_vapr2)
vapr3 <- cbind(sp_xy, xy_vapr3)
vapr4 <- cbind(sp_xy, xy_vapr4)
vapr5 <- cbind(sp_xy, xy_vapr5)
vapr6 <- cbind(sp_xy, xy_vapr6)
vapr7 <- cbind(sp_xy, xy_vapr7)
vapr8 <- cbind(sp_xy, xy_vapr8)
vapr9 <- cbind(sp_xy, xy_vapr9)
vapr10 <- cbind(sp_xy, xy_vapr10)
vapr11 <- cbind(sp_xy, xy_vapr11)
vapr12 <- cbind(sp_xy, xy_vapr12)
wind1 <- cbind(sp_xy, xy_wind1)
wind2 <- cbind(sp_xy, xy_wind2)
wind3 <- cbind(sp_xy, xy_wind3)
wind4 <- cbind(sp_xy, xy_wind4)
wind5 <- cbind(sp_xy, xy_wind5)
wind6 <- cbind(sp_xy, xy_wind6)
wind7 <- cbind(sp_xy, xy_wind7)
wind8 <- cbind(sp_xy, xy_wind8)
wind9 <- cbind(sp_xy, xy_wind9)
wind10 <- cbind(sp_xy, xy_wind10)
wind11 <- cbind(sp_xy, xy_wind11)
wind12 <- cbind(sp_xy, xy_wind12)

write.csv(prec1,"bry_prec1.csv")
write.csv(prec2,"bry_prec2.csv")
write.csv(prec3,"bry_prec3.csv")
write.csv(prec4,"bry_prec4.csv")
write.csv(prec5,"bry_prec5.csv")
write.csv(prec6,"bry_prec6.csv")
write.csv(prec7,"bry_prec7.csv")
write.csv(prec8,"bry_prec8.csv")
write.csv(prec9,"bry_prec9.csv")
write.csv(prec10,"bry_prec10.csv")
write.csv(prec11,"bry_prec11.csv")
write.csv(prec12,"bry_prec12.csv")
write.csv(srad1,"bry_srad1.csv")
write.csv(srad2,"bry_srad2.csv")
write.csv(srad3,"bry_srad3.csv")
write.csv(srad4,"bry_srad4.csv")
write.csv(srad5,"bry_srad5.csv")
write.csv(srad6,"bry_srad6.csv")
write.csv(srad7,"bry_srad7.csv")
write.csv(srad8,"bry_srad8.csv")
write.csv(srad9,"bry_srad9.csv")
write.csv(srad10,"bry_srad10.csv")
write.csv(srad11,"bry_srad11.csv")
write.csv(srad12,"bry_srad12.csv")
write.csv(tavg1,"bry_tavg1.csv")
write.csv(tavg2,"bry_tavg2.csv")
write.csv(tavg3,"bry_tavg3.csv")
write.csv(tavg4,"bry_tavg4.csv")
write.csv(tavg5,"bry_tavg5.csv")
write.csv(tavg6,"bry_tavg6.csv")
write.csv(tavg7,"bry_tavg7.csv")
write.csv(tavg8,"bry_tavg8.csv")
write.csv(tavg9,"bry_tavg9.csv")
write.csv(tavg10,"bry_tavg10.csv")
write.csv(tavg11,"bry_tavg11.csv")
write.csv(tavg12,"bry_tavg12.csv")
write.csv(tmax1,"bry_tmax1.csv")
write.csv(tmax2,"bry_tmax2.csv")
write.csv(tmax3,"bry_tmax3.csv")
write.csv(tmax4,"bry_tmax4.csv")
write.csv(tmax5,"bry_tmax5.csv")
write.csv(tmax6,"bry_tmax6.csv")
write.csv(tmax7,"bry_tmax7.csv")
write.csv(tmax8,"bry_tmax8.csv")
write.csv(tmax9,"bry_tmax9.csv")
write.csv(tmax10,"bry_tmax10.csv")
write.csv(tmax11,"bry_tmax11.csv")
write.csv(tmax12,"bry_tmax12.csv")
write.csv(tmin1,"bry_tmin1.csv")
write.csv(tmin2,"bry_tmin2.csv")
write.csv(tmin3,"bry_tmin3.csv")
write.csv(tmin4,"bry_tmin4.csv")
write.csv(tmin5,"bry_tmin5.csv")
write.csv(tmin6,"bry_tmin6.csv")
write.csv(tmin7,"bry_tmin7.csv")
write.csv(tmin8,"bry_tmin8.csv")
write.csv(tmin9,"bry_tmin9.csv")
write.csv(tmin10,"bry_tmin10.csv")
write.csv(tmin11,"bry_tmin11.csv")
write.csv(tmin12,"bry_tmin12.csv")
write.csv(vapr1,"bry_vapr1.csv")
write.csv(vapr2,"bry_vapr2.csv")
write.csv(vapr3,"bry_vapr3.csv")
write.csv(vapr4,"bry_vapr4.csv")
write.csv(vapr5,"bry_vapr5.csv")
write.csv(vapr6,"bry_vapr6.csv")
write.csv(vapr7,"bry_vapr7.csv")
write.csv(vapr8,"bry_vapr8.csv")
write.csv(vapr9,"bry_vapr9.csv")
write.csv(vapr10,"bry_vapr10.csv")
write.csv(vapr11,"bry_vapr11.csv")
write.csv(vapr12,"bry_vapr12.csv")
write.csv(wind1,"bry_wind1.csv")
write.csv(wind2,"bry_wind2.csv")
write.csv(wind3,"bry_wind3.csv")
write.csv(wind4,"bry_wind4.csv")
write.csv(wind5,"bry_wind5.csv")
write.csv(wind6,"bry_wind6.csv")
write.csv(wind7,"bry_wind7.csv")
write.csv(wind8,"bry_wind8.csv")
write.csv(wind9,"bry_wind9.csv")
write.csv(wind10,"bry_wind10.csv")
write.csv(wind11,"bry_wind11.csv")
write.csv(wind12,"bry_wind12.csv")

all_data = list(prec1=prec1, prec2=prec2, prec3=prec3, prec4=prec4, prec5=prec5, prec6=prec6,
                prec7=prec7, prec8=prec8, prec9=prec9, prec10=prec10, prec11=prec11, prec12=prec12,
                tavg1=tavg1, tavg2=tavg2, tavg3=tavg3, tavg4=tavg4, tavg5=tavg5, tavg6=tavg6,
                tavg7=tavg7, tavg8=tavg8, tavg9=tavg9, tavg10=tavg10, tavg11=tavg11, tavg12=tavg12,
                tmax1=tmax1, tmax2=tmax2, tmax3=tmax3, tmax4=tmax4, tmax5=tmax5, tmax6=tmax6,
                tmax7=tmax7, tmax8=tmax8, tmax9=tmax9, tmax10=tmax10, tmax11=tmax11, tmax12=tmax12,
                tmin1=tmin1, tmin2=tmin2, tmin3=tmin3, tmin4=tmin4, tmin5=tmin5, tmin6=tmin6,
                tmin7=tmin7, tmin8=tmin8, tmin9=tmin9, tmin10=tmin10, tmin11=tmin11, tmin12=tmin12,
                srad1=srad1, srad2=srad2, srad3=srad3, srad4=srad4, srad5=srad5, srad6=srad6,
                srad7=srad7, srad8=srad8, srad9=srad9, srad10=srad10, srad11=srad11, srad12=srad12,
                vapr1=vapr1, vapr2=vapr2, vapr3=vapr3, vapr4=vapr4, vapr5=vapr5, vapr6=vapr6,
                vapr7=vapr7, vapr8=vapr8, vapr9=vapr9, vapr10=vapr10, vapr11=vapr11, vapr12=vapr12,
                wind1=wind1, wind2=wind2, wind3=wind3, wind4=wind4, wind5=wind5, wind6=wind6,
                wind7=wind7, wind8=wind8, wind9=wind9, wind10=wind10, wind11=wind11, wind12=wind12,
                pet=pet, ai=ai)

all_data = list(srad1=srad1, srad2=srad2, srad3=srad3, srad4=srad4, srad5=srad5, srad6=srad6,
                srad7=srad7, srad8=srad8, srad9=srad9, srad10=srad10, srad11=srad11, srad12=srad12,
                vapr1=vapr1, vapr2=vapr2, vapr3=vapr3, vapr4=vapr4, vapr5=vapr5, vapr6=vapr6,
                vapr7=vapr7, vapr8=vapr8, vapr9=vapr9, vapr10=vapr10, vapr11=vapr11, vapr12=vapr12,
                wind1=wind1, wind2=wind2, wind3=wind3, wind4=wind4, wind5=wind5, wind6=wind6,
                wind7=wind7, wind8=wind8, wind9=wind9, wind10=wind10, wind11=wind11, wind12=wind12
                )
median_results = map2(all_data, names(all_data), ~{
  .x %>% 
    group_by(species) %>%
    summarize(median_val = median(`...4`, na.rm = TRUE)) %>%
    rename(!!paste0(.y) := median_val)
})

combined_result = reduce(median_results, function(df1, df2) {
  left_join(df1, df2, by = "species")
})

median_latitude = srad1 %>%
  group_by(species) %>%
  summarize(latitude = median(decimalLatitude, na.rm = TRUE))

median_chelsa = left_join(combined_result, median_latitude, by = "species")

# 第一列是物种名，sra1为第二列
columns_to_average <- c(2:13)
subset_data <- median_chelsa[, columns_to_average]
row_averages <- rowMeans(subset_data)
median_chelsa$srad <- row_averages

columns_to_average <- c(14:25)
subset_data <- median_chelsa[, columns_to_average]
row_averages <- rowMeans(subset_data)
median_chelsa$vapr <- row_averages

columns_to_average <- c(26:37)
subset_data <- median_chelsa[, columns_to_average]
row_averages <- rowMeans(subset_data)
median_chelsa$wind <- row_averages


write.csv(median_chelsa,"bry_median_CHELSA1.csv")

#合并所有的csv文件
library(dplyr)
library(readr)

# 假设所有CSV文件都在当前工作目录的一个名为"csv_files"的文件夹中
file_list <- list.files(path = "other", pattern = "\\.csv$", full.names = TRUE)
# 读取第一个CSV文件，作为合并的基础
combined_data <- read_csv(file_list[1])
# 从第二个CSV文件开始，依次合并
for (file in file_list[-1]) {
  # 读取当前CSV文件
  current_data <- read_csv(file)
  
  # 基于species列进行合并，这里使用left_join，如果你需要保留所有数据行，可以使用full_join
  combined_data <- left_join(combined_data, current_data, by = "species")
}
# 保存合并后的数据到一个新的CSV文件中
write_csv(combined_data, "combined_data.csv")

