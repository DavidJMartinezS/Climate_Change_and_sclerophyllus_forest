library(tidyverse)
library(sf)
library(ncdf4)
library(SPEI)
library(zoo)
library(scico)
library(devEMF)
library(raster)

# download study area
url <- "http://www.ide.cl/descargas/capas/MMA/PisosVegetacionalesPliscoff2017.zip"
download.file(url,destfile = 'data/Pisos_LyP2017.zip')
unzip('data/Pisos_LyP2017.zip',exdir = 'data/')
read_sf('data/PisosVegetacionalesPliscoff2017.shp') %>%
  filter(codigo=="P41") %>%
  st_transform(4326) %>%
  write_sf('data/bemaql.shp')

# determining critical thresholds
BEMAQL <- read_sf("data/bemaql.shp")
limits <- st_bbox(BEMAQL)

pp.nc <- nc_open("data/CR2MET_pr_v2.0_mon_1979_2019_005deg.nc")
tm.nc <- nc_open("data/CR2MET_t2m_v2.0_mon_1979_2019_005deg.nc")

pp <- ncvar_get(pp.nc, "pr")
tm <- ncvar_get(tm.nc, "t2m")
lat <- ncvar_get(pp.nc, "lat")
lon <- ncvar_get(pp.nc, "lon")
time <- ncvar_get(pp.nc, "time") 

lat.sa<-lat[443:485]
lon.sa<-lon[108:142]
pp.sa<-pp[c(108:142),c(443:485),]
tm.sa<-tm[c(108:142),c(443:485),]
Nlon <- length(lon.sa)
Nlat <- length(lat.sa)
Ntime <- length(time)

um.cri.10 <- matrix(NA, nrow = Nlon, ncol = Nlat)
dimnames(um.cri.10) <- list(lon.sa, lat.sa)
for (k in 1:Nlat){
  for (j in 1:Nlon) {
    pet <- thornthwaite(tm.sa[j, k, ], lat.sa[k], na.rm = T)
    bhs <- pp.sa[j, k, ] - pet
    bhs.sm <- rollsum(bhs, 12, fill = NA, align = "right")
    um.cri.10[j, k] <- quantile(bhs.sm, .1, na.rm = T)
  }
}
save(um.cri.10, file = "data/um_cri_grid_10.RData")

um.cri.25 <- matrix(NA, nrow = Nlon, ncol = Nlat)
dimnames(um.cri.25) <- list(lon.sa, lat.sa)
for (k in 1:Nlat){
  for (j in 1:Nlon) {
    pet <- thornthwaite(tm.sa[j, k, ], lat.sa[k], na.rm = T)
    bhs <- pp.sa[j, k, ] - pet
    bhs.sm <- rollsum(bhs, 12, fill = NA, align = "right")
    um.cri.25[j, k] <- quantile(bhs.sm, .25, na.rm = T)
  }
}
save(um.cri.25, file = "data/um_cri_grid_25.RData")

# Time of Emergence ----
pr.can<-"~/Memoria/Memoria/R_work/GCMs y downscaling/CanESM5_pr_SSP585_r2i1p1f1_UQM.nc"
tm.can<-"~/Memoria/Memoria/R_work/GCMs y downscaling/CanESM5_tas_SSP585_r2i1p1f1_UQM.nc"
pr.ipsl<-"~/Memoria/Memoria/R_work/GCMs y downscaling/IPSL_CM6A_LR_pr_SSP585_r1i1p1f1_UQM.nc"
tm.ipsl<-"~/Memoria/Memoria/R_work/GCMs y downscaling/IPSL_CM6A_LR_tas_SSP585_r1i1p1f1_UQM.nc"
pr.inm<-"~/Memoria/Memoria/R_work/GCMs y downscaling/INM_CM4_8_pr_SSP585_r1i1p1f1_UQM.nc"
tm.inm<-"~/Memoria/Memoria/R_work/GCMs y downscaling/INM_CM4_8_tas_SSP585_r1i1p1f1_UQM.nc"
pr.cesm<-"~/Memoria/Memoria/R_work/GCMs y downscaling/CESM2_WACCM_pr_SSP245_r3i1p1f1_UQM.nc"
tm.cesm<-"~/Memoria/Memoria/R_work/GCMs y downscaling/CESM2_WACCM_tas_SSP245_r3i1p1f1_UQM.nc"
pr.nor<-"~/Memoria/Memoria/R_work/GCMs y downscaling/NorESM2_LM_pr_SSP245_r2i1p1f1_UQM.nc"
tm.nor<-"~/Memoria/Memoria/R_work/GCMs y downscaling/NorESM2_LM_tas_SSP245_r2i1p1f1_UQM.nc"

pr.gcm <- '~/Memoria/Memoria/R_work/GCMs y downscaling/CanESM5_pr_SSP585_r2i1p1f1_UQM.nc'
tm.gcm <- '~/Memoria/Memoria/R_work/GCMs y downscaling/CanESM5_tas_SSP585_r2i1p1f1_UQM.nc'
ct <- um.cri.10
percentile <- "P10%"

toe.gcm.grid.fun <- function(pr.gcm, tm.gcm, ct, percentile) {
  nc.pr <- nc_open(pr.gcm)
  nc.tm <- nc_open(tm.gcm)
  
  pr <- ncvar_get(nc.pr, "pr")
  tm <- ncvar_get(nc.tm, "t2m")
  lon <- ncvar_get(nc.pr, "lon")
  lat <- ncvar_get(nc.pr, "lat")
  time <- ncvar_get(nc.pr, "time")
  Nlat <- length(lat)
  Nlon <- length(lon)
  Ntime <- length(time)
  
  date <- c(seq(as.Date("1979-01-01"), as.Date("2100-12-01"), by = "1 month"))
  years <- as.numeric(format(date, '%Y'))
  
  bhs.gcm <- array(NA, c(Nlon, Nlat, Ntime))
  for (k in 1:Nlat){
    for (j in 1:Nlon) {
      pet <- thornthwaite(tm[j, k, ], lat[k], na.rm = T)
      bhs <- pr[j, k, ] - pet
      bhs.sm <- rollsum(bhs, 12, fill = NA, align = "right")
      bhs.mm <- rollmean(bhs.sm, 360, fill = NA, align = "right", na.rm = T)
      bhs.gcm[j, k, ] <- bhs.mm
    }
  }
  
  toe.matrix <- matrix(NA, ncol = Nlat, nrow = Nlon)
  dimnames(toe.matrix) <- list(lon,lat)
  for (k in 1:Nlat){
    for (j in 1:Nlon) {
      toe <- c(NA)
      for (l in 1:Ntime) {
        if (is.na(bhs.gcm[j, k, l]) == TRUE) {
          toe <- c(toe, NA)
        }
        else if (bhs.gcm[j, k, l] > ct[j, k]){
          toe <- c(toe, NA)
        }
        else{
          toe <- c(toe, years[l])
        }
        toe1 <- toe[tail(which(is.na(toe)), 1) + 1]
        toe.matrix[j, k] <- toe1
      }
    }
  }
  
  gcm.name <- pr.gcm %>% basename() %>% str_extract("(?<=).+(?=_pr)")
  toe.df <- reshape2::melt(toe.matrix, varnames = c('x', 'y'), value.name = 'ToE') %>% 
    mutate(umbral = percentile, gcm = gcm.name)
  return(toe.df)
}

toe.10.can.df <- toe.gcm.grid.fun(pr.can, tm.can, um.cri.10, "P10%") 
toe.10.ipsl.df <- toe.gcm.grid.fun(pr.ipsl, tm.ipsl, um.cri.10, "P10%") 
toe.10.inm.df <- toe.gcm.grid.fun(pr.inm, tm.inm, um.cri.10, "P10%") 
toe.10.cesm.df <- toe.gcm.grid.fun(pr.cesm, tm.cesm, um.cri.10, "P10%") 
toe.10.nor.df <- toe.gcm.grid.fun(pr.nor, tm.nor, um.cri.10, "P10%") 

toe.25.can.df <- toe.gcm.grid.fun(pr.can, tm.can, um.cri.25, "P25%") 
toe.25.ipsl.df <- toe.gcm.grid.fun(pr.ipsl, tm.ipsl, um.cri.25, "P25%") 
toe.25.inm.df <- toe.gcm.grid.fun(pr.inm, tm.inm, um.cri.25, "P25%") 
toe.25.cesm.df <- toe.gcm.grid.fun(pr.cesm, tm.cesm, um.cri.25, "P25%") 
toe.25.nor.df <- toe.gcm.grid.fun(pr.nor, tm.nor, um.cri.25, "P25%") 

toe.uni.df <- rbind(
  toe.10.can.df,
  toe.10.ipsl.df,
  toe.10.inm.df,
  toe.10.cesm.df,
  toe.10.nor.df,
  toe.25.can.df,
  toe.25.ipsl.df,
  toe.25.inm.df,
  toe.25.cesm.df,
  toe.25.nor.df
) %>% 
  mutate(gcm = str_replace_all(gcm, '_', '-'),
         gcm = factor(gcm, levels = c("NorESM2-LM", "CESM2", "CESM2-WACCM", "INM-CM4-8", "IPSL-CM6A-LR", "CanESM5")),
         umbral = factor(umbral, levels = c("P25%", "P10%")))

map.toe <- toe.uni.df %>% 
  filter(x > limits[1]-0.05, x < limits[3]+0.05,
         y > limits[2]-0.05, y < limits[4]+0.05) %>% 
  ggplot() +
  geom_raster(aes(x = x, y = y, fill = ToE),
              show.legend = T) +
  scale_fill_stepsn(
    colours = c("#9c6644", "#b08968", "#ddb892", "#e6ccb2", "#ede0d4"),
    na.value = "white",
    n.breaks = 5,
    limits = c(2050,2100),
    guide=guide_colorsteps(show.limits = T)
  )+
  # scale_fill_binned(
  #   low = "red",
  #   high = "yellow",
  #   breaks = seq(2050, 2100, 10),
  #   na.value = "white"
  # ) +
  geom_sf(
    data = BEMAQL,
    color = 'black',
    fill = NA,
    size = 0.1
  ) +
  scale_x_continuous(breaks = seq(-71.4, -70.2, 0.6), expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0))+
  facet_grid(umbral ~ gcm) +
  labs(x = "", y = "") +
  theme_bw() +
  theme(
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.position = "bottom",
    panel.spacing.x = unit(1.2, "lines"),
    strip.background = element_blank(),
    axis.text = element_text(size = 7)
  ) +
  guides(
    fill = guide_colourbar(
      direction = 'horizontal',
      title = 'ToE',
      title.position = 'top',
      title.hjust = 0.5,
      ticks.colour = '#f5f5f2',
      ticks.linewidth = 2,
      barwidth = 25,
      barheight = 0.5
    )
  )
map.toe
#ggsave("map_toe.emf", map.toe, width = 8, height = 6, device = {function(filename, ...) devEMF::emf(file = filename, ...)})
#ggsave("map_toe.png", map.toe, width = 8, height = 6, dpi = 700)
ggsave("map_toe.pdf", map.toe, width = 8, height = 6)

