library(tidyverse)
library(sf)
library(ncdf4)
library(SPEI)
library(zoo)
library(ggh4x)
library(raster)
library(reshape2)

pp.nc <- nc_open("data/CR2MET_pr_v2.0_mon_1979_2019_005deg.nc")
tm.nc <- nc_open("data/CR2MET_t2m_v2.0_mon_1979_2019_005deg.nc")

pp <- ncvar_get(pp.nc, "pr")
tm <- ncvar_get(tm.nc, "t2m")
lat <- ncvar_get(pp.nc, "lat")
lon <- ncvar_get(pp.nc, "lon")
time <- ncvar_get(pp.nc, "time") 

lat.sa <- lat[443:485]
lon.sa <- lon[108:142]
pp.sa <- pp[c(108:142), c(443:485),]
tm.sa <- tm[c(108:142), c(443:485),]
Nlon <- length(lon.sa)
Nlat <- length(lat.sa)
Ntime <- length(time)
volt.lat <- sort(1:Nlat, decreasing = T)
crs<-CRS('+proj=longlat +datum=WGS84 +no_defs')
date <- c(seq(as.Date("1979-01-01"), as.Date("2019-12-01"), by = "1 month")) %>% as.character()

bhs.cr2 <- array(NA, c(Nlon, Nlat, Ntime))
dimnames(bhs.cr2) <- list(lon.sa, lat.sa, date)
for (k in 1:Nlat){
  for (j in 1:Nlon) {
    pet <- thornthwaite(tm.sa[j, k, ], lat.sa[k], na.rm = T)
    bhs <- pp.sa[j, k, ] - pet
    bhs.sm <- rollsum(bhs, 12, fill = NA, align =  "right")
    bhs.cr2[j, k, ] <- bhs.sm
  }
}

BEMAQL <- read_sf("data/bemaql.shp")
ptos_sf <- read.csv("data/pt_method_2.csv",header = T,sep = ";") %>% 
  st_as_sf(coords = c("Lon", "Lat"), crs = 4326, remove = F) 

stack.list<-c()
for (i in 1:Ntime) {
  ras.cr2 <-
    raster::raster(
      t(bhs.cr2[, volt.lat, i]),
      xmn = min(lon.sa),
      xmx = max(lon.sa),
      ymn = min(lat.sa),
      ymx = max(lat.sa),
      crs = crs
    ) %>% raster::crop(BEMAQL, snap = 'out')
  stack.list <- c(stack.list, ras.cr2)
}
stack.bhs.cr2 <- raster::stack(stack.list)
coordinates(ptos)= ~ Lon + Lat
crs(ptos)=crs
stack.bhs.pts<-raster::extract(stack.bhs.cr2,ptos_sf)

apply(stack.bhs.pts, 2, mean) %>% hist(breaks=50)
apply(stack.bhs.pts, 2, mean) %>% as.data.frame() %>% `colnames<-`("bhs") %>% cbind(fecha=seq(as.Date("1979-01-01"),as.Date("2019-12-01"),by="1 month")) %>% ggplot()+geom_line(aes(x=fecha,y=bhs))
apply(stack.bhs.pts, 2, mean)[373:492] %>% quantile(c(0.5,0.1),na.rm=T)#-383.3944 y -472.7023 desde 2010 ***

CT_P50<-(-383)
CT_P10<-(-473)

toe.gcm.grid.fun<-function(pr.gcm,tm.gcm,hum,percentil,m.average){
  nc.pr<-nc_open(pr.gcm)
  nc.tm<-nc_open(tm.gcm)
  
  pr<-ncvar_get(nc.pr,"pr")
  tm<-ncvar_get(nc.tm,"t2m")
  lon<-ncvar_get(nc.pr,"lon")
  lat<-ncvar_get(nc.pr,"lat")
  time<-ncvar_get(nc.pr,"time")
  Nlat<-length(lat)
  Nlon<-length(lon)
  Ntime<-length(time)
  
  date<-c(seq(as.Date("1979-01-01"),as.Date("2100-12-01"),by="1 month"))
  years<-as.numeric(format(fecha,'%Y'))
  
  bhs.gcm<-array(NA,c(Nlon,Nlat,Ntime))
  for (k in 1:Nlat){
    for (j in 1:Nlon) {
      pet <- thornthwaite(tm[j, k, ], lat[k], na.rm = T)
      bhs <- pr[j, k, ] - pet
      bhs.sm <- rollsum(bhs, 12, fill = NA, align =  "right")
      bhs.mm <- rollmean(bhs.sm, (m.average * 12), fill = NA, align = "right", na.rm=T)
      bhs.gcm[j, k, ] <- bhs.mm
    }
  }
  
  toe.matrix <- matrix(NA, ncol = Nlat, nrow = Nlon)
  dimnames(toe.matrix) <- list(lon, lat)
  for (k in 1:Nlat){
    for (j in 1:Nlon) {
      toe <- c(NA)
      for (l in 1:Ntime) {
        if (is.na(bhs.gcm[j, k, l]) == TRUE) {
          toe <- c(toe, NA)
        }
        else if (bhs.gcm[j, k, l] > ct){
          toe <- c(toe, NA)
        }
        else{
          toe<-c(toe, years[l])
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

pr.can<-"data/CanESM5_pr_SSP585_r2i1p1f1_UQM.nc"
tm.can<-"data/CanESM5_tas_SSP585_r2i1p1f1_UQM.nc"
pr.ipsl<-"data/IPSL_CM6A_LR_pr_SSP585_r1i1p1f1_UQM.nc"
tm.ipsl<-"data/IPSL_CM6A_LR_tas_SSP585_r1i1p1f1_UQM.nc"
pr.inm<-"data/INM_CM4_8_pr_SSP585_r1i1p1f1_UQM.nc"
tm.inm<-"data/INM_CM4_8_tas_SSP585_r1i1p1f1_UQM.nc"
pr.cesm<-"data/CESM2_WACCM_pr_SSP245_r3i1p1f1_UQM.nc"
tm.cesm<-"data/CESM2_WACCM_tas_SSP245_r3i1p1f1_UQM.nc"
pr.nor<-"data/NorESM2_LM_pr_SSP245_r2i1p1f1_UQM.nc"
tm.nor<-"data/NorESM2_LM_tas_SSP245_r2i1p1f1_UQM.nc"

toe.hc2.can.df <- toe.gcm.grid.fun(pr.can, tm.can, CT_P10, "P10%", m.average = 30)
toe.hc2.ipsl.df <- toe.gcm.grid.fun(pr.ipsl, tm.ipsl, CT_P10, "P10%", m.average = 30)
toe.hc2.inm.df <- toe.gcm.grid.fun(pr.inm, tm.inm, CT_P10, "P10%", m.average = 30)
toe.hc2.cesm.df <- toe.gcm.grid.fun(pr.cesm, tm.cesm, CT_P10, "P10%", m.average = 30)
toe.hc2.nor.df <- toe.gcm.grid.fun(pr.nor, tm.nor, CT_P10, "P10%", m.average = 30)

toe.hc1.can.df <- toe.gcm.grid.fun(pr.can, tm.can, CT_P50, "P50%", m.average = 30)
toe.hc1.ipsl.df <- toe.gcm.grid.fun(pr.ipsl, tm.ipsl, CT_P50, "P50%", m.average = 30)
toe.hc1.inm.df <- toe.gcm.grid.fun(pr.inm, tm.inm, CT_P50, "P50%", m.average = 30)
toe.hc1.cesm.df <- toe.gcm.grid.fun(pr.cesm, tm.cesm, CT_P50, "P50%", m.average = 30)
toe.hc1.nor.df <- toe.gcm.grid.fun(pr.nor, tm.nor, CT_P50, "P50%", m.average = 30)

toe.uni.df.30.10 <-
  rbind(
    toe.hc2.can.df,
    toe.hc2.ipsl.df,
    toe.hc2.inm.df,
    toe.hc2.cesm.df,
    toe.hc2.nor.df,
    toe.hc1.can.df,
    toe.hc1.ipsl.df,
    toe.hc1.inm.df,
    toe.hc1.cesm.df,
    toe.hc1.nor.df
  ) %>% 
  mutate(gcm = str_replace_all(gcm, '_', '-'),
         gcm = factor(gcm, levels = c("NorESM2-LM", "CESM2", "CESM2-WACCM", "INM-CM4-8", "IPSL-CM6A-LR", "CanESM5")),
         umbral = factor(umbral, levels = c("P50%", "P10%")))

map.toe.30.10 <- ggplot() +
  geom_tile(data = toe.uni.df.30.10,
            aes(x = x, y = y, fill = ToE),
            show.legend = T) +
  scale_fill_stepsn(
    colours = viridis::viridis(10),
    na.value = "white",
    breaks = seq(2010, 2100, 10)
  ) +
  geom_sf(
    data = BEMAQL,
    color = 'black',
    fill = NA,
    size = 0.1
  ) +
  scale_x_continuous(breaks = seq(-71.4, -70.2, 0.6)) +
  facet_grid(umbral ~ gcm) +
  labs(x = "", y = "") +
  theme_bw() +
  theme(
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.position = "bottom",
    legend.background = element_rect(fill = "transparent", color = NA),
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
      barwidth = 30,
      barheight = 0.5
    )
  )
map.toe.30.10
# ggsave("map_toe_MM30_hc2010.emf", map.toe.30.10, width = 8, height = 6, device = {function(filename, ...) devEMF::emf(file = filename, ...)})
# ggsave("map_toe_MM30_hc2010.png", map.toe.30.10, width = 8, height = 6, dpi = 700)
ggsave("map_toe_MM30_hc2010.pdf", map.toe.30.10, width = 8, height = 6)


