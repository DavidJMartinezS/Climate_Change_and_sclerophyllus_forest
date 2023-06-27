# Figure 1. Study area map ----
library(tidyverse)
library(sf)
library(ggmap)
library(ggsn)
library(ggspatial)

BEMAQL <- read_sf("data/bemaql.shp")

# Study area map 
snaspe <- read_sf("data/SNASPE_junio_2022.shp") %>% # from https://www.ide.cl/
  st_transform(4326) %>%
  .[st_as_sfc(st_bbox(BEMAQL)), ] %>%
  mutate(X = st_coordinates(st_centroid(geometry))[, 1],
         Y = st_coordinates(st_centroid(geometry))[, 2])

snaspe$Nombre_tot <- str_replace(snaspe$Nombre_tot, "Parque Nacional", "N.P.")
snaspe$Nombre_tot <- str_replace(snaspe$Nombre_tot, "Reserva Nacional", "N.R.")

sbbox <-make_bbox(lon = c(-70.17764,-71.40232),
                  lat = c(-32.84243,-34.76047))
SA <- get_stamenmap(sbbox %>% unname(), zoom = 9, maptype = "terrain")

SA_ggmap <- ggmap(SA) +
  geom_sf(
    data = BEMAQL,
    fill = "slategray2",
    alpha = 0.3,
    size = .7,
    inherit.aes = FALSE
  ) +
  geom_sf(
    data = snaspe,
    fill = "green",
    alpha = 0.3,
    inherit.aes = FALSE,
    show.legend = T
  ) +
  geom_text_repel(
    data = snaspe,
    aes(X, Y, label = Unidad),
    color = "darkblue",
    size = 3.5,
    check_overlap = FALSE
  ) +
  labs(x = "", y = "") +
  theme_bw() +
  theme(rect = element_rect(fill = "transparent")) +
  scalebar(
    BEMAQL,
    st.size = 3.5,
    height = 0.012,
    dist_unit = "km",
    dist = 20,
    transform = TRUE,
    model = 'WGS84'
  ) +
  annotation_north_arrow(
    location = "tr",
    which_north = "true",
    pad_x = unit(0.2, "in"),
    pad_y = unit(1, "in"),
    style = north_arrow_fancy_orienteering
  ) +
  annotate(
    "text",
    x = -71.15,
    y = -34.80,
    label = "© OpenStreetMap contributors",
    size = 3
  )

# Secondary maps 
library(rworldxtra)
data("countriesHigh")
class(countriesHigh)
mundo<-st_as_sf(countriesHigh)

chile_reg<-raster::getData("GADM", country= "CHL", level=1)
chile_reg<-st_as_sf(chile_reg)

south_america_map <- mundo %>%
  filter(continent == "South America and the Caribbean") %>%
  ggplot() +
  geom_sf() +
  geom_sf(data = chile_reg, fill = "honeydew3") +
  geom_sf(data = st_as_sfc(st_bbox(BEMAQL)),
          fill = "red",
          alpha = 0.5) +
  scale_x_continuous(limits = c(-78,-57)) +
  scale_y_continuous(limits = c(-55,-17), 
                     breaks = seq(-50, -20, 10)) +
  theme_bw() +
  theme(plot.background = element_rect(fill = "transparent", color = NA))

map_bemaql<- ggplot() +
  geom_sf(data = chile_reg,
          fill = NA,
          inherit.aes = FALSE) +
  geom_sf(data = BEMAQL,
          fill = "honeydew3",
          inherit.aes = FALSE) +
  coord_sf(xlim = c(-71.6,-70), 
           ylim = c(-35,-32.5),
           expand = F) +
  scale_x_continuous(limits = c(-71.8, -69.8)) +
  scale_y_continuous(limits = c(-34.9, -32.7),
                     breaks = seq(-35, -33, 1)) +
  theme_bw() +
  theme(plot.background = element_rect(fill = "transparent", color = NA))

# Unit maps 
study_area <- ggplot() +
  coord_equal(xlim = c(0, 23.5),
              ylim = c(0, 26),
              expand = FALSE) +
  annotation_custom(
    ggplotGrob(south_america_map),
    xmin = 0,
    xmax = 9,
    ymin = 10,
    ymax = 26
  ) +
  annotation_custom(
    ggplotGrob(map_bemaql),
    xmin = 0,
    xmax = 9,
    ymin = 0,
    ymax = 10
  ) +
  annotation_custom(
    ggplotGrob(SA_ggmap),
    xmin = 8,
    xmax = 23.5,
    ymin = -0.5,
    ymax = 26
  ) +
  theme(plot.background = element_rect(fill = "transparent", colour = NA))

# save
#ggsave("map_SA_paper.emf",study_area, width = 8, height = 9, device = {function(filename, ...) devEMF::emf(file = filename, ...)})
#ggsave("map_SA_paper.png",study_area, width = 8, height = 9, dpi = 700)
ggsave("map_SA_paper.pdf",study_area, width = 8, height = 9)

# Figure 3. Quadrants according to historical climate ----
library(tidyverse)
library(ncdf4)
library(sf)
library(zoo)
library(reshape2)
library(RColorBrewer)
library(wesanderson)
library(SPEI)

BEMAQL <- read_sf("data/bemaql.shp")
limits <- st_bbox(BEMAQL)

pp.nc <- nc_open("data/CR2MET_pr_v2.0_mon_1979_2019_005deg.nc")
tm.nc <- nc_open("data/CR2MET_t2m_v2.0_mon_1979_2019_005deg.nc")

pp <- ncvar_get(pp.nc, "pr")
tm <- ncvar_get(tm.nc, "t2m")
lat <- ncvar_get(pp.nc, "lat")
lon <- ncvar_get(pp.nc, "lon")
time <- ncvar_get(pp.nc, "time") 

# limit to the study area
lat.sa <- lat[(which.min(abs(lat - limits[2])) - 2):(which.min(abs(lat - limits[4])) + 2)]
lon.sa <- lon[(which.min(abs(lon - limits[1])) - 2):(which.min(abs(lon - limits[3])) + 2)]
pp.sa <- pp[(which.min(abs(lon - limits[1])) - 2):(which.min(abs(lon - limits[3])) + 2), (which.min(abs(lat - limits[2])) - 2):(which.min(abs(lat - limits[4])) + 2), ]
tm.sa <- tm[(which.min(abs(lon - limits[1])) - 2):(which.min(abs(lon - limits[3])) + 2), (which.min(abs(lat - limits[2])) - 2):(which.min(abs(lat - limits[4])) + 2), ]
Nlon <- length(lon.sa)
Nlat <- length(lat.sa)
Ntime <- length(time)

# Anual Precipitation (AP) 
pp.cr2 <- matrix(NA, nrow = Nlon, ncol = Nlat)
dimnames(pp.cr2) <- list(lon.sa, lat.sa)
for (k in 1:Nlat){
  for (j in 1:Nlon) {
    pp.sm <- rollsum(pp.sa[j, k, ], 12, fill = NA, na.rm = T)
    pp.mean <- mean(na.omit(pp.sm))
    pp.cr2[j, k] <- pp.mean
  }
}
pp.cr2.df <- reshape2::melt(pp.cr2, varnames = c('x', 'y'), value.name = 'val')

pp.cr2.map <- ggplot() +
  geom_tile(data = pp.cr2.df, aes(x = x, y = y, fill = val), alpha = 0.9) +
  geom_sf(
    data = BEMAQL,
    color = "black",
    fill = NA,
    size = 0.2
  ) +
  scale_x_continuous(breaks = seq(-71.4, -70.2, 0.6), expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_gradientn(
    trans = 'log10',
    colors = c(brewer.pal(n = 8, name = "YlGnBu")),
    breaks = c(300, 500, 800, 1200, 1800)
  ) +
  labs(subtitle = "a)") +
  theme_bw() +
  theme(legend.position = "right",
        axis.title = element_blank()) +
  guides(
    fill = guide_colourbar(
      direction = 'vertical',
      title = 'AP\n(mm)',
      title.position = 'top',
      title.hjust = 0,
      ticks.colour = 'gray40',
      frame.colour = 'gray40',
      ticks.linewidth = 1,
      barheight = 10,
      barwidth = 0.5
    )
  )

# Average Annual Temperature (AAT) 
tm.cr2 <- matrix(NA, nrow = Nlon, ncol = Nlat)
dimnames(tm.cr2) <- list(lon.sa, lat.sa)
for (k in 1:Nlat){
  for (j in 1:Nlon) {
    tm.mm <- rollmean(tm.sa[j, k, ], 12, fill = NA, na.rm = T)
    tm.mean <- mean(na.omit(tm.mm))
    tm.cr2[j, k] <- tm.mean
  }
}
tm.cr2.df <- reshape2::melt(tm.cr2, varnames = c('x', 'y'), value.name = 'val')

pal <- wes_palette("Zissou1", type = "continuous")
tm.cr2.map <- ggplot() +
  geom_tile(data = tm.cr2.df, aes(x = x, y = y, fill = val), alpha = 0.8) +
  geom_sf(
    data = BEMAQL,
    color = "black",
    fill = NA,
    size = 0.2
  ) +
  scale_x_continuous(breaks = seq(-71.4,-70.2, 0.6), expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_gradientn(colours = pal) +
  labs(subtitle = "b)") +
  theme_bw() +
  theme(legend.position = "right",
        axis.title = element_blank()) +
  guides(
    fill = guide_colourbar(
      direction = 'vertical',
      title = 'AAT\n(°C)',
      title.position = 'top',
      title.hjust = 0,
      ticks.colour = 'gray40',
      frame.colour = 'gray40',
      ticks.linewidth = 1,
      barheight = 10,
      barwidth = 0.5
    )
  )

# Simplified Water Balance (SWB) 
bhs.cr2 <- matrix(NA, nrow = Nlon, ncol = Nlat)
dimnames(bhs.cr2) <- list(lon.sa, lat.sa)
for (k in 1:Nlat){
  for (j in 1:Nlon) {
    pet <- thornthwaite(tm.sa[j, k, ], lat.sa[k], na.rm =  T)
    bhs <- pp.sa[j, k, ] - pet
    bhs.sm <- rollsum(bhs, 12, fill = NA, align = "right")
    bhs.mm <- mean(bhs.sm, na.rm = T)
    bhs.cr2[j, k] <- bhs.mm
  }
}
bhs.cr2.df <- reshape2::melt(bhs.cr2, varnames = c('x', 'y'), value.name = 'val')

bhs.cr2.map <- ggplot() +
  geom_tile(data = bhs.cr2.df, aes(x = x, y = y, fill = val), alpha = 0.8) +
  geom_sf(
    data = BEMAQL,
    color = "black",
    fill = NA,
    size = 0.2
  ) +
  scale_x_continuous(breaks = seq(-71.4,-70.2, 0.6), expand =  c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_viridis_c(
    option = "viridis",
    direction = -1,
    breaks = seq(-400, 1600, 400)
  ) +
  labs(subtitle = "c)") +
  theme_bw() +
  theme(legend.position = "right",
        axis.title = element_blank()) +
  guides(
    fill = guide_colourbar(
      direction = 'vertical',
      title = 'SWB\n(mm)',
      title.position = 'top',
      title.hjust = 0,
      ticks.colour = 'gray40',
      frame.colour = 'gray40',
      ticks.linewidth = 1,
      barheight = 10,
      barwidth = 0.5
    )
  )

# Quadrants
cuadrantes <- read_sf("data/cuadrantes.shp") %>% st_transform(4326)
cuadr.labels <- data.frame(
  long = c(-71.3, -70.3, -71.3, -70.3, -71.3, -70.3),
  lati = c(-33, -33, -33.9, -33.9, -34.1, -34.1),
  label = c("Q1", "Q2", "Q3", "Q4", "C5", "C6")
)

cuadr.map <- ggplot() +
  geom_sf(
    data = BEMAQL,
    fill = 'darkolivegreen3',
    color = 'black'
  ) +
  geom_sf(data = cuadrantes,
          fill = NA,
          color = 'black') +
  scale_x_continuous(breaks = seq(-71.4,-70.2, 0.6), expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  geom_text(data = cuadr.labels,
            aes(x = long, y = lati, label = label),
            size = 4) +
  labs(subtitle = "d)") +
  theme_bw() +
  theme(axis.title = element_blank(),
        plot.margin = margin(0, 1.8, 0, 0, "cm"))

# Union 
library(patchwork)
plot_pp_tm_bhs_cuadr <- (pp.cr2.map + tm.cr2.map) / (bhs.cr2.map + cuadr.map)

#ggsave("pp_tm_bhs_cuadr.emf",cuadrantes.paper, width = 7, height = 8, device = {function(filename, ...) devEMF::emf(file = filename, ...)})
#ggsave("pp_tm_bhs_cuadr.png",cuadrantes.paper, width = 7, height = 8, dpi = 700)
ggsave("pp_tm_bhs_cuadr.pdf",cuadrantes.paper, width = 7, height = 8)

# Figure 5,6 and 7 of main article and figures S2-S13 of Supplementary data. ----
# Precipitationn, Temperature and Simplified water balance future 
library(ncdf4)
library(dplyr)
library(ggplot2)
library(zoo)
library(sp)
library(raster)
library(rgdal)
library(sf)
library(wesanderson)
library(RColorBrewer)
library(devEMF)
library(stringr)
library(rgeos)
library(cowplot)
library(SPEI)

BEMAQL <- read_sf("data/bemaql.shp")
BEMAQL_E <- BEMAQL %>%
  st_cast('POLYGON') %>%
  .[2, ]
BEMAQL_O <- BEMAQL %>%
  st_cast('POLYGON') %>%
  .[-2, ] %>%
  group_by(objectid) %>%
  summarise(across(all_of(names(BEMAQL)[-c(1, length(BEMAQL))]), first))
crs <- CRS('+proj=longlat +datum=WGS84 +no_defs')

# Plot precipitation and temperature
plot.var.fun <- function(nc) {
  nc.var <- nc_open(nc)
  if (str_detect(nc, "pr") == TRUE) {
    var <- "pp"
    var.nc <- ncvar_get(nc.var, "pr")
  }
  else {
    var <- "tm"
    var.nc <- ncvar_get(nc.var, "t2m")
  }
  
  lon <- ncvar_get(nc.var, "lon")
  lat <- ncvar_get(nc.var, "lat")
  time <- ncvar_get(nc.var, "time")
  Nlat <- length(lat)
  Nlon <- length(lon)
  Ntime <- length(time)
  volt.lat <- sort(1:Nlat, decreasing = T)
  
  for (tm in 1:Ntime) {
    var.ras<-raster(t(var.nc[,volt.lat,tm]),xmn=min(lon), xmx=max(lon), ymn=min(lat), ymx=max(lat), crs=crs)
    var.este.df<-mask(var.ras,BEMAQL_E) %>% as.data.frame(xy=T)
    var.oeste.df<-mask(var.ras,BEMAQL_O) %>% as.data.frame(xy=T)
    if (tm==1) {
      var.C1.df<-var.oeste.df[var.oeste.df$y > -33.5,3] %>% as.data.frame()
      var.C3.df<-var.oeste.df[var.oeste.df$y < -33.5 & var.oeste.df$y > -34,3] %>% as.data.frame()
      var.C5.df<-var.oeste.df[var.oeste.df$y < -34,3] %>% as.data.frame()
      var.C2.df<-var.este.df[var.este.df$y > -33.5,3] %>% as.data.frame()
      var.C4.df<-var.este.df[var.este.df$y < -33.5 & var.este.df$y > -34,3] %>% as.data.frame()
      var.C6.df<-var.este.df[var.este.df$y < -34,3] %>% as.data.frame()
    }
    else  {
      var.C1.df<-cbind(var.C1.df,var.oeste.df[var.oeste.df$y > -33.5,3]) %>% `colnames<-`(c()) 
      var.C3.df<-cbind(var.C3.df,var.oeste.df[var.oeste.df$y < -33.5 & var.oeste.df$y > -34,3]) %>% `colnames<-`(c())
      var.C5.df<-cbind(var.C5.df,var.oeste.df[var.oeste.df$y < -34,3]) %>% `colnames<-`(c())
      var.C2.df<-cbind(var.C2.df,var.este.df[var.este.df$y > -33.5,3]) %>% `colnames<-`(c())
      var.C4.df<-cbind(var.C4.df,var.este.df[var.este.df$y < -33.5 & var.este.df$y > -34,3]) %>% `colnames<-`(c())
      var.C6.df<-cbind(var.C6.df,var.este.df[var.este.df$y < -34,3]) %>% `colnames<-`(c())
    }
  }
  
  if (var=="pp") {
    c1<-apply(var.C1.df,2, mean,na.rm=T) %>% rollsum(12,fill = NA,align = 'right') %>% rollmean(10,fill = NA,align = 'right',na.rm=T)%>% as.data.frame() %>% `colnames<-`("pp") %>% cbind(fecha=seq(as.Date("1979-01-01"),as.Date("2100-12-01"),by="1 month"),punto="Q1")
    c2<-apply(var.C2.df,2, mean,na.rm=T) %>% rollsum(12,fill = NA,align = 'right') %>% rollmean(10,fill = NA,align = 'right',na.rm=T)%>% as.data.frame() %>% `colnames<-`("pp") %>% cbind(fecha=seq(as.Date("1979-01-01"),as.Date("2100-12-01"),by="1 month"),punto="Q2")
    c3<-apply(var.C3.df,2, mean,na.rm=T) %>% rollsum(12,fill = NA,align = 'right') %>% rollmean(10,fill = NA,align = 'right',na.rm=T)%>% as.data.frame() %>% `colnames<-`("pp") %>% cbind(fecha=seq(as.Date("1979-01-01"),as.Date("2100-12-01"),by="1 month"),punto="Q3")
    c4<-apply(var.C4.df,2, mean,na.rm=T) %>% rollsum(12,fill = NA,align = 'right') %>% rollmean(10,fill = NA,align = 'right',na.rm=T)%>% as.data.frame() %>% `colnames<-`("pp") %>% cbind(fecha=seq(as.Date("1979-01-01"),as.Date("2100-12-01"),by="1 month"),punto="Q4")
    c5<-apply(var.C5.df,2, mean,na.rm=T) %>% rollsum(12,fill = NA,align = 'right') %>% rollmean(10,fill = NA,align = 'right',na.rm=T)%>% as.data.frame() %>% `colnames<-`("pp") %>% cbind(fecha=seq(as.Date("1979-01-01"),as.Date("2100-12-01"),by="1 month"),punto="Q5")
    c6<-apply(var.C6.df,2, mean,na.rm=T) %>% rollsum(12,fill = NA,align = 'right') %>% rollmean(10,fill = NA,align = 'right',na.rm=T)%>% as.data.frame() %>% `colnames<-`("pp") %>% cbind(fecha=seq(as.Date("1979-01-01"),as.Date("2100-12-01"),by="1 month"),punto="Q6")
    C1<-apply(var.C1.df,2, mean,na.rm=T) %>% rollsum(12,fill = NA,align = 'right') %>% as.data.frame() %>% `colnames<-`(var) %>% cbind(fecha=seq(1,1464,1),punto="Q1")
    C2<-apply(var.C2.df,2, mean,na.rm=T) %>% rollsum(12,fill = NA,align = 'right') %>% as.data.frame() %>% `colnames<-`(var) %>% cbind(fecha=seq(1,1464,1),punto="Q2")
    C3<-apply(var.C3.df,2, mean,na.rm=T) %>% rollsum(12,fill = NA,align = 'right') %>% as.data.frame() %>% `colnames<-`(var) %>% cbind(fecha=seq(1,1464,1),punto="C3")
    C4<-apply(var.C4.df,2, mean,na.rm=T) %>% rollsum(12,fill = NA,align = 'right') %>% as.data.frame() %>% `colnames<-`(var) %>% cbind(fecha=seq(1,1464,1),punto="Q4")
    C5<-apply(var.C5.df,2, mean,na.rm=T) %>% rollsum(12,fill = NA,align = 'right') %>% as.data.frame() %>% `colnames<-`(var) %>% cbind(fecha=seq(1,1464,1),punto="Q5")
    C6<-apply(var.C6.df,2, mean,na.rm=T) %>% rollsum(12,fill = NA,align = 'right') %>% as.data.frame() %>% `colnames<-`(var) %>% cbind(fecha=seq(1,1464,1),punto="Q6")
  }
  else  {
    c1<-apply(var.C1.df,2, mean,na.rm=T) %>% rollmean(12,fill = NA,align = 'right',na.rm=T)%>% as.data.frame() %>% `colnames<-`("tm") %>% cbind(fecha=seq(as.Date("1979-01-01"),as.Date("2100-12-01"),by="1 month"),punto="Q1")
    c2<-apply(var.C2.df,2, mean,na.rm=T) %>% rollmean(12,fill = NA,align = 'right',na.rm=T)%>% as.data.frame() %>% `colnames<-`("tm") %>% cbind(fecha=seq(as.Date("1979-01-01"),as.Date("2100-12-01"),by="1 month"),punto="Q2")
    c3<-apply(var.C3.df,2, mean,na.rm=T) %>% rollmean(12,fill = NA,align = 'right',na.rm=T)%>% as.data.frame() %>% `colnames<-`("tm") %>% cbind(fecha=seq(as.Date("1979-01-01"),as.Date("2100-12-01"),by="1 month"),punto="Q3")
    c4<-apply(var.C4.df,2, mean,na.rm=T) %>% rollmean(12,fill = NA,align = 'right',na.rm=T)%>% as.data.frame() %>% `colnames<-`("tm") %>% cbind(fecha=seq(as.Date("1979-01-01"),as.Date("2100-12-01"),by="1 month"),punto="Q4")
    c5<-apply(var.C5.df,2, mean,na.rm=T) %>% rollmean(12,fill = NA,align = 'right',na.rm=T)%>% as.data.frame() %>% `colnames<-`("tm") %>% cbind(fecha=seq(as.Date("1979-01-01"),as.Date("2100-12-01"),by="1 month"),punto="Q5")
    c6<-apply(var.C6.df,2, mean,na.rm=T) %>% rollmean(12,fill = NA,align = 'right',na.rm=T)%>% as.data.frame() %>% `colnames<-`("tm") %>% cbind(fecha=seq(as.Date("1979-01-01"),as.Date("2100-12-01"),by="1 month"),punto="Q6")
    C1<-apply(var.C1.df,2, mean,na.rm=T) %>% as.data.frame() %>% `colnames<-`(var) %>% cbind(fecha=seq(1,1464,1),punto="Q1")
    C2<-apply(var.C2.df,2, mean,na.rm=T) %>% as.data.frame() %>% `colnames<-`(var) %>% cbind(fecha=seq(1,1464,1),punto="Q2")
    C3<-apply(var.C3.df,2, mean,na.rm=T) %>% as.data.frame() %>% `colnames<-`(var) %>% cbind(fecha=seq(1,1464,1),punto="Q3")
    C4<-apply(var.C4.df,2, mean,na.rm=T) %>% as.data.frame() %>% `colnames<-`(var) %>% cbind(fecha=seq(1,1464,1),punto="Q4")
    C5<-apply(var.C5.df,2, mean,na.rm=T) %>% as.data.frame() %>% `colnames<-`(var) %>% cbind(fecha=seq(1,1464,1),punto="Q5")
    C6<-apply(var.C6.df,2, mean,na.rm=T) %>% as.data.frame() %>% `colnames<-`(var) %>% cbind(fecha=seq(1,1464,1),punto="Q6")
  }
  var.list <- list(C1, C2, C3, C4, C5, C6)
  var.union <- rbind(c1, c2, c3, c4, c5, c6)
  var.union.men <- rbind(C1, C2, C3, C4, C5, C6)
  
  if (var=="pp") {
    pend <- sapply(var.list, function(x) {
      round(coefficients(lm(data = x, pp ~ fecha))[2], 5)
    }) %>% as.data.frame() %>% 
      cbind(c("Q1", "Q2", "Q3", "Q4", "Q5", "Q6")) %>% 
      `colnames <- `(c('m', 'punto'))
    titulo <- nc %>% str_extract_all("(?<=).+(?=_pr)") %>% str_replace_all("_", "-")
    plt.var <- ggplot(var.union, aes(x = fecha, y = pp))+
      geom_line()+
      geom_smooth(method = lm, se = F, colour = "red")+
      geom_text(
        data = pend,
        x = Inf,
        y = Inf,
        aes(label = paste('m= ', round(m * 120, 2), " mm/decade")),
        vjust = "inward",
        hjust = "inward")+
      facet_wrap( ~ punto, ncol = 3)+
      labs(x = "Year", y = "Anual Precipitation (mm)")+
      theme_bw()+
      theme(plot.title = element_text(hjust = 0.5),
            plot.background = element_rect(colour = NA, fill = 'transparent'),
            strip.background = element_blank(),
            strip.text = element_text(hjust =  0))
  }
  else {
    pend <-
      sapply(var.list, function(x) {
        round(coefficients(lm(data = x, tm ~ fecha))[2], 5)
      }) %>% as.data.frame() %>% cbind(c("Q1", "Q2", "Q3", "Q4", "Q5", "Q6")) %>% `colnames<-`(c('m', 'punto'))  
    titulo <- nc %>% str_extract_all("(?<=).+(?=_tas)") %>% str_replace_all("_", "-")
    plt.var <- ggplot(var.union, aes(x = fecha, y = tm)) +
      geom_line() +
      geom_smooth(method = lm, se = F, colour = "red")+
      geom_text(
        data = pend,
        x = Inf,
        y = Inf,
        aes(label = paste('m= ', round(m * 120, 2), " °C/decade")),
        vjust = "inward",
        hjust = "inward"
      ) +
      facet_wrap( ~ punto, ncol = 3)+
      labs(x = "Year", y = "Anual Average Temperature (°C)") +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5),
            plot.background = element_rect(colour = NA, fill = 'transparent'),
            strip.background = element_blank(),
            strip.text = element_text(hjust = 0)) 
  }
  
  ############ MAP ###########
  if (var == "pp") {
    pp.mn.array <- array(0, c(Nlon, Nlat, Ntime))
    
    for (i in 1:Nlat) {
      for (j in 1:Nlon) {
        sm <- rollsum(var.nc[j, i, ], 12, fill = NA, align = 'right') 
        mm<-rollmean(sm, 240, fill = NA, align = 'right', na.rm=T)
        pp.mn.array[j, i, ] <- mm
      }
    }
    pp.20.40 <- pp.mn.array[, volt.lat, 732]
    pp.ras.20.40 <- raster(t(pp.20.40), xmn=min(lon), xmx=max(lon), ymn=min(lat), ymx=max(lat), crs=crs)
    pp.ras.20.40 <- crop(pp.ras.20.40, BEMAQL, snap = 'out')
    pp.ras.20.40.mask <- mask(pp.ras.20.40, BEMAQL)
    pp.ras.20.40.df <- as.data.frame(pp.ras.20.40, xy = T) %>% cbind(periodo = "2020-2040")
    pp.ras.20.40.df.mask <- as.data.frame(pp.ras.20.40.mask, xy = T) %>% cbind(periodo = "2020-2040")
    
    pp.40.60 <- pp.mn.array[, volt.lat, 972]
    pp.ras.40.60 <- raster(t(pp.40.60), xmn=min(lon), xmx = max(lon), ymn = min(lat), ymx = max(lat), crs = crs)
    pp.ras.40.60 <- crop(pp.ras.40.60, BEMAQL, snap = 'out')
    pp.ras.40.60.mask <- mask(pp.ras.40.60, BEMAQL)
    pp.ras.40.60.df <- as.data.frame(pp.ras.40.60, xy = T) %>% cbind(periodo = "2040-2060")
    pp.ras.40.60.df.mask <- as.data.frame(pp.ras.40.60.mask, xy = T) %>% cbind(periodo = "2040-2060")
    
    pp.60.80 <- pp.mn.array[, volt.lat, 1212]
    pp.ras.60.80 <- raster(t(pp.60.80), xmn = min(lon), xmx = max(lon), ymn = min(lat), ymx = max(lat), crs = crs)
    pp.ras.60.80 <- crop(pp.ras.60.80, BEMAQL, snap = 'out')
    pp.ras.60.80.mask <- mask(pp.ras.60.80, BEMAQL)
    pp.ras.60.80.df <- as.data.frame(pp.ras.60.80, xy = T) %>% cbind(periodo = "2060-2080")
    pp.ras.60.80.df.mask <- as.data.frame(pp.ras.60.80.mask, xy = T) %>% cbind(periodo = "2060-2080")
    
    pp.80.100 <- pp.mn.array[, volt.lat, 1452]
    pp.ras.80.100 <- raster(t(pp.80.100), xmn = min(lon), xmx = max(lon), ymn = min(lat), ymx = max(lat), crs = crs)
    pp.ras.80.100 <- crop(pp.ras.80.100, BEMAQL, snap = 'out')
    pp.ras.80.100.mask <- mask(pp.ras.80.100, BEMAQL)
    pp.ras.80.100.df <- as.data.frame(pp.ras.80.100, xy = T) %>% cbind(periodo = "2080-2100")
    pp.ras.80.100.df.mask <- as.data.frame(pp.ras.80.100.mask, xy = T) %>% cbind(periodo = "2080-2100")
    
    union.pp.fut <- rbind(pp.ras.20.40.df, pp.ras.40.60.df, pp.ras.60.80.df, pp.ras.80.100.df)
    union.pp.fut.mask <- rbind(pp.ras.20.40.df.mask, pp.ras.40.60.df.mask, pp.ras.60.80.df.mask, pp.ras.80.100.df.mask)
    
    map.pp <- ggplot()+
      geom_raster(data = union.pp.fut, aes(x = x, y = y, fill = layer), alpha = 0.9)+
      geom_sf(
        data = BEMAQL,
        color = "black",
        fill = NA,
        size = 0.2
      ) +
      scale_x_continuous(breaks = seq(-72, -70, 0.5), expand = c(0, 0)) +
      scale_y_continuous(expand = c(0, 0)) +
      facet_wrap( ~ periodo, ncol = 4) +
      labs(x = "", y = "") +
      theme_bw()+
      theme(strip.background = element_blank(),
            strip.text = element_text(hjust = 0.05),
            strip.text.x = element_text(size = 12),
            plot.background = element_rect(colour = NA, fill = 'transparent'),
            panel.spacing.x = unit(1.5, "lines"),
            legend.position = "right") +
      scale_fill_gradientn(trans = 'log10', colors = c(brewer.pal(n = 8, name = "YlGnBu"))) +
      guides(
        fill = guide_colourbar(
          direction = 'vertical',
          title = 'P\n(mm)',
          title.position = 'top',
          title.hjust = 0,
          ticks.colour = 'gray40',
          frame.colour = 'gray40',
          ticks.linewidth = 2,
          barwidth = 0.5,
          barheight = 8
        ))
    
    gcmname <- nc %>% str_extract_all("(?<=).+(?=_pr)") %>% str_replace_all("_","-")
    titulo <- ggdraw() + draw_label(gcmname, fontface = 'plain')
    
    pp.plot.uni <- plot_grid(
      titulo,
      plt.var,
      map.pp,
      ncol = 1,
      labels = c('', 'a)', 'b)'),
      label_fontface = 'plain',
      rel_heights = c(0.1, 1.2, 1)
    ) +
      theme(plot.margin = margin(.2, .2, .2, .2, "cm"))
    
    #ggsave(paste0("pp_fut_",substr(nc,1,3),".emf"),pp.plot.uni, width = 7, height = 8, device = {function(filename, ...) devEMF::emf(file = filename, ...)})
    #ggsave(paste0("pp_fut_",substr(nc,1,3),".png"),pp.plot.uni, width = 7, height = 8, dpi = 700)
    ggsave(paste0("pp_fut_",substr(nc,1,3),".pdf"),pp.plot.uni, width = 7, height = 8)
    return(pp.plot.uni)
    
  } else {
    tm.mn.array <- array(0, c(Nlon, Nlat, Ntime))
    
    for (i in 1:Nlat) {
      for (j in 1:Nlon) {
        tm.mm <- rollmean(var.nc[j, i, ], 240, fill = NA, align = 'right', na.rm=T)
        tm.mn.array[j, i, ] <- tm.mm
      }
    }
    tm.20.40 <- tm.mn.array[, volt.lat, 732]
    tm.ras.20.40 <- raster(t(tm.20.40), xmn = min(lon), xmx = max(lon), ymn = min(lat), ymx = max(lat), crs = crs)
    tm.ras.20.40 <- crop(tm.ras.20.40, BEMAQL, snap = 'out')
    tm.ras.20.40.mask <- mask(tm.ras.20.40, BEMAQL)
    tm.ras.20.40.df <- as.data.frame(tm.ras.20.40, xy = T) %>% cbind(periodo = "2020-2040")
    tm.ras.20.40.mask.df <- as.data.frame(tm.ras.20.40.mask, xy = T) %>% cbind(periodo = "2020-2040")
    
    tm.40.60 <- tm.mn.array[, volt.lat, 972]
    tm.ras.40.60 <- raster(t(tm.40.60), xmn = min(lon), xmx = max(lon), ymn = min(lat), ymx = max(lat), crs = crs)
    tm.ras.40.60 <- crop(tm.ras.40.60, BEMAQL, snap = 'out')
    tm.ras.40.60.df <- as.data.frame(tm.ras.40.60, xy = T) %>% cbind(periodo = "2040-2060")
    tm.ras.40.60.mask <- mask(tm.ras.40.60, BEMAQL)
    tm.ras.40.60.mask.df <- as.data.frame(tm.ras.40.60.mask, xy = T) %>% cbind(periodo = "2040-2060")
    
    tm.60.80 <- tm.mn.array[, volt.lat, 1212]
    tm.ras.60.80 <- raster(t(tm.60.80), xmn = min(lon), xmx = max(lon), ymn = min(lat), ymx = max(lat), crs = crs)
    tm.ras.60.80 <- crop(tm.ras.60.80, BEMAQL, snap = 'out')
    tm.ras.60.80.df <- as.data.frame(tm.ras.60.80, xy = T) %>% cbind(periodo = "2060-2080")
    tm.ras.60.80.mask <- mask(tm.ras.60.80, BEMAQL)
    tm.ras.60.80.mask.df <- as.data.frame(tm.ras.60.80.mask, xy = T) %>% cbind(periodo = "2060-2080")
    
    tm.80.100 <- tm.mn.array[, volt.lat, 1452]
    tm.ras.80.100 <- raster(t(tm.80.100), xmn = min(lon), xmx = max(lon), ymn = min(lat), ymx = max(lat), crs = crs)
    tm.ras.80.100 <- crop(tm.ras.80.100, BEMAQL, snap = 'out')
    tm.ras.80.100.df <- as.data.frame(tm.ras.80.100, xy = T) %>% cbind(periodo = "2080-2100")
    tm.ras.80.100.mask <- mask(tm.ras.80.100, BEMAQL)
    tm.ras.80.100.mask.df <- as.data.frame(tm.ras.80.100.mask, xy = T) %>% cbind(periodo = "2080-2100")
    
    union.tm.fut <- rbind(tm.ras.20.40.df, tm.ras.40.60.df, tm.ras.60.80.df, tm.ras.80.100.df)
    union.tm.fut.mask <- rbind(tm.ras.20.40.mask.df, tm.ras.40.60.mask.df, tm.ras.60.80.mask.df, tm.ras.80.100.mask.df)
    pal <- wes_palette("Zissou1", type = "continuous")
    
    map.tm <- ggplot()+
      geom_raster(data = union.tm.fut, aes(x = x, y = y, fill = layer), alpha=0.9)+
      geom_sf(
        data = BEMAQL,
        color = "black",
        fill = NA,
        size = 0.2)+
      scale_x_continuous(breaks = seq(-72, -70, 0.5), expand = c(0,0))+
      scale_y_continuous(expand = c(0, 0))+
      facet_wrap( ~ periodo, ncol =  4)+
      labs(x = "", y = "") +
      theme_bw() +
      theme(
        strip.background = element_blank(),
        strip.text = element_text(hjust = 0.05),
        strip.text.x = element_text(size = 12),
        plot.background = element_rect(colour = NA, fill = 'transparent'),
        panel.spacing.x = unit(1.5, "lines"),
        legend.position = "right"
      ) +
      scale_fill_gradientn(colours = pal) +
      guides(
        fill = guide_colourbar(
          direction = 'vertical',
          title = 'AAT\n(°C)',
          title.position = 'top',
          title.hjust = 0,
          ticks.colour = 'gray40',
          frame.colour = 'gray40',
          ticks.linewidth = 1,
          barwidth = 0.5,
          barheight = 8
        )
      )
    gcmname <- nc %>% str_extract_all("(?<=).+(?=_tas)") %>% str_replace_all("_", "-")
    titulo <- ggdraw() + draw_label(gcmname, fontface = 'plain')
    
    tm.plot.uni <- plot_grid(
      titulo,
      plt.var,
      map.tm,
      ncol = 1,
      labels = c('', 'a)', 'b)'),
      label_fontface = 'plain',
      rel_heights = c(0.1, 1.2, 1)
    ) +
      theme(plot.margin = margin(.2, .2, .2, .2, "cm"))
    
    #ggsave(paste0("tm_fut_",substr(nc,1,3),".emf"),tm.plot.uni, width = 7, height = 8, device = {function(filename, ...) devEMF::emf(file = filename, ...)})
    #ggsave(paste0("tm_fut_",substr(nc,1,3),".png"),tm.plot.uni, width = 7, height = 8, dpi = 700)
    ggsave(paste0("tm_fut_",substr(nc,1,3),".pdf"),tm.plot.uni, width = 7, height = 8)
    return(tm.plot.uni)
  }
}

can.pr <- "data/CanESM5_pr_SSP585_r2i1p1f1_UQM.nc"
can.tm <- "data/CanESM5_tas_SSP585_r2i1p1f1_UQM.nc"
plt.can.pr <- plot.var.fun(can.pr)
plt.can.tm <- plot.var.fun(can.tm)
ipsl.pr <- "data/IPSL_CM6A_LR_pr_SSP585_r1i1p1f1_UQM.nc"
ipsl.tm <- "data/IPSL_CM6A_LR_tas_SSP585_r1i1p1f1_UQM.nc"
plt.ipsl.pr <- plot.var.fun(ipsl.pr)
plt.ipsl.tm <- plot.var.fun(ipsl.tm)
inm.pr <- "data/INM_CM4_8_pr_SSP585_r1i1p1f1_UQM.nc"
inm.tm <- "data/INM_CM4_8_tas_SSP585_r1i1p1f1_UQM.nc"
plt.inm.pr <- plot.var.fun(inm.pr)
plt.inm.tm <- plot.var.fun(inm.tm)
cesm.pr <- "data/CESM2_WACCM_pr_SSP245_r3i1p1f1_UQM.nc"
cesm.tm <- "data/CESM2_WACCM_tas_SSP245_r3i1p1f1_UQM.nc"
plt.cesm.pr <- plot.var.fun(cesm.pr)
plt.cesm.tm <- plot.var.fun(cesm.tm)
nor.pr <- "data/NorESM2_LM_pr_SSP245_r2i1p1f1_UQM.nc"
nor.tm <- "data/NorESM2_LM_tas_SSP245_r2i1p1f1_UQM.nc"
plt.nor.pr <- plot.var.fun(nor.pr)
plt.nor.tm <- plot.var.fun(nor.tm)

# plot simplified water balance 
bhs.gcm.plt <- function(pr.gcm, tm.gcm){
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
  volt.lat <- sort(1:Nlat, decreasing = T)
  
  bhs.gcm <- array(NA, c(Nlon, Nlat, Ntime))
  bhs.gcm.mm <- array(NA, c(Nlon, Nlat, Ntime))
  for (k in 1:Nlat){
    for (j in 1:Nlon) {
      pet <- thornthwaite(tm[j, k, ], lat[k], na.rm = T)
      bhs <- pr[j, k, ] - pet
      bhs.sm <- rollsum(bhs, 12, fill = NA, align = "right")
      bhs.mm <- rollmean(bhs.sm, 240, fill = NA, align = "right")
      bhs.gcm[j, k, ] <- bhs.sm
      bhs.gcm.mm[j, k, ] <- bhs.mm
    }
  }
  
  for (j in 1:Ntime) {
    var.ras <- raster(t(bhs.gcm[, volt.lat, j]), xmn = min(lon), xmx = max(lon), ymn = min(lat), ymx = max(lat), crs = crs)
    var.este.df <- mask(var.ras, BEMAQL_E) %>% as.data.frame(xy = T)
    var.oeste.df <- mask(var.ras, BEMAQL_O) %>% as.data.frame(xy = T)
    if (j==1) {
      var.C1.df <- var.oeste.df[var.oeste.df$y > -33.5, 3] %>% as.data.frame()
      var.C3.df <- var.oeste.df[var.oeste.df$y < -33.5 & var.oeste.df$y > -34, 3] %>% as.data.frame()
      var.C5.df <- var.oeste.df[var.oeste.df$y < -34, 3] %>% as.data.frame()
      var.C2.df <- var.este.df[var.este.df$y > -33.5, 3] %>% as.data.frame()
      var.C4.df <- var.este.df[var.este.df$y < -33.5 & var.este.df$y > -34, 3] %>% as.data.frame()
      var.C6.df <- var.este.df[var.este.df$y < -34, 3] %>% as.data.frame()
    }
    else  {
      var.C1.df <- cbind(var.C1.df, var.oeste.df[var.oeste.df$y > -33.5, 3]) %>% `colnames<-`(c()) 
      var.C3.df <- cbind(var.C3.df, var.oeste.df[var.oeste.df$y < -33.5 & var.oeste.df$y > -34, 3]) %>%  `colnames<-`(c())
      var.C5.df <- cbind(var.C5.df, var.oeste.df[var.oeste.df$y < -34, 3]) %>%  `colnames<-`(c())
      var.C2.df <- cbind(var.C2.df, var.este.df[var.este.df$y > -33.5, 3]) %>% `colnames<-`(c())
      var.C4.df <- cbind(var.C4.df, var.este.df[var.este.df$y < -33.5 & var.este.df$y > -34, 3]) %>%  `colnames<-`(c())
      var.C6.df <- cbind(var.C6.df, var.este.df[var.este.df$y < -34, 3]) %>% `colnames<-`(c())
    }
  }  
  c1 <- apply(var.C1.df, 2, mean, na.rm = T) %>% as.data.frame() %>% `colnames<-`("bhs") %>% cbind(fecha = seq(as.Date("1979-01-01"), as.Date("2100-12-01"), by = "1 month"), date = seq(1, 1464, 1), punto = "Q1")
  c2 <- apply(var.C2.df, 2, mean, na.rm = T) %>% as.data.frame() %>% `colnames<-`("bhs") %>% cbind(fecha = seq(as.Date("1979-01-01"), as.Date("2100-12-01"), by = "1 month"), date = seq(1, 1464, 1), punto = "Q2")
  c3 <- apply(var.C3.df, 2, mean, na.rm = T) %>% as.data.frame() %>% `colnames<-`("bhs") %>% cbind(fecha = seq(as.Date("1979-01-01"), as.Date("2100-12-01"), by = "1 month"), date = seq(1, 1464, 1), punto = "Q3")
  c4 <- apply(var.C4.df, 2, mean, na.rm = T) %>% as.data.frame() %>% `colnames<-`("bhs") %>% cbind(fecha = seq(as.Date("1979-01-01"), as.Date("2100-12-01"), by = "1 month"), date = seq(1, 1464, 1), punto = "Q4")
  c5 <- apply(var.C5.df, 2, mean, na.rm = T) %>% as.data.frame() %>% `colnames<-`("bhs") %>% cbind(fecha = seq(as.Date("1979-01-01"), as.Date("2100-12-01"), by = "1 month"), date = seq(1, 1464, 1), punto = "Q5")
  c6 <- apply(var.C6.df, 2, mean, na.rm = T) %>% as.data.frame() %>% `colnames<-`("bhs") %>% cbind(fecha = seq(as.Date("1979-01-01"), as.Date("2100-12-01"), by = "1 month"), date = seq(1, 1464, 1), punto = "Q6")
  
  var.list <- list(c1, c2, c3, c4, c5, c6)
  var.union <- rbind(c1, c2, c3, c4, c5, c6)
  
  pend <- sapply(var.list, function(x) {
    round(coefficients(lm(data = x, bhs ~ date))[2], 5)
  }) %>% as.data.frame() %>% cbind(c("Q1", "Q2", "Q3", "Q4", "Q5", "Q6")) %>% `colnames<-`(c('m', 'punto'))
  plt.bhs <- ggplot(var.union, aes(x = fecha, y = bhs))+
    geom_line() +
    geom_smooth(method = lm,
                se = F,
                colour = "red") +
    geom_text(
      data = pend,
      x = Inf,
      y = Inf,
      aes(label = paste('m= ', round(m * 120, 2), " mm/decade")),
      vjust = "inward",
      hjust = "inward"
    ) +
    facet_wrap(~ punto, ncol = 3) +
    labs(x = "Year", y = "SWB (mm)") +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5),
      plot.background = element_rect(fill = 'transparent', colour = NA),
      strip.background = element_blank(),
      strip.text = element_text(hjust = 0)
    )
  
  ########## MAP ########
  bhs.20.40 <- bhs.gcm.mm[, volt.lat, 732]
  bhs.ras.20.40 <- raster(t(bhs.20.40), xmn = min(lon), xmx = max(lon), ymn = min(lat), ymx = max(lat), crs = crs)
  bhs.ras.20.40 <- crop(bhs.ras.20.40, BEMAQL, snap = 'out')
  bhs.ras.20.40.mask <- mask(bhs.ras.20.40, BEMAQL)
  bhs.ras.20.40.df <- as.data.frame(bhs.ras.20.40, xy = T) %>% cbind(periodo = "2020-2040")
  bhs.ras.20.40.df.mask <- as.data.frame(bhs.ras.20.40.mask, xy = T) %>% cbind(periodo = "2020-2040")
  
  bhs.40.60 <- bhs.gcm.mm[, volt.lat, 972]
  bhs.ras.40.60 <- raster(t(bhs.40.60), xmn = min(lon), xmx = max(lon), ymn = min(lat), ymx = max(lat), crs = crs)
  bhs.ras.40.60 <- crop(bhs.ras.40.60, BEMAQL, snap = 'out')
  bhs.ras.40.60.mask <- mask(bhs.ras.40.60, BEMAQL)
  bhs.ras.40.60.df <- as.data.frame(bhs.ras.40.60, xy = T) %>% cbind(periodo = "2040-2060")
  bhs.ras.40.60.df.mask <- as.data.frame(bhs.ras.40.60.mask, xy = T) %>% cbind(periodo = "2040-2060")
  
  bhs.60.80 <- bhs.gcm.mm[, volt.lat, 1212]
  bhs.ras.60.80 <- raster(t(bhs.60.80), xmn = min(lon), xmx = max(lon), ymn = min(lat), ymx = max(lat), crs = crs)
  bhs.ras.60.80 <- crop(bhs.ras.60.80, BEMAQL, snap = 'out')
  bhs.ras.60.80.mask <- mask(bhs.ras.60.80, BEMAQL)
  bhs.ras.60.80.df <- as.data.frame(bhs.ras.60.80, xy = T) %>% cbind(periodo = "2060-2080")
  bhs.ras.60.80.df.mask <- as.data.frame(bhs.ras.60.80.mask, xy = T) %>% cbind(periodo = "2060-2080")
  
  bhs.80.100 <- bhs.gcm.mm[, volt.lat, 1452]
  bhs.ras.80.100 <- raster(t(bhs.80.100), xmn = min(lon), xmx = max(lon), ymn = min(lat), ymx = max(lat), crs = crs)
  bhs.ras.80.100 <- crop(bhs.ras.80.100, BEMAQL, snap = 'out')
  bhs.ras.80.100.mask <- mask(bhs.ras.80.100, BEMAQL)
  bhs.ras.80.100.df <- as.data.frame(bhs.ras.80.100, xy = T) %>% cbind(periodo = "2080-2100")
  bhs.ras.80.100.df.mask <- as.data.frame(bhs.ras.80.100.mask, xy = T) %>% cbind(periodo = "2080-2100")
  
  union.bhs.fut <- rbind(bhs.ras.20.40.df, bhs.ras.40.60.df, bhs.ras.60.80.df, bhs.ras.80.100.df)
  pal <- wes_palette("Zissou1", type = "continuous")
  
  map.bhs <- ggplot() +
    geom_raster(data = union.bhs.fut, aes(x = x, y = y, fill = layer), alpha = 0.9) +
    geom_sf(
      data = BEMAQL,
      color = "black",
      fill = NA,
      size = 0.2
    ) +
    scale_x_continuous(breaks = seq(-72, -70, 0.5), expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    facet_wrap( ~ periodo, ncol = 4)+
    labs(x="",y="")+
    theme_bw() +
    theme(
      strip.background = element_blank(),
      strip.text = element_text(hjust = 0.05),
      strip.text.x = element_text(size = 12),
      plot.background = element_rect(colour = NA, fill = 'transparent'),
      panel.spacing.x = unit(1.5, "lines"),
      legend.position = "right"
    ) +
    scico::scale_fill_scico(palette = "roma") +
    guides(
      fill = guide_colourbar(
        direction = 'vertical',
        title = 'SWB\n(mm)',
        title.position = 'top',
        title.hjust = 0,
        ticks.colour = 'gray40',
        frame.colour = 'gray40',
        ticks.linewidth = 1,
        barwidth = 0.5,
        barheight = 8
      )
    )
  ######## UNION DE PLOTS #######
  gcmname <- pr.gcm %>% str_extract_all("(?<=).+(?=_pr)") %>% str_replace_all("_","-")
  titulo <- ggdraw() + draw_label(gcmname, fontface = 'plain')
  
  bhs.plot.uni <- plot_grid(
    titulo,
    plt.bhs,
    map.bhs,
    ncol = 1,
    labels = c('', 'a)', 'b)'),
    label_fontface = 'plain',
    rel_heights = c(0.1, 1.2, 1)
  ) +
    theme(plot.margin = margin(.2, .2, .2, .2, "cm"))
  
  #ggsave(paste0("bhs_fut_",substr(pr.gcm,1,3),".emf"),bhs.plot.uni, width = 7, height = 8, device = {function(filename, ...) devEMF::emf(file = filename, ...)})
  #ggsave(paste0("bhs_fut_",substr(pr.gcm,1,3),".png"),bhs.plot.uni, width = 7, height = 8, dpi = 700)
  ggsave(paste0("bhs_fut_",substr(pr.gcm,1,3),".pdf"),bhs.plot.uni, width = 7, height = 8)
  return(bhs.plot.uni)
}

can.pr <- "data/CanESM5_pr_SSP585_r2i1p1f1_UQM.nc"
can.tm <- "data/CanESM5_tas_SSP585_r2i1p1f1_UQM.nc"
plt.can.bhs <- bhs.gcm.plt(can.pr, can.tm)
ipsl.pr <- "data/IPSL_CM6A_LR_pr_SSP585_r1i1p1f1_UQM.nc"
ipsl.tm <- "data/IPSL_CM6A_LR_tas_SSP585_r1i1p1f1_UQM.nc"
plt.ipsl.bhs <- bhs.gcm.plt(ipsl.pr, ipsl.tm)
inm.pr <- "data/INM_CM4_8_pr_SSP585_r1i1p1f1_UQM.nc"
inm.tm <- "data/INM_CM4_8_tas_SSP585_r1i1p1f1_UQM.nc"
plt.inm.bhs <- bhs.gcm.plt(inm.pr, inm.tm)
cesm.pr <- "data/CESM2_WACCM_pr_SSP245_r3i1p1f1_UQM.nc"
cesm.tm <- "data/CESM2_WACCM_tas_SSP245_r3i1p1f1_UQM.nc"
plt.cesm.bhs <- bhs.gcm.plt(cesm.pr, cesm.tm)
nor.pr <- "data/NorESM2_LM_pr_SSP245_r2i1p1f1_UQM.nc"
nor.tm <- "data/NorESM2_LM_tas_SSP245_r2i1p1f1_UQM.nc"
plt.nor.bhs <- bhs.gcm.plt(nor.pr, nor.tm)

# Figure 8. critical thresholds ----
library(ncdf)
library(tidyverse)
library(sf)
library(scico)

load("data/um_cri_grid_10.RData")
load("data/um_cri_grid_25.RData")

BEMAQL <- read_sf("data/bemaql.shp")

pr.can <- "data/CanESM5_pr_SSP585_r2i1p1f1_UQM.nc"
nc.pr <- nc_open(pr.can)
lon <- ncvar_get(nc.pr, "lon")
lat <- ncvar_get(nc.pr, "lat")

p10.df <- raster::raster(t(um.cri.10[, 43:1]), xmn = min(lon), xmx = max(lon), ymn = min(lat), ymx = max(lat), crs = crs) %>% 
  raster::crop(BEMAQL, snap = 'out') %>% 
  raster::as.data.frame(xy = T) %>% 
  mutate(umbral = "P10%")
p25.df <- raster::raster(t(um.cri.25[, 43:1]), xmn = min(lon), xmx = max(lon), ymn = min(lat), ymx = max(lat), crs = crs) %>% 
  raster::crop(BEMAQL, snap = 'out') %>% 
  raster::as.data.frame(xy = T) %>% 
  mutate(umbral = "P25%")
um.uni <- bind_rows(p25.df, p10.df) %>% mutate(umbral = factor(umbral, levels = c("P25%","P10%")))

um.cri.map <- ggplot() +
  geom_tile(data = um.uni, aes(x = x, y = y, fill = layer)) +
  geom_sf(data = BEMAQL,
          fill = NA,
          color = "black") +
  scale_x_continuous(breaks = seq(-71.4, -70.2, 0.6), expand = c(0, 0)) +
  scico::scale_fill_scico(palette = "roma",
                          breaks = seq(-1000, 1000, 200),) +
  facet_wrap( ~ umbral) +
  labs(x = "", y = "") +
  theme_bw() +
  theme(
    plot.background = element_rect(fill = "transparent", color = NA),
    strip.background = element_blank(),
    panel.spacing.x = unit(1.2, "lines"),
    legend.position = "bottom",
    strip.text = element_text(size = 12)
  ) +
  scale_x_continuous(expand = c(0, 0), breaks = seq(-72, -70, 0.5)) +
  scale_y_continuous(expand = c(0, 0)) +
  guides(
    fill = guide_colourbar(
      direction = 'horizontal',
      title = 'SWB threshold (mm)',
      title.position = 'top',
      title.hjust = 0.5,
      ticks.colour = 'gray40',
      frame.colour = 'gray40',
      ticks.linewidth = 1.5,
      barwidth = 20,
      barheight = 0.7
    )
  )

#ggsave("um_cri_map.emf", um.cri.map, width = 6, height = 6, device = {function(filename, ...) devEMF::emf(file = filename, ...)})
#ggsave("um_cri_map.png", um.cri.map, width = 5.5, height = 6, dpi = 700)
ggsave("um_cri_map.pdf", um.cri.map, width = 5.5, height = 6)

pr.can<-"data/CanESM5_pr_SSP585_r2i1p1f1_UQM.nc"
nc.pr<-nc_open(pr.can)
lon<-ncvar_get(nc.pr,"lon")
lat<-ncvar_get(nc.pr,"lat")

um.cri.map<-ggplot()+
  geom_tile(data=um.uni,aes(x=x,y=y,fill=layer))+
  geom_sf(data = BEMAQL,fill=NA,color="black")+
  scale_x_continuous(breaks = seq(-71.4,-70.2,0.6),expand = c(0,0))+
  scico::scale_fill_scico(palette = "roma",
                          breaks=seq(-1000,1000,200)
  )+
  facet_wrap(~umbral)+
  labs(x="",y="")+
  theme_bw()+
  theme(plot.background = element_rect(fill = "transparent", color = NA),
        strip.background = element_blank(),
        panel.spacing.x = unit(1.2, "lines"),
        legend.position = "bottom",
        strip.text = element_text(size=12))+
  scale_x_continuous(expand = c(0,0),breaks = seq(-72,-70,0.5))+
  scale_y_continuous(expand = c(0,0))+
  guides(fill = guide_colourbar(direction = 'horizontal',
                                title='SWB threshold (mm)',  
                                title.position='top',
                                title.hjust=0.5,
                                ticks.colour='gray40',frame.colour = 'gray40',
                                ticks.linewidth=1.5,
                                barwidth = 20,
                                barheight = 0.7))

# Figure S1. Points second methodology ----
library(OpenStreetMap)
library(cowplot)
library(sf)
library(tidyverse)
library(magick)

BEMAQL <- read_sf("data/bemaql.shp")
puntos <- read.csv("C:/Users/Lenovo/Documents/Trabajo con profe/puntos_metodo_2.csv", header = T, sep = ";")
limits <- st_bbox(BEMAQL)

map<-openmap(c(unname(limits)[c(2,1)]),c(unname(limits)[c(4,3)]),minNumTiles=9,type = "bing")
plot(map)
map2 <- openproj(map)

map2_plt <- autoplot.OpenStreetMap(map2)
map2_plt <- map2_plt +
  geom_sf(
    data = BEMAQL,
    fill = NA,
    color = 'black',
    size = 0.7,
    inherit.aes = F
  ) +
  geom_point(
    data = puntos,
    aes(x = Lon, y = Lat),
    color = 'darkturquoise',
    size = 2.5
  ) +
  theme(axis.title = element_blank()) +
  annotate(
    "text",
    x = -70.50,
    y = -34.70,
    label = "? OpenStreetMap contributors",
    size = 3,
    color = 'white'
  )

C12019 <- magick::image_read("data/C1_2019.png")
C62019 <- magick::image_read("data/C6_2019.png")
C12019.plt <- ggdraw() + draw_image(C12019)
C62019.plt <- ggdraw() + draw_image(C62019)

pts.y.GEP <-
  plot_grid(
    map2_plt,
    plot_grid(
      C12019.plt,
      C62019.plt,
      labels = c('b)', 'c)'),
      label_fontface = "plain",
      nrow = 2
    ),
    ncol = 2,
    labels = "a)",
    label_fontface = "plain")

# ggsave("map_ptos_image.emf",pts.y.GEP, width = 11, height = 8, device = {function(filename, ...) devEMF::emf(file = filename, ...)})
# ggsave("map_ptos_image.png",pts.y.GEP, width = 11, height = 8, dpi = 700)
ggsave("map_ptos_image.pdf",pts.y.GEP, width = 11, height = 8)

# Figure S15. Vegetation unit of Luebert and Pliscoff (2017) ----
library(tidyverse)
library(sf)
library(scales)
library(RColorBrewer)
library(ggsci)
library(wesanderson)
library(ggthemes)
library(cowplot)
library(ggsn)
library(ggspatial)
library(raster)
library(rworldxtra)

pisos_sf <- read_sf("data/PisosVegetacionalesPliscoff2017.shp") %>% st_transform(4326)

#getdata Chile
chile_reg <- raster::getData("GADM", country = "CHL", level = 1) %>% st_as_sf()
sa <- c("Valparaíso", "Región Metropolitana de Santiago", "Libertador General Bernardo O'Higgins")

reg_sa <- chile_reg %>%
  st_cast("POLYGON") %>% 
  mutate(
    X = st_coordinates(st_centroid(geometry))[, 1],
    Y = st_coordinates(st_centroid(geometry))[, 2]
  ) %>%
  filter(NAME_1 %in% sa & X > -75)

reg_sa2 <- chile_reg %>%
  st_cast("POLYGON") %>% 
  mutate(
    X = st_coordinates(st_centroid(geometry))[, 1],
    Y = st_coordinates(st_centroid(geometry))[, 2]
  ) %>%
  filter(CC_1 %in% c("IV", "V", "RM", "VI", "VII") & X > -75)

sf::sf_use_s2(F)
pisos_sa <- st_crop(pisos_sf, reg_sa)

pal <- stata_pal("s2color")(15)[-c(6, 7)]
palettes <- ggthemes_data[["tableau"]][["color-palettes"]][["regular"]]
pal2<-pull(palettes$`Miller Stone`,value)[1:9]

pisos_sa <- pisos_sa %>% 
  filter(formacion %in% c("Bosque esclerofilo","Bosque caducifolio")) %>% 
  mutate(X=st_coordinates(st_centroid(geometry))[,1],
         Y=st_coordinates(st_centroid(geometry))[,2]) %>% 
  filter(!row_number() %in% c(10)) # remove bosque caducifolio mediterraneo-templado 

sp <- pisos_sa$piso %>%
  str_replace_all("/", "-") %>%
  str_replace_all("-", "and") %>%
  str_extract_all("(?<=de  ?).+") %>%
  unlist()
spp <- sp %>% str_split(" and ", n = 2, simplify = TRUE)
sp1 <- spp[, 1]
sp2 <- spp[, 2]

veget.español<-c("Bosque esclerofilo mediterráneo andino de ",
                 "Bosque esclerofilo mediterráneo costero de ",
                 "Bosque esclerofilo mediterráneo interior de ",
                 "Bosque caducifolio mediterráneo costero de ",
                 "Bosque caducifolio mediterráneo interior de ")
veget.ingles<-c("Andean mediterranean sclerophyllous forest of ",
                "Coastal mediterranean sclerophyllous forest of ",
                "Interior mediterranean sclerophyllous forest of ",
                "Coastal mediterranean deciduous forest of ",
                "Interior mediterranean deciduous forest of ")

pisos_sa$PV <- ''
for (i in 1:nrow(pisos_sa)) {
  for (j in 1:length(veget.ingles)) {
    if (str_detect(pisos_sa$piso[i],veget.español[j])==TRUE) {
      pisos_sa$PV[i]<-paste0(veget.ingles[j],sp[i])
    }
  }
}

labels.pv <- pisos_sa$PV %>% 
  str_extract_all("(?<=).+(?<= of)") %>%
  str_c("\n", sp)

data("countriesHigh")
class(countriesHigh)
mundo <- st_as_sf(countriesHigh)
limites.reg <- st_bbox(reg_sa)
south_america_map <- mundo %>% 
  filter(continent == "South America and the Caribbean") %>% 
  ggplot() +
  geom_sf()+
  geom_sf(data=chile_reg, fill="honeydew3")+
  geom_sf(data = st_as_sfc(limites.reg), fill = "red",alpha=0.5) +
  scale_x_continuous(limits = c(-78, -64),breaks = seq(-76,-64,4))+
  scale_y_continuous(limits = c(-55, -17),breaks = seq(-50,-20,10))+
  theme_bw()+
  theme(plot.background = element_rect(fill = "transparent", color = NA),
        axis.ticks = element_blank(),
        axis.text = element_blank())

plot.pisos<-ggplot()+
  geom_sf(data=reg_sa2,color='black',fill="antiquewhite1",alpha=0.5)+
  geom_sf(data=pisos_sa,aes(fill=PV))+
  scale_fill_manual(values = pal2, labels=labels.pv)+
  labs(fill="Vegetation Units\n  ")+
  coord_sf(xlim = c(-72.2,-68.3),ylim = c(-34.88,-32.2))+
  theme_bw()+
  theme(legend.spacing = unit(4, "cm"),
        legend.background = element_rect(color = 'black'),
        legend.text = element_text(size = 9),
        legend.title = element_text(size = 10,face = "bold"),
        legend.key.height = unit(1.1,"cm"),
        plot.background = element_blank(),
        panel.grid.major = element_line(color = "gray70", linetype = "dashed", size = 0.5),
        panel.background = element_rect(color = 'black',fill = "aliceblue"),
        legend.position = c(0.83,0.67)
  )+
  guides(fill = guide_legend(byrow = TRUE))+
  annotation_scale(location = "br", width_hint = 0.4,style = "ticks",height = unit(0.5, "cm"),line_width = 2,pad_x = unit(1, "cm"))+
  annotation_north_arrow(location = "br", which_north = "true", 
                         height = unit(2, "cm"), width = unit(2, "cm"), 
                         pad_x = unit(3, "cm"), pad_y = unit(2, "cm"),
                         style = north_arrow_fancy_orienteering(text_size = unit(14,"cm")))
map.LP <- ggdraw() +
  draw_plot(plot.pisos) +
  draw_plot(
    plot = south_america_map,
    x = -0.13,
    y = 0.48,
    height = .495,
    width = .495
  )
# ggsave("map_LP_R.emf",map.LP, width = 10, height = 8.3, device = {function(filename, ...) devEMF::emf(file = filename, ...)})
# ggsave("map_LP_R.png",map.LP, width = 10, height = 8.3, dpi = 700)
ggsave("map_LP_R.pdf",map.LP, width = 10, height = 8.3)






