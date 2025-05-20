library(ggplot2)
library(ggmap)
library(sp)
library(maptools)
library(maps)

####LOAD WORLD MAP
world  <- map_data("world")

attach(map_sites)
sample_map <- ggplot()+
  geom_map(
    data = world, map = world,
    aes(long, lat, map_id = region),
    fill="lightgray"
  ) +
  borders("world",
          colour = "white", size=0.1) +
  ylim(-90,90)+
  theme_void() +
  geom_point(data=map_sites,
             aes(x=longitude,
                 y=latitude  ),size=2.5,shape=22,
             color="white",fill="black"
  )+
  #scale_color_manual(name="",values = habitant_cols[dis$Var1 %>% unlist() %>% as.vector()])+
  theme(legend.position = "bottom")+
  theme(aspect.ratio = 0.5)+
  xlab(NULL)+ylab(NULL)+
  coord_fixed()
detach(map_sites)
