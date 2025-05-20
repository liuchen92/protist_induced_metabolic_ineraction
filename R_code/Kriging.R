kringh_input_df<-  data %>% 
  group_by(Latitude,Longitude ) %>% 
  summarise(mean_protist_diversity=mean(Protist_diversity),
            mean_bacdiversity=mean(Bacteria_diversity),
            mean_bp=mean(Protist_Bacteria_Diversity_Ratios  ),
            mean_mip_result=mean(MIP_SCORE ),
            mean_mro_result=mean(MRO_SCORE ))



sp::coordinates(kringh_input_df) = ~Longitude+Latitude
sp::proj4string(kringh_input_df) <- CRS("+proj=longlat +datum=WGS84")

#### gridding global map 
world_sf <- ne_countries(scale = 110, type = "countries", continent = NULL,
                         country = NULL, geounit = NULL, sovereignty = NULL,
                         returnclass = c("sf"))
world_sf <- world_sf %>% filter(sov_a3 != 'ATA')
grid_global <- expand.grid(x=seq(from =  sf::st_bbox(world_sf)[1],to =  sf::st_bbox(world_sf)[3],length.out = 800),
                    y=seq(from =  sf::st_bbox(world_sf)[2],to =  sf::st_bbox(world_sf)[4],length.out = 400))


sp::coordinates(grid_global) <- ~x+y

sp::proj4string(grid_global) <- CRS("+proj=longlat +datum=WGS84")
sp::gridded(grid_global) <- TRUE


#####variogram
semivariog_mro_div<-gstat::variogram(mean_mro_result  ~1, locations=kringh_input_df, data=kringh_input_df)
plot(semivariog_mro_div)
model.variog_mro<-gstat::vgm(psill=1.5e-03, model=c( "Mat"),  ### compare "Mat","Exp", "Sph", "Gau",  select MAT model
                             nugget=1e-04, range=4000) #### based on plot(semivariog_mro_div)
fit.variog_mro<-gstat::fit.variogram(semivariog_mro_div, model.variog_mro)
plot(semivariog_mro_div, model.variog_mro)
#### Cross validaiton
set.seed(100010)
kriging_cv<- krige.cv(formula = mean_mro_result  ~  1, 
                      locations=kringh_input_df,  
                      model=model.variog_mro ,nfold= 10)

 
kriging_cv_df<-as.data.frame(kriging_cv)
lm(observed~var1.pred, data =kriging_cv_df  ) %>% summary()

result <- cor.test(kriging_cv_df$observed, kriging_cv_df$var1.pre, method = "pearson")

print(result )


mro_krige_cv<-ggplot(kriging_cv_df  ,
       aes(x= observed   , y= var1.pred))+
  geom_point(size=4, alpha=1 , shape=22, color="white", fill="#795548")+
  stat_smooth(
    alpha=0.2 , method="lm",se=T, size=1, color="#536DFE",fill="#536DFE",
    formula =y~ poly(x,1) ,show.legend = FALSE)+
  ggpmisc::stat_poly_eq(size=5, family = "sans",color="black",
                        aes(label =  paste(..rr.label..,..p.value.label..,sep = "~~~~") ), 
                        formula = y ~ poly(x,1))+
  theme_bw() +
  ggsci::scale_color_aaas(name="")+
  theme(panel.grid=element_blank(),
        legend.position ="null") +
  theme(text = element_text(size = 20,family = "sans", color="black"),
        axis.text = element_text(size = 20,family = "sans", color="black"),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        aspect.ratio =1)+
  labs(y="Predicted MRO Scores", 
       x="Computed MRO Scores")+ theme (
         panel.background = element_rect(fill = alpha("grey",0.2)  ,
                                         colour =alpha("grey",0.2) )

    
#### predation and ploting
set.seed(100010)
krig_mro<-krige(formula=mean_mro_result  ~  1,  locations=kringh_input_df,  newdata=grid_global, model=model.variog_mro)
raster_result_mro <-raster::raster(krig_mro)
land_mask_mro <- rasterize(world_sf, raster_result_mro, field = 1, background = NA)
raster_mro_masked <- mask(raster_result_mro, land_mask_mro)
raster_mro_df <- as.data.frame(raster_mro_masked, xy = TRUE)
colnames(raster_mro_df) <- c("lon", "lat", "value")
raster_mro_df <- na.omit(raster_mro_df)


mro_krige_map<-ggplot(data = raster_mro_df, aes(x = lon, y = lat, fill = value)) +
  geom_raster() +
  scale_fill_gradientn(colours = my_colormap,  
                       guide = guide_colorbar(barwidth = 0.5, barheight = 1, 
                                              ticks.colour = 'black', 
                                              ticks.linewidth = 1,     
                                              frame.colour = "black") )+
  # spatial-aware automagic north arrow
  labs(
    title = "MRO score"
  ) +
  #theme_ipsum(base_family = "Roboto Condensed") +
  theme(aspect.ratio = 0.5,
        plot.title = element_text(hjust = 0.5))+
  coord_fixed()+ 
  guides(fill = guide_colorbar("MRO score",title.position = "bottom",
                               title.hjust = 0.5,
                               barwidth =unit(1, "null"),
                               barheight = unit(0.5, "cm"))) + 
  theme(legend.position = "bottom")+
  theme(text = element_text(size = 20,family = "sans", color="black"),
        axis.text = element_text(size = 20,family = "sans", color="black"))+ theme (
          panel.background = element_rect(fill = alpha("grey",0.0)  ,
                                          colour =alpha("grey",0.0) ),
          panel.grid.major=element_line(color='gray',linetype="dashed")
        )

       )
