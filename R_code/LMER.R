
#### lmer Model constructions 
model_MRO<-lmerTest::lmer(MRO_SCORE ~ Protist_Bacteria_Diversity_Ratios  + (1| Soils),  data )
summary(model_MRO);MuMIn::r.squaredGLMM(model_MRO )
model_MIP<-lmerTest::lmer(MIP_SCORE ~ Protist_Bacteria_Diversity_Ratios   +(1| Soils), data ) 
summary(model_MIP);MuMIn::r.squaredGLMM(model_MIP)
#### lmer Modelplot 
summary_model_MRO<-  summary(model_MRO)
summary_model_MRO_coeff<-summary_model_MRO$coefficients %>% as.data.frame()
model_1_plot<-ggplot(data,
                     aes(x= Protist_Bacteria_Diversity_Ratios   , y=MRO_SCORE))+
  # geom_point(size=1, alpha=0.5 , shape=16)+
  # geom_bin2d(bins = 25) +
  stat_density_2d(aes(fill = after_stat(density) ),
                  geom = "raster", contour = FALSE) +
  scale_fill_gradientn(
    colors = c('#ffffff',"#CFDFDF","#A0C0C0","#70A0A0","#408080"),
    breaks = c(0, 20, 40, 60, 80),
    labels = c("0", "20", "40", "60","80"))+
  guides(fill=guide_colourbar(title = "Density",barwidth=25,direction="horizontal",position="bottom"))+
  scale_y_continuous(expand = c(0, 0))+scale_x_continuous(expand = c(0, 0))+
  theme_bw() +
  geom_abline(intercept = as.vector(summary_model_MRO_coeff[1,1])  ,
              slope = as.vector(summary_model_MRO_coeff[2,1]), color="#FF5252",  
              linetype="solid", size=1.5)+
  ggsci::scale_color_aaas(name="")+
  theme(panel.grid=element_blank(),
        legend.position ="null") +
  theme(text = element_text(size = 20,family = "sans", color="black"),
        axis.text = element_text(size = 20,family = "sans", color="black"),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        aspect.ratio =0.5)+
  labs(x="Richness Ratio Between Protist and Bacterial Communities", 
       y="Metabolic Competition")

summary_model_MIP<-  summary(model_MIP)
summary_model_MIP_coeff<-summary_model_MIP$coefficients %>% as.data.frame()


model_2_plot<-ggplot(data,
                     aes(x= Protist_Bacteria_Diversity_Ratios   , y=MIP_SCORE))+  
  stat_density_2d(aes(fill = after_stat(density) ),
                  geom = "raster", contour = FALSE) +
  scale_fill_gradientn( 
    colors = c('#ffffff',"#CFDFDF","#A0C0C0","#70A0A0","#408080"),
    breaks = c(0, 0.1, 0.2, 0.3,0.4),
    labels = c("0", "0.1", "0.2", "0.3","0.4") )+
  guides(fill=guide_colourbar(title = "Density",barwidth=25,direction="horizontal",position="bottom"))+
  scale_y_continuous(expand = c(0, 0))+scale_x_continuous(expand = c(0, 0))+
  geom_abline(intercept = as.vector(summary_model_MIP[1,1])  ,
              slope = as.vector(summary_model_MIP[2,1]), color="#536DFE",  
              linetype="solid", size=1.5)+
  theme_bw() +
  ggsci::scale_color_aaas(name="")+
  theme(panel.grid=element_blank(),
        legend.position ="null") +
  theme(text = element_text(size = 20,family = "sans", color="black"),
        axis.text = element_text(size = 20,family = "sans", color="black"),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        aspect.ratio =0.5)+
  labs(x="Richness Ratio Between Protist and Bacterial Communities", 
       y="Metabolic Cooperation")
figure2_remend<-model_2_plot/model_1_plot

figure<-model_2_plot/model_1_plot

ggsave(filename = "Figure.pdf", figure,height = 12, width =8)



