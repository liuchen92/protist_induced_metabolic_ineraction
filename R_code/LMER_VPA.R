
model_mro_vap_analysis<-lmerTest::lmer(MRO_SCORE ~ Protist_diversity+Bacteria_diversity+Latitude +Longitude  + (1| Soils),  data )
require("glmm.hp")
glmm.hp(model_mro_vap_analysis)

model_mip_vap_analysis<-lmerTest::lmer(MIP_SCORE ~ Protist_diversity+Bacteria_diversity+Latitude +Longitude  + (1| Soils),  data )
require("glmm.hp")
glmm.hp(model_mip_vap_analysis)

