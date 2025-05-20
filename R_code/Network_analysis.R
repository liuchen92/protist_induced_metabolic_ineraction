#### R package loading
require(microeco)
require(tidyverse)
require(magrittr)
require(vegan)
require(igraph)
require(MASS )
require(furrr)
offspring.tbl_tree_item <- getFromNamespace("offspring", "tidytree")
assign("offspring.tbl_tree_item", offspring.tbl_tree_item, envir = .GlobalEnv)
#### Functions write to make network and also calcuate networ property
## Function 1: extract network and property
 network_plot_data<-function( microtable_class=network_class,
                              Treatment="treament",
                              cor_method = "spearman",
                              filter_thres = 0,
                              use_WGCNA_pearson_spearman=T,
                              COR_p_thres=0.01,
                              COR_cut = 0.7){
   
   ck <- clone(microtable_class)
   Treatment33=Treatment
   cor_method33 = cor_method
   filter_thres33 =filter_thres
   use_WGCNA_pearson_spearman33=use_WGCNA_pearson_spearman
   COR_p_thres33=COR_p_thres
   COR_cut33 = COR_cut
   ck$sample_table %<>% subset(Treatment == Treatment33)
   print(unique(ck$sample_table$Treatment))
   ck$tidy_dataset()
   ck <- trans_network$new(dataset = ck, cor_method = cor_method33, filter_thres = filter_thres33, use_WGCNA_pearson_spearman=use_WGCNA_pearson_spearman33)
   ck$cal_network(COR_p_thres = COR_p_thres33, COR_cut = COR_cut33)
   ck_net = return_networ(ck)
   edge_distribution=ck_net %>%
     activate(edges) %>%
     data.frame() %>%dplyr::select(edges) %>% table() %>% as.data.frame()
   negative_edga=edge_distribution %>% filter(edges =="negative") %>%  dplyr::select(Freq) %>% unlist() %>% as.vector()
   positive_edga=edge_distribution %>% filter(edges =="positive") %>%  dplyr::select(Freq) %>% unlist() %>% as.vector()
   coefficient=igraph::transitivity(ck_net)  
   p<-ggraph(ck_net, layout = 'kk') + 
     geom_edge_fan(aes(color=factor(edges) ),show.legend=T, alpha=0.5,width =0.5) + 
     geom_node_point(aes(fill=factor(taxonomic_annotation   ) ),size=4,shape=21)+ 
     ggsci::scale_fill_lancet()+
     scale_edge_color_manual(values  = c( "#FF5252","#388E3C"),
                             labels = c(   paste("Negative interactions:", negative_edga), 
                                           paste("Positive interactions:", positive_edga) ) )+
     guides(fill=F, 
            edge_colour=guide_legend("Correlation", keywidth = unit(1, 'cm'),
                                     position = "bottom", nrow = 2))+
     theme_graph()+theme(aspect.ratio = 1)+
     ggtitle(Treatment) +
     theme( text = element_text(size = 20,family = "sans", color="black"),
            panel.background = element_blank(),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank())
   list(ck_net, edge_distribution,p,coefficient)
 }

##Function 2 : reture network data table
 return_networ<-function(x) {
   x_net <- as_tbl_graph(x$res_network) 
   x_net=tbl_graph(edges =x_net %>%
                     activate(edges) %>%
                     data.frame() %>% mutate(edges=ifelse(label=="+", "positive", "negative")), 
                   nodes =x_net %>%
                     activate(nodes) %>%
                     data.frame()%>% 
                     mutate( taxonomic_annotation=ifelse(Phylum%in%color_for_phylum,Phylum, "other" ))  )
   
   nodes =x_net %>% activate(nodes)
   
   if("other" %in%  nodes$name ) {
     x_net=x_net %>%activate(nodes)   %>%  filter(!name%in%c("other"))
     
   }else (x_net=x_net)
   
   return(x_net)
 }
## Function 3 selecct sampels for network construction
cal_net_dege<-function(x) {
  
  seed=x
  print(seed)
  
  
  set.seed(seed)   
  g1<- sample_data  %>% filter(pro_diversity<=40) %>% dplyr::select(NewID_18S)  %>% unlist() %>% as.vector() %>% sample(25)
  g2<-  sample_data  %>% filter(pro_diversity>40,pro_diversity<=80) %>% dplyr::select(NewID_18S)  %>% unlist() %>% as.vector() %>% sample(25)
  g3<-  sample_data  %>%  filter(pro_diversity>80,pro_diversity<=120) %>% dplyr::select(NewID_18S)  %>% unlist() %>% as.vector() %>% sample(25)
  g4<-  sample_data  %>% filter(pro_diversity>120,pro_diversity<=160) %>% dplyr::select(NewID_18S)  %>% unlist() %>% as.vector() %>% sample(25)
  g5<-  sample_data  %>% filter(pro_diversity>160,pro_diversity<=200) %>% dplyr::select(NewID_18S)  %>% unlist() %>% as.vector() %>% sample(25)
  g6<-  sample_data  %>% filter(pro_diversity>200,pro_diversity<=240) %>% dplyr::select(NewID_18S)  %>% unlist() %>% as.vector() %>% sample(25)
  g7<-  sample_data  %>% filter(pro_diversity>240,pro_diversity<=280) %>% dplyr::select(NewID_18S)  %>% unlist() %>% as.vector() %>% sample(25)
  g8<-  sample_data  %>% filter(pro_diversity>280,pro_diversity<=320) %>% dplyr::select(NewID_18S)  %>% unlist() %>% as.vector() %>% sample(25)
  g9<-  sample_data  %>% filter(pro_diversity>320,pro_diversity<=360) %>% dplyr::select(NewID_18S)  %>% unlist() %>% as.vector() %>% sample(25)
  g10<-  sample_data  %>% filter(pro_diversity>360,pro_diversity<=400) %>% dplyr::select(NewID_18S)  %>% unlist() %>% as.vector() %>% sample(25)
  g11<-  sample_data  %>% filter(pro_diversity>400,pro_diversity<=440) %>% dplyr::select(NewID_18S)  %>% unlist() %>% as.vector() %>% sample(25)
  g12<-  sample_data  %>% filter(pro_diversity>440) %>% dplyr::select(NewID_18S)  %>% unlist() %>% as.vector() %>% sample(25)
  
  
   select_sample =data.frame(
    sample=c(g1,g2,g3,g4,g5,g6,g7,g8,g9,g10,g11,g12),
    Treatment=rep(c("g1","g2","g3","g4","g5","g6","g7","g8","g9","g10","g11","g12"), each=25)
    
  )
  
  
  select_sample %<>% 
    left_join(samples_information[c( "NewID_18S", "community","pro_diversity" ,"bac_diversity","pb_ratio","fungi_pro_ratio","prot_relative_abundance_to_18s","mean_mip","mean_mro" )],by=c("sample"   ="NewID_18S") )
  
  select_sample$NewID_16S=select_sample$community

  
  #### gather all samples 
  selcted_otu_table_16S <- all_table_16S[, ( select_sample$NewID_16S)  ]
  
  #### select top 150 asvs across samples selected in this simulaiton
  top_asv=rowSums(selcted_otu_table_16S) %>% sort(decreasing = T) %>% .[1:150] %>% names()
  #### all asvs not included were set as other
  new_out = selcted_otu_table_16S[top_asv,] %>% as.data.frame() %>% t() %>% as.data.frame() %>% mutate(other=rowSums(t(rhizo_otu_table_16S))- rowSums(.)) %>% t() %>% as.data.frame()

  #### taxonomy table were adjusted accordingly
  taxonomy_eco2=taxonomy_eco
  taxonomy_eco2["other",]=rep("other",6)
   #### give color for different phylum 
  color_for_phylum=  c("Proteobacteria"  ,
                       "Chloroflexi"  ,    "Actinobacteriota" ,
                       "Bacteroidota"  ,   "Acidobacteriota" ,
                       "Cyanobacteria" , "other")
  
  
  
  
  g1_plot_results<-network_plot_data(microtable_class=microtable$new(sample_table = select_sample %>%dplyr::filter(Treatment=="g1") %T>% {rownames(.)=.$NewID_16S} ,
                                                                     otu_table = new_out, 
                                                                     tax_table =taxonomy_eco2[c(top_asv,"other" ),]
  ),Treatment="g1")
  g2_plot_results<-network_plot_data(microtable_class=microtable$new(sample_table = select_sample %>%dplyr::filter(Treatment=="g2") %T>% {rownames(.)=.$NewID_16S} ,
                                                                     otu_table = new_out, 
                                                                     tax_table =taxonomy_eco2[c(top_asv,"other" ),]
  ),Treatment="g2")
  g3_plot_results<-network_plot_data(microtable_class=microtable$new(sample_table = select_sample %>%dplyr::filter(Treatment=="g3") %T>% {rownames(.)=.$NewID_16S} ,
                                                                     otu_table = new_out, 
                                                                     tax_table =taxonomy_eco2[c(top_asv,"other" ),]
  ),Treatment="g3")
  g4_plot_results<-network_plot_data(microtable_class=microtable$new(sample_table = select_sample %>%dplyr::filter(Treatment=="g4") %T>% {rownames(.)=.$NewID_16S} ,
                                                                     otu_table = new_out, 
                                                                     tax_table =taxonomy_eco2[c(top_asv,"other" ),]
  ),Treatment="g4")
  g5_plot_results<-network_plot_data(microtable_class=microtable$new(sample_table = select_sample %>%dplyr::filter(Treatment=="g5") %T>% {rownames(.)=.$NewID_16S} ,
                                                                     otu_table = new_out, 
                                                                     tax_table =taxonomy_eco2[c(top_asv,"other" ),]
  ),Treatment="g5")
  g6_plot_results<-network_plot_data(microtable_class=microtable$new(sample_table = select_sample %>%dplyr::filter(Treatment=="g6") %T>% {rownames(.)=.$NewID_16S} ,
                                                                     otu_table = new_out, 
                                                                     tax_table =taxonomy_eco2[c(top_asv,"other" ),]
  ),Treatment="g6")
  g7_plot_results<-network_plot_data(microtable_class=microtable$new(sample_table = select_sample %>%dplyr::filter(Treatment=="g7") %T>% {rownames(.)=.$NewID_16S} ,
                                                                     otu_table = new_out, 
                                                                     tax_table =taxonomy_eco2[c(top_asv,"other" ),]
  ),Treatment="g7")
  g8_plot_results<-network_plot_data(microtable_class=microtable$new(sample_table = select_sample %>%dplyr::filter(Treatment=="g8") %T>% {rownames(.)=.$NewID_16S} ,
                                                                     otu_table = new_out, 
                                                                     tax_table =taxonomy_eco2[c(top_asv,"other" ),]
  ),Treatment="g8")
  g9_plot_results<-network_plot_data(microtable_class=microtable$new(sample_table = select_sample %>%dplyr::filter(Treatment=="g9") %T>% {rownames(.)=.$NewID_16S} ,
                                                                     otu_table = new_out, 
                                                                     tax_table =taxonomy_eco2[c(top_asv,"other" ),]
  ),Treatment="g9")
  g10_plot_results<-network_plot_data(microtable_class=microtable$new(sample_table = select_sample %>%dplyr::filter(Treatment=="g10") %T>% {rownames(.)=.$NewID_16S} ,
                                                                     otu_table = new_out, 
                                                                     tax_table =taxonomy_eco2[c(top_asv,"other" ),]
  ),Treatment="g10")
  g11_plot_results<-network_plot_data(microtable_class=microtable$new(sample_table = select_sample %>%dplyr::filter(Treatment=="g11") %T>% {rownames(.)=.$NewID_16S} ,
                                                                     otu_table = new_out, 
                                                                     tax_table =taxonomy_eco2[c(top_asv,"other" ),]
  ),Treatment="g11")
  g12_plot_results<-network_plot_data(microtable_class=microtable$new(sample_table = select_sample %>%dplyr::filter(Treatment=="g12") %T>% {rownames(.)=.$NewID_16S} ,
                                                                     otu_table = new_out, 
                                                                     tax_table =taxonomy_eco2[c(top_asv,"other" ),]
  ),Treatment="g12")
  
 
  
  select_sample_sum=select_sample %>% group_by(Treatment ) %>% 
    summarise(mean_pro_diversity=mean(pro_diversity),
              mean_bac_diversity=mean(bac_diversity),
              fungi_protist_ratio_mean=mean(fungi_pro_ratio),
              prot_relative_abundance_to_18s_mean=mean(prot_relative_abundance_to_18s),
              mean_mip_mean=mean(mean_mip),
              mean_mro_mean=mean(mean_mro),
              mean_pb_diversity=mean(pb_ratio))
  
  
  select_sample_sum$postive=c(
    g1_plot_results[[2]][2,2],
    g2_plot_results[[2]][2,2],
    g3_plot_results[[2]][2,2],
    g4_plot_results[[2]][2,2],
    g5_plot_results[[2]][2,2],
    g6_plot_results[[2]][2,2],
    g7_plot_results[[2]][2,2],
    g8_plot_results[[2]][2,2],
    g9_plot_results[[2]][2,2],
    g10_plot_results[[2]][2,2],
    g11_plot_results[[2]][2,2],
    g12_plot_results[[2]][2,2]
  )
  
  select_sample_sum$negative=c(
    g1_plot_results[[2]][1,2],
    g2_plot_results[[2]][1,2],
    g3_plot_results[[2]][1,2],
    g4_plot_results[[2]][1,2],
    g5_plot_results[[2]][1,2],
    g6_plot_results[[2]][1,2],
    g7_plot_results[[2]][1,2],
    g8_plot_results[[2]][1,2],
    g9_plot_results[[2]][1,2],
    g10_plot_results[[2]][1,2],
    g11_plot_results[[2]][1,2],
    g12_plot_results[[2]][1,2]
    
  )
  
  select_sample_sum$coefficient =c(
    g1_plot_results[[4]],
    g2_plot_results[[4]],
    g3_plot_results[[4]],
    g4_plot_results[[4]],
    g5_plot_results[[4]],
    g6_plot_results[[4]],
    g7_plot_results[[4]],
    g8_plot_results[[4]],
    g9_plot_results[[4]],
    g10_plot_results[[4]],
    g11_plot_results[[4]],
    g12_plot_results[[4]]
    
    
  )
  
  select_sample_sum$all_edge=rhizo_select_sample_sum$postive+rhizo_select_sample_sum$negative
  select_sample_sum$seed=seed
  return(select_sample_sum)
}

#### make simulations 
require(furrr)
future::plan(multisession , workers = 60)
network_results_diversity_less_conserve=furrr::future_map_dfr(seq(400), function(x) {cal_net_dege(x) })
