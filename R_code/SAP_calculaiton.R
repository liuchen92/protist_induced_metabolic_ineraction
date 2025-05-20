
function(ko_list_bacteria){

n=length(ko_list_bacteria)

kegg_sugar <-   ko_list_bacteria[ko_list_bacteria%in%sugar]  
kegg_acid <-  ko_list_bacteria[ko_list_bacteria%in%acid]  
sugar_rel=length(kegg_sugar)/n
acid_rel=length(kegg_acid)/n
SAP=tanh(60.76* sugar_rel - 20.21*acid_rel)
} 
