#https://github.com/fbreitwieser/sankeyD3
#devtools::install_github("fbreitwieser/sankeyD3")
library(sankeyD3)
library(dplyr)



nodes = data.frame(name = c("init_PE","init_PP","init_EE", "multi_PE","multi_PP","multi_EE", paste0("PE_",1:10), 
                            paste0("PP_",1:4), paste0("EE_",1:15), 
                            paste0("combined_",1:29)), 
                   group=as.character(c("blue","red","green",
                                        "blue","red","green",
                                        rep("blue",10),rep("red",4),rep("green",15),
                                        rep("blue",10),rep("red",4),rep("green",15))))



#init pp pe ee and multi pp pe ee
left=c(rep(1,3617),rep(2,800),rep(3,5267),rep(4,3406),rep(5,687),rep(6,5462))-1



pe.init=read.table(paste0("/fast/AG_Ohler/henriette/PANCREAS_final/output_hg19/open_regions_subset/open_regions_pairs/timeless_init_prom-enh/2_30/classes-10.txt"))
right.pe.init=pe.init$V1+6-1

pp.init=read.table(paste0("/fast/AG_Ohler/henriette/PANCREAS_final/output_hg19/open_regions_subset/open_regions_pairs/timeless_init_prom-prom/2_30/classes-4.txt"))
right.pp.init=pp.init$V1+6+10-1

ee.init=read.table(paste0("/fast/AG_Ohler/henriette/PANCREAS_final/output_hg19/open_regions_subset/open_regions_pairs/timeless_init_enh-enh/2_30/classes-15.txt"))
right.ee.init=ee.init$V1+6+10+4-1






pe=read.table(paste0("/fast/AG_Ohler/henriette/PANCREAS_final/output_hg19/open_regions_subset/open_regions_pairs/timeless_multi_prom-enh/10/classes-10_afterEM.txt"))
right.pe=pe$V1+6+-1

pp=read.table(paste0("/fast/AG_Ohler/henriette/PANCREAS_final/output_hg19/open_regions_subset/open_regions_pairs/timeless_multi_prom-prom/4/classes-4_afterEM.txt"))
right.pp=pp$V1+6+10-1

ee=read.table(paste0("/fast/AG_Ohler/henriette/PANCREAS_final/output_hg19/open_regions_subset/open_regions_pairs/timeless_multi_enh-enh/15/classes-15_afterEM.txt"))
right.ee=ee$V1+6+10+4-1



right1=c(right.pe.init, right.pp.init, right.ee.init, right.pe, right.pp, right.ee)



n=length(right1)
n





l2=read.table(paste0("/fast/AG_Ohler/henriette/PANCREAS_final/output_hg19/open_regions_subset/open_regions_pairs/timeless_combined_initmulti/10_4_15/classes-29_afterEM.txt"))
right2=l2$V1+6+29-1



links1=as.data.frame(matrix(c(left,right1,rep(1,n)),byrow=F,ncol=3))



links2=as.data.frame(matrix(c(right1,right2,rep(1,n)),byrow=F,ncol=3))
links=rbind(links1,links2)



names(links) = c("source", "target", "value")

links$group=c(rep("blue",3617),rep("red",800),rep("green",5267),rep("blue",3406),rep("red",687),rep("green",5462))

#links$group[1] <- "blue"
#links$group[2] <- "green"
#mycolor <- 'd3.scaleOrdinal() .domain(["blue", "red", "green"]) .range(["blue", "red", green"])'

#colourScale = JS('d3.scaleOrdinal().range(["#7d3945","#e0677b", "#244457"])')

#colourScale = JS('d3.scaleOrdinal().domain(["blue", "red", "green"]).
#                 range(["blue", "red", "green"])')


sn <- sankeyNetwork(Links = links, Nodes = nodes,
                    Source = "source", Target = "target",
                    Value = "value", NodeID = "name",
                    fontSize= 12, nodeWidth = 30,NodeGroup = "group",LinkGroup = "group", iterations = 0, 
                    colourScale = JS('d3.scaleOrdinal().domain(["blue", "red", "green"]).
                 range(["blue", "red", "green"])'),
                    linkType = 'path1', linkOpacity = 0.3, zoom=T) 
# nodeLabelMargin = 4)

sn

sn %>% saveNetwork(file = paste0("/fast/AG_Ohler/henriette/PANCREAS_final/output_hg19/open_regions_subset/open_regions_pairs/timeless_combined_initmulti/10_4_15/sankey.html"))



sn <- sankeyNetwork(Links = links, Nodes = nodes,
                    Source = "source", Target = "target",
                    Value = "value", NodeID = "name",
                    fontSize= 12, nodeWidth = 30,NodeGroup = "group",LinkGroup = "group", 
                    colourScale = JS('d3.scaleOrdinal().domain(["blue", "red", "green"]).
                 range(["blue", "red", "green"])'),
                    linkType = 'path1', linkOpacity = 0.3, zoom=T) 
#nodeLabelMargin = 4)
sn

sn %>% saveNetwork(file = paste0("/fast/AG_Ohler/henriette/PANCREAS_final/output_hg19/open_regions_subset/open_regions_pairs/timeless_combined_initmulti/10_4_15/sankey_rearranged.html"))





q()


#######################









nodes = data.frame(name = c("init_PE","init_PP","init_EE", "multi_PE","multi_PP","multi_EE", paste0("PE_",1:10), 
                            paste0("PP_",1:5), paste0("EE_",1:17), 
                            paste0("combined_",1:32)), 
                   group=as.character(c("blue","red","green",
                                        "blue","red","green",
                                        rep("blue",10),rep("red",5),rep("green",17),
                                        rep("blue",10),rep("red",5),rep("green",17))))



#init pp pe ee and multi pp pe ee
left=c(rep(1,3617),rep(2,800),rep(3,5267),rep(4,3406),rep(5,687),rep(6,5462))-1



pe.init=read.table(paste0("/fast/AG_Ohler/henriette/PANCREAS_final/output_hg19/open_regions_subset/open_regions_pairs/timeless_init_prom-enh/2_30/classes-10.txt"))
right.pe.init=pe.init$V1+6-1

pp.init=read.table(paste0("/fast/AG_Ohler/henriette/PANCREAS_final/output_hg19/open_regions_subset/open_regions_pairs/timeless_init_prom-prom/2_30/classes-5.txt"))
right.pp.init=pp.init$V1+6+10-1

ee.init=read.table(paste0("/fast/AG_Ohler/henriette/PANCREAS_final/output_hg19/open_regions_subset/open_regions_pairs/timeless_init_enh-enh/2_30/classes-17.txt"))
right.ee.init=ee.init$V1+6+10+5-1






pe=read.table(paste0("/fast/AG_Ohler/henriette/PANCREAS_final/output_hg19/open_regions_subset/open_regions_pairs/timeless_multi_prom-enh/10/classes-10_afterEM.txt"))
right.pe=pe$V1+6+-1

pp=read.table(paste0("/fast/AG_Ohler/henriette/PANCREAS_final/output_hg19/open_regions_subset/open_regions_pairs/timeless_multi_prom-prom/5/classes-5_afterEM.txt"))
right.pp=pp$V1+6+10-1

ee=read.table(paste0("/fast/AG_Ohler/henriette/PANCREAS_final/output_hg19/open_regions_subset/open_regions_pairs/timeless_multi_enh-enh/17/classes-17_afterEM.txt"))
right.ee=ee$V1+6+10+5-1



right1=c(right.pe.init, right.pp.init, right.ee.init, right.pe, right.pp, right.ee)



n=length(right1)
n





l2=read.table(paste0("/fast/AG_Ohler/henriette/PANCREAS_final/output_hg19/open_regions_subset/open_regions_pairs/timeless_combined_initmulti/10_5_17/classes-32_afterEM.txt"))
right2=l2$V1+6+32-1



links1=as.data.frame(matrix(c(left,right1,rep(1,n)),byrow=F,ncol=3))



links2=as.data.frame(matrix(c(right1,right2,rep(1,n)),byrow=F,ncol=3))
links=rbind(links1,links2)



names(links) = c("source", "target", "value")

links$group=c(rep("blue",3617),rep("red",800),rep("green",5267),rep("blue",3406),rep("red",687),rep("green",5462))

#links$group[1] <- "blue"
#links$group[2] <- "green"
#mycolor <- 'd3.scaleOrdinal() .domain(["blue", "red", "green"]) .range(["blue", "red", green"])'

#colourScale = JS('d3.scaleOrdinal().range(["#7d3945","#e0677b", "#244457"])')

#colourScale = JS('d3.scaleOrdinal().domain(["blue", "red", "green"]).
#                 range(["blue", "red", "green"])')


sn <- sankeyNetwork(Links = links, Nodes = nodes,
                    Source = "source", Target = "target",
                    Value = "value", NodeID = "name",
                    fontSize= 12, nodeWidth = 30,NodeGroup = "group",LinkGroup = "group", iterations = 0, 
                    colourScale = JS('d3.scaleOrdinal().domain(["blue", "red", "green"]).
                 range(["blue", "red", "green"])'),
                    linkType = 'path1', linkOpacity = 0.3, zoom=T) 
# nodeLabelMargin = 4)

sn

sn %>% saveNetwork(file = paste0("/fast/AG_Ohler/henriette/PANCREAS_final/output_hg19/open_regions_subset/open_regions_pairs/timeless_combined_initmulti/10_5_17/sankey.html"))



sn <- sankeyNetwork(Links = links, Nodes = nodes,
                    Source = "source", Target = "target",
                    Value = "value", NodeID = "name",
                    fontSize= 12, nodeWidth = 30,NodeGroup = "group",LinkGroup = "group", 
                    colourScale = JS('d3.scaleOrdinal().domain(["blue", "red", "green"]).
                 range(["blue", "red", "green"])'),
                    linkType = 'path1', linkOpacity = 0.3, zoom=T) 
#nodeLabelMargin = 4)
sn

sn %>% saveNetwork(file = paste0("/fast/AG_Ohler/henriette/PANCREAS_final/output_hg19/open_regions_subset/open_regions_pairs/timeless_combined_initmulti/10_5_17/sankey_rearranged.html"))
