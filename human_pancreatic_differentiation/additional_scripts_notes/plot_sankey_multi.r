#https://github.com/fbreitwieser/sankeyD3
#devtools::install_github("fbreitwieser/sankeyD3")
library(sankeyD3)
library(dplyr)



nodes = data.frame(name = c("multi_PE","multi_PP","multi_EE", paste0("PE_",1:10), 
                            paste0("PP_",1:4), paste0("EE_",1:15), 
                            paste0("combined_",1:29)), 
                   group=as.character(c("blue","red","green",
                                        rep("blue",10),rep("red",4),rep("green",15),
                                        rep("blue",10),rep("red",4),rep("green",15))))



#multi pp pe ee
left=c(rep(1,3406),rep(2,687),rep(3,5462))-1



pe=read.table(paste0("/fast/AG_Ohler/henriette/PANCREAS_final/output_hg19/open_regions_subset/open_regions_pairs/timeless_multi_prom-enh/10/classes-10_afterEM.txt"))
right.pe=pe$V1+3-1

pp=read.table(paste0("/fast/AG_Ohler/henriette/PANCREAS_final/output_hg19/open_regions_subset/open_regions_pairs/timeless_multi_prom-prom/4/classes-4_afterEM.txt"))
right.pp=pp$V1+3+10-1

ee=read.table(paste0("/fast/AG_Ohler/henriette/PANCREAS_final/output_hg19/open_regions_subset/open_regions_pairs/timeless_multi_enh-enh/15/classes-15_afterEM.txt"))
right.ee=ee$V1+3+10+4-1



right1=c(right.pe, right.pp, right.ee)



n=length(right1)
n





l2=read.table(paste0("/fast/AG_Ohler/henriette/PANCREAS_final/output_hg19/open_regions_subset/open_regions_pairs/timeless_combined_multi/10_4_15/classes-29_afterEM.txt"))
right2=l2$V1+3+29-1



links1=as.data.frame(matrix(c(left,right1,rep(1,n)),byrow=F,ncol=3))






  
links2=as.data.frame(matrix(c(right1,right2,rep(1,n)),byrow=F,ncol=3))
links=rbind(links1,links2)



names(links) = c("source", "target", "value")

links$group=c(rep("blue",3406),rep("red",687),rep("green",5462))

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

sn %>% saveNetwork(file = paste0("/fast/AG_Ohler/henriette/PANCREAS_final/output_hg19/open_regions_subset/open_regions_pairs/timeless_combined_multi/10_4_15/sankey_multi.html"))



sn <- sankeyNetwork(Links = links, Nodes = nodes,
                    Source = "source", Target = "target",
                    Value = "value", NodeID = "name",
                    fontSize= 12, nodeWidth = 30,NodeGroup = "group",LinkGroup = "group", 
                    colourScale = JS('d3.scaleOrdinal().domain(["blue", "red", "green"]).
                 range(["blue", "red", "green"])'),
                    linkType = 'path1', linkOpacity = 0.3, zoom=T) 
                    #nodeLabelMargin = 4)
sn

sn %>% saveNetwork(file = paste0("/fast/AG_Ohler/henriette/PANCREAS_final/output_hg19/open_regions_subset/open_regions_pairs/timeless_combined_multi/10_4_15/sankey_multi_rearranged.html"))

q()

############################


#https://github.com/fbreitwieser/sankeyD3
#devtools::install_github("fbreitwieser/sankeyD3")
library(sankeyD3)
library(dplyr)



nodes = data.frame(name = c("multi_PE","multi_PP","multi_EE", paste0("PE_",1:10), 
                            paste0("PP_",1:5), paste0("EE_",1:17), 
                            paste0("combined_",1:32)), 
                   group=as.character(c("blue","red","green",
                                        rep("blue",10),rep("red",5),rep("green",17),
                                        rep("blue",10),rep("red",5),rep("green",17))))



#multi pp pe ee
left=c(rep(1,3406),rep(2,687),rep(3,5462))-1



pe=read.table(paste0("/fast/AG_Ohler/henriette/PANCREAS_final/output_hg19/open_regions_subset/open_regions_pairs/timeless_multi_prom-enh/10/classes-10_afterEM.txt"))
right.pe=pe$V1+3-1

pp=read.table(paste0("/fast/AG_Ohler/henriette/PANCREAS_final/output_hg19/open_regions_subset/open_regions_pairs/timeless_multi_prom-prom/5/classes-5_afterEM.txt"))
right.pp=pp$V1+3+10-1

ee=read.table(paste0("/fast/AG_Ohler/henriette/PANCREAS_final/output_hg19/open_regions_subset/open_regions_pairs/timeless_multi_enh-enh/17/classes-17_afterEM.txt"))
right.ee=ee$V1+3+10+5-1



right1=c(right.pe, right.pp, right.ee)



n=length(right1)
n





l2=read.table(paste0("/fast/AG_Ohler/henriette/PANCREAS_final/output_hg19/open_regions_subset/open_regions_pairs/timeless_combined_multi/10_5_17/classes-32_afterEM.txt"))
right2=l2$V1+3+32-1



links1=as.data.frame(matrix(c(left,right1,rep(1,n)),byrow=F,ncol=3))







links2=as.data.frame(matrix(c(right1,right2,rep(1,n)),byrow=F,ncol=3))
links=rbind(links1,links2)



names(links) = c("source", "target", "value")

links$group=c(rep("blue",3406),rep("red",687),rep("green",5462))

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

sn %>% saveNetwork(file = paste0("/fast/AG_Ohler/henriette/PANCREAS_final/output_hg19/open_regions_subset/open_regions_pairs/timeless_combined_multi/10_5_17/sankey_multi.html"))



sn <- sankeyNetwork(Links = links, Nodes = nodes,
                    Source = "source", Target = "target",
                    Value = "value", NodeID = "name",
                    fontSize= 12, nodeWidth = 30,NodeGroup = "group",LinkGroup = "group", 
                    colourScale = JS('d3.scaleOrdinal().domain(["blue", "red", "green"]).
                 range(["blue", "red", "green"])'),
                    linkType = 'path1', linkOpacity = 0.3, zoom=T) 
#nodeLabelMargin = 4)
sn

sn %>% saveNetwork(file = paste0("/fast/AG_Ohler/henriette/PANCREAS_final/output_hg19/open_regions_subset/open_regions_pairs/timeless_combined_multi/10_5_17/sankey_multi_rearranged.html"))





