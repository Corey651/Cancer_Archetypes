library(ComplexHeatmap)
library(circlize)
library(matrixStats)
library(RColorBrewer)
library(GetoptLong)
library(GenomicRanges)
library(ggplot2)

tmat=read.table("CNA_Table_Bin.txt",header=FALSE,sep=',')
tmat=as.matrix(tmat)

data_col=colorRamp2(c(-0.25,-0.07,0.07,0.25),c("#377EB8","#FFFFFF","#FFFFFF","#E41A1C"))
bin_chr=read.table("bins.txt")
bin_chr=bin_chr[1:20000,]

chr_col=structure(names=c(1:24),rep(c("#AEAEAE","#E6E6E6"),times=12))
#################################

# V2

tmat3=tmat
F=dim(tmat3)
rownames(tmat3)=c("Survival","Proliferation","Fibroblastic","Energy","Biomass","Senescence")
Q=apply(abs(tmat3),2,max)
pdf(file = "Plot_CNA.pdf",width=9,height=4)

column_ha=HeatmapAnnotation(Max.Corr.=anno_lines(Q,add_points=TRUE,pt_gp=gpar(col=ifelse(1:20000 %in% c(829,832,866,871,874,910,6354,7008,7009,8453,8474,8477,8791,8825,12238,12248,12975,13055,13101,17027,17344,17368,18674,18677,18691,18725,18727,18732),"red","NA"))))
ha=HeatmapAnnotation(
  Chr=bin_chr, col=list(Chr=chr_col),annotation_legend_param=list(
    Chr=list(title="Chr")
  ), show_legend=FALSE,annotation_name_side="right",show_annotation_name=FALSE
)



ht_list = Heatmap(tmat3, col = data_col, name = "Spearman Correlation", row_names_side="left",
                  show_row_dend = FALSE, show_column_dend = FALSE,row_order=1:F[1],column_order=1:F[2],
                  show_column_names = FALSE,top_annotation=column_ha, 
                  row_title_gp = gpar(col = "#FFFFFF00"),border=TRUE,
                  use_raster=TRUE,raster_quality=5,heatmap_legend_param=list(at=c(-0.25,0,0.25),labels=c("-0.25","NS","0.25"),
                                                                             nrow=1,border=TRUE,title_position="leftcenter-rot"),bottom_annotation=ha) %v%
  HeatmapAnnotation(label=anno_mark(at=c(12975,8453,13055,13101,18725,18674,7008,8825,866,866,12238,12248,8791,18674,6354,18677,871,18732,7009,8474,8477,910,17027,874,17344,17368,18691,832,829,18727),labels=c("ATM","CARD11","CBL","CHEK1","CHEK2","CRKL","CSF1R","EGFR","H3C13","H3C14","HRAS","IGF2","IKZF1","LZTR1","MAP3K1","MAPK1","MCL1","NF2","PDGFRB","PMS2","RAC1","RIT1","RPTOR","SETDB1","SMAD2","SMAD4","SMARCB1","TENT5C","VTCN1","ZNRF3"),side="bottom"))

draw(ht_list,
     annotation_legend_side = "bottom",heatmap_legend_side = "right")
decorate_annotation("Chr",{
  grid.text("Chr.",.537,-.8)
})

dev.off()
