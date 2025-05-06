setwd("path/to/working/directory")
rm (list = ls())



# Load Libraries                    ####
library(openxlsx)
library(dplyr)
library(tidyr)
library(ComplexUpset)
library(ggplot2)
library(data.table)
library(org.Hs.eg.db)
library(clusterProfiler)
library(ggpubr)

### Supplementary Figure 1 ###
colors <- RColorBrewer::brewer.pal(n = 10, name = "Paired")
colors_final <- c(colors[8],
                  colorRampPalette(c(colors[3], colors[4], "black"))(6)[4],
                  colorRampPalette(c(colors[3], colors[4], "black"))(6)[2],
                  colorRampPalette(c(colors[1], colors[2], "black"))(6)[4],
                  colorRampPalette(c(colors[1], colors[2], "black"))(6)[2],
                  colorRampPalette(c(colors[9], colors[10], "black"))(6)[4],
                  colorRampPalette(c(colors[9], colors[10], "black"))(6)[2],
                  colorRampPalette(c(colors[5], colors[6], "black"))(6)[4],
                  colorRampPalette(c(colors[5], colors[6], "black"))(6)[2])

names <- c("Human1",
           "Human1-O1", "Human1-O2",
           "Human1-D1", "Human1-D2",
           "Human1-T1", "Human1-T2",
           "Human1-S1", "Human1-S2", 
           "Human1-O1-SDL", "Human1-O2-SDL", 
           "Human1-D1-SDL", "Human1-D2-SDL",
           "Human1-T1-SDL", "Human1-T2-SDL", 
           "Human1-S1-SDL", "Human1-S2-SDL")

gMCSs.list <- list()

for (model in names){
  gMCSs.list[[model]] <- as.data.frame(read.xlsx('calculated_gMCSs_gMIS.xlsx', sheet = model))
  gMCSs.list[[model]] <- t(apply(gMCSs.list[[model]], 1, function(y){
    sub("_KO", "", y)
  }))
  gMCSs.list[[model]] <- unique(apply(gMCSs.list[[model]],1,function(y){paste(sort(y[y!=""]),collapse = "--")}))
}

names <- c("Human1-gMCS",
           "Human1-O1-gMCS", 
           "Human1-O2-gMCS",
           "Human1-D1-gMCS","Human1-D2-gMCS",
           "Human1-T1-gMCS", "Human1-T2-gMCS",
           "Human1-S1-gMCS", "Human1-S2-gMCS", 
           "Human1-O1-gMIS", "Human1-O2-gMIS", 
           "Human1-D1-gMIS", "Human1-D2-gMIS",
           "Human1-T1-gMIS","Human1-T2-gMIS", 
           "Human1-S1-gMIS", "Human1-S2-gMIS")

nameVector <- unlist(mapply(function(x,y){ rep(y, length(x)) }, gMCSs.list, names))

resultDF <- cbind.data.frame(unlist(gMCSs.list), nameVector)
colnames(resultDF) <- c("ENSEMBL", "layer")

model_array <- c("Human1-O", "Human1-D", "Human1-T", "Human1-S")

custom_geom_point <- function(data) {
  data$shape <- ifelse(grepl("SDL", data$category), 21, 22)  # Cambia la forma dependiendo de la condición
  ggplot2::geom_point(data=data, aes(shape=shape))
}

a <- 0
pos <- 1
plot_final_gmcs <- list()
for (model in model_array){
  resultDF_filtered <- resultDF[resultDF$layer %in% c(paste0(model, c("1-gMCS", "2-gMCS"), sep = ""), paste0(model, c("1-gMIS", "2-gMIS"), sep = "")),]
  
  data <- resultDF_filtered %>% 
    pivot_wider(id_cols = ENSEMBL,
                names_from = layer, 
                values_from = layer, 
                values_fn = list(x = length), 
                values_fill = list(x = 0))
  data <- t(data)
  colnames(data) <- data[1,]
  
  data <- data[-1,]
  data[!is.na(data)] <- TRUE
  data[is.na(data)] <- FALSE
  
  data <- as.data.frame(t(data))
  
  layer <- c(paste0(model, c("1-gMCS","1-gMIS", "2-gMCS", "2-gMIS"), sep = ""))
  names_colors <- setNames(rep(colors_final[(2+a):(3+a)], each = 2), layer)
  
  plot_final_gmcs[[pos]] <- ComplexUpset::upset(
    data,
    rev(layer), sort_sets = F,  stripes='white', name='Model',
    base_annotations=list(
      'gMCSs Intersection'=intersection_size(
        counts=TRUE,
        mapping=aes(fill='bars_color')
      ) + scale_fill_manual(values=c('bars_color'='black'), guide = 'none')
    ),matrix=(
      intersection_matrix(geom=geom_point(shape='circle filled', size=3))
      + scale_color_manual(
        values=names_colors,
        guide='none'
      )
    ),
    queries=list(
      upset_query(set=layer[1], fill= colors_final[2+a]),
      upset_query(set=layer[2], shape = 'square', fill= colors_final[2+a]),
      upset_query(set=layer[3], fill= colors_final[3+a]),
      upset_query(set=layer[4], shape = 'square', fill=colors_final[3+a])
    ),
    width_ratio=0.2, wrap=TRUE,  
    set_sizes=FALSE, 
    themes=upset_modify_themes(
      list(
        'gMISs Intersection'=theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_blank(),
                                   axis.title.x = element_blank(),
                                   panel.background = element_blank(), axis.line = element_line(colour = "black"))))
  )
  
  pos <- pos+1
  a<- a+2
}


final_plot1 <- (
  plot_final_gmcs[[1]] + ggtitle("Omnipath") + plot_final_gmcs[[2]] + ggtitle("Dorothea")  +
    plot_final_gmcs[[3]] + ggtitle("TRRUST")  + plot_final_gmcs[[4]] + ggtitle("Signor")  
)

### SL models ###
resultDF_filtered <- resultDF

layer <- c("Human1-gMCS",
           "Human1-O1-gMCS", "Human1-O1-gMIS",
           "Human1-O2-gMCS", "Human1-O2-gMIS",
           "Human1-D1-gMCS", "Human1-D1-gMIS",
           "Human1-D2-gMCS", "Human1-D2-gMIS",
           "Human1-T1-gMCS", "Human1-T1-gMIS",
           "Human1-T2-gMCS", "Human1-T2-gMIS", 
           "Human1-S1-gMCS", "Human1-S1-gMIS",
           "Human1-S2-gMCS", "Human1-S2-gMIS")

data <- resultDF_filtered %>% 
  pivot_wider(id_cols = ENSEMBL,
              names_from = layer, 
              values_from = layer, 
              values_fn = list(x = length), 
              values_fill = list(x = 0))
data <- t(data)
colnames(data) <- data[1,]

data <- data[-1,]
data[!is.na(data)] <- TRUE
data[is.na(data)] <- FALSE

data <- as.data.frame(t(data))

colors_final2 <- c(colors_final[1], rep(colors_final[2:length(colors_final)], each = 2))

names_colors <- setNames(colors_final2, layer)

plot_final_SL <- ComplexUpset::upset(
  data,
  rev(layer), sort_sets = F,  stripes='white', name='Model',
  base_annotations=list(
    'gMCSs Intersection'=intersection_size(
      counts=TRUE,
      mapping=aes(fill='bars_color')
      
    ) + scale_fill_manual(values=c('bars_color'='black'), guide = 'none')
  ),matrix=(
    intersection_matrix(geom=geom_point(shape='circle filled', size=3))
    + scale_color_manual(
      values=names_colors,
      guide='none'
    )
  ),
  queries=list(
    
    upset_query(set=layer[1], fill= colors_final2[1]),
    upset_query(set=layer[2], fill= colors_final2[2]),
    upset_query(set=layer[3], fill= colors_final2[3]),
    upset_query(set=layer[4], fill=colors_final2[4]),
    upset_query(set=layer[5], fill=colors_final2[5]),
    upset_query(set=layer[6], fill=colors_final2[6]),
    
    upset_query(set=layer[7], fill= colors_final2[7]),
    upset_query(set=layer[8], fill= colors_final2[8]),
    upset_query(set=layer[9], fill= colors_final2[9]),
    upset_query(set=layer[10], fill=colors_final2[10]),
    upset_query(set=layer[11], fill=colors_final2[11]),
    upset_query(set=layer[12], fill=colors_final2[12]),
    upset_query(set=layer[13], fill=colors_final2[13]),
    upset_query(set=layer[14], fill=colors_final2[14]),
    upset_query(set=layer[15], fill=colors_final2[15]),
    upset_query(set=layer[16], fill=colors_final2[16]),
    upset_query(set=layer[17], fill=colors_final2[17])
    
  ),
  width_ratio=0.2, wrap=TRUE,  
  set_sizes=FALSE, 
  themes=upset_modify_themes(
    list(
      'gMISs Intersection'=theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_blank(),
                                 axis.title.x = element_blank(),
                                 panel.background = element_blank(), axis.line = element_line(colour = "black")))))

final_plot3 <- (plot_final_SL + ggtitle("All integrated models")) / final_plot1 
ggsave(paste0("./Plots/Supp_Fig_1.png"), plot = final_plot3,
       device = "png", width = 14, height = 14, units = "in", dpi = 300, bg = "white")


### Supplementary Figure 2 ###

model <- c("Human1-O1","Human1-O2",
           "Human1-D1","Human1-D2",
           "Human1-T1","Human1-T2",
           "Human1-S1", "Human1-S2")

model_name <- c("Omnipath", "Omnipath", 
                "Dorothea", "Dorothea",
                "TRRUST", "TRRUST", 
                "Signor", "Signor")

layers <- rep(1:2, 4)

data_final_final <- data.frame()

names <- model[1]
final <- list()
cont <- 1
for (names in model){
  
  layer <- paste0(names, "-SDL")
  
  gMCSs <- read.xlsx("calculated_gMCSs_gMIS.xlsx", sheet = layer)
  
  gMCSs.ENSEMBL <- as.matrix(gMCSs)
  rownames(gMCSs.ENSEMBL) <- as.character(seq(1,nrow(gMCSs.ENSEMBL)))
  gMCSs.ENSEMBL <- unique(gMCSs.ENSEMBL)
  gMCSs.ENSEMBL[is.na(gMCSs.ENSEMBL)] <- ""
  
  lengths <- apply(gMCSs.ENSEMBL, 2, function(x){
    pos <- which(x == "")
    pos[length(pos)]
  })
  
  lengths <- do.call(cbind,lengths)
  lengths[length(lengths)+1] <- nrow(gMCSs.ENSEMBL)

  pos_off <- t(apply(gMCSs, 1, function(y){
    x <- grepl("_off", y)
    as.numeric(x)
  }))
  
  pos_on <- t(apply(gMCSs, 1, function(y){
    x <- !grepl("_off", y) & grepl("ENS", y) 
    as.numeric(x)
  }))
  
  
  pos_off <- rowSums(pos_off)
  pos_on <- rowSums(pos_on)
  
  times <- array()
  times[1] <- lengths[1]
  
  for (i in 2:length(lengths)){
    times[i] <- lengths[i]-lengths[i-1]
  }
  
  ideal <- rep(c(1:ncol(gMCSs)), times = times)
  
  pos_ON_final <- which(pos_on == ideal)
  pos_Off_final <- which(pos_off == ideal)
  
  all <- c(1:nrow(gMCSs))
  pos_rest <- all[c(pos_ON_final, pos_Off_final)]
  pos_comb_FINAL <- setdiff(all, pos_rest)
  
  counts_on <- array()
  counts_off <- array()
  counts_comb <- array()
  a <- 0
  
  for (i in 1:length(lengths)){
    counts_on[i] <- length(which(pos_ON_final <= lengths[i] & pos_ON_final > a))
    counts_off[i] <- length(which(pos_Off_final <= lengths[i] & pos_Off_final > a))
    counts_comb[i] <- length(which(pos_comb_FINAL <= lengths[i] & pos_comb_FINAL > a))
    a <- lengths[i]
  }
  
  
  data_final <- rbind(counts_on, counts_off, counts_comb)
  data <- reshape2::melt(as.matrix(data_final ) )
  data <- data[-3,]
  
  data$model <- model_name[cont]
  data$layers <- layers[cont]
  
  data_final_final <- rbind(data_final_final, data)
  cont <- cont + 1
  
}

data <- data_final_final
levels(data$Var1) <- c("-", "+", "-/+")
data$layers <- factor(data$layers)

names <- c("Omnipath", "Dorothea", "TRRUST", "Signor")


colors <- RColorBrewer::brewer.pal(n = 10, name = "Paired")
colors_final <- c(colorRampPalette(c(colors[3], colors[4], "black"))(6)[4],
                  colorRampPalette(c(colors[3], colors[4], "black"))(6)[2],
                  
                  colorRampPalette(c(colors[1], colors[2], "black"))(6)[4],
                  colorRampPalette(c(colors[1], colors[2], "black"))(6)[2],
                  
                  colorRampPalette(c(colors[9], colors[10], "black"))(6)[4],
                  colorRampPalette(c(colors[9], colors[10], "black"))(6)[2],
                  
                  colorRampPalette(c(colors[5], colors[6], "black"))(6)[4],
                  colorRampPalette(c(colors[5], colors[6], "black"))(6)[2])
cont <-   1
model <- names[1]

for (model in names){
  c3 <- list()
  a <- data[data$Var2 == 1 & data$model %in% model,]
  c3[[1]] <- ggbarplot(a, "Var1", "value",
                       fill = "layers", color = "layers", 
                       label = T, title = paste0("length = 1"), 
                       position = position_dodge(0.9),
                       xlab = "gMCS type", ylab = "number of gMCSs", lab.vjust = 0, lab.size = 3, 
                       palette = colors_final[cont:(cont+1)]) + 
    geom_vline(aes(xintercept = c(1.5,0,1.5, 0)), linetype = "dashed",  alpha = 0.5) + 
    theme(plot.title = element_text(hjust = 0.5, size = 12), axis.title.x = element_blank()) + labs(fill = "number of layers", color = "number of layers")
  
  for (i in 2:length(lengths)){
    a <- data[data$Var2 == i & data$model %in% model,]
    c3[[i]] <- ggbarplot(a, "Var1", "value",
                         fill = "layers", color = "layers", 
                         label = T, title = paste0("length = ", i), #lab.pos = "in",
                         position = position_dodge(0.9),
                         xlab = "gMCS type", ylab = "number of gMCSs", lab.vjust = 0, 
                         lab.size = 3, 
                         palette = colors_final[cont:(cont+1)]) +  geom_vline(aes(xintercept = rep(c(1.5, 2.5, 0), 2)), linetype = "dashed",  alpha = 0.5) + 
      theme(plot.title = element_text(hjust = 0.5, size = 12), axis.title.x = element_blank(), axis.title.y = element_blank()) + labs(fill = "number of layers", color = "number of layers")
    
  }
  cont <- cont+2
  final[[model]]<- annotate_figure(ggarrange(plotlist = c3, ncol = length(lengths), common.legend = T, 
                                             legend = "bottom", widths = c(0.8,1,1,1,1)), top = text_grob(paste0(model), face = "bold", size = 15))
  
}

plot <- ggarrange( plotlist = final, nrow = 4)
ggsave(paste0("Plots/Supp_Fig_2.png"), plot, bg = "white", height = 12, width = 14)


### Supplementary Figure 3 ###
CCLE_Exp <- as.data.frame(fread("DepMap_23Q4/OmicsExpressionProteinCodingGenesTPMLogp1.csv"))
colnames(CCLE_Exp) <- gsub("[\\(\\)]", "", regmatches(colnames(CCLE_Exp), gregexpr("\\(.*?\\)", colnames(CCLE_Exp))))
genes <- colnames(CCLE_Exp)
CCLE_Exp <- as.data.frame(t(CCLE_Exp))

colnames(CCLE_Exp) <- CCLE_Exp[1,]
CCLE_Exp$ENTREZID <- genes
toENSEMBL <- clusterProfiler::bitr(CCLE_Exp[,"ENTREZID"], "ENTREZID", "ENSEMBL", org.Hs.eg.db)
CCLE_Exp <- merge(toENSEMBL, CCLE_Exp, by = "ENTREZID")
CCLE_Exp <- CCLE_Exp[!duplicated(CCLE_Exp$ENTREZID),]
CCLE_Exp <- CCLE_Exp[!duplicated(CCLE_Exp$ENSEMBL),]
rownames(CCLE_Exp) <- CCLE_Exp[,2]
CCLE_Exp <- CCLE_Exp[,-c(1,2)]
CCLE_genes <- rownames(CCLE_Exp)
CCLE_Exp<- apply(CCLE_Exp,2,as.numeric)
rownames(CCLE_Exp) <- CCLE_genes

Achilles <- as.data.frame(fread("DepMap_23Q4/CRISPRGeneEffect.csv"))
colnames(Achilles) <- gsub("[\\(\\)]", "", regmatches(colnames(Achilles), gregexpr("\\(.*?\\)", colnames(Achilles))))
genes <- colnames(Achilles)
Achilles <- as.data.frame(t(Achilles))

colnames(Achilles) <- Achilles[1,]
Achilles$ENTREZID <- genes
toENSEMBL <- clusterProfiler::bitr(Achilles[,"ENTREZID"], "ENTREZID", "ENSEMBL", org.Hs.eg.db)
Achilles <- merge(toENSEMBL, Achilles, by = "ENTREZID")
Achilles <- Achilles[!duplicated(Achilles$ENTREZID),]
Achilles <- Achilles[!duplicated(Achilles$ENSEMBL),]
rownames(Achilles) <- Achilles[,2]
Achilles <- Achilles[,-c(1,2)]
Achilles_genes <- rownames(Achilles)
Achilles<- apply(Achilles,2,as.numeric)
rownames(Achilles) <- Achilles_genes

Achilles_binary <- t(apply(Achilles, 1, function(x){
  ifelse(x < -0.6, 1, 0) 
}))

epsilon <- epsilon <- 1e-7
threshold <- fread("th_log2TPM_CCLE.txt")
gene_essential <- "ENSG00000106462"
gMCS <- c("ENSG00000106462", "ENSG00000177606", "ENSG00000173039_off")  #EZH2, JUN, RELA
gMCS_on <- c("ENSG00000106462", "ENSG00000177606")
gMCS_off <- c("ENSG00000173039")

exp_on <- CCLE_Exp[gMCS_on, threshold$V1, drop = FALSE]
rownames(exp_on) <- gMCS_on
exp_on_binary <- as.data.table(t(apply(exp_on, 1, function(x) (x - threshold$threshold) > epsilon)))


exp_off <- CCLE_Exp[gMCS_off, , drop = FALSE]
rownames(exp_off) <- gMCS_off
exp_off_binary <-  as.data.table(t(apply(exp_off, 1, function(x) (x - threshold$threshold) < epsilon)))
final_binary <- rbind(exp_on_binary, exp_off_binary)

final_binary <- as.data.frame(final_binary)
rownames(final_binary) <- gMCS
one_hit <- which(colSums(final_binary) == 1)
essential_cells <- c()
for (i in 1:length(one_hit)) {
  essential_cells <- c(essential_cells, names(one_hit)[i])
  # col <- one_hit[i]
  # act_gene <- rownames(final_binary)[which(final_binary[,col])]
  # 
  # model_result$binary_genes[[act_gene]] <- c(model_result$binary_genes[[act_gene]], cell_line)
  # model_result$gmcs_pos[[act_gene]] <- c(model_result$gmcs_pos[[act_gene]], pos)
  # model_result$mat_update[act_gene, cell_line] <- 1
}

# essential_cells <- binary_genes_MM[[model]][[gene_essential]][gmcs.pos_MM[[model]][[gene_essential]] %in% pos_pos]

# essential_cells <- binary_genes_MM[[model]][[gene_essential]][gmcs.pos_MM[[model]][[gene_essential]] %in% pos_pos]

annot <- data.frame(essential = rep(0,ncol(CCLE_Exp)), row.names = colnames(CCLE_Exp))

annot$essential[colnames(CCLE_Exp) %in% essential_cells] <- 1

# gene_essential_no_off <- gsub(pattern = "_off", replacement = "", gene_essential)

gene_no_off <- gsub(pattern = "_off", replacement = "", gMCS)

genes_symbol <- bitr(gsub(pattern = "_off", replacement = "", gMCS), "ENSEMBL", "SYMBOL", drop = F, org.Hs.eg.db)

off.pos <- grep(pattern = "_off", gMCS)
genes_2 <- gsub(pattern = "_off", replacement = "", gMCS)
mat <- match(genes_2, genes_symbol$ENSEMBL)

final <- genes_symbol$SYMBOL[mat]
final[off.pos] <-  paste(final[off.pos], "_off", sep = "")

gene.ratio <- CCLE_Exp[gene_no_off,]
for (i in 1:ncol(CCLE_Exp)){
  gene.ratio[gene_no_off,i] <- CCLE_Exp[gene_no_off,i]/threshold$threshold[i]
}
gc()
final_exp <- as.matrix(log2(gene.ratio[gene_no_off,]))
final_exp[final_exp>10] <- 10
final_exp[final_exp<(-10)] <- (-10)
rownames(final_exp) <- final


Achilles_binary2 <- Achilles_binary[,colnames(Achilles_binary) %in% colnames(CCLE_Exp)]
samples <- colnames(CCLE_Exp)[!colnames(CCLE_Exp) %in% colnames(Achilles_binary2)]

# null_dataframe <- data.frame(patients = samples)
null_dataframe <- data.frame(matrix(NA, nrow = nrow(Achilles_binary2), ncol = length(samples)))
colnames(null_dataframe) <- samples

Achilles_binary <- cbind(Achilles_binary2, null_dataframe)
Achilles_binary <- Achilles_binary[,colnames(CCLE_Exp)]

clinical_data <- as.data.frame(fread("DepMap_23Q4/Model.csv"))


cell_lines <- colnames(CCLE_Exp)
clinical_data <- clinical_data[clinical_data$ModelID %in% cell_lines,]
rownames(clinical_data) <- clinical_data$ModelID
clinical_data <- clinical_data[cell_lines,]
# sample_class <- factor(clinical_data$OncotreePrimaryDisease)
sample_class <- clinical_data$OncotreeLineage
sample_class[sample_class == ""] <- "Non_cancerous"
sample_class <- factor(sample_class)


colors_heatmap = "blue_white_red"

Achilles_binary <- Achilles_binary[gene_essential,]
col_anot <- data.frame(lineage = sample_class, DepMap = factor(t(Achilles_binary)))
colnames(col_anot)[2] <- "DepMap_essentiality"

annot$essential <- factor(annot$essential)
rownames(col_anot) <- colnames(final_exp)
col_anot <- merge(col_anot, annot, by = 0)
rownames(col_anot) <- col_anot$Row.names
col_anot <- col_anot[,-1]
colnames(col_anot) <- c("lineage", "DepMap essentiality", "gMIS essentiality")

col_anot$`DepMap essentiality` <- factor(col_anot$`DepMap essentiality`)
col_anot$`gMIS essentiality` <- factor(col_anot$`gMIS essentiality`)

gene.set.list <- lapply(levels(sample_class), function(x){final_exp[,sample_class == x]})
names(gene.set.list) <- levels(sample_class)
sample.class.list <- lapply(levels(sample_class), function(ch){
  z <- sample_class[sample_class == ch]
  z <- factor(as.character(z),levels=unique(as.character(z), fromLast = T))
  return(z)
})
names(sample.class.list) <- levels(sample_class)
plot.list <- list()
for (ch in levels(sample_class)) {

  # order by gene essentiality
  if(!is.null(nrow(gene.set.list[[ch]]))){
    gene.set.rest <- gene.set.list[[ch]][-1,,drop = F]
    idx <- order(gene.set.list[[ch]][1,] - apply(gene.set.rest,2,max), decreasing = F)

    sample.class.list[[ch]] <- sample.class.list[[ch]][idx]
    gene.set.list[[ch]] <- gene.set.list[[ch]][,idx]
  }
  
  
}
gene.set.final <- do.call(cbind, gene.set.list)

anot.colors <- scales::hue_pal()(length(levels(sample_class)))
names(anot.colors) <- levels(as.factor(col_anot$lineage))

set.seed(2021)
prueba <- as.list(c('k', pals::polychrome(ncol(col_anot)-1)))

prueba <- as.list(c('k', a = "#009999", b = "#9900CC"))
names(prueba) <- colnames(col_anot)
prueba[[1]] <- anot.colors
prueba[2:ncol(col_anot)] <- lapply(prueba[2:ncol(col_anot)], function(x){zz <- c("#F5F5F5",x); names(zz) <- c("0","1"); return(zz)})


rownames(gene.set.final)[3] <- "RELA"
plot <- pheatmap::pheatmap(gene.set.final,
                 cluster_rows = F, cluster_cols = F, annotation_col = col_anot,
                 gaps_col = cumsum(table(sample_class)),
                 show_rownames = T, show_colnames = F,  
                 color = colorRampPalette(unlist(strsplit(colors_heatmap,'_')))(100), 
                 annotation_colors = prueba,
                 breaks = seq(-1,+1, length.out = 100),
                 scale = "none", 
                 annotation_legend = T,
                 legend = T,
                 silent = T)[[4]]

ggsave(paste0("./Plots/Supp_Fig_3.png"),
       plot = plot, device = "png", width = 11, height = 8.5, units = "in", dpi = 300, bg = "white")

gc()

### Supplementary Figure 4 ###
genes <- c("ENSG00000106462", "ENSG00000177606", "ENSG00000173039")

cells <- intersect(colnames(CCLE_Exp), colnames(Achilles))
CCLE_Exp <- CCLE_Exp[,cells]
ess <- Achilles[,cells]

data <- data.frame(exp = (CCLE_Exp[genes[2],]-CCLE_Exp[genes[3],]), ess = ess[genes[1],])
colnames(data) <- c("exp", "ess")

correlation_plot3 <- ggpubr::ggscatter(data, x = "exp", y = 'ess', 
                                       add = "reg.line", conf.int = T, 
                                       title = paste0("{EZH2-; JUN-; RELA+}"),
                                       ylab = paste0("Essentiality of EZH2 [Achilles score]"),
                                       xlab = paste0("log ratio expression of JUN and RELA [log2(TPM+1)]"),
                                       cor.coef = T, cor.method = "pearson")# + font("title", size = 10.5, face = "bold") + 

ggsave(paste0("Plots/Supp_Fig_4.png"),correlation_plot3,
       device = "png", width = 7, height = 5, units = "in", dpi = 300, bg = "white")


### Supplementary Figure 5 ###
ess1 <- Achilles["ENSG00000120875",] #DUSP4
exp_final <- CCLE_Exp["ENSG00000099860",] #GADD45B
cell_lines <- intersect(names(ess1), names(exp_final))

ess <-  ess1[cell_lines]
exp <- as.data.frame(exp_final[cell_lines])
data <- as.data.frame(cbind(ess, exp))

colnames(data) <- c("ess", "exp")


gene <- "DUSP4"
gMCS_genes_symbol <- "GADD45B"
correlation_plot<- ggpubr::ggscatter(data, x = "exp", y = 'ess', 
                                     add = "reg.line", conf.int = T, 
                                     # font.label = c(8), 
                                     # cor.coef.coord = c(3, 0.2),
                                     # title = expression("{ "*IQGAP[1]^"+" "),
                                     title = paste0("{DUSP4+; GADD45B-}"),
                                     # title = expression("IQGAP1^{+};GADD45B^{-}"),
                                     ylab = paste0("Essentiality of DUSP4 [Achilles score]"),
                                     xlab = paste0("Expression of GADD45B [log2(TPM+1)]"),
                                     cor.coef = T, cor.method = "pearson") #+ font("title", size = 10.5, face = "bold") + 

ggsave(paste0("Plots/Supp_Fig_5.png"),correlation_plot,
       device = "png", width = 7, height = 5, units = "in", dpi = 300, bg = "white")

### Supplementary Figure 6 ###
ess1 <- Achilles["ENSG00000108465",] #CDK5RAP3
exp_final <- CCLE_Exp["ENSG00000160691",] #SHC1
cell_lines <- intersect(names(ess1), names(exp_final))

ess <-  ess1[cell_lines]
exp <- as.data.frame(exp_final[cell_lines])
data <- as.data.frame(cbind(ess, exp))

colnames(data) <- c("ess", "exp")

correlation_plot2 <- ggpubr::ggscatter(data, x = "exp", y = 'ess', 
                                       add = "reg.line", conf.int = T, 
                                       # font.label = c(8), 
                                       cor.coef.coord = c(3.3, 0.2),
                                       title = paste0("{CDK5RAP3+; SHC1+}"),
                                       ylab = paste0("Essentiality of CDK5RAP3 [Achilles score]"),
                                       xlab = paste0("Expression of SHC1 [log2(TPM+1)]"),
                                       cor.coef = T, cor.method = "pearson")

ggsave(paste0("Plots/Supp_Fig_6.png"),correlation_plot2,
       device = "png", width = 7, height = 5, units = "in", dpi = 300, bg = "white")

gc()
