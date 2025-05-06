setwd("path/to/working/directory")
rm(list = ls())

library(openxlsx)
library(data.table)
library(org.Hs.eg.db)
library(readxl)
# library(tidyr)
library(clusterProfiler)
library(dplyr)

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

threshold <- fread("th_log2TPM_CCLE.txt")

databases <- excel_sheets("calculated_gMCSs_gMIS.xlsx")

mat_final_list <- list()
binary_genes <- list()
gmcs.pos <- list()

model <- databases
results <- list()
metrics_list <- list()

epsilon <- 1e-7
# databases <- databases[5]
# model <- databases
gMCSs.ENSEMBL.list <- list()

# databases <- databases[11]
for (model in databases){
  gMCSs <- read.xlsx("calculated_gMCSs_gMIS.xlsx", sheet = model)

  gMCSs.ENSEMBL <- as.matrix(gMCSs)
  rownames(gMCSs.ENSEMBL) <- as.character(seq(1,nrow(gMCSs.ENSEMBL)))
  gMCSs.ENSEMBL <- unique(gMCSs.ENSEMBL)
  gMCSs.ENSEMBL[is.na(gMCSs.ENSEMBL)] <- ""
  
  genes_in <- apply(gMCSs.ENSEMBL, 1, function(x){
    x <- gsub(pattern = "_off", replacement = "", x)
    x %in% CCLE_genes
  })
  
  genes_in <- colSums(genes_in) 
  
  lengths <- apply(gMCSs.ENSEMBL, 2, function(x){
    pos <- which(x == "")
    pos[length(pos)]
  })
  
  lengths <- do.call(cbind,lengths)
  lengths[length(lengths)+1] <- nrow(gMCSs.ENSEMBL)
  
  a <- 0
  keep <- array()
  for (i in 1:length(lengths)){
    pos <- lengths[i]
    keep[(a+1):(pos+a)] <- genes_in[(a+1):(pos+a)] == i
    a <- pos
  }
  
  gMCSs.ENSEMBL <- gMCSs.ENSEMBL[which(keep),]
  gMCSs.ENSEMBL.4.aux <- apply(gMCSs.ENSEMBL,1,function(x){paste0(x[x!=""],collapse = "--")})
  
  gMCSs.ENSEMBL.list[[model]] <- strsplit(gMCSs.ENSEMBL.4.aux, '--')

  CCLE_Exp[is.na(CCLE_Exp)] <- 0 
  genes_total <- unique(as.vector(gMCSs.ENSEMBL))
  mat_final <- matrix(0,nrow = length(genes_total), ncol = ncol(CCLE_Exp),
                      dimnames = list(genes_total, colnames(CCLE_Exp)))
  
  off_pos <- t(apply(gMCSs.ENSEMBL, 1, function(x){
    pos <- grepl("_off",x)
    pos[pos == T] <- 1
    pos[pos == F] <- 0
    empty <- !grepl("ENSG",x)
    pos[empty==T] <- -2
    empty <- x == ""
    pos[empty==T] <- -1
    return(pos)
  }))
  
  lengths <- apply(gMCSs.ENSEMBL, 1, function(x){
    sum(x >= 0)
  })

  ccle_cols <- colnames(CCLE_Exp)
  ccle_rows <- rownames(CCLE_Exp)
  
  gMCSs.ENSEMBL_DT <- as.data.table(gMCSs.ENSEMBL)
  model_result <- list(binary_genes = list(), gmcs_pos = list(), 
                       mat_update = matrix(0, nrow = length(genes_total), 
                                           ncol = ncol(CCLE_Exp),dimnames = list(genes_total, ccle_cols)))
  
  for (pos in 1:nrow(gMCSs.ENSEMBL_DT)){
    
    gene1 <- gMCSs.ENSEMBL_DT[pos, 1][[1]]
    gene2 <- gMCSs.ENSEMBL_DT[pos, 2][[1]]
    
    if (gene2 == "") {
      model_result$binary_genes[[gene1]] <- ccle_cols
      model_result$gmcs_pos[[gene1]] <- rep(pos, length(ccle_cols))
      model_result$mat_update[gene1,] <- rep(1, length(ccle_cols))
    } else {
      len_gmcs <- lengths[pos]
      gMCS <- unlist(gMCSs.ENSEMBL_DT[pos, 1:len_gmcs, with = FALSE])
      
      if (sum(off_pos[pos, ] == -2) < 1) {
        pos_on <- which(off_pos[pos, ] == 0)
        pos_off <- which(off_pos[pos, ] == 1)
        
        gMCS_off_clean <- gsub("_off", "", gMCS[pos_off])
        gMCS_on <- gMCS[pos_on]
        
        valid_off <- length(pos_off) == 0 || all(gMCS_off_clean %in% ccle_rows)
        valid_on <- length(pos_on) == 0 || all(gMCS_on %in% ccle_rows)
        
        if (valid_off && valid_on) {
          
          exp_off_binary <- if (length(pos_off) > 0) {
            exp_off <- CCLE_Exp[gMCS_off_clean, , drop = FALSE]
            rownames(exp_off) <- gMCS[pos_off]
            as.data.table(t(apply(exp_off, 1, function(x) (x - threshold$threshold) < epsilon)))
          } else NULL
          
          exp_on_binary <- if (length(pos_on) > 0) {
            exp_on <- CCLE_Exp[gMCS_on, threshold$V1, drop = FALSE]
            rownames(exp_on) <- gMCS[pos_on]
            as.data.table(t(apply(exp_on, 1, function(x) (x - threshold$threshold) > epsilon)))
          } else NULL
          
          if (!is.null(exp_on_binary) && !is.null(exp_off_binary)) {
            exp_on_binary[, gene := rownames(exp_on)]
            exp_off_binary[, gene := rownames(exp_off)]
            final_binary <- rbind(exp_on_binary, exp_off_binary)
            setkey(final_binary, gene)
          } else if (!is.null(exp_on_binary)) {
            final_binary <- copy(exp_on_binary)
            final_binary[, gene := rownames(exp_on)]
            setkey(final_binary, gene)
          } else if (!is.null(exp_off_binary)) {
            final_binary <- copy(exp_off_binary)
            final_binary[, gene := rownames(exp_off)]
            setkey(final_binary, gene)
          }
          final_binary <- as.data.frame(final_binary)
          rownames(final_binary) <- final_binary$gene
          pos_gene <- which(colnames(final_binary) == "gene")
          final_binary <- final_binary[,-pos_gene]
 
          one_hit <- which(colSums(final_binary) == 1)
          if (length(one_hit)>0){
            for (i in 1:length(one_hit)) {
              cell_line <- names(one_hit)[i]
              col <- one_hit[i]
              act_gene <- rownames(final_binary)[which(final_binary[,col])]
              
              model_result$binary_genes[[act_gene]] <- c(model_result$binary_genes[[act_gene]], cell_line)
              model_result$gmcs_pos[[act_gene]] <- c(model_result$gmcs_pos[[act_gene]], pos)
              model_result$mat_update[act_gene, cell_line] <- 1
            }
          }
          
        }
      }
    }
  }
  gc()
  
  if (file.exists(paste0("adaptation/Results_KO_", model, ".txt"))){
    KO_genes <- as.data.frame(fread(paste0("adaptation/Results_KO_", model, ".txt"), header = F))
    gMCSs_problem_genes <- fread(paste0("adaptation/Results_gMCSs_", model, ".txt"), header = F)
    KO_genes <- as.data.frame(KO_genes)
    gMCSs_problem_genes <- gMCSs_problem_genes
    
    gMCSs_problem_genes <- t(apply(gMCSs_problem_genes, 1, function(x){
      x <- gsub(pattern = "_auxKO_act", replacement = "_off", x)
      x <- gsub(pattern = "_auxKO", replacement = "", x)
    }))
    KO_genes <- as.data.frame(apply(KO_genes, 1, function(x){
      x <- gsub(pattern = "_auxKO_act", replacement = "_off", x)
      x <- gsub(pattern = "_auxKO", replacement = "", x)
    }))
    colnames(KO_genes) <- "V1"
    
    sdf <- reshape2::melt(model_result$mat_update)
    colnames(sdf) <- c("ENSEMBL", "CellLine", "Prediction")
    sdf$Prediction <- factor(as.character(sdf$Prediction), levels = c("1", "0"), labels = c("Essential", "Not essential"))
    
    ess_pred <- sdf[sdf$Prediction == "Essential",]
    ess_pred <- ess_pred[ess_pred$ENSEMBL %in% KO_genes$V1,]
    
    gMCSs_problem_genes <- apply(gMCSs_problem_genes,1,function(x){paste0(x[x!=""],collapse = "--")})
    
    KO_genes_final <- KO_genes[which(gMCSs_problem_genes %in% gMCSs.ENSEMBL.4.aux),]
    ess_pred <- ess_pred[ess_pred$ENSEMBL %in% KO_genes_final,]
    essential_genes <- unique(ess_pred[ess_pred$Prediction == "Essential","ENSEMBL"])
    
    if (!isEmpty(essential_genes)){
      for (i in 1:length(essential_genes)){
        list_pos <- model_result$gmcs_pos[as.character(essential_genes[i])]
        list_cell_line <- model_result$binary_genes[as.character(essential_genes[i])]
        
        pos_KO <- which(KO_genes$V1 %in% essential_genes[i])
        pos_pos <- which(gMCSs.ENSEMBL.4.aux %in% gMCSs_problem_genes[pos_KO])
        
        for (j in 1:length(list_pos[[1]])){
          if (!identical(list_pos[j], integer(0))){
            if (list_pos[[1]][j] %in% pos_pos){
              model_result$gmcs_pos[as.character(essential_genes[i])][[1]][j] <- 0
              model_result$binary_genes[as.character(essential_genes[i])][[1]][j] <- "NA"
            }
          }
        }
        model_result$gmcs_pos <- lapply(model_result$gmcs_pos, function(x) {x[x!=0]})
        model_result$binary_genes <- lapply(model_result$binary_genes, function(x) {x[x!="NA"]})
        
        if (identical(model_result$gmcs_pos[[as.character(essential_genes[i])]], numeric(0))){
          model_result$gmcs_pos[[as.character(essential_genes[i])]] <- NULL
          model_result$binary_genes[[as.character(essential_genes[i])]] <- NULL
          
          if (essential_genes[i] %in% rownames(model_result$mat_update)){
            pos_new <- which(rownames(model_result$mat_update) %in% essential_genes[i])
            model_result$mat_update[pos_new, ] <- 0
          }
        }
        if (essential_genes[i] %in% rownames(model_result$mat_update)){
          pos_new <- which( rownames(model_result$mat_update) %in% essential_genes[i])
          pos_essential <- which(model_result$mat_update[pos_new, ] == 1)
          essential_cell_lines <- names(pos_essential)[which(!names(pos_essential) %in% model_result$binary_genes[[as.character(essential_genes[i])]])]
          if (length(essential_cell_lines)>0){
            model_result$mat_update[pos_new, essential_cell_lines] <- 0
          }
        }
        
      }
      
    }
  }
  results[[model]] <- model_result
  
  genes_final <- intersect(rownames(Achilles), genes_total)
  cell_lines <- intersect(colnames(Achilles_binary), colnames(model_result$mat_update))
  Binary_Matrix_pred <- model_result$mat_update[genes_final,cell_lines]
  Binary_Matrix_DepMap <- Achilles_binary[genes_final, cell_lines]
  metrics = data.frame('cellLine' = cell_lines,
                       'TP' = NA, 'TN' = NA, 'FP' = NA, 'FN' = NA)

  for (x in 1:ncol(Binary_Matrix_pred)){

    genes_final_NArm <-  names(which(!is.na(Binary_Matrix_DepMap[,x])))
    
    Binary_Matrix_pred_NArm <- Binary_Matrix_pred[genes_final_NArm,]
    Binary_Matrix_DepMap_NArm <- Binary_Matrix_DepMap[genes_final_NArm, ]
    
    
    modelEssential <- rownames(Binary_Matrix_pred_NArm)[Binary_Matrix_pred_NArm[,x]>0]
    modelNonEssential <- setdiff(genes_final_NArm, modelEssential)
    
    DepMap_essential <- genes_final_NArm[Binary_Matrix_DepMap_NArm[,x]==1]
    DepMapNonEssential <- setdiff(genes_final_NArm,DepMap_essential)
   
    metrics$TP[x] = length(intersect(modelEssential, DepMap_essential))        
    metrics$TN[x] = length(intersect(modelNonEssential, DepMapNonEssential)) 
    metrics$FP[x] = length(intersect(modelEssential, DepMapNonEssential))     
    metrics$FN[x] = length(intersect(modelNonEssential, DepMap_essential))
    
  }
  metrics$PPV <-  metrics$TP/(metrics$TP+metrics$FP)
  metrics$Model <- model
  metrics_list[[model]] <- metrics
  gc()
  
}

metrics_df <- do.call(rbind, metrics_list)
SDL_results <- metrics_df

a <- aggregate(TP~Model,data=SDL_results,mean)
b <- aggregate(FP~Model,data=SDL_results,mean)
c <- aggregate(PPV~Model,data=SDL_results,mean)
g <-  aggregate(TN~Model,data=SDL_results,mean)
h <-  aggregate(FN~Model,data=SDL_results,mean)
final <- merge(merge(merge(merge(a,b),c),g),h)

final
# saveRDS(results, "results.RDS")

##Figure 3 ##
model_array <- c("Human1-O1", 
                 "Human1-O2",
                 "Human1-D1","Human1-D2",
                 "Human1-T1", "Human1-T2",
                 "Human1-S1", "Human1-S2", 
                 "Human1-O1-SDL", "Human1-O2-SDL", 
                 "Human1-D1-SDL", "Human1-D2-SDL", 
                 "Human1-T1-SDL","Human1-T2-SDL", 
                 "Human1-S1-SDL", "Human1-S2-SDL")

mat_Human1 <- results$Human1$mat_update
final_final <- data.frame()
final_list <- list()
final_list_all <- list()

for (j in 1:length(model_array)){
  model <- results[[model_array[j]]]
  mat <- model$mat_update
  gmcs.pos <- model$gmcs_pos
  binary_genes <- model$binary_genes
  genes <- rownames(mat)[rownames(mat) %in% rownames(mat_Human1)]
  genes <- genes[genes != ""]

  mat_mat <- mat[genes, ] + mat_Human1[genes,]
  
  mat[genes, ][mat_mat == 2] <- 0
  
  cell_lines <- colnames(mat)[colnames(mat) %in% colnames(Achilles_binary)]
  genes <- rownames(Achilles_binary)[rownames(Achilles_binary) %in% rownames(mat)]
  
  mat <- mat[!grepl("_off", rownames(mat)),]
  
  final <- list()
  final <- data.frame('1' = 0, '2' = 0, '3'= 0, '4'= 0, '5'=0)
  
  for (n in 1:ncol(mat)){
    pos <- lapply(binary_genes[names(which(mat[,n]==1))], function(x){
      which(x == colnames(mat)[n])
    })
    pos_in_list <- lapply(names(pos), function(name) {
      position <- pos[[name]]
      elements_in_position <- gmcs.pos[[name]][position]
      unique(lengths(gMCSs.ENSEMBL.list[[model_array[j]]][elements_in_position]))
    })
    if (length(table(unlist(pos_in_list)))<5){
      pre_final <- rep(0, 5)
      table_table <- table(unlist(pos_in_list))
      for (a in names(table_table)){
        pre_final[as.numeric(a)] <- table_table[a]
      }
    } else {
      pre_final <- table(unlist(pos_in_list))
    }
    final <- rbind(final, pre_final)
  }
  final <- final[-1,]
  
  rownames(final) <- colnames(mat)
  final_all <- as.data.frame(rowSums(final))
  final_all$model <- model_array[j]
  pre_mat <- reshape2::melt(as.matrix(final))
  pre_mat$model <- model_array[j]
  final_list[[model_array[j]]] <- pre_mat
  final_list_all[[model_array[j]]] <- final_all
  
  final_final[j,1:length(colMeans(final))] <- colMeans(final)
  
}
rownames(final_final) <- model_array
final_final_all <- as.data.frame(rowSums(final_final))
final_final_all$model <- model_array
colnames(final_final_all)[1] <- "value"

final_final_all$type <-  factor(ifelse(grepl("SDL", final_final_all$model), "SDL", "SL"), levels = c("SL", "SDL"))
gc()

final_to_boxplot_all <- do.call(rbind, final_list_all)
colnames(final_to_boxplot_all)[1] <- "value"

final_to_boxplot_all$type <- factor(ifelse(grepl("SDL", final_to_boxplot_all$model), "SDL", "SL"), levels = c("SL", "SDL"))
final_to_boxplot_all$model <- gsub("-SDL", "", final_to_boxplot_all$model)
final_to_boxplot_all$model <- factor(final_to_boxplot_all$model, 
                                     levels = c("Human1-O1", "Human1-O2", "Human1-O3",
                                                "Human1-D1", "Human1-D2","Human1-D3",
                                                                            "Human1-T1", "Human1-T2","Human1-T3",
                                                                            "Human1-S1", "Human1-S2", "Human1-S3"),
                                     labels = c("Human1-O1", "Human1-O2", "Human1-O3",
                                                "Human1-D1", "Human1-D2","Human1-D3",
                                                "Human1-T1", "Human1-T2","Human1-T3",
                                                "Human1-S1", "Human1-S2", "Human1-S3"))

colors <- RColorBrewer::brewer.pal(n = 10, name = "Paired")
colors_final <- c(colorRampPalette(c(colors[3], colors[4], "black"))(6)[4],
                  colorRampPalette(c(colors[3], colors[4], "black"))(6)[2],
                  colorRampPalette(c(colors[1], colors[2], "black"))(6)[4],
                  colorRampPalette(c(colors[1], colors[2], "black"))(6)[2],
                  colorRampPalette(c(colors[9], colors[10], "black"))(6)[4],
                  colorRampPalette(c(colors[9], colors[10], "black"))(6)[2],
                  colorRampPalette(c(colors[5], colors[6], "black"))(6)[4],
                  colorRampPalette(c(colors[5], colors[6], "black"))(6)[2])

data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

df3_all <- data_summary(final_to_boxplot_all, varname="value",
                        groupnames=c("type", "model"))

df3_all$model=as.factor(df3_all$model)


library(tidyverse)
library(RColorBrewer)
library(ggpubr)
library(rstatix)
library(ggpattern)

stats <- compare_means(value ~ type, group.by = c("model"), data = final_to_boxplot_all, method = "wilcox.test") 
stats$type <- "SL-SDL"

boxplot_all <- ggplot(df3_all, aes(x = model,y = value, fill = model, pattern = type))+
  geom_errorbar(aes(ymin=value-sd, ymax=value+sd), width=.2, position=position_dodge(.9))+
  theme_bw() +
  geom_bar_pattern(position = "dodge", stat = "identity",
                   color = "black", pattern_fill = "black", pattern_angle = 45, 
                   pattern_density = 0.1, pattern_spacing = 0.025, pattern_key_scale_factor = 0.6) +
  scale_pattern_manual(values = c(SL = "none", SDL = "stripe", "SL-SDL" = "stripe")) + 
  scale_fill_manual(values = colors_final, guide = guide_legend(override.aes = list(pattern = "none"), order = 1)) +
  scale_y_continuous(limits = c(0, NA),trans = "log1p", breaks = c(0,1,2,3,4,5,7,10,15,25,50,75,100,150,200,300,500)) +
  theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5), axis.text.x  = element_text(angle = 45, vjust = 1, hjust = 1),
        strip.text.x =element_blank(), axis.title.x = element_blank()) +
  ylab("New essential KO genes") +


ggtitle("New essential KO genes") 

##PPV 
Melted_gMCS_All_Statistics_th <- reshape2::melt(data = SDL_results,
                                      id.vars = c("Model", "cellLine"),
                                      measure.vars = c("TP", "FP",
                                                       "PPV"))

names2 <- c("Human1",
            "Human1-O1", "Human1-O1-SDL",
            "Human1-O2","Human1-O2-SDL", 
            "Human1-D1","Human1-D1-SDL",
            "Human1-D2","Human1-D2-SDL",
            "Human1-T1", "Human1-T1-SDL",
            "Human1-T2","Human1-T2-SDL",
            "Human1-S1", "Human1-S1-SDL", 
            "Human1-S2","Human1-S2-SDL"  
)

Melted_gMCS_All_Statistics_th <- Melted_gMCS_All_Statistics_th[Melted_gMCS_All_Statistics_th$Model %in% names2,]
Melted_gMCS_All_Statistics_th$Model <- factor(Melted_gMCS_All_Statistics_th$Model, names2)
Melted_gMCS_All_Statistics_th$Type <- ifelse(grepl(pattern = "-SDL", x = Melted_gMCS_All_Statistics_th$Model), "SDL", "SL")
Melted_gMCS_All_Statistics_th$Type <- factor(Melted_gMCS_All_Statistics_th$Type, levels = c("SL", "SDL"))
Melted_gMCS_All_Statistics_th$Model <- sub("-SDL", "",Melted_gMCS_All_Statistics_th$Model )
names2 <- c("Human1",
            "Human1-O1",
            "Human1-O2",
            "Human1-D1","Human1-D2",
            "Human1-T1",  "Human1-T2",
            "Human1-S1", 
            "Human1-S2")

Melted_gMCS_All_Statistics_th$Model <- factor(Melted_gMCS_All_Statistics_th$Model, names2)

levels(Melted_gMCS_All_Statistics_th$Model) <- c("Human1",
                                                 "Human1-O1",
                                                 "Human1-O2", 
                                                 "Human1-D1",
                                                 "Human1-D2", 
                                                 "Human1-T1",  
                                                 "Human1-T2",
                                                 "Human1-S1", 
                                                 "Human1-S2")
colors_final <- c(colors[8],
                  colorRampPalette(c(colors[3], colors[4], "black"))(6)[4],
                  colorRampPalette(c(colors[3], colors[4], "black"))(6)[2],
                  colorRampPalette(c(colors[1], colors[2], "black"))(6)[4],
                  colorRampPalette(c(colors[1], colors[2], "black"))(6)[2],
                  colorRampPalette(c(colors[9], colors[10], "black"))(6)[4],
                  colorRampPalette(c(colors[9], colors[10], "black"))(6)[2],
                  colorRampPalette(c(colors[5], colors[6], "black"))(6)[4],
                  colorRampPalette(c(colors[5], colors[6], "black"))(6)[2])

PPV <- ggplot(Melted_gMCS_All_Statistics_th[Melted_gMCS_All_Statistics_th$variable %in% "PPV",], aes(x = Model,y = value, fill = Model, pattern = Type))+
  geom_boxplot(aes(fill = Model)) +  #scale_fill_grey()+
  theme_bw() + 
 geom_boxplot_pattern(#position = position_dodge(preserve = "single"), 
    color = "black", pattern_fill = "black", pattern_angle = 45, pattern_density = 0.1, pattern_spacing = 0.025, pattern_key_scale_factor = 0.6) +
  scale_pattern_manual(values = c(SL = "none", SDL = "stripe")) + 
  scale_fill_manual(values = colors_final, guide = guide_legend(override.aes = list(pattern = "none"), order = 1)) +
  theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        strip.text.x =element_blank(), axis.title.x = element_blank(), axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1) ) + 
  ylab("PPV") +
  ggtitle("Positive Predictive Value (PPV)") 


#EZH2, JUN, RELA
genes <- c("ENSG00000106462", "ENSG00000177606", "ENSG00000173039")

CCLE_Exp <- t(CCLE_Exp)
cells <- intersect(rownames(CCLE_Exp), colnames(Achilles))
CCLE_Exp <- CCLE_Exp[cells,]
Achilles <- Achilles[,cells]

data <- data.frame(exp = (CCLE_Exp[,genes[2]]-CCLE_Exp[,genes[3]]), ess = Achilles[genes[1],])
expr_EHZ2 <- CCLE_Exp[, "ENSG00000106462"]
data$KO_level = "Q2-Q3"
data$KO_level[expr_EHZ2 >= quantile(expr_EHZ2,0.75)] = "Q4"
data$KO_level[expr_EHZ2 <= quantile(expr_EHZ2,0.25)] = "Q1"
data$KO_level <- factor(data$KO_level, levels = c("Q1", "Q2-Q3", "Q4"))
data$log_expr_EHZ2 <- expr_EHZ2

colnames(data) <- c("exp", "ess", "KO_level", "log_exp_EHZ2")

correlation_plot2 <- ggarrange(ggscatter(data,
                                         x = "exp", y = "ess", color = "KO_level",
                                         facet.by = "KO_level") + 
                                 geom_smooth(method = "lm", formula = y ~ x, col = "red") +
                                 scale_y_continuous(expand = c(0.1,0,0.1,0)) + 
                                 stat_cor(method="pearson", label.x.npc="right", label.y.npc="top", hjust=1, vjust = 0) + 
                                 geom_hline(yintercept = c(0, -0.6), colour = "grey", alpha = 0.9, linetype = "dashed") + 
                                 xlab(paste0("log ratio expression of JUN and RELA [log2(TPM+1)]")) +
                                 ylab(paste0("EZH2 Achilles score")) + 
                                 facet_wrap(~KO_level, nrow = 1) + ggtitle("EZH2-, JUN-, RELA+"),
                               ggviolin(data, x = "KO_level", y = "log_exp_EHZ2", fill = "KO_level",
                                        add = "jitter", add.params = list(alpha = 0.5, fill = "KO_level")) + 
                                 ylab("log2[TPM+1]") + 
                                 ggtitle("EHZ2", subtitle = "essential gene"),
                               nrow = 1, common.legend = T, legend = "right", widths = c(3,1))

correlation_plot<- ggpubr::ggscatter(data, x = "exp", y = 'ess',
                                     add = "reg.line",
                                     conf.int = T,
                                     cor.coeff.args = list(label.y = 0.6),
                                     cor.coef = T, cor.method = "pearson",
                                     xlab = paste0("log ratio expression of JUN and RELA [log2(TPM+1)]"), 
                                     ylab = paste0("Essentiality of EZH2 [Achilles score]"),title = "{EZH2-; JUN-; RELA+}")

plot_fig3 <- ggarrange(ggarrange(ggarrange(boxplot_all, NULL, widths = c(1,0), legend = F),PPV, common.legend = T,
                             legend = "bottom",
                             labels = c("A", "B")), correlation_plot2, nrow = 2, common.legend = T,labels = c("", "C"), heights = c(1,0.8))

ggsave(paste0("Plots/Figure3.pdf"),
       plot = plot_fig3, device = "pdf", width = 12, height = 10, units = "in", dpi = 300, bg = "white")


##Figure 4
##########################
#Dusp4+ & GADD45B-
gmcs <- c("ENSG00000099860", "ENSG00000160691_off")

ess1 <- Achilles["ENSG00000120875",]
exp_final <- CCLE_Exp[,"ENSG00000099860"]
cell_lines <- intersect(names(ess1), names(exp_final))



ess <-  as.data.frame(ess1[cell_lines])
exp <- as.data.frame(exp_final[cell_lines])
data <- as.data.frame(cbind(ess, exp))

expr_DUSP7 <- CCLE_Exp[cell_lines,"ENSG00000120875"]
data$KO_level = "Q2-Q3"
data$KO_level[expr_DUSP7 >= quantile(expr_DUSP7,0.75)] = "Q4"
data$KO_level[expr_DUSP7 <= quantile(expr_DUSP7,0.25)] = "Q1"
data$KO_level <- factor(data$KO_level, levels = c("Q1", "Q2-Q3", "Q4"))

data$log_exp_DUSP7 <- expr_DUSP7

colnames(data) <- c("ess", "exp", "KI_level", "log_exp_DUSP4")

correlation_plot <- ggarrange(ggscatter(data,
                                        x = "exp", y = "ess", color = "KI_level",
                                        facet.by = "KI_level") + 
                                geom_smooth(method = "lm", formula = y ~ x, col = "red") +
                                scale_y_continuous(expand = c(0.1,0,0.1,0)) + 
                                stat_cor(method="pearson", label.x.npc="right", label.y.npc="top", hjust=1, vjust = 0) + 
                                geom_hline(yintercept = c(0, -0.6), colour = "grey", alpha = 0.9, linetype = "dashed") + 
                                xlab(paste0("GADD45B expression [log2(TPM+1)]")) + 
                                ylab(paste0("DUSP4 Achilles score")) + 
                                facet_wrap(~KI_level, nrow = 1) + ggtitle("DUSP4+, GADD45B-"),

                              ggviolin(data, x = "KI_level", y = "log_exp_DUSP4", fill = "KI_level",
                                       add = "jitter", add.params = list(alpha = 0.5, fill = "KI_level")) + ylab("log2[TPM+1]") + 
                                ggtitle("DUSP4", subtitle = "tumor suppressor"),
                              
                              nrow = 1, common.legend = T, legend = "right", widths = c(3,1,1,1))

##### 
#CDK5RAP3+ & SHC1+
##########

ess1 <- Achilles["ENSG00000108465",]
exp_final <- CCLE_Exp[,"ENSG00000160691"]
cell_lines <- intersect(names(ess1), names(exp_final))

ess <-  as.data.frame(ess1[cell_lines])
exp <- as.data.frame(exp_final[cell_lines])
data <- as.data.frame(cbind(ess, exp))


expr_DUSP7 <- CCLE_Exp[cell_lines, "ENSG00000108465"]
data$KO_level = "Q2-Q3"
data$KO_level[expr_DUSP7 >= quantile(expr_DUSP7,0.75)] = "Q4"
data$KO_level[expr_DUSP7 <= quantile(expr_DUSP7,0.25)] = "Q1"
data$KO_level <- factor(data$KO_level, levels = c("Q1", "Q2-Q3", "Q4"))

data$log_exp_DUSP7 <- expr_DUSP7

colnames(data) <- c("ess", "exp", "KI_level", "log_exp_CDK5RAP3")

correlation_plot2 <- ggarrange(ggscatter(data,
                                         x = "exp", y = "ess", color = "KI_level",
                                         facet.by = "KI_level") + 
                                 geom_smooth(method = "lm", formula = y ~ x, col = "red") +
                                 scale_y_continuous(expand = c(0.1,0,0.1,0)) + 
                                 stat_cor(method="pearson", label.x.npc="right", label.y.npc="top", hjust=1, vjust = 0) + 
                                 geom_hline(yintercept = c(0, -0.6), colour = "grey", alpha = 0.9, linetype = "dashed") + 
                                 xlab(paste0("SHC1 expression [log2(TPM+1)]")) +
                                 ylab(paste0("CDK5RAP3 Achilles score")) + 
                                 facet_wrap(~KI_level, nrow = 1) + ggtitle("CDK5RAP3+, SHC1+"),
                               ggviolin(data, x = "KI_level", y = "log_exp_CDK5RAP3", fill = "KI_level",
                                        add = "jitter", add.params = list(alpha = 0.5, fill = "KI_level")) + ylab("log2[TPM+1]") + 
                                 ggtitle("CDK5RAP3", subtitle = "tumor suppressor"),
                               nrow = 1, common.legend = T, legend = "right", widths = c(3,1,1,1))

model_array <- c("Human1-O1-SDL", "Human1-O2-SDL", 
                 "Human1-D1-SDL", "Human1-D2-SDL", 
                 "Human1-T1-SDL","Human1-T2-SDL", 
                 "Human1-S1-SDL", "Human1-S2-SDL")
final_final <- data.frame()
final_list <- list()
for (j in 1:length(model_array)){
  model <- model_array[j]
  data <- results[[model]]
  mat <- data$mat_update
  gmcs.pos <- data$gmcs_pos
  binary_genes <- data$binary_genes

  cell_lines <- colnames(mat)[colnames(mat) %in% colnames(Achilles_binary)]
  genes <- rownames(Achilles_binary)[rownames(Achilles_binary) %in% rownames(mat)]

  mat <- mat[grepl("_off", rownames(mat)),]
  
  final <- list()
  final <- data.frame('1' = 0, '2' = 0, '3'= 0, '4'= 0, '5'=0)
  
  for (n in 1:ncol(mat)){
    pos <- lapply(binary_genes[names(which(mat[,n]==1))], function(x){
      which(x == colnames(mat)[n])
    })

    positions_in_gMCs <- lapply(names(pos), function(name) {
      position <- pos[[name]]
      elements_in_position <- gmcs.pos[[name]][position]
      unique(lengths(gMCSs.ENSEMBL.list[[model]][elements_in_position]))

    })

    if (length(table(unlist(positions_in_gMCs)))<5 ){
      pre_final <- rep(0, 5)
      table_table <- table(unlist(positions_in_gMCs))
      for (a in names(table_table)){
        pre_final[as.numeric(a)] <- table_table[a]
      }

    } else {
      pre_final <- table(unlist(posiciones_en_lista2))
    }
    final <- rbind(final, pre_final)
  }
  final <- final[-1,]
  pre_mat <- reshape2::melt(as.matrix(final))
  pre_mat$model <- model
  final_list[[model]] <- pre_mat
  
  final_final[j,1:length(colMeans(final))] <- colMeans(final)
}

rownames(final_final) <- model_array

gc()

final_to_boxplot <- do.call(rbind, final_list)

final_to_boxplot$type <- factor(ifelse(grepl("SDL", final_to_boxplot$model), "SDL", "SL"), levels = c("SL", "SDL"))
final_to_boxplot$model <- gsub("-SDL", "", final_to_boxplot$model)
final_to_boxplot$model <- factor(final_to_boxplot$model, 
                                 levels = c("Human1-O1", "Human1-O2",
                                            "Human1-D1", "Human1-D2",
                                            "Human1-T1", "Human1-T2",
                                            "Human1-S1", "Human1-S2"), 
                                 labels = c("Human1-O1", "Human1-O2", 
                                            "Human1-D1", "Human1-D2",
                                            "Human1-T1", "Human1-T2",
                                            "Human1-S1", "Human1-S2"))

colors <- RColorBrewer::brewer.pal(n = 10, name = "Paired")
colors_final <- c(colorRampPalette(c(colors[3], colors[4], "black"))(6)[4],
                  
                  colorRampPalette(c(colors[3], colors[4], "black"))(6)[2],
                  colorRampPalette(c(colors[1], colors[2], "black"))(6)[4],
                  
                  colorRampPalette(c(colors[1], colors[2], "black"))(6)[2],
                  colorRampPalette(c(colors[9], colors[10], "black"))(6)[4],
                  
                  colorRampPalette(c(colors[9], colors[10], "black"))(6)[2],
                  colorRampPalette(c(colors[5], colors[6], "black"))(6)[4],
                  colorRampPalette(c(colors[5], colors[6], "black"))(6)[2])

data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}
df3 <- data_summary(final_to_boxplot, varname="value", 
                    groupnames=c("Var2","type", "model"))

df3$model=as.factor(df3$model)
# head(df3)


barplot <- ggplot(df3, aes(x = Var2,y = value, fill = model, pattern = type)) +
  geom_errorbar(aes(ymin=value-sd, ymax=value+sd), width=.2, position=position_dodge(.9)) +
  theme_bw() +
  geom_bar_pattern(position = "dodge", stat = "identity",
                   color = "black", pattern_fill = "black", pattern_angle = 45, pattern_density = 0.1, pattern_spacing = 0.025, pattern_key_scale_factor = 0.6) +
  geom_vline(xintercept = c(1.5, 2.5, 3.5, 4.5), linetype = "dashed", alpha = 0.5) +
  scale_pattern_manual(values = c(SL = "none", SDL = "stripe")) + 
  scale_fill_manual(values = colors_final, guide = guide_legend(override.aes = list(pattern = "none"), order = 1)) +
  scale_x_discrete(labels=c("X1" = "length = 1", "X2" = "length = 2", "X3" = "length = 3", "X4" = "length = 4", "X5" = "length = 5")) +
  scale_y_continuous(limits = c(0, NA), trans = "log1p",breaks = c(0,1,2,3,4,5,7,10,15,25, 35,50,75,100,150,200,300,500)) +
  theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        strip.text.x =element_blank(), axis.title.x = element_blank(), 
        axis.text.y = element_text( size = 12), 
        axis.text.x = element_text( size = 12),
        axis.title.y = element_text( size = 12),
        legend.position = "top", panel.border = element_blank(), axis.line = element_line(colour = "black")) +
  ylab("KI genes")


plotfig4 <- ggarrange(correlation_plot, correlation_plot2, barplot, labels = c("A", "B", "C"),
                   nrow = 3,  heights = c(0.8,0.8, 1))

ggsave(paste0("Plots/Figure4.pdf"),plotfig4,
       device = "pdf", width = 13, height = 10, units = "in", dpi = 300, bg = "white")
