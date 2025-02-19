# set working dir
working_dir <- ""
if (!dir.exists(working_dir)) {
  dir.create(working_dir, recursive = TRUE)
}
setwd(working_dir)

#######COE####
library(dplyr)
library(tidyverse)
data_ori <- read.csv("", row.names=1)
data <- data_ori
# subdir
subdir <- "COE"
subdir_path <- file.path(working_dir, subdir)
if (!dir.exists(subdir_path)) {
  dir.create(subdir_path)
}

groupnames <- colnames(data)

groups <- gsub("(.*?)_\\d+", "\\1", groupnames)

colnames(data) <- groups


groups_names <- unique(groups)

frequency <- as.data.frame(table(groups))


num_rows <- nrow(frequency)

for (i in 1:num_rows) {

  current_column <- frequency[i, 2]
  

  assign(paste("df_", i, sep = ""), data[, (4*(i-1)+1):(4*(i-1)+current_column)])
}

COE <- data.frame(matrix(ncol = ncol(data), nrow = nrow(data)))



k <-num_rows


for (j in 1:k) {

  df_name <- paste0("df_", j)
  

  df_data <- get(df_name)
  

  
  num_columns <- ncol(df_data)

  for (i in 1:num_columns) {
    
    data_df <- df_data %>%
      filter(rowSums(!is.na(.)) == i)
    

    result_df <- apply(data_df, 1, function(x) sd(x, na.rm = TRUE)/mean(x, na.rm = TRUE)*100)  
    
    COE[,(4*j+i-4)] <- rep(c(result_df, rep(NA, nrow(COE) - length(result_df))), length.out = nrow(COE))
    
    
  }
}


groupnames_rep <- paste0(groupnames,"_rep")

colnames(COE) <- groupnames_rep



non_nas <- colSums(!is.na(COE)) > 0

COE <- COE[, non_nas, drop = FALSE]



cv_table_df  <- gather(COE, key=sample, value="quantification_COE")

cv_table_df$sample <- factor(cv_table_df$sample,
                             levels = groupnames_rep)  

p=ggplot(cv_table_df, aes(x = sample, y = quantification_COE, fill = sample)) +
  geom_boxplot() +
  theme(legend.position = "top")+scale_y_continuous(breaks = seq(0, 150, 20))

output_file_png <- file.path(subdir_path, "COE_Proteomics.png")
output_file_CSV <- file.path(subdir_path, "COE_Proteomics.csv")

ggsave(filename = output_file_png, plot = p, height = 6, width = 12, dpi = 300)

write.csv(COE, file = output_file_CSV, row.names = FALSE)

objects_to_save <- c("data_ori","working_dir")
rm(list = setdiff(ls(), objects_to_save))

#######box_plot_log10###
# Load the necessary packages
library(ggplot2)
library(tidyverse)

subdir <- "quality conrol"
subdir_path <- file.path(working_dir, subdir)
if (!dir.exists(subdir_path)) {
  dir.create(subdir_path)
}

output_folder <- subdir_path


output_file  <- file.path(output_folder, "log10_distributuion.png")

data <- data_ori


groupnames <- colnames(data)


log10_transformed <- log10(data)


log10_transformed  <- gather(log10_transformed , key=sample, value="log10_quantification")

log10_transformed$sample <- factor( log10_transformed$sample,
                                    levels = groupnames)  

p=ggplot(log10_transformed, aes(x = sample, y = log10_quantification, fill = sample)) +
  geom_boxplot() +
  theme(legend.position = "top")

ggsave(filename = output_file, plot = p, height = 8, width = 10, dpi = 300)

objects_to_save <- c("data_ori","working_dir")

rm(list = setdiff(ls(), objects_to_save))

########t_test#####################################
# Load the data from the CSV file
proteomics_data <- data_ori

# Perform a log2 transform on the data
proteomics_data_log2 <- log2(proteomics_data)

# Get the unique group names
group_names <- unique(gsub("_[0-9]+$", "", colnames(proteomics_data_log2)))

subdir <- "twobytwo"
subdir_path <- file.path(working_dir, subdir)
if (!dir.exists(subdir_path)) {
  dir.create(subdir_path)
}

output_folder <- subdir_path

# Set the output directory path
output_dir <- output_folder

# Check if the directory exists, and create it if it doesn't
if (!file.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Loop through each pair of groups and save them to a separate CSV file in the output directory
for (i in 1:(length(group_names)-1)) {
  for (j in (i+1):length(group_names)) {
    # Get the names of the current pair of groups
    group1_name <- group_names[i]
    group2_name <- group_names[j]
    
    # Subset the data for the current pair of groups
    group1_data_log2 <- proteomics_data_log2[, grep(paste0("^", group1_name, "_"), colnames(proteomics_data_log2))]
    group2_data_log2 <- proteomics_data_log2[, grep(paste0("^", group2_name, "_"), colnames(proteomics_data_log2))]
    
    # Combine the pair of groups into a single data frame
    output_data <- cbind(group1_data_log2, group2_data_log2)
    
    # Add row names as a column to the data frame
    output_data$Protein_ID <- rownames(output_data)
    
    # Save the pair of groups to a CSV file in the output directory, including row names
    write.csv(output_data, file.path(output_dir, paste0(group1_name, "_vs_", group2_name, "_log2.csv")), row.names = FALSE)
  }
}
objects_to_save <- c("subdir_path","working_dir")
rm(list = setdiff(ls(), objects_to_save))   

input_dir <- subdir_path

# 创建子目录
subdir_clean <- "twobytwo_clean"
subdir_clean_path <- file.path(working_dir, subdir_clean)
if (!dir.exists(subdir_path)) {
  dir.create(subdir_path)
}

output_folder <- subdir_clean_path

output_dir <- output_folder

# Create the output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# Get a list of all CSV files in the input directory
csv_files <- list.files(input_dir, pattern = "\\.csv$", full.names = TRUE)

# Loop through each CSV file and clean the data
for (csv_file in csv_files) {
  # Load the CSV file
  csv_data <- read.csv(csv_file, header = TRUE, check.names = FALSE)
  
  # 删除带有空值的行
  csv_data_cleaned <- csv_data[complete.cases(csv_data), ]
  
  # Get the row names from the last column of the data
  row_names <- csv_data_cleaned[, ncol(csv_data_cleaned)]
  
  # Remove the last column from the data
  csv_data_cleaned <- csv_data_cleaned[, -ncol(csv_data_cleaned)]
  
  # Assign the row names to the data frame
  rownames(csv_data_cleaned) <- row_names
  
  
  # Save the cleaned data to a new CSV file in the output directory, including row names
  output_file <- file.path(output_dir, basename(csv_file))
  write.csv(csv_data_cleaned, file = output_file, row.names = TRUE)
}

objects_to_save <- c("subdir_clean_path","working_dir")
rm(list = setdiff(ls(), objects_to_save))   


input_folder <- subdir_clean_path
subdir <- "t_test"
subdir_path_t_test <- file.path(working_dir, subdir)
if (!dir.exists(subdir_path_t_test)) {
  dir.create(subdir_path_t_test)
}


output_folder <- subdir_path_t_test

# Create the output directory if it doesn't exist
if (!dir.exists(output_folder )) {
  dir.create(output_folder)
}

file_paths <- list.files(path = input_folder, pattern = "\\.csv$", full.names = TRUE)


for (file_path in file_paths) {

  proteomics.csv <- read.csv(file_path)
  
  
  
  new_first_row <- colnames(proteomics.csv) # Copy the column names from the original data
  new_first_row[1] <- "Protein_ID" # Replace the first column name
  new_first_row <- gsub("_\\d+$", "", new_first_row) # Remove the _1, _2, _3 suffixes
  colnames(proteomics.csv) <- new_first_row # Replace the column names with the modified first row
  
  

  protein_names <- proteomics.csv[,1]
  proteomics.csv <- proteomics.csv[,2:9]
  

  group1 <- proteomics.csv[,1:4]
  group2 <- proteomics.csv[,5:8]

  group_1_title <- names(group1)[1]
  group_2_title <- names(group2)[1]
  

  group1_mean <- rowMeans(group1)
  group2_mean <- rowMeans(group2)
  
  results <- data.frame(Protein_ID = protein_names, stringsAsFactors = FALSE)
  for (i in 1:nrow(proteomics.csv)) {
    result <- t.test(group1[i,], group2[i,])
    results[i, "t_statistic"] <- result$statistic
    results[i, "p_value"] <- result$p.value
    results[i, group_1_title] <- group1_mean[i]
    results[i, group_2_title] <- group2_mean[i]
    results[i, "difference"] <- group2_mean[i] - group1_mean[i]
    results[i, "minus_log10_p_value"] <- -log10(result$p.value)
  }
  

  rownames(results) <- results$Protein_ID
  results$Protein_ID <- NULL

  results$Protein_ID <- rownames(results)
  

  output_file_name <- paste0("t_testresults_", gsub(".csv", "", basename(file_path)), ".csv")
  output_file_path <- file.path(output_folder, output_file_name)
  write.csv(results, file = output_file_path, row.names = TRUE)
}


objects_to_save <- c("subdir_path_t_test","working_dir")

rm(list = setdiff(ls(), objects_to_save))   

library(ggplot2)

input_folder <- subdir_path_t_test

subdir <- "t_test_vol"
subdir_path_t_test_vol <- file.path(working_dir, subdir)
if (!dir.exists(subdir_path_t_test_vol)) {
  dir.create(subdir_path_t_test_vol)
}

output_folder <- subdir_path_t_test_vol


if (!dir.exists(output_folder )) {
  dir.create(output_folder)
}  


dir.create(output_folder, showWarnings = FALSE, recursive = TRUE)


file_list <- list.files(input_folder, pattern = "\\.csv$")

for (file in file_list) {

  data <- read.csv(file.path(input_folder, file))
  

  filename <- tools::file_path_sans_ext(file)
  output_file <- file.path(output_folder, paste0(filename, ".png"))
  

  group_1_title <- names(data)[4]
  group_2_title <- names(data)[5]
  
  #Protein_up<-subset(Protein,significant=="TRUE" & log2FoldChange>=2)
  Protein_up<-subset(data,p_value <0.05 & difference>=log2(2))
  Protein_up_file <- cbind(rownames(Protein_up),Protein_up)
  
  
  #Protein_down<-subset(Protein,significant=="TRUE" & log2FoldChange<=(-2))
  Protein_down<-subset(data, p_value <0.05 & difference<=-log2(2)) 
  Protein_down_file <- cbind(rownames(Protein_down),Protein_down)
  
  #Protein_none<-subset(Protein,significant=="FALSE")
  Protein_none<-subset(data,p_value>0.05 | abs(difference) <log2(2))
  Protein_none_file <- cbind(rownames(Protein_none),Protein_none)
  

  up<-dim(Protein_up)[1]
  down<-dim(Protein_down)[1]
  total<-up+down
  uplable=paste("up :",up)
  downlable=paste("down :",down)
  Protein_up$sig<-uplable
  Protein_down$sig<-downlable
  Protein_none$sig<-"FALSE"
  data<-rbind(Protein_up,Protein_down,Protein_none)
  

  p<-ggplot(data)+ 
    geom_point(aes(x=difference,y=minus_log10_p_value,color=sig),size=0.8)+
    ggtitle(paste(group_1_title, "vs", group_2_title))+ theme(plot.title = element_text(hjust = 0.5)) +
    xlab(bquote(paste(log[2],"(fold change)",sep="")))+
    labs(x = "log2(Fold Change)", y = "-log10(p-value)") 
  
  p<-p+ scale_color_manual(paste("Protein_difference","(",total,")"),breaks=c(uplable,downlable,NA),values=c("#d71345","#7fb80e","#77787b"))
  p<-p+geom_hline(yintercept=-log10(0.05),linetype="dotdash",size=0.5)+geom_vline(xintercept=c(log2(2),-log2(2)),linetype="dotdash",size=0.5)
  

  p<-p+theme(panel.border=element_rect(fill=NA,colour="black"))
  p<-p+theme(
    panel.background = element_rect(fill = "transparent",colour =NA),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    plot.background = element_rect(fill  = "transparent",colour =NA)
  )
  

  ggsave(filename = output_file, plot = p, height = 6, width = 6, dpi = 1000)
}

##############t-test end############

objects_to_save <- c("subdir_path_t_test","working_dir")

rm(list = setdiff(ls(), objects_to_save))    


####Pathway analysis########

#####GO analysis#######
# Set the directory containing the CSV files
csv_dir <- subdir_path_t_test


subdir <- "higher_lower"
subdir_path_higher_lower <- file.path(working_dir, subdir)
if (!dir.exists( subdir_path_higher_lower)) {
  dir.create( subdir_path_higher_lower)
}

output_folder <-  subdir_path_higher_lower


# Get a list of all CSV files in the directory
csv_files <- list.files(path = csv_dir, pattern = "*.csv", full.names = TRUE)

# Loop through each CSV file
for (csv_file in csv_files) {
  # Read in the CSV file
  protein_data <- read.csv(csv_file, header = TRUE)
  
  # Rename the first column to "Protein_ID"
  colnames(protein_data)[1] <- "Protein_ID"
  
  # Filter for proteins with higher and lower abundance
  higher_abundance_list <- protein_data$Protein_ID[protein_data$p_value < 0.05 & protein_data$difference > log2(1.5)]
  lower_abundance_list <- protein_data$Protein_ID[protein_data$p_value < 0.05 & protein_data$difference < -log2(1.5)]
  
  # Output results to a new directory with the same filename as the input file
  output_dir <- output_folder
  output_file_higher <- paste0(output_dir, "/higher_abundance_", tools::file_path_sans_ext(basename(csv_file)), ".txt")
  output_file_lower <- paste0(output_dir, "/lower_abundance_", tools::file_path_sans_ext(basename(csv_file)), ".txt")
  
  # Check if there are any proteins with higher abundance
  if (length(higher_abundance_list) == 0) {
    cat(paste0("No proteins with higher abundance found in ", csv_file, "\n"))
  } else {
    # Output the list of proteins with higher abundance to a new file
    write.table(higher_abundance_list, file = output_file_higher, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
    cat(paste0("Proteins with higher abundance in ", csv_file, " written to ", output_file_higher, "\n"))
  }
  
  # Check if there are any proteins with lower abundance
  if (length(lower_abundance_list) == 0) {
    cat(paste0("No proteins with lower abundance found in ", csv_file, "\n"))
  } else {
    # Output the list of proteins with lower abundance to a new file
    write.table(lower_abundance_list, file = output_file_lower, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
    cat(paste0("Proteins with lower abundance in ", csv_file, " written to ", output_file_lower, "\n"))
  }
}


library(clusterProfiler)
library(org.Hs.eg.db)

# Set input and output paths
input_path <- output_folder

subdir <- "GO_BP"
subdir_path_GO_BP <- file.path(working_dir, subdir)
if (!dir.exists( subdir_path_GO_BP)) {
  dir.create( subdir_path_GO_BP)
}
output_folder <-  subdir_path_GO_BP


# Get list of text files in input folder
file_list <- list.files(input_path, pattern = "\\.txt$")

# Loop through each file and perform GO analysis
for (file in file_list) {
  # Load protein list using read.table
  uniprot_list <- read.table(file.path(input_path, file), header = FALSE, stringsAsFactors = FALSE)
  # Convert data.frame to vector
  uniprot_list <- as.vector(uniprot_list$V1)
  
  # Convert UniProt IDs to Entrez gene IDs
  entrez_ids <- bitr(uniprot_list,
                     fromType = "UNIPROT",
                     toType = "ENTREZID",
                     OrgDb = org.Hs.eg.db)
  
  # Remove any NA values introduced during conversion (if your IDs contain unrecognized ID mappings)
  entrez_ids <- entrez_ids[complete.cases(entrez_ids),]
  
  # Perform gene ontology analysis
  go_results_bp <- enrichGO(gene = entrez_ids[, 2], # Use the column with Entrez IDs 
                            OrgDb = org.Hs.eg.db,
                            keyType = "ENTREZID",
                            ont = "BP",
                            pAdjustMethod = "BH",
                            pvalueCutoff = 0.05,)
  
  # Set output file paths
  go_output_file <- file.path(output_folder, paste0(gsub("\\.txt$", "", file), "_go_results_bp.csv"))
  
  # Write GO results to CSV file
  write.csv(as.data.frame(go_results_bp), go_output_file, row.names = FALSE)
  
  
}  

library(ggplot2)


input_path <- output_folder

subdir <- "GO_BP_plot/"
subdir_path_GO_BP_plot <- file.path(working_dir, subdir)
if (!dir.exists( subdir_path_GO_BP_plot)) {
  dir.create( subdir_path_GO_BP_plot)
}

output_folder <-  subdir_path_GO_BP_plot

output_path <- paste0(output_folder,"/")


file_list <- list.files(input_path, pattern = "\\.csv$")


for (file in file_list) {
 
  go_results <- read.csv(file.path(input_path, file), header = TRUE, stringsAsFactors = FALSE)
  

  go_results_top10 <- head(go_results[order(go_results$pvalue), ], 10)
  

  base_size <- 12
  

  go_plot <- ggplot(data = go_results_top10, aes(x = GeneRatio, y = Description)) +
    geom_point(aes(size = Count, color = p.adjust)) +
    scale_size(range = c(2, 8)) + 
    labs(title = gsub("_go_results_bp.csv", "", file), x = "GeneRatio", y = "Description") +
    theme(axis.text.y = element_text(size = base_size))
  

  go_output_file <- file.path(output_path, gsub("_bp.csv", "_top10.png", file))
  

  ggsave(go_output_file, go_plot, width = 16, height = 8, dpi = 300)
  
}




objects_to_save <- c("subdir_path_t_test","working_dir")

rm(list = setdiff(ls(), objects_to_save))    


# Load required libraries
library(clusterProfiler)
library(org.Hs.eg.db)
library(biomaRt)
library(dplyr)
library(limma)
library(ggplot2)


#####KEGG-GSEA analysis#######
# Set the directory containing the CSV files
csv_dir <- subdir_path_t_test

subdir <- "GSEA_KEGG"
subdir_path_GSEA_KEGG <- file.path(working_dir, subdir)
if (!dir.exists( subdir_path_GSEA_KEGG)) {
  dir.create( subdir_path_GSEA_KEGG)
}
output_folder <-  subdir_path_GSEA_KEGG


# Get a list of all CSV files in the directory
csv_files <- list.files(path = csv_dir, pattern = "*.csv", full.names = TRUE)

# Loop through each CSV file
for (csv_file in csv_files) {
  # Read in the CSV file
  result <- read.csv(csv_file, header = TRUE)
  result <- result[!grepl("_2", result$Protein_ID), ]
  
  
  fold_change_pvalue <- data.frame(matrix(ncol = 3, nrow = nrow(result)))
  fold_change_pvalue[,1] <- result$Protein_ID
  fold_change_pvalue[,2] <- result$difference
  fold_change_pvalue[,3] <- result$p_value
  column_name_fold <- c("Protein_ID","logFC","P.Value")
  names(fold_change_pvalue) <- column_name_fold
  
  

  genelist <- fold_change_pvalue %>%
    filter(!is.na(Protein_ID), !is.na(logFC)) %>% 
    mutate(ranking_metric = -log10(P.Value)*sign(logFC)) %>% 
    group_by(Protein_ID) %>% 
    summarise(ranking_metric = mean(ranking_metric, na.rm = TRUE)) %>% 
    arrange(-ranking_metric) %>% 
    tibble::deframe() 
  
  
  res <- gseKEGG(
    genelist,    
    organism = "hsa",    
    keyType = "uniprot",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    minGSSize = 15, 
    maxGSSize = 500, 
    eps = 0, 
    nPermSimple = 10000,
    seed = FALSE)
  
  df <- as.data.frame(res@result)
  
  # Set output file paths
  kegg_output_file <- file.path(output_folder, paste0(gsub("\\.csv$", "", basename(csv_file)), "_GSEA_KEGG.csv"))
  
  write.csv(df, file = kegg_output_file, row.names = FALSE)
  
  
  plot <- dotplot(res, showCategory = 10, split = ".sign") + facet_grid(. ~ .sign)
  
  kegg_output_png <- file.path(output_folder, paste0(gsub("\\.csv$", "", basename(csv_file)), "_GSEA_KEGG.png"))

  ggsave(kegg_output_png, plot, width = 16, height = 12, dpi = 300)
  
  
  
  
}

