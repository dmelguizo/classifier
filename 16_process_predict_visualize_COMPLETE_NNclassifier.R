# setwd("C:/Users/dariu/Dropbox/0_Bioinformatics_group_Skovde/8_vinnova/scripts/stemcells/")

# Best model
version<-"v26"

libraries <- c("ggplot2","keras","dplyr","ggpubr")
check.libraries <- is.element(libraries, installed.packages()[, 1])==FALSE
libraries.to.install <- libraries[check.libraries]
if (length(libraries.to.install!=0)) {
  install.packages(libraries.to.install,repos='http://cran.us.r-project.org')
}

success <- sapply(libraries,require, quietly = FALSE,  character.only = TRUE)
if(length(success) != length(libraries)) {stop("A package failed to return a success in require() function.")}

# Notes:

# Requirements for single-cell input data:

# 1. QC -> remove empty/low quality/dead cells
# 2. Samples should have been normalized (using Seurat  SC transform or standard LogNormalization, ....)
  # This normalization only applies to the input data (no leakage since the model has been already trained)

# 1. Process data

#####################################  1.1 Import input data ##################################### 
# Data must be a txt file with the following requirements:
# tab separation (sep='\\t')
# 1st row as column names (header=TRUE)
# 1st column as row names (row.names=1)
data<- readRDS("./output/Round2/6_test_lognorm_completeset.rds") # introduce as input in the workflow
# read.table("$expData", as.is=TRUE, header=TRUE, sep='\\t', row.names=1)
# data<-data[1:3800,1:1000]

###### 1.2 Import parameters used to standarize the model
# From 5_Train_Test_best_model_v26._COMPLETE_SET.R
# Locally
# zstand_parameters<-paste0("./output/6_z_stand_param_train_complete_set_",version,".rds")

# GitHub
param_url<-paste0("https://raw.githubusercontent.com/dmelguizo/classifier/main/6_z_stand_param_train_completeset_",version,".rds")
zstand_parameters <- readRDS(url(param_url, "rb"))

####### 1.3 Scale (z-standarization) input data
  
data_scaled<-as.matrix((data - zstand_parameters[["means_train"]]) / zstand_parameters[["sds_train"]])

####### 1.4 Transpose input data

data_scaled<-t(data_scaled)

################################ 2. Predict using complete NN classifier ################################## 

####### 2.1 Import classifier
# From 5_Train_Test_best_model_v26._COMPLETE_SET.R
classifier<-load_model_hdf5(paste0("./output/6_best_model_completeset_",version,".h5")) # import classifier from GitHub as h5 object

# GitHub
h5_url <- "https://raw.githubusercontent.com/dmelguizo/classifier/main/6_best_model_completeset_v26.h5"
temp_file <- tempfile(fileext = ".h5")
download.file(h5_url, temp_file, mode = "wb")
classifier <- load_model_hdf5(temp_file)

######## 2.2 Predict cell state 

prediction<-predict(classifier, data_scaled) #  early_diff=[1], undiff=[2] 

df<-data.frame(predicted=as.factor(ifelse(max.col(prediction[,1:2])==1,"early_diff","undiff")))

######################################  3. Visualize results #####################################

df2 <- df %>% 
  group_by(predicted) %>%
  summarise(count = n()) %>%
  mutate(perc = count/sum(count))

df2$class<-"prediction"

plot<-ggplot(data=df2, aes(x=predicted, y=perc*100, fill=predicted)) +
  geom_bar(stat="identity") +
  geom_text(aes(label=round(perc*100,2)), vjust=1.6, color="black", size=5)+
  labs(title="PSC state prediction",
       x ="Predicted state", 
       y = "Percentage of cells") +
  theme_light()+
  theme(axis.text.x = element_text(color="black",size=14),
        axis.text.y = element_text(color="black",size=12),
        axis.title.x = element_text(face="bold", color="black",size=15),
        axis.title.y = element_text(face="bold", color="black",size=15),
        plot.title = element_text(color="black", size=20,hjust = 0.5))+
  guides(fill="none")


# Report

sample<-"ID"
date<-as.character(Sys.Date())
num_cells<-ncol(data)
num_genes<-nrow(data)
perc_undiff<-round(df2$perc[2],2)*100
perc_earlydiff<-round(df2$perc[1],2)*100

report<-as.data.frame(rbind(sample,date,version,num_cells,num_genes,perc_undiff,perc_earlydiff))

colnames(report)<-NULL
rownames(report)<-c("SampleID",
                    "Date and Time",
                    "Version classifier",
                    "Number of cells predicted",
                    "Number of genes used for prediction",
                    "% PSC Undifferentiated",
                    "% PSC Early differentiated")

report<-ggtexttable(report)

# Export results
pdf("Results_prediction_PSC_state.pdf")
ggarrange(report,plot
          ncol = 1, nrow = 2)
dev.off()





plot_grid(ncol=1,
          nrow=2,
          rel_widths = c(1, 0.5),
          plot,
          grid.table(report))

# Stacked bar plot
ggplot(df2, aes(x = factor(class), y = perc*100,fill=factor(predicted))) +
  geom_bar(width = 0.2,stat="identity") +
  coord_flip() +
  labs(title="PSC state prediction",
       x ="", 
       y = "Percentage of cells",
       fill = "PSC state") +
  geom_text(aes(label = round(perc,2)), vjust = 1.5, colour = "white")
  

