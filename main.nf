#!/usr/bin/env nextflow

process data_processing {
	input:
	path expData
	
	output:
	path "data_scaled.rds"

	script:
	"""
	#!/usr/bin/env Rscript

	# Define version of classifier
	version<-"v26"
	# Import normalized single-cell input data
	data <- readRDS("$expData")
	data<-as.matrix(data)
	# Import zstand parameters for scaling from Github
	param_url<-paste0("https://raw.githubusercontent.com/dmelguizo/classifier/main/6_z_stand_param_train_completeset_",version,".rds")
	zstand_parameters <- readRDS(url(param_url, "rb"))
	# Scale input data (z-standarization)
	data_scaled<-(data - zstand_parameters[["means_train"]]) / zstand_parameters[["sds_train"]]
	# Transpose scaled data
	data_scaled<-t(data_scaled)
	# Export scaled data for prediciton with classifier
	saveRDS(data_scaled, file = "data_scaled.rds")
	"""
}

process prediction_PSC_state{
	input:
	path expData
	
	output:
	path "df_prediction_results.rds"

	script:
	"""
	#!/usr/bin/env Rscript

	# Define version classifier
	version<-"v26"
	# Load libraries required
	library(keras)
	library(tensorflow)
	# Import scaled data
	data_scaled<-readRDS("$expData")
	# Import classifier from GitHub
	h5_url <- "https://raw.githubusercontent.com/dmelguizo/classifier/main/6_best_model_completeset_v26.h5"
	temp_file <- tempfile(fileext = ".h5")
	download.file(h5_url, temp_file, mode = "wb")
	classifier <- load_model_hdf5(temp_file)
	# Predict cell state
	prediction<-predict(classifier, data_scaled)
	# Create data frame with results (labels)
	df<-data.frame(predicted=as.factor(ifelse(max.col(prediction[,1:2])==1,"early_diff","undiff")))
	# Export data frame with prediction results
	saveRDS(data_scaled, file = "data_scaled.rds")
	saveRDS(df, file = "df_prediction_results.rds")
	"""
}

process visualization_results {
	input:
	path expData1
	path expData2

	output:
	path "Results_prediction_PSC_state.pdf"
	
	script:
	"""
        #!/usr/bin/env Rscript

        # Define version classifier
        version<-"v26"
        # Load libraries required
        library("dplyr")
	library("ggpubr")
	library("ggplot2")
	# Modify data frame with results for plotting
	df<-readRDS("$expData2")
	df2 <- as.data.frame(df %>% group_by(predicted) %>% summarise(count = n()) %>% mutate(perc = count/sum(count)))
	plot<-ggplot(data=df2, aes(x=predicted, y=perc*100, fill=predicted)) +
	geom_bar(stat="identity") +
	geom_text(aes(label=round(perc*100,2)), vjust=1.6, color="black", size=5)+
	labs(title="PSC state prediction", x ="Predicted state",y = "Percentage of cells") +
	theme_light()+
	theme(axis.text.x = element_text(color="black",size=14),
		axis.text.y = element_text(color="black",size=12),
	        axis.title.x = element_text(face="bold", color="black",size=15),
	        axis.title.y = element_text(face="bold", color="black",size=15),
	        plot.title = element_text(color="black", size=20,hjust = 0.5)) +
	guides(fill="none")
	# Report
	data_scaled<-readRDS("$expData1")
	sample<-"ID"
	date<-as.character(Sys.Date())
	num_cells<-nrow(data_scaled)
	num_genes<-ncol(data_scaled)
	perc<-df2["perc"]
	perc_undiff<-round(perc[2,1],2)*100
	perc_earlydiff<-round(perc[1,1],2)*100
	report<-as.data.frame(rbind(sample,date,version,num_cells,num_genes,perc_undiff,perc_earlydiff))
	colnames(report)<-NULL
	rownames(report)<-c("SampleID","Date and Time","Version classifier","Number of cells predicted",
		"Number of genes used for prediction","% PSC Undifferentiated","% PSC Early differentiated")
	report<-ggtexttable(report)
	# Export results
	pdf("Results_prediction_PSC_state.pdf")
		ggarrange(report,plot, ncol = 1, nrow = 2)
	dev.off()
	"""
}

workflow {
	input_data = Channel.fromPath(params.input)
	data_scaled = data_processing(input_data)
	prediction_results = prediction_PSC_state(data_scaled)
	plots = visualization_results(data_scaled,prediction_results)
}
