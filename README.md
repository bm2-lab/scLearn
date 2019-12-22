# **scLearn: Learning for single cell assignment**

* **scLearn** is a learning-based model for single cell assignment. We applied scLearn to four main tasks of single cell assignment with different perspectives and annotation levels cross 35 datasets and proved that scLearn outperforms all existing methods.
* **scLearn** intuitively carrys out a search like scmap-cluster by measuring the similarity between query cells and each reference cluster centroid but with metric and similarity threshold learned from trained reference datasets. It consists of three steps: data preprocessing, model learning and assignment:
  * **Data preprocessing**: Besides the routine normalization and quality control, scLearn removes rare cell types whose cell number are less than 10 from reference datasets. After the two steps above, scLearn performs feature selection with M3Drop which is based on dropout rate and has been proved suitable for single cell assignment.
  * **Model learning**: scLearn builds a learning-based model. In this model, scLearn first randomly selects a part of samples from each class of the labeled reference with the information (similar or dissimilar) of these selected samples. By bootstrapping ten times, scLearn finally gets a stable classifier with learned metric specific to the trained rederence datasets.
  * **Assignment**: With trained learning-based classifier, scLearn performs assignment for the query cells with learned metric and learned threshold.
  


* For illustration purpose, we took the dataset **data_example/baron_human.rds**. as an example.
    * **Install**: You can install the **scLearn** package from Github using **devtools** packages with R>=3.6.1.
    ```
    library(devtools)
    install.github("bm2-lab/scLearn")
    ```
    * **Data preprocessing**:
    ```
    # loading the reference dataset
    data<-readRDS('example_data/baron-human.rds')
    rawcounts=assays(data)[[1]]
    ann<-as.character(data$cell_type1)
    names(ann)<-colnames(data)
    
    # cell quality control and rare cell type filtered and feature selection
    data_qc<-Cell_qc(rawcounts,ann,species="Hs")
    data_type_filtered<-Cell_type_filter(data_qc$expression_profile,data_qc$sample_information,min_cell_number = 10)
    high_varGene_names <- Feature_selection_M3Drop(data_type_filtered$expression_profile)
    ```
    
    * **model learning**:
    ```
    # training the classifier
    scLearn_classifier_result<-scLearn_classifier(high_varGene_names,data_type_filtered$expression_profile,data_type_filtered$sample_information)
    ```
    
    * **Assignment**:
    ```
    # loading the quary cell and performing cell quality control
    data2<-readRDS('example_data/xin-human.rds')
    rawcounts2=assays(data2)[[1]]
    ann2<-as.character(data2$cell_type1)
    names(ann2)<-colnames(data2)
    ann2<-ann2[ann2 %in% c("alpha","beta","delta","gamma")]
    rawcounts2<-rawcounts2[,names(ann2)]
    data_qc_query<-Cell_qc(rawcounts2,ann2,species="Hs")
    data_type_filtered2<-data_qc_query
    
    # assignment with trained classifier
    scLearn_predict_result<-scLearn_predict(scLearn_classifier_result,data_qc_query$expression_profile)
    
    ```
