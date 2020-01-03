# **scLearn: Learning for single cell assignment**

* **scLearn** is a learning-based framework for single cell assignment. We applied scLearn to four main tasks of single cell assignment with different perspectives and annotation levels across 30 datasets and proved that scLearn outperforms all existing methods.
* **scLearn** intuitively carrys out a search like scmap-cluster by measuring the similarity between query cells and each reference cluster centroid but with metric and similarity threshold learned from trained reference datasets. It consists of three steps: data preprocessing, model learning and assignment:
  * **Data preprocessing**: Besides the routine normalization and quality control, scLearn removes rare cell types whose cell number are less than 10 from reference datasets. After the two steps above, scLearn performs feature selection with M3Drop which is based on dropout rate and has been proved suitable for single cell assignment.
  * **Model learning**: scLearn builds a learning-based model. In this model, scLearn first randomly selects a part of samples from each class of the labeled reference with the information (similar or dissimilar) of these selected samples. By bootstrapping ten times, scLearn finally gets a stable classifier with learned metric specific to the trained rederence datasets.
  * **Assignment**: With trained learning-based classifier, scLearn performs assignment for the query cells with learned metric and learned threshold.
  


* For illustration purpose, we took the dataset **data_example/baron_human.rds**. as an example.
    * **Install**: You can install the **scLearn** package from Github using **devtools** packages with R>=3.6.1.
    ```
    library(devtools)
    install_github("bm2-lab/scLearn")
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

| Trained model names | Description | No. of cell types | corresponding dataset | PMID |
| :------: | :------: | :------: | :------: | :------: |
| pancreas_mouse_baron.rds | Mouse pancreas |  | Baron |  |
| pancreas_human_baron.rds | Human pancreas |  | Baron |  |
| pancreas_human_muraro.rds | Human pancreas |  | Muraro |  |
| pancreas_human_segerstolpe.rds | Human pancreas |  | Segerstolpe |  |
| pancreas_human_xin.rds | Human pancreas |  | Xin |  |
| embryo_development_mouse_deng.rds | Mouse embryo development |  | Deng |  |
| cerebral_cortex_human_pollen.rds | Human cerebral cortex |  | Pollen |  |
| colorectal_tumor_human_li.rds | Human colorectal tumors |  | Li |  |
| brain_mouse_usoskin.rds | Mouse brain |  | Usoskin |  |
| cortex_mouse_tasic.rds | Mouse cortex |  | Tasic |  |
| embryo_stem_cell_mouse_klein.rds | Mouse embryo stem cells |  | Klein |  |
| brain_mouse_zeisel.rds | Mouse brain |  | Zeisel |  |
| retina_mouse_shekhar_shallow_annotation.rds | Mouse retina |  | Shekhar |  |
| retina_mouse_shekhar_deep_annotation.rds | Mouse retina |  | Shekhar |  |
| retina_mouse_macosko.rds | Mouse retina |  | Macosko |  |
| lung_cancer_cell_lines_human_cellbench10X.rds | Mixture of five human lung cancer cell lines |  | CellBench10X |  |
| lung_cancer_cell_lines_human_cellbenchCelSeq.rds | Mixture of five human lung cancer cell lines |  | CellBenchCelSeq |  |
| whole_mus_musculus_mouse_TM.rds | Whole Mus musculus |  | TM |  |
| primary_visual_cortex_mouse_AMB_shallow_annotation_3.rds | Primary mouse visual cortex |  | AMB |  |
| primary_visual_cortex_mouse_AMB_shallow_annotation_16.rds | Primary mouse visual cortex |  | AMB |  |
| primary_visual_cortex_mouse_AMB_shallow_annotation_92.rds | Primary mouse visual cortex |  | AMB |  |
| PBMC_human_zheng_sorted.rds | FACS-sorted PBMC |  | Zheng sorted |  |
| PBMC_human_zheng_68K.rds | PBMC |  | Zheng 68K |  |
| primary_visual_cortex_mouse_VISP.rds | Mouse primary visual cortex |  | VISP |  |


