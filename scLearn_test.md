# **scLearn: Learning for single cell assignment**


```R

```





data<-readRDS('~/my_project/assignment_project/data/human/baron-human.rds')
data
rawcounts=assays(data)[[1]]
ann<-as.character(data$cell_type1)
names(ann)<-colnames(data)
table(ann)



```R
data_qc<-Cell_qc(rawcounts,ann,species="Hs")
data_type_filtered<-Cell_type_filter(data_qc$expression_profile,data_qc$sample_information,min_cell_number = 10)
high_varGene_names <- Feature_selection_M3Drop(data_type_filtered$expression_profile)
length(high_varGene_names)
```

    Loading required package: stringr
    Loading required package: numDeriv
    Warning message in bg__calc_variables(expr_mat):
    “Warning: Removing 2642 undetected genes.”


1889



![png](output_3_2.png)



```R
metric_learned<-runDCA(high_varGene_names,data_type_filtered$expression_profile,data_type_filtered$sample_information)
```

    Loading required package: dml
    Loading required package: MASS


    [1] "building chunks ..."
    [1] "building negtive links ..."
    [1] "performing dca ..."
    [1] "Done!"



```R
dim(metric_learned$expression_profile_trans)
```


<ol class=list-inline>
	<li>25</li>
	<li>8529</li>
</ol>




```R
#pdf("~/metric_learning_result_new/metric_learning_effect_tsne.pdf")
library(gridExtra)
ex_ml_tsne<-DrawCluster(metric_learned$expression_profile_trans,label=metric_learned$sample_information,pca = FALSE,method="tsne",calculated = TRUE,plot=FALSE)
ex_tsne<-DrawCluster(metric_learned$expression_profile_origin,label=metric_learned$sample_information,pca = FALSE,method="tsne",calculated = TRUE,plot=FALSE)
grid.arrange(ex_tsne$p,ex_ml_tsne$p, nrow=2, ncol=1)
#dev.off()
```

    Loading required package: Rtsne



![png](output_6_1.png)



```R
data2<-readRDS('~/my_project/assignment_project/data/human/xin-human.rds')
data2
rawcounts2=assays(data2)[[1]]
ann2<-as.character(data2$cell_type1)
names(ann2)<-colnames(data2)
table(ann2)
```


    class: SingleCellExperiment 
    dim: 39851 1600 
    metadata(0):
    assays(2): normcounts logcounts
    rownames(39851): A1BG A2M ... LOC102724004 LOC102724238
    rowData names(1): feature_symbol
    colnames(1600): Sample_1 Sample_2 ... Sample_1599 Sample_1600
    colData names(6): donor.id condition ... gender cell_type1
    reducedDimNames(0):
    spikeNames(1): ERCC
    altExpNames(0):



    ann2
                 alpha alpha.contaminated               beta  beta.contaminated 
                   886                 60                472                 31 
                 delta delta.contaminated              gamma gamma.contaminated 
                    49                  9                 85                  8 



```R
ann2<-ann2[ann2 %in% c("alpha","beta","delta","gamma")]
rawcounts2<-rawcounts2[,names(ann2)]
```


```R
data_qc_query<-Cell_qc(rawcounts2,ann2,species="Hs")
data_type_filtered2<-data_qc_query
table(data_type_filtered$sample_information)
table(data_type_filtered2$sample_information)
```


    
                acinar activated_stellate              alpha               beta 
                   958                283               2317               2520 
                 delta             ductal        endothelial            epsilon 
                   599               1072                246                 18 
                 gamma         macrophage               mast quiescent_stellate 
                   255                 55                 22                171 
               schwann 
                    13 



    
    alpha  beta delta gamma 
      885   472    49    85 



```R
scLearn_classifier_result<-scLearn_classifier(high_varGene_names,data_type_filtered$expression_profile,data_type_filtered$sample_information)
```

    [1] "Calculating 1"
    [1] "Calculating 2"
    [1] "Calculating 3"
    [1] "Calculating 4"
    [1] "Calculating 5"
    [1] "Calculating 6"
    [1] "Calculating 7"
    [1] "Calculating 8"
    [1] "Calculating 9"
    [1] "Calculating 10"



```R

scLearn_predict_result<-scLearn_predict(scLearn_classifier_result,data_qc_query$expression_profile)
```

    [1] "The number of missing features in the query data is  39 "
    [1] "The rate of missing features in the query data is  0.0206458443620963 "



```R
head(scLearn_predict_result)
head(data_qc_query$sample_information,30)
```


<table>
<caption>A data.frame: 6 × 2</caption>
<thead>
	<tr><th scope=col>predict_cell_type</th><th scope=col>sample_name</th></tr>
	<tr><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><td>beta</td><td>Sample_1</td></tr>
	<tr><td>beta</td><td>Sample_2</td></tr>
	<tr><td>beta</td><td>Sample_3</td></tr>
	<tr><td>beta</td><td>Sample_4</td></tr>
	<tr><td>beta</td><td>Sample_5</td></tr>
	<tr><td>beta</td><td>Sample_6</td></tr>
</tbody>
</table>




<dl class=dl-horizontal>
	<dt>Sample_1</dt>
		<dd>'beta'</dd>
	<dt>Sample_2</dt>
		<dd>'beta'</dd>
	<dt>Sample_3</dt>
		<dd>'beta'</dd>
	<dt>Sample_4</dt>
		<dd>'beta'</dd>
	<dt>Sample_5</dt>
		<dd>'beta'</dd>
	<dt>Sample_6</dt>
		<dd>'beta'</dd>
	<dt>Sample_7</dt>
		<dd>'beta'</dd>
	<dt>Sample_8</dt>
		<dd>'beta'</dd>
	<dt>Sample_9</dt>
		<dd>'beta'</dd>
	<dt>Sample_10</dt>
		<dd>'beta'</dd>
	<dt>Sample_11</dt>
		<dd>'beta'</dd>
	<dt>Sample_12</dt>
		<dd>'beta'</dd>
	<dt>Sample_13</dt>
		<dd>'beta'</dd>
	<dt>Sample_14</dt>
		<dd>'beta'</dd>
	<dt>Sample_15</dt>
		<dd>'beta'</dd>
	<dt>Sample_16</dt>
		<dd>'beta'</dd>
	<dt>Sample_17</dt>
		<dd>'beta'</dd>
	<dt>Sample_18</dt>
		<dd>'beta'</dd>
	<dt>Sample_19</dt>
		<dd>'beta'</dd>
	<dt>Sample_20</dt>
		<dd>'beta'</dd>
	<dt>Sample_21</dt>
		<dd>'beta'</dd>
	<dt>Sample_22</dt>
		<dd>'beta'</dd>
	<dt>Sample_23</dt>
		<dd>'beta'</dd>
	<dt>Sample_24</dt>
		<dd>'beta'</dd>
	<dt>Sample_25</dt>
		<dd>'beta'</dd>
	<dt>Sample_26</dt>
		<dd>'beta'</dd>
	<dt>Sample_27</dt>
		<dd>'beta'</dd>
	<dt>Sample_28</dt>
		<dd>'beta'</dd>
	<dt>Sample_29</dt>
		<dd>'beta'</dd>
	<dt>Sample_30</dt>
		<dd>'beta'</dd>
</dl>




```R
estimate_result<-Estimation_result(scLearn_predict_result,data_qc_query$sample_information,data_type_filtered$sample_information)
```


```R
estimate_result
```


<dl>
	<dt>$F1_meature</dt>
		<dd>0.99394754539341</dd>
	<dt>$accuracy</dt>
		<dd>0.991281019450034</dd>
	<dt>$precision</dt>
		<dd>0.996628455832771</dd>
	<dt>$recall</dt>
		<dd>0.991281019450034</dd>
	<dt>$true_negative_ratio</dt>
		<dd>0</dd>
	<dt>$false_unassigned_ratio</dt>
		<dd>0.00536552649228706</dd>
	<dt>$true_unassigned_ratio</dt>
		<dd>NaN</dd>
	<dt>$cluster_original_projection</dt>
		<dd><table>
<caption>A matrix: 1491 × 3 of type chr</caption>
<thead>
	<tr><th scope=col>orign_label</th><th scope=col>projection_label</th><th scope=col>sample_name</th></tr>
</thead>
<tbody>
	<tr><td>beta</td><td>beta</td><td>Sample_1 </td></tr>
	<tr><td>beta</td><td>beta</td><td>Sample_2 </td></tr>
	<tr><td>beta</td><td>beta</td><td>Sample_3 </td></tr>
	<tr><td>beta</td><td>beta</td><td>Sample_4 </td></tr>
	<tr><td>beta</td><td>beta</td><td>Sample_5 </td></tr>
	<tr><td>beta</td><td>beta</td><td>Sample_6 </td></tr>
	<tr><td>beta</td><td>beta</td><td>Sample_7 </td></tr>
	<tr><td>beta</td><td>beta</td><td>Sample_8 </td></tr>
	<tr><td>beta</td><td>beta</td><td>Sample_9 </td></tr>
	<tr><td>beta</td><td>beta</td><td>Sample_10</td></tr>
	<tr><td>beta</td><td>beta</td><td>Sample_11</td></tr>
	<tr><td>beta</td><td>beta</td><td>Sample_12</td></tr>
	<tr><td>beta</td><td>beta</td><td>Sample_13</td></tr>
	<tr><td>beta</td><td>beta</td><td>Sample_14</td></tr>
	<tr><td>beta</td><td>beta</td><td>Sample_15</td></tr>
	<tr><td>beta</td><td>beta</td><td>Sample_16</td></tr>
	<tr><td>beta</td><td>beta</td><td>Sample_17</td></tr>
	<tr><td>beta</td><td>beta</td><td>Sample_18</td></tr>
	<tr><td>beta</td><td>beta</td><td>Sample_19</td></tr>
	<tr><td>beta</td><td>beta</td><td>Sample_20</td></tr>
	<tr><td>beta</td><td>beta</td><td>Sample_21</td></tr>
	<tr><td>beta</td><td>beta</td><td>Sample_22</td></tr>
	<tr><td>beta</td><td>beta</td><td>Sample_23</td></tr>
	<tr><td>beta</td><td>beta</td><td>Sample_24</td></tr>
	<tr><td>beta</td><td>beta</td><td>Sample_25</td></tr>
	<tr><td>beta</td><td>beta</td><td>Sample_26</td></tr>
	<tr><td>beta</td><td>beta</td><td>Sample_27</td></tr>
	<tr><td>beta</td><td>beta</td><td>Sample_28</td></tr>
	<tr><td>beta</td><td>beta</td><td>Sample_29</td></tr>
	<tr><td>beta</td><td>beta</td><td>Sample_30</td></tr>
	<tr><td>⋮</td><td>⋮</td><td>⋮</td></tr>
	<tr><td>gamma</td><td>gamma</td><td>Sample_1569</td></tr>
	<tr><td>gamma</td><td>gamma</td><td>Sample_1570</td></tr>
	<tr><td>gamma</td><td>gamma</td><td>Sample_1571</td></tr>
	<tr><td>gamma</td><td>gamma</td><td>Sample_1572</td></tr>
	<tr><td>gamma</td><td>gamma</td><td>Sample_1573</td></tr>
	<tr><td>gamma</td><td>gamma</td><td>Sample_1574</td></tr>
	<tr><td>gamma</td><td>gamma</td><td>Sample_1575</td></tr>
	<tr><td>gamma</td><td>gamma</td><td>Sample_1576</td></tr>
	<tr><td>gamma</td><td>gamma</td><td>Sample_1577</td></tr>
	<tr><td>gamma</td><td>gamma</td><td>Sample_1578</td></tr>
	<tr><td>gamma</td><td>gamma</td><td>Sample_1579</td></tr>
	<tr><td>gamma</td><td>gamma</td><td>Sample_1580</td></tr>
	<tr><td>gamma</td><td>gamma</td><td>Sample_1581</td></tr>
	<tr><td>gamma</td><td>gamma</td><td>Sample_1582</td></tr>
	<tr><td>gamma</td><td>gamma</td><td>Sample_1583</td></tr>
	<tr><td>gamma</td><td>gamma</td><td>Sample_1584</td></tr>
	<tr><td>gamma</td><td>gamma</td><td>Sample_1585</td></tr>
	<tr><td>gamma</td><td>gamma</td><td>Sample_1586</td></tr>
	<tr><td>gamma</td><td>gamma</td><td>Sample_1587</td></tr>
	<tr><td>gamma</td><td>gamma</td><td>Sample_1588</td></tr>
	<tr><td>gamma</td><td>gamma</td><td>Sample_1589</td></tr>
	<tr><td>gamma</td><td>gamma</td><td>Sample_1590</td></tr>
	<tr><td>gamma</td><td>gamma</td><td>Sample_1591</td></tr>
	<tr><td>gamma</td><td>gamma</td><td>Sample_1592</td></tr>
	<tr><td>gamma</td><td>gamma</td><td>Sample_1593</td></tr>
	<tr><td>gamma</td><td>gamma</td><td>Sample_1594</td></tr>
	<tr><td>gamma</td><td>gamma</td><td>Sample_1595</td></tr>
	<tr><td>gamma</td><td>gamma</td><td>Sample_1597</td></tr>
	<tr><td>gamma</td><td>gamma</td><td>Sample_1598</td></tr>
	<tr><td>gamma</td><td>gamma</td><td>Sample_1600</td></tr>
</tbody>
</table>
</dd>
</dl>




```R
library(scmap)
```

    Creating a generic function for ‘toJSON’ from package ‘jsonlite’ in package ‘googleVis’
    
    Attaching package: ‘scmap’
    
    The following object is masked _by_ ‘.GlobalEnv’:
    
        ann
    



```R
labels<-estimate_result$cluster_original_projection

```


```R
head(labels)
sankey_prepare<-function(predict_result,sample_information_reference){
  require(dplyr)
  cell_type_reference<-names(table(sample_information_reference))
  other_cell_type<-setdiff(cell_type_reference,intersect(cell_type_reference,unique(predict_result[,2])))
  if(length(other_cell_type)>0){
    other_cell_type<-as.data.frame(other_cell_type)
    other_cell_type$orign_label<-"NULL"
    other_cell_type$sample_name<-"for_sankey_plot"
    colnames(other_cell_type)<-c("projection_label","orign_label","sample_name")
    predict_result<-rbind(predict_result,other_cell_type)
  }
  return(predict_result)
}

```


<table>
<caption>A matrix: 6 × 3 of type chr</caption>
<thead>
	<tr><th scope=col>orign_label</th><th scope=col>projection_label</th><th scope=col>sample_name</th></tr>
</thead>
<tbody>
	<tr><td>beta</td><td>beta</td><td>Sample_1</td></tr>
	<tr><td>beta</td><td>beta</td><td>Sample_2</td></tr>
	<tr><td>beta</td><td>beta</td><td>Sample_3</td></tr>
	<tr><td>beta</td><td>beta</td><td>Sample_4</td></tr>
	<tr><td>beta</td><td>beta</td><td>Sample_5</td></tr>
	<tr><td>beta</td><td>beta</td><td>Sample_6</td></tr>
</tbody>
</table>




```R
labels_prepare<-sankey_prepare(labels,data_type_filtered$sample_information)
```


```R
tail(labels_prepare)
```


<table>
<caption>A data.frame: 6 × 3</caption>
<thead>
	<tr><th></th><th scope=col>orign_label</th><th scope=col>projection_label</th><th scope=col>sample_name</th></tr>
	<tr><th></th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>1494</th><td>NULL</td><td>ductal            </td><td>for_sankey_plot</td></tr>
	<tr><th scope=row>1495</th><td>NULL</td><td>endothelial       </td><td>for_sankey_plot</td></tr>
	<tr><th scope=row>1496</th><td>NULL</td><td>epsilon           </td><td>for_sankey_plot</td></tr>
	<tr><th scope=row>1497</th><td>NULL</td><td>macrophage        </td><td>for_sankey_plot</td></tr>
	<tr><th scope=row>1498</th><td>NULL</td><td>quiescent_stellate</td><td>for_sankey_plot</td></tr>
	<tr><th scope=row>1499</th><td>NULL</td><td>schwann           </td><td>for_sankey_plot</td></tr>
</tbody>
</table>




```R
save(labels_prepare,file="~/metric_learning_result_new/labels_example_sankey.RData")
```


```R
plot(getSankey(labels[,1], labels[,2], plot_width = 300, plot_height = 300))

```


```R
a
```


'/tmp/RtmpQz5BoA/SankeyIDe41937a9ea2e.html'



```R

```
