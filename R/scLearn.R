### cell quality control
Cell_qc<-function(expression_profile,sample_information_cellType=NULL,sample_information_timePoint=NULL,species="Hs",gene_low=500,gene_high=10000,mito_high=0.1,umi_low=1500,umi_high=Inf,logNormalize=TRUE,plot=FALSE,plot_path="./quality_control.pdf"){
  require(stringr)
  if(species=="Hs"){
    mito.genes <- grep("^MT-", rownames(expression_profile), value = FALSE)
  }else if(species=="Mm"){
    mito.genes <- grep("^mt-", rownames(expression_profile), value = FALSE)
  }else{
    stop("species should be 'Mm' or 'Hs'")
  }
  mito_percent<-function(x,mito.genes){
    return(sum(x[mito.genes])/sum(x))
  }    
  percent.mito<-apply(expression_profile,2,mito_percent,mito.genes=mito.genes)
  nUMI<-apply(expression_profile,2,sum)
  nGene<-apply(expression_profile,2,function(x){length(x[x>0])})
  expression_profile<-apply(expression_profile,2,function(x){x/(sum(x)/10000)})
  if(logNormalize==TRUE){
    expression_profile<-log(expression_profile+1)    
  }
  if(species=="Hs"){
    expression_profile<-expression_profile[!str_detect(row.names(expression_profile),"^MT-"),]
  }else if(species=="Mm"){
    expression_profile<-expression_profile[!str_detect(row.names(expression_profile),"^mt-"),]
  }
  filter_retain<-rep("filter",ncol(expression_profile))
  for(i in 1:length(filter_retain)){
    if(nGene[i]>gene_low & nGene[i]<gene_high & percent.mito[i]<mito_high & nUMI[i]>umi_low & nUMI[i]<umi_high){
      filter_retain[i]<-"retain"
    }
  }
  if(plot){
    pdf(file=plot_path)
    par(mfrow=c(1,3))
    hist(nGene,breaks = 20,freq=FALSE,xlab = "Gene numbers",ylab = "Density",main = "Gene numbers distribution")
    lines(nGene,col="red",lwd=1)
    hist(nUMI,breaks = 20,freq=FALSE,xlab = "UMI numbers",ylab = "Density",main = "UMI numbers distribution")
    lines(density(nUMI),col="red",lwd=1)
    hist(percent.mito,breaks = 20,freq=FALSE,xlab = "Percent of mito",ylab = "Density",main = "Percent of mito distribution")
    lines(density(percent.mito),col="red",lwd=1)
    dev.off()
  }
  SQ_filter<-as.matrix(expression_profile[,which(filter_retain=="retain")])
  SQ_data_qc<-SQ_filter#have adopted total_expr=10000 and log
  if(is.null(sample_information_cellType)){
    if(is.null(sample_information_timePoint)){
      return(list("expression_profile"=SQ_data_qc))
    }
  }
 if(is.null(sample_information_timePoint)){
    sample_information_qc<-sample_information_cellType[colnames(SQ_data_qc)]
    return(list("expression_profile"=SQ_data_qc,"sample_information_cellType"=sample_information_qc))
  }
  else{
    sample_information_qc<-gsub("_","-",sample_information_cellType[colnames(SQ_data_qc)])
    sample_information2_qc<-gsub("_","-",sample_information_timePoint[colnames(SQ_data_qc)])
    return(list("expression_profile"=SQ_data_qc,"sample_information_cellType"=sample_information_qc,"sample_information_timePoint"=sample_information2_qc))
  }
}

### filter out cell type without the minimal cell number
Cell_type_filter<-function(expression_profile,sample_information_cellType,sample_information_timePoint=NULL,min_cell_number=10){
  if(is.null(sample_information_timePoint)){
    sample_information<-sample_information_cellType
  }else{
    sample_information<-paste(sample_information_cellType,sample_information_timePoint,sep="_")
  }
  names(sample_information)<-names(sample_information_cellType)
  cell_number_type<-table(sample_information)
  cell_type_filtered<-names(cell_number_type[cell_number_type<min_cell_number])
  sample_information_filtered<-sample_information[!(sample_information %in% cell_type_filtered)]
  expression_profile_filtered<-expression_profile[,names(sample_information_filtered)]
  sample_information_cellType<-sample_information_cellType[names(sample_information_filtered)]
  if(!is.null(sample_information_timePoint)){
    sample_information_timePoint<-sample_information_timePoint[names(sample_information_filtered)]
    return(list("expression_profile"=expression_profile_filtered,"sample_information_cellType"=sample_information_cellType,"sample_information_timePoint"=sample_information_timePoint))
  }else{
    return(list("expression_profile"=expression_profile_filtered,"sample_information_cellType"=sample_information_cellType))
  }
}

### feature selection
Feature_selection_M3Drop<-function(expression_profile,log_normalized=TRUE,threshold=0.05){
  library(M3Drop)
  if(log_normalized){
    high_varGenes <- M3DropFeatureSelection(exp(expression_profile)-1,mt_method="fdr", mt_threshold=threshold)
    high_varGene_names<-row.names(high_varGenes)
    return(high_varGene_names)
  }else{
    high_varGenes <- M3DropFeatureSelection(expression_profile,mt_method="fdr", mt_threshold=threshold)
    high_varGene_names<-row.names(high_varGenes)
    return(high_varGene_names)
  }
}

### train the reference 
scLearn_model_learning<-function(high_varGene_names,expression_profile,sample_information_cellType,sample_information_timePoint=NULL,bootstrap_times=10,cutoff=0.01,dim_para=0.999){
  qiu_ji<-function(x,y){
        return(x %*% y)
    }
  Feature_cluster<-function(expression_profile,sample_information){
    sample_information<-sort(sample_information)
    expression_profile<-expression_profile[,names(sample_information)]
    num_each_class<-table(as.character(sample_information))
    feature_matrix<-matrix(0,length(num_each_class),nrow(expression_profile))
    row.names(feature_matrix)<-names(num_each_class)
    a=1
    for(i in 1:(length(num_each_class))){
      #print(i)
      expression_profile_choose<-expression_profile[,a:(a+num_each_class[i]-1)]
      a=a+num_each_class[i]
      feature_matrix[i,]<-apply(expression_profile_choose,1,mean)
    }
    return(feature_matrix)
  }
  Threshold_similarity<-function(expression_profile,sample_information,cutoff=0.01){
    sample_information<-sort(sample_information)
    expression_profile<-expression_profile[,names(sample_information)]
    num_each_class<-table(as.character(sample_information))
    simi_mean_list<-list()
    km<-function(vec,k){
      km<-mean(sort(vec,decreasing = TRUE)[2:(k+1)])
      return(km)
    }
    a=1
    thre<-c()
    for(i in 1:(length(num_each_class))){
      #print(i)
      expression_profile_choose<-expression_profile[,a:(a+num_each_class[i]-1)]
      a=a+num_each_class[i]
      feature_vecter<-apply(expression_profile_choose,1,mean)
      simi_mean_list[[i]]<-apply(expression_profile_choose,2,cor,y=feature_vecter)
      thre[i]<-sort(simi_mean_list[[i]],decreasing = F)[ceiling(length(simi_mean_list[[i]])*cutoff)+1]
    }
    names(thre)<-names(num_each_class)
    names(simi_mean_list)<-names(num_each_class)
    return(list("simi_mean_list"=simi_mean_list,"threshold"=thre))
  }
  runDCA<-function(high_varGenes,expression_profile,sample_information,strength=0.1,seed=1){
    require(dml)
    sample_information<-sort(sample_information)
    expression_profile<-expression_profile[high_varGenes,names(sample_information)]
    num_eachclass<-table(sample_information)
    #print("building chunks ...")
    chks<-list()
    a=1
    k=1
    for(i in 1:(length(num_eachclass))){
      s<-ceiling(num_eachclass[i]*strength)+1
      set.seed(seed)
      chks[[k]]=sample(a:(a+num_eachclass[i]-1),s)
      k=k+1
      set.seed(seed)
      chks[[k]]=sample(setdiff(a:(a+num_eachclass[i]-1),chks[[k-1]]),s)
      k=k+1
      a=a+num_eachclass[i]
    }
    chunks = rep(-1, sum(num_eachclass))
    for (i in 1:(2*(length(num_eachclass)))) {
      for (j in chks[[i]]) {
        chunks[j] = i
      }
    }
    #print("building negtive links ...")
    neglinks = matrix(rep(1,4*(length(num_eachclass))*(length(num_eachclass))),ncol = 2*(length(num_eachclass)), byrow = TRUE)
    for(i in seq(1,ncol(neglinks),2)){
      neglinks[i,i]=0
      neglinks[i,i+1]=0
      neglinks[i+1,i]=0
      neglinks[i+1,i+1]=0
    }
    #print("performing dca ...")
    dca_result = dca(data = t(expression_profile), chunks = chunks, neglinks = neglinks)
    newData<-dca_result$newData
    colnames(newData)<-paste("dc",1:ncol(newData),sep="")
    trans_matrix<-dca_result$DCA
    ex<-t(newData)
    #print("Done!")
    return(list("expression_profile_trans"=ex,"expression_profile_origin"=expression_profile,"trans_matrix"=trans_matrix,"sample_information"=sample_information))
  }
  if(is.null(sample_information_timePoint)){
    threshold_cluster_trans<-list()
    feature_matrix_trans<-list()
    trans_matrix<-list()
    if(bootstrap_times<2){
      warning("The bootstrap_times should be at least larger than 1.")
    }
    for(r in 1:bootstrap_times){
      print(paste("Bootstrapying",r))
      trans_result<-runDCA(high_varGene_names,expression_profile,sample_information_cellType,strength = 0.1,seed=r)
      trans_matrix[[r]]<-trans_result$trans_matrix
      thre_result_trans_cluster<-Threshold_similarity(trans_result$expression_profile_trans,trans_result$sample_information,cutoff=cutoff)
      threshold_cluster_trans[[r]]<-thre_result_trans_cluster$threshold
      feature_matrix_trans[[r]]<-Feature_cluster(trans_result$expression_profile_trans,trans_result$sample_information)
    }
    return(list("high_varGene_names"=high_varGene_names,"simi_threshold_learned"=threshold_cluster_trans,"feature_matrix_learned"=feature_matrix_trans,"trans_matrix_learned"=trans_matrix))
  }else{
    label_len<-length(table(sample_information_cellType))+length(table(sample_information_timePoint))
    label_matrix<-matrix(0,label_len,ncol(expression_profile))
    
    row.names(label_matrix)<-c(names(table(sample_information_cellType)),names(table(sample_information_timePoint)))
    
    colnames(label_matrix)<-colnames(expression_profile)
    
    sample_information_all<-c(sample_information_cellType,sample_information_timePoint)
    label_all<-c(names(table(sample_information_cellType)),names(table(sample_information_timePoint)))
    for (i in 1:nrow(label_matrix)) {
        label_matrix[i,names(sample_information_all[sample_information_all==label_all[i]])] <- 1
    }
    
    L_kernel<-matrix(0,ncol(label_matrix),ncol(label_matrix))
    colnames(L_kernel)<-colnames(expression_profile)
    row.names(L_kernel)<-colnames(expression_profile)
    L_kernel<-t(label_matrix) %*% label_matrix
    getProperDim<-function(lambda,dim_para=0.999){
        thr<-dim_para
        sum_lambda<-sum(lambda)
        lambda_num<-length(lambda)
        tmp_lambda<-0
        for(lind in 1:lambda_num){
          tmp_lambda<-tmp_lambda+lambda[lind]
          if(tmp_lambda>=thr*sum(lambda)){
            proper_dim<-lind
            return(proper_dim)  
          }
        }
    }
    mddm_linear<-function(X,L,dim_para=0.999){
      D<-nrow(X)
      N<-ncol(X)
      nor_mean<-function(vec){
        vec<-vec-mean(vec)
        return(vec)
      }
      tmpL<-apply(L,2,nor_mean)
      HLH<-apply(tmpL,1,nor_mean)
      S=X %*% HLH %*% t(X)
      B=diag(D)
      B_tr=B
      eig_result<-eigen(B_tr %*% S)
      tmp_lambda<-diag(eig_result$values)
      tmp_P<-eig_result$vectors
      tmp_P<-Re(tmp_P)
      tmp_lambda <- Re(tmp_lambda[row(tmp_lambda)==col(tmp_lambda)])
      lambda<-sort(tmp_lambda,decreasing=T)
      order<-order(tmp_lambda,decreasing = T)
      P<-tmp_P[,order]
      proper_dim<- getProperDim(lambda, dim_para);
      P<- P[,1:proper_dim]
      lambda<- lambda[1:proper_dim]
      return(list("Projection_matrix"=P,"eigenvalues"=lambda))        
    }
    
    mddm_result<-mddm_linear(expression_profile[high_varGene_names,],L_kernel,dim_para = dim_para)
    trans_result<-list()
    trans_result$trans_matrix<-t(mddm_result$Projection_matrix)
    trans_result$expression_profile_orign<-expression_profile[high_varGene_names,]
    trans_result$expression_profile_trans<-t(mddm_result$Projection_matrix) %*% expression_profile[high_varGene_names,]
    sample_information_cb<-paste(sample_information_cellType,sample_information_timePoint,sep="_")
    names(sample_information_cb)<-names(sample_information_cellType)
    trans_result$sample_information_combine<-sample_information_cb
    
    thre_result_trans_cluster<-Threshold_similarity(trans_result$expression_profile_trans,trans_result$sample_information_combine,cutoff=cutoff)
    
    threshold_cluster_trans<-thre_result_trans_cluster$threshold
    feature_matrix_trans<-Feature_cluster(trans_result$expression_profile_trans,trans_result$sample_information_combine)
    
    return(list("high_varGene_names"=high_varGene_names,"simi_threshold_learned"=threshold_cluster_trans,"feature_matrix_learned"=feature_matrix_trans,"trans_matrix_learned"=trans_result$trans_matrix))
  }
}

### draw graph with tsne or umap
DrawCluster <- function(data, label = NULL, point_size = 1,method=c("tsne","umap"),draw_cluster_text=TRUE,calculated=TRUE,pca=TRUE,perplexity=100,plot=TRUE,seed=1){
  require(ggplot2)
  require(dplyr)
  set.seed((seed))
  method=method[1]
  if(calculated){
    if(method=="tsne"){
      require(Rtsne)
      tsneresult2 <- Rtsne(t(data), perplexity = perplexity, pca = pca)
      X <- as.data.frame(tsneresult2$Y)
    }else if(method=="umap"){
      require(umap)
      umapresult1<-umap(t(data))
      X<-as.data.frame(umapresult1$layout)
    }else{
      print("method must be tsne or umap.")
      break
    }
  }else{
    X=data
  }
  if(length(label) == 0){
    label <- array(1, dim(X)[1])
    labelname = c(1)
  }
  labelname<-names(table(label))
  p <- ggplot(X, aes(x=X[,1], y=X[,2]))
  cell_group=factor(label)
  if(method=="tsne"){
    p <- p + geom_point(aes(color=cell_group), size = point_size)+ xlab("tSNE1") + ylab("tSNE2")
  }else{
    p <- p + geom_point(aes(color=cell_group), size = point_size)+ xlab("umap1") + ylab("umap2")
  }
  if(draw_cluster_text){
    Label_cal<-X
    Label_cal$cluster<-label
    cluster_x_y<-Label_cal %>% group_by(cluster) %>% summarise(x_median=median(V1),y_median=median(V2))
    p<-p+annotate("text",x=cluster_x_y$x_median,y=cluster_x_y$y_median,label=cluster_x_y$cluster)
  }
  mytheme <- theme_bw() +
    theme(plot.title=element_text(size=rel(1.5),hjust=0.5),
          axis.title=element_text(size=rel(1)),
          axis.text=element_text(size=rel(1)),
          panel.grid.major=element_line(color="white"),
          panel.grid.minor=element_line(color="white"),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 15)
    )
  p<-p + mytheme + guides(colour = guide_legend(override.aes = list(size = 4)))
  if(plot){
    print(p)
  }
  #print(p + mytheme)
  return(list("p"=p,"x"=X,"cell_group"=cell_group))
}

### the column is sample names, draw the similarity graph
correlation<-function(matrix,method=c("pearson","spearman","cosin","euclidean"),cpu_num=8){
  cosdist <- function(x1,x2){
    n1 <- sqrt(sum(x1^2))
    n2 <- sqrt(sum(x2^2))
    d <- as.numeric(x1 %*% x2) / n1 / n2
    d
  }
  euclidean<-function(x1,x2){
    return(sqrt(t(x1-x2) %*% (x1-x2)))
  }
  method<-method[1]
  if(!(method %in% c("pearson","spearman","cosin","euclidean"))){
    print("method must be 'pearson','spearman','cosin' or 'euclidean'.")
    break
  }
  require(parallel)
  cpu_num_set <- makeCluster(getOption("cluster.cores", cpu_num))
  sample_num<-ncol(matrix)
  cor_matrix<-matrix(rep(1,sample_num^2),sample_num)
  colnames(cor_matrix)<-colnames(matrix)
  row.names(cor_matrix)<-colnames(matrix)
  simi_cor<-function(vec,matrix,method){
    if(method=="cosin"){
      cor_matrix<-apply(matrix,2,cosdist,x2=vec)
    }else if(method=="euclidean"){
      cor_matrix<-apply(matrix,2,euclidean,x2=vec)
    }else{
      cor_matrix<-apply(matrix,2,cor,y=vec,method=method)
    }
    return(cor_matrix)
  }
  cor_matrix<-parApply(cpu_num_set,matrix,2,simi_cor,matrix=matrix,method=method)
  stopCluster(cpu_num_set)
  return(cor_matrix)
}

### predict cell type
scLearn_cell_assignment<-function(scLearn_model_learning_result,expression_profile_query,vote_rate=0.6,diff=0.05){
  Vote_class<-function(vec,vote_rate=0.6){
    vec_len<-length(vec)
    num<-length(table(vec)[table(vec)==max(table(vec))])
    if(num>1){
      return("unassigned")
    }else if(max(table(vec))>vec_len*vote_rate){
      return(names(table(vec)[table(vec)==max(table(vec))]))
    }else{
      return("unassigned")
    }
  }
  Get_query_hvg<-function(expression_profile,high_varGene_names){
    missing_num<-length(high_varGene_names)-length(intersect(row.names(expression_profile),high_varGene_names))
    missing_features<-setdiff(high_varGene_names,intersect(row.names(expression_profile),high_varGene_names))
    missing_rate<-missing_num/length(high_varGene_names)
    print(paste("The number of missing features in the query data is ",missing_num,seq=""))
    print(paste("The rate of missing features in the query data is ",missing_rate,seq=""))
    if(missing_num>0){
      missing_data<-matrix(0,missing_num,ncol(expression_profile))
      row.names(missing_data)<-missing_features
      expression_profile<-rbind(expression_profile,missing_data)
    }
    expression_profile_hvg<-expression_profile[high_varGene_names,]
    return(expression_profile_hvg)
  }
  Assignment_result<-function(expression_profile_query_hvg,feature_matrix,threshold,diff=0.05){
    options(warn=-1)
    feature_matrix<-t(feature_matrix)
    result<-data.frame(cluster_lab=1:ncol(expression_profile_query_hvg),unassigned_case=rep(NA,ncol(expression_profile_query_hvg)),cluster_cor=rep(0,ncol(expression_profile_query_hvg)))
    row.names(result)<-colnames(expression_profile_query_hvg)
    for(i in 1:ncol(expression_profile_query_hvg)){
      cor_result<-rep(0,ncol(feature_matrix))
      names(cor_result)<-colnames(feature_matrix)
      for(j in 1:ncol(feature_matrix)){
        if(is.na(cor(expression_profile_query_hvg[,i],feature_matrix[,j]))){
          cor_result[j]<-0
        }else{
          cor_result[j]<-cor(expression_profile_query_hvg[,i],feature_matrix[,j])
          }
      }
      cor_compare<-cor_result-threshold 
      if(length(cor_compare[cor_compare>0])==0){
        result[i,1]<-"unassigned"
        result[i,2]<-"Dissimilar to all reference cell types"
        result[i,3]<-NA
      }
      if(length(cor_compare[cor_compare>0])==1){
        result[i,1]<-names(cor_compare[cor_compare>0])
        result[i,3]<-cor_result[result[i,1]]
      }
      if(length(cor_compare[cor_compare>0])>=2){
        if(sort(cor_result[names(cor_compare[cor_compare>0])],decreasing=T)[1]-sort(cor_result[names(cor_compare[cor_compare>0])],decreasing=T)[2]>=diff){
          result[i,1]<-names(cor_result[names(cor_compare[cor_compare>0])])[which(cor_result[names(cor_compare[cor_compare>0])]==max(cor_result[names(cor_compare[cor_compare>0])]))][1]
          result[i,3]<-cor_result[result[i,1]]            
        }else{
          result[i,1]<-"unassigned"
          result[i,2]<-"Similar to multiple existing labels"
          result[i,3]<-NA
        }
      }
    }
    return(result)
  }
  expression_profile_query_hvg<-Get_query_hvg(expression_profile_query,scLearn_model_learning_result$high_varGene_names)
  if(is.list(scLearn_model_learning_result$trans_matrix_learned)){
    predict_result<-matrix(0,ncol(expression_profile_query),length(scLearn_model_learning_result$trans_matrix_learned))
    for(r in 1:length(scLearn_model_learning_result$trans_matrix_learned)){
      expression_profile_query_hvg_ml<-scLearn_model_learning_result$trans_matrix_learned[[r]] %*% expression_profile_query_hvg
      assignment_result<-Assignment_result(expression_profile_query_hvg_ml,scLearn_model_learning_result$feature_matrix_learned[[r]],threshold = scLearn_model_learning_result$simi_threshold_learned[[r]],diff=diff)
      predict_result[,r]<-assignment_result[,1]
    }
    predict_result_final<-apply(predict_result,1,Vote_class,vote_rate=vote_rate)
    predict_result_final<-as.data.frame(predict_result_final)
    predict_result_final$sample<-colnames(expression_profile_query)
    predict_result_final$predict_result_final<-as.character(predict_result_final$predict_result_final)
    colnames(predict_result_final)<-c("Predict_cell_type","Query_cell_id")
    predict_result_final<-predict_result_final[,c("Query_cell_id","Predict_cell_type")]
    return(predict_result_final)
  }else{
    expression_profile_query_hvg_ml <- scLearn_model_learning_result$trans_matrix_learned %*% expression_profile_query_hvg
    assignment_result <- Assignment_result(expression_profile_query_hvg_ml, scLearn_model_learning_result$feature_matrix_learned, threshold = scLearn_model_learning_result$simi_threshold_learned, diff = diff)
    assignment_result$Query_cell_id <- row.names(assignment_result)
    assignment_result$Predict_cell_type <- as.character(assignment_result$cluster_lab)
    assignment_result<- assignment_result[, c("Query_cell_id","Predict_cell_type")]
    return(assignment_result)
  }
}

### sankey graph
sankey_plot<-function(predict_result,sample_information_reference,plot=FALSE){
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
  if(plot){
    require(scmap)
    plot(getSankey(predict_result[,1], predict_result[,2], plot_width = 300, plot_height = 300))
  }
  return(predict_result)
}

