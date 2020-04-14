ge.na.ratio <- function(x){
  sum(is.na(x))/dim(x)[1]/dim(x)[2]
}

ge.split <- function(data,split,which=1){
  sapply(data,function(v){strsplit(v,split)[[1]][which]})
}

# ge.setwd() <- function(){
#   setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# }

ge.readtable <- function(data,sep = "\t",header = T){
  read.table(data,sep = sep,header = header,stringsAsFactors = F)
}

ge.writetable <- function(data,filename ,sep = "\t",col.names = T,row.names = T,quote = F){
  write.table(data,filename,sep=sep,col.names = col.names,row.names = row.names,quote = quote)
}

ge.plot.density <- function(data){
  plot(density(na.omit(unlist(data))),main="density default")
}

ge.remove.techrep <- function(data,pattern="_repB",method="mean"){
  repB <- names(data)[grepl(pattern, names(data))]
  for (i in repB) {
    repa <- str_split(i,pattern)[[1]][1]
    df1 <- data[,grepl(repa, names(data))]
    data <- data[,!grepl(repa, names(data))]
    new_mean <- apply(df1, 1, function(x){ifelse(sum(is.na(x))==2,NA, mean(as.numeric(x),na.rm=T))} )
    data <- cbind(data,new_mean)
    names(data)[ncol(data)] <- repa
  }
  return(data)
}

ge.plot.techrep.correlation <- function(cor1,cor2,name="pearson_correlation"){
  pdf(paste0(name,".pdf"))
  r <- cor(cor1, cor2, use = "pairwise.complete.obs")   
  smoothScatter(cor1, cor2, nrpoints = 100,cex = 2,
                colramp = colorRampPalette(c(blues9,"orange", "red")),
                main = name, xlab = "repA", ylab = "repB")
  abline(lm(cor1 ~ cor2), col="red", lwd=2, lty=2)
  text(min(cor1,na.rm = T)*1.2,max(cor2,na.rm = T)*0.9,labels =paste0( "r =", as.character(round(r,4))),cex = 1.2)
  dev.off()
}

ge.plot.pool.correlation <- function(data,name="bio_cor",method="circle"){
  library(corrplot)
  df_cor <- data.frame(data)
  pdf(paste0(name,".pdf"))
  mycor=cor(df_cor, use = "pairwise.complete.obs")
  corrplot(mycor, method=method,type = "upper",tl.col = "black",tl.srt = 45, tl.cex = 0.5)
  dev.off()
}




ge.plot.boxplot <- function(data,x,y,type,filename,title="boxplot"){
  a <- ggplot(data=data, aes(x =x, y =y ,color=type,group=type)) +
    geom_jitter(alpha = 0.3,size=3) +
    geom_boxplot(alpha = .5,size=1)+
    labs(x="sample",y="value",fill= "type")+
    ggtitle(title)+
    theme_bw() + 
    theme(panel.border = element_blank())+
    theme(axis.line = element_line(size=1, colour = "black")) +
    theme(panel.grid =element_blank())+  
    theme(axis.text = element_text(size = 15,colour = "black"),text = element_text(size = 15,colour = "black"))+
    theme(axis.text.x = element_text( hjust = 1,angle = 45))
  ggsave(paste0(filename, ".pdf"),plot=a,width=8,height=8)
}

# p <- ggboxplot(df1, x="dose", y="len", color = "dose", 
#                palette = c("#00AFBB", "#E7B800", "#FC4E07"), 
#                add = "jitter", shape="dose")
# my_comparisons <- list(c("0.5", "1"), c("1", "2"), c("0.5", "2"))
# p+stat_compare_means(comparisons = my_comparisons)+#不同组间的比较
#   stat_compare_means(label.y = 50)


ge.plot.valcano <- function(data, title,fd=1,pvalue=0.05){
  df.t <- data
  cut.fd <- fd
  pdf(paste0(title, "_volcano.pdf"))
  plot(df.t$fd, -log10(df.t$P_value_adjust), col="#00000033", pch=19,
       xlab=paste("log2 (fold change)"),
       ylab="-log10 (P_value_adjust)",
       main= title)
  
  up <- subset(df.t, df.t$P_value_adjust < pvalue & df.t$fd > cut.fd)
  down <- subset(df.t, df.t$P_value_adjust< pvalue & df.t$fd< as.numeric(cut.fd*(-1)))
  write.csv(up,file = paste0(title, "_up.csv"))
  write.csv(down,file = paste0(title, "_down.csv"))
  points(up$fd, -log10(up$P_value_adjust), col=1, bg = brewer.pal(9, "YlOrRd")[6], pch=21, cex=1.5)
  points(down$fd, -log10(down$P_value_adjust), col = 1, bg = brewer.pal(11,"RdBu")[9], pch = 21,cex=1.5)
  abline(h=-log10(pvalue),v=c(-1*fd,fd),lty=2,lwd=1)
  
  dev.off()
}

ge.plot.pca <- function(data,type,title=""){
  df10 <- data
  df10[is.na(df10)] <- 0
  names <-type
  df10 <- t(apply(df10, 1, scale))
  colnames(df10) <- names
  df.pr <- prcomp(t(df10))
  a<- ggbiplot(df.pr, obs.scale = 1, var.scale = 10, groups =names,alpha = 0,varname.size= 1, ellipse =F, circle = F,var.axes = F)+
    geom_point(aes(colour=names),size = 3,alpha=1)+
    # geom_point(aes(shape=df1$column),size = 3,alpha=1/2)+
    #scale_color_manual(name="type",values=c("#537e35","#e17832","#f5b819","#5992c6","#282f89"))+
    theme(legend.direction = 'horizontal',legend.position = 'top',legend.text = element_text(size = 15,color = "black"), legend.title = element_text(size=15,color="black") ,panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"))+ theme(panel.grid =element_blank())+
    theme(axis.text = element_text(size = 15,color = "black"))+
    theme(plot.subtitle=element_text(size=30, hjust=0, color="black"))+
    theme(axis.title.x=element_text(size=17, hjust=0.5, color="black"))+
    theme(axis.title.y=element_text(size=17, hjust=0.5, color="black"))             #geom_text(aes(label=type,vjust = -0.8, hjust = 0.5,size=0.5),show.legend = FALSE)
  ggsave(paste0(title,"_pca.pdf"),plot =a ,width=12,height=8,device = NULL)
}
ge.plot.pca.label <- function(data,type,label,title=""){
  df10 <- data
  df10[is.na(df10)] <- 0
  names <-type
  df10 <- t(apply(df10, 1, scale))
  colnames(df10) <- names
  df.pr <- prcomp(t(df10))
  a<- ggbiplot(df.pr, obs.scale = 1, var.scale = 10, groups =names,alpha = 0,varname.size= 1, ellipse =F, circle = F,var.axes = F)+
    geom_point(aes(colour=names),size = 3,alpha=1)+
    # geom_point(aes(shape=df1$column),size = 3,alpha=1/2)+
    #scale_color_manual(name="type",values=c("#537e35","#e17832","#f5b819","#5992c6","#282f89"))+
    theme(legend.direction = 'horizontal',legend.position = 'top',legend.text = element_text(size = 15,color = "black"), legend.title = element_text(size=15,color="black") ,panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"))+ theme(panel.grid =element_blank())+
    theme(axis.text = element_text(size = 15,color = "black"))+
    theme(plot.subtitle=element_text(size=30, hjust=0, color="black"))+
    theme(axis.title.x=element_text(size=17, hjust=0.5, color="black"))+
    theme(axis.title.y=element_text(size=17, hjust=0.5, color="black"))+         geom_text(aes(label=label,vjust = -0.8, hjust = 0.5,size=0.2),show.legend = FALSE)
  ggsave(paste0(title,"_pca.pdf"),plot =a ,width=12,height=8,device = NULL)
}

drawPCA<- function(data,type,title="",ptColors=NULL,label=NULL){ 
  M <- t(data)
  M <- apply(M,2,function(v){(v-mean(v,na.rm=T))/sd(v,na.rm=T)})
  clnames <- row.names(data)
  M[is.na(M)] <- 0
  m1 <- prcomp(M);
  Y  <- scale(M, m1$center, m1$scale) %*% m1$rotation 
  Y  <- Y[,c(1,2)]
  
  Y <- data.frame(Y,type);
  colnames(Y) <- c("PC1","PC2","label")
  eigs <- m1$sdev^2
  percentages <- eigs[1:2] / sum(eigs)
  p <- ggplot(Y, aes(x=PC1, y=PC2, colour=label)) + geom_point(size=4)
  p <- p + theme(  panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.border = element_blank(),
                   axis.line.x = element_line(color="black", size = 0.25),
                   axis.line.y = element_line(color="black", size = 0.25),
                   plot.title   = element_text(size=16),
                   panel.background = element_blank())
  
  strLabx <- sprintf("PC1(%4.2f%%)",percentages[1]*100)
  p <- p +  labs(x =strLabx,y = sprintf("PC2(%4.2f%%)",percentages[2]*100),
                 title =sprintf("PCA:%d features",length(clnames)))
  if(!is.null(ptColors)){
    p <- p +   scale_colour_manual(values=ptColors)
  }
  if(!is.null(label)){
    p <- p +   geom_text(aes(label=type,vjust = -0.8, hjust = 0.5,size=0.5),show.legend = FALSE)
  }
  
  ggsave(paste0(title,"_pca.pdf"),plot =p ,device = NULL)
}
ge.plot.tsne <- function(data,type,title=""){
  df10 <- data
  df10[is.na(df10)] <- 0
  names <-type
  df10 <- t(apply(df10, 1, scale))
  colnames(df10) <- names
  
  color <- factor(names) #,levels = c("red","#74A9CF")
  pdf(paste0(title,"_TNSE.pdf"))
  df11.tsne <- Rtsne(t(df10), dims = 2, perplexity = (ncol(data)-1)/3-1, verbose = T , check_duplicates = FALSE)
  plot(df11.tsne$Y,col=color, main = "tsne", pch = 20,cex=2,cex.axis=2,cex.lab=2)
  
  plot(df11.tsne$Y, type = "n", main = "tsne", pch = 20)
  text(df11.tsne$Y, labels = type, col= "DimGrey",cex = 0.8)
  dev.off()
}

ge.plot.umap<- function(data,type,title=""){
  df10 <- data
  df10[is.na(df10)] <- 0
  names <-type
  df10 <- t(apply(df10, 1, scale))
  colnames(df10) <- names
  
  color <- factor(names) #,levels = c("red","#74A9CF")
  pdf(paste0(title,"_UMAP.pdf"))
  df.umap <- umap(t(df10))
  plot(df.umap$layout,col = color, main = "umap", pch = 20,cex=2,cex.axis=2,cex.lab=2)
  
  plot(df.umap$layout, type = "n", main = "umap", pch = 20)
  text(df.umap$layout, labels = type, col= "DimGrey",cex = 0.8)
  dev.off()
}

ge.plot.bar <- function(data,sample,value,group=group,title="",xlab="sample",ylab="value"){
  a <- ggplot(data,aes(x=sample,y=value,group=group))+ 
    geom_bar(position = "dodge",stat = "identity",width =0.8,alpha=0.8,aes(fill=group))+
    ggtitle(paste0(title,"_barplot"))+
    xlab(xlab)+
    ylab(ylab)+
    theme(legend.text = element_text(size = 15,color = "black"),legend.position = 'top',
          legend.title = element_text(size=15,color="black") ,
          panel.grid.major =element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"))+
    theme(panel.grid =element_blank())+
    theme(axis.text = element_text(size = 10,color = "black"))+
    theme(axis.text.x = element_text( hjust = 1,angle = 45))+
    theme(plot.subtitle=element_text(size=30, hjust=0, color="black"))+
    theme(axis.title.x=element_text(size=17, hjust=0.5, color="black"))+
    theme(axis.title.y=element_text(size=17, hjust=0.5, color="black"))  #geom_text(aes(label=value,vjust = -0.8, hjust = 0.5),show.legend = FALSE)
  ggsave(paste0(title,"_barplot.png"),plot=a,width=10,height=8)
}


ge.plot.line <- function(data,sample,value,group=group,title="",xlab="sample",ylab="value"){
  a <- ggplot(data,aes(x=sample,y=value,group=group,color=group))+ 
    geom_line()+
    geom_point()+
    ggtitle(paste0(title,"_lineplot"))+
    xlab(xlab)+
    ylab(ylab)+
    theme(legend.text = element_text(size = 15,color = "black"),legend.position = 'top',
          legend.title = element_text(size=15,color="black") ,
          panel.grid.major =element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"))+
    theme(panel.grid =element_blank())+
    theme(axis.text = element_text(size = 10,color = "black"))+
    theme(axis.text.x = element_text( hjust = 1,angle = 45))+
    theme(plot.subtitle=element_text(size=30, hjust=0, color="black"))+
    theme(axis.title.x=element_text(size=17, hjust=0.5, color="black"))+
    theme(axis.title.y=element_text(size=17, hjust=0.5, color="black"))+ geom_text(aes(label=group,vjust = -0.8, hjust = 0.5),show.legend = FALSE)
  ggsave(paste0(title,"_lineplot.png"),plot=a,width=10,height=8)
}


# ge.plot.vioplot <- function(sample1,sample2,title="",xlab="sample",ylab="value"){
# pdf(paste0(title, "_violin.pdf"))
# vioplot(sample1 ,sample2  ,
#          areaEqual=FALSE, 
#         # rectCol= color, col= color,
#         lineCol=c("black", "black"),
#         border=c("black","black"),
#         names=c("DIANN","DIANN_quantile"),
#         main="biological replicates", xlab=xlab, ylab=ylab,plotCentre = "point")
# dev.off()
# }


ge.mfuzz.cselection <- function(data,range=seq(5,50,5),repeats = 5){
  df3a<-as.matrix(data)
  df3Ex<- ExpressionSet(assayData = df3a)
  if(interactive()){
    df3F <- filter.NA(df3Ex)
    df3F <- fill.NA(df3F)
    df3F <- standardise(df3F)
  }
  
  df3F <- filter.NA(df3F)
  m<-mestimate(df3F)
  cselection(df3F,m=m,crange = range,repeats = repeats,visu = T)
  return(df3F)
}

ge.mfuzz.getresult <- function(data, pic,filename){
  cl <- mfuzz(data,c=pic,m=1.25)
  dir.create(path=filename,recursive = TRUE)
  pdf(paste0(filename,".pdf"))
  mfuzz.plot2(data, cl=cl,mfrow=c(4,4),centre=TRUE,x11=F,centre.lwd=0.2)#min.mem=0.99
  dev.off()
  
  for(i in 1:pic){
    potname<-names(cl$cluster[unname(cl$cluster)==i])
    write.csv(cl[[4]][potname,i],paste0(filename,"/mfuzz_",i,".csv"))
  }
}



