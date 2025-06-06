
setwd(dir = "D:/R_project/UCR_project/")
rm(list = ls())

#引用包
library(Matrix)
library(tidyverse)
library(pROC)
library(ggplot2)
library(survival)
library(regplot)
library(ggsci)
library(survminer)
library(timeROC)
library(ggDCA)
library(limma)
library(rms)

inputFile="GSE104954.txt"       #表达矩阵
hub="LASSO.txt"        #核心基因

#读取输入文件
rt=read.table(file = "05-scripts/43-ncRNA_ANN_model/Nomogram/GSE104954.txt", 
              header=T, sep="\t", check.names=F)

rt <- read_tsv(file = "05-scripts/43-ncRNA_ANN_model/Nomogram/GSE104954.txt", col_names = TRUE)

rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=t(data)
sample=read.table("05-scripts/43-ncRNA_ANN_model/Nomogram/sample.txt", 
                  sep="\t", header=F, check.names=F)
colnames(sample)=c("ID","Type")
data=data[sample$ID,]
aSAH1=data[,read.table("05-scripts/43-ncRNA_ANN_model/Nomogram/LASSO.txt", 
                       header=F, sep="\t", check.names=F)[,1]]
aSAH=cbind(sample,aSAH1)
#简单看一下ROC曲线AUC的情况
aflist=roc(Type~HERC6+MTNR1A+MYLIP+PLSCR1, data = aSAH)
g3 <- ggroc(aflist, size = 1.2,alpha=.6,)
g5=g3+ggsci::scale_color_lancet()
print(g5)
#诺曼图，高度比：8比10
############逻辑回归模型
dd <- datadist(aSAH)
options(datadist="dd")
fit <- lrm(formula = Type ~ HERC6+MTNR1A+MYLIP+PLSCR1, data =aSAH)
print(fit)
coef=as.data.frame(fit$coefficients)[-1,,drop=F]
coefout=cbind(ID=rownames(coef),coef)
write.table(coefout,file="coefficients.txt",sep="\t",quote=F,row.names = F)
#绘图
pdf(file="nomogram.pdf", width=9, height=7.5)
plot(nomogram(fit,fun.at = seq(0.05,0.95,0.05)),funlabel = "nomogram model")
dev.off()

plot(regplot(fit,plots=c("density","boxes"), observation=T, title="Prediction Nomogram", clickable=T, points=TRUE,droplines=TRUE))

nomoscore=predict(fit, data=t(aSAH))
aSAH$nomoscore=nomoscore
write.table(aSAH,file="nomoscore.txt",sep="\t",quote=F,row.names = F)
