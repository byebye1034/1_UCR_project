## Singh AK et al: Proteins with amino acid repeats constitute a rapidly evolvable and human-specific essentialome
##This file contains all scripts for statistical analysis.

#1. Computing enrichment using permutation testing: 10000 random sampling
##Note: The following example is to test whether human essential genes are enriched for having homorepeat containing proteins.
##Modify the input lists based on the question.
##1.1 Storing the background and test lists
background=as.vector(unlist(read.table(pipe('pbpaste')))) #All proteins in the dataset
background=unique(background)
length(background)

ess=as.vector(unlist(read.table(pipe('pbpaste')))) #All essential proteins
ess=intersect(ess,background)
length(ess)

hr=as.vector(unlist(read.table(pipe('pbpaste')))) #All homorepeat containing proteins
hr=intersect(hr,background)
length(hr)

hr.ess=intersect(hr,ess) #Overap of essential proteins and homorepeat containing proteins in real dataset
length(hr.ess)

##1.2 Generating random sampling
random=matrix(ncol=10000,nrow=length(ess)) 
for (i in 1:10000){
  sam=sample(background,size=length(ess),replace=FALSE,prob=NULL) #Taking random set (size=length(ess)) from background
  random[,i]=sam
}

##1.3 Checking overlap of random sampling with homorepeat containing proteins
result=vector(length=10000)
for(i in 1:ncol(random)){
  result[i]=length(intersect(random[,i],hr))
}

##1.4 Computing statistical measures
l=length(hr.ess) #Value in real dataset
m=mean(result) #Computing mean of random expectation
s=sd(result) #Computing standard deviation of random expectation
z=(length(hr.ess)-mean(result))/sd(result) #Computing Z-score
pval=length(which(result>=length(hr.ess)))/length(result) #Computing p-value of enrichment
ipval=length(which(result<=length(hr.ess)))/length(result)	#Computing p-value of depletion

##1.5 Output
res=matrix(ncol=6, nrow=2)
res[1,]=c("Value in real dataset","Mean of random expectation", "SD of random expectation", "Z-score","Pval_enrich","Pval_depletion")
res[2,]=c(l,m,s,z,pval,ipval)
write.table(res,pipe('pbcopy'),sep="\t",quote=F,row.names=F,col.names=F)

##1.6 Getting list of random expectation values for plotting as histogram
write.table(result,pipe('pbcopy'),sep="\t",quote=F,row.names=F,col.names=F)

##end----------------------end


#2. Fishers Exact test
df=read.table(pipe("pbpaste"),sep="\t",header=TRUE) #To input contingency matrix
x=fisher.test(df,workspace=2e8,simulate.p.value=TRUE)$p.value #Computing p-value
print(x)

##end----------------------end


#3. Chi-squared test
df=read.table(pipe("pbpaste"),sep="\t",header=TRUE) #To input contingency table with counts
y=chisq.test(df)$p.value #Computing p-value
print(y)

##end----------------------end


#4. Wilcoxon rank sum test
df=read.table(pipe("pbpaste"),sep="\t",header=TRUE) #To input file with Classes and Values as two different columns
df1=subset(df,df$Class=="Class1") #Subset class 1 in a dataframe
df2=subset(df,df$Class=="Class2") #Subset class 2 in a dataframe
z=wilcox.test(df1$Value,df2$Value)$p.value #Computing p-value between two classes
print(z)

##end----------------------end


#5. Multiple hypothesis testing: FDR correction
library(fuzzySim)
df=read.table(pipe("pbpaste"),sep="\t",header=TRUE) #To input table with Test name and P-values as twodifferent columns
a=FDR(pvalues=df$P-values) #Stores all excluded and included p-values after FDR correction

##end----------------------end
