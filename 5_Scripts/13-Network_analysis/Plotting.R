## Singh AK et al: Proteins with amino acid repeats constitute a rapidly evolvable and human-specific essentialome
##This file contains all scripts for plotting.

#1. Plotting all ggplot2 plots
##Note: Add hash to the lines which are not required when plotting any specific plot.
df=read.table(pipe("pbpaste"),header=TRUE,sep="\t") #Input file
library(ggplot2)
library(ggridges)
ggplot(df,aes(df$X_axis,df$Y_axis,fill=df$Legend))+
  geom_bar(stat="identity",mapping=NULL,position="dodge")+ #For bar plot
  geom_histogram()+ #For histogram
  geom_density(na.rm=TRUE)+ #For density plot
  geom_line(stat="identity",mapping=NULL,group=0)+ #For line plot
  geom_smooth(method="lm")+ #For smoothened line plot
  geom_point(aes(size=df$Size))+ #For bubble plot
  geom_boxplot(outlier.shape=NA,na.rm=TRUE)+ #For boxplot
  stat_density_ridges(jittered_points=TRUE,position=position_points_jitter(width=0,height=0),point_shape='|',point_size=3,point_alpha=1,alpha=0.7)+ #For ridge plot
  geom_tile()+ #For heatmap
  scale_fill_continuous(high="#000066",low="#A0A0E0")+ #For continuous colour gradient
  geom_violin(trim=TRUE)+ #For violin plot
  coord_cartesian(ylim=c(0,100),xlim=c(0,500))+ #Setting x and y-axis limits
  scale_x_discrete(limits=c("Sample 1","Sample 2","Sample 3"))+ #For assigning order of x-axis labels
  labs(title="Title of plot",x="X-axis label",y="Y-axis label")+ #To set plot title and axes labels
  coord_flip()+ #To flip x and y-axis
  theme_classic()

##end----------------------end



#2. Plotting a circos plot
##2.1. Input and processing the file
df=read.table(pipe("pbpaste"),header=TRUE,sep="\t") #Input file with processes 1 and 2 and number of intersecting genes in three distinct columns
df$Process1=factor(df$Process1,levels=unique(df$Process1),ordered=TRUE) #Freezing the sequence of process 1
df$Process2=factor(df$Process2,levels=unique(df$Process2),ordered=TRUE) #Freezing the sequence of process 2

library(reshape2)
df2=dcast(df,df$Process1 ~ df$Process2)
df3=df2[,c(2:5)] #Removing row names for final input
processes=c("Signaling","Development and differentiation","Cell-cell communication","Transcription") #Assign names to all the processes

##2.2. Plotting
dimnames(df3)=list(have=processes,prefer=processes)
groupColors=c("#FFDD89", "#957244","#800000","#A0A0A0","#CCCC00","#009999","#009900","#FFFF00","#CC00CC","#994C00") #Assigning colours to each process

library(circlize)
chordDiagram(df3,transparency=0.5)

##end----------------------end
