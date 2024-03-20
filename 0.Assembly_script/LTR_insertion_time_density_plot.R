###Input file source from EDTA/S1/S1_hifi.asm.bp.p_ctg.fa.mod.EDTA.raw/LTR/S1_hifi.asm.bp.p_ctg.fa.mod.pass.list & EDTA/S2/S2_hifi.asm.bp.p_ctg.fa.mod.EDTA.raw/LTR/S2_hifi.asm.bp.p_ctg.fa.mod.pass.list
###due to unsetted parametar -u, so using time_convert.py modified the insertion time
library(ggplot2)
data<-read.table("S1_hifi.asm.bp.p_ctg.fa.mod.pass.list2",header = F)
den<-density(data$V12/1000000) ###calculate insertion time density, transfor the unit to MYA, Million Year Age
denx=den$x ##x axis
deny=den$y ##y axis
ggplot()+geom_density(mapping=aes(x=data$V12/1000000),fill="dodgerblue2")+ 
  geom_vline(xintercept = denx[which.max(deny)])+ # Add a vertical line at the x-intercept corresponding to the maximum density
  geom_segment(aes(xend=denx[which.max(deny)],yend=0.715,x=1,y=0.75),arrow=arrow(length=unit(0.2,"cm")))+ # Add an arrow segment pointing from (1, 0.75) to the x-intercept of the maximum density
  annotate("text",x=1.6,y=0.75,label=paste(round(denx[which.max(deny)],4),"MYA"),size=4)+labs(x="LTR insertion time (Million Year Ago)",y="Density")+ # Add text annotation at (1.6, 0.75) with the x-intercept value (in MYA) corresponding to the maximum density
  theme_classic()




