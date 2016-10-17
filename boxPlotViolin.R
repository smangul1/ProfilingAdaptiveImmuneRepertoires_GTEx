#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)<2) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} 


library("ggplot2")
dataIGH=read.csv(args[1])

png(height=3000, 10000, pointsize=5, file=args[2])

#ggplot(dataIGH, aes(x=tissue, y= d, fill=tissue))+geom_dotplot(binaxis='y', stackdir='center')

#ggplot(dataIGH, aes(x=tissue, y= d, fill=tissue))+geom_boxplot(aes(fill = d))+theme(axis.text=element_text(size=65),axis.title=element_text(size=100,face="bold"))+theme(axis.text.x = element_text(angle=90))

#ggplot(dataIGH, aes(x=tissue, y= d, fill=tissue))+geom_violin(trim = FALSE)  +theme(axis.text=element_text(size=55),axis.title=element_text(size=55,face="bold"))+theme(axis.text.x = element_text(angle=90))

#ggplot(dataIGH, aes(x=tissue, y= d, fill=tissue))+geom_violin(trim = FALSE)  +theme(axis.text=element_text(size=75),axis.title=element_text(size=75,face="bold"))+theme(axis.text.x = element_text(angle=90))+geom_point(stat = "summary", fun.y = "mean", size = I(3), color = I("black")) + geom_point(stat = "summary", fun.y = "mean", size = I(2.2), color = I("orange"))


ggplot(dataIGH, aes(x=tissue, y= d, fill=tissue))+geom_violin(trim = FALSE)  +theme(axis.text=element_text(size=75),axis.title=element_text(size=75,face="bold"))+theme(axis.text.x = element_text(angle=90))+geom_point(stat = "summary", fun.y = "mean", size = I(3), color = I("black")) + geom_point(stat = "summary", fun.y = "mean", size = I(2.2), color = I("orange"))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))

dev.off()

