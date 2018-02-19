#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  args[2] = "out.pdf"
}

## program...
library(ggplot2)
library(gridExtra)
#print(args[2])
df = read.table(args[1], header=FALSE,sep="\t")
colnames(df)<-c("Chr","Position","Coverage")
maxcoverage<-max(df$Coverage)
#print(maxcoverage)
c1<-ggplot(df, aes(x=Position,y=Coverage)) + geom_line() + labs(x = "Position", y = "Depth of Coverage") + theme_bw() + scale_y_continuous(limits = c(1, maxcoverage)) + facet_grid(Chr~.)
ggsave(file=args[2], width=8, height=4)
  

#gene<-levels(df$Chr)
#p <- list()
#for (i in gene) {
# # print(i)
#  flist<-c(i)
#  subset<-df[df$Chr %in% flist,]
#  lab<-i
# # print (lab)
#  p[[i]]<-ggplot(subset, aes(x=Position,y=Coverage)) + geom_line() + labs(x = "Position", y = "Depth of Coverage") + scale_y_continuous(limits = c(1, maxcoverage)) + labs(title = lab) + theme_bw() 
#}
  
#pdf(file=args[2],width=8.3, height=11.7)
#do.call(grid.arrange, c(p, ncol=2))
#dev.off()

