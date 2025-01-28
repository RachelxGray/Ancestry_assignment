ldepth<-read.table("missing10.ldepth.mean", sep="\t", header = TRUE)
mean(ldepth$MEAN_DEPTH)
sd(ldepth$MEAN_DEPTH)