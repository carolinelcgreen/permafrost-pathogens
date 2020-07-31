# Function to get the probability into a whole matrix not half, 
# here is Spearman you can change it to Kendall or Pearson
cor.prob.all <- function (X, dfr = nrow(X) - 2) {
  R <- cor(X, use="pairwise.complete.obs",method="spearman")
  r2 <- R^2
  Fstat <- r2 * dfr/(1 - r2)
  R<- 1 - pf(Fstat, 1, dfr)
  R[row(R) == col(R)] <- NA
  R
}
# Change matrices to dataframes
nbar<- as.data.frame(cor(nba[2:ncol(nba)]),method="spearman") # to a dataframe for r^2
nbap<- as.data.frame(cor.prob.all(nba[2:ncol(nba)])) # to a dataframe for p values
# Reset rownames
nbar <- data.frame(row=rownames(nbar),nbar) # create a column called "row" 
rownames(nbar) <- NULL
nbap <- data.frame(row=rownames(nbap),nbap) # create a column called "row" 
rownames(nbap) <- NULL
# Melt
nbar.m <- melt(nbar)
nbap.m <- melt(nbap)
# Classify (you can classify differently for nbar and for nbap also)         
nbar.m$value2<-cut(nbar.m$value,breaks=c(-1,-0.75,-0.5,-0.25,0,0.25,0.5,0.75,1),
                   include.lowest=TRUE, 
                   label=c("(-0.75,-1)","(-0.5,-0.75)","(-0.25,-0.5)",
                           "(0,-0.25)","(0,0.25)","(0.25,0.5)","(0.5,0.75)",
                           "(0.75,1)")) # the label for the legend
nbap.m$value2<-cut(nbap.m$value,
                   breaks=c(-Inf, 0.001, 0.01, 0.05),
                   label=c("***", "** ", "*  ")) 
nbar.m<-cbind.data.frame(nbar.m,nbap.m$value,nbap.m$value2) # adding the p value and its cut to the first dataset of R coefficients
names(nbar.m)[5]<-paste("valuep") # change the column names of the dataframe 
names(nbar.m)[6]<-paste("signif.")
nbar.m$row <- factor(nbar.m$row,
                     levels=rev(unique(as.character(nbar.m$variable)))) # reorder the variable factor
# Plotting the matrix correlation heatmap
# Set options for a blank panel
po.nopanel <-list(opts(panel.background=theme_blank(),
                       panel.grid.minor=theme_blank(),
                       panel.grid.major=theme_blank()))
pa<-ggplot(nbar.m, aes(row, variable)) +
  geom_tile(aes(fill=value2),colour="white") +
  scale_fill_brewer(palette = "RdYlGn",name="Correlation")+ # RColorBrewer package
  opts(axis.text.x=theme_text(angle=-90))+
  po.nopanel
pa # check the first plot
# Adding the significance level stars using geom_text 
pp<- pa +
  geom_text(aes(label=signif.),size=2,na.rm=TRUE) # you can play with the size
# Workaround for the alpha aesthetics if it is good to represent significance level, the same workaround can be applied for size aesthetics in ggplot2 as well. Applying the alpha aesthetics to show significance is a little bit problematic, because we want the alpha to be low while the p value is high, and vice verse which can't be done without a workaround
nbar.m$signif.<-rescale(as.numeric(nbar.m$signif.),to=c(0.1,0.9)) # I tried to use to=c(0.1,0.9) argument as you might expect, but to avoid problems with the next step of reciprocal values when dividing over one, this is needed for the alpha aesthetics as a workaround
nbar.m$signif.<-as.factor(0.09/nbar.m$signif.) # the alpha now behaves as wanted  except for the NAs values stil show as if with three stars level, how to fix that?
# Adding the alpha aesthetics in geom_point in a shape of squares (you can improve here)
pp<- pa +
  geom_point(data=nbar.m,aes(alpha=signif.),
             shape=22,size=5,colour="darkgreen",
             na.rm=TRUE,legend=FALSE) # you can remove this step, the result of this step is seen in one of the layers in the above green heatmap, the shape used is 22 which is again a square but the size you can play with it accordingly  

