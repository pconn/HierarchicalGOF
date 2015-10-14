### create plot from Web of Science search data
require(ggplot2)
Count=c(6,20,2,2,1,1)
Type=c("NA","None","Pp.pval","Pp.gc","Cross.val","Non-Bayes")
Plot.data=data.frame(Titles=Titles,COunts=Counts)

p <- ggplot(Plot.data,aes(x=Type,y=Count))+geom_bar(stat="identity")+labs(x="Model checking procedure",y="Number of articles")
p <- p + theme(text=element_text(size=16))
p

ggsave(filename="WOSsearch.jpeg")