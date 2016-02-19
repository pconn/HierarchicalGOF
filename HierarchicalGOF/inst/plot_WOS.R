### create plot from Web of Science search data
require(ggplot2)
Count=c(20,6,2,2,1,1)
Type=c("None","NA","Pp.pval","Pp.gc","Cross.val","Non-Bayes")
Titles=factor(Type,levels=c("None","NA","Pp.pval","Pp.gc","Cross.val","Non-Bayes"))
Plot.data=data.frame(Titles=Titles,Count=Count)

p <- ggplot(Plot.data,aes(x=Titles,y=Count))+geom_bar(stat="identity")+labs(x="Model checking procedure",y="Number of articles")
p <- p + theme(text=element_text(size=16))
p

ggsave(filename="WOSsearch.jpeg")