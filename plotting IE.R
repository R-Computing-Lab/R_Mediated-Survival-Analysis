data=basic_w
data$IE_2.5<-data$IE_2.5*1000
data$IE_97.5<-data$IE_97.5*1000
data$IE<-data$IE*1000
ggplot(data=data) +  geom_ribbon(aes(x=Time, ymin=IE_2.5, ymax=IE_97.5), fill= 'lightskyblue', alpha=.35) + geom_step(aes(x=Time, y=IE_97.5), directions="hv", linetype=2,alpha=0.5) +geom_step(aes(x=Time, y=IE), direction="hv") + geom_step(aes(x=Time,y=IE_2.5), direction="hv", linetype=2, alpha=0.5) + theme_bw() + theme(axis.text=element_text(size=16), axis.title=element_text(size=16,face="bold"), plot.title=element_text(size=16,face="bold")) +  theme(legend.position="top",  legend.justification=c(1, 0.5)) +   xlab("Years") + ylab("Saved Lives per 1000 People") +  ggtitle("TITLE") + coord_cartesian(xlim=c(45,100))



