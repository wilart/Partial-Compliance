library(ggplot2)



pdf("Alc_Coc_histogram4.pdf",width=10,height=10)
ggplot2::ggplot()+geom_histogram(aes(x=dat_for_analysis_4$y),color="gray",fill="black",bins=30)+theme_bw()+ggtitle("log(0.5 + Alcohol Days + Cocaine Days ) Histogram")+ theme(plot.title = element_text(size=24),
                                                                                                                                                                               axis.title = element_text(size=22),
                                                                                                                                                                               axis.text = element_text(size=22))+xlab("Outcome")+ylab("Count")
dev.off()
