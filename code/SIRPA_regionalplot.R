# regional plot as part of the Figure 5a

library(data.table)
library(dplyr)
library(ggplot2)


SIRPA <- fread('/path/to/SIRPA_region.txt')

SIRPA$FINEMAP <- as.factor(SIRPA$FINEMAP)


p1 <- ggplot(SIRPA, aes(x=BEG, y=-log10(as.double(PVAL)), color=FINEMAP)) + geom_point( size=1.5) + scale_color_manual(values = c("darkgrey","red1")) + theme_classic() + theme(legend.position = "none") + scale_y_continuous(limits = c(0,120))  + geom_hline(yintercept = -log10(5e-8),linetype="dashed", col="darkgrey")

p1  + scale_y_continuous(breaks = c(0,10,25,50,75,100,110)) + scale_x_continuous(breaks = c(1500000,1544167,1561029,1567964,1600000,1700000,1800000,1875154,1900000,1920543,2000000),labels = c("1,500,000","","","","1,600,000","1,700,000","1,800,000","","1,900,000","","2,000,000")) + ylab("-log 10 P-value") + xlab("")

