##########
#name:          plot_AIC_BIC.r
#description:   plots AIC, BIC and likelihood ratio (LR) test
#author:        Henriette Miko (henriette.miko@mdc-berlin.de)
#date:          March 28, 2020
##########


library(ggplot2)
library(reshape2)

args = commandArgs(trailingOnly=TRUE)

curdir <- getwd()
setwd(curdir)

AIC = read.table(args[1], header=F)
BIC = read.table(args[2], header=F)
timeless.dir=args[3]
likelihood = read.table(args[4], header=F)
numPar = read.table(args[5], header=F)


#AIC and BIC

df <- as.data.frame(cbind(AIC, BIC))
colnames(df)=c("AIC", "BIC")
df["numbers"] <- 2:30

df_melted <- melt(df, id.vars="numbers")

ggplot(df_melted, aes(x=numbers,y=value, color=variable)) + 
    geom_point() + geom_line() + labs(x="cluster number", y = "criterion") + 
    scale_x_continuous(breaks=2:30, minor_breaks=2:30) +
    ggtitle("Bayesian and Akaike information criterion for model selection") +
    theme(plot.title = element_text(hjust = 0.5)) + 
    theme(legend.title=element_blank())

ggsave(paste0(timeless.dir,"/AIC_BIC.pdf"), height=6, width=8)


#likelihood

df2 <- as.data.frame(likelihood)
colnames(df2)=c("likelihood")
df2["numbers"] <- 2:30

df2_melted <- melt(df2, id.vars="numbers")

ggplot(df2_melted, aes(x=numbers,y=value))  + geom_point() + geom_line() +
    labs(x="cluster number", y = "log likelihood") + 
    scale_x_continuous(breaks=2:30, minor_breaks=2:30) +
    ggtitle("log likelihood") + 
    theme(plot.title = element_text(hjust = 0.5)) + 
    theme(legend.title=element_blank())

ggsave(paste0(timeless.dir,"/likelihood.pdf"),height=6,width=8)


#likelihood ratio test
#2*(LL(L_B)-LL(L_A))

LR = 2*diff(as.vector(t(likelihood))) 

#likelihood is not always decreasing increasing
#set negative numbers to 0
LR[LR<0]=0

#differences of number of parameters
#for example between model11 and model10
#numPar(model11) - numPar(model10)
diffnumPar=diff(as.vector(t(numPar)))


#chi-square distribution with numPar(model11)-numPar(model10), 
#95th quantile
qchisq(0.95,df=diffnumPar)
print(qchisq(0.95,df=diffnumPar))

#threshold stays the same for all cluster numbers
threshold=qchisq(0.95,df=diffnumPar)[1]

print("LR test")
print(which(LR > threshold))
print(which(LR < threshold))

df3 <- as.data.frame(LR)
colnames(df3)=c("LR")
df3["numbers"] <- 2:29 #ratios

df3_melted <- melt(df3, id.vars="numbers")

ggplot(df3_melted, aes(x=numbers,y=value))  + geom_point() + geom_line() +
    labs(x="cluster number and next for ratio", y = "likelihood ratio") + 
    scale_x_continuous(breaks=2:29, minor_breaks=2:29) +
    ggtitle("likelihood ratio test")  + 
    theme(plot.title = element_text(hjust = 0.5)) + 
    theme(legend.title=element_blank()) +
    geom_hline(yintercept=threshold,color="red")

ggsave(paste0(timeless.dir,"/LR_nonneg.pdf"),height=6,width=8)


file.exists("Rplots.pdf")
file.remove("Rplots.pdf")


q()

