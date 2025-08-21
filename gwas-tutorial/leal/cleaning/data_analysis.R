# Sex Check
sexcheck = read.table("GWAS_sex_checking.sexcheck", header=T)
names(sexcheck)
sex_problem = sexcheck[which(sexcheck$STATUS=="PROBLEM"),]
sex_problem

# Duplicates
dups = read.table("duplicates.genome", header = T)
problem_pairs = dups[which(dups$PI_HAT > 0.4),]
problem_pairs

problem_pairs = dups[which(dups$PI_HAT > 0.05),]
myvars = c("FID1", "IID1", "FID2", "IID2", "PI_HAT")
problem_pairs[myvars]

# Excess Zygosity
Dataset <- read.table("plink.het", header=TRUE, sep="", na.strings="NA",
                      dec=".", strip.white=TRUE)
mean(Dataset$F)
sd(Dataset$F)

jpeg("hist.jpeg", height=1000, width=1000)
hist(scale(Dataset$F), xlim=c(-4,4))
dev.off()

# Hardy-Weinberg Equilibrium
hardy = read.table("plink.hwe", header = T)
names(hardy)
hwe_prob = hardy[which(hardy$P < 0.0000009),]
hwe_prob


# Multidimensional Scaling
mydata = read.table("mds_components.txt", header=T)
mydata$pch[mydata$Group==1 ] <-15
mydata$pch[mydata$Group==2 ] <-16
mydata$pch[mydata$Group==3 ] <-2
jpeg("mds.jpeg", height=1000, width=1000)
plot(mydata$C1, mydata$C2 ,pch=mydata$pch)
dev.off()

# QQ Plot
broadqq <-function(pvals, title)
{
  observed <- sort(pvals)
  lobs <- -(log10(observed))
  expected <- c(1:length(observed))
  lexp <- -(log10(expected / (length(expected)+1)))
  plot(c(0,7), c(0,7), col="red", lwd=3, type="l", xlab="Expected (-logP)", ylab="Observed (-logP)",
       xlim=c(0,max(lobs)), ylim=c(0,max(lobs)), las=1, xaxs="i", yaxs="i", bty="l", main = title)
  points(lexp, lobs, pch=23, cex=.4, bg="black") }
jpeg("qqplot_compare.jpeg", height=1000, width=1000)
par(mfrow=c(2,1))
aff_unadj<-read.table("unadj.assoc.logistic", header=TRUE)
aff_unadj.add.p<-aff_unadj[aff_unadj$TEST==c("ADD"),]$P
broadqq(aff_unadj.add.p,"Some Trait Unadjusted")
aff_C1C2<-read.table("C1-C2.assoc.logistic", header=TRUE)
aff_C1C2.add.p<-aff_C1C2[aff_C1C2$TEST==c("ADD"),]$P
broadqq(aff_C1C2.add.p, "Some Trait Adjusted")
dev.off()

# SNP Search
gws_unadj = aff_unadj[which(aff_unadj$P < 0.0000001),]
gws_unadj
gws_adjusted = aff_C1C2[which(aff_C1C2$P < 0.0000001),]
gws_adjusted