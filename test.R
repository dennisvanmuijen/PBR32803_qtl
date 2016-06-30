
library(qtl)
library(dplyr)
cross <- read.cross(
  format = "csv",
  file = "data/cross_pheno.csv",
  genotypes = c("A","H","B"),
  alleles = c("A","B"),
  estimate.map = TRUE,
  F.gen = 8,
  BC.gen = 0
)
cross <- calc.genoprob(cross, step = 5, map.function = "kosambi")
cross <- convert2riself(cross)



cross <- fill.geno(cross)
cross <- jittermap(cross)
cofactors <- mqmautocofactors(cross ,50)	# Set 15 Cofactors



result <- mqmscan(cross,cofactors, pheno.col = 2)	# Backward model selection

(cross$pheno %>% names)[8]
phen <- 8
result <- mqmscan(cross,cofactors, pheno.col = phen)
s1 <- scanone(cross, pheno.col = phen, method = "hk")

plot(s1, result, col = c("black"," orange"), lty = 1, lwd = 2)


mqmres <- mqmgetmodel(result)
mqmres 

qtl <- makeqtl(cross %>% calc.genoprob(), mqmres$chr, mqmres$pos, what = "prob")


phen <- 16
sw.out <- stepwiseqtl(cross %>% calc.genoprob(), pheno.col = phen, max.qtl = 5, method = "hk", additive.only = T,
         penalties = 3)
plot(sw.out)
summary(sw.out)
plotLodProfile(sw.out)


names(cross$pheno)










qtl <- makeqtl(cross %>% calc.genoprob(), c(5,8,9), c(95,56.3,73.8), what = "prob")


fit <- fitqtl(cross, pheno.col = phen, qtl = qtl, formula ="y~Q1+Q2+Q3")
summary(fit)




fit <- fitqtl(cross, pheno.col = phen, qtl = qtl, formula = paste0("y~",paste(mqmres$altname, collapse = "+")))
summary(fit)




plot(scanone(cross, pheno.col = 1))

plotdata <- scanone(cross, pheno.col = 1, n.perm = 1000, n.cluster = 4)
plotdata2 <- scanone(cross, pheno.col = 1)

summary(plotdata , alpha=c(0.05, 0.10))
# Results above the 0.05 threshold
res <- summary(plotdata2, perms=plotdata, alpha=0.5)
res[1,1]
res

mqtl
mqtl <- makeqtl(cross,res[,1], res[,2], what="prob")
form <- 

fit <- fitqtl(cross, pheno.col = 8, qtl = mqtl, method = "hk", model = "normal",
              formula = form, get.ests = T)
summary(fit)


names(cross$pheno)

?stepwiseqtl


cross <- jittermap(cross)
cross <- calc.genoprob(cross)
pens <- c(3.913854, 5.836726, 3.194303 )



mqtl <- makeqtl(cross,test$chr, test$pos, what="prob")

test <- stepwiseqtl(cross, pheno.col = 9, max.qtl = 3, 
                    method = "hk", additive.only = F,
                    penalties = pens)
attributes(test)$formula



fit <- fitqtl(cross, pheno.col = 9, qtl = mqtl, method = "hk", model = "normal",
              formula = attributes(test)$formula, get.ests = T)
summary(fit)



?stepwiseqtl


pens <- calc.penalties(out.2dim, alpha = 0.1)

out.2dim <- scantwo(cross, method="hk", pheno.col = 8, n.perm = 1000, n.cluster = 4)





QTLnames <- unlist(mqtl$altname)
form <- paste(QTLnames, collapse="+")
form <- paste("y~",form, sep="")
