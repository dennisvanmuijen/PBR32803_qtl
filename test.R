

cross <- read.cross(
  format = "csv",
  file = "data/cross_pheno.csv",
  genotypes = c("A","H","B"),
  alleles = c("A","B"),
  estimate.map = TRUE,
  F.gen = 2,
  BC.gen = 0
)
cross <- calc.genoprob(cross, step = 5, map.function = "kosambi")
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
