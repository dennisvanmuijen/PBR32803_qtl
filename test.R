

cross <- read.cross(
  format = "csv",
  file = "data/cross.csv",
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

fit <- fitqtl(cross, pheno.col = 1, qtl = mqtl, method = "hk", model = "normal",
              formula = form, get.ests = T)
summary(fit)




QTLnames <- unlist(mqtl$altname)
form <- paste(QTLnames, collapse="+")
form <- paste("y~",form, sep="")
