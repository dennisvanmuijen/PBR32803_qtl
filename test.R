
library(qtl)
library(dplyr)
cross <- read.cross(
  format = "csv",
  file = "data/WageningenWeek.csv",
  genotypes = c("A","H","B"),
  alleles = c("A","B"),
  estimate.map = TRUE,
  F.gen = 7,
  BC.gen = 0
)
cross <- convert2riself(cross)
cross <- calc.genoprob(cross, step = 5, map.function = "kosambi")



names(cross$pheno)[44]

ptm <- proc.time()
out.s1perm <- scanone(cross, pheno.col = 44, n.perm = 5000, n.cluster = 3, method = "hk")
proc.time() - ptm
summary(out.s1perm, alpha = c(0.05,0.01))
out.s1perm <- scanone(cross, pheno.col = 44, n.perm = 5000, n.cluster = 3, method = "hk")

out.s1 <- scanone(cross, pheno.col = 44, method = "hk")


plotdata <- out.s1
plotdata$data_id <- tolower(row.names(plotdata))
plotdata$tooltip <- row.names(plotdata)
thr <- summary(out.s1perm, alpha = c(0.05,0.01))

plotdata$marker <- row.names(plotdata)
library(ggplot2)
p <- ggplot(plotdata, aes(x=pos, y=lod)) +
  geom_line(aes(group = 1))+geom_rug(data = plotdata[which(!grepl("loc",plotdata$marker)),], sides = "b")+
  facet_wrap(~chr, nrow=3)+
  geom_point(color="orange", size=0.1)+
  theme_dark() +
  geom_hline(yintercept = thr[1], lwd = 0.5, lty = 2, col = "white") +
  geom_hline(yintercept = thr[2], lwd = 0.5, lty = 2, col = "orange") + 
  # geom_text(aes( 0, thr[1], label = thr[1], vjust = -1), size = 3)
  geom_text(aes(x, y, label=lab),
            data=data.frame(x=10, y=Inf, lab=rep(thr[1],18) %>% round(digits = 2)), vjust=1)

p


cross <- cross %>% jittermap 
cross <- cross %>% calc.genoprob()
out.s1perm <- scanone(cross, pheno.col = 19, n.perm = 1000, n.cluster = 3, method = "hk")

out.s1 <- scanone(cross, pheno.col = 19, method = "hk")
res2 <- summary(out.s1, perms=out.s1perm, alpha=0.05)



# res2 <- summary(out.s2)
# if(attributes(out.s2)$pLOD != 0){

  mqtl <- makeqtl(cross,res2[,1], res2[,2], what="prob")
  QTLnames <- unlist(mqtl$altname)
  form <- paste(QTLnames, collapse="+")
  form <- paste("y~",form, sep="")
  fit <- fitqtl(cross, pheno.col = 19, qtl = mqtl, method = "hk", model = "normal",
                formula = form, get.ests = T)
  fit %>% summary
  
  
  
  
  