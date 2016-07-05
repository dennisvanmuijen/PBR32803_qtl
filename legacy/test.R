
library(qtl)
library(dplyr)
cross <- read.cross(
  format = "csv",
  file = "data/Kamzi_2012.csv",
  genotypes = c("A","H","B"),
  alleles = c("A","B"),
  estimate.map = TRUE,
  F.gen = 7,
  BC.gen = 0
)
# cross <- convert2riself(cross)
cross <- calc.genoprob(cross, step = 5, map.function = "kosambi")


phen <- 5
out.s1perm <- scanone(cross, pheno.col = phen, n.perm = 2500, n.cluster = 4, method = "hk")
out.s1  <- scanone(cross, pheno.col = phen, method = "hk")###############################
thrs <- summary(out.s1perm , alpha=c(0.05, 0.01))
# QTL detect at alpha = 0.05
res <- summary(out.s1, perms=out.s1perm, alpha=0.05)

out.s2 <- cim(cross, pheno.col = phen, method = "hk", n.marcovar = 0)
res2 <- summary(out.s2, perms=out.s1perm, alpha=0.05)
res2




library(dplyr)
library(ggplot2)



  # cross <- geno() %>% convert2riself()
  cross <- geno()
  cross <- cross %>% jittermap()
  cross <- calc.genoprob(cross, step = 5, map.function = "kosambi")

  plotdata <- scanone(cross, pheno.col = 2)
  plotdata$data_id <- tolower(row.names(plotdata))
  plotdata$tooltip <- row.names(plotdata)
  thr <- IMapping()[[2]]
  plotdata$marker <- row.names(IMapping()[[1]])
  plotdata_set <- filter(plotdata, chr %in% 3:6)
  p <- ggplot(data = plotdata_set, aes(x=pos, y=lod, tooltip = tooltip, data_id = data_id))+
    geom_line(aes(group = 1))+geom_rug(data = plotdata_set, sides = "b")+
    facet_wrap(~chr, nrow=3, scales = "free_x") +
    geom_point_interactive(color="orange", size=0.1)+
    theme_dark()
  p
  
 
  