setwd("/home/sonia/Documentos/Covid19/shiny-app/toy/SYBR-Biorad/")

a <- read.table("20200610_SYBR_JNG_P1 -  Melt Curve RFU Results_MOD_H30.csv", header = T, sep = ",")
a$X <-NULL
b <- read.table("20200610_SYBR_JNG_P1 -  Melt Curve RFU Results_MOD_N.csv", header = T, sep = ",")
b$X <- NULL
c <- read.table("20200610_SYBR_JNG_P1 -  Melt Curve RFU Results_MOD_RdRp.csv", header = T, sep = ",")
c$X <- NULL
d <- read.table("20200610_SYBR_JNG_P1 -  Melt Curve RFU Results_MOD_S.csv", header = T, sep = ",")
d$X <- NULL

df <- merge(a,b,by.x = "Cycle")
df2 <- merge(df,c, by.x = "Cycle")
df3 <- merge(df2,d, by.x = "Cycle")

df3$Cycle <- rownames(df3)
df3$Cycle <- NULL
at <- t(df3)
atm <- melt(at)

atm2 <- atm %>%
  group_split(Var2)

def <- lapply(atm2, function(x){
  x <- as.data.frame(x)
  x$let <- substr(x$Var1, 0,1)
  x$num <- substr(x$Var1, 2, length(x$Var1))
  x <- x[order(x$let, as.numeric(as.character(x$num))),]
})

def2 <- do.call(rbind, def)

write.table(def2$value, "forApplied.csv",sep = ",", dec = ".", row.names = FALSE)
  