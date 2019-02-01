library(reshape)
library(vegan)
library(wesanderson)
# Read data
LTcover2013 <- read.csv("C:/Users/lgherar1/GitHub/Microbial-vss-Plant-diversity-Resp/data/LTcover2013.csv")
treatments<-read.csv("C:/Users/lgherar1/GitHub/Microbial-vss-Plant-diversity-Resp/data/trt.csv")

# Format data so species are in rows and plots are in columns
df<-melt(LTcover2013,id="plot")
df<-cast(df,variable~plot)
rownames(df)<-df[,1]
df<-df[,-1]
trans.df <- t(df)
#Site x Species matrix; now begin playing in Vegan
#Make a vector for coloring points
colvec_nitrogen <- c(wes_palette("Darjeeling1")[5], wes_palette("Darjeeling1")[3], wes_palette("Darjeeling1")[3],
                     wes_palette("Darjeeling1")[3], wes_palette("Darjeeling1")[5], wes_palette("Darjeeling1")[3],
                     wes_palette("Darjeeling1")[5], wes_palette("Darjeeling1")[3], wes_palette("Darjeeling1")[5],
                     wes_palette("Darjeeling1")[5], wes_palette("Darjeeling1")[3], wes_palette("Darjeeling1")[5],
                     wes_palette("Darjeeling1")[5], wes_palette("Darjeeling1")[3], wes_palette("Darjeeling1")[5],
                     wes_palette("Darjeeling1")[5], wes_palette("Darjeeling1")[5], wes_palette("Darjeeling1")[3],
                     wes_palette("Darjeeling1")[5], wes_palette("Darjeeling1")[3], wes_palette("Darjeeling1")[3],
                     wes_palette("Darjeeling1")[5], wes_palette("Darjeeling1")[3], wes_palette("Darjeeling1")[3],
                     wes_palette("Darjeeling1")[3], wes_palette("Darjeeling1")[5], wes_palette("Darjeeling1")[3],
                     wes_palette("Darjeeling1")[5], wes_palette("Darjeeling1")[5], wes_palette("Darjeeling1")[5],
                     wes_palette("Darjeeling1")[3], wes_palette("Darjeeling1")[3], wes_palette("Darjeeling1")[5],
                     wes_palette("Darjeeling1")[3], wes_palette("Darjeeling1")[5], wes_palette("Darjeeling1")[3])

nitrogen<-treatments[,2]
#PCoA with bray curtis dissimilarity
pcoa.res <- dbrda(trans.df ~ 1, distance = "bray")
summary(pcoa.res)

png("figures/Plant_Nitrogen_all_water_treatments.png",res = 300, pointsize = 12,height=5.5,width=8,units="in")
ordiplot(pcoa.res, type = "none", las = 1, main = "PLANT - PCoA with Bray-Curtis - Nitrogen (95%)",xlab = "MDS1 (26%)", ylab = "MDS2 (21%)")

#Ellipses are put onto the image before points, so points sit on top of confidence intervals and are more visable
ordiellipse(pcoa.res, nitrogen, draw="polygon", border = NA, kind = "sd", conf = 0.95, label = FALSE, col = c(wes_palette("Darjeeling1")[5],wes_palette("Darjeeling1")[3]), lty = c(2))

points(pcoa.res, display = "sites", col = colvec_nitrogen, pch = 20, cex = 2)

legend(x="bottomleft", pch = c(16, 16, 16), col = c(wes_palette("Darjeeling1")[5], wes_palette("Darjeeling1")[3]), legend = c("No Nitrogen", "Nitrogen"), bty = "n")
dev.off()

#Site x Species matrix; now begin playing in Vegan
colvec_water <- c(wes_palette("Darjeeling1")[3], wes_palette("Darjeeling1")[3], wes_palette("Darjeeling1")[5], 
                  wes_palette("Darjeeling1")[2], wes_palette("Darjeeling1")[2], wes_palette("Darjeeling1")[5], 
                  wes_palette("Darjeeling1")[3], wes_palette("Darjeeling1")[3], wes_palette("Darjeeling1")[2], 
                  wes_palette("Darjeeling1")[5], wes_palette("Darjeeling1")[2], wes_palette("Darjeeling1")[5],
                  wes_palette("Darjeeling1")[3], wes_palette("Darjeeling1")[2], wes_palette("Darjeeling1")[2], 
                  wes_palette("Darjeeling1")[3], wes_palette("Darjeeling1")[5], wes_palette("Darjeeling1")[2],
                  wes_palette("Darjeeling1")[2], wes_palette("Darjeeling1")[5], wes_palette("Darjeeling1")[3], 
                  wes_palette("Darjeeling1")[5], wes_palette("Darjeeling1")[3], wes_palette("Darjeeling1")[5],
                  wes_palette("Darjeeling1")[3], wes_palette("Darjeeling1")[2], wes_palette("Darjeeling1")[5], 
                  wes_palette("Darjeeling1")[2], wes_palette("Darjeeling1")[3], wes_palette("Darjeeling1")[5],
                  wes_palette("Darjeeling1")[5], wes_palette("Darjeeling1")[2], wes_palette("Darjeeling1")[3], 
                  wes_palette("Darjeeling1")[2], wes_palette("Darjeeling1")[5], wes_palette("Darjeeling1")[3])


water<-treatments[,3]
#PCoA with bray curtis dissimilarity
pcoa.res <- dbrda(trans.df ~ 1, distance = "bray")

summary(pcoa.res)
png("figures/Plant_Water.png",res = 300, pointsize = 12,height=5.5,width=8,units="in")
ordiplot(pcoa.res, type = "none", las = 1, main = "PCoA with Bray-Curtis - Water (95%)", xlab = "MDS1 (8.3%)", ylab = "MDS2 (4.5%)")
mtext(text ="ALL NITROGEN TREATMENTS",side = 3)
ordiellipse(pcoa.res, water, draw="polygon", border = NA, kind = "sd", conf = 0.95, label = FALSE, col = c(wes_palette("Darjeeling1")[2],wes_palette("Darjeeling1")[3], wes_palette("Darjeeling1")[5]), lty = c(2))

points(pcoa.res, display = "sites", col = colvec_water, pch = 20, cex = 2)

legend(x="bottomleft", pch = c(16, 16, 16), col = c(wes_palette("Darjeeling1")[5], wes_palette("Darjeeling1")[3], wes_palette("Darjeeling1")[2]), legend = c("Irrigated", "Drought", "Control"), bty = "n")
dev.off()



