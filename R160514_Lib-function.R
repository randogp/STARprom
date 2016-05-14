# GR 21.08.2015 Lib8 analysis
# This script does the following:
# import promoter - bc association
# import RNAseq libraries
#       lib3 = vinblastine duplicate 1, marked as "a"
#       lib4 = jasplakinolide duplicate 1, marked as "b", bc are in reverse complement.
#       lib8 = vinblastine and jasplakinolide duplicate 2, marked as "c"
#
# import RNAseq libraries
#       lib3 = vinblastine duplicate 1, marked as "a"
#       lib4 = jasplakinolide duplicate 1, marked as "b", bc are in reverse complement.
#       lib8 = vinblastine and jasplakinolide duplicate 2, marked as "c"
#
# adapt bcrevc in lib4 with function str_rev
#
# prune libs to have same number of barcodes (rows)
# merge the read tables together
#
# normalize reads of each sample by their median and scale each sample to a global median of 500
# --> using function median.norm
# 
#        la media globale 500 ?? inventata, per farla bene, implementare:
#        cmed <- sapply(lib8[,2:ncol(lib8)], median, na.rm=TRUE) #median of each column
#        lib8.norm[,2:ncol(lib8.norm)] <- lib8.norm[,2:ncol(lib8.norm)] * median(cmed) #rescale normalized values to global median
#
# plot normalization boxplots
#
# define groups
# ctrl.a <- c("D0a", "D1a", "D2a", "D8a", "V0a") #library 3
# vinb.a <- c("V1a", "V2a", "V4a", "V8a")
# ctrl.b <- c("D0b", "D1b", "D2b", "D8b", "J0b") #library 4
# jas.b <- c("J1b", "J2b", "J4b", "J8b")
# ctrl.c <- c("D0c", "D1c", "D2c", "D4c", "D8c", "J0c", "V0c") #library 8
# vinb.c <- c("V1c", "V2c", "V4c", "V8c")
# jas.c <- c("J1c", "J2c", "J4c", "J8c")
#
# calculate fold changes by dividing each read by the mean of their control reads.
# --> using function fold.change
#
# merge tables of the fold changes calculated on the corresponding controls
#
# sort barcodes by decrescent order of mean drug / mean controls
# --> using function rank.fc
#
# transform fold changes into log 2 values to have equidistant values on activation (fc >1) and repression (fc < 1).
#
# plot a heatmap of all values
#
# plot a heatmap of tail values (top 20, bottom 20)
# --> using function cut.tails

# test if DMSO is unducing something by by calculating the fold changefor dmso / time 0
# --> a slightly different rank.fc.dmso is required.
# For this control, the groups are
#  time0.a <- c("D0a", "V0a") #library 3
#  time0.b <- c("D0b", "J0b") #library 4
#  time0.c <- c("D0c", "J0c", "V0c") #library 8
#  dmso.a <- c("D1a", "D2a", "D8a") #library 3
#  dmso.b <- c("D1b", "D2b", "D8b") #library 4
#  dmso.c <- c("D1c", "D2c", "D4c", "D8c") #library 8


# at the end of the script, ls() gives the following items 
#
#[1] "bks"              "cols"             "ctrl.a"           "ctrl.b"           "ctrl.c"           "cut.tails"       
#[7] "dmso.a"           "dmso.b"           "dmso.c"           "fc"               "fc.a"             "fc.b"            
#[13] "fc.c"             "fc.dmso.a"        "fc.dmso.b"        "fc.dmso.c"        "fold.change"      "jas.b"           
#[19] "jas.c"            "lib"              "lib.norm"         "lib3"             "lib4"             "lib8"            
#[25] "list"             "log.dmso.fc"      "log.fc"           "median.norm"      "prom"             "rank.a"          
#[31] "rank.b"           "rank.b.back"      "rank.c"           "rank.d"           "rank.dmso.a"      "rank.dmso.b"     
#[37] "rank.dmso.c"      "rank.dmso.d"      "rank.fc"          "rank.fc.dmso"     "table.dmso.fc"    "table.fc"        
#[43] "tails.dmso.fc"    "tails.fc"         "tails.table"      "tails.table.dmso" "time0.a"          "time0.b"         
#[49] "time0.c"          "vinb.a"           "vinb.c" 



rm(list=ls())

#setwd('/Users/Santa/Google Drive/Schibler Lab/STARpromBC/R-scripts_STARprom') #windows
setwd('/Users/randogp/Google Drive/Schibler Lab/STARpromBC/R-scripts_STARprom') #MacBook Air
setwd('/Users/admin/Google Drive/Schibler Lab/STARpromBC/R-scripts_STARprom') #iMac21

#import promoter - bc association
prom<-read.delim("Galaxy-[First-5000-library-BC-prom_SUPP].tabular", header = FALSE, as.is = TRUE, quote = "\'")
colnames(prom)<-c("ID", "DNAreads", "bc", "prom-miseq", "reads_pacbio", "bc_pacbio", "prom_pacbio")
prom<-prom[c("bc", "ID", "DNAreads", "prom_pacbio")]



lib3<-read.delim("Galaxy-[Lib3-all-reads].tabular", header = FALSE, as.is = TRUE, quote = "\'")
colnames(lib3)<-c("bc","V0a", "Vh", "V1a", "V2a", "V4a", "V8a", "D0a", "D1a", "D2a", "D8a") #Vh is half-hour will not be considered
lib3 <- lib3[c("bc","D0a", "D1a", "D2a", "D8a", "V0a", "V1a", "V2a", "V4a", "V8a")] #change column order
lib3 <- lib3[lib3$bc %in% prom$bc, ] # the 3110 rows in lib3 whose bc has a corresponding promoter

lib4<-read.delim("Galaxy-[Lib4-all-reads].tabular", header = FALSE, as.is = TRUE, quote = "\'")
colnames(lib4)<-c("bcrevc", "D8b", "J0b", "Jh", "J1b", "J2b", "J4b", "J8b", "D0b", "D1b", "D2b")
lib4 <- lib4[c("bcrevc","D0b", "D1b", "D2b", "D8b", "J0b", "J1b", "J2b", "J4b", "J8b")] #change column order
lib4 <- cbind(bc = 0, lib4) #add an empty column named BC
lib4$bc <- chartr("ATGC","TACG",lib4$bcrevc) #write complement

str_rev <- function(x) {
  sapply( x, function(xx) { 
    intToUtf8( rev( utf8ToInt( xx ) ) )
  } )
}
lib4$bc <- str_rev(lib4$bc) #write reverse using str-rev function
lib4 <- subset(lib4, select = -c(bcrevc)) #remove column with BCrevc
lib4 <- lib4[lib4$bc %in% prom$bc, ] # the 3156 rows

rm(str_rev)

lib8<-read.delim("Galaxy-[Lib8-all-reads].tabular", header = FALSE, as.is = TRUE, quote = "\'")
colnames(lib8)<-c("bc","D0c", "D1c", "D2c", "D4c", "D8c", "J0c", "J1c", "J2c", "J4c", "J8c", "V0c", "V1c", "V2c", "V4c", "V8c")
lib8 <- lib8[lib8$bc %in% prom$bc, ] # the 3182 rows

lib8 <- lib8[lib8$bc %in% lib3$bc, ] # the 3106 rows
lib8 <- lib8[lib8$bc %in% lib4$bc, ] # the 3065 rows
lib4 <- lib4[lib4$bc %in% lib3$bc, ] # the 3069 rows
lib4 <- lib4[lib4$bc %in% lib8$bc, ] # the 3065 rows
lib3 <- lib3[lib3$bc %in% lib4$bc, ] # the 3065 rows
lib3 <- lib3[lib3$bc %in% lib8$bc, ] # the 3065 rows


# Group definitions
ctrl.a <- c("D0a", "D1a", "D2a", "D8a", "V0a") #library 3
vinb.a <- c("V1a", "V2a", "V4a", "V8a")

ctrl.b <- c("D0b", "D1b", "D2b", "D8b", "J0b") #library 4
jas.b <- c("J1b", "J2b", "J4b", "J8b")

ctrl.c <- c("D0c", "D1c", "D2c", "D4c", "D8c", "J0c", "V0c") #library 8
vinb.c <- c("V1c", "V2c", "V4c", "V8c")
jas.c <- c("J1c", "J2c", "J4c", "J8c")

# I try experiment: merge 3 tables

lib <- merge(lib3, lib4, by="bc", all=FALSE)
lib <- merge(lib, lib8, by="bc", all=FALSE)
lib[lib=="."]<-1 #replace . in few cells with 1.
#lib <- lib[complete.cases(lib),] #alternative way: replace . with NA and then remove the raw.

lib[,2:ncol(lib)] <- as.numeric(as.character(unlist(lib[,2:ncol(lib)]))) #required for some graphs


# median normalization
#        funziona ma la media globale 500 ?? inventata, per farla bene, implementare:
#        cmed <- sapply(lib8[,2:ncol(lib8)], median, na.rm=TRUE) #median of each column
#        lib8.norm[,2:ncol(lib8.norm)] <- lib8.norm[,2:ncol(lib8.norm)] * median(cmed) #rescale normalized values to global median

median.norm <- function (x) {
  x[,2:ncol(x)] <- as.numeric(as.character(unlist(x[,2:ncol(x)]))) #change bc columns from character to numeric (R faq 7.10) OK
  x[,2:ncol(x)] <- sweep(x[,2:ncol(x)], 2,sapply(x[,2:ncol(x)], median, na.rm=TRUE), `/`) #divide each bc column by its median
  x[,2:ncol(x)] <- x[,2:ncol(x)] * 500 #rescale normalized values to global median
}

lib.norm <- lib
lib.norm[,2:ncol(lib)] <- median.norm(lib)

#plot normalization OK
par(mfrow=c(1,2))
boxplot(lib[,2:ncol(lib)], outline = FALSE, main = "Before normalization", xlab = "samples", ylab = "reads") #boxplot before normalization
boxplot(lib.norm[,2:ncol(lib.norm)], outline = FALSE, main = "After median norm", xlab = "samples", ylab = "reads") #boxplot after normalization
par(mfrow=c(1,1))

#correlation dmso

lib.dmso<-lib.norm[c("D0b","D0c","D1b","D1c","D2b","D2c","D8b","D8c")]
log.lib.dmso<-log(lib.dmso,10)
par(mfrow=c(2,2))
plot(log.lib.dmso[c("D0b","D0c")])
plot(log.lib.dmso[c("D1b","D1c")])
plot(log.lib.dmso[c("D2b","D2c")])
plot(log.lib.dmso[c("D8b","D8c")])
par(mfrow=c(1,1))

plot(log.lib.dmso[,1])

write.table(log.lib.dmso, "160307_Fig-S1C-RAW-DATA.txt", sep="\t")

#Fold Induction

fold.change <- function (lib.norm,ctrl) {
  lib.norm[,2:ncol(lib.norm)] <- lib.norm[,2:ncol(lib.norm)] / rowMeans(lib.norm[,ctrl]) #divide each bc column by the control mean
  lib.norm[,2:ncol(lib.norm)]<-signif(lib.norm[,2:ncol(lib.norm)], digits = 2)
}

fc.a <- lib.norm
fc.b <- lib.norm
fc.c <- lib.norm
fc.a[,2:ncol(lib.norm)] <- fold.change(lib.norm,ctrl.a)
fc.b[,2:ncol(lib.norm)] <- fold.change(lib.norm,ctrl.b)
fc.c[,2:ncol(lib.norm)] <- fold.change(lib.norm,ctrl.c)



#Ranking

rank.fc <- function (lib.fc,drug,ctrl) {
  rank <- data.frame(rowMeans(lib.fc[,drug])/rowMeans(lib.fc[,ctrl])) #the rank is higher when drug>controls
  colnames(rank) <- paste("rank.", substr(ctrl[1], 3, nchar(ctrl[1])), sep ="") #identify rank column with last letter of ctrl (a,b,c)
  lib.fc <- cbind(rank,lib.fc)
  lib.fc <- lib.fc[order(-rank),] #sort decrescent order by rank
  lib.fc <- cbind(lib.fc[,1:2],lib.fc[ctrl],lib.fc[drug])
}

rank.a <- rank.fc(fc.a,vinb.a,ctrl.a) #rank fc.a according to ratio of the row means drug/ctrl
rank.b <- rank.fc(fc.b,jas.b,ctrl.b)
rank.c <- rank.fc(fc.c,jas.c,ctrl.c)
rank.d <- rank.fc(fc.c,vinb.c,ctrl.c)

#Create table of all fold changes, calculated on the corresponding controls

table.fc <- merge(rank.b, rank.a, by="bc", all=FALSE)
table.fc <- merge(table.fc, rank.c, by="bc", all=FALSE)
table.fc <- merge(table.fc, rank.d, by="bc", all=FALSE)

table.fc <- table.fc[c("bc", "rank.b", "rank.a", "rank.c.x", "rank.c.y",  "D0b",  "D1b", "D2b", "D8b",  "J0b", "J1b", "J2b",  "J4b", "J8b", "D0a", "D1a", "D2a", "D8a", "V0a",  "V1a", "V2a", "V4a",  "V8a", "D0c.x", "D1c.x", "D2c.x", "D4c.x", "D8c.x", "J0c.x", "J1c", "J2c", "J4c", "J8c", "D0c.y", "D1c.y", "D2c.y", "D4c.y", "D8c.y", "V0c.y", "V1c",  "V2c", "V4c", "V8c")]
table.fc[,2:ncol(table.fc)] <- as.numeric(as.character(unlist(table.fc[,2:ncol(table.fc)]))) #change bc columns from character to numeric


#Heatmap preparation
table.fc <- table.fc[ order(-table.fc[,"rank.c.x"]), ] #order induced by jasp lib8
#write.table(table.fc, "150825_FC-All-libraries-ranked-by-Ueli.txt", sep="\t")
fc <- c("D0b",  "D1b", "D2b", "D8b",  "J0b", "J1b", "J2b", "J4b", "J8b", "D0a", "D1a", "D2a", "D8a", "V0a",  "V1a", "V2a", "V4a",  "V8a", "D0c.x", "D1c.x", "D2c.x", "D4c.x", "D8c.x", "J0c.x", "J1c", "J2c", "J4c", "J8c", "D0c.y", "D1c.y", "D2c.y", "D4c.y", "D8c.y", "V0c.y", "V1c",  "V2c", "V4c", "V8c")
#fc <- c(ctrl.a, vinb.a, ctrl.b, jas.b, ctrl.c, jas.c, vinb.c)


log.fc <- table.fc
log.fc[,fc] <- log(log.fc[,fc],2) #uso il vettore fc per dire tutto


library(gplots)

draw_heatmap <- function(log.fc,group,scale) {
  cols <- c(colorRampPalette(c("cornflowerblue", "white"))(30), colorRampPalette(c("white", "red"))(30)) #color palette 30 between blue and white, 30 between white and red
  bks <- c(seq(-scale,0,length=30),seq(0,scale,length=31)) #set the breaks (breaks should be 1 more than colors)
  heatmap.2(as.matrix(log.fc[,group]), Rowv=NULL, Colv=NULL, dendrogram="none", trace="none",
            key=FALSE, lhei=c(0.1,1), lwid=c(0.1,1), margins=c(5,12), col=cols, breaks=bks) #ignore error of figures too large  
}

draw_heatmap.dendro <- function(log.fc,group,scale) {
  cols <- c(colorRampPalette(c("cornflowerblue", "white"))(30), colorRampPalette(c("white", "red"))(30)) #color palette 30 between blue and white, 30 between white and red
  bks <- c(seq(-scale,0,length=30),seq(0,scale,length=31)) #set the breaks (breaks should be 1 more than colors)
  heatmap.2(as.matrix(log.fc[,group]), Rowv=TRUE, Colv=NULL, dendrogram="row", trace="none",
            key=FALSE, lhei=c(0.1,1), lwid=c(0.1,1), margins=c(5,12), col=cols, breaks=bks) #ignore error of figures too large  
}

# I select the best palette (logFC 8 as max looks good)


pdf(file="pdf150827_heatmap-prova_scale4.pdf",width=4.86,height=7.01)
draw_heatmap(tails.fc,fc,4)
dev.off()

pdf(file="pdf150827_heatmap-prova_scale6.pdf",width=4.86,height=7.01)
draw_heatmap(tails.fc,fc,6)
dev.off()

pdf(file="pdf150827_heatmap-prova_scale8.pdf",width=4.86,height=7.01)
draw_heatmap(tails.fc,fc,8)
dev.off()

pdf(file="pdf150827_heatmap-prova_scale10.pdf",width=4.86,height=7.01)
draw_heatmap(tails.fc,fc,10)
dev.off()

pdf(file="pdf150831_heatmap-completa_scale8.pdf",width=4.86,height=7.01)
draw_heatmap(log.fc,fc,8)
dev.off()

pdf(file="pdf150831_heatmap-completa_dendro8.pdf",width=4.86,height=7.01)
draw_heatmap.dendro(log.fc,fc,8)
dev.off()

#Top-low 20 heatmap ranked by Jas lib8

cut.tails <- function (x,y) {
  top <- x[1:y,]
  bottom <- x[(nrow(x)-y+1):(nrow(x)),]
  x<-rbind(top,bottom)
}

cut.top <- function (x,y) {
  top <- x[1:y,]
  x<-rbind(top)
}

tails.table <- cut.tails(table.fc,20)
tails.fc <- tails.table
tails.fc[,fc] <- log(tails.table[,fc],2)


heatmap.2(as.matrix(tails.fc[,fc]), Rowv=NULL, Colv=NULL, dendrogram="none", trace="none",
          key=FALSE, lhei=c(0.1,1), lwid=c(0.1,1), margins=c(5,12), col=cols, breaks=bks) #ignore error of figures too large

# Heatmap ranked by Vinb 2
rm(log.fc)
rm(tails.table)
rm(tails.fc)

table.fc <- table.fc[ order(-table.fc[,"rank.c.y"]), ] #order induced by vinb lib8
log.fc <- table.fc
log.fc[,fc] <- log(log.fc[,fc],2) #uso il vettore fc per dire tutto
heatmap.2(as.matrix(log.fc[,fc]), Rowv=NULL, Colv=NULL, dendrogram="none", trace="none",
          key=FALSE, lhei=c(0.1,1), lwid=c(0.1,1), margins=c(5,12), col=cols, breaks=bks) #ignore error of figures too large


tails.table <- cut.tails(table.fc,20)
tails.fc <- tails.table
tails.fc[,fc] <- log(tails.table[,fc],2)


heatmap.2(as.matrix(tails.fc[,fc]), Rowv=NULL, Colv=NULL, dendrogram="none", trace="none",
          key=FALSE, lhei=c(0.1,1), lwid=c(0.1,1), margins=c(5,12), col=cols, breaks=bks) #ignore error of figures too large


# Control: heatmap ranked on DMSO induction vs time 0.

time0.a <- c("D0a", "V0a") #library 3
time0.b <- c("D0b", "J0b") #library 4
time0.c <- c("D0c", "J0c", "V0c") #library 8

dmso.a <- c("D1a", "D2a", "D8a") #library 3
dmso.b <- c("D1b", "D2b", "D8b") #library 4
dmso.c <- c("D1c", "D2c", "D4c", "D8c") #library 8

fc.dmso.a <- lib.norm
fc.dmso.b <- lib.norm
fc.dmso.c <- lib.norm

fc.dmso.a[,2:ncol(lib.norm)] <- fold.change(lib.norm,time0.a)
fc.dmso.b[,2:ncol(lib.norm)] <- fold.change(lib.norm,time0.b)
fc.dmso.c[,2:ncol(lib.norm)] <- fold.change(lib.norm,time0.c)


rank.fc.dmso <- function (lib.fc,dmso,t0,drug,ctrl) {
  rank <- data.frame(rowMeans(lib.fc[,dmso])/rowMeans(lib.fc[,t0])) #the rank is higher when drug>controls
  colnames(rank) <- paste("rank.", substr(t0[1], 3, nchar(ctrl[1])), sep ="") #identify rank column with last letter of ctrl (a,b,c)
  lib.fc <- cbind(rank,lib.fc)
  lib.fc <- lib.fc[order(-rank),] #sort decrescent order by rank
  lib.fc <- cbind(lib.fc[,1:2],lib.fc[ctrl],lib.fc[drug])
}

rank.dmso.a <- rank.fc.dmso(fc.dmso.a,dmso.a,time0.a,vinb.a,ctrl.a) #rank fc.a according to ratio of the row means drug/ctrl
rank.dmso.b <- rank.fc.dmso(fc.dmso.b,dmso.b,time0.b,jas.b,ctrl.b)
rank.dmso.c <- rank.fc.dmso(fc.dmso.c,dmso.c,time0.c,vinb.c,ctrl.c)
rank.dmso.d <- rank.fc.dmso(fc.dmso.c,dmso.c,time0.c,jas.c,ctrl.c)

table.dmso.fc <- merge(rank.dmso.b, rank.dmso.a, by="bc", all=FALSE)
table.dmso.fc <- merge(table.dmso.fc, rank.dmso.c, by="bc", all=FALSE)
table.dmso.fc <- merge(table.dmso.fc, rank.dmso.d, by="bc", all=FALSE)

table.dmso.fc <- table.dmso.fc[c("bc", "rank.b", "rank.a", "rank.c.x", "rank.c.y",  "D0b",  "D1b", "D2b", "D8b",  "J0b", "J1b", "J2b",  "J4b", "J8b", "D0a", "D1a", "D2a", "D8a", "V0a",  "V1a", "V2a", "V4a",  "V8a", "D0c.x", "D1c.x", "D2c.x", "D4c.x", "D8c.x", "J0c.x", "J1c", "J2c", "J4c", "J8c", "D0c.y", "D1c.y", "D2c.y", "D4c.y", "D8c.y", "V0c.y", "V1c",  "V2c", "V4c", "V8c")]
table.dmso.fc[,2:ncol(table.dmso.fc)] <- as.numeric(as.character(unlist(table.dmso.fc[,2:ncol(table.dmso.fc)]))) #change bc columns from character to numeric

table.dmso.fc <- table.dmso.fc[ order(-table.dmso.fc[,"rank.c.y"]), ] #order induced by vinb lib8


log.dmso.fc <- table.dmso.fc
log.dmso.fc[,fc] <- log(log.dmso.fc[,fc],2) #uso il vettore fc per dire tutto
heatmap.2(as.matrix(log.dmso.fc[,fc]), Rowv=NULL, Colv=NULL, dendrogram="none", trace="none",
          key=FALSE, lhei=c(0.1,1), lwid=c(0.1,1), margins=c(5,12), col=cols, breaks=bks) #ignore error of figures too large


tails.table.dmso <- cut.tails(table.dmso.fc,20)
tails.dmso.fc <- tails.table.dmso
tails.dmso.fc[,fc] <- log(tails.table.dmso[,fc],2)


heatmap.2(as.matrix(tails.dmso.fc[,fc]), Rowv=NULL, Colv=NULL, dendrogram="none", trace="none",
          key=FALSE, lhei=c(0.1,1), lwid=c(0.1,1), margins=c(5,12), col=cols, breaks=bks) #ignore error of figures too large

#heatmaps for Ueli 1.9.15

table.fc<-table.fc[ order(-table.fc[,"rank.c.x"]), ]
jas.table <- cut.top(table.fc,20)
jas.fc <- jas.table
jas.fc[,fc] <- log(jas.table[,fc],2)

pdf(file="pdf150901_heatmap-ranked-Jasp.pdf",width=7.01,height=4.86)
draw_heatmap(jas.fc,fc,8)
dev.off()

table.fc<-table.fc[ order(-table.fc[,"rank.c.y"]), ]
vb.table <- cut.top(table.fc,20)
vb.fc <- vb.table
vb.fc[,fc] <- log(vb.table[,fc],2)

pdf(file="pdf150901_heatmap-ranked-Vb.pdf",width=7.01,height=4.86)
draw_heatmap(vb.fc,fc,8)
dev.off()

top.fc <- rbind(jas.fc,vb.fc)
top.fc <- unique(top.fc)
pdf(file="pdf150901_heatmap-dendrogram-Jas-OR-Vb.pdf",width=7.01,height=4.86)
draw_heatmap.dendro(top.fc,fc,8)
dev.off()


#Work on promoters

#add rolling circle tail to promoter sequence
prom.5p <- "CAGAGCCA"
prom.3p <- "GAAGGCTG"
prom$prom <- paste(prom.5p, prom$prom, prom.3p, sep="")
rm(prom.5p)
rm(prom.3p)

plot(log(prom$DNAreads,2)) #the curve makes me realize after prom 3500 there are sequencing errors
#on a to-find analysis I measured 3363 promoters, see email 16 03 2015
ab<-prom[1:3363,]
ab <- ab[ab$bc %in% lib$bc, ]
ab <- merge(prom[1:3363,], lib, by="bc", all=FALSE) #3199 barcodes
write.table(ab$prom_pacbio, "151104_3199-prom-sequences-with-some-duplicates.txt", sep="\t")

prom.fas <- ab$prom_pacbio
prom01 <- prom.fas[1:200]
prom02 <- prom.fas[201:400]
prom03 <- prom.fas[401:600]
prom04 <- prom.fas[601:800]
prom05 <- prom.fas[801:1000]
prom06 <- prom.fas[1001:1200]
prom07 <- prom.fas[1201:1400]
prom08 <- prom.fas[1401:1600]
prom09 <- prom.fas[1601:1800]
prom10 <- prom.fas[1801:2000]
prom11 <- prom.fas[2001:2200]
prom12 <- prom.fas[2201:2400]
prom13 <- prom.fas[2401:2600]
prom14 <- prom.fas[2601:2800]
prom15 <- prom.fas[2801:3000]
prom16 <- prom.fas[3001:3199]

prom01<-paste(prom01, collapse="") #sep="" is not what we need
prom02<-paste(prom02, collapse="")
prom03<-paste(prom03, collapse="")
prom04<-paste(prom04, collapse="")
prom05<-paste(prom05, collapse="")
prom06<-paste(prom06, collapse="")
prom07<-paste(prom07, collapse="")
prom08<-paste(prom08, collapse="")
prom09<-paste(prom09, collapse="")
prom10<-paste(prom10, collapse="")
prom11<-paste(prom11, collapse="")
prom12<-paste(prom12, collapse="")
prom13<-paste(prom13, collapse="")
prom14<-paste(prom14, collapse="")
prom15<-paste(prom15, collapse="")
prom16<-paste(prom16, collapse="")

write(prom01, file = "151105_prom01.txt", ncolumns = 1, append = FALSE, sep = "")
write(prom02, file = "151105_prom02.txt", ncolumns = 1, append = FALSE, sep = "")
write(prom03, file = "151105_prom03.txt", ncolumns = 1, append = FALSE, sep = "")
write(prom04, file = "151105_prom04.txt", ncolumns = 1, append = FALSE, sep = "")
write(prom05, file = "151105_prom05.txt", ncolumns = 1, append = FALSE, sep = "")
write(prom06, file = "151105_prom06.txt", ncolumns = 1, append = FALSE, sep = "")
write(prom07, file = "151105_prom07.txt", ncolumns = 1, append = FALSE, sep = "")
write(prom08, file = "151105_prom08.txt", ncolumns = 1, append = FALSE, sep = "")
write(prom09, file = "151105_prom09.txt", ncolumns = 1, append = FALSE, sep = "")
write(prom10, file = "151105_prom10.txt", ncolumns = 1, append = FALSE, sep = "")
write(prom11, file = "151105_prom11.txt", ncolumns = 1, append = FALSE, sep = "")
write(prom12, file = "151105_prom12.txt", ncolumns = 1, append = FALSE, sep = "")
write(prom13, file = "151105_prom13.txt", ncolumns = 1, append = FALSE, sep = "")
write(prom14, file = "151105_prom14.txt", ncolumns = 1, append = FALSE, sep = "")
write(prom15, file = "151105_prom15.txt", ncolumns = 1, append = FALSE, sep = "")
write(prom16, file = "151105_prom16.txt", ncolumns = 1, append = FALSE, sep = "")

rm(prom01,prom02,prom03,prom04,prom05,prom06,prom07,prom08,prom09,prom10,prom11,prom12,prom13,prom14,prom15,prom16)

#import motifs read with Jaspar
motif.all<-read.delim("151105_All_Jaspar_Motifs.txt", header = TRUE, as.is = TRUE, quote = "\'")
motif.agg <-aggregate(motif.all, by=list(motif.all$Model.ID), FUN=NROW)
motif.agg <- motif.agg[,1:2]
colnames(motif.agg)<-c("Model.ID", "counts")

#prepare table to associate prom expression reads to lumicycler cts

# important:
# import workspace W151105_Workspace-reads-prom-motifs.RData

basal.expr <- ab[c("ID","bc","DNAreads","prom_pacbio","D0a","V0a","D0b","J0b","D0c","J0c","V0c")]
prom.norm <- function (x) {
  x[,5:ncol(x)] <- as.numeric(as.character(unlist(x[,5:ncol(x)]))) #change bc columns from character to numeric (R faq 7.10) OK
  x[,5:ncol(x)] <- sweep(x[,5:ncol(x)], 2,sapply(x[,5:ncol(x)], median, na.rm=TRUE), `/`) #divide each bc column by its median
  x[,5:ncol(x)] <- x[,5:ncol(x)] * 500 #rescale normalized values to global median
}

basal.norm <- basal.expr
basal.norm[,5:ncol(basal.expr)] <- prom.norm(basal.expr)

#plot normalization OK
par(mfrow=c(1,2))
boxplot(basal.expr[,5:ncol(basal.expr)], outline = FALSE, main = "Before normalization", xlab = "samples", ylab = "reads") #boxplot before normalization
boxplot(basal.norm[,5:ncol(basal.norm)], outline = FALSE, main = "After median norm", xlab = "samples", ylab = "reads") #boxplot after normalization
par(mfrow=c(1,1))

write.table(basal.norm, "151110_prom-sequences-with-normalised-t0-reads.txt", sep="\t")

#export FC candidates for heatmap
write.table(top.fc, "151111_top-fold-changes-Jas-OR-Vb.txt", sep="\t")

#import top20 candidates for heatmap

top20<-read.delim("151111_top20-fold-changes-Jas-OR-Vb.tsv", header = TRUE, as.is = TRUE, quote = "\'")
				
stuff<-c("ID",  "Single.clone.tested",  "bc",	"prom_pacbio")
jas1<-c("D0b", "D1b",	"D2b",	"D8b",	"J0b",	"J1b", "J2b",	"J4b",	"J8b")
vb1<-c("D0a",  "D1a",	"D2a",	"D8a",	"V0a",	"V1a",	"V2a",	"V4a",	"V8a")
jas2<-c("D0c.x",  "D1c.x",	"D2c.x",	"D4c.x",	"D8c.x",	"J0c.x",	"J1c",	"J2c",	"J4c",	"J8c")
vb2<-c("D0c.y", "D1c.y",	"D2c.y",	"D4c.y",	"D8c.y",	"V0c.y",	"V1c",	"V2c",	"V4c",	"V8c")

column.order<-c(jas1,jas2,vb1,vb2)

draw_heatmap.dendro <- function(log.fc,group,scale) {
  cols <- c(colorRampPalette(c("cornflowerblue", "white"))(30), colorRampPalette(c("white", "red"))(30)) #color palette 30 between blue and white, 30 between white and red
  bks <- c(seq(-scale,0,length=30),seq(0,scale,length=31)) #set the breaks (breaks should be 1 more than colors)
  heatmap.2(as.matrix(log.fc[,group]), Rowv=TRUE, Colv=NULL, dendrogram="row", trace="none",
            key=FALSE, lhei=c(0.1,1), lwid=c(0.1,1), margins=c(5,12), col=cols, breaks=bks) #ignore error of figures too large  
}

pdf(file="pdf151111_heatmap-dendrogram-top20-Jas-OR-Vb.pdf",width=7.01,height=4.86)
draw_heatmap.dendro(top20,column.order,8)
dev.off()

#sposto le righe in excel in base al dendrogramma
dendro20<-read.delim("151112_dendrogram.tsv", header = TRUE, as.is = TRUE, quote = "\'")
draw_heatmap <- function(log.fc,group,scale) {
  cols <- c(colorRampPalette(c("cornflowerblue", "white"))(30), colorRampPalette(c("white", "red"))(30)) #color palette 30 between blue and white, 30 between white and red
  bks <- c(seq(-scale,0,length=30),seq(0,scale,length=31)) #set the breaks (breaks should be 1 more than colors)
  heatmap.2(as.matrix(log.fc[,group]), Rowv=NULL, Colv=NULL, dendrogram="none", trace="none",
            key=FALSE, lhei=c(0.1,1), lwid=c(0.1,1), margins=c(5,12), col=cols, breaks=bks) #ignore error of figures too large  
}
pdf(file="pdf151112_ordered-heatmap-based-on-dendrogram-top20-Jas-OR-Vb.pdf",width=7.01,height=4.86)
draw_heatmap(dendro20,column.order,8)
dev.off()





#write.table(table.fc, "150825_FC-All-libraries-ranked-by-Ueli.txt", sep="\t")

#concatenate a fasta header, two columns one is header one is dna seq
load("W151113_Workspace_R_heatmap.RData")
rm(list=setdiff(ls(), "basal.norm"))
prom.fas<-basal.norm[c("ID","bc","prom_pacbio")]
rm(basal.norm)
#add rolling circle tail to promoter sequence
prom.5p <- "CAGAGCCA"
prom.3p <- "GAAGGCTG"
prom.fas$prom_pacbio <- paste(prom.5p, prom.fas$prom_pacbio, prom.3p, sep="")
rm(prom.5p)
rm(prom.3p)

prom.fas$ID <- paste(">",prom.fas$ID, "-", prom.fas$bc, sep="")
prom.fas <- prom.fas[,c("ID","prom_pacbio")] #remove unwanted columns
prom.fas <- unique(prom.fas)
prom.fas <- do.call(rbind, lapply(seq(nrow(prom.fas)), function(i) t(prom.fas[i, ]))) #write a fasta format
write.table(prom.fas, "151131_3188-promoters.fasta", row.names = FALSE, col.names = FALSE, quote = FALSE)
#note: FIMO found 6 duplicated ID in lines 401, 2863, 2865, 3539, 6291, 6293 that were removed by hand



#print distribution of known barcodes
load("W151113_Workspace_R_heatmap.RData")
rm(list=setdiff(ls(), "lib"))
a<-rowMeans(lib[2:34])
a<-sort(a, decreasing=TRUE)
plot(log(a,10))



#Draw complete heatmaps for supp figure
rm(list=ls())
setwd('/Users/randogp/Google Drive/Schibler Lab/STARpromBC/R-scripts_STARprom') #MacBook Air
setwd('/Users/admin/Google Drive/Schibler Lab/STARpromBC/R-scripts_STARprom') #iMac21

tax_lat <- read.delim("150419_Cytoskeletal-drugs-reads_AND_FC.txt", header = TRUE, as.is = TRUE)

tax_lat <- tax_lat[c("bc",
                     "l.rank","l6.fc","l7.fc","l8.fc","l9.fc", "l0.fc","l1.fc","l2.fc","l3.fc","l4.fc","l5.fc",
                     "t.rank","t6.fc","t7.fc","t8.fc","t9.fc","t0.fc","t1.fc","t2.fc","t3.fc","t4.fc", "t5.fc"
)]
colnames(tax_lat)<-c("bc",
                     "l1.rank","l1.dmso.0h","l1.dmso.1h","l1.dmso.2h","l1.dmso.8h", "l1.lat.0h","l1.lat.05h","l1.lat.1h","l1.lat.2h","l1.lat.4h","l1.lat.8h",
                     "t1.rank","t1.dmso.0h","t1.dmso.1h","t1.dmso.2h","t1.dmso.8h", "t1.tax.0h","t1.tax.05h","t1.tax.1h","t1.tax.2h","t1.tax.4h","t1.tax.8h"
)

jas_vb <- read.delim("150825_FC-All-libraries-ranked-by-Ueli.txt", header = TRUE, as.is = TRUE)
colnames(jas_vb)<-c("bc",
                    "j1.rank","v1.rank","j2.rank","v2.rank",
                    "j1.dmso.0h","j1.dmso.1h","j1.dmso.2h","j1.dmso.8h", "j1.jas.0h","j1.jas.1h","j1.jas.2h","j1.jas.4h","j1.jas.8h",
                    "v1.dmso.0h","v1.dmso.1h","v1.dmso.2h","v1.dmso.8h", "v1.vb.0h","v1.vb.1h","v1.vb.2h","v1.vb.4h","v1.vb.8h",
                    "j2.dmso.0h","j2.dmso.1h","j2.dmso.2h","j2.dmso.4h", "j2.dmso.8h", "j2.jas.0h","j2.jas.1h","j2.jas.2h","j2.jas.4h","j2.jas.8h",
                    "v2.dmso.0h","v2.dmso.1h","v2.dmso.2h","v2.dmso.4h", "v2.dmso.8h", "v2.vb.0h","v2.vb.1h","v2.vb.2h","v2.vb.4h","v2.vb.8h"
)

jas_vb <- jas_vb[c("bc",
                   "j1.rank","j1.dmso.0h","j1.dmso.1h","j1.dmso.2h","j1.dmso.8h", "j1.jas.0h","j1.jas.1h","j1.jas.2h","j1.jas.4h","j1.jas.8h",
                   "v1.rank","v1.dmso.0h","v1.dmso.1h","v1.dmso.2h","v1.dmso.8h", "v1.vb.0h","v1.vb.1h","v1.vb.2h","v1.vb.4h","v1.vb.8h",
                   "j2.rank","j2.dmso.0h","j2.dmso.1h","j2.dmso.2h","j2.dmso.4h", "j2.dmso.8h", "j2.jas.0h","j2.jas.1h","j2.jas.2h","j2.jas.4h","j2.jas.8h",
                   "v2.rank","v2.dmso.0h","v2.dmso.1h","v2.dmso.2h","v2.dmso.4h", "v2.dmso.8h", "v2.vb.0h","v2.vb.1h","v2.vb.2h","v2.vb.4h","v2.vb.8h"
)]

all.drugs<-merge(jas_vb, tax_lat, by="bc", all=FALSE)
rm(jas_vb)
rm(tax_lat)

jasplakinolide1<-all.drugs[c("bc","j1.rank","j1.dmso.0h","j1.dmso.1h","j1.dmso.2h","j1.dmso.8h", "j1.jas.0h","j1.jas.1h","j1.jas.2h","j1.jas.4h","j1.jas.8h")]
jasplakinolide1<-jasplakinolide1[ order(-jasplakinolide1[,"j1.rank"]), ]

jasplakinolide2<-all.drugs[c("bc","j2.rank","j2.dmso.0h","j2.dmso.1h","j2.dmso.2h","j2.dmso.4h", "j2.dmso.8h", "j2.jas.0h","j2.jas.1h","j2.jas.2h","j2.jas.4h","j2.jas.8h")]
jasplakinolide2<-jasplakinolide2[ order(-jasplakinolide2[,"j2.rank"]), ]

vinblastine1<-all.drugs[c("bc","v1.rank","v1.dmso.0h","v1.dmso.1h","v1.dmso.2h","v1.dmso.8h", "v1.vb.0h","v1.vb.1h","v1.vb.2h","v1.vb.4h","v1.vb.8h")]
vinblastine1<-vinblastine1[ order(-vinblastine1[,"v1.rank"]), ]

vinblastine2<-all.drugs[c("bc","v2.rank","v2.dmso.0h","v2.dmso.1h","v2.dmso.2h","v2.dmso.4h", "v2.dmso.8h", "v2.vb.0h","v2.vb.1h","v2.vb.2h","v2.vb.4h","v2.vb.8h")]
vinblastine2<-vinblastine2[ order(-vinblastine2[,"v2.rank"]), ]

taxol1<-all.drugs[c("bc","t1.rank","t1.dmso.0h","t1.dmso.1h","t1.dmso.2h","t1.dmso.8h", "t1.tax.0h","t1.tax.05h","t1.tax.1h","t1.tax.2h","t1.tax.4h","t1.tax.8h")]
taxol1<-taxol1[ order(-taxol1[,"t1.rank"]), ]

latrunculinB1<-all.drugs[c("bc","l1.rank","l1.dmso.0h","l1.dmso.1h","l1.dmso.2h","l1.dmso.8h", "l1.lat.0h","l1.lat.05h","l1.lat.1h","l1.lat.2h","l1.lat.4h","l1.lat.8h")]
latrunculinB1<-latrunculinB1[ order(-latrunculinB1[,"l1.rank"]), ]

j1.fc <- c("j1.dmso.0h","j1.dmso.1h","j1.dmso.2h","j1.dmso.8h", "j1.jas.0h","j1.jas.1h","j1.jas.2h","j1.jas.4h","j1.jas.8h")
j2.fc <- c("j2.dmso.0h","j2.dmso.1h","j2.dmso.2h","j2.dmso.4h", "j2.dmso.8h", "j2.jas.0h","j2.jas.1h","j2.jas.2h","j2.jas.4h","j2.jas.8h")
v1.fc <- c("v1.dmso.0h","v1.dmso.1h","v1.dmso.2h","v1.dmso.8h", "v1.vb.0h","v1.vb.1h","v1.vb.2h","v1.vb.4h","v1.vb.8h")
v2.fc <- c("v2.dmso.0h","v2.dmso.1h","v2.dmso.2h","v2.dmso.4h", "v2.dmso.8h", "v2.vb.0h","v2.vb.1h","v2.vb.2h","v2.vb.4h","v2.vb.8h")
t1.fc <- c("t1.dmso.0h","t1.dmso.1h","t1.dmso.2h","t1.dmso.8h", "t1.tax.0h","t1.tax.05h","t1.tax.1h","t1.tax.2h","t1.tax.4h","t1.tax.8h")
l1.fc <- c("l1.dmso.0h","l1.dmso.1h","l1.dmso.2h","l1.dmso.8h", "l1.lat.0h","l1.lat.05h","l1.lat.1h","l1.lat.2h","l1.lat.4h","l1.lat.8h")

log.j1.fc <- jasplakinolide1
log.j2.fc <- jasplakinolide2
log.v1.fc <- vinblastine1
log.v2.fc <- vinblastine2
log.t1.fc <- taxol1
log.l1.fc <- latrunculinB1

write.table(log.j1.fc, "160514_jasplakinolide1_FC.tab", row.names = FALSE, col.names = TRUE, quote = FALSE)
write.table(log.j2.fc, "160514_jasplakinolide2_FC.tab", row.names = FALSE, col.names = TRUE, quote = FALSE)
write.table(log.v1.fc, "160514_vinblastine1_FC.tab", row.names = FALSE, col.names = TRUE, quote = FALSE)
write.table(log.v2.fc, "160514_vinblastine2_FC.tab", row.names = FALSE, col.names = TRUE, quote = FALSE)
write.table(log.t1.fc, "160514_paclitaxel1_FC.tab", row.names = FALSE, col.names = TRUE, quote = FALSE)
write.table(log.l1.fc, "160514_latrunculin1_FC.tab", row.names = FALSE, col.names = TRUE, quote = FALSE)

log.j1.fc[,j1.fc] <- log( log.j1.fc[,j1.fc],2) 
log.j2.fc[,j2.fc] <- log( log.j2.fc[,j2.fc],2) 
log.v1.fc[,v1.fc] <- log( log.v1.fc[,v1.fc],2) 
log.v2.fc[,v2.fc] <- log( log.v2.fc[,v2.fc],2) 
log.t1.fc[,t1.fc] <- log( log.t1.fc[,t1.fc],2) 
log.l1.fc[,l1.fc] <- log( log.l1.fc[,l1.fc],2) 

library(gplots)

draw_heatmap <- function(log.fc,group,scale) {
  cols <- c(colorRampPalette(c("cornflowerblue", "white"))(30), colorRampPalette(c("white", "red"))(30)) #color palette 30 between blue and white, 30 between white and red
  bks <- c(seq(-scale,0,length=30),seq(0,scale,length=31)) #set the breaks (breaks should be 1 more than colors)
  heatmap.2(as.matrix(log.fc[,group]), Rowv=NULL, Colv=NULL, labRow = NULL, dendrogram="none", trace="none",
            key=FALSE, lhei=c(0.1,1), lwid=c(0.1,1), margins=c(5,12), col=cols, breaks=bks) #ignore error of figures too large  
}


pdf(file="pdf151203_heatmap-j1_scale4.pdf",width=7.01,height=4.86)
draw_heatmap(log.j1.fc,j1.fc,4)
dev.off()
pdf(file="pdf151203_heatmap-j2_scale4.pdf",width=7.01,height=4.86)
draw_heatmap(log.j2.fc,j2.fc,4)
dev.off()
pdf(file="pdf151203_heatmap-v1_scale4.pdf",width=7.01,height=4.86)
draw_heatmap(log.v1.fc,v1.fc,4)
dev.off()
pdf(file="pdf151203_heatmap-v2_scale4.pdf",width=7.01,height=4.86)
draw_heatmap(log.v2.fc,v2.fc,4)
dev.off()
pdf(file="pdf151203_heatmap-l1_scale4.pdf",width=7.01,height=4.86)
draw_heatmap(log.l1.fc,l1.fc,4)
dev.off()
pdf(file="pdf151203_heatmap-t1_scale4.pdf",width=7.01,height=4.86)
draw_heatmap(log.t1.fc,t1.fc,4)
dev.off()

#Export log fold changes for Promoter Set Enrichment analysis

load("W151113_Workspace_R_heatmap.RData")
basal.norm<-basal.norm[c("ID","bc","prom_pacbio")]
rm(list=setdiff(ls(), c("basal.norm","log.j1.fc","log.j2.fc","log.v1.fc","log.v2.fc","log.l1.fc","log.t1.fc")))
#save.image("W151210_worskpace_input_PSEA.RData")



#Promoter Set Enrichment Analysis, build the gct tables
load("W151210_worskpace_input_PSEA.RData")

j1<-merge(basal.norm, log.j1.fc, by="bc", all=FALSE)
j1 <- j1[order(-j1$j1.rank),]
j1$ID <- paste("BC",j1$ID, sep="")            
j1<-j1[,c("ID","ID",          
          "j1.dmso.0h", "j1.dmso.1h",	"j1.dmso.2h",	"j1.dmso.8h",
          "j1.jas.1h",  "j1.jas.2h",	"j1.jas.4h",	"j1.jas.8h"
)]
j1 <- unique(j1)
write.table(j1, "151210_3118BC-j1_gct-for-GSEA.tab", row.names = FALSE, col.names = TRUE, quote = FALSE) #I will add the required 2 rows by hand on excel

j2<-merge(basal.norm, log.j2.fc, by="bc", all=FALSE)
j2 <- j2[order(-j2$j2.rank),]
j2$ID <- paste("BC",j2$ID, sep="")            
j2<-j2[,c("ID","ID",
          "j2.dmso.1h",  "j2.dmso.2h",	"j2.dmso.4h",	"j2.dmso.8h",
          "j2.jas.1h",  "j2.jas.2h",	"j2.jas.4h",	"j2.jas.8h"
)]
j2 <- unique(j2)
write.table(j2, "151210_3118BC-j2_gct-for-GSEA.tab", row.names = FALSE, col.names = TRUE, quote = FALSE)

v1<-merge(basal.norm, log.v1.fc, by="bc", all=FALSE)
v1 <- v1[order(-v1$v1.rank),]
v1$ID <- paste("BC",v1$ID, sep="")            
v1<-v1[,c("ID","ID",
          "v1.dmso.0h", "v1.dmso.1h",	"v1.dmso.2h",	"v1.dmso.8h",
          "v1.vb.1h", "v1.vb.2h",	"v1.vb.4h",	"v1.vb.8h"
)]
v1 <- unique(v1)
write.table(v1, "151210_3118BC-v1_gct-for-GSEA.tab", row.names = FALSE, col.names = TRUE, quote = FALSE)

v2<-merge(basal.norm, log.v2.fc, by="bc", all=FALSE)
v2 <- v2[order(-v2$v2.rank),]
v2$ID <- paste("BC",v2$ID, sep="")            
v2<-v2[,c("ID","ID",
          "v2.dmso.1h", "v2.dmso.2h",	"v2.dmso.4h",	"v2.dmso.8h",
          "v2.vb.1h",  "v2.vb.2h",	"v2.vb.4h",	"v2.vb.8h"
)]
v2 <- unique(v2)
write.table(v2, "151210_3118BC-v2_gct-for-GSEA.tab", row.names = FALSE, col.names = TRUE, quote = FALSE)


t1<-merge(basal.norm, log.t1.fc, by="bc", all=FALSE)
t1 <- t1[order(-t1$t1.rank),]
t1$ID <- paste("BC",t1$ID, sep="")            
t1<-t1[,c("ID","ID",      
          "t1.tax.1h",  "t1.tax.2h",	"t1.tax.4h",	"t1.tax.8h",
          "t1.dmso.0h",  "t1.dmso.1h",	"t1.dmso.2h",	"t1.dmso.8h"
)]
t1 <- unique(t1)
write.table(t1, "151210_3118BC-t1_gct-for-GSEA.tab", row.names = FALSE, col.names = TRUE, quote = FALSE)

l1<-merge(basal.norm, log.l1.fc, by="bc", all=FALSE)
l1 <- l1[order(-l1$l1.rank),]
l1$ID <- paste("BC",l1$ID, sep="")            
l1<-l1[,c("ID","ID",      
          "l1.dmso.0h", "l1.dmso.1h",	"l1.dmso.2h",	"l1.dmso.8h",
          "l1.lat.1h",  "l1.lat.2h",	"l1.lat.4h",	"l1.lat.8h"
)]
l1 <- unique(l1)
write.table(l1, "151210_3118BC-l1_gct-for-GSEA.tab", row.names = FALSE, col.names = TRUE, quote = FALSE)

rm(list=setdiff(ls(), c("j1","j2","v1","v2","t1","l1")))
#save.image("W151210_worskpace_output_PSEA-prep.RData")


#old stuff deprecated, merge jasp and vinb duplicates to have a stronger PSEA table
# jasp<-merge(log.j1.fc, log.j2.fc, by="bc", all=FALSE)
# jasp <- cbind(rank = 0, jasp) #add an empty column named rank
# jasp$rank <- rowMeans(jasp[c("j1.rank","j2.rank")])
# jasp <- subset(jasp, select = -c(j1.rank,j2.rank) )
# jasp <- jasp[order(-jasp$rank),]
# 
# load("W151113_Workspace_R_heatmap.RData")
# rm(list=setdiff(ls(), c("basal.norm","jasp")))
# basal.norm<-basal.norm[c("ID","bc","prom_pacbio")]
# 
# jasp<-merge(basal.norm, jasp, by="bc", all=FALSE)
# rm(basal.norm)
# jasp <- jasp[order(-jasp$rank),]
# jasp$ID <- paste("BC",jasp$ID, sep="")            
# jasp<-jasp[,c("ID",
#               "j1.dmso.0h","j1.dmso.2h","j1.dmso.8h","j1.jas.2h","j1.jas.4h","j1.jas.8h",
#               "j2.dmso.0h","j2.dmso.2h","j2.dmso.4h","j2.jas.2h","j2.jas.4h","j2.jas.8h"
#               )]
# jasp <- unique(jasp)
# write.table(jasp, "151203_3118BC-Jasplakinolide_gct-for-GSEA.tab", row.names = FALSE, col.names = FALSE, quote = FALSE)
# 
# 
# 
# #build the vinb gct table for motif enrichment analysis
# vinb<-merge(log.v1.fc, log.v2.fc, by="bc", all=FALSE)
# vinb <- cbind(rank = 0, vinb) #add an empty column named BC
# vinb$rank <- rowMeans(vinb[c("v1.rank","v2.rank")])
# vinb <- subset(vinb, select = -c(v1.rank,v2.rank) )
# vinb <- vinb[order(-vinb$rank),]
# 
# load("W151113_Workspace_R_heatmap.RData")
# rm(list=setdiff(ls(), c("basal.norm","vinb")))
# basal.norm<-basal.norm[c("ID","bc","prom_pacbio")]
# 
# vinb<-merge(basal.norm, vinb, by="bc", all=FALSE)
# rm(basal.norm)
# vinb <- vinb[order(-vinb$rank),]
# vinb$ID <- paste("BC",vinb$ID, sep="")            
# vinb<-vinb[,c("ID",
#               "v1.dmso.0h","v1.dmso.2h","v1.dmso.8h","v1.vb.2h","v1.vb.4h","v1.vb.8h",
#               "v2.dmso.0h","v2.dmso.2h","v2.dmso.4h","v2.vb.2h","v2.vb.4h","v2.vb.8h"
# )]
# vinb <- unique(vinb)
# write.table(vinb, "151207_3118BC-vinblastine_gct-for-GSEA.tab", row.names = FALSE, col.names = FALSE, quote = FALSE)





#draw heatmap for supp1, including top 100.


load("W151210_worskpace_input_PSEA.RData")

j1.fc <- c("j1.dmso.0h","j1.dmso.1h","j1.dmso.2h","j1.dmso.8h", "j1.jas.0h","j1.jas.1h","j1.jas.2h","j1.jas.4h","j1.jas.8h")
j2.fc <- c("j2.dmso.0h","j2.dmso.1h","j2.dmso.2h","j2.dmso.4h", "j2.dmso.8h", "j2.jas.0h","j2.jas.1h","j2.jas.2h","j2.jas.4h","j2.jas.8h")
v1.fc <- c("v1.dmso.0h","v1.dmso.1h","v1.dmso.2h","v1.dmso.8h", "v1.vb.0h","v1.vb.1h","v1.vb.2h","v1.vb.4h","v1.vb.8h")
v2.fc <- c("v2.dmso.0h","v2.dmso.1h","v2.dmso.2h","v2.dmso.4h", "v2.dmso.8h", "v2.vb.0h","v2.vb.1h","v2.vb.2h","v2.vb.4h","v2.vb.8h")
t1.fc <- c("t1.dmso.0h","t1.dmso.1h","t1.dmso.2h","t1.dmso.8h", "t1.tax.0h","t1.tax.05h","t1.tax.1h","t1.tax.2h","t1.tax.4h","t1.tax.8h")
l1.fc <- c("l1.dmso.0h","l1.dmso.1h","l1.dmso.2h","l1.dmso.8h", "l1.lat.0h","l1.lat.05h","l1.lat.1h","l1.lat.2h","l1.lat.4h","l1.lat.8h")

top100.j1.fc<-log.j1.fc[1:100,]
top100.j2.fc<-log.j2.fc[1:100,]
top100.v1.fc<-log.v1.fc[1:100,]
top100.v2.fc<-log.v2.fc[1:100,]
top100.l1.fc<-log.l1.fc[1:100,]
top100.t1.fc<-log.t1.fc[1:100,]

library(gplots)
draw_heatmap <- function(log.fc,group,scale) {
  cols <- c(colorRampPalette(c("cornflowerblue", "white"))(30), colorRampPalette(c("white", "red"))(30)) #color palette 30 between blue and white, 30 between white and red
  bks <- c(seq(-scale,0,length=30),seq(0,scale,length=31)) #set the breaks (breaks should be 1 more than colors)
  heatmap.2(as.matrix(log.fc[,group]), Rowv=NULL, Colv=NULL, labRow = NULL, dendrogram="none", trace="none",
            key=FALSE, lhei=c(0.1,1), lwid=c(0.1,1), margins=c(5,12), col=cols, breaks=bks) #ignore error of figures too large  
}

pdf(file="pdf151203_heatmap-j1_scale4.pdf",width=7.01,height=4.86)
draw_heatmap(log.j1.fc,j1.fc,4)
dev.off()

draw_heatmap(top100.j1.fc,j1.fc,4)
draw_heatmap(top100.j2.fc,j2.fc,4)
draw_heatmap(top100.v1.fc,v1.fc,4)
draw_heatmap(top100.v2.fc,v2.fc,4)
draw_heatmap(top100.l1.fc,l1.fc,4)
draw_heatmap(top100.t1.fc,t1.fc,4)

draw_heatmap(log.j1.fc,j1.fc,4)
draw_heatmap(log.j2.fc,j2.fc,4)
draw_heatmap(log.v1.fc,v1.fc,4)
draw_heatmap(log.v2.fc,v2.fc,4)
draw_heatmap(log.l1.fc,l1.fc,4)
draw_heatmap(log.t1.fc,t1.fc,4)
#save.image("~/Google Drive/Schibler Lab/STARpromBC/R-scripts_STARprom/W151217_Heatmaps.RData")


#count the median motif for each
rm(list=ls())
fimo<-read.delim("151131_fimo_Transfac_2012_motifs_mapped_on_STAR-PROM.txt", header = TRUE, as.is = TRUE, quote = "\'")
fimo <- cbind(motif.count = 1, fimo)
aggdata <-aggregate(fimo$motif.count, by=list(fimo$sequence.name), FUN=sum, na.rm=TRUE)
summary(aggdata) #median 29 motifs, min 2 max 204


