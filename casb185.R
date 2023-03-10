library(GEOquery)
library(limma)

# load data from GEO database
gset <- getGEO("GSE42890", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL8759", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group membership for all samples
gsms <- paste0("1101101111111211121XX210X1100110112100101110221110",
               "12211220000111111102202212022X221111100011111110X2",
               "0011110111010022222002110022X120001011110010211002",
               "100XX000110111122011111210011001022")
sml <- strsplit(gsms, split="")[[1]]

# filter out excluded samples (marked as "X")
sel <- which(sml != "X")
sml <- sml[sel]
gset <- gset[ ,sel]

# log2 transformation
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }

# assign samples to groups and set up design matrix
gs <- factor(sml)
groups <- make.names(c("0-10","10-20","20-30"))
levels(gs) <- groups
gset$group <- gs
design <- model.matrix(~group + 0, gset)
colnames(design) <- levels(gs)

# fit linear model to be used in ebayes 
fit <- lmFit(gset, design)

# set up contrasts of interest and recalculate model coefficients
cts <- paste(groups, c(tail(groups, -1), head(groups, 1)), sep="-")
cont.matrix <- makeContrasts(contrasts=cts, levels=design)
model <- contrasts.fit(fit, cont.matrix)

# compute t statistics and add the 250 most significant genes to a table
model <- eBayes(model, 0.01)
tT <- topTable(model, adjust="fdr", sort.by="B", number=250)

tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","F","Gene.symbol","Gene.title"))
write.table(tT, file=stdout(), row.names=F, sep="\t")

# histogram of observed p values
tT2 <- topTable(model, adjust="fdr", sort.by="B", number=Inf)
hist(tT2$adj.P.Val, col = "grey", border = "white", xlab = "P-adj",
     ylab = "Number of genes", main = "P-adj value distribution")

# summarize test results as "up", "down" or "not expressed"
dT <- decideTests(model, adjust.method="fdr", p.value=0.05)

# volcano plot
colnames(model)
ct <- 1 
volcanoplot(model, coef=ct, main=colnames(model)[ct], pch=20,
            highlight=length(which(dT[,ct]!=0)), names=rep('+', nrow(model)))

# Mean difference plot
plotMD(model, column=ct, status=dT[,ct], legend=F, pch=20, cex=1)
abline(h=0)

