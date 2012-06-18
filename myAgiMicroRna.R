# LGMB   Laboratory of Genetics and Molecular Biology
# BiT    Bioinformatics Team
#
# Project...: myArray
# Date......: 22/07/2011
# Author(s).: daniel, matheus
###############################################################################



my.qcPlots <- function (dd, offset = 5, MeanSignal = TRUE, ProcessedSignal = FALSE, 
		TotalProbeSignal = FALSE, TotalGeneSignal = FALSE, BGMedianSignal = FALSE, 
		BGUsed = FALSE, targets) 
{
	if (!is(dd, "uRNAList")) {
		stop("'input' must be a uRNAList")
		if (is.null(dim(dd)[1])) {
			stop("'input' is empty")
		}
	}
	if (MeanSignal) {
		MMM = dd$meanS
		min = min(MMM)
		for (i in 1:dim(MMM)[2]) {
			MMM[, i] = MMM[, i] + (abs(min) + offset)
		}
		MMM = log2(MMM)
		maintitle = "MeanSignal"
		colorfill = "orange"
		dev.set(dev.next())
		boxplotMicroRna(MMM, maintitle, colorfill)
		dev.set(dev.next())
		plotDensityMicroRna(MMM, maintitle)
		dev.set(dev.next())
		ddaux = dd
		ddaux$meanS = MMM
		mvaMicroRna(ddaux, maintitle, verbose = FALSE)
		rm(ddaux)
		maintitle = "MeanSignal- RLE "
		dev.set(dev.next())
		RleMicroRna(MMM, maintitle, colorfill)
		dev.set(dev.next())
		hierclusMicroRna(MMM, targets$GErep, methdis = "euclidean", 
				methclu = "complete", sel = FALSE, 100)
	}
	if (ProcessedSignal) {
		MMM = dd$procS
		min = min(MMM)
		for (i in 1:dim(MMM)[2]) {
			MMM[, i] = MMM[, i] + (abs(min) + offset)
		}
		MMM = log2(MMM)
		maintitle = "ProcessedSignal"
		colorfill = "blue"
		dev.set(dev.next())
		boxplotMicroRna(MMM, maintitle, colorfill)
		dev.set(dev.next())
		plotDensityMicroRna(MMM, maintitle)
		dev.set(dev.next())
		ddaux = dd
		ddaux$TPS = MMM
		mvaMicroRna(ddaux, maintitle, verbose = FALSE)
		rm(ddaux)
		maintitle = "ProcessedSignal - RLE "
		dev.set(dev.next())
		RleMicroRna(MMM, maintitle, colorfill)
	}
	if (TotalProbeSignal) {
		up = which(duplicated(dd$genes$ProbeName) == FALSE)
		ddaux = dd[up, ]
		MMM = ddaux$TPS
		min = min(MMM)
		for (i in 1:dim(MMM)[2]) {
			MMM[, i] = MMM[, i] + (abs(min) + offset)
		}
		MMM = log2(MMM)
		maintitle = "TotalProbeSignal"
		colorfill = "red"
		dev.set(dev.next())
		boxplotMicroRna(MMM, maintitle, colorfill)
		dev.set(dev.next())
		plotDensityMicroRna(MMM, maintitle)
		maintitle = " TotalProbeSignal - RLE "
		dev.set(dev.next())
		RleMicroRna(MMM, maintitle, colorfill)
	}
	if (TotalGeneSignal) {
		ddTGS = tgsMicroRna(dd, offset, half = FALSE, makePLOT = FALSE, 
				verbose = FALSE)
		MMM = log2(ddTGS$TGS)
		maintitle = "TotalGeneSignal"
		colorfill = "green"
		dev.set(dev.next())
		boxplotMicroRna(MMM, maintitle, colorfill)
		dev.set(dev.next())
		plotDensityMicroRna(MMM, maintitle)
		maintitle = " TotalGeneSignal - RLE "
		dev.set(dev.next())
		RleMicroRna(MMM, maintitle, colorfill)
	}
	if (BGMedianSignal) {
		MMM = log2(dd$other$gBGMedianSignal)
		maintitle = "BGMedianSignal"
		colorfill = "yellow"
		dev.set(dev.next())
		boxplotMicroRna(MMM, maintitle, colorfill)
	}
	if (BGUsed) {
		MMM = log2(dd$other$gBGUsed)
		maintitle = "BGused"
		colorfill = "cyan"
		dev.set(dev.next())
		boxplotMicroRna(MMM, maintitle, colorfill)
	}
}
#
#
#
#qualityPlots <- function (dd, offset = 5, MeanSignal = TRUE, ProcessedSignal = FALSE, 
#		TotalProbeSignal = FALSE, TotalGeneSignal = FALSE, BGMedianSignal = FALSE, 
#		BGUsed = FALSE, targets) 
#{
#	if (!is(dd, "uRNAList")) {
#		stop("'input' must be a uRNAList")
#		if (is.null(dim(dd)[1])) {
#			stop("'input' is empty")
#		}
#	}
#	if (MeanSignal) {
#		MMM = dd$meanS
#		min = min(MMM)
#		for (i in 1:dim(MMM)[2]) {
#			MMM[, i] = MMM[, i] + (abs(min) + offset)
#		}
#		MMM = log2(MMM)
#		maintitle = "MeanSignal"
#		colorfill = "orange"
#		dev.set(dev.next())
#		boxplotMicroRna(MMM, maintitle, colorfill)
#		dev.set(dev.next())
#		plotDensityMicroRna(MMM, maintitle)
#		dev.set(dev.next())
#		ddaux = dd
#		ddaux$meanS = MMM
#		mvaMicroRna(ddaux, maintitle, verbose = FALSE)
#		rm(ddaux)
#		maintitle = "MeanSignal- RLE "
#		dev.set(dev.next())
#		RleMicroRna(MMM, maintitle, colorfill)
#		dev.set(dev.next())
#		hierclusMicroRna(MMM, targets$GErep, methdis = "euclidean", 
#				methclu = "complete", sel = FALSE, 100)
#	}
#	if (ProcessedSignal) {
#		MMM = dd$procS
#		min = min(MMM)
#		for (i in 1:dim(MMM)[2]) {
#			MMM[, i] = MMM[, i] + (abs(min) + offset)
#		}
#		MMM = log2(MMM)
#		maintitle = "ProcessedSignal"
#		colorfill = "blue"
#		dev.set(dev.next())
#		boxplotMicroRna(MMM, maintitle, colorfill)
#		dev.set(dev.next())
#		plotDensityMicroRna(MMM, maintitle)
#		dev.set(dev.next())
#		ddaux = dd
#		ddaux$TPS = MMM
#		mvaMicroRna(ddaux, maintitle, verbose = FALSE)
#		rm(ddaux)
#		maintitle = "ProcessedSignal - RLE "
#		dev.set(dev.next())
#		RleMicroRna(MMM, maintitle, colorfill)
#	}
#	if (TotalProbeSignal) {
#		up = which(duplicated(dd$genes$ProbeName) == FALSE)
#		ddaux = dd[up, ]
#		MMM = ddaux$TPS
#		min = min(MMM)
#		for (i in 1:dim(MMM)[2]) {
#			MMM[, i] = MMM[, i] + (abs(min) + offset)
#		}
#		MMM = log2(MMM)
#		maintitle = "TotalProbeSignal"
#		colorfill = "red"
#		dev.set(dev.next())
#		boxplotMicroRna(MMM, maintitle, colorfill)
#		dev.set(dev.next())
#		plotDensityMicroRna(MMM, maintitle)
#		maintitle = " TotalProbeSignal - RLE "
#		dev.set(dev.next())
#		RleMicroRna(MMM, maintitle, colorfill)
#	}
#	if (TotalGeneSignal) {
#		ddTGS = tgsMicroRna(dd, offset, half = FALSE, makePLOT = FALSE, 
#				verbose = FALSE)
#		MMM = log2(ddTGS$TGS)
#		maintitle = "TotalGeneSignal"
#		colorfill = "green"
#		dev.set(dev.next())
#		boxplotMicroRna(MMM, maintitle, colorfill)
#		dev.set(dev.next())
#		plotDensityMicroRna(MMM, maintitle)
#		maintitle = " TotalGeneSignal - RLE "
#		dev.set(dev.next())
#		RleMicroRna(MMM, maintitle, colorfill)
#	}
#	if (BGMedianSignal) {
#		MMM = log2(dd$other$gBGMedianSignal)
#		maintitle = "BGMedianSignal"
#		colorfill = "yellow"
#		dev.set(dev.next())
#		boxplotMicroRna(MMM, maintitle, colorfill)
#	}
#	if (BGUsed) {
#		MMM = log2(dd$other$gBGUsed)
#		maintitle = "BGused"
#		colorfill = "cyan"
#		dev.set(dev.next())
#		boxplotMicroRna(MMM, maintitle, colorfill)
#	}
#}

# old getTGS()
my.tgsMicroRna <- function (dd, offset = 0, half = TRUE, makePLOT = FALSE, verbose = FALSE) 
{
	if (!is(dd, "uRNAList")) {
		stop("'input' must be a uRNAList")
		if (is.null(dim(dd)[1])) {
			stop("'input' is empty")
		}
	}
	TT = dim(dd)[1]
	uniqueProbe = unique(dd$genes$ProbeName)
	LUP = length(uniqueProbe)
	uniqueGene = unique(dd$genes$GeneName)
	LUG = length(uniqueGene)
	if (verbose) {
		cat("\n")
		cat("GETTING Agilent Feature Extraction TotalGeneSignal", 
				"\n")
		cat("\n")
		cat("\tTotal Probes:\t", TT, "\n")
		cat("\tUnique Probe: \t", LUP, "\n")
		cat("\tUnique Gene: \t", LUG, "\n")
		cat("\n")
	}
	ug = which(duplicated(dd$genes$GeneName) == FALSE)
	nGEN = length(ug)
	ddTGS = dd[ug, ]
	if (half) {
		for (i in 1:dim(ddTGS)[2]) {
			index = which(ddTGS$TGS[, i] < 0.5)
			ddTGS$TGS[index, i] = 0.5
		}
	}
	else {
		min = min(ddTGS$TGS)
		for (i in 1:dim(ddTGS)[2]) {
			ddTGS$TGS[, i] = ddTGS$TGS[, i] + (abs(min) + offset)
		}
	}
	ddTGS$TGS = ddTGS$TGS
	ddTGS$TGS = ddTGS$TGS
	ddTGS$TGS = ddTGS$TGS
	ddTGS$TGS = ddTGS$TGS
	nARR = dim(ddTGS)[2]
	geneNames = list(c(dd$genes$GeneName[ug]), c(1:nARR))
	if (!missing(makePLOT)) {
		if (makePLOT) {
			MMM = log2(ddTGS$TGS)
			maintitle = "TotalGeneSignal"
			colorfill = "green"
			dev.set(dev.next())
			boxplotMicroRna(MMM, maintitle, colorfill)
			dev.set(dev.next())
			plotDensityMicroRna(MMM, maintitle)
			dev.set(dev.next())
			ddaux = ddTGS
			ddaux$meanS = MMM
			mvaMicroRna(ddaux, maintitle, TRUE)
			rm(ddaux)
		}
	}
	return(ddTGS)
}


my.tgsNormalization <- function (ddTGS, NORMmethod = "quantile", makePLOTpre = FALSE, 
		makePLOTpost = FALSE, targets, verbose = FALSE) 
{
	if (!is(ddTGS, "uRNAList")) {
		stop("'input' must be a uRNAList")
		if (is.null(dim(ddTGS)[1])) {
			stop("'input' is empty")
		}
	}
	if (NORMmethod != "none" && NORMmethod != "quantile" && NORMmethod != 
			"scale" && NORMmethod != "Tquantile" && NORMmethod != "scale75thp" && NORMmethod != "Tscale75thp") {
		stop("NORMmethod should be one of 'none', 'quantile','scale', 'Tquantile', 'scale75thp' or 'Tscale75thp' ")
	}
	if (!missing(makePLOTpre)) {
		if (makePLOTpre) {
			MMM = log2(ddTGS$TGS)
			colorfill = "blue"
			maintitle = "NOT NORM."
			dev.set(dev.next())
			plotDensityMicroRna(MMM, maintitle)
			dev.set(dev.next())
			boxplotMicroRna(MMM, maintitle, colorfill)
			dev.set(dev.next())
			mvaBASIC(MMM, colorfill, maintitle)
			dev.set(dev.next())
			hierclusMicroRna(MMM, targets$GErep, methdis = "euclidean", 
					methclu = "complete", sel = FALSE, 100)
			maintitle = "NOT NORM - RLE "
			dev.set(dev.next())
			RleMicroRna(MMM, maintitle, colorfill)
			rm(MMM)
		}
	}
	if (NORMmethod != "none") {
		
		source("/work/projects/manalysis/myArray/myAgi4x44PreProcess.R")
		
		ddNORM = ddTGS
		
		if (NORMmethod == "quantile") {
			exprsNORM = normalizeBetweenArrays(ddTGS$TGS, method = NORMmethod)
			ddNORM$TGS = log2(exprsNORM)
		}
		else if (NORMmethod == "scale75thp") {
			exprsNORM = my.normalizeBetweenArrays(ddTGS$TGS, method = NORMmethod)
			ddNORM$TGS = exprsNORM
		}
		else if (NORMmethod == "Tquantile") {
			
			if (!is.null(ddTGS$targets$GErep)) {
				exprsNORM <- ddTGS$TGS
				
				for (u in unique(ddTGS$targets$GErep)) {
					j <- ddTGS$targets$GErep == u
					if (length(j[j==TRUE]) > 1) {
						print(paste("Applying", NORMmethod,"normalization between arrays",paste(colnames(exprsNORM[,j]),collapse=","), sep=" "))
						exprsNORM[,j] <- my.normalizeBetweenArrays(ddTGS$TGS[,j], method = 'quantile')
					}
				}
				ddNORM$TGS = log2(exprsNORM)
			}
			else {
				stop("GErep in targets of uRNAList is null")
			}
			
		}
		else if (NORMmethod == "Tscale75thp") {
			if (!is.null(ddTGS$targets$GErep)) {
				exprsNORM <- ddTGS$TGS
				
				for (u in unique(ddTGS$targets$GErep)) {
					j <- ddTGS$targets$GErep == u
					if (length(j[j==TRUE]) > 1) {
						print(paste("Applying", NORMmethod,"normalization between arrays",paste(colnames(exprsNORM[,j]),collapse=","), sep=" "))
						exprsNORM[,j] <- my.normalizeBetweenArrays(ddTGS$TGS[,j], method = 'scale75thp')
					}
				}
				ddNORM$TGS = exprsNORM
			}
			else {
				stop("GErep in targets of uRNAList is null")
			}			
		}				
		else {
			exprsNORM = normalizeBetweenArrays(ddTGS$TGS, method = NORMmethod)
			ddNORM$TGS = exprsNORM
		}
		rm(exprsNORM)
	}
	else {
		ddNORM = ddTGS
		ddNORM$TGS = log2(ddTGS$TGS)
	}
	if (verbose) {
		cat("------------------------------------------------------", 
				"\n")
		cat("\tNORMMALIZATION:\t", NORMmethod, "\n")
		cat("\tOUTPUT in log-2 scale", "\n")
		cat("------------------------------------------------------", 
				"\n")
	}
	if (!missing(makePLOTpost)) {
		if (makePLOTpost) {
			colorfill = "red"
			maintitle = "NORMALIZED SIGNAL"
			MMM = ddNORM$TGS
			dev.set(dev.next())
			plotDensityMicroRna(MMM, maintitle)
			dev.set(dev.next())
			boxplotMicroRna(MMM, maintitle, colorfill)
			dev.set(dev.next())
			mvaBASIC(MMM, colorfill, maintitle)
			dev.set(dev.next())
			hierclusMicroRna(MMM, targets$GErep, methdis = "euclidean", 
					methclu = "complete", sel = FALSE, 100)
			maintitle = "NORM DATA - RLE "
			dev.set(dev.next())
			RleMicroRna(MMM, maintitle, colorfill)
			rm(MMM)
		}
	}
	ddNORM$TGS = ddNORM$TGS
	ddNORM$TGS = ddNORM$TGS
	ddNORM$TGS = ddNORM$TGS
	return(ddNORM)
}

my.esetMicroRna <- function (uRNAList, targets, makePLOT = FALSE, verbose = FALSE, unique.rownames, signal="TGS") {
	if (!is(uRNAList, "uRNAList")) {
		stop("'input' must be a uRNAList", call. = FALSE)
		if (is.null(dim(uRNAList)[1])) {
			stop("'input' is empty", call. = FALSE)
		}
	}
	if (missing(targets)) {
		stop("'targets' is missing ", call. = FALSE)
	}
	if ("GErep" %in% colnames(targets)) {
		GErep = targets$GErep
		nGE = sum(table(table(GErep)))
		g1 = targets$GErep
		g2 = rownames(targets)
	}
	else {
		stop("'targets' needs 'GErep' field")
	}
	goON = all(rownames(targets) == colnames(uRNAList[[signal]]))
	if (!goON) {
		stop("rownames in pData(targets) different from colnames uRNAList[[",signal,"]]", 
				call. = FALSE)
	}
	phenoData = new("AnnotatedDataFrame", data = targets)
	TMP = uRNAList[[signal]]
	if (!missing(unique.rownames)) {
		rownames(TMP) = unique.rownames
	} 
	else {
		rownames(TMP) = uRNAList$genes$ProbeName
		#rownames(TMP) = uRNAList$genes$GeneName
	}	
	colnames(TMP) = g2
	nGEN = dim(TMP)[1]
	nARR = dim(TMP)[2]
	esetPROC = new("ExpressionSet", exprs = TMP, phenoData = phenoData)
	if (verbose) {
		cat("outPUT DATA: esetPROC", "\n")
		print(dim(esetPROC))
		cat("------------------------------------------------------", 
				"\n")
	}
	if (!missing(makePLOT)) {
		if (makePLOT) {
			dev.set(dev.next())
			hierclusMicroRna(exprs(esetPROC), GErep, methdis = "euclidean", 
					methclu = "complete", sel = FALSE, 100)
			dev.set(dev.next())
			size = 100
			maintitle = "50 High Variance"
			HeatMapMicroRna(exprs(esetPROC), 50, maintitle)
			dev.set(dev.next())
			PCAplotMicroRna(esetPROC, targets)
		}
	}
	return(esetPROC)
}


#' run.miRNA.TGS.extraction - Run signal extraction based on f.sigprep
#' @param f.sigprep, f.targets, f.dd
#' @returnType RGList
#' @return f.dd
#' @author daniel
#' @export
run.miRNA.TGS.extraction <- function(f.targets, f.dd.in, f.vNORMmethod, do.rma=FALSE,f.half=FALSE,f.offset=0 ) {
	if (do.rma) {
		opt <- 'RMA'
		print(paste("Run analysis: TGS from ",opt, "; normalized by ",f.vNORMmethod,sep=""))
	} else {
		opt <- 'AFE'
		print(paste("Run analysis: TGS from ",opt, "; normalized by ",f.vNORMmethod,sep=""))
	}
	###
	# Agilent Feature Extraction
	###
	
	if (is.null(f.dd.in)) {
		f.dd.in <- readMicroRnaAFE(f.targets, verbose=TRUE)
	}
	
	x11()
	plotDensityMicroRna(log2(f.dd.in$TGS), maintitle = paste("log2 TGS (",opt,") Signal",sep=""))
	graphics.off()
	
	###
	# qcPlots
	###
	
	if (  length(list.files("quality",pattern="^qcPlots")) == 0 ) {
		png(filename=paste("quality",paste("qcPlots","_%d.png",sep=""),sep="/"), bg="white", res=300, width=3000, height=3000)
		my.qcPlots(f.dd.in, offset=5, MeanSignal=TRUE, ProcessedSignal=TRUE, 
				TotalProbeSignal=TRUE, TotalGeneSignal=TRUE, BGMedianSignal= TRUE, 
				BGUsed=TRUE, targets)
		graphics.off()		
	}
	###
	# cvArray
	###
	
	if ( length(list.files("quality",pattern="^cvArray")) == 0  ) {
		f.dd.tmp <- f.dd.in
		for (i in 1:dim(f.dd.tmp)[2]) {
			index = which(f.dd.tmp$procS[, i] < 0.1)
			f.dd.tmp$procS[index, i] <- 0.1	
		}
		
		for (s in c('MeanSignal','ProcessedSignal')) {
			png(filename=paste("quality",paste("cvArray",s,"_%d.png",sep=""),sep="/"), bg="white", res=300, width=3000, height=3000)
			cvArray(f.dd.tmp, s, f.targets, verbose=TRUE)	
			graphics.off()
		}			
		rm(f.dd.tmp)
	}
	
	###
	# Pre-processing
	###
	
	## Change rownames based on the other columns
	
	con <- textConnection(as.character(f.targets$Subject))	# text connection
	subject <- read.table(con,sep='.')
	close(con)
	
	strReverse <- function(x) sapply(lapply(strsplit(x, NULL), rev), paste, collapse="") 
	
	aux <- gsub(".txt", "", as.character(f.targets$FileName)) 
	con <- textConnection(as.character(strReverse(aux)))
	filename <- read.table(con,sep='_',fill=TRUE)
	close(con)
	rm(aux)
	
	aux <- paste(f.targets$Treatment, subject$V1, f.targets$Array, filename$V2, filename$V1,sep = '.')	# criando rownames ('condicao'.'nro da microdissec��o')
	rm(subject, filename)
	
	rownames(f.targets) <- aux
	rm(aux)
	
	# Add column 'Name' in f.targets 
	f.targets$Name <- row.names(f.targets)
	
	# Replicate f.targets in f.dd.in
	f.dd.in <- rep.targets(f.dd.in, f.targets)
	
	if (! file.exists(paste('quality/aQM',sep=""))) {
		library('matlab')
		detach("package:matlab")
		
		## Creating the Expression Set Object (before background correction and normalization)		
		esetRAW=my.esetMicroRna(f.dd.in,f.targets,makePLOT=FALSE, verbose=TRUE,
				unique.rownames=paste(f.dd.in$genes$ProbeName, seq(1:length(f.dd.in$genes$ProbeName)), sep="/"))
		dim(esetRAW)
		library("matlab")
		
		## Quality analysis (before background correction and normalization)
		#featureData(esetRAW)$X <- f.dd.in$other$Row[,1]
		#featureData(esetRAW)$Y <- f.dd.in$other$Col[,1]
		
		aqm.raw <- arrayQualityMetrics(	expressionset = esetRAW,
				outdir = paste('quality/aQM',sep=""),
				force = TRUE,
				do.logtransform = TRUE,
				intgroup=c("Treatment","Subject","Array")
		)
		dev.off()
		rm(esetRAW)
		rm(aqm.raw)
	}
	
	if (do.rma) {		
		f.dd <- rmaMicroRna(f.dd.in, normalize=FALSE, background=TRUE)		
	} else {
		f.dd <- my.tgsMicroRna(f.dd.in, half=f.half, offset=f.offset, verbose=TRUE)
	}
	
	f.dd <- rep.targets(f.dd, f.targets)
	f.dd$targets <- f.targets
	
	if (! file.exists(paste('quality/aQM_',opt,sep=""))) {
		library('matlab')
		detach("package:matlab")
		
		## Creating the Expression Set Object (before background correction and normalization)		
		esetRAW.2=my.esetMicroRna(f.dd,f.targets,makePLOT=FALSE, verbose=TRUE,
				unique.rownames=paste(f.dd$genes$GeneName, seq(1:length(f.dd$genes$GeneName)), sep="/"))
		dim(esetRAW.2)
		library("matlab")
		
		## Quality analysis (before background correction and normalization)
		#featureData(esetRAW.2)$X <- f.dd$other$Row[,1]
		#featureData(esetRAW.2)$Y <- f.dd$other$Col[,1]
		
		aqm.raw.2 <- arrayQualityMetrics(	expressionset = esetRAW.2,
				outdir = paste('quality/aQM_',opt,sep=""),
				force = TRUE,
				do.logtransform = TRUE,
				intgroup=c("Treatment","Subject","Array")
		)
		dev.off()
		rm(esetRAW.2)
		rm(aqm.raw.2)
	}
	
	f.dd$targets <- f.targets 
	
	# Background Correction and Normalization
	f.ddNORM=my.tgsNormalization(f.dd,NORMmethod=f.vNORMmethod, targets=f.targets,	makePLOTpre=FALSE,makePLOTpost=FALSE)
	
	x11()
	plotDensityMicroRna(f.ddNORM$TGS, maintitle = paste("log2 Normalized TGS (",opt,") Signal",sep=""))
	graphics.off()

	if (! file.exists(paste('quality/aQM_',opt,'_',f.vNORMmethod,sep=""))) {
		library('matlab')
		detach("package:matlab")
		
		## Creating the Expression Set Object (before background correction and normalization)		
		esetNORM=my.esetMicroRna(f.ddNORM,f.targets,makePLOT=FALSE, verbose=TRUE,
				unique.rownames=paste(f.ddNORM$genes$GeneName, seq(1:length(f.ddNORM$genes$GeneName)), sep="/"))
		dim(esetNORM)
		library("matlab")
		
		## Quality analysis (before background correction and normalization)
		#featureData(esetNORM)$X <- f.ddNORM$other$Row[,1]
		#featureData(esetNORM)$Y <- f.ddNORM$other$Col[,1]
		
		aqm.norm <- arrayQualityMetrics(	expressionset = esetNORM,
				outdir = paste('quality/aQM_',opt,'_',f.vNORMmethod,sep=""),
				force = TRUE,
				do.logtransform = FALSE,
				intgroup=c("Treatment","Subject","Array")
		)
		dev.off()
		rm(esetNORM)
		rm(aqm.norm)
	}
	
	###
	# Filtering probes
	###
	f.ddFILT=filterMicroRna(	ddNORM=f.ddNORM, 
								dd=f.dd.in,
								control=TRUE,
								IsGeneDetected=FALSE,
								wellaboveNEG=FALSE,
								limIsGeneDetected=25,
								limNEG=25,
								makePLOT=FALSE,
								targets=f.targets,
								verbose=TRUE,
								writeout=TRUE) 
	

	
	for (f in c('IsNOTGeneDetected.txt', 'NOCtrl_exprs.txt', 'NOCtrl_FlagIsGeneDetected.txt',
				'IsNOTWellAboveNEG.txt','')) {
				
		if (file.exists(f)) {
			system(paste("mv ",f," results/tables/ 2>&1"), intern=TRUE)		
		}
	}
	
	if (! file.exists(paste('quality/aQM_',opt,'_',f.vNORMmethod,'_filtered',sep=""))) {
		library('matlab')
		detach("package:matlab")
		
		## Creating the Expression Set Object (before background correction and normalization)		
		esetFILT=my.esetMicroRna(f.ddFILT,f.targets,makePLOT=FALSE, verbose=TRUE,
				unique.rownames=paste(f.ddFILT$genes$GeneName, seq(1:length(f.ddFILT$genes$GeneName)), sep="/"))
		dim(esetFILT)
		library("matlab")
		
		## Quality analysis (before background correction and normalization)
		#featureData(esetFILT)$X <- f.ddFILT$other$Row[,1]
		#featureData(esetFILT)$Y <- f.ddFILT$other$Col[,1]
		
		aqm.filt <- arrayQualityMetrics(	expressionset = esetFILT,
				outdir = paste('quality/aQM_',opt,'_',f.vNORMmethod,'_filtered',sep=""),
				force = TRUE,
				do.logtransform = FALSE,
				intgroup=c("Treatment","Subject","Array")
		)
		dev.off()
		rm(esetFILT)
		rm(aqm.filt)
	}
	
	esetFILT = esetMicroRna(f.ddFILT, f.targets, verbose = TRUE)
		
	## Writing the Expression Set Object: Processeddata.txt
	writeEset(esetFILT, f.ddFILT, f.targets,verbose=TRUE)
	
	if (file.exists("ProcessedData.txt")) {
		system("mv ProcessedData.txt results/tables/ 2>&1", intern=TRUE)
	}
	
	## Principal Component Analysis
	#png(filename=paste("quality",paste("PCA_",opt,"_",f.vNORMmethod,'_filtered',"_%d.png",sep=""),sep="/"), bg="white", res=300, width=3000, height=3000)
	#PCAplotMicroRna(esetFILT,f.targets)
	#graphics.off()
		
	return(list(dd=f.dd.in,ddNORM=f.ddNORM,ddFILT=f.ddFILT))
}




###############################################################################
# Date: Jul 18, 2011
# USP-Ribeirao Preto, Wilson's Lab
# 
# Author: matheus
###############################################################################


#
#
#`filter.IsGeneDetected` <-
#		function(ddFILT,limIsGeneDetected,targets,verbose,writeout){
#	
#	
#	minFLAGisf=1 	# gIsFound: FLAG ok: 1 = feature FOUND (58)
#	
#	if (!is(ddFILT, "uRNAList")){
#		stop("'input' must be a uRNAList")
#		if (is.null(dim(ddFILT)[1])) {
#			stop("'input' is empty")
#		}
#	}
#	
#	if(missing(targets)){
#		stop("'targets' is missing ")
#	}
#	
#	if("GErep" %in% colnames(targets)){
#		g1=targets$GErep  # g1 must be numeric, from 1:n
#		g2=rownames(targets)
#		GErep=targets$GErep 
#		nGE=sum(table(table(GErep)))
#	}else{
#		stop("'targets' needs 'GErep' field")
#	}
#	
#	FLAG=ddFILT$other$gIsGeneDetected
#	indexSNR=apply(FLAG,1,filterFLAG.micro,GErep,nGE,minFLAGisf,limIsGeneDetected)
#	selSNR=which(unlist(indexSNR) == 1) # posiciones que pasan filtro 
#	
## WRITING OUT THE REMOVED DATA only geneIDs and selected flags
#	
#	if(writeout){
#		outfile="IsNOTGeneDetected.txt" 
#		write.filt.out.miRNA(ddFILT,selSNR,outfile,FLAG,targets) 
#	}
#	if(verbose){
#		cat("FILTERING BY IsGeneDetected FLAG","\n") 
#		cat("\n")
#		cat("	FLAG FILTERING OPTIONS - FLAG OK = 1 - limIsGeneDetected: ",limIsGeneDetected,"%","\n")
#		cat("	FEATURES AFTER IsGeneDetected FILTERING: ",sum(indexSNR),"\n")
#		cat("	NON Gene Detected :",length(ddFILT$genes$ProbeName[-selSNR]),"\n")
#		cat("------------------------------------------------------","\n")
#	}
#	
## EXTRACTING THE SELECTED PROBES  
#	
#	ddFILT=ddFILT[selSNR,]	 
#	return(ddFILT)
#	
#}  # end function 
#
#`filter.control.miRNA` <-
#		function(ddNORM,targets,verbose,writeout){
#	
#	if (!is(ddNORM, "uRNAList")){
#		stop("'input' must be a uRNAList")
#		if (is.null(dim(ddNORM)[1])) {
#			stop("'input' is empty")
#		}
#	}
#	
#	if(missing(targets)){
#		stop("'targets' is missing ")
#	}
#	
#	if("GErep" %in% colnames(targets)){
#		g1=targets$GErep  # g1 must be numeric, from 1:n
#		g2=rownames(targets)
#	}else{
#		stop("'targets' needs 'GErep' field")
#	}
#	
#	cat("\n")
#	cat("FILTERING BY ControlType","\n")
#	
#	selSNR=which(ddNORM$genes$ControlType==0)
#	
#	ddFILT=ddNORM[selSNR,]
#	
#	
## WRITING OUT THE RAW DATA WITHOUT CONTROLS 
#	
#	if(writeout){
#		write.control.out.miRNA(ddFILT,selSNR,targets)
#	}	
#	if(verbose){
#		cat("\n")
#		cat("   FEATURES BEFORE FILTERING: ",dim(ddNORM)[1],"\n")
#		cat(" 	FEATURES AFTER ControlType FILTERING: ",dim(ddFILT)[1],"\n")
#		cat("------------------------------------------------------","\n")
#	}
#	
#	return(ddFILT)
#} # end function 
#
#`filter.wellaboveNEG.miRNA` <-
#		function(ddFILT,dd,limNEG,SDtimes,targets,verbose,writeout){
#	
#	
#	if (!is(ddFILT, "uRNAList")){
#		stop("'input ddFILT' must be a uRNAList")
#		if (is.null(dim(ddFILT)[1])) {
#			stop("'input' is empty")
#		}
#	}
#	
#	if (!is(dd, "uRNAList")){
#		stop("'input dd' must be a uRNAList")
#		if (is.null(dim(dd)[1])) {
#			stop("'input' is empty")
#		}
#	}
#	
#	if(missing(targets)){
#		stop("'targets' is missing ")
#	}
#	
#	if("GErep" %in% colnames(targets)){
#		g1=targets$GErep  # g1 must be numeric, from 1:n
#		g2=rownames(targets)
#		GErep=targets$GErep 
#		nGE=sum(table(table(GErep)))
#	}else{
#		stop("'targets' needs 'GErep' field")
#	}
#	
#	indexneg=which(dd$genes$ControlType == -1)
#	NEG=log2(dd$meanS[indexneg,])  #  MeanSignal in original dd 
#	MeanNeg=apply(NEG,2,mean)
#	SdNeg=apply(NEG,2,sd)
#	Limit=MeanNeg + SDtimes*(SdNeg)
#	
#	
#	FLAG=ddFILT$meanS 
#	indexSNR=apply(FLAG,1,filterWellAboveSIGNALv2,GErep,nGE,Limit,limNEG)
#	selSNR=which(unlist(indexSNR) == 1) #posiciones que pasan filtro 
#	
#	
## WRITING OUT THE WellAboveNeg.out filtered DATA 
#	
#	if(writeout){
#		outfile="IsNOTWellAboveNEG.txt"
#		VALUES=round(ddFILT$meanS,3)
#		write.filt.out.miRNA(ddFILT,selSNR,outfile,VALUES,targets) 
#	}
#	
#	
#	if(verbose){
#		cat("FILTERING BY WellAboveNeg filterWellAboveSIGNALv2 ~ FLAG","\n")
#		cat("\n")
#		cat("	FLAG FILTERING OPTIONS - limNEG: ",limNEG,"%","\n")
#		cat("	Limit computed as MeanNeg + ",SDtimes," x (SDNeg) ","\n")
#		
#		cat("	Limit in each array: ",round(Limit,2),"\n")
#		cat("\n")
#		cat("PROBES AFTER WellAboveNeg FILTERING: ",sum(indexSNR),"\n")
#		cat("WellAboveNeg OUT :",length(ddFILT$genes$ProbeName[-selSNR]),"\n")
#		cat("------------------------------------------------------","\n")
#	}
#	
## EXTRACTING THE SELECTED PROBES  
#	
#	ddFILT=ddFILT[selSNR,]
#	return(ddFILT)
#	
#} # end function 
#
#`filterFLAG.micro` <-
#		function(flag,GErep,nGE,minFLAG,limSNR) {
#	
#	index=0
#	ii=0
#	while(ii < nGE) {
#		ii=ii+1
#		minSize=length(which(GErep == unique(GErep)[ii]))*(limSNR/100)  
#		pos=which(GErep == ii)  
#		
#		aux=flag[pos]
#		nan=!is.na(aux)         
#		aux=aux[nan]           
#		lp=length(aux)         
#		
#		sumSNR=0
#		jj=0
#		while(jj < lp ){
#			jj=jj+1
#			
#			if(aux[jj] == minFLAG){
#				sumSNR=sumSNR+1
#				if(sumSNR >= minSize){
#					index=1
#					jj=lp + 1  # exit the while(jj)
#					ii=nGE +1  # exit the while(ii)
#				}
#			}
#		}
#	}
#	
#	return(index)
#} # end of function 
#
#
#`filterWellAboveSIGNALv2` <-
#		function(flag,GErep,nGE,Limit,limNEG) {
#	
#	index=0
#	ii=0
#	while(ii < nGE) {
#		ii=ii+1
#		minSize=length(which(GErep == unique(GErep)[ii]))*(limNEG/100)  
#		pos=which(GErep == ii) 
#		
#		aux=flag[pos]
#		nan=!is.na(aux)        
#		aux=aux[nan]            
#		lp=length(aux)          
#		
#		sumSNR=0
#		jj=0
#		while(jj < lp ){
#			jj=jj+1
#			
#			if(aux[jj] >= Limit[pos[jj]]){
#				sumSNR=sumSNR+1
#				if(sumSNR >= minSize){
#					index=1
#					jj=lp + 1  # exit the while(jj)
#					ii=nGE +1  # exit the while(ii)
#				}
#			}
#		}
#	}    
#	return(index)
#} # end of function
#
#densMiR <- function (object, maintitle) 
#{
#	samples = colnames(object)
#	nARR = dim(object)[2]
#	colors <- rainbow(nARR, s = 1, v = 1, start = 0, end = max(1, 
#					nARR - 1)/nARR, gamma = 1)
#	y.max = c()
#	x.max = c()
#	for (n in 1:nARR) {
#		y.max[n] = max(density(object[, n], na.rm = TRUE)$y)
#		x.max[n] = max(density(object[, n], na.rm = TRUE)$x)
#	}
#	y.pos = order(y.max, decreasing = TRUE, na.last = NA)
#	x.pos = order(x.max, decreasing = TRUE, na.last = NA)
#	for (n in y.pos) {
#		k = which(y.pos == n)
#		if (n == y.pos[1]) 
#			plot(density(object[, n], na.rm = TRUE), col = colors[n], 
#					main = "", asp = 0.7 * x.max[x.pos[1]]/y.max[y.pos[1]])
#		else lines(density(object[, n], na.rm = TRUE), col = colors[n])
#	}
#	title(main = maintitle)
#	if(length(samples) > 0 ){
#		legend(x = "topright", legend = samples, cex = 0.8, fill = colors, inset = 0.05)
#	}
#}
#
#filterMiR <- function (ddNORM, dd, control, IsGeneDetected, wellaboveNEG, 
#		limIsGeneDetected, limNEG, makePLOT, targets, verbose, writeout) 
#{
#	if (verbose) {
#		cat("FILTERING PROBES BY FLAGS", "\n")
#		cat("\n")
#	}
#	if (missing(targets)) {
#		stop("'targets' is missing ")
#	}
#	if ("GErep" %in% colnames(targets)) {
#		g1 = targets$GErep
#		g2 = rownames(targets)
#		GErep = targets$GErep
#		nGE = sum(table(table(GErep)))		
#	}
#	else {
#		stop("'targets' needs 'GErep' field")
#	}
#	if (!is(ddNORM, "uRNAList")) {
#		stop("'input' must be a uRNAList")
#		if (is.null(dim(ddNORM)[1])) {
#			stop("'input' is empty")
#		}
#	}
#	if (!missing(control)) {
#		if (control) {
#			ddFILT = filter.control.miRNA(ddNORM, targets, verbose, 
#					writeout)
#		}
#	}
#	if (!missing(IsGeneDetected)) {
#		if (IsGeneDetected) {
#			ddFILT = filter.IsGeneDetected(ddFILT, limIsGeneDetected, 
#					targets, verbose, writeout)
#		}
#	}
#	if (!missing(wellaboveNEG)) {
#		if (wellaboveNEG) {
#			ddFILT = filter.wellaboveNEG.miRNA(ddFILT, dd, limNEG, 
#					SDtimes = 1.5, targets, verbose, writeout)
#		}
#	}
#	if(!is.null(colnames(ddNORM$TGS)))
#		colnames(ddFILT$TGS) <- colnames(ddNORM$TGS)
#	if (!missing(makePLOT)) {
#		if (makePLOT) {
#			colorfill = "yellow"
#			maintitle = "FILTERED SIGNAL"
#			dev.set(dev.next())
#			densMiR(ddFILT$TGS, maintitle)
#			dev.set(dev.next())
#			boxplotMicroRna(ddFILT$TGS, maintitle, colorfill)
#			dev.set(dev.next())
#			hierclusMicroRna(ddFILT$TGS, paste(targets$Subject, targets$GErep, sep=" - "), methdis = "euclidean", 
#					methclu = "complete", sel = FALSE, 100)
#			maintitle = "FILTERED SIGNAL - RLE "
#			dev.set(dev.next())
#			RleMicroRna(ddFILT$TGS, maintitle, colorfill)
#		}
#	}
#	return(ddFILT)
#}
#
#library(marray)
#
#
#
#getEset <- function (uRNAList, targets, verbose = FALSE) 
#{
#	if (!is(uRNAList, "uRNAList")) {
#		stop("'input' must be a uRNAList", call. = FALSE)
#		if (is.null(dim(uRNAList)[1])) {
#			stop("'input' is empty", call. = FALSE)
#		}
#	}
#	if (missing(targets)) {
#		stop("'targets' is missing ", call. = FALSE)
#	}
#	if ("GErep" %in% colnames(targets)) {
#		GErep = targets$GErep
#		nGE = base::sum(table(table(GErep)))
#		g1 = targets$GErep
#		g2 = rownames(targets)
#	}
#	else {
#		stop("'targets' needs 'GErep' field")
#	}
#	goON = all(rownames(targets) == colnames(uRNAList$iTGS))
#	if (!goON) {
#		stop("rownames in pData(targets) different from colnames uRNAList$TGS", 
#				call. = FALSE)
#	}
#	phenoData = new("AnnotatedDataFrame", data = targets)
#	TMP = uRNAList$TGS
#	rownames(TMP) = uRNAList$genes$GeneName
#	colnames(TMP) = g2
#	nGEN = dim(TMP)[1]
#	nARR = dim(TMP)[2]
#	esetPROC = new("ExpressionSet", exprs = TMP, phenoData = phenoData)
#	if (verbose) {
#		cat("outPUT DATA: esetPROC", "\n")
#		print(dim(esetPROC))
#		cat("------------------------------------------------------", 
#				"\n")
#	}
#	return(esetPROC)
#}
#
#pvalHisto <- function (fit2, DE, PVcut, DEmethod, MTestmethod, CM, verbose = FALSE) 
#{
#	nCON = dim(DE)[2]
#	if (nCON > 1 && DEmethod == "nestedF") {
#		ord = which(p.adjust(fit2$F.p.value, method = MTestmethod) <= 
#						PVcut)
#		if (verbose) {
#			cat("num F Contrasts:", nCON, "\n")
#			cat("DEG Significants by F.p.value <= ", PVcut, "-", 
#					length(ord), MTestmethod, "\n")
#			cat("------------------------------------------------------", 
#					"\n")
#		}
#		aux = paste.character(colnames(CM))
#		maintitle = paste("F-test ", aux, DEmethod, sep = " - ")
#		dev.set(dev.next())
#		hist(fit2$F.p.value, freq = TRUE, col = "yellow", main = maintitle, 
#				xlab = "F test p.value", ylab = " freq ")
#	}
#	if (nCON == 1 || DEmethod == "separate") {
#		for (i in 1:nCON) {
#			maintitle = paste(colnames(CM)[i], DEmethod, sep = " - ")
#			dev.set(dev.next())
#			hist(fit2$p.value[, i], freq = TRUE, col = "green", 
#					main = maintitle, xlab = "t test p.value", ylab = " freq ")
#		}
#	}
#}
#
#
##' getDecideTests - A interface for decideTests from limma package
##' @param fit MArrayLM object
##' @param DEmethod method for decideTests, only 'separate' or 'nestedF' are implemented. see decideTests in limma package.
##' @param MTestmethod method for multiple test, choices are 'none','BH', 'BY', ... see p.adjust
##' @param PVcut p-value threshold to declare significant features
##' @param lfc minimum log2-fold-change required 
##' @returnType NULL
##' @return none
##' @author matheus
##' @export
#getDecideTests <- function (fit2, DEmethod, MTestmethod, PVcut, lfc = 0, verbose = FALSE) 
#{
#	if (!is(fit2, "MArrayLM")) {
#		stop("'design' must be a 'MArrayLM")
#		if (is.null(length(fit2$coefficients))) {
#			stop("fit2' is empty")
#		}
#	}
#	if (missing(DEmethod)) {
#		stop(" method for decideTests 'separate' or 'nestedF' is needed")
#	}
#	if (missing(MTestmethod)) {
#		stop(" method for multiple test 'none','BH' or 'BY' is needed")
#	}
#	if (missing(PVcut)) {
#		stop("'PVcut' is missing")
#	}
#	DE = decideTests(fit2, method = DEmethod, adjust.method = MTestmethod, 
#			p.value = PVcut, lfc = lfc)
#	sumDE = summary(DE)
#	rownames(sumDE)[1] = "DOWN"
#	rownames(sumDE)[3] = "UP"
#	sumDE = sumDE[c(3, 1), ]
#	if (verbose) {
#		cat("\n")
#		cat("------------------------------------------------------", 
#				"\n")
#		cat(" Method for Selecting DEGs:", DEmethod, "\n")
#		cat(" Multiple Testing  method: ", MTestmethod, "- pval", 
#				PVcut, " - minimum log2(foldc) ", lfc, " \n")
#		cat("\n")
#		print(sumDE)
#		cat("------------------------------------------------------", 
#				"\n")
#	}
#	return(DE)
#}
#
### Rascunho
#
#
#MA.plot.miRNA <- function (fit2, DE, CM, i) 
#{
#	posDE = which(DE[, i] > 0)
#	negDE = which(DE[, i] < 0)
#	M = fit2$coef[, i]
#	A = fit2$Amean
#	dev.set(dev.next())
#	plot(A, M, cex = 1, col = "black", pch = 20)
#	points(A, M, cex = 0.6, col = "cyan2", pch = 19)
#	points(A[posDE], M[posDE], cex = 0.6, col = "red", pch = 19)
#	points(A[negDE], M[negDE], cex = 0.6, col = "blue", pch = 19)
#	smooth.fit = fitted(loess(M ~ A))
#	points(A, smooth.fit, col = "black", cex = 0.5, pch = 19)
#	colors = c("red", "blue", "black")
#	samples = c("Up", "Down", "M~A smooth fit")
#	legend(x = "topright", legend = samples, cex = 0.8, fill = colors, 
#			inset = 0.05)
#	nn = paste(i, " - ", colnames(CM)[i])
#	title(main = nn)
#}
#
##'
##' @param MTestMethod2 Parametro utilizado para escolher os genes que terao seus nomes plotados
#significantMir <- function (eset, ddset, targets, fit2, CM, DE, DEmethod, MTestmethod, 
#		PVcut, lfc = 0, verbose = FALSE, MTestmethod2= "none") 
#{
#	require("geneplotter")
#	require("codelink")
#	
#	if (!is(eset, "ExpressionSet")) {
#		stop("'eset' must be a ExpressionSet", call. = FALSE)
#		if (is.null(nrow(exprs(eset)))) {
#			stop("'eset' is empty", call. = FALSE)
#		}
#	}
#	if (!is(ddset, "uRNAList")) {
#		stop("'ddset' must be a uRNAList", call. = FALSE)
#		if (is.null(dim(ddset)[1])) {
#			stop("'ddset' is empty", call. = FALSE)
#		}
#	}
#	if (missing(targets)) {
#		stop("'targets' is missing ", call. = FALSE)
#	}
#	else {
#		if ("GErep" %in% colnames(targets)) {
#			GErep = targets$GErep
#			nGE = sum(as.numeric(table(table(GErep))))
#			g1 = targets$GErep
#			g2 = rownames(targets)
#		}
#		else {
#			stop("'targets' needs 'GErep' field", call. = FALSE)
#		}
#	}
#	if (!is(fit2, "MArrayLM")) {
#		stop("'fit2' must be a 'MArrayLM", call. = FALSE)
#		if (is.null(length(fit2$coefficients))) {
#			stop("fit2' is empty", call. = FALSE)
#		}
#	}
#	if (!is(CM, "matrix")) {
#		stop("'CM' must be a 'Contrast Matrix'", call. = FALSE)
#		if (is.null(dim(CM))) {
#			stop("'CM' is empty", call. = FALSE)
#		}
#	}
#	if (!is(DE, "TestResults")) {
#		stop("'DE' must be a 'TestResults'", call. = FALSE)
#		if (is.null(dim(DE))) {
#			stop("'DE' is empty", call. = FALSE)
#		}
#	}
#	if (missing(DEmethod)) {
#		stop(" method for decideTests 'separate' or 'nestedF' is needed", 
#				call. = FALSE)
#	}
#	if (missing(MTestmethod)) {
#		stop(" method for multiple test 'none','BH' or 'BY' is needed", 
#				call. = FALSE)
#	}
#	if (missing(PVcut)) {
#		stop("'PVcut' is missing", call. = FALSE)
#	}
#	if (verbose) {
#		cat("------------------------------------------------------", 
#				"\n")
#	}
#	nGEN = dim(eset)[1]
#	nARR = dim(eset)[2]
#	nCON = dim(CM)[2]
#	for (i in 1:nCON) {
#		if (verbose) {
#			cat("CONTRAST: ", i, " - ", colnames(CM)[i], "\n")
#			cat("\n")
#		}
#		notDE = which(DE[, i] == 0)
#		nDDEE = which(DE[, i] != 0)
#		if (length(nDDEE) > 0) {
#			AgiMicroRna:::DEG.print.info(eset, DE, i, verbose)
#		}
#		method = match.arg(DEmethod, c("separate", "nestedF"))
#		switch(method, separate = {
#					adj.pval = p.adjust(fit2$p.value[, i], method = MTestmethod)
#					fdr = p.adjust(fit2$p.value[, i], method = "BH")
#					ord.pval = order(fit2$p.value[, i], decreasing = F)
#					eset.ord = eset[ord.pval, ]
#					fit2.ord = fit2[ord.pval, ]
#					ddset.ord = ddset[ord.pval, ]
#					adj.pval = adj.pval[ord.pval]
#					fdr = fdr[ord.pval]
#				}, nestedF = {
#					adj.Fpval = p.adjust(fit2$F.p.value, method = MTestmethod)
#					fdr.Fpval = p.adjust(fit2$F.p.value, method = "BH")
#					ord.pval.DE = nDDEE[order(fit2$F.p.value[nDDEE], 
#									decreasing = F)]
#					ord.pval.notDE = notDE[order(fit2$F.p.value[notDE], 
#									decreasing = F)]
#					eset.ord.DE = eset[ord.pval.DE, ]
#					fit2.ord.DE = fit2[ord.pval.DE, ]
#					ddset.ord.DE = ddset[ord.pval.DE, ]
#					adj.Fpval.DE = adj.Fpval[ord.pval.DE]
#					fdr.Fpval.DE = fdr.Fpval[ord.pval.DE]
#					eset.ord.notDE = eset[ord.pval.notDE, ]
#					fit2.ord.notDE = fit2[ord.pval.notDE, ]
#					ddset.ord.notDE = ddset[ord.pval.notDE, ]
#					adj.Fpval.notDE = adj.Fpval[ord.pval.notDE]
#					fdr.Fpval.notDE = fdr.Fpval[ord.pval.notDE]
#				})
#		nDDEE.ord = seq(1:length(nDDEE))
#		notDE.ord = seq((length(nDDEE) + 1):dim(eset)[1])
#		method = match.arg(DEmethod, c("separate", "nestedF"))
#		switch(method, separate = {
#					DEGLIST.ALL = paste("LIST.ALL", DEmethod, MTestmethod, 
#							colnames(CM)[i], "txt", sep = ".")
#					index.ALL = seq(1:dim(eset)[1])
#					AgiMicroRna:::write.LIST.miRNA(fit2.ord, index.ALL, DEGLIST.ALL, 
#							DEmethod, i, adj.pval, fdr, ddset.ord)
#					filename = paste(DEGLIST.ALL, "html", sep = ".")
#					genelist = fit2.ord$genes[adj.pval <= PVcut, ]
#					title = paste("DEGs", DEmethod, MTestmethod, PVcut, 
#							colnames(CM)[i], sep = ".")
#					head <- c("miRNA ID")
#					AgiMicroRna:::miRNA.htmlpage(genelist, filename, title, table.head = head, 
#							table.center = TRUE)
#					MA.plot.miRNA(fit2, DE, CM, i)
#					#Volcano
#					dev.set(dev.next())
#					x<-data.frame(fold=fit2$coef[,i],p=fit2$p.value[,i],p.adj=p.adjust(fit2$p.value[,i],method = MTestmethod2 ))
#					title.plot = paste(i, " - ", colnames(CM)[i])
#					posDE = which(DE[, i] > 0)
#					posDE = intersect(posDE,which(x[, "p.adj"] <= PVcut))
#					negDE = which(DE[, i] < 0)
#					negDE = intersect(negDE,which(x[, "p.adj"] <= PVcut))
#					func.list$volcano.plus(x,"fold","p", title.plot, PVcut, lfc, -lfc, "p.adj",ncolors=5,text=as.numeric(c(posDE,negDE)))
#					
#					#if(length(posDE) > 0) text(x$fold[posDE], -log10(x$p.adj[posDE]), labels=rownames(x)[posDE])
#				}, nestedF = {
#					DEGLIST = paste("DEGs", DEmethod, MTestmethod, PVcut, 
#							colnames(CM)[i], "txt", sep = ".")
#					NOTDE = paste("notDEG", Mcut, DEmethod, MTestmethod, 
#							PVcut, colnames(CM)[i], "txt", sep = ".")
#					AgiMicroRna:::write.LIST.miRNA(fit2.ord.DE, nDDEE.ord, DEGLIST, 
#							DEmethod, i, adj.Fpval.DE, fdr.Fpval.DE, ddset.ord.DE)
#					filename = paste(DEGLIST, "html", sep = ".")
#					genelist = fit2.ord.DE$genes[nDDEE.ord, ]
#					title = paste("DEGs", DEmethod, MTestmethod, PVcut, 
#							colnames(CM)[i], sep = ".")
#					head <- c("miRNA ID")
#					AgiMicroRna:::miRNA.htmlpage(genelist, filename, title, table.head = head, 
#							table.center = TRUE)
#					AgiMicroRna:::write.LIST.miRNA(fit2.ord.notDE, notDE.ord, NOTDE, 
#							DEmethod, i, adj.Fpval.notDE, fdr.Fpval.notDE, 
#							ddset.ord.notDE)
#					MA.plot.miRNA(fit2, DE, CM, i)
#					#Volcano
#					dev.set(dev.next())
#					x<-data.frame(fold=fit2$coef[,i],p=fit2$p.value[,i],p.adj=p.adjust(fit2$p.value[,i],method = MTestmethod2 ))
#					title.plot = paste(i, " - ", colnames(CM)[i])
#					posDE = which(DE[, i] > 0)
#					posDE = intersect(posDE,which(x[, "p.adj"] <= PVcut))
#					negDE = which(DE[, i] < 0)
#					negDE = intersect(negDE,which(x[, "p.adj"] <= PVcut))
#					func.list$volcano.plus(x,"fold","p", title.plot, PVcut, lfc, -lfc, "p.adj",ncolors=5,text=as.numeric(c(posDE,negDE)))
#					#if(length(posDE) > 0) text(x$fold[posDE], -log10(x$p.adj[posDE]), labels=rownames(x)[posDE])
#				})
#	}
#}

###############################################################################
# History:
# o - Initial creation, date: Jul 18, 2011
###############################################################################
