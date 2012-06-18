# LGMB   Laboratory of Genetics and Molecular Biology
# BiT    Bioinformatics Team
#
# Project...: myArray
# Date......: 22/07/2011
# Author(s).: daniel, matheus
###############################################################################

filter.isvalid <- function (ddNORM, ManuelaGO, targets, annotation.package) 
{
	if (!is(ddNORM, "RGList")) {
		stop("'input' must be a RGList")
		if (is.null(dim(ddNORM)[1])) {
			stop("'input' is empty")
		}
	}
	if (missing(targets)) {
		stop("'targets' is missing ")
	}
	if ("GErep" %in% colnames(targets)) {
		g1 = targets$GErep
		g2 = rownames(targets)
	}
	else {
		stop("'targets' needs 'GErep' field")
	}
	cat("\n")
	cat("FILTERING BY IsValid FLAG", "\n")
	selSNR = which(ddNORM$genes$gIsValid == 1)
	ddFILT = ddNORM[selSNR, ]
	#write.control.out(ddFILT, selSNR, annotation.package, ManuelaGO, 
	#		targets)
	cat("------------------------------------------------------", 
			"\n")
	cat(" PROBES BEFORE FILTERING: ", dim(ddNORM)[1], "\n")
	cat(" PROBES AFTER IsValid FILTERING: ", dim(ddFILT)[1], 
			"\n")
	cat(" RAW DATA WITHOUT INVALID OUT :", length(ddFILT$genes$ProbeName), 
			"\n")
	cat("------------------------------------------------------", 
			"\n")
	return(ddFILT)
}

my.filter.probes <- function (ddNORM, control = NULL, wellaboveBG = NULL, isfound = NULL, 
		wellaboveNEG = NULL, sat = NULL, PopnOL = NULL, NonUnifOL = NULL, 
		nas = NULL, limWellAbove, limISF, limNEG, limSAT, limPopnOL, 
		limNonUnifOL, limNAS, makePLOT, annotation.package, flag.counts = NULL, 
		targets) 
{
	cat("FILTERING PROBES BY FLAGS", "\n")
	cat("\n")
	if (!is(ddNORM, "RGList")) {
		stop("'input' must be a RGList")
		if (is.null(dim(ddNORM)[1])) {
			stop("'input' is empty")
		}
	}
	if (missing(annotation.package)) {
		stop("'annotation.package' is needed ")
	}
	else {
		if (annotation.package == "notAnnPack") {
			cat("NOT ANNOTATION PACKAGE AVAILABLE AT THIS MOMENT", 
					"\n")
			ManuelaGO = FALSE
		}
		else {
			if (annotation.package %in% .packages(all = TRUE)) {
			}
			else {
				cat("\tR will proceed to install the annotation package:", 
						annotation.package, "\n")
				install.packages(annotation.package)
			}
			ManuelaGO = TRUE
		}
	}
	if (missing(targets)) {
		stop("'targets' is missing ")
	}
	if (missing(limWellAbove)) {
		limWellAbove = 75
	}
	if (missing(limISF)) {
		limISF = 75
	}
	if (missing(limNEG)) {
		limNEG = 75
	}
	if (missing(limSAT)) {
		limSAT = 75
	}
	if (missing(limPopnOL)) {
		limPopnOL = 25
	}
	if (missing(limNonUnifOL)) {
		limNonUnifOL = 25
	}
	if (missing(limNAS)) {
		limNAS = 100
	}
	if (!missing(control)) {
		if (control) {
			ddFILT = my.filter.control(ddNORM, ManuelaGO, targets, 
					annotation.package)
		}
	}
	if (!missing(wellaboveBG)) {
		if (wellaboveBG) {
			ddFILT = filter.wellaboveBG(ddFILT, limWellAbove, 
					ManuelaGO, targets, annotation.package)
		}
	}
	if (!missing(isfound)) {
		if (isfound) {
			ddFILT = filter.isfound(ddFILT, limISF, ManuelaGO, 
					targets, annotation.package)
		}
	}
	if (!missing(wellaboveNEG)) {
		if (wellaboveNEG) {
			ddFILT = filter.wellaboveNEG(ddFILT, ddNORM, limNEG, 
					SDtimes = 1.5, ManuelaGO, targets, annotation.package)
		}
	}
	if (!missing(sat)) {
		if (sat) {
			ddFILT = filter.saturated(ddFILT, limSAT, ManuelaGO, 
					targets, annotation.package)
		}
	}
	if (!missing(PopnOL)) {
		if (PopnOL) {
			ddFILT = filter.PopnOL(ddFILT, limPopnOL, ManuelaGO, 
					targets, annotation.package)
		}
	}
	if (!missing(NonUnifOL)) {
		if (NonUnifOL) {
			ddFILT = filter.NonUnifOL(ddFILT, limNonUnifOL, ManuelaGO, 
					targets, annotation.package)
		}
	}
	if (!missing(nas)) {
		if (nas) {
			ddFILT = Agi4x44PreProcess:::filter.nas(ddFILT, limNAS, targets)
		}
	}
	if (!missing(flag.counts)) {
		if (flag.counts) {
			flag = 1
			FLAG = ddFILT$other$gIsWellAboveBG
			probeFLAG = apply(FLAG, 1, countFLAG, flag)
			rm(FLAG)
			cat("COUNT FLAG gIsWellAboveBG", "\n")
			cat("PROBES DISTRIBUTION across exp.cond. WITH (FLAG OK) gIsWellAboveBG = ", 
					flag, "\n")
			print(table(probeFLAG))
			cat("------------------------------------------------------", 
					"\n")
			flag = 1
			FLAG = ddFILT$other$gIsFound
			probeFLAG = apply(FLAG, 1, countFLAG, flag)
			rm(FLAG)
			cat("COUNT FLAG gIsFound", "\n")
			cat("PROBES DISTRIBUTION across exp.cond. WITH (FLAG OK) gIsFound = ", 
					flag, "\n")
			print(table(probeFLAG))
			cat("------------------------------------------------------", 
					"\n")
			flag = 0
			FLAG = ddFILT$other$gIsSaturated
			probeFLAG = apply(FLAG, 1, countFLAG, flag)
			rm(FLAG)
			cat("COUNT FLAG gIsSaturated", "\n")
			cat("PROBES DISTRIBUTION across exp.cond. WITH (FLAG OK) gIsSaturated = ", 
					flag, "\n")
			print(table(probeFLAG))
			cat("------------------------------------------------------", 
					"\n")
			flag = 0
			FLAG = ddFILT$other$gIsFeatPopnOL
			probeFLAG = apply(FLAG, 1, countFLAG, flag)
			rm(FLAG)
			cat("COUNT FLAG gIsFeatPopnOL", "\n")
			cat("PROBES DISTRIBUTION across exp.cond. WITH (FLAG OK) gIsFeatPopnOL = ", 
					flag, "\n")
			print(table(probeFLAG))
			cat("------------------------------------------------------", 
					"\n")
			flag = 0
			FLAG = ddFILT$other$gIsFeatNonUnifOL
			probeFLAG = apply(FLAG, 1, countFLAG, flag)
			rm(FLAG)
			cat("COUNT FLAG gIsFeatNonUnifOL", "\n")
			cat("PROBES DISTRIBUTION across exp.cond. WITH (FLAG OK) gIsFeatNonUnifOL = ", 
					flag, "\n")
			print(table(probeFLAG))
			cat("------------------------------------------------------", 
					"\n")
		}
	}
	if (!missing(makePLOT)) {
		if (makePLOT) {
			colorfill = "orange"
			maintitle = "FILTERED SIGNAL"
			X11()
			par(mfrow = c(1, 1))
			plotDensity(ddFILT$G, maintitle)
			X11()
			par(mfrow = c(1, 1))
			BoxPlot(ddFILT$G, maintitle, colorfill)
			maintitle = "FILTERED DATA - RLE "
			X11()
			par(mfrow = c(1, 1))
			RLE(ddFILT$G, maintitle, colorfill)
		}
	}
	return(ddFILT)
}

my.filter.control <- function (ddNORM, ManuelaGO, targets, annotation.package) 
{
	if (!is(ddNORM, "RGList")) {
		stop("'input' must be a RGList")
		if (is.null(dim(ddNORM)[1])) {
			stop("'input' is empty")
		}
	}
	if (missing(targets)) {
		stop("'targets' is missing ")
	}
	if ("GErep" %in% colnames(targets)) {
		g1 = targets$GErep
		g2 = rownames(targets)
	} else {
		stop("'targets' needs 'GErep' field")
	}
	cat("\n")
	cat("FILTERING BY ControlType FLAG", "\n")
	selSNR = which(ddNORM$genes$ControlType == 0)
	ddFILT = ddNORM[selSNR, ]
	#write.control.out(ddFILT, selSNR, annotation.package, ManuelaGO, 
	#		targets)
	cat("------------------------------------------------------", 
			"\n")
	cat(" PROBES BEFORE FILTERING: ", dim(ddNORM)[1], "\n")
	cat(" PROBES AFTER ControlType FILTERING: ", dim(ddFILT)[1], 
			"\n")
	cat(" RAW DATA WITHOUT CONTROLS OUT :", length(ddFILT$genes$ProbeName), 
			"\n")
	cat("------------------------------------------------------", 
			"\n")
	return(ddFILT)
}



#' run.signal.extraction - Run signal extraction based on f.sigprep
#' @param f.sigprep, f.targets, f.dd
#' @returnType RGListq
#' @return f.dd
#' @author daniel
#' @export
run.signal.extraction <- function(f.vGf,f.vGb,f.vBGmethod,f.vNORMmethod, f.targets, f.dd, f.annpkg, chip.validate=TRUE,AFE.2=FALSE,add.col=NULL) {
	
	print(paste("Run analysis: ",f.vGf,f.vGb,f.vBGmethod,f.vNORMmethod,sep=" "))
	
	###
	# Agilent Feature Extraction
	###
	
	if (is.null(f.dd)) {
		if (AFE.2) {
			f.dd=my.read.AgilentFE.2(f.targets, makePLOT=FALSE, Gf=f.vGf, Gb=f.vGb)
		} else {
			f.dd=my.read.AgilentFE(f.targets, makePLOT=FALSE, Gf=f.vGf, Gb=f.vGb)
		}
	}
	
	x11()
	plotDensity(log2(f.dd$G), "Density Plot Log2 Signal")
	dev.off()
	
	###
	# Pre-processing
	###
	
	## Change rownames based on the other columns
	
	con <- textConnection(as.character(f.targets$Subject))	# text connection
	subject <- read.table(con,sep='.')
	close(con)
	
	strReverse <- function(x) sapply(lapply(strsplit(x, NULL), rev), paste, collapse="") 
	
	aux <- gsub(".txt", "", as.character(f.targets$FileName)) 
	loc <- gsub("^(\\d+_\\d+).+","\\1",as.character(strReverse(aux)),perl=TRUE)
	con <- textConnection(as.character(loc))
	filename <- read.table(con, sep='_', fill=TRUE)
	close(con)
	rm(aux)
	
	if (!is.null(add.col)) {
		aux <- paste(f.targets[[add.col]],f.targets$Treatment, subject$V1, f.targets$Array, filename$V2, filename$V1,sep = '.')
	}
	else {
		aux <- paste(f.targets$Treatment, subject$V1, f.targets$Array, filename$V2, filename$V1,sep = '.')	# criando rownames ('condicao'.'nro da microdissec��o')
	}	
	rm(subject, filename)
	
	rownames(f.targets) <- aux
	rm(aux)
	
	# Add column 'Name' in f.targets 
	f.targets$Name <- row.names(f.targets)
	
	# Replicate f.targets in f.dd
	f.dd <- rep.targets(f.dd, f.targets)
	
	if (! file.exists(paste('quality/aQM_',f.vGf,'_',f.vGb,sep=""))) {
		## Creating the Expression Set Object (before background correction and normalization)
		esetRAW=my.build.eset(f.dd,f.targets,makePLOT=FALSE, annotation.package=f.annpkg, 
				unique.rownames=paste(f.dd$genes$ProbeName, seq(1:length(f.dd$genes$ProbeName)), sep="/"))
		dim(esetRAW)
		
		## Quality analysis (before background correction and normalization)
		#featureData(esetRAW)$X <- f.dd$other$Row[,1]
		#featureData(esetRAW)$Y <- f.dd$other$Col[,1]
		
		aqm.raw <- arrayQualityMetrics(	expressionset = esetRAW,
				outdir = paste('quality/aQM_',f.vGf,'_',f.vGb,sep=""),
				force = TRUE,
				do.logtransform = TRUE,
				intgroup=c("Treatment","Subject","Array")
		)
		dev.off()
		rm(esetRAW)
		rm(aqm.raw)
	}

	if (chip.validate) {
	
		all.annotation.data.chip.validate <-
				read.delim(file='/work/projects/microarray_pipeline/chip-validate/all-annotation-data-chip-validate',
						comment.char = "#",
						sep="\t",
						quote = ""
				)
		
		#colnames(all.annotation.data.chip.validate)
		#sum(summary(as.factor(all.annotation.data.chip.validate[,'C.V_STATUS'])))
		
		#ls(f.dd$other)
		f.dd$other$gIsValid <- c()
		cv.good <- subset(all.annotation.data.chip.validate, C.V_STATUS != 5 & is.na(C.V_STATUS)==FALSE)
		cv.good.names <- unique( cv.good$NAME )
		length(cv.good.names)
		
		#cv.bad <- subset(all.annotation.data.chip.validate, C.V_STATUS == 5 | is.na(C.V_STATUS)==TRUE)
		#dim(cv.bad)
		#cv.bad.names <- unique( cv.bad$NAME )
		#length(cv.bad.names)
		
		#summary(as.factor(all.annotation.data.chip.validate[,'C.V_STATUS']))
		
		f.dd$genes$gIsValid[ 1:length(f.dd$genes$ProbeName) ] <- 0
		f.dd$genes$gIsValid[ f.dd$genes$ProbeName %in% cv.good.names ] <- 1

		summary(as.factor(f.dd$genes$gIsValid))
		
		f.dd.filtered <- filter.isvalid(f.dd, targets=f.targets, ManuelaGO = FALSE, annotation.package=f.annpkg)
		summary(as.factor(f.dd.filtered$genes$gIsValid))
		
		f.dd <- f.dd.filtered
	}
	
	# Background Correction and Normalization
	f.ddNORM=my.BGandNorm(f.dd,BGmethod=f.vBGmethod,NORMmethod=f.vNORMmethod,
			foreground='Green',background='Green',
			offset=50,makePLOTpre=FALSE,makePLOTpost=FALSE,log.transform=TRUE)
	
	x11()
	plotDensity(f.ddNORM$G, "Density Plot Log2 Background/Normalized Signal")
	dev.off()
	
	if (! file.exists(paste('quality/aQM_',f.vGf,'_',f.vGb,'_',f.vBGmethod,'_',f.vNORMmethod,'_',chip.validate,sep=""))) {
		## Creating the Expression Set Object (after background correction and normalization)
		esetNORM=my.build.eset(f.ddNORM,f.targets,makePLOT=FALSE, annotation.package=f.annpkg, 
				unique.rownames=paste(f.ddNORM$genes$ProbeName, seq(1:length(f.ddNORM$genes$ProbeName)), sep="/"))
		dim(esetNORM)
		
		## Quality analysis (after background correction and normalization)
		#featureData(esetNORM)$X <- f.ddNORM$other$Row[,1]
		#featureData(esetNORM)$Y <- f.ddNORM$other$Col[,1]
		
		aqm.norm <- arrayQualityMetrics(	expressionset = esetNORM,
				outdir = paste('quality/aQM_',f.vGf,'_',f.vGb,'_',f.vBGmethod,'_',f.vNORMmethod,'_',chip.validate,sep=""),
				force = TRUE,
				do.logtransform = FALSE,
				intgroup=c("Treatment","Subject","Array")
		)
		dev.off()
		rm(esetNORM)
		rm(aqm.norm)
	}
	
	if (AFE.2) {
		## Filtering Probes
		f.ddFILT <- my.filter.probes(f.ddNORM,
				control=TRUE,                   # remoção de sondas controles
				wellaboveBG=TRUE,               # remoção de sondas não distinguíveis do background
				isfound=FALSE,                   # remoção de sondas não encontradas (1 - diferença:  sinal é 1.5 vezes o ruído do background local / 2 - o diâmetro do spot é ao menos 0.30 vezes o tamanho do diâmetro	nominal do spot)
				wellaboveNEG=TRUE,              # remoção de sondas abaixo dos valores dos controles negativos=  > Mean Negative Controls + 1.5*(Std. dev.Negative Controls).
				sat=TRUE,                       # remoção de sondas saturadas
				PopnOL=TRUE,                    # remoção de sondas outiliers
				NonUnifOL=TRUE,                 # remoção de sondas se pixels ruídos excedem determinado valor para um atributo uniforme
				nas=TRUE,                       # remove NAs
				makePLOT=FALSE,annotation.package=f.annpkg,flag.counts=FALSE,targets=f.ddNORM$targets)
	} else {
		## Filtering Probes
		f.ddFILT <- filter.probes(f.ddNORM,
				control=TRUE,                   # remoção de sondas controles
				wellaboveBG=TRUE,               # remoção de sondas não distinguíveis do background
				isfound=TRUE,                   # remoção de sondas não encontradas (1 - diferença:  sinal é 1.5 vezes o ruído do background local / 2 - o diâmetro do spot é ao menos 0.30 vezes o tamanho do diâmetro	nominal do spot)
				wellaboveNEG=TRUE,              # remoção de sondas abaixo dos valores dos controles negativos=  > Mean Negative Controls + 1.5*(Std. dev.Negative Controls).
				sat=TRUE,                       # remoção de sondas saturadas
				PopnOL=TRUE,                    # remoção de sondas outiliers
				NonUnifOL=TRUE,                 # remoção de sondas se pixels ruídos excedem determinado valor para um atributo uniforme
				nas=TRUE,                       # remove NAs
				limWellAbove=75,
				limISF=75,
				limNEG=75,
				limSAT=75,
				limPopnOL=75,
				limNonUnifOL=75,
				limNAS=100,
				makePLOT=FALSE,annotation.package=f.annpkg,flag.counts=TRUE,f.ddNORM$targets)
	}
	
	for (f in c('IsFeatPopnOL.txt',	'IsNOTWellAboveBG.txt',	'IsNOTFound.txt',
			'IsNOTWellAboveNEG.txt', 'RawDataNOCtrl.txt', 'RawDataNOCtrlWABKGandISF.txt',
			'IsSaturated.txt','IsFeatNonUnifOL.txt') ) {
		
		if (file.exists(f)) {
			system(paste("mv ",f," results/tables/ 2>&1"), intern=TRUE)		
		}
	}
	
	if (! file.exists(paste('quality/aQM_',f.vGf,'_',f.vGb,'_',f.vBGmethod,'_',f.vNORMmethod,'_',chip.validate,'_','filtered',sep=""))) {
		## Creating the Expression Set Object (after background correction and normalization)
		esetFILT=my.build.eset(f.ddFILT,f.targets,makePLOT=FALSE, annotation.package=f.annpkg, 
				unique.rownames=paste(f.ddFILT$genes$ProbeName, seq(1:length(f.ddFILT$genes$ProbeName)), sep="/"))
		dim(esetFILT)
		
		## Quality analysis (after background correction, normalization and filtering)
		#featureData(esetFILT)$X <- f.ddFILT$other$Row[,1]
		#featureData(esetFILT)$Y <- f.ddFILT$other$Col[,1]
		
		aqm.filt <- arrayQualityMetrics(	expressionset = esetFILT,
				outdir = paste('quality/aQM_',f.vGf,'_',f.vGb,'_',f.vBGmethod,'_',f.vNORMmethod,'_',chip.validate,'_','filtered',sep=""),
				force = TRUE,
				do.logtransform = FALSE,
				intgroup=c("Treatment","Subject","Array")
		)
		dev.off()
		rm(esetFILT)
		rm(aqm.filt)
	}
	
	#f.ddSUM <- summarize.tech.reps(f.ddFILT)
	#f.targets <- f.ddSUM$targets
	## Summarizing Probes
	#f.ddPROC=summarize.probe(f.ddSUM,makePLOT=FALSE,f.targets)
	
	## Summarizing Probes
	f.ddPROC=summarize.probe(f.ddFILT,makePLOT=FALSE,f.targets)
	
	## Creating the Expression Set Object (after background correction and normalization)
	esetPROC=build.eset(f.ddPROC,f.targets,makePLOT=FALSE, annotation.package=f.annpkg)
	dim(esetPROC)
	
	
	## Writing the Expression Set Object: Processeddata.txt
	write.eset(esetPROC,f.ddPROC,f.annpkg,f.ddPROC$targets)
	
	if (file.exists("ProcessedData.txt")) {
		system("mv ProcessedData.txt results/tables/ 2>&1", intern=TRUE)
	}
	
	if (! file.exists(paste('quality/aQM_',f.vGf,'_',f.vGb,'_',f.vBGmethod,'_',f.vNORMmethod,'_',chip.validate,'_','filtered','_','summarized',sep=""))) {
		## Quality analysis (after background correction, normalization, filtering and summarization)
		#featureData(esetPROC)$X <- f.ddPROC$other$Row[,1]
		#featureData(esetPROC)$Y <- f.ddPROC$other$Col[,1]
		
		aqm.processed <- arrayQualityMetrics(	expressionset = esetPROC,
				outdir = paste('quality/aQM_',f.vGf,'_',f.vGb,'_',f.vBGmethod,'_',f.vNORMmethod,'_',chip.validate,'_','filtered','_','summarized',sep=""),
				force = TRUE,
				do.logtransform = FALSE,
				intgroup=c("Treatment","Subject","Array")
		)
		dev.off()
		rm(aqm.processed)
	}		

	return(list(dd=f.dd,ddNORM=f.ddNORM,ddFILT=f.ddFILT,ddPROC=f.ddPROC))
}


#' normalize75thPercentileAbsValues - Normalization by 75th percentile
#' @param raw matrix
#' @returnType matrix
#' @return normalized matrix
#' @author daniel
#' @export
normalize75thPercentileAbsValues <- function (x) {
	
	percentile75 <- function(w, ...) {
		as.numeric(quantile(w, probs=0.75, ...))
	}
	
	narrays <- NCOL(x)
	if (narrays == 1) 
		return(x)
	cmed <- log(apply(abs(x), 2, percentile75, na.rm = TRUE))
	cmed <- exp(cmed - mean(cmed))
	t(t(x)/cmed)
}

#' my.normalizeBetweenArrays - Added normalization by 75th percentile (scale75thp)
#' @param matrix
#' @returnType matrix
#' @return normalized matrix
#' @author daniel
#' @export
my.normalizeBetweenArrays <- function (object, method = NULL, targets = NULL, ...) {
	if (is.null(method)) {
		if (is(object, "matrix")) {
			method = "quantile"
		}
		else if (is(object, "EListRaw")) {
			method = "quantile"
		}
		else {
			method = "Aquantile"
		}
	}
	choices <- c("none", "scale", "scale75thp", "quantile", "Aquantile", "Gquantile", 
			"Rquantile", "Tquantile", "vsn")
	method <- match.arg(method, choices)
	if (method == "vsn") 
		stop("vsn method no longer supported. Please use normalizeVSN instead.")
	if (is(object, "matrix")) {
		if (!(method %in% c("none", "scale", "scale75thp", "quantile", "vsn"))) 
			stop("method not applicable to matrix objects")
		return(switch(method, none = object, scale = normalizeMedianAbsValues(object), 
						scale75thp = normalize75thPercentileAbsValues(object),
						quantile = normalizeQuantiles(object, ...), vsn = exprs(vsnMatrix(x = object, 
										...))))
	}
	if (is(object, "EListRaw")) {
		object$E <- log2(Recall(object$E, method = method, ...))
		object <- new("EList", unclass(object))
		return(object)
	}
	if (is(object, "RGList")) 
		object <- MA.RG(object)
	if (is.null(object$M) || is.null(object$A)) 
		stop("object doesn't appear to be RGList or MAList object")
	switch(method, scale = {
				object$M <- normalizeMedianAbsValues(object$M)
				object$A <- normalizeMedianAbsValues(object$A)
			}, scale75thp = {
				object$M <- normalize75thPercentileAbsValues(object$M)
				object$A <- normalize75thPercentileAbsValues(object$A)
			}, quantile = {
				narrays <- NCOL(object$M)
				Z <- normalizeQuantiles(cbind(object$A - object$M/2, 
								object$A + object$M/2), ...)
				G <- Z[, 1:narrays]
				R <- Z[, narrays + (1:narrays)]
				object$M <- R - G
				object$A <- (R + G)/2
			}, Aquantile = {
				object$A <- normalizeQuantiles(object$A, ...)
			}, Gquantile = {
				G <- object$A - object$M/2
				E <- normalizeQuantiles(G, ...) - G
				object$A <- object$A + E
			}, Rquantile = {
				R <- object$A + object$M/2
				E <- normalizeQuantiles(R, ...) - R
				object$A <- object$A + E
			}, Tquantile = {
				narrays <- NCOL(object$M)
				if (NCOL(targets) > 2) targets <- targets[, c("Cy3", 
									"Cy5")]
				targets <- as.vector(targets)
				Z <- cbind(object$A - object$M/2, object$A + object$M/2)
				for (u in unique(targets)) {
					j <- targets == u
					Z[, j] <- normalizeQuantiles(Z[, j], ...)
				}
				G <- Z[, 1:narrays]
				R <- Z[, narrays + (1:narrays)]
				object$M <- R - G
				object$A <- (R + G)/2
			})
	object
}


# Added parameters:  Rf="gProcessedSignal", Gf="gMeanSignal", Rb="gBGMedianSignal", Gb="gBGUsed"
# Capture columns Row and Col
my.read.AgilentFE <- function (targets, makePLOT,	Rf="gProcessedSignal", Gf="gMeanSignal", 
													Rb="gBGMedianSignal", Gb="gBGUsed") {
	if (!is(targets, "data.frame")) {
		stop("'targets' must be a data.frame")
	}
	ddaux = read.maimages(files = targets$FileName, source = "agilent", 
			other.columns = list(IsFound = "gIsFound", IsWellAboveBG = "gIsWellAboveBG", 
					IsSaturated = "gIsSaturated", IsFeatNonUnifOF = "gIsFeatNonUnifOL", 
					IsFeatPopnOL = "gIsFeatPopnOL", ChrCoord = "chr_coord", Row = "Row", Col = "Col"), 
			columns = list(Rf = Rf, Gf = Gf, 
					Rb = Rb, Gb = Gb), verbose = T, 
			sep = "\t", quote = "")
	if (length(ddaux$genes$Sequence) == 0) {
		cat(" INPUT DATA DOES NOT CONTAIN - Sequence and chr_coord", 
				"\n")
		cat(" SCANN THE DATA USING AFE 9.5.3.1", "\n")
		stop(" the script will stop now")
	}
	dd = new("RGList")
	dd$R = ddaux$R
	dd$G = ddaux$G
	dd$Rb = ddaux$Rb
	dd$Gb = ddaux$Gb
	dd$targets = ddaux$targets
	dd$genes = ddaux$genes[, c(4, 6, 7, 8, 9, 10)]
	dd$other = ddaux$other
	rm(ddaux)
	cat("", "\n")
	cat("  RGList:", "\n")
	cat("\tdd$R:\t'",Rf,"' ", "\n")
	cat("\tdd$G:\t'",Gf,"' ", "\n")
	cat("\tdd$Rb:\t'",Rb,"' ", "\n")
	cat("\tdd$Gb:\t'",Gb,"' ", "\n")
	cat("", "\n")
	if (!missing(makePLOT)) {
		if (makePLOT) {
			MMM = log2(dd$R)
			maintitle = Rf
			colorfill = "red"
			X11()
			par(mfrow = c(1, 1))
			BoxPlot(MMM, maintitle, colorfill)
			X11()
			par(mfrow = c(1, 1))
			plotDensity(MMM, maintitle)
			MMM = log2(dd$G)
			maintitle = Gf
			colorfill = "green"
			X11()
			par(mfrow = c(1, 1))
			BoxPlot(MMM, maintitle, colorfill)
			X11()
			par(mfrow = c(1, 1))
			plotDensity(MMM, maintitle)
			MMM = log2(dd$Gb)
			maintitle = Gb
			colorfill = "orange"
			X11()
			par(mfrow = c(1, 1))
			BoxPlot(MMM, maintitle, colorfill)
			MMM = log2(dd$Rb)
			maintitle = Rb
			colorfill = "blue"
			X11()
			par(mfrow = c(1, 1))
			BoxPlot(MMM, maintitle, colorfill)
		}
	}
	return(dd)
}

# Added parameters:  Rf="gProcessedSignal", Gf="gMeanSignal", Rb="gBGMedianSignal", Gb="gBGUsed"
# Capture columns Row and Col
my.read.AgilentFE.2 <- function (targets, makePLOT,	Rf="gProcessedSignal", Gf="gMeanSignal", 
		Rb="gBGMedianSignal", Gb="gBGUsed") {
	if (!is(targets, "data.frame")) {
		stop("'targets' must be a data.frame")
	}
	ddaux = read.maimages(files = targets$FileName, source = "agilent", 
			other.columns = list(IsWellAboveBG = "gIsWellAboveBG", 
					IsSaturated = "gIsSaturated", IsFeatNonUnifOF = "gIsFeatNonUnifOL", 
					IsFeatPopnOL = "gIsFeatPopnOL", Row = "Row", Col = "Col"), 
			columns = list(Rf = Rf, Gf = Gf, 
					Rb = Rb, Gb = Gb), verbose = T, 
			sep = "\t", quote = "")
#	if (length(ddaux$genes$Sequence) == 0) {
#		cat(" INPUT DATA DOES NOT CONTAIN - Sequence and chr_coord", 
#				"\n")
#		cat(" SCANN THE DATA USING AFE 9.5.3.1", "\n")
#		stop(" the script will stop now")
#	}
	dd = new("RGList")
	dd$R = ddaux$R
	dd$G = ddaux$G
	dd$Rb = ddaux$Rb
	dd$Gb = ddaux$Gb
	dd$targets = ddaux$targets
	dd$genes = ddaux$genes
	dd$other = ddaux$other
	rm(ddaux)
	cat("", "\n")
	cat("  RGList:", "\n")
	cat("\tdd$R:\t'",Rf,"' ", "\n")
	cat("\tdd$G:\t'",Gf,"' ", "\n")
	cat("\tdd$Rb:\t'",Rb,"' ", "\n")
	cat("\tdd$Gb:\t'",Gb,"' ", "\n")
	cat("", "\n")
	if (!missing(makePLOT)) {
		if (makePLOT) {
			MMM = log2(dd$R)
			maintitle = Rf
			colorfill = "red"
			X11()
			par(mfrow = c(1, 1))
			BoxPlot(MMM, maintitle, colorfill)
			X11()
			par(mfrow = c(1, 1))
			plotDensity(MMM, maintitle)
			MMM = log2(dd$G)
			maintitle = Gf
			colorfill = "green"
			X11()
			par(mfrow = c(1, 1))
			BoxPlot(MMM, maintitle, colorfill)
			X11()
			par(mfrow = c(1, 1))
			plotDensity(MMM, maintitle)
			MMM = log2(dd$Gb)
			maintitle = Gb
			colorfill = "orange"
			X11()
			par(mfrow = c(1, 1))
			BoxPlot(MMM, maintitle, colorfill)
			MMM = log2(dd$Rb)
			maintitle = Rb
			colorfill = "blue"
			X11()
			par(mfrow = c(1, 1))
			BoxPlot(MMM, maintitle, colorfill)
		}
	}
	return(dd)
}






# Added options:  'Tquantile', 'scale75thp', 'Tscale75thp'
# Added parameter: log.transform (TRUE/FALSE)
# Changed options and parameters for background (Green/Red) and foreground (Green/Red)
my.BGandNorm <- function (	RGlist, BGmethod = NULL, NORMmethod = NULL, foreground = NULL, 
							background = NULL, offset = NULL, makePLOTpre = NULL, makePLOTpost = NULL, 
							log.transform = FALSE) 
{
	require(vsn)
	cat("BACKGROUND CORRECTION AND NORMALIZATION ", "\n")
	cat("\n")
	if (!is(RGlist, "RGList")) {
		stop("'input' must be a RGList")
		if (is.null(dim(RGlist)[1])) {
			stop("'input' is empty")
		}
	}
	if (missing(BGmethod)) {
		BGmethod = "none"
	}
	else {
		if (BGmethod != "none" && BGmethod != "half" && BGmethod != 
				"normexp") {
			stop("BGmethod should be one of 'none', 'half','normexp' ")
		}
	}
	if (missing(NORMmethod)) {
		NORMmethod = "none"
	}
	else {
		if (NORMmethod != "none" && NORMmethod != "quantile" && 
				NORMmethod != "vsn" && NORMmethod != "Tquantile" && NORMmethod != "scale75thp" && NORMmethod != "Tscale75thp") {
			stop("NORMmethod should be one of 'none', 'quantile', 'vsn', 'Tquantile', 'scale75thp' or 'Tscale75thp' ")
		}
	}
	
	if (missing(foreground)) {
		cat("'foreground signal is missing => foreground used is 'Green'", 
				"\n")
		cat("\tNORMmethod is set to 'quantile' ", "\n")
		NORMmethod = "quantile"
	}
	else {
		if (foreground != "Red" && foreground != 
				"Green") {
			cat("foreground: ", foreground, "\n")
			stop("'foreground' must be either 'Red' or 'Green'")
		}
		else {
			if (foreground == "Red") {
				RGlist$G = RGlist$R
				cat("\tforeground: Red", "\n")
			}
			else {
				cat("\tforeground: Green", "\n")
			}
		}
	}
	if (missing(background)) {
		cat("'background signal is missing => BGmethod is set to 'none'")
		BGmethod = "none"
	}
	else {
		if (background != "Red" && background != "Green") {
			cat("background: ", background, "\n")
			stop("'background' must be either 'Red' or 'Green'")
		}
		else {
			if (background == "Red") {
				RGlist$Gb = RGlist$Rb
				cat("\tbackground: Red", "\n")
			}
			else {
				cat("\tbackground: Green", "\n")
			}
		}
	}
	if (missing(offset)) {
		offset = 0
	}
	p = c(1:dim(RGlist)[2])
	controlNA = length(which(is.na(RGlist$G[, p]) == TRUE))
	if (controlNA > 0) {
		cat("# NAs in RGlist$G - ", controlNA, "\n")
		stop("program will stop now")
	}
	if (!missing(makePLOTpre)) {
		if (makePLOTpre) {
			ddaux = RGlist
			ddaux$G = log2(ddaux$G)
			colorfill = "blue"
			maintitle = "RAW SIGNAL"
			MMM = ddaux$G
			X11()
			par(mfrow = c(1, 1))
			plotDensity(MMM, maintitle)
			X11()
			par(mfrow = c(1, 1))
			boxplotNegCtrl(ddaux, Log2 = TRUE)
			X11()
			par(mfrow = c(1, 1), ask = T)
			MVAplotMEDctrl(ddaux, maintitle)
			maintitle = "RAW DATA - RLE "
			X11()
			par(mfrow = c(1, 1))
			RLE(MMM, maintitle, colorfill)
			rm(MMM)
			rm(ddaux)
		}
	}
	if (BGmethod == "half") {
		ddBG = backgroundCorrect(RGlist, method = BGmethod, offset = offset)
	}
	if (BGmethod == "normexp") {
		ddBG = RGlist
		MMM = RGlist$G - RGlist$Gb
		MMM = backgroundCorrect(MMM, method = BGmethod, offset = offset)
		ddBG$G = MMM
		rm(MMM)
	}
	if (BGmethod == "none") {
		ddBG = RGlist
	}
	if (NORMmethod != "none") {
		ddNORM = RGlist
		
			if (NORMmethod == "Tquantile") {
				
				if (!is.null(ddNORM$targets$GErep)) {
					exprsNORM <- ddBG$G
					
					for (u in unique(ddNORM$targets$GErep)) {
						j <- ddNORM$targets$GErep == u
						if (length(j[j==TRUE]) > 1) {
							print(paste("Applying", NORMmethod,"normalization between arrays",paste(colnames(exprsNORM[,j]),collapse=","), sep=" "))
							exprsNORM[,j] <- my.normalizeBetweenArrays(ddBG$G[,j], method = 'quantile')
						}
					}
					
				}
				else {
					stop("GErep in targets of RGlist is null")
				}
				
			}
			else if (NORMmethod == "Tscale75thp") {
				if (!is.null(ddNORM$targets$GErep)) {
					exprsNORM <- ddBG$G
					
					for (u in unique(ddNORM$targets$GErep)) {
						j <- ddNORM$targets$GErep == u
						if (length(j[j==TRUE]) > 1) {
							print(paste("Applying", NORMmethod,"normalization between arrays",paste(colnames(exprsNORM[,j]),collapse=","), sep=" "))
							exprsNORM[,j] <- my.normalizeBetweenArrays(ddBG$G[,j], method = 'scale75thp')
						}
					}
				}
				else {
					stop("GErep in targets of RGlist is null")
				}			
			}
			else if (NORMmethod == "vsn") {
				exprsNORM = normalizeVSN(ddBG$G)
			}
			else {
				exprsNORM = my.normalizeBetweenArrays(ddBG$G, method = NORMmethod)
			}
		
#		if (	(NORMmethod == "quantile")||
#				(NORMmethod == "Tquantile")) {
#			ddNORM$G = exprsNORM
#			ddNORM$G = log2(exprsNORM)
#		}
#		else {
			ddNORM$G = exprsNORM
#		}
		rm(exprsNORM)
		
	}
	else {
		ddNORM = ddBG
#		ddNORM$G = log2(ddBG$G)
	}
	cat("\n")
	cat("\tBGmethod:\t", BGmethod, "\n")
	cat("\tNORMmethod:\t", NORMmethod, "\n")
	
	if (log.transform) {
		if (	(NORMmethod == "quantile")||
				(NORMmethod == "Tquantile")
				) {
		
					ddNORM$G = log2(ddNORM$G)
					cat("\tOUTPUT in log-2 scale", "\n")					
		}
		else {
			cat("\tNOT PERFORMING log-2 transformation")
		}
	}
	else {
		cat("\tOUTPUT not-scaled", "\n") 
	}
	
	cat("\n")
	cat("------------------------------------------------------", 
	"\n")

	ddNORM$R = ddNORM$G
	ddNORM$Gb = ddNORM$G
	ddNORM$Rb = ddNORM$G

	if (!missing(makePLOTpost)) {
		if (makePLOTpost) {
			colorfill = "red"
			maintitle = "NORMALIZED SIGNAL"
			MMM = ddNORM$G
			X11()
			par(mfrow = c(1, 1))
			plotDensity(MMM, maintitle)
			X11()
			par(mfrow = c(1, 1))
			boxplotNegCtrl(ddNORM, Log2 = TRUE)
			X11()
			par(mfrow = c(1, 1), ask = T)
			MVAplotMEDctrl(ddNORM, maintitle)
			maintitle = "NORMALIZED DATA - RLE "
			X11()
			par(mfrow = c(1, 1))
			RLE(MMM, maintitle, colorfill)
			rm(MMM)
		}
	}
	return(ddNORM)
}

# Added parameter:  unique.rownames
# Enable the creation of the Expression Set Object to use in aQM
my.build.eset <- function (RGlist, targets, makePLOT, annotation.package, unique.rownames) 
{
	if (!is(RGlist, "RGList")) {
		stop("'input' must be a RGList")
		if (is.null(dim(RGlist)[1])) {
			stop("'input' is empty")
		}
	}
	if (missing(targets)) {
		stop("'targets' is missing ")
	}
	if ("GErep" %in% colnames(targets)) {
		GErep = targets$GErep
		nGE = sum(as.integer(table(table(GErep))))
		g1 = targets$GErep
		g2 = rownames(targets)
	}
	else {
		stop("'targets' needs 'GErep' field")
	}
	goON = all(rownames(targets) == colnames(RGlist$G))
	if (!goON) {
		stop("rownames in pData(targets) different from colnames RGlist$G", 
				"\n")
	}
	phenoData = new("AnnotatedDataFrame", data = targets)
	TMP = RGlist$G
	
	if (!missing(unique.rownames)) {
		rownames(TMP) = unique.rownames
	} 
	else {
		rownames(TMP) = RGlist$genes$ProbeName	
	}
		
	colnames(TMP) = g2
	nGEN = dim(TMP)[1]
	nARR = dim(TMP)[2]
	esetPROC = new("ExpressionSet", exprs = TMP, phenoData = phenoData, 
			annotation = annotation.package)
	if (!missing(makePLOT)) {
		if (makePLOT) {
			X11()
			par(mfrow = c(1, 1))
			hierclus(exprs(esetPROC), GErep, methdis = "euclidean", 
					methclu = "complete", sel = FALSE, 100)
			X11()
			par(mfrow = c(1, 1))
			size = 100
			maintitle = "100 High Var"
			HeatMap(exprs(esetPROC), size, maintitle)
			X11()
			par(mfrow = c(1, 1))
			PCAplot(esetPROC, targets)
		}
	}
	return(esetPROC)
}


# Enable the creation of the Expression Set Object to use in aQM from matrx
my.build.eset.matrix <- function (m, targets, annotation.package, unique.rownames) 
{
	if (!is(m, "matrix")) {
		stop("'input' must be a matrix")
		if (is.null(dim(m)[1])) {
			stop("'input' is empty")
		}
	}
	if (missing(targets)) {
		stop("'targets' is missing ")
	}
	if ("GErep" %in% colnames(targets)) {
		GErep = targets$GErep
		nGE = sum(as.integer(table(table(GErep))))
		g1 = targets$GErep
		g2 = rownames(targets)
	}
	else {
		stop("'targets' needs 'GErep' field")
	}
	goON = all(rownames(targets) == colnames(m))
	if (!goON) {
		stop("rownames in pData(targets) different from matrix colnames", 
				"\n")
	}
	phenoData = new("AnnotatedDataFrame", data = targets)
	TMP = m
	
	if (!missing(unique.rownames)) {
		rownames(TMP) = unique.rownames
	} 
	
	colnames(TMP) = g2
	nGEN = dim(TMP)[1]
	nARR = dim(TMP)[2]
	esetPROC = new("ExpressionSet", exprs = TMP, phenoData = phenoData, 
			annotation = annotation.package)
	return(esetPROC)
}
