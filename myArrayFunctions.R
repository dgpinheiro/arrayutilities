# LGMB   Laboratory of Genetics and Molecular Biology
# BiT    Bioinformatics Team
#
# Project...: myArray
# Date......: 22/07/2011
# Author(s).: daniel, matheus
###############################################################################



#' saveArray - Save all objects in a archive, named with pattern name.
#' @param dir Directory in which the objects will be saved
#' @returnType NULL
#' @return none
#' @author matheus, daniel
#' @export
saveArray <- function(dir,prefix='R')
{
		
	if( !file.exists(dir) )
		stop("No such file or directory.")
	
	date = format(Sys.Date(), "%d_%m_%Y")
	
	i = 0
	
	file = paste(prefix, date, i, sep="_")
	file = paste(file, "RData", sep=".")
		
	while( file.exists(paste(dir,file,sep="/")) )
	{
		i = i + 1
		file = paste(prefix, date, i, sep="_")
		file = paste(file, "RData", sep=".")
	}
		
	save(list = ls(all=TRUE,envir=globalenv()), file = paste(dir,file, sep="/"))
	print(paste("GlobalEnv saved in '",paste(dir,file, sep="/"),"'",sep=""))
	
}


#' loadArray - Load all objects saved in the most recent archive
#' @param dir Directory in which the objects will be loaded
#' @returnType NULL
#' @return none
#' @author matheus
#' @export
loadArray <- function(dir, prefix='R')
{
	if( !file.exists(dir) )
		stop("No such file or directory.")
	files = list.files(path = dir, pattern = paste(prefix,"[^.]*\\.RData",sep=""))
	
	info <- file.info(paste(dir,files,sep="/"))	# obtem informacoes dos arquivos
	file <- ""
	file <- dimnames(info[which(info$mtime==max(info$mtime)),])[[1]]	# encontra o arquivo mais recente
	
	print(paste("Loading archive ",file," ...",sep=""))
	
	load(file=file,envir = globalenv())		
}

#' create.structure - Create structure of directories
#' @param NULL
#' @returnType NULL
#' @return none
#' @author daniel
#' @export
create.structure <- function()
{
	dirs <- c('quality','saves','results/tables','results/plots', 'results/plots/probes', 'results/tables/probes'
			
			)
	
	for (d in dirs) {
		if (!file.exists(d)) {
			print(paste("Creating directory:",d))
			system(paste("mkdir -p ",d," 2>&1"), intern=TRUE)
		}
	}
}

#' rep.targets - Replicate changes in rownames and colnames of RGlist based on targets
#' @param RGlist, targets
#' @returnType RGlist
#' @return RGlist
#' @author daniel
#' @export
rep.targets <- function(lst, rtargets) {
	
	if (as.character(class(lst)) == 'uRNAList') {
		colnames(lst$meanS) <- rtargets$Name
		colnames(lst$procS) <- rtargets$Name
		colnames(lst$TGS) <- rtargets$Name
		colnames(lst$TPS) <- rtargets$Name
		colnames(lst$other$gBGMedianSignal) <- rtargets$Name
		colnames(lst$other$gBGUsed) <- rtargets$Name
		colnames(lst$other$gIsFeatNonUnifOL) <- rtargets$Name 
		colnames(lst$other$gIsFeatPopnOL) <- rtargets$Name
		colnames(lst$other$gIsGeneDetected) <- rtargets$Name
		colnames(lst$other$gIsSaturated) <- rtargets$Name
	}
	else {
		lst$targets$FileName <- rtargets$Name
		colnames(lst$R) <- rtargets$Name
		colnames(lst$G) <- rtargets$Name
		colnames(lst$Rb) <- rtargets$Name
		colnames(lst$Gb) <- rtargets$Name		
		colnames(lst$other$gIsWellAboveBG) <- rtargets$Name
		if (!is.null(lst$other$gIsFound))
			colnames(lst$other$gIsFound) <- rtargets$Name
		
		colnames(lst$other$gIsSaturated) <- rtargets$Name
		colnames(lst$other$gIsFeatPopnOL) <- rtargets$Name
		colnames(lst$other$gIsFeatNonUnifOL) <- rtargets$Name
		
			if (!is.null(lst$other$chr_coord))
		colnames(lst$other$chr_coord) <- rtargets$Name
	
		colnames(lst$other$Col) <- rtargets$Name
		colnames(lst$other$Row) <- rtargets$Name
	}
	
	lst$targets <- rtargets
	return(lst)
}

#' summarize.tech.reps - Summarize technical replicates into RGlist
#' @param RGlist
#' @returnType RGlist
#' @return RGlist
#' @author daniel
#' @export
summarize.tech.reps <- function (RGlist) {
	
	RGlist$targets$SubjectTreatment <- paste(RGlist$targets$Subject,RGlist$targets$Treatment,sep="/")
	
	for (st in levels(as.factor(RGlist$targets$SubjectTreatment)) ) {
		
		nr <- dim(RGlist$targets[which(RGlist$targets$SubjectTreatment == st),])[1]
		#print(paste(st,nr,sep="="))
		
		if (nr > 1) {
			s.names <- RGlist$targets[which(RGlist$targets$SubjectTreatment == st),'Name']
			
			# [1] "FileName"         "Treatment"        "GErep"            "Subject"         
			# [5] "Array"            "Name"             "SubjectTreatment"		
			rtargets <- c()
			for (col in colnames(RGlist$targets)) {
				rtargets[col] <- paste(unique(RGlist$targets[s.names, col]), collapse='|')
			}
			
			s.new.name <- s.names[1]	
			for (i in 2:nr) {
				s.new.name <- paste(s.new.name, 
						paste(strsplit( s.names[i] , "\\.")[[1]][c(4,5)], collapse='.'), sep='/')
			}
			
			rtargets['Name'] <- s.new.name
			
			print(paste("Deriving another sample (",s.new.name,") from: ",paste(s.names,collapse=";"), sep=""))
			
			for (k in c('G','Gb')) {				
				RGlist[[k]] <- cbind(RGlist[[k]], apply(RGlist[[k]][,s.names], 1, median))
				colnames(RGlist[[k]])[length(colnames(RGlist[[k]]))] <- s.new.name
				RGlist[[k]] <- subset(RGlist[[k]], select=!(colnames(RGlist[[k]]) %in% s.names))				
			}
			
			RGlist[['R']] <-  RGlist[['G']]
			RGlist[['Rb']] <-  RGlist[['Gb']]
			
			for (k in c('gIsWellAboveBG', 'gIsFound', 'gIsSaturated',	'gIsFeatPopnOL', 'gIsFeatNonUnifOL'	)) {
				RGlist$other[[k]] <- cbind(RGlist$other[[k]], apply(RGlist$other[[k]][,s.names], 1, median))
				colnames(RGlist$other[[k]])[length(colnames(RGlist$other[[k]]))] <- s.new.name
				RGlist$other[[k]] <- subset(RGlist$other[[k]], select=!(colnames(RGlist$other[[k]]) %in% s.names))				
			}	
			
			for (k in c('chr_coord','Col','Row') ) {
				RGlist$other[[k]] <- cbind(RGlist$other[[k]], RGlist$other[[k]][,s.names[1]])
				colnames(RGlist$other[[k]])[length(colnames(RGlist$other[[k]]))] <- s.new.name
				RGlist$other[[k]] <- subset(RGlist$other[[k]], select=!(colnames(RGlist$other[[k]]) %in% s.names))
			}
			
			RGlist$targets <- RGlist$targets[!(rownames(RGlist$targets) %in% s.names),]		
						
			#print(as.character(rtargets))
			
			RGlist$targets <- rbind(RGlist$targets, as.character(rtargets))
			
			rownames(RGlist$targets)[length(rownames(RGlist$targets))] <- rtargets['Name']
		}
	}
	rep.targets(RGlist, RGlist$targets)
	return(RGlist)
}

#' summarize.gene.reps - Summarize gene replicates into a data.frame
#' @param data.frame
#' @param samples
#' @returnType data.frame
#' @return data.frame summarized
#' @author daniel
#' @export
summarize.gene.reps <- function(temp.s, samps) {
	
	temp.s$Symbol <- as.factor(temp.s$Symbol)
	unique.sym <- as.character(levels(temp.s$Symbol))
	
	nas <- which(is.na(temp.s$Symbol))
	
	lnas <- length(nas)
	
	temp.sym <- as.data.frame(matrix(rep(0,length(samps)*(length(setdiff(unique.sym,NA))+lnas)),			
			ncol=length(samps),
			nrow=(length(setdiff(unique.sym,NA))+lnas)))

	colnames(temp.sym) <- samps
					
	temp.sym$Symbol <- rep(NA,dim(temp.sym)[1])
		
	rownames(temp.sym) <- c(setdiff(unique.sym,NA), temp.s$ProbeName[nas])
	
	for (symbol in unique.sym ) {
		#print(symbol)			
		if (is.na(symbol)) {
			
			greps <- nas
			
			for (ena in greps) {
				tmp <- temp.s[ena, samps]	
				temp.sym[temp.s$ProbeName[ena], samps] <- as.numeric(tmp)
				temp.sym[temp.s$ProbeName[ena], 'Symbol'] <- NA
			}
		} else {
			
			greps <- which(temp.s$Symbol==symbol)
			
			if (length(greps)>1) {
				#print(paste(symbol, paste(greps, collapse=";")))
				tmp <- apply(temp.s[greps, samps], 2, median)
				
			} else {
				tmp <- temp.s[greps, samps]	
			}
			temp.sym[symbol, samps] <- as.numeric(tmp)
			temp.sym[symbol, 'Symbol'] <- symbol
		}
		
	}
	
	return(temp.sym)
}


#' create.contrast.matrix - Create contrast matrix
#' @param targets
#' @returnType data.frame
#' @return data.frame
#' @author daniel
#' @export
create.contrast.matrix <- function(targets) {
	addvar<-targets
	
	tmp <- addvar

	unique.factor <- c()
	for (t in rownames(tmp)) {
		z <- 1
		for (f in (unlist(strsplit(as.character(tmp[t,'Factor']),'.',fixed=TRUE)))) {
			fn <- as.character(paste(c('factor',z),collapse='.'))
			tmp[t,fn] <- f
			z<-z+1
			unique.factor <- unique(c(unique.factor,fn))
		}
	}
	
	for (fn in unique.factor) {
		unique.treatment <- unique(tmp[[fn]])
		length.unique.treatment <-  length(unique.treatment)	
		print(fn,quote=FALSE)
		for (i in c(1:length.unique.treatment)) {
			if (i < length.unique.treatment) {
				for (j in c((i+1):length.unique.treatment) ) {
					cmp.name <- paste(unique.treatment[i],'x',unique.treatment[j],sep="")
					
					
					print(paste("", cmp.name, collapse="\t"),quote=FALSE)
					addvar[[cmp.name]] <- rep(0,length(addvar$Name))
					addvar[which(tmp[[fn]] == unique.treatment[i]),cmp.name] <- 1
					addvar[which(tmp[[fn]] == unique.treatment[j]),cmp.name] <- 2
					
					# check if there is a paired comparison
					# unique.treatment[i] x unique.treatment[j]
					
					paired.sbjct.tmp <-
					apply(
					table(
						data.frame(
							Factor=as.factor(c(
								tmp[[fn]][which(tmp[[fn]] == unique.treatment[i])],
								tmp[[fn]][which(tmp[[fn]] == unique.treatment[j])]
								)
							),
							Subject=as.factor(c(
								as.character(addvar[which(tmp[[fn]] == unique.treatment[i]),'Subject']),
								as.character(addvar[which(tmp[[fn]] == unique.treatment[j]),'Subject']))
							)
						)
					)>=1, 2, function(x) { x[1] & x[2] } )
					
					paired.sbjct <- names(paired.sbjct.tmp[paired.sbjct.tmp==TRUE])
					#unpaired.sbjct <- names(paired.sbjct.tmp[paired.sbjct.tmp==FALSE])
					
					if (length(paired.sbjct)>1) {
						addvar[[paste('Paired:',cmp.name,sep="")]] <- rep(0,length(addvar$Name))
						p.id <- 1
						for (s in paired.sbjct) {
							addvar[ addvar$Subject == s & tmp[[fn]] == unique.treatment[i], paste('Paired:',cmp.name,sep="") ] <- p.id
							addvar[ addvar$Subject == s & tmp[[fn]] == unique.treatment[j], paste('Paired:',cmp.name,sep="") ] <- p.id
							p.id <- p.id+1
						}

					}
				}
			}
		}
	}
	
	return(addvar)
}


#welch.df(	as.numeric(dat['CD86/A_24_P131589',grep('^(EM-PRE)',colnames(dat))]),
#		as.numeric(dat['CD86/A_24_P131589',grep('^(EM-POS)',colnames(dat))])
#)


#' welch.df - Welch (or Satterthwaite) approximation to the degrees of freedom
#' @param ...
#' @returnType ...
#' @return ...
#' @author daniel
#' @export
welch.df <- function(x,y) {
	n<-length(x)
	m<-length(y)
	return(
			(	(	(sd(x)^2/n)	+
							(	 sd(y)^2/m)
							)^2)/
					(	(	(sd(x)^4/((n^2)*(n-1)))	+
							(	 sd(y)^4/((m^2)*(m-1)))
							)))
}

#' my.html3D - html3D of made4 using jmol
#' @param ...
#' @returnType ...
#' @return ...
#' @author daniel
#' @export
my.html3D <- function (df, classvec = NULL, writepdb = FALSE, filenamebase = "output", 
		writehtml = FALSE, title = NULL, scaled = TRUE, xyz.axes = c(1:3), jmolpath= ".", 
		...) 
{
	
	require(made4)
	
	if (ncol(df) < 3) 
		stop("need 3 columns to create 3D plot")
	df <- df[, xyz.axes]
	btt <- function(x) {
		for (i in c(1e-05, 1e-04, 0.001, 0.01, 0.1, 1, 10, 100, 
				1000, 10000, 1e+05)) {
			if (min(x) >= min(i * (-1)) && max(x) <= max(i)) 
				return(x * 10/i)
		}
	}
	if (!is.null(classvec)) 
		classvec = made4:::checkfac(classvec)
	addaxes <- function() {
		XX = matrix(c(c(2, 4, -2, -4), c(0, 0, 0, 0), c(0, 0, 
								0, 0)), ncol = 3, dimnames = list(NULL, LETTERS[24:26]))
		axescord = rbind(c(0, 0, 0), XX, XX[, c(3, 1, 2)], XX[, 
						c(3, 2, 1)])
		axescord = cbind(Ind = 1:nrow(axescord), axescord)
		axesvec = as.factor(c("AXA", unlist(lapply(c("XXX", "YYY", 
												"ZZZ"), function(x) {
											rep(c("AXA", x), 2)
										}))))
		formataxes <- function(x) {
			return(sprintf(paste("ATOM  %5.0f  %2s  %3s X %3.0f     %7.3f %7.3f %7.3f  1.00 %5.2f", 
									sep = ""), x[1], c("CA", "X", "Y", "Z")[x[5]], 
							substr(levels(axesvec)[x[5]], 1, 3), x[1], x[2], 
							x[3], x[4], x[4]))
		}
		axescord = cbind(axescord, as.numeric(axesvec))
		axes = apply(axescord, 1, formataxes)
		axes[14] = "TER"
		return(axes)
	}
	formatline <- function(x, classvec, ...) {
		if (!is.null(classvec)) {
			return(sprintf(paste("ATOM  %5.0f  ID  %3s %1s %3.0f     %7.3f %7.3f %7.3f  1.00 %5.2f", 
									sep = ""), x[1], substr(levels(classvec)[x[5]], 
									1, 3), LETTERS[x[5]], x[1], x[2], x[3], x[4], 
							x[4]))
		}
		if (is.null(classvec)) {
			return(sprintf(paste("ATOM  %5.0f  ID  SAM %5.0f     %7.3f %7.3f %7.3f  1.00 %5.2f", 
									sep = ""), x[1], x[1], x[2], x[3], x[4], x[4]))
		}
	}
	bels <- df
	if (scaled) {
		bels <- btt(bels)
	}
	rn <- row.names(bels)
	bels <- cbind(Ind = 1:nrow(bels), bels)
	if (!is.null(classvec)) 
		bels = cbind(bels, vec = as.numeric(classvec))
	pdb <- c(addaxes(), apply(bels, 1, formatline, classvec = classvec))
	if (writehtml) {
		pdbfilename = paste(filenamebase, ".pdb", sep = "")
		write(pdb, pdbfilename)
		htmlfilename = paste(filenamebase, ".html", sep = "")
		my.jmol3D(pdbfilename, classvec = classvec, title = title, 
				filename = htmlfilename, jmolpath=jmolpath, ...)
	}
	if (writepdb) {
		pdbfilename = paste(filenamebase, ".pdb", sep = "")
		write(pdb, pdbfilename)
	}
	return(as.matrix(pdb))
}

#' my.jmol3D - chime3D of made4 using jmol
#' @param ...
#' @returnType ...
#' @return ...
#' @author daniel
#' @export
my.jmol3D <- function (pdbfilename, classvec = NULL, title = NULL, filename = "output.html", 
		point.size = 40, cols = NULL, jmolpath, ...) 
{
	outfile <- file(filename, "w")
	if (!is.null(classvec)) 
		classvec = made4:::checkfac(classvec)
	if (is.null(title)) 
		title = filename
	cat("<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.01 Transitional//EN\"
            \"http://www.w3.org/TR/html4/loose.dtd\">
<html>\n <head><title>", title, "</title>

        <script src=\"",jmolpath,"/Jmol.js\" type=\"text/javascript\"></script>
        <script type=\"text/javascript\">
                jmolInitialize(\"",jmolpath,"\", \"JmolAppletSigned0.jar\");    
                _jmol.menuCssText = \"style='width:180px'\";
        </script>
		</head>\n", 
			"<body bgcolor=\"#FFFFFF\">\n", "<h1 align=\"center\">", 
			title, "</h1>\n", file = outfile, sep = "")
	cat(
			"<br><b>Rotate</b> using your <font color =\"green\">left</font>", 
			"mouse button. <b>Zoom</b> by pressing the <font color=\"green\">shift key </font> while using your <font color='green'>", 
			"left</font> mouse button.\n", file = outfile, sep = "\n")
	
	cat("<table cellpadding=\"5\" cellspacing=\"6\" border=\"0\" width=\"850\" align=\"center\">\n", 
			file = outfile, sep = "\n")
	
	cat("<tr> <td colspan=\"4\"> <hr align=\"center\" width=\"100%\" style=\"background-color:green;\" size=\"5\"></td></tr>\n", 
			file = outfile, sep = "\n")
	cat("<tr> <td width=\"300\" valign=\"top\">\n", file = outfile, 
			sep = "\n")
	cat("\n<p><b>Spin Graph</b></p>\n<ul><li>Start continuous spinning", 
			"<script type=\"text/javascript\">jmolButton(\"spin on\",\"Spin on\",\"\",\"Start continuos spinning\")</script>",
			"</li>\n", 
			"<li> Stop continuous spinning ", 
			"<script type=\"text/javascript\">jmolButton(\"spin off\",\"Spin off\",\"\",\"Stop continuos spinning\")</script>",
			"</li></ul>\n", 
			file = outfile, sep = "\n")
	cat("\n<p><b>Restore</b></p><ul><li>original rotation and zoom", 
			"<script type=\"text/javascript\">jmolButton(\"reset; rotate x 180\",\"Restore\",\"\",\"Restore rotation and zoom\")</script>", 
			"</li></ul>\n", file = outfile, 
			sep = "\n")
	doclasscol <- function(classvec) {
		cat("\n<p><b>Colour classses</b>\n<ul>", file = outfile, 
				sep = "/n")
		nclass <- length(levels(classvec))
		letts <- LETTERS[1:nclass]
		graphcols = cols
		if (is.null(graphcols)) 
			graphcols <- getcol(nc = c(1:nclass), palette = "colours1")
		for (i in c(1:nclass)) {
			cat("\n<li>", levels(classvec)[i], " samples (", 
					graphcols[i], ")\n", "
<script type=\"text/javascript\">jmolButton(\"select *; colour atoms grey; select *",letts[i], "; colour atoms ", graphcols[i],"\",\"Highlight\",\"\",\"Colour\")</script></li>\n", 
					file = outfile, sep = "")
		}
		cat("\n</ul><p><b>Restore</b></p><ul><li>Original Colours\n", "<script type=\"text/javascript\">jmolButton(\"",	
				file = outfile, sep = "")
		for (i in c(1:nclass)) {
			cat(" select *", letts[i], "; colour atoms ", graphcols[i],	";", 
				file = outfile, sep = "")
		}
		cat("\",\"Colour\",\"\",\"Colour\")</script></li></ul>\n", file = outfile, sep = "")

	}
	if (!is.null(classvec)) 
		doclasscol(classvec)
	cat("</td>\n\n<td width=\"5\" bgcolor=\"green\"><br></td> <td width=\"10\"><br></td>\n", 
			file = outfile, sep = "\n")
	cat("<td valign=\"top\" width=\"450\" bgcolor=\"white\">\n", file = outfile, 
			sep = " ")
	cat("
<script type=\"text/javascript\">
   jmolApplet([\"100%\",\"100%\"],\"set echo top left;echo loading...;refresh; load ", pdbfilename, "; select *; connect DELETE; colour atoms grey; spacefill ",point.size,"; set ambient ",point.size,"; select *X; colour red; spacefill off; select none; set axes on; rotate 180;", 
file = outfile,	sep = " ")
	cat("select XXX3; label F1; colour labels red; select YYY7; label F2; colour labels green; select ZZZ11; label F3; colour labels blue; select *B; colour atoms blue; zoom 80; rotate X 20; rotate Y -10;", 
			file = outfile, sep = " ")

	if (!is.null(classvec)) {
		nclass <- length(levels(classvec))
		letts <- LETTERS[1:nclass]
		graphcols = cols
		if (is.null(graphcols)) 
			graphcols <- getcol(nc = c(1:nclass), palette = "colours1")
		for (i in c(1:nclass)) {
			cat("select *", letts[i], "; colour atoms ", graphcols[i], 
					";", file = outfile, sep = " ")
		}		
	}
	if (is.null(classvec)) 
		cat("select SAM; colour atoms red'", file = outfile, 
				sep = " ")
	cat("echo\")</script>\n</td>\n", 
			file = outfile, sep = "\n")
	cat("</tr>\n</table>\n</body></html>", file = outfile, sep = "\n")
	close(outfile)
}

#' hclust.ave - hclust(d,"average",...)
#' @param hclust params except method="average"
#' @returnType hclust
#' @return hclust
#' @author daniel
#' @export
hclust.ave <- function (d,...) {
	return(hclust(d=d,method="average",...))
}


#' hclust.diana - as.hclust(diana(x,diss=TRUE,...))
#' @param hclust params except diss=TRUE
#' @returnType hclust
#' @return hclust
#' @author daniel
#' @export
hclust.diana <- function (x,...) {
	return(as.hclust(diana(x=x,diss=TRUE,...)))
}

#' dianaHook - diana Hook
#' @param data matrix or data frame, k value
#' @returnType cutree returns
#' @return cutree returns
#' @author daniel
#' @export
dianaHook = function(this_dist,k){
	tmp = diana(this_dist,diss=TRUE)
	assignment = cutree(tmp,k)
	return(assignment)
}

#' my.central.coord - Central coordinate for a class after ordination
#' @param eigenarrays, array classification
#' @returnType data.frame
#' @return central point coordinates
#' @author daniel
#' @export
my.central.coord <- function(dfxy, classvec) {
	
	f1 <- function(cl) {
		n <- length(cl)
		cl <- as.factor(cl)
		x <- matrix(0, n, length(levels(cl)))
		x[(1:n) + n * (unclass(cl) - 1)] <- 1
		dimnames(x) <- list(names(cl), levels(cl))
		data.frame(x)
	}
	
	dfxy <- data.frame(dfxy)
	if (!is.data.frame(dfxy)) 
		stop("Non convenient selection for dfxy")
	if (any(is.na(dfxy))) 
		stop("NA non implemented")
	if (!is.factor(classvec)) 
		stop("factor expected for classvec")
	
	wt = rep(1, length(classvec))
	#xax = 1 
	#yax = 2
	
	dfdistri <- f1(classvec) * wt
	
	w1 <- unlist(lapply(dfdistri, sum))
	
	dfdistri <- t(t(dfdistri)/w1)
	
	df <- data.frame()
	
	for (col in colnames(dfxy)) {
		#print(col)
		mat <- as.matrix(t(dfdistri)) %*% dfxy[,col]
		#print(mat)
		for (r in  rownames(mat)) {
			df[r,col] <- mat[r,1]
		}		
	}
	
	#print(paste(df[,1],collapse=","))
	#print(paste(df[,2],collapse=","))
	
	#coox <- as.matrix(t(dfdistri)) %*% dfxy[, xax]
	#cooy <- as.matrix(t(dfdistri)) %*% dfxy[, yax]
	
	#print(paste(coox,collapse=";"))
	#print(paste(cooy,collapse=";"))
	return(df)
}
	


# ----- Define a function for plotting a matrix ----- #
my.CorrPlot <- function(x, ...){
	min <- min(x)
	max <- max(x)
	yLabels <- rownames(x)
	xLabels <- colnames(x)
	title <-c("")
	# check for additional function arguments
	if( length(list(...)) ){
		Lst <- list(...)
		if( !is.null(Lst$zlim) ){
			min <- Lst$zlim[1]
			max <- Lst$zlim[2]
		}
		if( !is.null(Lst$yLabels) ){
			yLabels <- c(Lst$yLabels)
		}
		if( !is.null(Lst$xLabels) ){
			xLabels <- c(Lst$xLabels)
		}
		if( !is.null(Lst$title) ){
			title <- Lst$title
		}
	}
	# check for null values
	if( is.null(xLabels) ){
		xLabels <- c(1:ncol(x))
	}
	if( is.null(yLabels) ){
		yLabels <- c(1:nrow(x))
	}
	
	layout(matrix(data=c(1,2), nrow=1, ncol=2), widths=c(4,1), heights=c(1,1))
	
	# Red and green range from 0 to 1 while Blue ranges from 1 to 0
	ColorRamp <- rgb( 	seq(0,1,length=256),  # Red
			seq(0,1,length=256),  # Green
			seq(1,0,length=256))  # Blue
	ColorLevels <- seq(min, max, length=length(ColorRamp))
	
	# Reverse Y axis
	reverse <- nrow(x) : 1
	yLabels <- yLabels[reverse]
	x <- x[reverse,]
	
	# Data Map
	par(mar = c(3,5,2.5,2))
	image(1:length(xLabels), 1:length(yLabels), t(x), col=ColorRamp, xlab="",
			ylab="", axes=FALSE, zlim=c(min,max))
	if( !is.null(title) ){
		title(main=title)
	}
	#axis(BELOW<-1, at=1:length(xLabels), labels=xLabels, cex.axis=0.7)
	
	#axis(BELOW<-1, at=1:length(xLabels), lab=FALSE, cex.axis=0.7)
	axis(LEFT <-2, at=1:length(yLabels), labels=yLabels, las= HORIZONTAL<-1, cex.axis=0.7)
	
	text(1:length(xLabels), par("usr")[3]-1.5, srt=45, adj=1,
			labels=xLabels,		
			xpd=T, cex=0.7)
	
	
	# Color Scale
	par(mar = c(3,2.5,2.5,2))
	image(1, ColorLevels,
			matrix(data=ColorLevels, ncol=length(ColorLevels),nrow=1),
			col=ColorRamp,
			xlab="",ylab="",
			xaxt="n")
	
	layout(1)
}
# ----- END plot function ----- #


## Houtan functions

func.list<-list()

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### make volcano plot plus
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
func.list$volcano.plus<-function(x,fold.change.col,pv.col, title.plot, cut.line, fold.cut1, fold.cut2, pv.adj.col,ncolors=1, text=NA, angle=-45) {
	
	op <- par(no.readonly = TRUE)
	p<-x[,pv.col]

	if (all(is.na(p))==FALSE) {
		
		M<-x[,fold.change.col]

		upr = p<=cut.line & M >=fold.cut1
		dwr = p<=cut.line & M<=fold.cut2
		
		q <- x[,pv.adj.col]
		q.f <- x[(upr|dwr),pv.adj.col]
		minq.f <- 0
		maxq.f <- 0
	
		if (length(q.f) > 0) {
			minq.f <- min(q.f)
			maxq.f <- max(q.f)
			if (is.na(minq.f)) {
				print(q.f)
			}
		}
	
		library('marray')
		
		Gcol <- maPalette(low = "#C8FFC8", high = "#006400", k = ncolors)
		Rcol <- maPalette(low = "#FFC8C8", high = "#640000", k = ncolors)

		xlim.r<-(range(x[,fold.change.col]))
		
		if ((length(q.f) > 0) && ( (maxq.f-minq.f) > 0 )) {
			par(fig=c(0, 0.8, 0, 1), mar=c(4, 4, 4, 1))
		}

		if(abs(xlim.r[1])>abs(xlim.r[2])){
			plot(
					x[,fold.change.col], #x-axis
					-1*log10(x[,pv.col]), #y-axis
					xlim=c(xlim.r[1], abs(xlim.r[1])), #x-axis limits
					main=title.plot,
					xlab="Gene Expression\nlog2(fold change)",
					ylab="-1 * log10 of the Significance",
					cex=.5,
					pch=20
			)
				
			if ((length(q.f) > 0) && ( (maxq.f-minq.f) > 0 )) {
				points(M[upr], -log10(p[upr]), col=Rcol[(as.integer(((q[upr]-minq.f)/((maxq.f-minq.f)/(ncolors-1))))+1)],cex=.5,pch=20)
				points(M[dwr], -log10(p[dwr]), col=Gcol[(as.integer(((q[dwr]-minq.f)/((maxq.f-minq.f)/(ncolors-1))))+1)],cex=.5,pch=20)
			}	
			else {
				points(M[upr], -log10(p[upr]), col=Rcol[ncolors],cex=.5,pch=20)
				points(M[dwr], -log10(p[dwr]), col=Gcol[ncolors],cex=.5,pch=20)			
			}
			
			#points(M[upr], -log10(p[upr]), col="red")
			#points(M[dwr], -log10(p[dwr]), col="green")
			
			abline(h= -log10(cut.line), lty=3, lwd=1)
			abline(v= fold.cut1, lty=3, lwd=1, col="black")
			abline(v= fold.cut2, lty=3, lwd=1, col="black")
		}else{
			plot(
					x[,fold.change.col], #x-axis
					-1*log10(x[,pv.col]), #y-axis
					xlim=c(-1*xlim.r[2], xlim.r[2]), #x-axis limits
					main=title.plot,
					xlab="Gene Expression\nlog2(fold change)",
					ylab="-1 * log10 of the Significance",
					cex=.5,
					pch=20
			)
		
			if ((length(q.f) > 0) && ( (maxq.f-minq.f) > 0 )) {
				points(M[upr], -log10(p[upr]), col=Rcol[(as.integer(((q[upr]-minq.f)/((maxq.f-minq.f)/(ncolors-1))))+1)],cex=.5,pch=20)
				points(M[dwr], -log10(p[dwr]), col=Gcol[(as.integer(((q[dwr]-minq.f)/((maxq.f-minq.f)/(ncolors-1))))+1)],cex=.5,pch=20)
			}
			else {
				points(M[upr], -log10(p[upr]), col=Rcol[ncolors],cex=.5,pch=20)
				points(M[dwr], -log10(p[dwr]), col=Gcol[ncolors],cex=.5,pch=20)			
			}
			
			#points(M[upr], -log10(p[upr]), col="red")
			#points(M[dwr], -log10(p[dwr]), col="green")
		
			abline(h= -log10(cut.line), lty=3, lwd=1)
			abline(v= fold.cut1, lty=3, lwd=1, col="black")
			abline(v= fold.cut2, lty=3, lwd=1, col="black")
		}
		if(!is.na(text[1])){
				text(M[text], -log10(p[text]), labels=rownames(x)[text],cex=0.5,srt=angle,offset=0.25,pos=4)		
		}
		
		if ((length(q.f) > 0) && ( (maxq.f-minq.f) > 0 )) {
			ColorLevels <- seq( minq.f, maxq.f, ((maxq.f-minq.f)/(ncolors-1)))
		
			#print(paste(minq.f,maxq.f))
			#print(ncolors)
			#print(ColorLevels)
			
			#Stay on same page and set up region and coordinates for legend.
			par(new=TRUE)
			par(fig=c(0.8, 0.9, 0.2, 0.8), mar=c(4, 0, 4, 3))
		
			#maColorBar(ColorLevels, col=Gcol, horizontal=FALSE, k=0, cex.axis=.8)
		
			image(1, seq(1,ncolors,1),
				matrix(data=seq(1,ncolors,1), ncol=ncolors, nrow=1),
				col=Gcol,
				xlab="",ylab="",
				axes=FALSE
			)
			
	
			
			axis(4, at = seq(0.5, (ncolors-0.5), 1), 
					labels = rep('',ncolors),cex.axis=.5,mgp=c(0, .2, 0)
			)
		
			par(new=TRUE)
			par(fig=c(0.9, 1, 0.2, 0.8), mar=c(4, 0, 4, 3))
			
			#maColorBar(ColorLevels, col=Rcol, horizontal=FALSE, k=ncolors, cex.axis=.8)
		
			image(1, seq(1,ncolors,1),
					matrix(data=seq(1,ncolors,1), ncol=ncolors, nrow=1),
					col=Rcol,
					xlab="",ylab="",
					axes=FALSE
			)
			
			axis(4, at = seq(0.5, (ncolors-0.5), 1), 
					labels = signif(ColorLevels,2),cex.axis=.5,mgp=c(0, .2, 0),
			)
		}
	}	
	par(op)
}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### make volcano plot
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
func.list$volcano<-function(x,fold.change.col,pv.col, title.plot, cut.line, fold.cut1, fold.cut2) {
	p<-x[,pv.col]
	M<-x[,fold.change.col]
	upr = p<=cut.line & M >=fold.cut1
	dwr = p<=cut.line & M<=fold.cut2
	
	xlim.r<-(range(x[,fold.change.col]))
	
	if(abs(xlim.r[1])>abs(xlim.r[2])){
		plot(
				x[,fold.change.col], #x-axis
				-1*log10(x[,pv.col]), #y-axis
				xlim=c(xlim.r[1], abs(xlim.r[1])), #x-axis limits
				main=title.plot,
				xlab="Gene Expression\nlog2(fold change)",
				ylab="-1 * log10 of the Significance",
				cex=.5,
				pch=20
		)
		
		points(M[upr], -log10(p[upr]), col="red",cex=.5,pch=20)
		points(M[dwr], -log10(p[dwr]), col="green",cex=.5,pch=20)			
				
		abline(h= -log10(cut.line), lty=3, lwd=1)
		abline(v= fold.cut1, lty=3, lwd=1, col="black")
		abline(v= fold.cut2, lty=3, lwd=1, col="black")
	}else{
		plot(
				x[,fold.change.col], #x-axis
				-1*log10(x[,pv.col]), #y-axis
				xlim=c(-1*xlim.r[2], xlim.r[2]), #x-axis limits
				main=title.plot,
				xlab="Gene Expression\nlog2(fold change)",
				ylab="-1 * log10 of the Significance",
				cex=.5,
				pch=20
		)
		
		points(M[upr], -log10(p[upr]), col="red",cex=.5,pch=20)
		points(M[dwr], -log10(p[dwr]), col="green",cex=.5,pch=20)			
				
		abline(h= -log10(cut.line), lty=3, lwd=1)
		abline(v= fold.cut1, lty=3, lwd=1, col="black")
		abline(v= fold.cut2, lty=3, lwd=1, col="black")
	}
	
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### student t.test
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
func.list$studentT<-function(x,s1,s2,var.pvalue=0.01,...) {
	
	x1 <- x[s1]
	x2 <- x[s2]
	x1 <- as.numeric(x1)
	x2 <- as.numeric(x2)
	
	my.t.test.p.value <- function(...) {
	    obj<-try(t.test(...), silent=TRUE)
	    if (is(obj, "try-error")) return(NA) else return(obj$p.value)
	}

	var.equal.bool = FALSE

	vt <- var.test(x1,x2)$p.value
	if (is.nan(vt)){
		x1[1] <- x1[1]+0.000000000000001
		x2[1] <- x2[1]+0.000000000000001
		vt <- var.test(x1,x2)$p.value
	}
	
	if (vt <= var.pvalue) {
		var.equal.bool = FALSE
	}
	else {
		var.equal.bool = TRUE
	}

	
	t.out.pvalue <- my.t.test.p.value(x1,x2, alternative="two.sided", na.rm=TRUE, var.equal=var.equal.bool, ...)
	
	return(as.numeric(t.out.pvalue))
}
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### vlookup method
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
func.list$vlookup<-function(val, df, col){
	
	df[df[1] == val, col][1]
}
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### wilcox exact test (deals with ties)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
func.list$wilcox.exact<-function(x,s1,s2) {
	x1 <- x[s1]
	x2 <- x[s2]
	x1 <- as.numeric(x1)
	x2 <- as.numeric(x2)
	t.out <- wilcox.exact(x1,x2, alternative="two.sided", exact=TRUE, na.rm=TRUE)
	out <- as.numeric(t.out$p.value)
	return(out)
}
#to sort data.frame matrix.  Use this function.  then use sort.data.frame(dat, key = "LOC")
sort.data.frame <- function(x, key, ...) {
	if (missing(key)) {
		rn <- rownames(x)
		if (all(rn %in% 1:nrow(x))) rn <- as.numeric(rn)
		x[order(rn, ...), , drop=FALSE]
	} else {
		x[do.call("order", c(x[key], ...)), , drop=FALSE]
	}
}
gg.volcano<-function(){
	#plot volcano
	red.axis<-c(1.3)
	func.list$volcano(data.1,fold.change="fold.change",pv.col="pvalues", title.plot="Volcano Plot\nNC vs. NC AZA",cut.line=red.axis)
	
	xlim.r<-(range(data.1[,"fold.change"]))
	p.1<-ggplot(zz, na.rm=T, aes(y=-1*log10(rawp.t), x=fold.change, size=(rawp.t<10^(-1*red.axis) & fold.change>2), colour=(rawp.t<10^(-1*red.axis) & fold.change>2)))
	p.1<-ggplot(zz, na.rm=T, aes(y=-1*log10(rawp.t), x=fold.change))
	p.1<-p.1 + 
			opts(title = "NC (Control) vs NC (AZA)") +
			geom_point(na.rm=T) + 
			#coord_equal() +
			#geom_text(size=2,hjust=-0.25, vjust=-0.5, colour="red") +
			scale_y_continuous("-1*Log10(q-value)\nSignificance") +
			#geom_point(x>=2, colour="red") +
			#facet_grid(. ~ CPG_ISLAND)	+
			#stat_smooth(method=lm, se=TRUE, fullrange=FALSE, colour=(values=c("NO"="BLUE")), linetype=2, size=0.25) +
			geom_hline(yintercept = c(red.axis), colour="black", size=0.5, linetype=2)+
			geom_vline(xintercept = c(-2), colour="black", size=0.5, linetype=2)+
			geom_vline(xintercept = c(2), colour="black", size=0.5, linetype=2)+
			geom_hline(yintercept = 0, colour="black", size=0.5, linetype=1)+
			geom_vline(xintercept = 0, colour="black", size=0.5, linetype=1)+
			#coord_equal() +
			scale_colour_gradientn(colour = rev(jet.colors(200)),  "-log10(FDR Methylation)")+
			opts(legend.position = "none")+
			scale_colour_manual(values = c("Black","RED"), "Significant?")   + 
			scale_x_continuous("log2(Gene Expression Fold Change)") 
	x11();p.1
}

color.map <- function(mol.biol) {
	if (mol.biol==1) "#24ac68" #black
	else if (mol.biol==2) "#88a8f9" #blue
	else if (mol.biol==3) "#e04cf0" #light blue
	else if (mol.biol==4) "#0000FF" #BLUE
	else if (mol.biol==5) "turquoise" #red
	else if (mol.biol==6) "purple" #red
	else if (mol.biol==7) "yellow" #red
	else if (mol.biol==8) "grey" #red
	else "#FFFFFF" #white
} 
color.map.labels <- function(mol.biol) {
	if (mol.biol=="TC2") "#24ac68" #black
	else if (mol.biol=="OI1") "#88a8f9" #blue
	else if (mol.biol=="OI3") "#e04cf0" #light blue
	else if (mol.biol==4) "#0000FF" #BLUE
	else if (mol.biol==5) "turquoise" #red
	else if (mol.biol==6) "purple" #red
	else if (mol.biol==7) "yellow" #red
	else if (mol.biol==8) "grey" #red
	else "#FFFFFF" #white
}

HeatMap.2<-function (object, size, maintitle, hc) 
{
	require(marray)
	require(gplots)
	require(gtools)
	require(gdata)
	samples = colnames(object)
	names = rownames(object)
	#variance using size
	genes.var = apply(object, 1, sd, na.rm = TRUE)
	genes.var.select = order(genes.var, decreasing = T)[1:size]
	DD.s = object[genes.var.select, ]
	samples = colnames(DD.s)
	names = rownames(DD.s)
	c <- rainbow(ncol(DD.s), start = 0, end = 0.3)
	rc <- rainbow(nrow(DD.s), start = 0, end = 0.3)
	rbg = maPalette(low = "green", high = "red", mid = "black", 
			k = 50)
	heatmap.2(
			DD.s, 
			labCol = samples, 
			labRow = names, 
			hclustfun = hc,
			scale = "none", 
			#RowSideColor=probe.cc,
			#ColSideColors=cc.col,
			col = rbg, 
			margin = c(12, 12), 
			key=T,
			#symkey=FALSE,
			#density.info="none",
			Rowv=T, # cluster by rows
			Colv=T, #cluster by columns
			cexRow=0.45,
			cexCol=1,
			sepwidth=c(0.02,2),
			dendrogram=c("both"),
			keysize=1,
			trace = "none"
	#tracecol = "cyan"
	)
	subtitle = paste(as.character(size), " high variance genes by SD")
	title(main = maintitle, sub = subtitle)
}

