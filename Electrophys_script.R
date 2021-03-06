###########################################################
#####
##### R script to analyze electrophys experiment data #####
#####
##### author: Claire Churchhouse
##### contact: cchurch@broadinstitute.org
#####
###########################################################


##### define user parameters
wkd <- "/Users/cchurch/Documents/AndrewAllen" # set the working directory
filename <- "16-4-4-de\ novo\ variants.csv" # set the file name containing the data
dir.out <- "/Users/cchurch/Documents/AndrewAllen" # set the output directory
file.volt <- "voltages.txt" # set the file containing voltages (path relative to working directory)
current.orient <- c("neg", "pos")
Loess <- FALSE
Rugc <-  8.314472 # R universal gas constant
Farad <- 9.64853399*10^4 # Faraday constant
Temp <- 293 # temperature in Kelvins to infer Boltzmann parameters


#### Set parameter defaults that are used if not defined by user
if( !exists("file.volt") ){ file.volt <- "voltages.txt" } 
if( !exists("current.orient") ){ current.orient <- rep("pos", length(PARAMS)) } 
PARAMS <- c("c1","c2")
LABELS <- c("C1", "C2")

#### Install required packages
#install.packages("ggplot2")
#library(ggplot2)
#install.packages("chemCal")
#library(chemCal)
#install.packages("MASS")
#library(MASS)
install.packages("fitdistrplus")
library(fitdistrplus)
install.packages("drc") # for Boltzman
library(drc)

#### Create a function to write comments to a logfile "output/logfile.txt"
output.log <- function( somestring ){
	sink( paste( dir.out, "/logfile.txt", sep=""), append=TRUE)
	cat( somestring )
	cat("\n")
	sink()
}


#### Create the output subdirectories
dir.create(dir.out, showWarnings = TRUE, recursive = TRUE, mode = "0777")
dir.tab.out <- paste(dir.out,"/Tables",sep="")
dir.create(dir.tab.out, showWarnings = TRUE, recursive = FALSE, mode = "0777")
dir.plot.out <- paste(dir.out,"/Plots",sep="")
dir.create(dir.plot.out, showWarnings = TRUE, recursive = FALSE, mode = "0777")


##### report user parameters to logfile
output.log("###########################################################\n")
output.log( format(Sys.time(), "%a %b %d %H:%M:%S %Y"))
output.log( "\n###### Reporting parameters: ######\n" )
output.log( paste("Working directory:", wkd) )
setwd(wkd)
output.log( paste("Output directory:", dir.out) )
output.log( paste("Voltages file:", file.volt) )
output.log( paste("Will analyze Function.C1 and Function.C2 data") )
output.log( paste(c("Respective current orientation is:", current.orient)) )
output.log( paste("Reading data from:", filename) )
##### read in data
data <- read.csv( filename, header=TRUE, as.is=TRUE, skip=1)
data.swp <- read.csv( filename, header=TRUE, as.is=TRUE, nrows=1)
swp.cols <- colnames(data.swp[ , grep("Sweep", colnames(data.swp))])

#### Read in the voltages from voltages.txt
volt <- read.table( file.volt, header=F)
volt <- volt[,1]
output.log( paste("Voltages are:") )
output.log( volt )



##### Format the data
output.log( "\n\n###### Formatting the data: ######\n" )
# write a function to format the Sweep variable name
trim.sweep <- function(Strng){
	a <- unlist(gregexpr(pattern="\\.", Strng ))
	if( length(a)>= 2 ){ return(substr(Strng, 1, a[2]-1)) }else{ return(Strng) }
	}
# use it for formatting
swp.edit <- rep(0, length(swp.cols))
for( i in 1:length(swp.cols)){
	swp.edit[i] <- trim.sweep(swp.cols[i])
}
swp.fac <- as.factor(swp.edit) # now swp.fac is the sweep variable label

# Identify the Capacitance columns:
cap.col <- grep( "Capacitance", data.swp[1,])
output.log( paste("Number of capacitance columns:", length(cap.col)) )
# Identify the C1 activation columns:
c1.col <- (1:ncol(data))[grepl( "Function.C1", data.swp[1,]) & !grepl( "C2", data.swp[1,])] 
output.log( paste("Number of C1 columns:", length(c1.col)) )
# Identify the C2 inactivation columns:
c2.col <- (1:ncol(data))[grepl( "Function.C2", data.swp[1,]) & !grepl( "C1", data.swp[1,])]
output.log( paste("Number of C2 columns:", length(c2.col)) )
# Identify the C1/C2 columns:
rat.col <- grep( "C1 / C2",  data.swp[1,])
output.log( paste("Number of C1/C2 columns:", length(rat.col)) )
# Identify the cell type column:
type.col <- which( colnames(data) == "Cell.Type")
output.log( paste("Cell type column is number:", type.col) )
# If cell type column is not found report error:
if( length(type.col)==0 ){
	output.log(paste("ERROR: Cell.Type column not found!") )
}
data[,type.col] <- as.factor(data[,type.col])

# Check that there are equal number of capacitance, C1, C2 and ratio columns identified:
if( length(cap.col) != length(c1.col) ){
	output.log(paste("ERROR: Unequal number of capacitance values and C1 values identified!") )
}
if( length(cap.col) != length(c2.col) ){
	output.log(paste("ERROR: Unequal number of capacitance values and C2 values identified!") )
}
if( length(cap.col) != length(rat.col) ){
	output.log(paste("ERROR: Unequal number of capacitance values and C1/C2 values identified!") )
}


# Remove cells in which there are any NaN values as entries:
nan.count <- apply( apply( data[,sort(c(c1.col,c2.col))], 1, is.nan), 1, sum)
bad.cell <- which( nan.count > 0 )
if( length(bad.cell) > 0 ){data <- data[ -bad.cell,]}
output.log( paste("Removed",  length(bad.cell), "cells for having at least one NaN entry.") )

# Count the number of Well IDs and unique ones:
no.IDs <- length(data$Well.ID)
no.IDs.uniq <- length( unique( data$Well.ID))
output.log( paste("Number of well IDs:", no.IDs) )
output.log( paste("Number of unique well IDs:", no.IDs.uniq) )
if(no.IDs.uniq != no.IDs){
	output.log( paste("ERROR: Not all the Well IDs are unique!") )
}
# Count the number of wells of each cell type:
output.log("Number of cells per well type:")
output.log( print(table(data$Cell.Type))[0] ) 
max.type <- max(table(data$Cell.Type))





################# Analysis ############
output.log( "\n###### Analysis: ######\n" )




##### Analysis of the capacitance data
output.log( "### Capacitance: ###" )
# Assign capacitance value of each well to be the first measured value, call it cap1:
cap1 <- data[,cap.col[1]]
names(cap1) <- data[,1]
# Write cap1 of each well to file:
filename <- paste(dir.tab.out, "/cap1_well.txt",sep="")
output.log( paste("Capacitance of each well i.e. cap1: ",filename,sep=""))
colnames(cap1) <- "cap1"
write.table( cap1, file=filename, row=TRUE, quote=F, sep="\t", col.names=NA)

# Calculate mean and SEM of cap1 per cell type:
OUT.ac <- NULL
OUT.sem <- NULL
for( i in unique(data$Cell.Type)){
	these <- which( data$Cell.Type == i)
	OUT.ac <- c( OUT.ac, mean(cap1[these], na.rm=TRUE))
	OUT.sem <- c( OUT.sem, sd(cap1[these])/sqrt(sum(!is.na(cap1[these]))) )
}
OUT <- cbind(OUT.ac, OUT.sem)
rownames(OUT) <- unique(data$Cell.Type)
colnames(OUT) <- c("mean_cap1", "SEM")
# Write mean cap1 and SEM per cell type to file:
filename <- paste(dir.tab.out, "/mean_cap1_celltype.txt",sep="")
output.log( paste("Mean cap1 and SEM of each cell type: ",filename,sep=""))
write.table( OUT, file=filename, row=TRUE, col.names=NA, quote=FALSE, sep="\t")

n <- ceiling(sqrt(length(unique(data$Cell.Type))))

# Plot distribution of cap1 per cell type in a single boxplot:
filename <- paste(dir.plot.out, "/cap1_celltype_boxplots.jpeg",sep="")
jpeg( file=filename, height=4200, width=4200, units="px", quality=75, res=600)
OUTcap <- NULL
for( i in unique(data$Cell.Type)){
	these <- which( data$Cell.Type == i)
	OUTcap <- cbind( OUTcap, c(cap1[these], rep(NA, max.type-length(these))))
}

ylm <- 1.2*max(apply(OUTcap,2,max, na.rm=T))
ymx <- apply(OUTcap,2,max, na.rm=T)
boxplot( OUTcap , names=rownames(OUT), ylim=c(0,ylm) , main="Cap1 for all cell types", ylab="Cap1 (pF)", xlab="Cell type" )
for(i in 1:length(unique(data$Cell.Type))){
	ctype <- unique(data$Cell.Type)[i]
	cnt <- length(which(data$Cell.Type == ctype))
	text( i, 1.1*ymx[i], paste(cnt))
}

dev.off()
output.log( paste("Boxplots of cap1 per cell type:", filename,sep="") )


# Barplot of mean and SEM of cap1 per cell type in a single barplot:
filename <- paste(dir.plot.out, "/mean_cap1_celltype_barplot.jpeg",sep="")
jpeg( file=filename, height=4200, width=4200, units="px", quality=75, res=600)
b1 <- OUT.ac + OUT.sem
b2 <- OUT.ac - OUT.sem
pos <- max(abs(b1)) > max(abs(b2))
if(pos){ blim <- max(b1)}else{ blim <- min(b2)}
Bars <- barplot ( OUT.ac , names=rownames(OUT), main="Mean cap1 for each cell type", ylab="mean cap1 (pF)", xlab="Cell type" , ylim=c(0, 1.2*blim) ) 
segments(Bars, b2, Bars, b1, lwd = 1.5)
arrows(Bars, b2, Bars, b1, lwd = 1.5, angle = 90, code = 3, length = 0.05)
for(i in 1:length(unique(data$Cell.Type))){
	ctype <- unique(data$Cell.Type)[i]
	cnt <- length(which(data$Cell.Type == ctype))
	text(Bars[i], min(c(b1[i],b2[i]))-3 , paste(cnt))
}
dev.off()
output.log( paste("Barplot of mean cap1 per cell type:", filename,sep="") )





##### Analysis of both the activation (C1) and inactivation (C2) currents
output.log( "\n### Current: ###\n" )

BparamsOUT <- NULL
for(p in 1:length(PARAMS)){
	
	param.c <- PARAMS[p]
	label.c <- LABELS[p]
	output.log( paste("Analyzing current", label.c,"\n") )
	c.col <- get(paste( param.c, ".col", sep=""))
	M <- data[,c.col]
	
	# set the current normalization factor
	if( param.c == "c1" ){ 
		cf <- apply( M, 1, min, na.rm=TRUE)
		Ih <- -0.5
		output.log( "Current is normalized to minimum value in each well." )
	}
	if( param.c == "c2" ){
		cf <- M[,1]
		Ih <- 0.5
		output.log( "Current is normalized to first recorded value in each well." )	
	}
	N <- matrix( cf, nrow(M), ncol(M))
	cnorm <- M/N # normalized the current
	cnorm <- as.matrix(cnorm)
	if(current.orient[p]=="neg"){ cnorm <- -1*cnorm} # make the normalized currents negative if user specified so.
	#write the normalized current to file
	filename <- paste(dir.tab.out, "/cnorm_well.txt",sep="")
	output.log( paste("Normalized current in each well i.e. cnorm: ",filename,sep=""))
	write.table( cbind(data$Well.ID, data$Cell.Type, signif(cnorm,digits=3)), file=filename, row=FALSE, col=c(colnames(data)[1:2], paste("cnorm_", volt,"V", sep="")), quote=F, sep="\t")
	
	# Plot normalized current against voltage
	for(fix.yaxis in c(FALSE)){
	#for(fix.yaxis in c(TRUE, FALSE)){	
		if( fix.yaxis ){strng <- "_fixed"}else{strng <- ""}
		filename <- paste(dir.plot.out, "/",label.c,"norm_vs_volt_celltype",strng, ".jpeg",sep="")
		output.log( paste("Plots of normalized current",label.c,"against voltage, one plot per cell type:", filename) )
		jpeg( file=filename, height=4200, width=4200, units="px", quality=75, res=600)
		par( mfrow =c(n,n) )
		for( i in unique(data$Cell.Type)){
			these <- which( data$Cell.Type == i)
			lt <- length(these)
			if( fix.yaxis ){
				Xvals <- matrix(volt, byrow=TRUE, nrow=nrow(cnorm), ncol=length(volt))
				Yvals <- cnorm
				}else{
					Xvals <- matrix(volt, byrow=TRUE, nrow=length(these), ncol=length(volt))
					Yvals <- cnorm[these,]
			 }
			plot( Xvals, Yvals, type="n", main=i, xlab="Voltage", ylab=paste(label.c, " normalized",sep=""))
			for(j in 1:length(these) ){
				points( volt, cnorm[these[j],], col=rainbow(lt)[j])
				if(Loess){
					lo <- loess(cnorm[these[j],]~volt)
					xl <- seq(min(volt), max(volt), (max(volt)-min(volt))/1000)
					lines(xl, predict(lo,xl), col=rainbow(lt)[j])		
				}else{
					lines( volt, cnorm[these[j],], col=rainbow(lt)[j])
				}

			}
		}
		dev.off()
	}
	
	# Plot the mean normalized current (with SEM) for each cell type, against voltage
	filename <- paste(dir.plot.out, "/mean_",label.c,"norm_vs_volt_celltype.jpeg",sep="")
	output.log( paste("Plot of mean of normalized current",label.c,"against voltage, all cell types in one plot:", filename) )
	jpeg( file=filename, height=4200, width=4200, units="px", quality=75, res=600)
	lt <- length(unique(data$Cell.Type))
	CM <- NULL; SEM <- NULL
	for( i in unique(data$Cell.Type)){
		these <- which( data$Cell.Type == i )
		cm <- apply( cnorm[these,], 2, mean, na.rm=TRUE)
		CM <- rbind(CM, cm)
		sem <- apply( cnorm[these,], 2, sd, na.rm=TRUE)/sqrt( apply(!is.na(cnorm[these,]),2,sum) )
		SEM <- rbind(SEM, sem)
	}
	plot( matrix(volt, byrow=TRUE, nrow=nrow(CM), ncol=ncol(CM)), CM, col=rainbow(lt), main=paste("Mean of normalized", label.c, "current"), ylab=paste("mean ", label.c, " normalized",sep=""), xlab="Voltage")

	for( i in 1:lt ){
		x <- volt
		y <- CM[i,]
		sd <- SEM[i,]

		if(Loess){
			lo <- loess(y~x)
			xl <- seq(min(x), max(x), (max(x)-min(x))/1000)
			lines(xl, predict(lo,xl), col=rainbow(lt)[i])
		}else{
			lines(matrix(volt, byrow=TRUE, nrow=nrow(CM), ncol=ncol(CM))[i,], CM[i,], col=rainbow(lt)[i])
		}
		segments(x, y-sd,x, y+sd, col=rainbow(lt)[i])
		epsilon = 1
		segments(x-epsilon,y-sd,x+epsilon,y-sd, col=rainbow(lt)[i])
		segments(x-epsilon,y+sd,x+epsilon,y+sd, col=rainbow(lt)[i])
	}
	legend("bottomleft", legend=unique(data$Cell.Type), lty=1, pch=21, col=rainbow(lt))
	dev.off()
		
	# write the mean normalized current and SEM to file
	filename <- paste(dir.tab.out, "/mean_",label.c,"norm_vs_volt_celltype.txt",sep="")
	output.log( paste("Mean of normalized", label.c, "current for each cell type with SEM:", filename) )
	OUT <- NULL; nms <- NULL
	for(i in 1:ncol(CM)){
		OUT <- cbind( OUT, CM[,i], SEM[,i] )
		nms <- c( nms, paste("mean_", volt[i],"V",sep=""), paste("SEM_", volt[i], "V", sep="") )
	}
	rownames(OUT) <- unique(data$Cell.Type)
	colnames(OUT) <- nms
	write.table( OUT, file=filename, row=T, col.names=NA, quote=F, sep="\t")

	##### Fitting the partial Boltzman function
	output.log(paste("\n### Fitting the partial Boltzman function for",param.c,"###"))
	X_norm <- NULL
	for(i in 1:nrow(CM)){
		output.log(paste("\nFor",unique(data$Cell.Type)[i]))
		x <- CM[i,]
		#calculate slope between consecutive points
		m <- abs(x[2:length(x)]-x[1:(length(x)-1)])/(volt[2:length(volt)]-volt[1:(length(volt)-1)])
		# find range to normalize to zero
		# For C1, start for the left
		if( param.c == "c1"){
			k=1
			while( (m[k] < 0.001) == TRUE & (k < length(m)) ){k <- k + 1}
			# normalize the first k points to zero
			output.log(paste("Normalizing the first",k,"points to zero."))
			Mod1 <- lm( x[1:k] ~ volt[1:k] )
		}

		# For C2, start from the right
		if( param.c == "c2"){
			k=length(volt)
			while( (m[k-1] < 0.001) == TRUE & (k > 1) ){k <- k - 1}
			# normalize the last k points to zero
			output.log(paste("Normalizing the last",k,"points to zero."))
			Mod1 <- lm( x[k:(1+length(m))] ~ volt[k:(1+length(m))] )
		}		
		
		output.log(paste("Normalization for slope", signif(Mod1$coefficients[2],3), "and intercept", signif(Mod1$coefficients[1],3)))
		# Normalize all data points to this fitted line
		x_norm <- x - Mod1$coefficients[1] - volt*Mod1$coefficients[2]
		X_norm <- rbind(X_norm, x_norm)		
	}
	rownames(X_norm) <- unique(data$Cell.Type)
	
	#identify which cell types to plot
	if( param.c == "c1"){
		good.plots <- which(abs(CM[,1]) < 0.2)
	}
	if( param.c == "c2"){
		good.plots <- which(abs(CM[,ncol(CM)]) < 0.2)
	}
	#Plot the normalized current versus voltage
	X_norm_plot <- X_norm[ good.plots,]
	SEM_plot <- SEM[ good.plots,]
	CLR <- rainbow(length(unique(data$Cell.Type)))[good.plots]
	filename <- paste(dir.plot.out, "/mean_",label.c,"norm_vs_volt_celltype_Transformed.jpeg",sep="")
output.log( paste("Plot of mean of normalized current",label.c,"against voltage, all cell types in one plot:", filename) )
	jpeg( file=filename, height=4200, width=4200, units="px", quality=75, res=600)
	lt <- length(unique(data$Cell.Type)[good.plots])
	plot( matrix(volt, byrow=TRUE, nrow=nrow(X_norm_plot), ncol=ncol(X_norm_plot)), X_norm_plot, col=CLR, main=paste("Transformed mean of normalized", label.c, "current"), ylab=paste("transformed mean ", label.c, " normalized",sep=""), xlab="Voltage")
	for( i in 1:lt ){
		x <- volt
		y <- X_norm_plot[i,]
		sd <- SEM_plot[i,]
		if(Loess){
			lo <- loess(y~x)
			xl <- seq(min(x), max(x), (max(x)-min(x))/1000)
			lines(xl, predict(lo,xl), col=CLR[i])
		}else{
			lines(matrix(volt, byrow=TRUE, nrow=nrow(X_norm_plot), ncol=ncol(X_norm_plot))[i,], X_norm_plot[i,], col=CLR[i])
		}
		segments(x, y-sd,x, y+sd, col=CLR[i])
		epsilon = 1
		segments(x-epsilon,y-sd,x+epsilon,y-sd, col=CLR[i])
		segments(x-epsilon,y+sd,x+epsilon,y+sd, col=CLR[i])
	}
	legend("bottomleft", legend=unique(data$Cell.Type)[good.plots], lty=1, pch=21, col=CLR)
	dev.off()
	
	if( param.c == "c1"){
		Gnorm <- NULL
		for(i in 1:nrow(X_norm_plot)){
			x <- X_norm_plot[i,]
			pts <- (length(volt)-2):length(volt) #take the last three points
			ModVr <- lm( x[pts] ~ volt[pts])
			Vr <- -(ModVr$coefficients[1])/(ModVr$coefficients[2]) # this is Vr
			G <- x/(volt-Vr)
			# Gmax is based on first derivative. Last point with pos derivative after the first maximum.
			spl <- smooth.spline( volt, G )
			pred <- predict(spl)
			Gprime <- diff(G)/diff(volt)
			mi <- order(Gprime)[length(Gprime)] # max index of Gprime
			Gprimedelta <- Gprime[mi:(length(Gprime)-1)] - Gprime[(mi+1):length(Gprime)]
			if( length(which(Gprimedelta < 0 )) > 0){
				 stp <- min(which(Gprimedelta < 0 ))
				 Gmax <- max(G[1:(mi+stp-1)])
				 }else{ 
				 Gmax <- max(G) }
			gnorm <- G/Gmax
			Gnorm <- rbind(Gnorm, gnorm) #these are the max-normalized G values
			rm(gnorm)
		}
	}
	
	#The Boltzman function
	Blt <- function( x, b, c, d, e){
		y <- c + ((d-c)/((1+exp(b*(x-e)))^1))
		return(y)
	}

	#Integral of Boltzman
	IntBlt <- function( x, b, c, d, e){
		y <- c*x + (d-c)*(x - (1/b)*log(exp(b*x) + exp(b*e)))
		return(y)
	}	

	#Estimate of AUC with trapezoid sums
	AUC <- function( x, y){
		auc <- 0
		for(i in 1:(length(x)-1)){
			auci <- (x[i+1]-x[i])*y[i+1] + 0.5*(y[i+1]-y[i])*(x[i+1]-x[i])
			auc <- c(auc , auci)
		}
		return(auc)
	}	

	#Specify what you are fitting the Boltzmann to
	if( param.c == "c1"){
		DUM <- Gnorm
		DUMlabel <- "Gnorm"
	}
	if( param.c == "c2"){
		DUM <- X_norm_plot 
		DUMlabel <- "C2norm"
	}
	
	#Plot Gnorms & fit Boltzmann
	filename <- paste(dir.plot.out, "/Boltzmann_",label.c,"_vs_volt_celltype_Transformed.jpeg",sep="")
	output.log( paste("Plot of mean of ",DUMlabel,label.c,"against voltage, all cell types in one plot, include Boltmann fit:", filename) )
		jpeg( file=filename, height=4200, width=4200, units="px", quality=75, res=600)
		lt <- length(unique(data$Cell.Type)[good.plots])
		plot( matrix(volt, byrow=TRUE, nrow=nrow(DUM), ncol=ncol(DUM)), DUM, col=CLR, main=paste(DUMlabel), ylab=paste(DUMlabel, label.c), xlab="Voltage")
		Bfit <- NULL
		Bparams <- NULL
		xseq <- seq(volt[1],volt[length(volt)],length.out=1000)
		for( i in 1:nrow(DUM) ){
			x <- volt
			y <- DUM[i,]
			if(Loess){
				lo <- loess(y~x)
				xl <- seq(min(x), max(x), (max(x)-min(x))/1000)
				lines(xl, predict(lo,xl), col=CLR[i])
			}else{
				lines(matrix(volt, byrow=TRUE, nrow=nrow(DUM), ncol=ncol(X_norm_plot))[i,], DUM[i,], col=CLR[i])
			}
			# Fit the Boltzmann function
			Bolt.m1 <- drm(DUM[i,] ~ volt, fct = L.5(fixed = c(NA, NA, NA, NA, 1), names = c("b", "c", "d", "e", "f")))
			Bparams <- rbind( Bparams, c(Bolt.m1$coefficients[1], Bolt.m1$coefficients[2], Bolt.m1$coefficients[3], Bolt.m1$coefficients[4]))
			rownames(Bparams)[i] <- paste(label.c, unique(data$Cell.Type)[good.plots])[i]
			Bfit <- rbind( Bfit, lapply(xseq , Blt, Bolt.m1$coefficients[1], Bolt.m1$coefficients[2], Bolt.m1$coefficients[3], Bolt.m1$coefficients[4]))		
			# show the Boltzmann fit as dotted line on plot
			lines( xseq, Bfit[i,], col=CLR[i] , lty=2)
			# return the estimates of Vhalf(e), and z (b*R*T/F)
			vhf <- signif(Bolt.m1$coefficients[4],3)
			zch <- signif(Bolt.m1$coefficients[1]*Rugc*Temp/Farad,3)
			output.log( paste("Boltmann fit for", unique(data$Cell.Type)[good.plots][i], ": Vhalf=", vhf, "  z=",zch) )
			# save fitted parameters of Boltzmann to file
			filename.tab <- paste(dir.tab.out, "/",label.c,"_BoltmannFit.txt",sep="")
			cat( paste("\n\n", unique(data$Cell.Type)[good.plots][i], "\n"), file=filename.tab, append=T)
			write.table(signif(coef(summary(Bolt.m1)),3), file=filename.tab, row.names=T, col.names=NA, quote=F, sep="\t", append=T)
		}
		leg1 <- paste(unique(data$Cell.Type)[good.plots])
		leg2 <- paste(unique(data$Cell.Type)[good.plots], "Boltzmann")
		if(DUM[1,ncol(DUM)] < DUM[1,1]){pos.leg <- "topright"}else{pos.leg <- "topleft"}
		legend(pos.leg, legend=c(leg1, leg2), lty=rep(c(1,2), each=length(good.plots)), pch=rep(c(21,NA), each=length(good.plots)), col=rep(CLR, 2))
		dev.off()
		
		# Save the Boltzman parameter fits 
		BparamsOUT <- rbind( BparamsOUT, Bparams)			
}

filename <- paste(dir.tab.out, "/BoltzmannFit.txt",sep="")
write.table(signif(BparamsOUT,3), file=filename, row.names=T, col.names=T, quote=F, sep="\t")


#Write a function that takes two Boltzmann fits and returns where they intersect:

CrossOver_points <- function( b1, c1, d1, e1, b2, c2, d2, e2, xl, xr){
	
	x <- seq( xl, xr, by=0.01)
	delta <- NULL
	for(i in 1:length(x)){
		y1 <- Blt( x[i], b1, c1, d1, e1)
		y2 <- Blt( x[i], b2, c2, d2, e2)
		delta[i] <- abs(y1-y2)	
	}
	n <- order(delta)[1]
	xi <- x[n]
	yi <- mean( c(Blt( x[n], b1, c1, d1, e1), Blt( x[n], b2, c2, d2, e2)) )
	return( c(xi, yi) )
}


#Write a function that takes two Boltzmann fits and returns the AUC of the region defined by them:

CrossOver_AUC <- function( b1, c1, d1, e1, b2, c2, d2, e2, xl, xr, xi, yi){

	Boltleft <- order(c(Blt( xl, b1, c1, d1, e1), Blt( xl, b2, c2, d2, e2)))[1]
	if( Boltleft==1){ 
		bl <- b1; br <- b2
		cl <- c1; cr <- c2
		dl <- d1; dr <- d2
		el <- e1; er <- e2
		}else{
			bl <- b2; br <- b1
			cl <- c2; cr <- c1
			dl <- d2; dr <- d1
			el <- e2; er <- e1		
			}
			
	aucl <- abs(IntBlt( xi, bl, cl, dl, el) - IntBlt( xl, bl, cl, dl, el))
	aucr <- abs(IntBlt( xi, br, cr, dr, er) - IntBlt( xr, br, cr, dr, er))
	auctot <- aucr + aucl
	return( auctot )
}


output.log("\nDetermining the intersection points and AUC of the intersection for the two curves for each cell type")
rwnms <- gsub("[-]", "minus", gsub("[+]", "plus", rownames(BparamsOUT)))
lbls <- unique( gsub( "C2 ", "", gsub("C1 ", "", rwnms)))
xl <- volt[1]
xr <- volt[length(volt)]
OUTauc <- NULL
OUTx <- NULL
OUTy <- NULL
for(i in lbls){
	here <- grep( i, rwnms)
	b1 <- BparamsOUT[here[1],1]
	c1 <- BparamsOUT[here[1],2]
	d1 <- BparamsOUT[here[1],3]
	e1 <- BparamsOUT[here[1],4]
	b2 <- BparamsOUT[here[2],1]
	c2 <- BparamsOUT[here[2],2]
	d2 <- BparamsOUT[here[2],3]
	e2 <- BparamsOUT[here[2],4]
	COP <- CrossOver_points( b1, c1, d1, e1, b2, c2, d2, e2, xl, xr)
	xi <- COP[1]
	OUTx <- c(OUTx, xi)
	yi <- COP[2]
	OUTy <- c(OUTy, yi)
	OUTauc <- c( OUTauc, CrossOver_AUC( b1, c1, d1, e1, b2, c2, d2, e2, xl, xr, xi, yi))
}
 
OUTx <- signif(OUTx,3) 
OUTy <- signif(OUTy,3) 
OUTauc <- signif(OUTauc,3) 
for(i in 1:length(OUTauc)){
 	output.log(paste( unique( gsub( "C2 ", "", gsub("C1 ", "", rownames(BparamsOUT))))[i],":", "x=",OUTx[i], "y=", OUTy[i], "AUC=", OUTauc[i]))
}

filename <- paste(dir.tab.out, "/Intersection_properties.txt",sep="")
write.table( cbind(OUTx, OUTy, OUTauc), file=filename, row=unique( gsub( "C2 ", "", gsub("C1 ", "", rownames(BparamsOUT)))) , col=c("x","y","AUC") , quote=F)





##### Analysis of the activation (C1) current density
output.log( "\n### Current density: ###\n" )

param.c <- PARAMS[1] # only do this for C1
label.c <- LABELS[1]
c.col <- get(paste( param.c, ".col", sep=""))
M <- data[,c.col]
cf <- apply( M, 1, min, na.rm=TRUE)
cden <- cf/cap1 # calculate curent density

# write the current density to file
filename <- paste(dir.tab.out, "/",label.c,"currentdensity_well.txt",sep="")
output.log( paste(label.c, "current density in each well:", filename) )
write.table( cden, file=filename, row=names(cden), col=paste(label.c, "_currentdensity",sep=""), sep="\t", quote=F)

# plot histograms of current density, one for each cell type, and fit a normal to each.
xmn <- floor( min( cden, na.rm=TRUE) * 1.1)
xmx <- ceiling( max( cden, na.rm=TRUE) * 1.1)
filename <- paste(dir.plot.out, "/currentdensity_celltype_hist.jpeg",sep="")
jpeg( file=filename, height=4200, width=4200, units="px", quality=75, res=600)
par( mfrow =c(n,n) )
for( i in unique(data$Cell.Type)){
	these <- which( data$Cell.Type == i )
	cnt <- length(these)
	cden.type <- cden[these]
	hist(cden.type, xlab="Current density", main=i, prob=TRUE, xlim=c(xmn, xmx),axes=F, ylab="")
	ft <- fitdistr( cden.type, "normal")
	xn <- seq( min(cden.type), max(cden.type), length.out=1000)
	yn <- dnorm(xn, mean=ft$estimate[1], sd=ft$estimate[2])
	points( xn, yn, yaxt="n" , type="l", col="red")
	par( new=TRUE )
	hist(cden.type, xlab="Current density", ylab="count", main=i, freq=TRUE, xlim=c(xmn, xmx), plot=TRUE)
	mn.rd <- signif(ft$estimate[1],3)
	sd.rd <- signif(ft$estimate[2],3)
	#legend("topleft", leg=as.expression(bquote("N( ",mu,"=",.(mn.rd), ", ",sigma,"=", .(sd.rd),")") ), col="red", lty=1)
	legend("topleft", leg=c( as.expression(bquote(mu == .(mn.rd))), as.expression(bquote(sigma == .(sd.rd))), paste(cnt) ), col="red", lty=c(1, NA, NA), cex=0.75) 
	rm(cden.type)
}
dev.off()

# Calculate mean current density and SEM for each cell type:
CDEN.MN <- NULL; CDEN.SEM <- NULL
for( i in unique(data$Cell.Type)){
	these <- which( data$Cell.Type == i )
	cden.mn <- mean( cden[these], na.rm=TRUE )
	cden.sem <- sd(cden[these], na.rm=TRUE) / sqrt( sum(!is.na(cden[these])) )
	CDEN.MN <- c(CDEN.MN, cden.mn)
	CDEN.SEM <- c(CDEN.SEM, cden.sem)
}

# Barplot of  mean current density and SEM for each cell type:
filename <- paste(dir.plot.out, "/",label.c,"currentdensity_celltype.jpeg",sep="")
output.log( paste("Barplot of mean", label.c, "current density in each cell type:", filename) )
jpeg( file=filename, height=4200, width=4200, units="px", quality=75, res=600)
b1 <- CDEN.MN + CDEN.SEM
b2 <- CDEN.MN - CDEN.SEM
Ylim <- 1.2*max( abs(c(b1,b2)) )
Bars <- barplot( CDEN.MN, ylim=c(0,-Ylim), names=unique(data$Cell.Type), ylab=paste(label.c,"current density"), xlab="Cell type", main=paste("Mean",label.c,"current density"))
segments(Bars, b2, Bars, b1, lwd = 1.5)
arrows(Bars, b2, Bars, b1, lwd = 1.5, angle = 90, code = 3, length = 0.05)
for(i in 1:length(unique(data$Cell.Type))){
	ctype <- unique(data$Cell.Type)[i]
	cnt <- length(which(data$Cell.Type == ctype))
	text(Bars[i], min(c(b1[i],b2[i]))-3 , paste(cnt))
}
dev.off()


# Test for normality of each current density distribution, assuming mean & variance are sample estimates:
# Use Shapiro Wilks test. Null hypothesis is normality.
PVAL <- NULL
for( i in unique(data$Cell.Type)){
	these <- which( data$Cell.Type == i )
	cden.type <- cden[these]
	PVAL<- c( PVAL, signif(shapiro.test(cden.type)$p.value, 3))
}
names(PVAL) <- unique(data$Cell.Type)
IS.NORM <- PVAL >= 0.05 
OUT <- data.frame( PVAL, IS.NORM)
filename <- paste(dir.tab.out, "/",label.c,"currentdensity_normality_test.txt",sep="")
output.log( paste(label.c, "current density normality test:", filename) )
write.table( OUT, file=filename, sep="\t", quote=F, col.names=NA)



# Pairwise comparisons of current density distributions for different cell types:
# Use normality assumption if both pass Shapiro Wilks test for normality
# Use non-parametric KS test if at least one distribution is not normal
lt <- length(unique(data$Cell.Type))
COMP <- matrix(NA, nrow=lt, ncol=lt)
MEAN <- matrix(NA, nrow=lt, ncol=lt)
VAR <- matrix(NA, nrow=lt, ncol=lt)
for( i in 1:lt){
	these <- which( data$Cell.Type == unique(data$Cell.Type)[i] )
	cden.type.i <- cden[these]
	j <- i + 1
	while( j <= lt){
		these <- which( data$Cell.Type == unique(data$Cell.Type)[j] )
		cden.type.j <- cden[these]
		use.Normality <- IS.NORM[i] && IS.NORM[j] 
		if( use.Normality ){
			output.log( paste("Normality assumption applied for comparison of",unique(data$Cell.Type)[i],"and",unique(data$Cell.Type)[j]) )
			#Kolmogorov-Smirnov test for complete distribution
			COMP[i,j] <- signif( ks.test( cden.type.i, cden.type.j, "two.sided")$p.value, 3)
			#Test of equal means
			MEAN[i,j] <- signif( t.test(cden.type.i, cden.type.j, "two.sided")$p.value, 3)
			#Test of equal variances 
			VAR[i,j] <- signif( var.test(cden.type.i, cden.type.j, ratio=1, "two.sided")$p.value, 3)
			
		}else{
			output.log( paste("Normality assumption is NOT applied for comparison of",unique(data$Cell.Type)[i],"and",unique(data$Cell.Type)[j]) )
			#Kolmogorov-Smirnov test for complete distribution
			COMP[i,j] <- signif( ks.test( cden.type.i, cden.type.j, "two.sided")$p.value, 3)
			#Wilcoxon Rank Sum Test for means 
			MEAN[i,j] <- signif( wilcox.test( cden.type.i, cden.type.j, alternative="two.sided")$p.value, 3)
			VAR[i,j] <- signif( fligner.test(list(cden.type.i, cden.type.j))$p.value, 3)			
		}
		j <- j + 1
	}
}
#write conclusions to log file
rownames(COMP) <- unique(data$Cell.Type)
colnames(COMP) <- unique(data$Cell.Type)
IS.SAME <- COMP >= 0.05 
filename <- paste(dir.tab.out, "/",label.c,"currentdensity_same_distr_test.txt",sep="")
output.log( paste(label.c, "current density test for distributions being the same between pairs of cell types:", filename) )
cat("\nP values:\n", file=filename, sep="\n")
write.table( COMP, file=filename, sep="\t", quote=F, append=TRUE, col.names=NA)
cat("\n \nAre the current density distributions the same:\n", file=filename, sep="\n", append=TRUE)
write.table( IS.SAME, file=filename, sep="\t", quote=F, append=TRUE, col.names=NA)

rownames(MEAN) <- unique(data$Cell.Type)
colnames(MEAN) <- unique(data$Cell.Type)
IS.SAME <- MEAN >= 0.05 
filename <- paste(dir.tab.out, "/",label.c,"currentdensity_same_mean_or_median_test.txt",sep="")
output.log( paste(label.c, "current density test for means or medians being the same between pairs of cell types:", filename) )
cat("\nP values:\n", file=filename, sep="\n")
write.table( MEAN, file=filename, sep="\t", quote=F, append=TRUE, col.names=NA)
cat("\n \nAre the current density means or medians the same:\n", file=filename, sep="\n", append=TRUE)
write.table( IS.SAME, file=filename, sep="\t", quote=F, append=TRUE, col.names=NA)

rownames(VAR) <- unique(data$Cell.Type)
colnames(VAR) <- unique(data$Cell.Type)
IS.SAME <- VAR >= 0.05 
filename <- paste(dir.tab.out, "/",label.c,"currentdensity_same_variance_test.txt",sep="")
output.log( paste(label.c, "current density test for variances being the same between pairs of cell types:", filename) )
cat("\nP values:\n", file=filename, sep="\n")
write.table( VAR, file=filename, sep="\t", quote=F, append=TRUE, col.names=NA)
cat("\n \nAre the current density variances the same:\n", file=filename, sep="\n", append=TRUE)
write.table( IS.SAME, file=filename, sep="\t", quote=F, append=TRUE, col.names=NA)


output.log( "\n\n###### Completed ######\n" )

