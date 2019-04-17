#########
#Misc Data functions
#########

########
#Rescaling to a particular mean/sd
########
rescale<-function(olddata,newmean,newsd){
	.oldmean<-mean(olddata)
	.oldsd<-sd(olddata)
	.t1<-newmean
	.t2<-olddata-.oldmean
	.t3<-newsd/.oldsd
	.y<-.t1+(.t2*.t3)
	return(.y)
}

###########
#Calculate euclidean distance
###########
cossim<-function(x,y){
	#########################
	#Calculate cosine simimilarity in N-dimensional space from X to Y,
	#	where x and y are dataframes or matrices
	#########################
	if( typeof(x)=="double" & length(x)>1 & is.null(dim(x)) ){
		warning("Treating first input as single vector")
		.matX<-t(matrix(x))
	}else{
		.matX<-as.matrix(x)
	}
	
	if( typeof(y)=="double" & length(y)>1 & is.null(dim(y)) ){
		warning("Treating second input as single vector")
		.matY<-t(matrix(y))
	}else{
		.matY<-as.matrix(y)
	}
	
	.num<-rowSums(.matX*.matY)
	.den<-sqrt(rowSums(.matX^2))*sqrt(rowSums(.matY^2))
	.out<-.num/.den
	names(.out)<-NULL
	return(.out)
}

eucdist<-function(x,y){
	################
	#Calculate euclidean distance in N-dimensional space from x to y, 
	#	where x and y are dataframes or matrices.
	################
	if( typeof(x)=="double" & length(x)>1 & is.null(dim(x)) ){
		warning("Treating first input as single vector")
		.matX<-t(matrix(x))
	}else{
		.matX<-as.matrix(x)
	}
	
	if( typeof(y)=="double" & length(y)>1 & is.null(dim(y)) ){
		warning("Treating second input as single vector")
		.matY<-t(matrix(y))
	}else{
		.matY<-as.matrix(y)
	}

	stopifnot(dim(.matX)==dim(.matY))
	.matDiff<-(.matX-.matY)^2
	.eucdistSQ<-rowSums(.matDiff)
	.eucdist<-sqrt(.eucdistSQ)
	return(as.vector(.eucdist))
}

#####################################
#Panel lagging (For use with "plyr" package)
#Syntax: "ddply([dataframe],~[panelvar],transform,[newvar]=lg([oldvar],[lvalue]))"
######################################
check<-function(h,lvalue){length(h)-lvalue}
lag<-function(f,lvalue,NAinit=T){if(NAinit){c(rep(NA,lvalue),f[1:check(f,lvalue)])}else{c(f[1:lvalue],f[1:check(f,lvalue)])}}
nolag<-function(g){g}
lg<-function(x,lvalue,NAinit=T){if(check(x,lvalue)>0){lag(x,lvalue,NAinit)}else{nolag(x)}}