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
eucdist<-function(x,y){
	#Calculate euclidean distance in N-dimensional space from x to y, where x and y are dataframes or matrices.
	.matX<-as.matrix(x)
	.matY<-as.matrix(y)
	.matDiff<-(.matX-.matY)^2
	.eucdistSQ<-rowSums(.matDiff)
	.eucdist<-sqrt(.eucdistSQ)
	return(as.vector(.eucdist))
}

#####################################
#Panel lagging (For use with "ddply")
#Syntax: "ddply([dataframe],~[panelvar],transform,[newvar]=lg([oldvar],[lvalue]))"
######################################
check<-function(h,lvalue){length(h)-lvalue}
lag<-function(f,lvalue){c(f[1:lvalue],f[1:check(f,lvalue)])}
nolag<-function(g){g}
lg<-function(x,lvalue){if(check(x,lvalue)>0){lag(x,lvalue)}else{nolag(x)}}
