library(devtools)
library(xtable)
library(texreg)
library(R2jags)
library(rjags)
library(runjags)
library(coda)
library(lattice)
library(ggmcmc)

#########
#Bayes functions
#########

post_pred<-function(sampleDF,postDF,invlinkfun){
	"This is a function which allows you to generate posterior predictions
		for all observations in the sample data. 
		It uses the dataset and creates an N*K matrix for the sample data, 
		where N is the num of observations and K is the number of predictions generated
		via the MCMC estimates (e.g. around 2,000).
		
		Inputs:
			sampleDF = data frame with columns per variable
			postDF = posterior estimates for each parameter.
			invlinkfun = link function to transform estimate to outcome. 
				Note you can simply supply a function which just returns 
				the estimate value (e.g. for a linear model).
			
		Notes:
			The order of how sampleDF lists variables must be
			identical to how they were inputted into the JAGS model.
			If theres an intercept, the first column in sampleDF should be 1's."
	
	#sample data
	sd.fin<-sampleDF
	samp.colnames<-colnames(sd.fin)
	sample.mat<-as.matrix(sd.fin)
	Csamp<-dim(sample.mat)[2]
	
	#posteriors
	post.mat<-as.matrix(postDF)
	post.mat.T<-t(post.mat)
	Rpost<-dim(post.mat.T)[1]
	if(Rpost!=Csamp){
		stop("Columns in sampleDF != Rows in postDF. Did you forget an intercept or an interaction term in your sampleDF? Did you leave the deviance in your postDF by accident?")
	}
	
	#Calculate
	post.pred<-sample.mat%*%post.mat.T
	post.pred.p<-invlinkfun(post.pred)
	post.pred.p.quantiles<-apply(post.pred.p, 1, quantile, probs = 
								c(0.025,0.1,0.25,0.5,0.75,0.9,0.975),  
								na.rm = F)
	post.pred.p.q.mat<-as.matrix(post.pred.p.quantiles)
	post.pred.p.q.mat.T<-t(post.pred.p.q.mat)
	
	#Return sample matrix with predictions tacked on at the end
	post.pred.df<-as.data.frame(post.pred.p.q.mat.T)
	sample.pred<-as.data.frame(sample.mat)
	sample.pred<-cbind(sample.pred,post.pred.df)
	colnames(sample.pred)<-c(samp.colnames,c("perc2p5","perc10","perc25","perc50","perc75","perc90","perc97p5"))
	return(sample.pred)	
}

bayes_outtab<-function(fit,caption="Posterior fit with 95% upper/lower",label="tab:output"){
	devtools::source_url("https://raw.githubusercontent.com/jkarreth/JKmisc/master/mcmctab.R")
	out.prep<-mcmctab(as.mcmc(fit))
	out.tab<-xtable(out.prep,caption=caption,label=label)
	return(out.tab)
}

bayes_postDF<-function(fit){
	"This function requires a fitted object that can be represented by the 
	as.mcmc() function."
	.mat<-as.matrix(as.mcmc(fit))
	.postdf<-as.data.frame(.mat)
	.postdf$deviance<-NULL
	return(.postdf)
}

#########
#Inverse links
#########
inv_logit<-function(x){
	"An inverse logit function to transform linear prediction to logit space"
	ret<-(exp(x))/(1+exp(x))
	return(ret)
}

inv_probit<-function(x){
	ret<-pnorm(x)
	return(ret)
}

##########
#JAGS hints
##########

jags_workflow<-function(){
	print("Model: see jags_model(type)")
	print("Data as list")
	print("Params to monitor")
	print("Starting vals as list(list(1),list(2),...list(#chains))")
	print("jags(data,inits,params,n.chains,n.iter,n.burnin,model)")
}

jags_model<-function(type){
	stopifnot(typeof(type)=="character")
	if(type=="ols"){
		cat("
		for(i in 1:N){
			y[i]~dnorm(mu[i],tau)
			mu[i]<-alpha+b*x[i]
		}
		alpha~dnorm(0,0.01)
		b~dnorm(0,0.01)
		tau~dnorm(0,0.01)
		")
	}else if(type=="multilevel"|type=="hierarchical"){
		cat("
		for(i in 1:N){
			y1[i]~dnorm(mu[i],tau)
			mu[i]<-b.1*x.1[i]+b.2[2ndlevel[i]]
		}
		b.1~dnorm(0,0.01)
		tau~dgamma(1,1)
		
		for(j in 1:N.level2){
			b.2[j]~dnorm(b.2.hat[j],tau.level2)
			b.2.hat[j]~a.2+b.2*x.2[j]
		}
		
		a.2~dnorm(0,0.01)
		b.2~dnorm(0.0.01)
		tau.level2~dgamma(1,1)
		")
	}else if(type=="irt"){
		cat("
		#b=discrimination parameters
		#a='difficulty' parameter
		#ipoint=estimate of relative position, e.g. ideal point or ability
		#obs_left=observation number of actor defining left-most extreme
		#obs_right=observation number of actor defining right-most extreme
		for(i in 1:actors){
			for(j in 1:votes){
				y[i,j]~dbern(p[i,j])
				logit(p[i,j])<-ipoint[i]*b[j] - a[j]
			}
		}
		
		ipoint[#obs_left]<- -2.0 #low on left
		ipoint[#obs_right]<- 2.0 #high on right
		
		for(i in 1:(#obs_left-1)){ipoint[i]~dnorm(0,1)} #exclude obs_left
		for(i in (#obs_left+1):(#obs_right-1)){ipoint[i]~dnorm(0,1)} #exclude obs_right
		for(i in (#obs_right+1):actors){ipoint[i]~dnorm(0,1)}
		
		for(j in 1:votes){
			b[j]~dnorm(0,.1)
			a[j]~dnorm(0,.1)
		}
		")
	}else if(type=="ologit"|type=="ordered"){
		cat("
		#Notice that there's CUTPOINTS and CATEGORIES.

		for(i in 1:N){
			for(j in 1:cutpoints){ 
				logit(gamma[i,j])<-theta1[j]-mu[i]
			}

		#Build your staggered distribution
		y[i]~dcat(p[i,1:categories])
		p[i,1]<-gamma[i,1]
		p[i,2]<-gamma[i,2]-gamma[i,1]
		p[i,3]<-gamma[i,3]-gamma[i,2]
		...
		p[i,categories]<-1-gamma[i,(categories-1)]

		#Get your linear prediction
		mu[i]<-b*x[i]

		#Spit out predictions for each category
		#pred[i,1]<-equals(p[i,1],max(p[i,1],p[i,2],p[i,3],...p[i,categories]))
		#pred[i,2]<-equals(p[i,2],max(p[i,1],p[i,2],p[i,3],...p[i,categories]))
		#pred[i,3]<-equals(p[i,3],max(p[i,1],p[i,2],p[i,3],...p[i,categories]))
		#pred[i,4]<-equals(p[i,4],max(p[i,1],p[i,2],p[i,3],...p[i,categories]))
		#predcat[i]<-pred[i,1]+2*pred[i,2]+3*pred[i,3]+...categories*pred[i,categories]
		}

		for(j in 1:cutpoints){
			theta[j]~dnorm(0,1.0E-3)
		}
		theta1[1:cutpoints]<-sort(theta)

		b~dnorm(0,0.01)
		")
	}else if(type=="logit"){
		cat("
		for(i in 1:N){
			y[i] ~ dbern(pi[i])
			logit(pi[i])<-mu[i]
			mu[i] <-a + b*x[i] 
		}

		a~dnorm(0,0.01)
		b~dnorm(0,0.01)
		")
	}else{
		stop("You provided an invalid type, or one I haven't coded up yet.")
	}
}