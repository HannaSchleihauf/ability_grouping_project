drop1p<-function(model.res, para=F, data=NULL, contr=NULL, n.cores=c("all-1", "all"), to.del=NULL, return.model.results=F){
	#written Nov. 2018
	if(is.null(contr)){contr=glmmTMBControl()}
	##determine terms that can be dropped:
	xx=as.character(model.res$call)
	names(xx)=names(model.res$call)
  model.fe.re=xx["formula"]
  model.zi=as.formula(xx["ziformula"])
	model.disp=as.formula(xx["dispformula"])
	xfamily=xx["family"]
	if(grepl(xfamily, pattern="(", fixed=T)){
		xfamily=unlist(strsplit(xfamily, split="(", fixed=T))
		xfamily[2]=gsub(x=xfamily[2], pattern=")", replacement="", fixed=T)
		xfamily[2]=gsub(x=xfamily[2], pattern="link = ", replacement="", fixed=T)
		xfamily[2]=gsub(x=xfamily[2], pattern="\"", replacement="", fixed=T)
		xfamily=get(xfamily[1])(xfamily[2])
	}
	##need to address how weights are recognized
	if(any(names(xx)=="weights")){
		data$XXXweights=data[, xx["weights"]]
	}else{
		data$XXXweights=1
	}
	
	xcall=as.character(model.fe.re)
	#xcall=gsub(x=xcall, pattern="\n", replacement=T)
	model.terms=attr(terms(as.formula(xcall)), "term.labels")
	model.terms=model.terms[!grepl(x=model.terms, pattern="|", fixed=T)]
	model.terms=strsplit(model.terms, split=":", fixed=T)
	if(length(to.del)==0){
		to.del=unlist(lapply(1:length(model.terms), function(m1){
			xx=unlist(lapply((1:length(model.terms))[-m1], function(m2){
				length(intersect(model.terms[[m1]], model.terms[[m2]]))==length(model.terms[[m1]])
			}))
			sum(xx)>0
		}))
		model.terms=model.terms[!to.del]
		model.terms=unlist(lapply(model.terms, paste, collapse=":"))
	}else{
		model.terms=to.del
	}
	#if data are handed over
	if(length(data)>0){
		#check whether all columns needed are in them
		#figure out names of the variables involved:
		model.terms2=attr(terms(as.formula(xcall)), "term.labels")
		model.terms2=gsub(x=model.terms2, pattern="||", replacement="|", fixed=T)
		model.terms2=gsub(x=model.terms2, pattern="|", replacement="+", fixed=T)
		model.terms2=gsub(x=model.terms2, pattern=":", replacement="+", fixed=T)
		model.terms2=gsub(x=model.terms2, pattern="*", replacement="+", fixed=T)
		model.terms2=gsub(x=model.terms2, pattern=" ", replacement="", fixed=T)
		model.terms2=gsub(x=model.terms2, pattern="I(", replacement="", fixed=T)
		model.terms2=gsub(x=model.terms2, pattern="^2)", replacement="", fixed=T)
		model.terms2=gsub(x=model.terms2, pattern="(", replacement="", fixed=T)
		model.terms2=gsub(x=model.terms2, pattern=")", replacement="", fixed=T)
		model.terms2=unique(unlist(strsplit(model.terms2, split="+",fixed=T)))
		model.terms2=model.terms2[!model.terms2%in%c("0", "1")]
		#deal with cbind response:
		#guess not needed anymore
# 		if(substr(x=as.character(xcall), start=1, stop=6)=="cbind("){
# 			xx=unlist(strsplit(as.character(xcall), split="~", fixed=T))[1]
# 			xx=gsub(xx, pattern="cbind(", replacement="", fixed=T)
# 			xx=gsub(xx, pattern=")", replacement="", fixed=T)
# 			xx=gsub(xx, pattern=" ", replacement="", fixed=T)
# 			model.terms2=c(model.terms2, unlist(strsplit(xx, split=",", fixed=T)))
# 		}
		#deal with sine/cosine included into model:
		#guess also not needed anymore:
# 		model.terms2=model.terms2[!substr(model.terms2, start=1, stop=4)=="sin("]
# 		model.terms2=model.terms2[!substr(model.terms2, start=1, stop=4)=="cos("]
		
		#still missing: offset term... (weights and squared terms should be fine)
		if(any(!model.terms2%in%names(data))){
			xx=model.terms2[!model.terms2%in%names(data)]
			stop(paste(c(paste(xx, collapse=", "), ifelse(length(xx)==1, "is ", "are "), "missing in data"), collape=""))
		}
	}else{
		stop("Error: no data frame handed over to argument 'data'")
	}
	model.terms2=model.terms[!grepl(x=model.terms, pattern="|", fixed=T)]
	##prepare data:
	#model=paste(model.terms, collapse="+")##get model wrt fixed effects
	ii.data=data
	#all.models=c(xcall, paste(xcall, model.terms, sep="-"))
	all.models=paste(as.character(model.fe.re), model.terms, sep="-")
	##define function doing the model comparison:
	model.fun<-function(i.model, idata, contr., xfamily, xweights){
		i.model=as.formula(i.model)
		i.res=try(glmmTMB(formula=i.model, ziformula=model.zi, dispformula=model.disp, data=idata, family=xfamily, weights=XXXweights, control=contr), silent=T)
		return(i.res)
	}
	##run all models:
  if(para){
    require(parallel)
    n.cores=n.cores[1]
    cl <- makeCluster(getOption("cl.cores", detectCores()))
		if(n.cores=="all"){
			n.cores=length(cl)
		}else if(n.cores=="all-1"){
			if(n.cores=="all-1"){n.cores=length(cl)-1}#else{n.cores=length(cl)}
		}
		n.cores=min(c(n.cores, length(all.models)))
		cl=cl[1:n.cores]
    parLapply(cl=cl, 1:length(cl), fun=function(x){
      return(invisible(library(glmmTMB)))
    })
    x.all.res=parLapply(cl=cl, X=all.models, fun=model.fun, idata=ii.data, contr.=contr, xfamily=xfamily, xweights="xweights")
    parLapply(cl=cl, X=1:length(cl), fun=function(x){invisible(rm(list=ls()))})
    stopCluster(cl)
  }else{
		x.all.res=lapply(X=all.models, FUN=model.fun, idata=ii.data, contr.=contr, xfamily=xfamily, xweights="weights")
  }
  xclass=unlist(lapply(x.all.res, function(x){class(x)[[1]]}))
  failed.models=model.terms[xclass=="try-error"]
  x.all.res=x.all.res[xclass!="try-error"]
	if(length(x.all.res)>0){
		all.tests=lapply(1:length(x.all.res), function(x){
			#xx=as.data.frame(anova(x.all.res[[x]]$value, x.all.res[[1]]$value))
			ires=c(Chisq=-2*(as.vector(logLik(x.all.res[[x]]))-as.vector(logLik(model.res))), "Chi Df"=length(fixef(model.res)$cond)-length(fixef(x.all.res[[x]])$cond))
			ires=c(ires, "Pr(>Chisq)"=as.vector(pchisq(q=ires["Chisq"], df=ires["Chi Df"], lower.tail=F)))
			#ires$AIC=xx[1, "AIC"]
			return(list(logLik=as.vector(logLik(x.all.res[[x]])), AIC=extractAIC(x.all.res[[x]])[2], Chisq=as.vector(ires["Chisq"]), "Chi Df"=as.vector(ires["Chi Df"]), 
				"Pr(>Chisq)"=ires["Pr(>Chisq)"], n=x.all.res[[x]]$modelInfo$nobs,
				conv=x.all.res[[x]]$sdr$pdHess))
		})
		aic=extractAIC(model.res)[2]
		####
		#browser()
		all.tests=data.frame(logLik=c(as.vector(logLik(model.res)), unlist(lapply(all.tests, function(x){x["logLik"]}))),
			AIC=c(aic, unlist(lapply(all.tests, function(x){x["AIC"]}))),
			Chisq=c(NA, unlist(lapply(all.tests, function(x){x["Chisq"]}))),
			"Chi Df"=c(NA, unlist(lapply(all.tests, function(x){x["Chi Df"]}))),
			"Pr(>Chisq)"=c(NA, unlist(lapply(all.tests, function(x){x["Pr(>Chisq)"]}))),
			n=c(model.res$modelInfo$nobs, unlist(lapply(all.tests, function(x){x["n"]}))),
			conv=c(model.res$sdr$pdHess, unlist(lapply(all.tests, function(x){x["conv"]}))))
		rownames(all.tests)=c("none", model.terms[xclass!="try-error"])
		if(length(failed.models)>0){
			xx=matrix(NA, ncol=ncol(all.tests), nrow=length(failed.models))
			colnames(xx)=names(all.tests)
			rownames(xx)=failed.models
			all.tests=rbind(all.tests, xx)
		}
		#all.tests[ncol(all.tests)]=as.character(all.tests[, ncol(all.tests)])
		if(return.model.results){
			to.return=list(drop1.res=all.tests, model.results=x.all.res)
		}else{
			to.return=list(drop1.res=all.tests, model.results=NULL)
		}
	}else{
		to.return=list(drop1.res=NULL, model.results=NULL)
	}
  return(to.return)
}

