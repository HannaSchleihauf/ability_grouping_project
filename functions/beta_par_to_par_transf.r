beta.par1topar2<-function(shape1=NULL, shape2=NULL, mu=NULL, phi=NULL){
	if(is.null(shape1) & is.null(shape2) & !is.null(mu) & !is.null(phi)){
			shape1=mu*phi
		shape2=(1-mu)*phi
		return(list(shape1=shape1, shape2=shape2))
	}else if(!is.null(shape1) & !is.null(shape2) & is.null(mu) & is.null(phi)){
		mu=shape1/(shape1+shape2)
		phi=shape1+shape2
		return(list(mu=mu, phi=phi))
	}else{
		stop("don't know what to do: either you give me shape1 and shape2 or mu an phi")
	}
}
