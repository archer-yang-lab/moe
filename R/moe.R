moe <- function(x, y, K, lam1 = NULL, lam2 = NULL, 
					maxit = 200, eps = 1e-2, option = TRUE){
    this.call <- match.call()
	# initialization of regularization parameter
    if (!is.matrix(x)) 
        stop("x has to be a matrix")	
    if (any(is.na(x))) 
        stop("Missing values in x not allowed!")
    if (is.null(lam1)) {
      stop("user must provide a lambda1 sequence")
    }
    if (is.null(lam2)) {
      stop("user must provide a lambda2 sequence")
    }
    y <- drop(y)
    np <- dim(x)
    nobs <- as.integer(np[1])
    nvars <- as.integer(np[2])
    vnames <- colnames(x)
    
    if (is.null(vnames)) 
        vnames <- paste("V", seq(nvars), sep = "")
    
    if (length(y) != nobs) 
        stop("x and y have different number of rows")
    
    if (!is.numeric(y)) 
        stop("The response y must be numeric. Factors must be converted to numeric")
	
	phi <- rep(1,K)
	beta <- matrix(rnorm((nvars+1)*K), nvars+1, K)
	alpha <- matrix(rnorm((nvars+1)*K), nvars+1, K)
	
	xx <- cbind(1,x) # design matrix with intercept
	obj <- rep(NA,maxit)
	npass = 0
	while(1){
		#################
		# E-step
		#################
		# prepare to compute the weight w_ik
		# first part: gk
		exp_xa <- exp(xx %*% alpha)
		rowsum_xa <- rowSums(exp_xa)
		gk <- exp_xa/rowsum_xa
		# second part: hk
		xb <- xx %*% beta
		hk <- matrix(NA, nobs, K)
		for(j in 1:K) hk[,j] <- dnorm(y,mean=xb[,j],sd=sqrt(phi[j]))
		# get weight w_ik
		wik <- (gk * hk) / rowSums(gk * hk)
	
		#################
		# M-step
		#################
		# Q_1: weighted lasso problem
		for(j in 1:K){
		# estimate beta	
			fit1 <- glmnet(x=x, y=y, lambda = lam1, alpha = 1, 
							standardize = FALSE,
							family="gaussian", weights=wik[,j])
			beta[,j] <- as.double(rbind2(fit1$a0,fit1$beta))
		# estimate phi
			phi[j] <- sum(wik[,j]*(y-predict(fit1,newx=x))^2) / sum(wik[,j])
		}
	
		if(option == TRUE){
		#################
		# E-step
		#################
		# prepare to compute the weight w_ik
		# first part: gk
		exp_xa <- exp(xx %*% alpha)
		rowsum_xa <- rowSums(exp_xa)
		gk <- exp_xa/rowsum_xa
		# second part: hk
		xb <- xx %*% beta
		hk <- matrix(NA, nobs, K)
		for(j in 1:K) hk[,j] <- dnorm(y,mean=xb[,j],sd=sqrt(phi[j]))
		# get weight w_ik
		wik <- (gk * hk) / rowSums(gk * hk)
		}
		#################
		# M-step
		#################
		# Q_2: multinomial grouped lasso problem with w_ik as working response
		fit2 <- glmnet(x=x, y=wik, lambda = lam2, alpha = 1, 
						standardize = FALSE,
						family="multinomial", type.multinomial="grouped")
		alpha[1,] <- fit2$a0
		for(j in 1:K) alpha[-1,j] <- as.double(fit2$beta[[j]])
	
		#################
		# Check convergence
		#################
		# compute the objective function
	    npass = npass + 1
	    if(npass > maxit) break
		obj[npass] <-  sum(log(rowSums(gk * hk))) 
					- nobs * lam1 * sum(abs(beta)) 
					- lam2 * sum(sqrt(rowSums(alpha * alpha)))
		if(npass>1 && abs(obj[npass]-obj[npass-1])<eps)	break	
	}
	outlist <- list(alpha=alpha, beta=beta, phi=phi, obj=obj[1:npass])
    class(outlist) <- "gglasso"
	outlist
}
