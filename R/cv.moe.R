cv.moe <-
  function(x, y, K, lam1 = NULL, lam2 = NULL, nfolds=5, ...) {
    if (is.null(lam1)) {
      stop("user must provide a lam1 sequence")
    }
    if (is.null(lam2)) {
      stop("user must provide a lam2 sequence")
    }
    y <- drop(y)
    y <- as.double(y)
    x <- as.matrix(x)
    nobs <- NROW(x)
    if (nfolds < 3)
      stop("nfolds must be bigger than 3; nfolds=5 recommended")
    foldid <- sample(rep(seq(nfolds), length = nobs))
	cvm_mat <- matrix(NA, length(lam1), length(lam2))
	cvsd_mat <- matrix(NA, length(lam1), length(lam2))
	lam_index <- expand.grid(lam1,lam2)
	lam_posi <- expand.grid(1:length(lam1),1:length(lam2))
	MonteCarlo <- function(i) {		 
	  predvec <- rep(NA, length(y))
	  for (k in seq(nfolds)) {
	    which <- foldid == k
	    y_sub <- y[!which]
	    fit <- moe(x=x[!which, , drop = FALSE], y=y_sub, K=K, lam1=lam_index[i,1], lam2=lam_index[i,2],...)
		xx <- cbind(1,x[which, , drop = FALSE]) # design matrix with intercept 
		exp_xa <- exp(xx %*% fit$alpha)
		rowsum_xa <- rowSums(exp_xa)
		gk <- exp_xa/rowsum_xa
		# second part: hk
		xb <- xx %*% fit$beta
		hk <- matrix(NA, NROW(xx), K)
		for(l in 1:K) hk[,l] <- dnorm(y[which],mean=xb[,l],sd=sqrt(fit$phi[l]))
		# get weight w_ik
		predvec[which] <- log(rowSums(gk * hk))
		cvm <- mean(predvec)
		cvsd <- sd(predvec)
	  }
	  res <- c(cvm, cvsd)	
	  return(res)
	}	
	sim_table <- mclapply(seq(NROW(lam_index)),MonteCarlo,mc.cores = getOption("mc.cores", 4L))

	#     cvm_mat[lam_posi[i,1],lam_posi[i,2]] <- cvm
	#     cvsd_mat[lam_posi[i,1],lam_posi[i,2]] <- cvsd
	# rownames(cvm_mat) <- paste("lam1", seq(length(lam1)), sep = "")
	# colnames(cvm_mat) <- paste("lam2", seq(length(lam2)), sep = "")
	# rownames(cvsd_mat) <- paste("lam1", seq(length(lam1)), sep = "")
	# colnames(cvsd_mat) <- paste("lam2", seq(length(lam2)), sep = "")
    # list(cvm_mat=cvm_mat, cvsd_mat = cvsd_mat)
	sim_table
}
