cv2d.moe <-
  function(x, y, K, lam1 = NULL, lam2 = NULL, ...) {
    if (is.null(lam1)) {
      stop("user must provide a lam1 sequence")
    }
    if (is.null(lam2)) {
      stop("user must provide a lam2 sequence")
    }
    y <- drop(y)
    y <- as.double(y)
    x <- as.matrix(x)
    N <- NROW(x)
    if (missing(foldid))
      foldid <-
      sample(rep(seq(nfolds), length = N))
    else
      nfolds <- max(foldid)
    if (nfolds < 3)
      stop("nfolds must be bigger than 3; nfolds=5 recommended")
    mm.cvm <- Inf
	out_mat <- matrix(NA, length(lam1), length(lam2))
    for (i in seq.int(length(lam1))) {
	  for (j in seq.int(length(lam2))) {
		  
		  
		  
		  
		  
	      cv_out <- moe(x, y, K, lam1 = lam1[i], lam2 = lam2[j], ...)
		  out_mat[i,j] <- cv_out$cvm
	      if (mm.cvm > cv_out$cvm.min) {
	        mm.cvm <- cv_out$cvm.min
	        mm.lambda <- cv_out$lambda.min
	        loc.lam1 <- i
	      }
	      cat("lam1 ", i, " completed.\n")
    	}
	}
	rownames(out_mat) <- paste("S", seq(length(lam1)), sep = "")
	colnames(out_mat) <- paste("L", seq(length(lambda)), sep = "")
    loc.lambda <- which(mm.lambda == lambda)
    list(out_mat=out_mat,
      mm.cvm = mm.cvm, loc.lambda = loc.lambda,
      loc.lam1 = loc.lam1, mm.lambda = mm.lambda,
      mm.lam1 = lam1[loc.lam1]
    )
  }
