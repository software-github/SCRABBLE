svt <- function(data, lambda = NaN, k = 5, incre = 5, m = NaN, n = NaN,
                tol = 1e-10, maxit = 300, method = 'deflation'){

  # define the method used in the svt
  if (strcmpi(method,'deflation')){
    def = 1;
  }else{
    def = 0;
  }

  # set the maximum succession number
  m <- nrow(data)
  n <- ncol(data)
  iter <- min(m,n)
  # check k
  if (k > iter){
    k <- iter
    warning('K is out of dimension, reset to maximum value')
  }

  # check lambda
  if (!(is.nan(lambda))){
    if (is.infinite(lambda)){
      return(0)
    }else{
      if (lambda < 0){
        stop('Lambda must be a nonnegative value.')
      }
    }
  }else{
    return(0)
  }

  # check the increment
  if (!(is.nan(lambda))){
    if (incre > iter - k) {
      incre <- iter - k
    }
  }

  # keep the eigenvectors
  w <- c()

  # keep the eigenvalues
  e <- c()

  # the number of fetching eigenvalues
  n_eig <- 1

  matvec <- function(x,args){
    mv <- c(mSparseT(args, x[(n+1):(m+n)]), mSparse(args, x[1:n]))
    if (length(e) > 0){
      mv <- mv - mMatrix(w, e*(mMatrixT(w, x)))
    }
    return(as.matrix(mv))
  }


  # main loop for computing singular values sequentially
  while(iter > 0){

    # nground <- length(e)
    nground <- 0
    value_svt_tmp <- rARPACK::eigs(matvec, n_eig*k + nground, n = m+n, args = data, which = "LR")
    value_svt_tmp$values <- Re(value_svt_tmp$values)
    value_svt_tmp$vectors <- Re(value_svt_tmp$vectors)
    index <- which(value_svt_tmp$values >= sort(value_svt_tmp$values, decreasing=T)[k], arr.ind=TRUE)
    eigvalues <- value_svt_tmp$values[index]


    if (!(is.nan(lambda))){
      if (value_svt_tmp$nconv < (n_eig*k + nground)){
        warning('The eigenvalues are not convergent and we need refresh with warm start.')
        k <- max(length(e) - 1,6)
        def <- 0
        value_svt_tmp <- rARPACK::eigs(matvec, n_eig*k + nground, n = m+n, args = data, which = "LR")
        value_svt_tmp$values <- Re(value_svt_tmp$values)
        index <- which(value_svt_tmp$values >= sort(value_svt_tmp$values, decreasing=T)[k], arr.ind=TRUE)
        while(value_svt_tmp$nconv < (n_eig*k + nground)){
          warning('The eigenvalues are not convergent and we need refresh with warm start.')
          k <- max(k - 1,1)
          value_svt_tmp <- rARPACK::eigs(matvec, n_eig*k + nground, n = m+n, which = "LR")
          value_svt_tmp$values <- Re(value_svt_tmp$values)
        }
        index <- which(value_svt_tmp$values >= sort(value_svt_tmp$values, decreasing=T)[k], arr.ind=TRUE)
        eigvalues <- value_svt_tmp$values[index]
      }

    }


    if (!(is.nan(lambda))){
      tmp <- which(eigvalues <= lambda)
      vectors <- value_svt_tmp$vectors[, index]
      vectors <- as.matrix(vectors)
      values <- value_svt_tmp$values[index]
      values <- as.vector(values)
      if (!(all(isempty(tmp)))){
        i <- tmp[1]
        if (def == 1){
          w <- cbind(w, vectors[, 1:(i-1)])
          e <- c(e,values[1:(i-1)])
          if (all(isempty(e))){
            return(0)
          }
        }else{
          w <- vectors[, c(1:i-1)]
          e <- values[c(1:i-1)]
          if (all(isempty(e))){
            return(0)
          }
        }
        break
      }
    }


    if(def == 1){
      w <- cbind(w, as.matrix(vectors))
      e <- c(e, as.vector(values))
    }else{
      w <- vectors
      e <- values
    }



    if (is.nan(lambda)){
      break
    }



    iter <- iter - k;


    if (def == 1){
      k <- min(incre, iter)
    }else{
      k <- min(k + incre, iter)
    }


    incre <- 2*incre

  }

  # prepare the output of the whole algorithm
  result <- list()
  result[[1]] <- e
  w <- w*sqrt(2)
  result[[2]] <- w[(n+1):(n+m),]
  result[[3]] <- w[1:n,]
  return(result)
}



