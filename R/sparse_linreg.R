# inputs: covariate matrix X, ideally in sparse matrix format, and outcome vector y
# outputs: dataframe of coefficient estimates, standard error, T statistic for
# each estimate, and its associated p-value

#'sparse_linreg
#'
#'Conduct linear regression on sparse matrix
#'@import tidyverse
#'@param X covariate matrix X
#'@param y outcome vector y
#'@return This function would return a data frame
#'@return \code{Estimate}: the estimated coefficients
#'@return \code{SE}: standard errors
#'@return \code{T}:  T statistics
#'@return \code{p_value}:  p_value of the corresponding statistics
#'@examples
#' #One dataset from LSQ
#' data(lsp)
#' X <- lsp[,3:714]
#' y_C <- lsp[,1]
#' lin_fit <- sparse_linreg(X,y_C)


#'@export


sparse_linreg <- function(X, y) {
  n <- nrow(X)
  p <- ncol(X)
  U <- chol(t(X) %*% X)
  Uinv <- backsolve(U, diag(1, nrow = p))
  XtX_inv <- Uinv %*% t(Uinv)
  B <- XtX_inv %*% t(X) %*% y
  XB <- X %*% B
  sq_err <- as.numeric(t(y-XB) %*% (y-XB))
  dof <- n-p
  se <- sqrt((sq_err/dof)*diag(XtX_inv))
  Ts <- B/se
  prs <- 2*(1-sapply(Ts, pt, dof))
  output <- data.frame("Estimate" = as.matrix(B), "SE" = as.matrix(se),
                       "T" = as.matrix(Ts), "p_value" = as.matrix(prs))
  return(output)
}
