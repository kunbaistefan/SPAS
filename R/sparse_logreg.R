## Function to perform logistic regression given a data frame containing a column
# corresponding to a binary outcome and other columns corresponding to covariates.
# Arguments are the (1) outcome_column, which is the column name of the binary outcome,
# (2) design_columns = the column names corresponding to the covariates, and (3)
# data_frame, which is the name of the data frame
# output is a matrix with same number of rows the length of the design_columns argument
# (corresponding to β estimates) and 4 columns. For each covariate, the estimate,
# standard error, z score, and p value are given

#'sparse_logreg
#'
#'Conduct logistic regression on sparse matrix
#'@import tidyverse
#'@param outcome_column the column name of the binary outcome
#'@param design_columns the column names corresponding to the covariates
#'@param data_frame the name of the data frame
#'@return A matrix with same number of rows the length of the design_columns argument (corresponding to β estimates) and 4 columns.
#'@return \code{Estimate}: the estimated coefficients
#'@return \code{SE}: standard errors
#'@return \code{T}:  T statistics
#'@return \code{p_value}:  p_value of the corresponding statistics
#'@examples
#'# simulate mateix
#' set.seed(12032020)
#' sample.space<-c(rep(0,40),seq(1,10))
#' mat<-matrix(sample(sample.space,n*m,replace=T),nrow=50000,ncol=100)
#' mat<-cbind(sample(c(0,1),n,replace=T),mat)
#' colnames(mat)<-c('environment',paste0('gene',seq(1,100)))
#' text<-as.data.frame(mat)
#' log_fit<-sparse_logreg(outcome_column = 'environment',design_columns = paste0('gene',seq(1,15)),data_frame = test)

#'@export

sparse_logreg<-function(outcome_column,design_columns,data_frame){
  n=nrow(data_frame)
  Y_star=data_frame[,which(colnames(data_frame)==outcome_column)]
  Y_star = matrix(Y_star, nrow = n, 1)
  X= model.matrix(as.formula(paste('~',paste(design_columns,collapse = '+'))),as.data.frame(data_frame))
  q=ncol(X)
  beta=rep(0,q)
  #implement iteratively reweighted least squares
  tol=0.00001
  epsilon=99
  ite_max=25
  ite=0
  while (epsilon > tol & ite <= ite_max){
    eta=X %*% beta
    mu=1/(1+exp(-1*eta))
    nu=mu*(1-mu)
    V = c(nu)
    Z = eta + 1/(nu) * (Y_star-mu)
    VX = V * X
    VZ = V * Z

    beta_new=solve(crossprod(X, VX), crossprod(X, VZ))
    epsilon = as.numeric(sqrt(crossprod(beta_new-beta,beta_new-beta)))
    beta=beta_new
    ite=ite+1
    beta_t=t(beta)
  }
  se<-sqrt(diag(solve(crossprod(X,VX))))
  z_scores<-beta_t/se
  p_values<-2*pnorm(abs(as.numeric(z_scores)),lower.tail = F)
  out<-cbind(as.vector(beta),as.vector(se),as.vector(z_scores),as.vector(p_values))
  colnames(out)<-c('Estimate','SE','Z','p_value')
  rownames(out)<-colnames(X)
  return(out)
}
