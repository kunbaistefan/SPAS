test_that("multiplication works", {
  simulate.sparse.data.mat.linear<-function(n,m){
    sample.space<-c(rep(0,40),seq(1,10))
    mat<-matrix(sample(sample.space,n*m,replace=T),nrow=n,ncol=m)
    mat<-cbind(sample(c(0,50),n,replace=T),mat)
    colnames(mat)<-c('timepoint',paste0('gene',seq(1,m)))
    mat<-as.data.frame(mat)
    sparse_X <- as.matrix(mat[,-1])
    sparse_coeffs <- sparse_linreg(sparse_X, mat[,1])$Estimate
    lm_coeffs <- lm(as.formula(paste("mat[,1] ~ -1 + ", paste0(colnames(mat[,-1]), collapse="+"))), data=mat)$coefficients
    coeffs <- cbind(sparse_coeffs,lm_coeffs)
    return(coeffs)
  }
  a <- simulate.sparse.data.mat.linear(50000,10)
  b <- simulate.sparse.data.mat.linear(50000,20)
  c <- simulate.sparse.data.mat.linear(50000,40)
  d <- simulate.sparse.data.mat.linear(50000,60)
  e <- simulate.sparse.data.mat.linear(50000,100)
  expect_equal(a[,1],a[,2])
  expect_equal(b[,1],b[,2])
  expect_equal(c[,1],c[,2])
  expect_equal(d[,1],d[,2])
  expect_equal(e[,1],e[,2])
})
