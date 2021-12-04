test_that("multiplication works", {
  simulate.sparse.data.mat<-function(n,m){
    #set.seed(12032020)
    sample.space<-c(rep(0,40),seq(1,10))
    mat<-matrix(sample(sample.space,n*m,replace=T),nrow=n,ncol=m)
    mat<-cbind(sample(c(0,1),n,replace=T),mat)
    colnames(mat)<-c('environment',paste0('gene',seq(1,m)))
    mat<-as.data.frame(mat)
    return(mat)
  }
  beta.diffs<-function(n,m){
    data.sim<-simulate.sparse.data.mat(n=n,m=m)
    glm_form<-as.formula(paste('environment ~ ',paste0('gene',seq(1,m),collapse = '+')))
    Sparse_logistic_reg = sparse_logreg(outcome_column = 'environment',design_columns = paste0('gene',seq(1,m)),data_frame = data.sim)
    GLM_fun = glm(glm_form, family = 'binomial',data=data.sim)
    sparse.coeffs<-cbind(as.vector(Sparse_logistic_reg[,1]),as.vector(GLM_fun$coefficients))
    colnames(sparse.coeffs)<-c('sparse_logistic_regression_estimate','glm_estimate')
    sparse.coeffs<-as.data.frame(sparse.coeffs)
    return(sparse.coeffs)
  }
  sparse.coeffs1 <- beta.diffs(5000,10)
  a1 <- sparse.coeffs1$sparse_logistic_regression_estimate
  b1 <- sparse.coeffs1$glm_estimate
  sparse.coeffs2 <- beta.diffs(5000,20)
  a2 <- sparse.coeffs2$sparse_logistic_regression_estimate
  b2 <- sparse.coeffs2$glm_estimate
  sparse.coeffs3 <- beta.diffs(5000,25)
  a3 <- sparse.coeffs3$sparse_logistic_regression_estimate
  b3 <- sparse.coeffs3$glm_estimate
  sparse.coeffs4 <- beta.diffs(10000,10)
  a4 <- sparse.coeffs4$sparse_logistic_regression_estimate
  b4 <- sparse.coeffs4$glm_estimate
  sparse.coeffs5 <- beta.diffs(10000,20)
  a5 <- sparse.coeffs5$sparse_logistic_regression_estimate
  b5 <- sparse.coeffs5$glm_estimate
  expect_equal(a1, b1)
  expect_equal(a2, b2)
  expect_equal(a3, b3)
  expect_equal(a4, b4)
  expect_equal(a5, b5)
})
