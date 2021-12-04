## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(Spas)
library(microbenchmark)
library(Matrix)
library(ggplot2)
library(gridExtra)

## -----------------------------------------------------------------------------
data(lsp)
X <- lsp[,3:714]
y_C <- lsp[,1]
y_B <- lsp[,2]
lin_fit <- sparse_linreg(X,y_C)
head(lin_fit)

## -----------------------------------------------------------------------------
# simulates a matrix of sparse covariate data and continuous outcome data,
# where covariate values are drawn from 0-10 with a 4/5
# chance of drawing a 0 and 1/5 chance of drawing a number 1-10,
# and outcome are drawn from 0-50 with uniform probability
simulate.sparse.data.mat.linear<-function(n,m){
  sample.space<-c(rep(0,40),seq(1,10))
  mat<-matrix(sample(sample.space,n*m,replace=T),nrow=n,ncol=m)
  mat<-cbind(sample(c(0,50),n,replace=T),mat)
  colnames(mat)<-c('timepoint',paste0('gene',seq(1,m)))
  mat<-as.data.frame(mat)
  return(mat)
}

# simulates a matrix of sparse covariate data and binary outcome data,
# where covariate values are drawn from 0-10 with a 4/5
# chance of drawing a 0 and 1/5 chance of drawing a number 1-10,
# and outcome are drawn from 1 or 0 with uniform probability
simulate.sparse.data.mat<-function(n,m){
  #set.seed(12032020)
  sample.space<-c(rep(0,40),seq(1,10))
  mat<-matrix(sample(sample.space,n*m,replace=T),nrow=n,ncol=m)
  mat<-cbind(sample(c(0,1),n,replace=T),mat)
  colnames(mat)<-c('environment',paste0('gene',seq(1,m)))
  mat<-as.data.frame(mat)
  return(mat)
}

## -----------------------------------------------------------------------------
accuracy_sim <- function(n, m, times=10) {
  diffs <- numeric(0)
  for(i in 1:times) {
    X <- simulate.sparse.data.mat(n, m)
    sparse_X <- as.matrix(X[,-1])
    sparse_coeffs <- sparse_linreg(sparse_X, X[,1])$Estimate
    lm_coeffs <- lm(as.formula(paste("X[,1] ~ -1 + ", paste0(colnames(X[,-1]), collapse="+"))), data=X)$coefficients
    diffs <- c(diffs, mean(abs(sparse_coeffs - lm_coeffs)))
  }
  return(diffs)
}
# linear accuracy
test <- simulate.sparse.data.mat.linear(50000, 100)
lm_results <- lm(as.formula(paste("test[,1] ~ -1 + ",
                                  paste0(colnames(test[,-1]), collapse="+"))), data=test)
sparse_linreg_results <- sparse_linreg(as.matrix(test[,-1]), test[,1])


accuracy_results <- data.frame("n_obs" = c(500, 5000, 10000, 50000, 50000, 50000, 50000, 50000),
                               "n_feats" = rep(c(10, 10, 10, 10, 15, 20, 30, 50)))

accu_mean_and_sd <- sapply(c(1:nrow(accuracy_results)), function(j) {
  accu <- accuracy_sim(accuracy_results[j,"n_obs"], accuracy_results[j,"n_feats"])
  accu_mean <- mean(accu)
  accu_sd <- sd(accu)
  return(c(accu_mean, accu_sd))
}) #should take several minutes to run
accuracy_results$accu_mean <- accu_mean_and_sd[1,]
accuracy_results$accu_sd <- accu_mean_and_sd[2,]
accuracy_results
coeffs <- data.frame("sparse_coefficients" = sparse_linreg_results$Estimate, "lm_coefficients" = lm_results$coefficients)
ggplot(data = coeffs, aes(x = sparse_coefficients, y = lm_coefficients)) + geom_point() + geom_abline(slope = 1, intercept = 0)

## -----------------------------------------------------------------------------
test <- simulate.sparse.data.mat(50000, 100)
log_fit<-sparse_logreg(outcome_column = 'environment',design_columns = paste0('gene',seq(1,15)),data_frame = test)
log_fit

## -----------------------------------------------------------------------------
#test similarity of estimates:
sparse.est <- sparse_logreg(outcome_column = 'environment',design_columns = paste0('gene',seq(1,15)),data_frame = test)
glm_form<-as.formula(paste('environment ~ ',paste0('gene',seq(1,15),collapse = '+')))
glm.est<-glm(glm_form, family = 'binomial',data=test)
sparse.coeffs<-cbind(as.vector(sparse.est[,1]),as.vector(glm.est$coefficients))
colnames(sparse.coeffs)<-c('sparse_logistic_regression_estimate','glm_estimate')
sparse.coeffs<-as.data.frame(sparse.coeffs)

ggplot(sparse.coeffs,aes(x=sparse_logistic_regression_estimate,y=glm_estimate))+geom_point()+ggtitle('n_obs=50,000; n_features=10')
norm2<-function(x,y){
  diff=x-y
  sum.diff<-0
  for(i in diff){
    sum.diff=sum.diff+i^2
  }
  return(sqrt(sum.diff))
}

beta.diffs<-function(n,m){
  data.sim<-simulate.sparse.data.mat(n=n,m=m)
  glm_form<-as.formula(paste('environment ~ ',paste0('gene',seq(1,m),collapse = '+')))
  Sparse_logistic_reg = sparse_logreg(outcome_column = 'environment',design_columns = paste0('gene',seq(1,m)),data_frame = data.sim)
  GLM_fun = glm(glm_form, family = 'binomial',data=data.sim)
  sparse.coeffs<-cbind(as.vector(Sparse_logistic_reg[,1]),as.vector(GLM_fun$coefficients))
  colnames(sparse.coeffs)<-c('sparse_logistic_regression_estimate','glm_estimate')
  sparse.coeffs<-as.data.frame(sparse.coeffs)
  mean.diff<-mean(abs(sparse.coeffs$sparse_logistic_regression_estimate-sparse.coeffs$glm_estimate))
  sd.diff<-sd(sparse.coeffs$sparse_logistic_regression_estimate-sparse.coeffs$glm_estimate)
  #print(ggplot(sparse.coeffs,aes(x=sparse_logistic_regression_estimate,y=glm_estimate))+geom_point()+ggtitle(paste0('n_obs=50,000; n_features=',m,'; 2-norm=',round(diffs,digits = 3)))+geom_abline(intercept = 0,slope = 1))
  return(c(mean.diff,sd.diff))
}

beta.diffs(50000,10)

diff.mat<-cbind(c(500,5000,10000,50000,10000,10000,10000,10000),c(10,10,10,10,15,20,30,50),rep(NA,8),rep(NA,8))

for(i in 1:nrow(diff.mat)){
  out<-beta.diffs(diff.mat[i,1],diff.mat[i,2])
  print(out)
  diff.mat[i,3]<-out[1]
  diff.mat[i,4]<-out[2]
}

write.table(diff.mat,file = 'beta_diffs.tab',row.names = F,col.names = F,quote = F,sep='\t')

## -----------------------------------------------------------------------------
# using microbenchmark, simulates an n x m sparse data matrix and calculates the computational
# time for running our sparse linear function versus R's lm() function, iterates 100 times
run.feature.sim.linear<-function(n,m){
  X <- simulate.sparse.data.mat.linear(n=n,m=m)
  sparse_X <- as.matrix(X[,-1])
  times<-microbenchmark(
    sparse_linear_regression_results = sparse_linreg(sparse_X, X[,1]),
    lm_results = lm(as.formula(paste("X[,1] ~ -1 + ", paste0(colnames(X[,-1]), collapse="+"))), data=X),
    times = 100
  )
  return(times)
}

# linear time
times_500_10 <- run.feature.sim.linear(500, 10)
times_5000_10 <- run.feature.sim.linear(5000,10)
times_50000_10 <- run.feature.sim.linear(10000,10)
times_500000_10 <- run.feature.sim.linear(50000,10)

times_5000_15 <- run.feature.sim.linear(5000, 15)
times_5000_20 <- run.feature.sim.linear(5000, 20)
times_5000_50 <- run.feature.sim.linear(5000, 50)
plot_it <- function(df, title) {
  df$time <- log(df$time)
  g <- ggplot(data = df, aes(x = expr, y = time))
  g <- g + geom_boxplot()
  g <- g + ggtitle(title)
  g <- g + ylab("log(Time (nanoseconds))")
  g <- g + xlab("function")
  g <- g + scale_x_discrete(labels = c("SLR", "LM"))
  return(g)
}

## -----------------------------------------------------------------------------
features10 <- plot_it(times_5000_10, "10 Features")
features15 <- plot_it(times_5000_15, "15 Features")
features20 <- plot_it(times_5000_20, "20 Features")
features50 <- plot_it(times_5000_50, "50 Features")
grid.arrange(features10, features15, features20, features50, nrow=2)

obs500 <- plot_it(times_500_10, "500 Observations")
obs5000 <- plot_it(times_5000_10, "5,000 Observations")
obs50000 <- plot_it(times_50000_10, "10,000 Observations")
obs500000 <- plot_it(times_500000_10, "500,00 Observations")
grid.arrange(obs500, obs5000, obs50000, obs500000, nrow=2)


## -----------------------------------------------------------------------------
# using microbenchmark, simulates an n x m sparse data matrix and calculates the computational
# time for running our logistic function versus R's glm() function, iterates 100 times
run.feature.sim<-function(n,m){
  data.sim<-simulate.sparse.data.mat(n=n,m=m)
  glm_form<-as.formula(paste('environment ~ ',paste0('gene',seq(1,m),collapse = '+')))
  times<-microbenchmark(
    Sparse_logistic_regression = sparse_logreg(outcome_column = 'environment',design_columns = paste0('gene',seq(1,m)),data_frame = data.sim),
    GLM = glm(glm_form, family = 'binomial',data=data.sim),
    times = 100
  )
  return(times)
}

test1<-run.feature.sim(n=10000,m=10)
test2<-run.feature.sim(n=10000,m=15)
test3<-run.feature.sim(n=10000,m=20)
test4<-run.feature.sim(n=10000,m=50)

test1.2<-test1 %>% as.data.frame() %>% mutate(feat='10 Features')
test2.2<-test2 %>% as.data.frame() %>% mutate(feat='15 Features')
test3.2<-test3 %>% as.data.frame() %>% mutate(feat='20 Features')
test4.2<-test4 %>% as.data.frame() %>% mutate(feat='50 Features')

vary.feat<-rbind(test1.2,test2.2,test3.2,test4.2)
vary.feat$logtime=log(vary.feat$time)
vary.feat$feat = factor(vary.feat$feat, levels=c("10 Features","15 Features","20 Features","50 Features"))

ggplot(vary.feat,aes(x=expr,y=logtime))+geom_boxplot() + ylab('log(Time (nanoseconds))')+xlab('Function') + facet_wrap(~feat) + scale_x_discrete(labels=c("Sparse_logistic_regression" = "SLR", "GLM" = "GLM"))

#modulate number of observations
test.obs1<-run.feature.sim(n=500,m=10)
test.obs2<-run.feature.sim(n=5000,m=10)
test.obs3<-run.feature.sim(n=10000,m=10)
test.obs4<-run.feature.sim(n=50000,m=10)

test.obs1.2<-test.obs1 %>% as.data.frame() %>% mutate(obs='500 Observations')
test.obs2.2<-test.obs2 %>% as.data.frame() %>% mutate(obs='5000 Observations')
test.obs3.2<-test.obs3 %>% as.data.frame() %>% mutate(obs='10000 Observations')
test.obs4.2<-test.obs4 %>% as.data.frame() %>% mutate(obs='50000 Observations')

vary.obs<-rbind(test.obs1.2,test.obs2.2,test.obs3.2,test.obs4.2)
vary.obs$logtime=log(vary.obs$time)
vary.obs$obs = factor(vary.obs$obs, levels=c("500 Observations","5000 Observations","10000 Observations","50000 Observations"))

ggplot(vary.obs,aes(x=expr,y=logtime))+geom_boxplot() + ylab('log(Time (nanoseconds))')+xlab('Function') + facet_wrap(~obs) + scale_x_discrete(labels=c("Sparse_logistic_regression" = "SLR", "GLM" = "GLM"))

