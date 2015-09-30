
##read pvalue data from GitHub
library(maxLik)
library(foreach)
library(sandwich)
library(lmtest)
library(car)
pvalue_data <- read.csv("pvalue_temp2.csv")
attach(pvalue_data)

## pvalue_all contains estimated p-value, pvalue without enough information to back-up is imputed by pvalue_alt (*, ***, <0.05, etc.)
pvalue_all = pval
pvalue_all[is.na(pval)] = pvalue_alt[is.na(pval)]

## code for online infomration
##depending on maxLik version 1.2.4
##set up constraints and loglikelihood of beta-uniform distribution
A4 <- matrix(c(1,-1,0,0,0,0),2,3)
B4 <- matrix(c(0,1),2,1)  
loglik4 <- function(param) {
  ll <- sum ( log( param[1]*dunif(d,beginning,range)
 + (1-param[1]) * dbeta(d,param[2],param[3])
/ (pbeta(range,param[2],param[3])-pbeta(beginning,param[2],param[3])) ) )
  
}
## code to generaet figure 1
png(file="pvalue_figure1_600.png", height=7, width=7, res=600, units='in')
rounding_random_pvalue_all = pvalue_all
rounding_random_pvalue_all[is.na(pvalue)==F]=pvalue[is.na(pvalue)==F]
min(rounding_random_pvalue_all[is.na(pvalue)==F])
beginning = 0;range = 1
d <- rounding_random_pvalue_all[test==1&rounding_random_pvalue_all>beginning &rounding_random_pvalue_all<=range&is.na(rounding_random_pvalue_all)==F]
res <- maxLik(loglik4,  start=c(0.3,0.15,2), constraints=list(ineqA=A4, ineqB=B4),method="BFGS" );print(res)
plot_test = curve( ( (coef(res)[1])*dunif(x,0,range)+(1-coef(res)[1])*dbeta(x,coef(res)[2],coef(res)[3])/pbeta(range,coef(res)[2],coef(res)[3]) ),col="blue",lwd=2, lty=2,n = 10000)
d <- rounding_random_pvalue_all[test==0&rounding_random_pvalue_all>beginning &rounding_random_pvalue_all<=range&is.na(rounding_random_pvalue_all)==F]
res <- maxLik(loglik4,  start=c(0.3,0.15,2), constraints=list(ineqA=A4, ineqB=B4),method="BFGS" );print(res)
plot_control = curve( ( (coef(res)[1])*dunif(x,0,range)+(1-coef(res)[1])*dbeta(x,coef(res)[2],coef(res)[3])/pbeta(range,coef(res)[2],coef(res)[3]) ),col="blue",lwd=2, lty=3,n = 10000)
beginning = 0;range = .15
hist(rounding_random_pvalue_all[test==1],breaks=200,col=rgb(0,0,0,0.5), border=rgb(0,0,0,0),xlim=c(beginning,range),freq=F,main="Distribution of p-values from accounting literature
                (2011 issues of TAR, JAE, JAR, CAR)", ylab="Density",xlab="p-values", ylim=c(0,20),right = T,sub="(assuming rounding-up of p-values)")
hist(rounding_random_pvalue_all[test==0],breaks=200, col=rgb(0,0,0,0.25), border=rgb(0,0,0,.25) ,xlim=c(beginning,range),freq=F, add=T)
abline(v=c(0.01,0.02,0.03,0.04,0.05,.06,.07,.08,.09,.1), col="red", lty=2, lwd=1)
legend("topright", legend=c("Test","Control","Overlap","Test-fitted","Control-fitted"), col=c(rgb(0,0,0,0.5),rgb(0,0,0,0.25),rgb(0,0,0,0.75),"black","black"), lwd=c(10,10,10,2,2), lty=c(1,1,1,2,3))
lines(plot_control,lwd=2, lty=3)
lines(plot_test,lwd=2, lty=2)
dev.off()

## code to generaet statsitics in Table 1

	## z-test
interval=0.005;beginning=0;range=.105
d <- rounding_random_pvalue_all[Audit==1&rounding_random_pvalue_all<=range&rounding_random_pvalue_all>beginning &rounding_random_pvalue_all<=range&is.na(rounding_random_pvalue_all)==F]
length(d)
res <- maxLik(loglik4,  start=c(0.3,0.15,2), constraints=list(ineqA=A4, ineqB=B4),method="BFGS" );print(res)
dist_f <- function(x) {(coef(res)[1])*dunif(x,0,range)+(1-coef(res)[1])*dbeta(x,coef(res)[2],coef(res)[3])/pbeta(range,coef(res)[2],coef(res)[3]) }
count_oe <- function(lower,upper) {
return(c(length(d[d>lower&d<=upper]),integrate(dist_f,lower,upper)$value*length(d)))}
	## estimate for observed,expected counts, and diff in percentage
justbelow = foreach (i=1:10, .combine=rbind) %do% {count_oe(0.005+0.01*(i-1), 0.01+0.01*(i-1))}
justabove = foreach (i=1:10, .combine=rbind) %do% {count_oe(0.01+0.01*(i-1), 0.015+0.01*(i-1))}
cbind(justbelow,justabove)
result = cbind(justbelow,justabove)
result = cbind(result,(result[,1] - result[,2])/result[,2],
			(result[,3] - result[,4])/result[,4],
(result[,1] - result[,2])/result[,2] - (result[,3] - result[,4])/result[,4])
result
justbelow_z = (justbelow[,1] - justbelow[,2])/
		sqrt(justbelow[,2]*(1-rowSums(cbind(justbelow[,1],justabove[,1]))/sum(justbelow,justabove))*
			 rep((1-sum(justbelow[,1])/sum(justbelow,justabove)),nrow(justbelow)))
justbelow_z
justabove_z = (justabove[,1] - justabove[,2])/
		sqrt(justabove[,2]*(1-rowSums(cbind(justbelow[,1],justabove[,1]))/sum(justbelow,justabove))*
			 rep((1-sum(justabove[,1])/sum(justbelow,justabove)),nrow(justabove)))
justabove_z
	##rearrange table
table1z = data.frame(O1 = result[,1], E1 = result[,2], excess1 = result[,5], zstat1 = justbelow_z,
			O2 = result[,3], E2 = result[,4], excess2 = result[,6], zstat2 = justabove_z,diff = result[,7])
table1z
	## t-test for difference in relative excess/shortage below and above critical values
test_twogroup <- function(data,indices) {
interval=0.005;beginning=0;range=0.105

A4 <- matrix(c(1,-1,0,0,0,0),2,3)
B4 <- matrix(c(0,1),2,1)  
loglik4 <- function(param) {
  ll <- sum ( log( param[1]*dunif(d,beginning,range) + (1-param[1]) * dbeta(d,param[2],param[3])/ (pbeta(range,param[2],param[3])-pbeta(beginning,param[2],param[3])) ) )
  ll
}
dist_f <- function(x) {(coef(res)[1])*dunif(x,0,range)+(1-coef(res)[1])*dbeta(x,coef(res)[2],coef(res)[3])/pbeta(range,coef(res)[2],coef(res)[3]) }
count_oe <- function(lower,upper) {
return(c(length(d[d>lower&d<=upper]),integrate(dist_f,lower,upper)$value*length(d)))}

d = data[indices]
res <- maxLik(loglik4,  start=c(0.3,0.15,2), constraints=list(ineqA=A4, ineqB=B4),method="BFGS" )
justbelow = foreach (i=1:10, .combine=rbind) %do% {count_oe(0.005+0.01*(i-1), 0.01+0.01*(i-1))}
justabove = foreach (i=1:10, .combine=rbind) %do% {count_oe(0.01+0.01*(i-1), 0.015+0.01*(i-1))}
result = cbind(justbelow,justabove)
return(
c( (result[,1] - result[,2])/result[,2] - (result[,3] - result[,4])/result[,4] )
)
}
d <- rounding_random_pvalue_all[rounding_random_pvalue_all<=range&rounding_random_pvalue_all>beginning &rounding_random_pvalue_all<=range&is.na(rounding_random_pvalue_all)==F]
all_twogroup_result <- boot(d,test_twogroup,R=10000, parallel="snow", ncpus=10)
all_twogroup_result

## code to estimate statisitcs in Table 2
## setup beta uniform mixture between 0 and 0.05
bum1 <- function(data, indices) {
beginning =0;range = 0.05
  d <-data[indices]
  library(maxLik)
  A4 <- matrix(c(1,-1,0,0,0,0),2,3)
  B4 <- matrix(c(0,1),2,1)  
  loglik4 <- function(param) {
    ll <- sum ( log( param[1]*dunif(d,beginning,range) + (1-param[1]) * dbeta(d,param[2],param[3])/ (pbeta(range,param[2],param[3])-pbeta(beginning,param[2],param[3])) ) )
    ll
  }
  res <- maxLik(loglik4,  start=c(0.3,0.15,2), constraints=list(ineqA=A4, ineqB=B4),method="BFGS" )
  return( res$estimate[1] )
}

		## all p-values
beginning =0;range = 0.05
d <- pvalue_all[pvalue_all>beginning&pvalue_all<range&is.na(pvalue_all)==F];length(d); ## 6015
res <- maxLik(loglik4,  start=c(0.3,0.15,2), constraints=list(ineqA=A4, ineqB=B4),method="BFGS" );print(res)
curve( ( (coef(res)[1])*dunif(x,0,range)+(1-coef(res)[1])*dbeta(x,coef(res)[2],coef(res)[3])/pbeta(range,coef(res)[2],coef(res)[3]) ),add=T,col="blue",lwd=2)
all_all = boot(data=d, statistic=bum1,R=10000, parallel="snow", ncpus=8)
mean(log(size[is.na(pvalue_all)==F&pvalue_all>beginning&pvalue_all<=range]), na.rm=T)
## code to estimate statistics in Table 3
## regression of FDR
bum2 <- function(data, indices) {
library(maxLik)
beginning=0;range=.05
  d <-data[indices,1]
  s <-data[indices,2]
  A4 <- matrix(c(1,-1,0,0,0,0),2,3)
  B4 <- matrix(c(0,1),2,1)  
  loglik4 <- function(param) {
    ll <- sum ( log( param[1]*dunif(d,beginning,range) + (1-param[1]) * dbeta(d,param[2],param[3])/ (pbeta(range,param[2],param[3])-pbeta(beginning,param[2],param[3])) ) )
    ll
  }
  res <- maxLik(loglik4,  start=c(0.3,0.15,2), constraints=list(ineqA=A4, ineqB=B4),method="BFGS" )
  return( c(coef(res),length(d), mean(log(s),na.rm=T) ))
}

bum2n <- function(data, indices) {
library(maxLik)
beginning=0.05;range=1
  d <-data[indices,1]
  s <-data[indices,2]
  A4 <- matrix(c(1,-1,0,0,0,0),2,3)
  B4 <- matrix(c(0,1),2,1)  
  loglik4 <- function(param) {
    ll <- sum ( log( param[1]*dunif(d,beginning,range) + (1-param[1]) * dbeta(d,param[2],param[3])/ (pbeta(range,param[2],param[3])-pbeta(beginning,param[2],param[3])) ) )
    ll
  }
  res <- maxLik(loglik4,  start=c(0.3,0.15,2), constraints=list(ineqA=A4, ineqB=B4),method="BFGS" )
  return( c(coef(res),length(d), mean(log(s),na.rm=T) ))
}
attach(pvalue_data)
beginning=0;range=.05
k=100
pvalue_all
##rounded = rep(0,length(pvalue_all))
##rounded[is.na(pvalue)==F]=1
pval_boot_size = cbind(pvalue_all,pvalue_data$size)

## positive / test

set.seed(4912392);result1 <- boot(data=pval_boot_size[Journal=="TAR"&is.na(pvalue_all)==F&pvalue_all>beginning&pvalue_all<=range,], statistic=bum2,R=k, parallel="snow", ncpus=8);result1
set.seed(4912392);result2 <- boot(data=pval_boot_size[Journal=="JAE"&is.na(pvalue_all)==F&pvalue_all>beginning&pvalue_all<=range,], statistic=bum2,R=k, parallel="snow", ncpus=8);result2
set.seed(4912392);result3 <- boot(data=pval_boot_size[Journal=="JAR"&is.na(pvalue_all)==F&pvalue_all>beginning&pvalue_all<=range,], statistic=bum2,R=k, parallel="snow", ncpus=8);result3
set.seed(4912392);result4 <- boot(data=pval_boot_size[Financial==1&is.na(pvalue_all)==F&pvalue_all>beginning&pvalue_all<=range,], statistic=bum2,R=k, parallel="snow", ncpus=8);result4
set.seed(4912392);result5 <- boot(data=pval_boot_size[Managerial==1&is.na(pvalue_all)==F&pvalue_all>beginning&pvalue_all<=range,], statistic=bum2,R=k, parallel="snow", ncpus=8);result5
set.seed(4912392);result6 <- boot(data=pval_boot_size[Tax==1&is.na(pvalue_all)==F&pvalue_all>beginning&pvalue_all<=range,], statistic=bum2,R=k, parallel="snow", ncpus=8);result6
set.seed(4912392);result7 <- boot(data=pval_boot_size[Audit==1&is.na(pvalue_all)==F&pvalue_all>beginning&pvalue_all<=range,], statistic=bum2,R=k, parallel="snow", ncpus=8);result7
set.seed(4912392);result8 <- boot(data=pval_boot_size[AIS==1&is.na(pvalue_all)==F&pvalue_all>beginning&pvalue_all<=range,], statistic=bum2,R=k, parallel="snow", ncpus=8);result8
set.seed(4912392);result9 <- boot(data=pval_boot_size[Archival==1&is.na(pvalue_all)==F&pvalue_all>beginning&pvalue_all<=range,], statistic=bum2,R=k, parallel="snow", ncpus=8);result9
set.seed(4912392);result10 <- boot(data=pval_boot_size[Experimental==1&is.na(pvalue_all)==F&pvalue_all>beginning&pvalue_all<=range,], statistic=bum2,R=k, parallel="snow", ncpus=8);result10
set.seed(4912392);result11 <- boot(data=pval_boot_size[is.na(pvalue_all)==F&pvalue_all>beginning&pvalue_all<=range&rounded==0,], statistic=bum2,R=k, parallel="snow", ncpus=8);result11
set.seed(4912392);result12 <- boot(data=pval_boot_size[is.na(pvalue_all)==F&pvalue_all>beginning&pvalue_all<=range&rounded==1,], statistic=bum2,R=k, parallel="snow", ncpus=8);result12
set.seed(4912392);result13 <- boot(data=pval_boot_size[is.na(pvalue_all)==F&pvalue_all>beginning&pvalue_all<=range&test==0,], statistic=bum2,R=k, parallel="snow", ncpus=8);result13
set.seed(4912392);result14 <- boot(data=pval_boot_size[is.na(pvalue_all)==F&pvalue_all>beginning&pvalue_all<=range&test==1,], statistic=bum2,R=k, parallel="snow", ncpus=8);result14
set.seed(4912392);result15 <- boot(data=pval_boot_size[is.na(pvalue_all)==F&pvalue_all>beginning&pvalue_all<=range,], statistic=bum2,R=k, parallel="snow", ncpus=8);result15
set.seed(4912392);result16 <- boot(data=pval_boot_size[is.na(pvalue_all)==F&pvalue_all>beginning&pvalue_all<=range&econ_sig==1,], statistic=bum2,R=k, parallel="snow", ncpus=8);result16
set.seed(4912392);result17 <- boot(data=pval_boot_size[is.na(pvalue_all)==F&pvalue_all>beginning&pvalue_all<=range&econ_sig==0,], statistic=bum2,R=k, parallel="snow", ncpus=8);result17
set.seed(4912392);result18 <- boot(data=pval_boot_size[is.na(pvalue_all)==F&pvalue_all>beginning&pvalue_all<=range&(regul==1|practitioners==1),], statistic=bum2,R=k, parallel="snow", ncpus=8);result18
set.seed(4912392);result19 <- boot(data=pval_boot_size[is.na(pvalue_all)==F&pvalue_all>beginning&pvalue_all<=range&(regul==0&practitioners==0),], statistic=bum2,R=k, parallel="snow", ncpus=8);result19
set.seed(4912392);result20 <- boot(data=pval_boot_size[Journal=="CAR"&is.na(pvalue_all)==F&pvalue_all>beginning&pvalue_all<=range&(regul==0&practitioners==0),], statistic=bum2,R=k, parallel="snow", ncpus=8);result20


dim(pval_boot_size[is.na(pvalue_all)==F&econ_sig==1,])
dim(pval_boot_size[is.na(pvalue_all)==F&econ_sig==0,])
dim(pval_boot_size[is.na(pvalue_all)==F,])
dim(pval_boot_size[is.na(pvalue_all)==F&(regul==1|practitioners==1),]) 
dim(pval_boot_size[is.na(pvalue_all)==F&(regul==0&practitioners==0),])

## negative / test
beginning=0.05;range=1
set.seed(4912392);result1n <- boot(data=pval_boot_size[Journal=="TAR"&is.na(pvalue_all)==F&pvalue_all>beginning&pvalue_all<=range,], statistic=bum2n,R=k, parallel="snow", ncpus=8);result1n
set.seed(4912392);result2n <- boot(data=pval_boot_size[Journal=="JAE"&is.na(pvalue_all)==F&pvalue_all>beginning&pvalue_all<=range,], statistic=bum2n,R=k, parallel="snow", ncpus=8);result2n
set.seed(4912392);result3n <- boot(data=pval_boot_size[Journal=="JAR"&is.na(pvalue_all)==F&pvalue_all>beginning&pvalue_all<=range,], statistic=bum2n,R=k, parallel="snow", ncpus=8);result3n
set.seed(4912392);result4n <- boot(data=pval_boot_size[Financial==1&is.na(pvalue_all)==F&pvalue_all>beginning&pvalue_all<=range,], statistic=bum2n,R=k, parallel="snow", ncpus=8);result4n
set.seed(4912392);result5n <- boot(data=pval_boot_size[Managerial==1&is.na(pvalue_all)==F&pvalue_all>beginning&pvalue_all<=range,], statistic=bum2n,R=k, parallel="snow", ncpus=8);result5n
set.seed(4912392);result6n <- boot(data=pval_boot_size[Tax==1&is.na(pvalue_all)==F&pvalue_all>beginning&pvalue_all<=range,], statistic=bum2n,R=k, parallel="snow", ncpus=8);result6n
set.seed(4912392);result7n <- boot(data=pval_boot_size[Audit==1&is.na(pvalue_all)==F&pvalue_all>beginning&pvalue_all<=range,], statistic=bum2n,R=k, parallel="snow", ncpus=8);result7n
set.seed(4912392);result8n <- boot(data=pval_boot_size[AIS==1&is.na(pvalue_all)==F&pvalue_all>beginning&pvalue_all<=range,], statistic=bum2n,R=k, parallel="snow", ncpus=8);result8n
set.seed(4912392);result9n <- boot(data=pval_boot_size[Archival==1&is.na(pvalue_all)==F&pvalue_all>beginning&pvalue_all<=range,], statistic=bum2n,R=k, parallel="snow", ncpus=8);result9n
set.seed(4912392);result10n <- boot(data=pval_boot_size[Experimental==1&is.na(pvalue_all)==F&pvalue_all>beginning&pvalue_all<=range,], statistic=bum2n,R=k, parallel="snow", ncpus=8);result10n
set.seed(4912392);result11n <- boot(data=pval_boot_size[is.na(pvalue_all)==F&pvalue_all>beginning&pvalue_all<=range&rounded==0,], statistic=bum2n,R=k, parallel="snow", ncpus=8);result11n
set.seed(4912392);result12n <- boot(data=pval_boot_size[is.na(pvalue_all)==F&pvalue_all>beginning&pvalue_all<=range&rounded==1,], statistic=bum2n,R=k, parallel="snow", ncpus=8);result12n
set.seed(4912392);result13n <- boot(data=pval_boot_size[is.na(pvalue_all)==F&pvalue_all>beginning&pvalue_all<=range&test==0,], statistic=bum2n,R=k, parallel="snow", ncpus=8);result13n
set.seed(4912392);result14n <- boot(data=pval_boot_size[is.na(pvalue_all)==F&pvalue_all>beginning&pvalue_all<=range&test==1,], statistic=bum2n,R=k, parallel="snow", ncpus=8);result14n
set.seed(4912392);result15n <- boot(data=pval_boot_size[is.na(pvalue_all)==F&pvalue_all>beginning&pvalue_all<=range,], statistic=bum2n,R=k, parallel="snow", ncpus=8);result15n
set.seed(4912392);result16n <- boot(data=pval_boot_size[is.na(pvalue_all)==F&pvalue_all>beginning&pvalue_all<=range&econ_sig==1,], statistic=bum2n,R=k, parallel="snow", ncpus=8);result16n
set.seed(4912392);result17n <- boot(data=pval_boot_size[is.na(pvalue_all)==F&pvalue_all>beginning&pvalue_all<=range&econ_sig==0,], statistic=bum2n,R=k, parallel="snow", ncpus=8);result17n
set.seed(4912392);result18n <- boot(data=pval_boot_size[is.na(pvalue_all)==F&pvalue_all>beginning&pvalue_all<=range&(regul==1|practitioners==1),], statistic=bum2n,R=k, parallel="snow", ncpus=8);result18n
set.seed(4912392);result19n <- boot(data=pval_boot_size[is.na(pvalue_all)==F&pvalue_all>beginning&pvalue_all<=range&(regul==0&practitioners==0),], statistic=bum2n,R=k, parallel="snow", ncpus=8);result19n
set.seed(4912392);result20n <- boot(data=pval_boot_size[Journal=="CAR"&is.na(pvalue_all)==F&pvalue_all>beginning&pvalue_all<=range,], statistic=bum2n,R=k, parallel="snow", ncpus=8);result20n


result_neg = rbind(result1n$t,result2n$t,result3n$t,result4n$t,result5n$t,result6n$t,result7n$t,result8n$t,result9n$t,result10n$t,result11n$t,result12n$t,result13n$t,result14n$t,result15n$t,result16n$t,result17n$t,result18n$t,result19n$t,result20n$t)
result_pos = rbind(result1$t,result2$t,result3$t,result4$t,result5$t,result6$t,result7$t,result8$t,result9$t,result10$t,result11$t,result12$t,result13$t,result14$t,result15$t,result16$t,result17$t,result18$t,result19$t,result20$t)


result_all = data.frame(cbind(result_pos,result_neg))
names(result_all)[1]="FDR";names(result_all)[4]="npos";names(result_all)[6]="NDR";names(result_all)[9]="nneg"
result_all$type1 = result_all$FDR*result_all$npos / (result_all$npos + result_all$nneg)
result_all$type2 = (1-result_all$NDR)*result_all$nneg / (result_all$npos + result_all$nneg)
result_all$power = - result_all$type2 / (1-(result_all$npos*result_all$FDR+result_all$nneg*result_all$NDR)/(result_all$npos+result_all$nneg))+1

summary(result_all)
names(result_all)[5] = "size1";names(result_all)[10] = "size2";

result_all$type = NA
result_all$type[1:k] = "TAR";result_all$type[(k+1):(2*k)] = "JAE";result_all$type[(2*k+1):(3*k)] = "JAR";
result_all$type[(3*k+1):(4*k)] = "Financial";result_all$type[(4*k+1):(5*k)] = "Managerial";result_all$type[(5*k+1):(6*k)] = "Tax";result_all$type[(6*k+1):(7*k)] = "Audit";result_all$type[(7*k+1):(8*k)] = "AIS"
result_all$type[(8*k+1):(9*k)] = "Archival";result_all$type[(k*9+1):(10*k)] = "Experimental"
result_all$type[(k*10+1):(11*k)] = "noround"
result_all$type[(k*11+1):(12*k)] = "round"
result_all$type[(k*12+1):(13*k)] = "control"
result_all$type[(k*13+1):(14*k)] = "test"
result_all$type[(k*14+1):(15*k)] = "All"
result_all$type[(k*15+1):(16*k)] = "sig1"
result_all$type[(k*16+1):(17*k)] = "sig0"
result_all$type[(k*17+1):(18*k)] = "prac"
result_all$type[(k*18+1):(19*k)] = "liter"
result_all$type[(k*19+1):(20*k)] = "CAR"



result_all$TAR = NA;result_all$JAE = NA;result_all$JAR = NA
result_all$Financial = NA;result_all$Managerial = NA;result_all$Tax = NA;result_all$Audit = NA;result_all$AIS = NA;
result_all$Archival = NA;result_all$Experimental = NA;result_all$noround = NA;result_all$round = NA;result_all$control = NA;
result_all$test = NA
result_all$econsig = NA;result_all$pracaud = NA;result_all$noeconsig = NA;result_all$nopracaud = NA;result_all$CAR = NA


result_all[(1:k),15:33] =           data.frame(1,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0)
result_all[(k+1):(2*k),15:33] =     data.frame(0,1,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0)
result_all[(2*k+1):(3*k),15:33] =   data.frame(0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0)
result_all[(3*k+1):(4*k),15:33] =   data.frame(1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1)
result_all[(4*k+1):(5*k),15:33] =   data.frame(1,1,1,0,1,0,0,0,1,1,1,1,1,1,1,1,1,1,1)
result_all[(5*k+1):(6*k),15:33] =   data.frame(1,1,1,0,0,1,0,0,1,1,1,1,1,1,1,1,1,1,1)
result_all[(6*k+1):(7*k),15:33] =   data.frame(1,1,1,0,0,0,1,0,1,1,1,1,1,1,1,1,1,1,1)
result_all[(7*k+1):(8*k),15:33] =   data.frame(1,1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1)
result_all[(8*k+1):(9*k),15:33] =   data.frame(1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1)
result_all[(9*k+1):(10*k),15:33] =  data.frame(1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1)
result_all[(10*k+1):(11*k),15:33] = data.frame(1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1)
result_all[(11*k+1):(12*k),15:33] = data.frame(1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1)
result_all[(12*k+1):(13*k),15:33] = data.frame(1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1)
result_all[(13*k+1):(14*k),15:33] = data.frame(1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1)
result_all[(14*k+1):(15*k),15:33] = data.frame(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)
result_all[(15*k+1):(16*k),15:33] = data.frame(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1)
result_all[(16*k+1):(17*k),15:33] = data.frame(1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1)
result_all[(17*k+1):(18*k),15:33] = data.frame(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1)
result_all[(18*k+1):(19*k),15:33] = data.frame(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1)
result_all[(19*k+1):(20*k),15:33] = data.frame(0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)


head(result_all)

result_all$w1 = NA
for ( i in 1:nrow(result_all)) {
	result_all$w1[i] = sd(result_all$FDR[result_all$type==result_all$type[i]])
}
summary(result_all$w1)


reg = lm(FDR ~ TAR+ JAE + JAR + CAR + Financial + Managerial + Audit + Tax + AIS
 + Archival + Experimental + noround + round + control + test + econsig + noeconsig
 + pracaud + nopracaud + size1 + power, data=result_all, weights=1/(w1^2))
summary(reg)
hc3_result = coeftest(reg,vcov=vcovHC(reg,"HC3"))
linearHypothesis(reg, "round - noround = 0",vcov=vcovHC(reg,"HC3"))
linearHypothesis(reg, "test - control = 0",vcov=vcovHC(reg,"HC3"))

reg = lm(FDR ~ TAR+ JAE + JAR + CAR + Financial + Managerial + Audit + Tax 
 + Archival + Experimental + noround + round + control + test + econsig + noeconsig
 + pracaud + nopracaud + size1 + power, data=result_all[c(1:700,801:2000),], weights=1/(w1^2))
summary(reg)
hc3_result = coeftest(reg,vcov=vcovHC(reg,"HC3"))
linearHypothesis(reg, "round - noround = 0",vcov=vcovHC(reg,"HC3"))
linearHypothesis(reg, "test - control = 0",vcov=vcovHC(reg,"HC3"))



 library(texreg)
htmlreg(list(reg), file = "E:/Project/pvalue/texreg.doc", digits = 3,
override.coef = hc3_result[,1], override.se = hc3_result[,2], 
    override.pval = hc3_result[,4],inline.css = FALSE, doctype = TRUE, html.tag = TRUE, 
    head.tag = TRUE, body.tag = TRUE)
unlink("E:/Project/pvalue/texreg.doc")
