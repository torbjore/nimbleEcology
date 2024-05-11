####  R Code to calculate likelihoods for Standard, Royle-Link, Multiple-State, and Multiple Method models

expit <- function(mu) {1/((1/exp(mu)) + 1)}
logit <- function(mu) {log(mu/(1-mu))}


###  nLL in each model is the negative log-likelihood
###  INPUT PARAMETERS
###  beta = vector of starting values for betas
###  data = encounter histories
###  counts = number of sites with each of the encounter histories (same length as data)
###  occ = number of occassions (2 values for the TwoSurvey model)


###  Likelihoods for each of the methods for analyses with no covariates for detection

Standard <- function(beta,data,counts,occ){
	nLL = sum(counts*-log((1-expit(beta[2]))*dbinom(data,occ,0)+expit(beta[2])*dbinom(data,occ,expit(beta[1]))))
}

Royle <- function(beta,data,counts,occ){
	if (beta[1] < beta[2]) {beta = beta[c(2,1,3)]}     #### this enforces the constraint that the false positive detection probability will be less than the true positive detection probability
	nLL = sum(counts*-log((1-expit(beta[3]))*dbinom(data,occ,expit(beta[2]))+expit(beta[3])*dbinom(data,occ,expit(beta[1]))))
}

MultipleState <- function(beta,data,counts){
	p = expit(beta[1:3])   ## vector ("true positive","proportion certain","false positive")
	psi = expit(beta[4])
	PI = t(matrix(c((1-p[3]),p[3],0,(1-p[1]),(p[1]*(1-p[2])),(p[1]*p[2])),3,2))
	nLL = c()
	for (a in 1:nrow(data)) nLL = rbind(nLL, counts[a]*-log((1-psi)*dmultinom(data[a,],prob = PI[1,])+(psi)*dmultinom(data[a,],prob = PI[2,])))
	nLL = sum(nLL)
	return(nLL)
}

MultipleMethod <- function(beta,data,counts,occ){
	p = expit(beta[1:3])   ## vector ("detection survey 1 = p11","detection survey 2 = r11","false positive = p01")
	psi = expit(beta[4])
	Y = data[,1]
	W = data[,2]
	PI = p[c(3,1)]
	TAU = c(0,p[2])
	nLL = sum(counts*-log((1-psi)*dbinom(Y,occ[1],PI[1])*dbinom(W,occ[2],TAU[1])+psi*dbinom(Y,occ[1],PI[2])*dbinom(W,occ[2],TAU[2])))
	return(nLL)
}

load("FPdatasets.Rdata")   ### this is a simulated data set with false positive rate of 0.05 and psi = 0.5
###  data for Standard and Royle models are the number of detections at the site (range from 0 to # of occassions)
###  data for multiple state is the number of 0's, 1's, and 2's
###  data for multiple method is the number of detections in the first survey and in the second


outS = optim(betaS,Standard,data = dataS,counts = countS,occ=occS)   # MacKenzie model
expit(outS$par)  ## parameters are p and psi

outR = optim(betaR,Royle,data = dataS,counts = countS,occ=occS)      # The Royle-Link model  (the false positive rate is whichever of the first two betas is smaller - the other is the true positive rate)
expit(outR$par)  ## parameters are p1, p01, and psi

outMS = optim(betaMS,MultipleState,data = dataMS,counts = countMS)
expit(outMS$par)  ## parameters are p1, b, p01, psi

outMM = optim(betaMM,MultipleMethod,data = dataMM,counts = countMM, occ =occMM)
expit(outMM$par)  ## parameters are p1, r, p01, psi


####  Examples of likelihoods with a covariate used for detection (dist) from the analysis in Appendix E

Standard <- function(beta,enc,dist){
	psi <- expit(beta[1])
	P1 <- expit(beta[2] + dist*beta[3])  	### logit link function incorporating distance as a covariate
	nLL = sum(-log((1-psi)*dbinom(enc[,1],enc[,2],0)+psi*dbinom(enc[,1],enc[,2],P1)))
}

Royle <- function(beta,enc,dist){
	psi <- expit(beta[1])
	P1 <- expit(beta[2] + dist*beta[3])		### logit link function incorporating distance as a covariate
	P01 <- expit(beta[4])
	nLL = sum(-log((1-psi)*dbinom(enc[,1],enc[,2],P01)+psi*dbinom(enc[,1],enc[,2],P1)))
}

MultipleMethod <- function(beta,enc,dist){
	psi <- expit(beta[1])
	P1 <- expit(beta[2] + dist*beta[3])		### logit link function incorporating distance as a covariate
	P01 <- expit(beta[4])
	TAU <- expit(beta[5])
	nLL = sum(-log((1-psi)*dbinom(enc[,1],enc[,2],P01)*dbinom(enc[,3],1,0)+psi*dbinom(enc[,1],enc[,2],P1)*dbinom(enc[,3],1,TAU)))
}

MultipleState <- function(beta,enc,dist){
	psi = expit(beta[1])
	P1 <- expit(beta[2] + dist*beta[3])		### logit link function incorporating distance as a covariate
	b <- expit(beta[4])
	P01 <- expit(beta[5])
	nLL = c()
	for (a in 1:nrow(data)) nLL = rbind(nLL,-log((1-psi)*dmultinom(enc[a,],prob = c((1-P01),P01,0))+(psi)*dmultinom(enc[a,],prob = c((1-P1[a]),P1[a]*(1-b),P1[a]*b))))
	nLL = sum(nLL)
}


