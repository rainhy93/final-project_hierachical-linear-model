# Edps/Stat/Psych 587
# c.j.andersion
#
#  Last modified:  3/26/2021
#
#  All of the functions for HLM that I have written
#  and that we use.  More information about functions are 
#  at the top of R for that functions
#
#  Contents:
#
#   icc(model)                                Simple version  
#                                             For more options, use icc R package
#
#   bic.hlm(model,cluster.id)                 Gives BIC that balences between and within
#
#   robust(inmodel,response,idvar,dftype)     Sandwich standard errors and stuff
#
#   contrast(model,L)                         Where L is c x p matrix of contrast 
#                                              coefficients
#    
#   hlmRsq(dataset,model)                     Computes R_1^2 and R_2^2
#
#
#####################################################################
#####################################################################
#####################################################################
icc <- function(modl) { 
   vars <- as.data.frame(VarCorr(modl))[4]
   total <- sum(vars)
   tau00 <- vars[1,1]
   icc <- tau00/total
   
   return(icc)  }
#####################################################################   
#####################################################################
#####################################################################

########################################################################
# Function bic.hlm                                                     #
# Carolyn J Anderson                                                   #
# Feb  25, 2020                                                        #
#                                                                      #
# Based off of results in                                              #
#   https://projecteuclid.org/download/pdfview_1/euclid.ejs/1399035845 #
#   Delattre, M., Lavielle, M., Poursat, M.A. (2014).  A note on BIC in#
#   mixed-effects models.  {\it Electronic Journal of Statistics},     #
#   \it{8}, 456--475.  DOI: 10.1214/140EJS890.                         #
#                                                                      #
# Gives an approximation for BIC for HLM that used 2 sample sizes      # 
#   The number of groups/clusters and the total number of observations #
#   Their derivation used balanced design                              #
#                                                                      #
#  This should be a good approximation for large number of groups and  #
#  large group sizes                                                   #
#                                                                      #
#  I transformed their result such that smaller is  better.  The f     #
#  function computes 2 varients: one based on total sample size and the#
#  other used the harmonical mean of n_j time number of clusters. If   #
#                                                                      #
#  Input lmer object                                                   # 
#                                                                      #
#   bic.hlm(model, cluster.id)                                         #
#                                                                      #
#  Studying the behavoir of this would be a nice project               #        #                                                                      #
########################################################################

bic.hlm <- function(model.fit, cluster.id) {
	N <- ngrps(model.fit)         # number of groups
	n <- nobs(model.fit)          # total number of observations
	varcov <- as.data.frame(VarCorr(model.fit))[4]    
	nrandom <- nrow(varcov)      # number of random effect parameters
	gamma <- fixef(model.fit)
	nfixed <- length(gamma)      # number of fixed effects
	dev <- deviance(model.fit)    # -2lnlikelihood
	
    nj <- table(cluster.id)          # number level 1 units per cluster
    harmonic.mean <- N/sum(1/nj)

    number.parm <- nfixed + nrandom
    aic <- dev + 2*number.parm


    bic.new   <- dev + nrandom*log(N) + nfixed*log(n)                 # New
	bic.harm  <- dev + nrandom*log(N) + nfixed*log(N*harmonic.mean)   # use harmonic mean of n_j
	bic.ngrps <- dev + number.parm*log(N)                             # SAS mixed formula
	bic.ntot  <- dev + number.parm*log(n)                             # lmer formula
	
    results <- cbind(N,n,harmonic.mean, nrandom,nfixed,dev,
	              bic.new,bic.harm,bic.ngrps,bic.ntot,aic)
	colnames(results) <- c('n groups', 'total obs','haroninc.mean.nj', 'n random','n fixed','deviance',
	'bic.new','bic.harm','bic.ngrps','bic.ntot','aic')
    return(results)
}
##############################################################################
##############################################################################
##############################################################################   

###############################################################################
###############################################################################
###############################################################################
# Function robust computed robust/sandwich/empirical standard errors for      #
#   in 2 level linear mixed model (HLM) fixed effects                         #
#                                                                             #
#  version 8                                                                  #
#  March 26,2021                                                              # 
#                                                                             #
#  C.J. Anderson                                                              #  
#  cja@illinois.edu                                                           #  
#                                                                             #
#  **********************IMPORTANT************************                    #
#     For this to run correctly requires that                                 #
#      o the random effects in the formula must be the first ones after "~"   #
#         (at least until I get better are R)                                 # 
#  *******************************************************                    #
#                                                                             #
#   Input:                                                                    #
#      o model0 is the results from lmer                                      #
#      o the response variable                                                #            
#      o the id variable that defines clusters                                # 
#      o dftype="residual" or "between/within" denominater degrees of freedom #
#                                                                             #
#   Output:  A table with                                                     # 
#      o estimtaed fixed effects                                              #
#      o degrees of freedom (residual, between/within, or satterthwaite)      #
#      o model based                                                          #
#           - s.e.                                                            #
#           - t-statistic                                                     #
#           - p-value from Student's t-distribution                           #
#      o robust (empirical, sandwiche estimators)  s.e.                       #
#           - s.e.                                                            #
#           - t-statistic                                                     #
#           - p-value from Student's t-distribution                           #
#                                                                             #       
#  Notes:                                                                     #
#     o For satterthwaite df, you must be using lmerTest                      #
#     o This function computes between/within df differently than SAS when    #
#         discrete variables are treated as.factor; however, when use dummy   #
#         or other coding the df are the same                                 #
#                                                                             #
###############################################################################
###############################################################################
robust <- function(model.fit,response,idvar,dftype) { 

model0 <- model.fit
y <- response
id <- idvar
##########################
# Set up                 #
##########################

varcov <- as.data.frame(VarCorr(model0))[4]   # extract variance covariance as vector
q <- -.5 + sqrt(1/4 + 2*(nrow(varcov)-1))     # number of random effects

nclusters <-length(unique(id))                # number of clusters

X <- model.matrix(model0)                     # Extract design matrix for fixed effects
n <- nrow(X)                                  # total sample size
p <- ncol(X)                                  # number of fixed effects
ncov <- q*(q-1)/2                             # number of covariances 

############################################################################
# This is general and works but perhaps not as efficient as could be       #
############################################################################
if(q==1) { 
   T <- varcov[1,1] 
   Z <- X[,1] 
 } else {
 Z <- X[, 1:q]
 T <- diag(varcov[1:q,])                       
ncov <- q*(q-1)/2
justcov <- varcov[(q+1):(q+ncov),]

T1 <- matrix(,(q), (q))
T1[lower.tri(T1, diag=FALSE)] <- justcov
T2 <- t(T1)
T2[lower.tri(T2, diag=FALSE)] <- justcov
T2 <- as.data.frame(T2)
T2[is.na(T2)] <- 0

T <- T + T2
}
T <- as.matrix(T)

nj <- table(id)                              # number level 1 units per cluster
csum <- cumsum(table(id))                    # cumulative frequencies

cut2 <- csum                                 # end index
cut1 <- c(0, csum) + 1                       # start index
cut1 <- cut1[1:(nclusters)]

sigma2 <- varcov[(q+ncov+1),]                 # put within variance estimated into sigma2

gamma <- as.matrix(fixef(model0), nrow=q, ncol=1)  # extract fixed effects and put into column vector
yhat <-  X %*% gamma                           # model based value of response variable)
yvec <- as.matrix(y, n)                        # turn y into a matrix (vector actually)

model.cov <-  vcov(model0)                     # These are model based ses from lmer
Robust.cov <- matrix(0,nrow=p, ncol=p)

# loop throught to get robust.cov

for (i in 1:nclusters){         
  # This is for getting model based covariance matrix              
    Zj <- X[cut1[i]:cut2[i],1:q]              # extract columns of X for group j  ***********
    Zj <- as.matrix(Zj)
    I <-diag(nrow(Zj))                        # create identity matirx of appropirate size
    Vj <- Zj %*% T %*% t(Zj) + sigma2*I       # compute V_j

    Xj <- X[cut1[i]:cut2[i],]
    iVj <- solve(Vj)
    A  <- model.cov %*% t(Xj) %*% iVj

  # This is for getting data based covaraiance matrix
    yj     <- yvec[cut1[i]:cut2[i],1]          # extract columns of y for group j
    yhatj <- yhat[cut1[i]:cut2[i],1]           # extract columns of yhat for group j
    ssresj <- (yj-yhatj) %*% t(yj - yhatj)     # compute sum sq res for group j

    Robust.cov <- Robust.cov + A %*% ssresj %*% t(A)
}

#################################################################
# Compute test statistics                                       #
#################################################################

model.se <-sqrt(diag(model.cov))
model.se
model.t <- gamma/model.se                      

robust.se <- sqrt(diag(Robust.cov))
robust.se
robust.t <- gamma/robust.se

################################################################
# Compute chosen type of (denominator) df                      #
#   if (var(X[cut1[1]:cut2[1],i]) < .000000001 (i.e., zero)    #
################################################################
rank.design <- rankMatrix(X)
df.residual = n -rank.design[1]

if (dftype=="residual"){ df <- df.residual }        # for residual df 
if  (dftype=="satterthwaite") {                     # for Satterthwaite
     sdf <- summary(model.fit)
	 df <- sdf$coefficients[, 3]
	 }
if (dftype=="between/within"){                      # for between/within df
    pbetween <- 0                                   # find number of between variables
    for (i in 1:p){
    if (var(X[cut1[i]:cut2[i],i]) < .0000001)       # if variance w/in=0, then it's a 
        pbetween <- pbetween + 1                    #   between cluster variable         
                  }

    df <- matrix(, nrow = p, ncol = 1)              # initalize matrix
    for (i in 1:p){
      if (var(X[cut1[1]:cut2[1],i]) >  0.0000){     # checking to see if variance >0
      tmp <- df.residual - nclusters + pbetween     # then this is a within cluster fixed effect
     }   else { tmp <- nclusters - pbetween         # else this is a between cluster fixed effect
     } 
    df[i] <- tmp                          
} 
}

################################################
# Compute p-values                             #
################################################
p.valueR <- 2*(1-pt(abs(robust.t),df))
p.valueM <- 2*(1-pt(abs(model.t),df))

################################################
# Output table of results                      #
################################################
fixed.table <- cbind(gamma, df, model.se, model.t, p.valueM, robust.se, robust.t, p.valueR)
colnames(fixed.table) <- c('Fixed Est.', dftype, 'Model se.','Model t', 'p-value', 'Robust se', 'Robust t', 'p-value')

return(fixed.table)
}
##########################################################################
################ End function robust #####################################
##########################################################################

##########################################################################
##########################################################################
##########################################################################
#  Edps/Psych/Stat 587
#  Carolyn J. Anderson
#
#  Function to do same constrasts as SAS  (used standard multivariate concepts
#    regarding linear combination of normaly distributed random variables)
#
#  Ho:    L*Gamma = 0               Ha:  L*Gamma <> 0
#
#  Where L is an (c X p) matirx of c contrast on the p fixed effects
#        Gamma is (p X 1) vector of fixed effects
#
#
#  Input: model
#         matrix (vector) with number of columns=number of fixed effects
#                   and elements are constrast coefficients (sum to 0)
#
#  Output:
#          F statistics
#          Guess as df
#          p-value for F
#          X2  (chi-square test statistic)
#          df
#          p-value for X2
#
#  Requirements:  Model should be fit by lmer with lmerTest (the latter is used
#                     as a guess at den df for F).
#
#  I am currently working on options for ddfm and empirical/robus/sandwich se's 
#
#####################################################################################

contrast <- function(model,L) {
 
  gamma <- fixef(model)
  nfixed <- length(gamma)                          # number of fixed effects
  ftable <- summary(model)[10]                     # use this for df (satterthwaite)
  ftable <- matrix(unlist(ftable),nrow=nfixed)     # so I can put out df for fixed effects
  cov.gamma <- vcov(model)                         # covariance matrix of parameter estimates
 
# check that number of columms of L = number of fixed effects
  nL <- ncol(L)
  if (nfixed != nL) { 
     print('The number of columns of L must equal number of fixed effects.')
     return( )
   }

# check for linearly dependent rows of L
  ck <- L %*% t(L)
  det.ck <- det(ck)
  if (det.ck <= .0000001) {
     print('Check for linearly dependent rows of L')
     return()
   }

  estimate<- L %*% gamma                           # value of constrast(s)
                                                   # F test statistic
  F <-  as.numeric((t(estimate) %*% solve(L %*% cov.gamma %*% t(L)) %*% estimate)/nrow(L))
                                                   # Wald test statistic
  X2 <-  as.numeric((t(estimate) %*% solve(L %*% cov.gamma %*% t(L)) %*% estimate))


  which.df <- which(L !=0, arr.ind=TRUE)           # Find elements of L that are non-zero

  ddf <- ftable[which.df[1,2],3]                   # I don't know what these should be                
  ndf <- nrow(L)                                   # these are fine

  pvalue <- pf(F,ndf,ddf,lower.tail=FALSE)         # p-value of F
  pchi  <- pchisq(X2,ndf,lower.tail=FALSE)         # p-value for Wald
                                                   # output results in nice table/format
  result <- c(F, ndf, ddf, pvalue, X2, ndf, pchi)
  names(result) <- c('F', 'num df', 'den df guess', 'p-value', 'Wald', 'df', 'p-value')
  return(result)

}
##########################################################################
################### End function contrast  ###############################
##########################################################################

   
#####################################################################################
#####################################################################################
#####################################################################################
#  EdPsy/Stat/Psych 587
#  Carolyn J. Anderson  
#
#####################################################################################
#####################################################################################
#    Function 'hlmRsq'           (good one)                                         #
#                                                                                   #
#    Rsq are computed by the code below equal the proportional reduction of         #
#      prediction errorvariance accounted for the level 1 and Level 2 Models.       #
#                                                                                   #
#    In population these will always be positive.  If they are negative, this       #
#      is indicates model miss-specfication.  If they get smaller as add variables, #
#      this also could be model miss-specification                                  #
#                                                                                   #
#    There were proposed by Snijders, TAB, Bosker, RJ (1994). "Modeled  variance    #
#      in  two-level  models." Sociological Methods & Research, 22(3), 342-363.     #
#                                                                                   #
#  Input (arguments):                                                               #
#    "dataset" to used when fitting the model                                       #
#    "model1"  is the model for which you want R1sq and R2sq                        #
#                                                                                   #
#  Requirements:                                                                    #
#    o variable identifying group/cluster should be "id"                            #
#    o first (set of fixed) variable(s) in formula are the random ones              #
#                                                                                   #
#####################################################################################
#####################################################################################

hlmRsq <- function(dataset,model1) {

X <- model.matrix(model1)                   # Extract design matrix for fixed effects
n <- nrow(X)                                # total sample size
varcov <- as.data.frame(VarCorr(model1))[4] # extract variance covariance of Random effects as vector
q <- -.5 + sqrt(1/4 + 2*(nrow(varcov)-1))   # number of random effects (solution to quadratic eq)
ncov <- q*(q-1)/2                           # number of covariances 
justcov <- varcov[(q+1):(q+ncov),]          # vector of just the covariances
nclusters <- length(unique(id))             # number of clusters

T <- if(q==1) { 
   Z <- X[,1]
   zbar <- mean(Z)
   T <- varcov[1,1]
 } else{ 
   Z <- X[, 1:q]
   zbar <- colMeans(Z)
   T <- diag(varcov[1:q,])                       
   T1 <- matrix(,(q), (q))
   T1[lower.tri(T1, diag=FALSE)] <- justcov
   T2 <- t(T1)
   T2[lower.tri(T2, diag=FALSE)] <- justcov
   T2 <- as.data.frame(T2)
   T2[is.na(T2)] <- 0
   T <- T + T2
}
T <- as.matrix(T)                       # Known at Psi in my glmm book

# Set up for loop

nj <- table(id)                         # number level 1 units per cluster
csum <- cumsum(table(id))               # cumulative frequencies
cut2 <- csum                            # end index
cut1 <- c(0, csum) + 1                  # start index
cut1 <- cut1[1:(nclusters)]

nj <- table(id)                         # number level 1 units
Sw <- matrix(0, nrow=q, ncol=q)         # initial Sum of Square between
Sb <- matrix(0, nrow=q, ncol=q)         # initial sum of squares within

# loop throught to get Sb and Sw

for (i in 1:nclusters){         
  # This is for getting model based covariance matrix              
    Zj <- X[cut1[i]:cut2[i],1:q]        # extract columns of X for random design for group j
    Zj <- as.matrix(Zj)
    onej <- matrix(1, nrow=1, ncol=nj[i])
    zbarj <- (t(onej %*% Zj)/nj[i])
    zbar  <- matrix(zbar, nrow=q, ncol=1)
    Sb <- Sb + nj[i] * (zbarj - zbar) %*% t(zbarj - zbar) 

    zbarj <- t(matrix(zbarj, nrow= q, ncol=nj[i]))
    Sw <- Sw +  t(Zj - zbarj) %*% (Zj - zbarj)
}
Sb <- Sb/(n-1)
Sw <- Sw/(n-1)

sigma2 <- varcov[(q+ncov+1),]            # put within variance estimated into sigma2

#
# Fit the null model using same estimation method at model input
#

# Determines whether ML or REML was used to fit the model

sum <- summary(model1)
if (sum[1]=="Linear mixed model fit by maximum likelihood "){
   method="FALSE"
} else {
   method='TRUE'
}

# Fits the null model by updating formula from model input to functon

f <- formula(model1)
all <- all.vars(f)
null <- update(f, .~ 1 -all + (1 | id))
model0 <- lmer(null, data=dataset, REML=method)
varcov0 <- as.data.frame(VarCorr(model0))[4] 

####################################################################
#  R1sq    & R2sq                                                  #
####################################################################

numerator1 <- t(zbar) %*% T %*% zbar + sum(diag( (T %*% (Sb + Sw)) )) + sigma2
denominator1 <- sum(varcov0)
R1sq <- 1 -numerator1/denominator1

harmonic.mean <- nclusters/sum(1/nj)

numerator2 <- t(zbar) %*% T %*% zbar + sum(diag((T %*% (Sb + Sw/harmonic.mean)))) + sigma2/harmonic.mean
denominator2 <- sum(varcov0[1,1]) + varcov0[2,1]/harmonic.mean

R2sq <- 1 - numerator2/denominator2

Rsq.values <- cbind(harmonic.mean, R1sq,R2sq)
colnames(Rsq.values) <- c('harmonic.mean', 'R1sq', 'R2sq')
return(Rsq.values)
}
####################################################################
##################### end function hlmRsq ##########################
####################################################################

   