#############################################################################
#Dependencies:
#############################################################################
library(dplyr)   #tidy code

#############################################################################
#Format initial state conditions
#############################################################################

###### TO DO - check this description  -BHR ###########################

#Input: init.df
#       -A data frame with columns being the disease states (SIR etc) and the 
#        rows being the spatial subdivisions (districts, grid cells etc)

#Output: a named vector that numbers the sub-states (S1, S2, ... Sn, I1....In, etc)
#       -This vector can be used as the initial state vector for the deSolve ODE solver

#----------------------------------------------------
# #Example: 
# init <- data.frame(S=c(100,100,100,100),
#                    E= c(0,0,0,0),
#                    I= c(10,0,0,0),
#                    V= c(0,0,0,0),
#                    H= c(0,0,0,0)) %>%
#   rfun.FormatInit()
# 
# #Return: 
# init
# > S1  S2  S3  S4  E1  E2  E3  E4  I1   I2  I3  I4  V1  V2  V3  V4  H1  H2  H3  H4 
# > 100 100 100 100 0   0   0   0   10   0   0   0   0   0   0   0   0   0   0   0 
#----------------------------------------------------

rfun.FormatInit <- function(init.df){
  
  #create flexible naming structure for vectors
  num_patch <<- nrow(init.df)
  names_states <- as.vector(colnames(init.df))
  num_state <<- length(names_states)
  num_comps <<- num_state*num_patch
  
  #Carrying capacity
  K <<- rowSums(init.df)
  
  varnames <- NULL
  for(i in 1:num_state){
    varnames <- cbind(varnames, paste(names_states[i], seq(1:num_patch), sep=""))
  }
  as.vector(varnames)
  
  #vectorize matrix input
  init.vector <- as.data.frame(t(as.vector(as.matrix(init.df))))
  colnames(init.vector) = varnames
  
  return(init.vector)
}



#############################################################################
#Format contact matrix 
#############################################################################
#Input: df.contact, beta.contact_constant
#       -df.contact: A data frame with both rows and columns being spatial 
#        subdivisions (districts, grid cells) and the entries being the intersections
#       -beta.contact_constant: transmission coefficient scalar

#Output: a matrix of the transmission coefficients between subdivisions

#----------------------------------------------------
# #Example:
# beta.contact_constant = 0.0001
# df.contact <- data.frame(q1 = c(0,1,1,0),
#                          q2 = c(1,0,0,1),
#                          q3 = c(1,0,0,1),
#                          q4 = c(0,1,1,0))
# beta.contact = rfun.FormatContact(df.contact, beta.contact_constant)
# 
# #Return:
# beta.contact
# >      q1    q2    q3    q4
# > [1,] 0e+00 1e-04 1e-04 0e+00
# > [2,] 1e-04 0e+00 0e+00 1e-04
# > [3,] 1e-04 0e+00 0e+00 1e-04
# > [4,] 0e+00 1e-04 1e-04 0e+00
#----------------------------------------------------

rfun.FormatContact = function(df.contact, beta.contact_constant){
  l = length(df.contact)
  matrix.contact  = as.matrix(df.contact, ncol=l)
  return(beta.contact_constant*matrix.contact)
}





#############################################################################
#Model 1: PatchSEIV
#############################################################################

###### TO DO - update this description for stochastic model -BHR ######## 

#Input: times, state, parameters_full
#       -times: a vector of time points to evaluate the model. Unit = days
#       -state: disease states, specified by initial conditions and model set up
#       -parameters_full: a list that includes the following fields:
#                -signal: a matrix that has 1 column for time series and n columns
#                 for n spatial subdivisions. The rows are time points. The first 
#                 column values are the time points, the others are the instantaneous
#                 vaccination rates at that time for that subdivision
#                -parameters: a named vector (c(beta1=xxx, beta2=xxx, ...., alpha1 = xxxx, ...)
#                 of parameter values
#                -beta.contact: a matrix of transmission coefficients between spatial subdivisions
#
#Output: SpatialPulseVax1
#        -a function that can be fed into deSolve that specifies canine rabies (SIEV) model with pulsed
#        -vaccination
#----------------------------------------------------
# #Example:
# source("Functions_InputFormatting.R")
# 
# #times
# times <- seq(0, 365, by=1)
# 
# #state (initial conditions)
# init <- data.frame(S=c(100,100), E= c(0,0), I= c(10,0), V= c(0,0), H= c(0,0)) %>%
#   rfun.FormatInit()
# 
# #parameters_full
# parms <- data.frame(beta = c(0.01,0.01), mu = c(0.001, 0.001), gamma = c(0.2, 0.2),
#                     alpha = c(0.5, 0.05), nu1   = c(0, 0),N = c(100, 100))%>%
#   rfun.FormatParams() #model parameters
# beta.contact = rfun.FormatContact(data.frame(q1 = c(0,1), q2 = c(1,0)), 0.0001) #spatial transmission matrix
# signal <- matrix(c(times, rep(0, (nrow(beta.contact)*length(times)))), ncol=(nrow(beta.contact)+1)) #signal for pulsevax
# parameters_full <- list(signal = signal, parameters = parms, beta.contact = beta.contact)
# 
# #Solve ode
# out <- ode(y=init, times=times, func=SpatialPulseVax1, parms=parameters_full)
# 
# #Return:
# head(out)
# >     time  S1        S2        E1         E2        I1          I2       V1 V2  H1       H2           11 12
# >[1,]    0 100.00000 100.00000  0.000000 0.00000000 10.000000 0.000000000  0  0  0.0000000 0.000000000  0  0
# >[2,]    1  96.10555  99.91626  7.152329 0.07566424  6.732125 0.008074672  0  0  0.7969166 0.008218785  0  0
# >[3,]    2  93.41096  99.83980 11.009819 0.13217176  5.559244 0.028029699  0  0  2.6469033 0.029051113  0  0
# >[4,]    3  91.13992  99.74612 13.506164 0.19548958  5.323956 0.058385758  0  0  5.1123594 0.061563715  0  0
# >[5,]    4  88.98693  99.61810 15.452465 0.28040111  5.520686 0.101498860  0  0  8.0136521 0.108691194  0  0
# >[6,]    5  86.83334  99.43874 17.191142 0.39911586  5.925642 0.162146190  0  0 11.2799163 0.175973877  0  0


#----------------------------------------------------


rfun.PatchSEIV <- function(init, parms, beta.contact, vax, Time, dt){
  #list frame to hold num individuals in each state in each sim
  
  T= Time/dt 
  
  X <- init%>%
    rfun.FormatInit()
  
  times <- c(seq(dt, Time, by=dt))
  
  #loop through each step of approximation
  for(q in 1:T){ #check to see if you want to write this as second step or based of init values
    
    #current time point
    current_time <- times[q]
    X1 <- X[nrow(X),]
    
    #current N
    N = (rep(0,num_patch)) #create empty frame
    for(j in 1:num_patch){
      for(i in 0:(num_state-1)){
        N[j] = N[j] + X1[j + (i*num_patch)][,1] 
      }
    } 
    
    
    #Pull apart states
    S= as.matrix(t(X1[1, 1:num_patch]))
    E= as.matrix(t(X1[1, (num_patch+1):(2*num_patch)]))
    I= as.matrix(t(X1[1, (2*num_patch+1):(3*num_patch)]))
    V= as.matrix(t(X1[1, (3*num_patch+1):(4*num_patch)]))
    
    #Pull apart parms
    beta = parms$beta
    b = parms$b
    mu =parms$mu
    gamma = parms$gamma
    alpha = parms$alpha
    nu2 <- vax[current_time,] #subset to current vax time
    #nu2 <- inst.vac(nu2)
    
    
    #Terms used
    bSI = beta*S*I
    bSI.contact = beta.contact%*%I*S
    birth.rate = b*N*(1-N/K) #carrying capacity birth rate
    birth.rate = ifelse(birth.rate > 0, birth.rate, 0) #no negative births
    
    #System of equations
    dS <- birth.rate -bSI - bSI.contact - nu2*S - mu*S
    dE <- bSI + bSI.contact -gamma*E - mu*E
    dI <- gamma*E - alpha*I - mu*I
    dV <- nu2*S-mu*V
    
    #abs value of all terms in each state
    S_terms = birth.rate  + bSI + bSI.contact + nu2*S + mu*S
    E_terms = bSI + bSI.contact + gamma*E + mu*E
    I_terms = gamma*E + alpha*I + mu*I
    V_terms = nu2*S + mu*V
    
    #Expected values
    Ex.X <- rbind(dS, dE, dI, dV)
    
    #Covariance matrix
    V11 = matrix(diag(S_terms[,1]), ncol=num_patch)
    V21 = beta.contact*matrix(rep(I, num_patch), ncol=num_patch)*
      matrix(rep(S, num_patch), ncol=num_patch, byrow=TRUE) + 
      matrix(diag(bSI[,1]), ncol=num_patch)
    V31 = matrix(rep(0, num_patch^2), nrow=num_patch)
    V41 = matrix(diag(nu2*S[,1]), ncol=num_patch) #FIX THIS WITH APPROX FUN
    
    V12 = t(V21)
    V22 = matrix(diag(E_terms[,1]), ncol=num_patch)
    V32 = matrix(diag(gamma*E[,1]), ncol=num_patch)
    V42 = matrix(rep(0, num_patch^2), nrow=num_patch)
    
    V13 = t(V31)
    V23 = t(V32)
    V33 = matrix(diag(I_terms[,1]), ncol=num_patch)
    V43 = matrix(rep(0, num_patch^2), nrow=num_patch)
    
    V14 = t(V41)
    V24 = t(V42)
    V34 = t(V43)
    V44 = matrix(diag(V_terms[,1]), ncol=num_patch)
    
    
    V = cbind(rbind(V11, V21, V31, V41),
              rbind(V12, V22, V32, V42),
              rbind(V13, V23, V33, V43),
              rbind(V14, V24, V34, V44))
    
    C=sqrt(V) 
    
    #############################################################################
    #Calculate population at next step
    #############################################################################
    #Create a random vector for each population
    ra = matrix(rnorm(num_comps^2, mean=0, sd=1), nrow = num_comps)
    
    #Sum cov for next step
    F_sum = rowSums((ra*C))
    
    #Caluclate the populations at the next step
    X2 <- X1 + dt*Ex.X + F_sum*sqrt(dt)
    X2 <- ifelse(X2 > 0.5, X2, 0)
    X <- rbind(X, X2)
    
    
  } #end of big loop
  out <- data.frame(time=seq(0, Time, by=dt), X)
  return(out)
}

