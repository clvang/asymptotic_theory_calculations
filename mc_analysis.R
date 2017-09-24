
rm(list=ls(all=TRUE))   #remove all variables in  workspace
graphics.off()  #close all graphics windows

calcRoot <- function(tau_vector, LHS_vector, K_vector, err_tol, N){
# SUBROUTINE mc_uncertainty( tau_vector, LHS_vector, err_tol, N, D_vector )
	retvals <- .Fortran("mc_uncertainty",
				tau_vector=as.double(tau_vector),
				LHS_vector=as.double(LHS_vector),
				K_vector = as.double(K_vector),
				err_tol = as.double(err_tol),
				N = as.integer(N),
				D_vector = double(N) )
	#there are the returned values
	list(dc_mcVals=retvals$D_vector)         #values for D
}

if (!is.loaded('mc_uncertainty')){
	dyn.load("mc_uncertainty.so")
}

# read data from fcprops.txt file
data_import <- read.table('fcprops.txt',skip=4,nrows=11, sep="!")
data_numeric <- as.numeric(as.character(data_import$V1[seq(1,10)]))

p 			<- data_numeric[1]    # chamber pressure [atm]
dc_sq 		<- data_numeric[2]    # drop diameter @ onset of fc squared [mm^2]
K 			<- data_numeric[3]    # burning rate constant prior to onset of fc [mm^2/s]
do_measured <- data_numeric[4]    # initial drop diameter measured [mm]
yo 			<- data_numeric[5]    # initial mass frac of low volitiliy component
err_tol 	<- data_numeric[6]    # error tolerrance for bisection method
Uk_sq       <- data_numeric[7]^2  # uncertainty in K squared (U_k ^2) [mm^2/s]^2
Udo_sq      <- data_numeric[8]^2  # uncertainty in do squared (U_do^2) [mm^2]
Udc_sq      <- data_numeric[9]^2  # uncertainty in dc squared (U_dc^2) [mm^2]
UYo_sq      <- (data_numeric[10]*yo)^2  #uncertainty in Yo squared
sol_id      <- data_import$V1[11] #solvent id (1) - Heptane (2)- Propanol


y_ofc <- 1.0  # (defined) mass fraction of low volatility comp at onset of fc

if (sol_id == 1){

	if ( yo == 0.05){
		# heptane95-hexadecane5 d_o corrections
		if (p == 1) {
		    percent_increase <- 0.055  #so far these numbers are only valid for hep-hex exp's
		} else {
		    percent_increase <- 0.086
		}	
	}
	if ( yo == 0.20){
		# heptane80-hexadecane20 d_o corrections
		if (p == 1) {
		    percent_increase <- 0.047  
		} else {
		    percent_increase <- 0.076
		}	
	}

}

if (sol_id == 2){
	# propanol-glycerol d_o corrections
	if (p == 1) {
	    percent_increase <- 0.03  #
	} else {
	    percent_increase <- 0.05  #at 3 atm
	}
}


d_o <- do_measured * (1.0 + percent_increase) #do accounting for droplet swelling
tau_o <- log(d_o / sqrt(dc_sq) )

# NOTE:
# 	tau = LOG(d_o / d_c)
# 	LHS = y_ofc / yo

N <- 1000000
P <- 0.95

print("----- BEGIN random number generation----- ")
#---- generate NORMAL random variables for 95% uncertianties --- ###
factor <- 1.25

set.seed(5)
do_temp <- rnorm( n=N*factor, mean=d_o, sd= sqrt(Udo_sq) )

set.seed(5)
Yo_mcN95 <- rnorm( n=N, mean=yo, sd=sqrt(UYo_sq) )

LHS_mcN95 <- y_ofc / Yo_mcN95
dc_mcN95 <- 1
## sample larger than needed values for dc
## then throw out-of-range values away, element-wise
i <- 5
while(length(dc_mcN95) < N){
	set.seed(i)
	dc_mcN95 <- rnorm(n=N*factor, mean = sqrt(dc_sq), sd=sqrt(Udc_sq) )
	index <- which(do_temp > dc_mcN95)
	dc_mcN95 <- dc_mcN95[index]
	i <- i + 1
}
dc_mcN95seed <- i

dc_mcN95 <- dc_mcN95[1:N]
do_mcN95 <- do_temp[index]
do_mcN95 <- do_mcN95[1:N]
tau_mcN95 <- log(do_mcN95 / dc_mcN95)
## sample larger than needed values for K
## then throw out-of-range values away
set.seed(5)
K_temp <- rnorm(n=N*factor, mean= K, sd= sqrt(Uk_sq))
K_mcN95 <- 1
i <- 5
while( length(K_mcN95)<N ){
	set.seed(i)	
	K_mcN95 <- rnorm(n=N*factor, mean= K, sd= sqrt(Uk_sq))
	K_mcN95 <- K_mcN95[which( K_mcN95 > 0)]
	i <- i + 1	
}
K_mcN95seed <- i
K_mcN95 <- K_mcN95[1:N]
#------------------------------------------------------------ ###



#---- generate NORMAL random variables for STANDARD uncertianties --- ###
## convert 95% uncertainties to standard deviations
cf <- 2.0
sigma_do <- sqrt(Udo_sq) / cf
sigma_dc <- sqrt(Udc_sq) / cf
sigma_K <- sqrt(Uk_sq) / cf
sigma_Y <- sqrt(UYo_sq) / cf

set.seed(5)
do_temp <- rnorm(n=N*factor, mean=d_o, sd= sigma_do )

set.seed(5)
Yo_mcNS <- rnorm( n=N, mean=yo, sd=sigma_Y )

LHS_mcNS <- y_ofc / Yo_mcNS
dc_mcNS <- 1
## sample larger than needed values for dc
## then throw out-of-range values away, element-wise
i <- 5
while(length(dc_mcNS) < N){
	set.seed(i)
	dc_mcNS <- rnorm(n=N*factor, mean = sqrt(dc_sq), sd=sigma_dc )
	index <- which(do_temp > dc_mcNS)
	dc_mcNS <- dc_mcNS[index]
	i <- i + 1
}
dc_mcNSseed <- i

dc_mcNS <- dc_mcNS[1:N]
do_mcNS <- do_temp[index]
do_mcNS <- do_mcNS[1:N]
tau_mcNS <- log(do_mcNS / dc_mcNS)
## sample larger than needed values for K
## then throw out-of-range values away
set.seed(5)
K_temp <- rnorm(n=N*factor, mean= K, sd= sigma_K)
K_mcNS <- 1
i <- 5
while( length(K_mcNS) < N ){
	set.seed(i)
	K_mcNS <- rnorm( n=N*factor, mean= K, sd= sigma_K )
	K_mcNS <- K_mcNS[ which( K_mcNS > 0) ]
	i <- i + 1
}
K_mcNSseed <- i 
K_mcNS <- K_mcNS[1:N]
# dev.new()
# hist(dc_mcNS, prob=TRUE, n=100)
# d <- density(rnorm(n=N, mean = sqrt(dc_sq), sd=sigma_dc ))
# lines(d,col="black",lwd=2)

#-------------------------------------------------------------------- ###


#---- generate UNIFORM random variables for STANDARD uncertianties --- ###
## convert 95% uncertainties to standard deviations
do_max <- d_o + sqrt(Udo_sq)
do_min <- d_o - sqrt(Udo_sq)
sigma_do <- (do_max - do_min) / sqrt(12)

dc_max <- sqrt(dc_sq) + sqrt(Udc_sq)
dc_min <- sqrt(dc_sq) - sqrt(Udc_sq)
sigma_dc <- (dc_max - dc_min) / sqrt(12)

K_max <- K + sqrt(Uk_sq)
K_min <- K - sqrt(Uk_sq)
sigma_K <- 2*sqrt(Uk_sq) / sqrt(12)

Y_max <- yo + sqrt(UYo_sq)
Y_min <- yo - sqrt(UYo_sq)
sigma_Y <- (Y_max - Y_min) / sqrt(12)

set.seed(5)
do_temp <- runif(n=N*factor, min=do_min , max=do_max )

set.seed(5)
Yo_mcUS <- runif( n=N, min=Y_min, max=Y_max )
LHS_mcUS <- y_ofc / Yo_mcUS
## sample larger than needed values for dc
## then throw out-of-range values away, element-wise
dc_mcUS <- 1
i <- 6
while( length(dc_mcUS) < N){
	set.seed(i)
	dc_mcUS <- runif(n=N*factor,min=dc_min,max=dc_max )
	index <- which(do_temp > dc_mcUS)
	dc_mcUS <- dc_mcUS[index]
	i <- i + 1
}
dc_mcUSseed <- i 
dc_mcUS <- dc_mcUS[1:N]
do_mcUS <- do_temp[index]
do_mcUS <- do_mcUS[1:N]
tau_mcUS <- log(do_mcUS / dc_mcUS)
## sample larger than needed values for K
## then throw out-of-range values away
if (K_min < 0){
	K_min <- 0
}
set.seed(5)
K_mcUS <- runif(n=N, min=K_min, max=K_max)

print("----- END random number generation----- ")

#---------------- calculate normal 95% uncertainties ------------------------ ###
D_N95 <- calcRoot(tau_vector=tau_mcN95, 
					LHS_vector=LHS_mcN95, 
					K_vector=K_mcN95,
					err_tol =err_tol, 
					N = N)$dc_mcVals
D_N95 <- D_N95* (1/1000)^2  #convert to m^2/s
DN95_lower <- quantile(D_N95,probs=c((1-P)/2,(1+P)/2))[[1]]
DN95_upper <- quantile(D_N95,probs=c((1-P)/2,(1+P)/2))[[2]] 
DN95_bar <- mean(D_N95 )
d <- density(D_N95)
DN95_most_probable <- d$x[which(d$y==max(d$y))]
print("********************************************************************")
print("--- results assume dependent variables are from NORMAL distribution ---")
print("MC 95% uncertainty LOWER D limit [m^2/s]: ")
print(DN95_lower)
print("MC 95% uncertainty MEAN value D [m^2/s]: ")
print(mean(D_N95))
print("MC 95% uncertainty UPPER D limit [m^2/s]: ")
print(DN95_upper)
print("MC 95% uncertainty MOST PROBABLE D [m^2/s]: ")
print(DN95_most_probable)

hist(D_N95,prob=TRUE,n=100,  
	main=paste0("95% N Histogram of D effective"),
	xlab="D [m^2/s]")
d <-density(D_N95)
lines(d,col="black",lwd=2)
abline(v=DN95_upper,col='red',lwd=1.5,lty="dashed")
abline(v=DN95_lower,col='red',lwd=1.5,lty="dotted")
abline(v=DN95_bar,col='red',lwd=2.3)
legend("topright", c("MC"), col=c("red"), lwd=2)


#---------------- calculate normal standard errors ------------------------ ###
D_NS <- calcRoot(tau_vector=tau_mcNS, 
					LHS_vector=LHS_mcNS, 
					K_vector=K_mcNS,
					err_tol =err_tol, 
					N = N)$dc_mcVals
D_NS <- D_NS * (1/1000)^2  #convert to m^2/s
DNS_bar <- mean(D_NS)
sigma_DNS <- sd(D_NS)
d <- density(D_NS)
DNS_upper <- DNS_bar + sigma_DNS
DNS_lower <- DNS_bar - sigma_DNS
D_NS_most_probable <- d$x[which(d$y==max(d$y))]

print("********************************************************************")
print("--- results assume dependent variables are from NORMAL distribution ---")
print("MC Standard ERROR -1SD D limit [m^2/s]: ")
print(DNS_lower)
print("MC Standard ERROR MEAN value D [m^2/s]: ")
print(DNS_bar )
print("MC Standard ERROR +1SD D limit [m^2/s]: ")
print(DNS_upper)
print("MC Standard ERROR MOST PROBABLE D [m^2/s]: ")
print(D_NS_most_probable)

dev.new()
hist(D_NS,prob=TRUE,n=100,  
	main=paste0("Normal Standard Error Histogram of D effective"),
	xlab="D [m^2/s]")
d <-density(D_NS)
lines(d,col="black",lwd=2)
abline(v=DNS_upper,col='red',lwd=1.5,lty="dashed")
abline(v=DNS_lower,col='red',lwd=1.5,lty="dotted")
abline(v=DNS_bar,col='red',lwd=2.3)
legend("topright", c("MC"), col=c("red"), lwd=2)


#---------------- calculate uniform standard error ------------------------ ###
D_US <- calcRoot(tau_vector=tau_mcUS, 
					LHS_vector=LHS_mcUS, 
					K_vector=K_mcUS,
					err_tol =err_tol, 
					N = N)$dc_mcVals
D_US <- D_US * (1/1000)^2  #convert to m^2/s
DUS_bar <- mean(D_US)
sigma_DUS <- sd(D_US)
d <- density(D_US)
DUS_upper <- DUS_bar + sigma_DUS
DUS_lower <- DUS_bar - sigma_DUS

print("********************************************************************")
print("--- results assume dependent variables are from UNIFORM distribution ---")
print("MC Standard ERROR -1SD D limit [m^2/s]: ")
print(DUS_lower)
print("MC Standard ERROR MEAN value D [m^2/s]: ")
print(DUS_bar )
print("MC Standard ERROR +1SD D limit [m^2/s]: ")
print(DUS_upper)


dev.new()
hist(D_US,prob=TRUE,n=100,  
	main=paste0("Uniform Standard Error Histogram of D effective"),
	xlab="D [m^2/s]")
d <-density(D_US)
lines(d,col="black",lwd=2)
abline(v=DUS_upper,col='red',lwd=1.5,lty="dashed")
abline(v=DUS_lower,col='red',lwd=1.5,lty="dotted")
abline(v=DUS_bar,col='red',lwd=2.3)
legend("topright", c("MC"), col=c("red"), lwd=2)






