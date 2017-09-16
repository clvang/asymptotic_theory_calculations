
rm(list=ls(all=TRUE))   #remove all variables in  workspace

calcRoot <- function(d_c, Udc, d_o, Udo, yo, UYo, y_ofc, err_tol, N, K, UK){
# SUBROUTINE mc_uncertainty( dc_mu, dc_sigma, do_mu, do_sigma, yo_mu, & 
# 						yo_sigma, y_ofc, err_tol, &
# 						dc_mc_out, do_mc_out, yo_mc_out, tau_mc_out, LHS_mc_out, eps_mc, & 
# 						NR, K_mu, dK_sigma, K_mc_out, D_mc )	
	retvals <- .Fortran("mc_uncertainty",dc_mu=as.double(d_c), 
				dc_sigma = as.double(Udc), 
				do_mu = as.double(d_o),
				do_sigma = as.double(Udo), 
				yo_mu = as.double(yo), 
				yo_sigma = as.double(UYo), 
				y_ofc = as.double(y_ofc), 
				err_tol = as.double(err_tol),
				dc_mc_out = double(N),
				do_mc_out = double(N),
				yo_mc_out = double(N),
				tau_mc_out = double(N),
				LHS_mc_out = double(N),
				eps_mc = double(N),
				NR = as.integer(N),
				K_mu = as.double(K),
				dK_sigma = as.double(UK),
				K_mc_out = double(N),
				D_mc = double(N) )
	#there are the returned values
	list(dc_mcVals=retvals$dc_mc_out,
		do_mcVals=retvals$do_mc_out,
		yo_mcVals=retvals$yo_mc_out,
		tau_mcVals=retvals$tau_mc_out,
		LHS_mcVals=retvals$LHS_mc_out,
		eps_mcVals=retvals$eps_mc,
		K_mcVals=retvals$K_mc_out,
		D_mcVals=retvals$D_mc)
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

returnedValues <- calcRoot(d_c =sqrt(dc_sq), Udc=sqrt(Udc_sq), 
		d_o=d_o, Udo=sqrt(Udo_sq), yo=yo, 
		UYo=sqrt(UYo_sq), y_ofc=y_ofc, 
		err_tol=err_tol, N=N, K = K, UK = sqrt(Uk_sq))

# hist(returnedValues$do_mcVals)
# hist(returnedValues$D_mcVals)
D_lower <- quantile(returnedValues$D_mcVals,probs=c((1-P)/2,(1+P)/2))[[1]] * (1/1000)^2
D_upper <- quantile(returnedValues$D_mcVals,probs=c((1-P)/2,(1+P)/2))[[2]] * (1/1000)^2
D_bar <- mean(returnedValues$D_mcVals * (1/1000)^2)
d <- density(returnedValues$D_mcVals * (1/1000)^2)

D_most_probable <- d$x[which(d$y==max(d$y))]

print("do percent_increase :")
print(percent_increase)

print("MC lower D limit [m^2/s]: ")
print(D_lower)
print("MC most probable D [m^2/s]: ")
print(D_most_probable)
print("MC upper D limit [m^2/s]: ")
print(D_upper)

# print(paste0("MC lower D limit: ", D_lower) )
# print(paste0("MC most probable D: ", D_most_probable))
# print(paste0("MC upper D limit: ", D_upper))


hist(returnedValues$D_mcVals*(1/1000)^2,prob=TRUE,n=100,  
	main=paste0("Histogram of D effective"),
	xlab="D [m^2/s]")
d <-density(returnedValues$D_mcVals*(1/1000)^2)
lines(d,col="black",lwd=2)
abline(v=D_upper,col='red',lwd=1.5,lty="dashed")
abline(v=D_lower,col='red',lwd=1.5,lty="dotted")
abline(v=D_most_probable,col='red',lwd=2.3)
# abline(v=D_upper_TSM,col='red',lwd=1.5,lty="dashed")
# abline(v=D_lower_TSM,col='red',lwd=1.5,lty="dotted")
# abline(v=D_TSM,col='red',lwd=2.3)
# legend("topright", c("MC", "TSM"), col=c("green", "red"), lwd=2)
legend("topright", c("MC"), col=c("red"), lwd=2)




















