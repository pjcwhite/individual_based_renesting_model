############################################################################################################
# R CODE TO ACCOMPANY:                                                                                     #
# Testing re-nesting: Time-variable re-nesting probability functions in avian seasonal productivity models #
# White et al.                                                                                             #
############################################################################################################ 

# NOTE - BY THE NATURE OF THE STOCHASTIC RE-NESTING MODEL, OUTPUTS WILL VARY, EVEN IF RE-NESTING MODEL IS RUN
# USING THE SAME INITIAL CONDITIONS. CODE HERE ACCOMPANIES SECTIONS 'A' AND 'C' OF THE PAPER

# CLEAN ENVIRONMENT IF WANTED
rm(list=ls())

# SET A WORKING DIRECTORY AS YOUR "NETLOGO OUTPUTS" FOLDER (as per SET UP INSTRUCTIONS document)
setwd(dir="C:/Users/40007266/OneDrive - Edinburgh Napier University/Documents/WORK/RESEARCH/PUBLICATIONS/RENESTING PAPER/IBIS SUBMISSION 2022/RENESTING MODEL TO PUBLISH/NETLOGO OUTPUTS")

#########################################################################################
# SECTION A - PATTERN-ORIENTATED MODELLING APPROACH (BLACK REDSTARTS AND YELLOWHAMMERS) #
#########################################################################################

# READ IN CUMULATIVE DISTRIBUTIONS OF FIRST-EGG DATES, SIMULATED AND OBSERVED, FROM "NETLOGO OUTPUTS" FOLDER
combined.simulated<-read.csv("Simulated nest initiations per day - combined.csv",skip=16)
combined.observed<-read.csv("Observed nest initiations - combined.csv",skip=16)

# NOTE IF YOUR WORKING DIRECTLY IS THE "NETLOGO OUTPUTS FOLDER, THEN THIS WILL BE DATA FROM THE LAST RE-NESTING 
# MODEL RUN YOU PERFORMED. EXAMPLE DATA ARE AVAILABLE IN THE GITHUB REPOSITORY AS AN ALTERNATIVE

  #######################################################################
  # CLUTCH INITIATION PERIOD CALCULATIONS (relating to Table 2b and 2g) #
  #######################################################################

  # "clutch initiation period" calculations - observed
  cum.observed<-cumsum(combined.observed$y)/max(cumsum(combined.observed$y))
  sum(cum.observed<=0.1);sum(cum.observed<=0.9) # 10th and 90th centiles
  sum(cum.observed<=0.9)-sum(cum.observed<=0.1) # "clutch inititaion period" sensu Weggler
  # "clutch initiation period" calculation - simulated
  cum.simulated<-cumsum(combined.simulated$y)/max(cumsum(combined.simulated$y))
  sum(cum.simulated<=0.1);sum(cum.simulated<=0.9) # 10th and 90th centiles
  sum(cum.simulated<=0.9)-sum(cum.simulated<=0.1) # "clutch inititaion period" sensu Weggler

  ##############################################################################################
  # KOLMOGOROV-SMIRNOFF TESTS OF NEST INITIATION DISTRIBUTIONS (relating to Table 2a and 2f) #
  ##############################################################################################

  # NOTE - THE PLOTS HERE ARE FOR EXPLANATORY PURPOSES, TO SHOW WHICH PARTS OF THE DISTRIBUTIONS ARE 
  # ACTUALLY COMPARED. THE PLOTS WERE NOT INCLUDED IN THE PAPER.
  
  # Plot of cumulative distribution of nest initiations per Julian day (not shown in paper)
  cum.observed<-cumsum(combined.observed$y)/max(cumsum(combined.observed$y))
  plot(cum.simulated,type="l",lty=2,lwd=2,xlab="Julian day",ylab="Cumulative proportion of initiations",xlim=c(100,300))
  points(cum.observed,type="l",lwd=2)

  # EXCLUDING DATA BEFORE PEAK OF FIRST INITIATIONS, AND AFTER THE POINT WHEN OBSERVED 
  # INITIATIONS FALL TO ZERO
  
  # The peak of nest initiations is indicated on the plot
  peak<-144 # This is the season-start-mean entered into NetLogo 
  abline(v=peak,lty=3) 
  
  # All data before the peak is removed
  cum.simulated<-cum.simulated[-(1:peak)]
  cum.observed<-cum.observed[-(1:peak)]
  
  # Data is removed from where observed initiations fall to zero  
  cum.simulated<-cum.simulated[cum.observed!=1];cum.simulated
  cum.observed<-cum.observed[cum.observed!=1];cum.observed
  
  # The point where observed nest initiations fall to zero is also indicated on the plot
  abline(v=(length(cum.observed)+peak),lty=3)

  # Kolmolgorov-Smirnoff test of the distribution only for included data
  ks.test(x=cum.simulated/max(cum.simulated),y=cum.observed/max(cum.observed))
  max(abs(cum.observed-cum.simulated)) # Maximum difference between the two distributions
  # For interest, scatter-plot of observed and simulated data with line at x=y (not included in paper)
  plot(x=cum.simulated,y=cum.observed,xlim=c(0,1),ylim=c(0,1))
  abline(a=0,b=1)

#########################################################
# SECTION C - POPULATION MODELLING (YELLOWHAMMERS ONLY) #
#########################################################

  ##VARIABLE PARAMETERS FOR DETERMINISTIC MODEL
  DSPE<-0.988 # DAILY SURVIVAL PROBABILTY AT FLEDGLING STAGE DURING PREDATOR REMOVAL
  #DSPE<-0.947 # DAILY SURVIVAL PROBABILTY AT FLEDGLING STAGE WHEN NO PREDATOR REMOVAL
  DSPN<-0.949 # DAILY SURVIVAL PROBABILTY AT FLEDGLING STAGE
  B<-0.51 # POST-FLEDGING SURVIVAL RATE
  n<-2.385 # NUMBER OF ATTEMPTS
  E<-16 # EGG PERIOD
  N<-12 # NESTLING PERIOD
  f<-2.37 # BROOD SIZE AT FLEDGING
  a<-1 # AGE AT FIRST REPRODUCTION
  L<-0.481 # FIRST YEAR SURVIVAL
  p<-0.503 # ADULT ANNUAL SURVIVAL


  # SEASONAL PRODUCTIVITY (NUMBER OF FEMALES PER FEMALE PER SEASON - see paper)
    # DETERMINISTIC FECUNDITY MODEL
      b<- 0.5 * n * f * (DSPE)^E * (DSPN)^N
    # STOCHASTIC FECUNDITY MODEL
      # SEASONAL.PRODUCTIVITY<-1.59  # this is the seasonal productivity value outputted from any run of the individual-based re-nesting model
      # b<- 0.5 * SEASONAL.PRODUCTIVITY

  # POPULATION MODELS AND OUTPUTS FOR Table 3 and Figure 5    
    print("SEASONAL PRODUCTIVITY =")
    round(b*2,3)
    w<-round((log(0.01)/log(p)),1) # age at last reproduction
    print("w =")
    print(w)
    slade<-function(y){
      (p * y^-1) + ((B*L) * b * y^-a) - ((B*L) * b * p^(w-a+1) * y^-(w+1))}
    y<-seq(0.00001,10,by=0.00001)
    iterative<-slade(y)
    iterative<-signif(iterative,5)
    answer<-y[iterative==1.00000e-00]
    answer
    print("POPN GROWTH =")
    round(max(answer),3)
