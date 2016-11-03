# This script belongs to the paper "Development of a pathway model to assess the exposure of European pine trees to pine wood nematode via the trade of wood" 
# written by Bob Douma, Wopke van der Werf, Lia Hemerik, Christer Magnusson and Christelle Robinet
# For details about the model we refer to the paper and its appendices.
# The script below is commented and runs with the files available in github. However, the file may not always be self-explanatory. 
# Please contact Bob Douma (bob.douma@wur.nl) or Christelle Robinet for any questions (christelle.robinet@orleans.inra.fr)
# Version 18
# Wagenigen, Orleans, Oct 31, 2016
##########################################################

# set location where model and files are stored in
persofolder = "D:/BobDouma/Paper_PWN/Model_Coniferous_PWN_China/pathwaymodel_PWN/"

# load libraries (if not installed yet, install them first)
library(reshape)
library(ggplot2)

#################
# FUNCTIONS NEEDED
####################
# give relative weight per country:
frequency.a = function(x,...)
{
  x = x/sum(x,...)
}
# probability to escape from at least one element
probacountry = function(x,...)
{
  1-prod(1-x,...)
}

# this functions calculates the probability that at least one event will take place. t events
# are stored in a list with t elements. The prob.name refers to the element in the list element
# this function is made to calculate the prob of at least one event at the end of the script
prob.t = function(t,x,prob.name){
  prob = 1
  for (i in 1:t){
    eval(parse(text=paste("prob = prob * (1-","x[[",i,"]]$",prob.name,")",sep=""))) 
  }
  prob = 1-prob
  return(prob)
}


# function to determine probability of detection at import
Q = function(Vtot,Vinf,sample.size,consignment.size,insp.per.cons,quarantine=T,type, extract.eff,freq.insp.lots,shortcut=FALSE,q=0,...){ # probability to detect the pest at inspection may depend on the total volume of product imported, the volume infested, and the quarantine status of the pest, and the type of wood (RW or SW)
  # if shortcut = false, the probability to detect is calculated from sampling theory, if true a value can be inserted directly
  if (shortcut == FALSE){
    # probability to detect an infested product (between 0 and 1) #+++++++++++ arbitrary data
    # for roundwood
    if (type=="RW"){
      # total imported volume divided by consignment size, gives potential number of inspections units
      k = floor(consignment.size[1]/sample.size[1]) # potential number of sampling units in consigments
      m = ceiling(k*(Vinf/Vtot)) # infested volume (l) expressed in number of sampling units
      n = k-m # non-infested volume (l) expressed in number of sampling units
      # probability to detect on an inspected lot
      prob.detect = (1-phyper(0,m=m,n=n,k=insp.per.cons))*extract.eff
      # probability to detect a pest given an imported volume of infested wood
      prop.detect = 1-((1-freq.insp.lots) + (freq.insp.lots) * (1-prob.detect))
    }
    # for sawnwood
    if (type=="SW"){
      # total imported volume divided by consignment size, gives potential number of inspections units
      k = floor(consignment.size[2]/sample.size[2]) # potential number of sampling units in consigments
      m = floor(k*(Vinf/Vtot)) # infested volume (l) expressed in number of sampling units
      n = k-m # non-infested volume (l) expressed in number of sampling units
      # probability to detect on an inspected lot
      prob.detect = (1-phyper(0,m=m,n=n,k=insp.per.cons))*extract.eff      # freq.insp.lots = is the frequency of lots that is not inspected!!!
      # prop.detected = the proportion of wood that not passes detection
      prop.detect = 1-((1-freq.insp.lots) + (freq.insp.lots) * (1-prob.detect))
    }
  } else {
    prop.detect = q # the proportion of wood that not passes inspection 
  }
  return(prop.detect)
}


####################
# DATA
####################

# Reference list of the European countries:
#-----------------------------------------
EUcountries.data = read.table(paste(persofolder,"countries_35.txt",sep=""),header=T)
####################


# Country level data:
#--------------------
# amount of RW imported in country j from CO over a given year (in TONS):
# col1: declarant, col2: 2001 data, col2: 2002 data, ..., col11: mean over 2001-2009, col12: sd over 2001-2009
#RW = read.table(paste(persofolder,"oaks.tc.yearly.output_TONS_01082013.csv",sep=""),header=T,sep=",")
RW = read.table(paste(persofolder,"rw.pines.tc.China.yearly.output_TONS_26092013.csv",sep=""),header=T,sep=",")
RW.data = RW[,c(1,11)] # selection of the year of interest, for instance 2009
#RW.data = data.frame(RW[,1],apply(RW[,c(2:10)],1,mean,na.rm=TRUE))
# amounrt of SW imported in country j from CO over a given year (in TONS): 
SW = read.table(paste(persofolder,"sw.pines.tc.China.yearly.output_TONS_26092013.csv",sep=""),header=T,sep=",")
SW.data = SW[,c(1,11)] # selection of the year of interest, for instance 2009
#SW.data = data.frame(SW[,1],apply(SW[,c(2:10)],1,mean,na.rm=TRUE))
#SW[,2] = 1000

# transition matrix for RW (oak) leaving each MSi and going to any MSj, for instance for 2009:
TRANSITIONRW = read.table(paste(persofolder,"transition_matrix_RW_pine_reverse_avg_08102013.csv",sep=""),header=T,sep=",")
# col1: name of departure MSi, col2.... : probability to go from MSi to MSj, the col names give the arrival MSj
TRANSITIONRW=TRANSITIONRW[,-1]
TRANSITIONRW[is.na(TRANSITIONRW)]=0
TRANSITIONRW.data = TRANSITIONRW

# transition matrix for SW (oak) leaving each MSi and going to any MSj, for instance for 2009:
TRANSITIONSW = read.table(paste(persofolder,"transition_matrix_SW_pine_reverse_avg_26092013.csv",sep=""),header=T,sep=",")
TRANSITIONSW=TRANSITIONSW[,-1]
TRANSITIONSW[is.na(TRANSITIONSW)]=0
TRANSITIONSW.data = TRANSITIONSW

# proportion of round wood that is staying in the 35 countries versus outside other countries
relative.flow.EU.RW = read.table(paste(persofolder,"rw.pine.relative.flow.toEU_06032013.csv",sep=""),header=T,sep=",")
# col1: name of MSi, col2....: proportion of wood staying in Europe for different years: 2001...2009
# select appropriate year
relative.flow.EU.RW.mean = apply(relative.flow.EU.RW[,-1],1,mean,na.rm=TRUE)
relative.flow.EU.RW.data = data.frame(relative.flow.EU.RW[,1],relative.flow.EU.RW.mean)

# proportion of sawn wood that is staying in the 35 countries versus outside other countries
relative.flow.EU.SW = read.table(paste(persofolder,"sw.pine.relative.flow.toEU_06032013.csv",sep=""),header=T,sep=",")
# col1: name of MSi, col2....: proportion of wood staying in Europe for different years: 2001...2009
# select appropriate year
relative.flow.EU.SW.mean = apply(relative.flow.EU.SW[,-1],1,mean,na.rm=TRUE)
relative.flow.EU.SW.data = data.frame(relative.flow.EU.SW[,1],relative.flow.EU.SW.mean)

# probability to disperse in each MS, depends on climate, time-window, vector presence in each MS, for now we assume = 1
probadispMS = read.csv(paste(persofolder,"Vector.activity_MS.csv",sep="")) 
probadispMS.data = probadispMS

# Transportation data (distance matrix) (depends on the proportion of factories in each NUTS2 region, as the gravity center in each MS is weighted by this proportion)
DISTANCE = read.table(paste(persofolder,"distancematrix_11032014.txt",sep=""),header=T) # distance between 2 MSs, col1: "country of departure"-"country of arrival", col2: MS1, ... Col: MS27 = distance spent in each MS, in meters 
DISTANCE = as.data.frame(DISTANCE)


# NUTS2-level data
#------------------
nuts2host =read.table(paste(persofolder,"Nuts2_pinus_spp.csv",sep=""),sep=",",header=T) 	# Host cover, col1: NUTS2 ID, col2: area (km?) of the NUTS2 region, col3: area of host (km?)
nuts2host.data = nuts2host
nuts2area = nuts2host[,2]

nuts2factoriestmp = read.table(paste(persofolder,"NUTS2_factories_relative_recalculated_27082014.csv",sep=""),header=T,sep=",") 
#nuts2factoriestmp = data.frame(nuts2factoriestmp[,-3])   	# Factories, col1: NUTS2 ID (2 letters), col2:number of wood processing factories 
nuts2factories = data.frame(nuts2host[,1],nuts2factoriestmp[][match(nuts2host[,1],nuts2factoriestmp[,1]),c(2,3)]) # to match the MS order with host data at NUTS2 level
nuts2factories$country = substr(nuts2factories[,1],start=1,stop=2)	# We add a column to give the name (2 letters) of the MS

# probability to disperse in each NUTS2 region, depends on climate, time-window, vector presence in each MS, for now we assume = 1
# probadispnuts2 = read.table("") # depends on climate, time-window, vector presence in each NUTS2 region, for now we assume = 1
probadispnuts2 = read.csv(paste(persofolder,"Vector.activity_NUTS2.csv",sep="")) 
probadispnuts2 = probadispnuts2[probadispnuts2[,1] %in% nuts2factories[,1],2]
probadispnuts2.data = data.frame(nuts2factories[,1],probadispnuts2)

# Initializing EfactRW and EfactSW
# EfactRW = 0.12, and for UK and IR 0.1
EfactRW = rep(0.12,nrow(nuts2factories))
UK.IR.index = which(nuts2factories[,4]=="UK" | nuts2factories[,4]=="IE")
EfactRW[UK.IR.index]=0.10
EfactRW.data = data.frame(nuts2factories[,1],EfactRW)
# EfactSW = 0.04
EfactSW =  rep(0.04,nrow(nuts2factories))
EfactSW.data = data.frame(nuts2factories[,1],EfactSW)

# Port data
#-----------------
portflow = read.table(paste(persofolder,"ports.import_China_03102013.csv",sep=""),header=T,sep=",") 	
#portflow=portflow[,c(1,2,4)]   # col1: port ID, col2: flow of forestry product from the USA for 2011 (1000 tons), col3:port name
portflow=portflow[,c(5,6)]   # col1: port ID, col2: flow of forestry product from the USA for 2011 (1000 tons), col3:port name
# add a column with country names:
portflow$country = substr(portflow[,1],start=1,stop=2)
# select unique countries for which port information is available
unique.country.ports = unique(portflow$country)
# add a column to data.frame where relative data can be stored in
portflow$relativeflow = NA
for (i in 1:length(unique.country.ports)){
  a = frequency.a(portflow[portflow$country == unique.country.ports[i],2],na.rm=TRUE)
  if (all(is.na(a))){a = rep(1/length(a),length(a))}
  portflow[portflow$country == unique.country.ports[i],4] = a
}
portflow.data = portflow

porthost = read.table(paste(persofolder,"Ports_pinus_spp.csv",sep=""),header=T,sep=",") 
porthost.data = porthost

# modify port import data such that for each country with a seaport import through a port is possible
portflow.data = data.frame(porthost[,1],portflow.data[match(porthost[,1],portflow.data[,1]),4])
portflow.data$country = substr(portflow.data[,1],start=1,stop=2)
sum.ports = tapply(portflow.data[,2],portflow.data$country,sum,na.rm=TRUE)
portflow.data[is.na(portflow.data[,2]),2] = 0
div.ports = tapply(portflow.data[,2],portflow.data$country,length)
add = (1-sum.ports)/div.ports
portflow.data[,2] = portflow.data[,2] + add[match(portflow.data$country,names(add))]
portflow.data$relativeflow = portflow.data[,2]
names(portflow.data) = c("portname","relativeflow","country","relativeflow")
# probability to disperse around each port, depends on climate, time-window, vector presence in each MS, for now we assume = 1
probadispport = read.csv(paste(persofolder,"Vector.activity_ports.csv",sep="")) 
probadispport.data = probadispport
# read population size
popsize.data = read.csv(paste(persofolder,"20150513_population_size_nuts2region.csv",sep="")) 
# match nuts2 regions of population size to nuts2host
popsize.data = data.frame(nuts2host[,1],popsize.data[(match(nuts2host[,1],popsize.data[,1])),4])

####################
# START Model 
####################

Conif_wood_PM = function(
  ####################
  # PARAMETERS
  ####################
  
	EUcountries = EUcountries.data, # vector of EUcountries two letter abbreviation
	RW = RW.data, # imported volume m3 of round wood. Dataframe with length nbMS and two columns: col 1 MS, col 2: imported volume m3 of round wood
	SW = SW.data, # imported volume m3 of sawn wood. Dataframe with length nbMS and two columns: col 1 MS, col 2: imported volume m3 of sawn wood
	TRANSITIONRW = TRANSITIONRW.data, #transition matrix of Roundwood of size nbMS * nbMS
	TRANSITIONSW = TRANSITIONSW.data, #transition matrix of Sawnwood of size nbMS * nbMS
	relative.flow.EU.RW = relative.flow.EU.RW.data, # proportion of round wood exported to extra EU. Col 1 MS, col 2: proportion of round wood exported to countries other than MS
	relative.flow.EU.SW = relative.flow.EU.SW.data, # proportion of sawn wood exported to extra EU. Col 1 MS, col 2: proportion of sawn wood exported to countries other than MS
	probadispMS = probadispMS.data, # vector with dispersal probabilities: vector of size nbMS
	nuts2host = nuts2host.data, # proportion host cover: dataframe: col1: MS, col2 total area of NUTS2region, col3: total area of host
	probadispnuts2  = probadispnuts2.data, # vector with dispersal probabilities: vector of size NUTS2
	EfactRW = EfactRW.data, # escape probabilities of pest at F1 nuts2level: vector of length NUTS2
	EfactSW = EfactSW.data, # escape probabilities of pest at F2 nuts2level: vector of length NUTS2
	portflow = portflow.data, # relative importance of ports. Dataframe with col1: port code, col2: absolute import, col3: country code, col4: relative amount
	porthost = porthost.data,# host cover at ports: dataframe: col1: port code, col2: host cover
	probadispport = probadispport.data, # probability to disperse at ports, vector with length NUSTsregions
  popsize = popsize.data, # number of inhabitants per nuts2 region
#   
	wooddensity = 428,  # kg/m3, to correct the importation data in tons
	ICO = 0.024,		# infestation rate in CO (between 0 and 1)
	pb = 0.05,		# proba to find pest in the bark (between 0 and 1), 1-pb gives the probability to find the pest in other parts of the tree
	rhob = 0.05,		# proportion of imported products with bark (between 0 and 1)
	sample.size = c(0.0003,0.0003), # sample size in m3 for RW and SW!!!
	consignment.size = c(55,55),# m3 (based on finish data for RW and SW)
	insp.per.cons = c(2,2), # number of samples per consignment for RW and SW
	extract.eff = 0.75,
	freq.insp.lots = 1, # percentage of inspected lots
#  
	d = 2,			# dispersal distance (km) of the pest or vector
	dres = 0.5, # dispersal distance (km) for non-vector pest transmission (bound to the pixel size of 1km)
	transportationspeedRW = 50, # speed of trucks in km/h
	transportationspeedSW = 70, # speed of trucks in km/h
	workinghoursperday = 7,	# number of working hours per day (for driving trucks for the transportation of wood)
	
# number of propagules per m3 product
	Nind.RW = 5000, # number of propagules per m3 roundwood 	
	Nind.SW = 5000, # number of propagules per m3 sawnwood
	Nind.res = 5000, # number of propagules per m3 of residues
	Nind.FP = 5000, # number of propagules per m3 of final product
#	
  ERW = 0.10,		# probability that one propagule escapes from 1 m3 of infested RW around ports (time implicit)
	ESW = 0.04,		# probability that one propagule escapes from 1 m3 of infested SW around ports (time implicit)
	EtransRW = 0.0075,	# probability that one propagule escapes from 1 m3 of infested RW when transported (per day)
	EtransSW = 0.002,	# probability that one propagule escapes from 1 m3 of infested SW when transported (per day)
	Eres= 0.00001,		# probability that one propagule escapes from 1 m3 of infested residues (time implicit)
	EFP	= 0.000001,		# probability that one propagule escapes from 1 m3 of infested FP (time implicit)
	## all these values should be 0 and 1
#	  
	ar = 0.5,		# proportion of wood ending in residues at 1st transformation and that could still transmit the pest
	arsw = 0.5,		# proportion of wood ending in residues at 2nd transformation and that could still transmit the pest
	### all these values should be 0 and 1
#  
	SCO = 0.05,		# survival to treatment in CO
	SF1CO = 1,    # survival to sawnwood in CO
	ST = 0.5,		# survival rate to treatment in EU
	SF1 = 1.0,		# survival rate to first transformation
	SF2 = 1.0,		# survival rate to second transformation
	Sr = 1.0,		# survival rate in wood residues
  b = 0.5,			# proportion of SW going to 2nd transformation (among SW not exported) => 1-b = proportion locally consumed
  ### all these values should be 0 and 1
  
	nsteps = 2,		# Number of exchanges between MSs
	sensitivity = FALSE, # default no sensitivity analysis
	sensitivity.value = 1.1, # percentage increase in parameter value
	sens.parameter = NA, # name of parameter
	par.type = NA, # does the parameter represents a probabiliy (P) or another value (NP)
	...)
{
  ###########################################################################
  # Do some additional data manipulation based on other data
  ###########################################################################
  
  # number of EU countries
  nbMS = nrow(EUcountries) 
  nbNUTS2 = 310
  nbPorts = 343

  relative.flow.EU.RW =  relative.flow.EU.RW[match(EUcountries[,1],relative.flow.EU.RW[,1]),] # makes that the data.frame matches the list of countries
  relative.flow.EU.RW[is.na(relative.flow.EU.RW)] = mean(relative.flow.EU.RW[,2],na.rm=TRUE) # replace values that are NA by mean of that year
  
  relative.flow.EU.SW =  relative.flow.EU.SW[match(EUcountries[,1],relative.flow.EU.SW[,1]),] # makes that the data.frame matches the list of countries
  relative.flow.EU.SW[is.na(relative.flow.EU.SW)] = mean(relative.flow.EU.SW[,2],na.rm=TRUE) # replace values that are NA by mean of that year
  
  # col1: name of ports, col2 - end: probability to encounter a host given host cover at various dispersal distance (1 km, 2 km, 5 km, 10 km, 25 km, 50 km)
  if(d==1){porthost = porthost[,c(1,2)]}
  if(d==2){porthost = porthost[,c(1,3)]}
  if(d==5){porthost = porthost[,c(1,4)]}
  if(d==10){porthost = porthost[,c(1,5)]}
  if(d==25){porthost = porthost[,c(1,6)]}
  if(d==50){porthost = porthost[,c(1,7)]}
  
  # select ports importing pine only:
  probhost = porthost[match(portflow[,1], porthost[,1]),]
  # select dispersal probability for ports that import pine
  probadispport = probadispport[match(portflow[,1], probadispport[,1]),c(1,4)]

  # match to EU countries
  probadispMS = probadispMS[as.character(probadispMS[,1]) %in% EUcountries[,1],]
    
  # calculate host cover at country level by nuts region in this country
  nuts2host$country = substr(nuts2host[,1],start=1,stop=2)
  nuts2host.orig = nuts2host
  countryhosttmp = tapply(nuts2host[,3],nuts2host$country, sum,na.rm=T)/tapply(nuts2host[,2],nuts2host$country, sum,na.rm=T)
  countryhost = data.frame(unique(nuts2host[,4]),countryhosttmp)
  countryhost = data.frame(RW[,1],countryhost[,2][match(RW[,1],countryhost[,1])])
  # col1: country names, col2: mean host proportion per country
  # This "nuts2host" file defines the reference list for the NUTS2regions considered:

  nuts2host = data.frame(nuts2host[,1],nuts2host[,3]/nuts2host[,2])


  
  ########################################################################################
  # Check whether data is of correct size
  ########################################################################################
  if(nrow(RW) != nbMS | ncol(RW) != 2) stop("dataframe dimensions RW not correct")
  if(nrow(SW) != nbMS | ncol(SW) != 2) stop("dataframe dimensions SW not correct")
  if(dim(TRANSITIONRW)[1] != nbMS | dim(TRANSITIONRW)[2] != nbMS) stop("dataframe dimensions TRANSITIONRW not correct")  
  if(dim(TRANSITIONSW)[1] != nbMS | dim(TRANSITIONSW)[2] != nbMS) stop("dataframe dimensions TRANSITIONSW not correct")  
  if(dim(relative.flow.EU.RW)[1] != nbMS | dim(relative.flow.EU.RW)[2] != 2) stop("dataframe dimensions relative.flow.EU.RW not correct")  
  if(dim(relative.flow.EU.SW)[1] != nbMS | dim(relative.flow.EU.SW)[2] != 2) stop("dataframe dimensions relative.flow.EU.RW not correct")  
  #if(length(probadispMS) != nbMS ) stop("vector length probabdispMS not correct")
  if(nrow(probadispMS) != nbMS  | ncol(probadispMS) != 2) stop("vector length probabdispMS not correct")
  if(dim(nuts2host)[1] != nbNUTS2 | dim(nuts2host)[2] != 2) stop("dataframe dimensions nuts2host not correct")  
  if(nrow(probadispnuts2) != nbNUTS2 | ncol(probadispnuts2) != 2) stop("vector probadispnuts2 not correct")
  if(nrow(EfactRW) != nbNUTS2 | ncol(EfactRW) != 2) stop("vector length EfactRW not correct")
  if(nrow(EfactSW) != nbNUTS2 | ncol(EfactSW) != 2) stop("vector length EfactSW not correct")
  if(dim(portflow)[1] < nbMS | dim(portflow)[2] != 4) stop("dataframe dimensions portflow not correct")  
  if(dim(porthost)[1] != nbPorts | dim(porthost)[2] != 2) stop("dataframe dimensions portflow not correct")  
  if(nrow(probadispport) != dim(portflow)[1] | ncol(probadispport) != 2) stop("vector length probadispport not correct")
  

  #############################################################################
  # ESTIMATING ALPHA (proportion of RW going to 1st transformation in each MS) 
  #        and BETA (proportion of SW going to 2nd transformation in each MS)
  #############################################################################
  
  # INPUT FILE
  inputfile = read.table(paste(persofolder,"wood_flow_for_alpha_beta_conif.csv",sep=""),header=T,na.strings="NA",sep=",")
  #attach(inputfile)
  
  # read species cover in EU
  species.cover = read.csv(paste(persofolder,"country_treecover_Pine_NonPine.csv",sep=""),header=T,sep=",")
  
  # adapt below to the "inputfile" that is read in
  SPEC_RW_IMP = inputfile$PIN_RW_IMP
  SPEC_RW_EXP = inputfile$PIN_RW_EXP
  SPEC_SW_IMP = inputfile$PIN_SW_EXP
  SPEC_SW_EXP = inputfile$PIN_SW_EXP
  
  # transforming NA to zero
  species.prod = inputfile$RW_PROD*(species.cover[,4])
  species.prod[is.na(species.prod)]<-0
  SPEC_RW_IMP[is.na(SPEC_RW_IMP)]<-0
  SPEC_RW_EXP[is.na(SPEC_RW_EXP)]<-0
  
  # calculating alpha:
  # proportion of roundwood that is transformed to sawnwood 
  # vector (value for each MS)
  TRANSFORM_RW_SW = SPEC_RW_IMP +  species.prod - SPEC_RW_EXP
  TRANSFORM_RW_SW[TRANSFORM_RW_SW<0] = 0
  alpha = TRANSFORM_RW_SW / (SPEC_RW_EXP+TRANSFORM_RW_SW)
  
  # calculating beta:
  # proportion of sawnwood that is leaving the sawnwood chain (transformed or consumed)  
  # vector (value for each MS)
  LOCALPROD_SW  = TRANSFORM_RW_SW *(1-ar)				
  species.prod.sw = LOCALPROD_SW *(species.cover[,4])
  species.prod.sw[is.na(species.prod.sw)]<-0
  SPEC_SW_IMP[is.na(SPEC_SW_IMP)]<-0
  SPEC_SW_EXP[is.na(SPEC_SW_EXP)]<-0   
  LEAVING_SW_CHAIN = SPEC_SW_IMP + species.prod.sw - SPEC_SW_EXP 	
  LEAVING_SW_CHAIN[LEAVING_SW_CHAIN<0] = 0
  beta = LEAVING_SW_CHAIN   / (SPEC_SW_EXP+LEAVING_SW_CHAIN )
  
  woodflow = inputfile[,-(5:14)]
  woodflow = cbind(woodflow,alpha,beta)
  
  woodflow = woodflow[match(RW[,1],woodflow[,1]),]
  
  alpha = data.frame(woodflow[,1],woodflow$alpha) # proportion of roundwood that is transformed to sawnwood, estimated with total wood trade 
  beta = 	data.frame(woodflow[,1],woodflow$beta) # proportion of sawnwood that is leaving the sawnwood chain (transformed or consumed)	

  ############### Sensitivity Analysis - Parameter change ################
  # Do sensitivity analysis if sensitivity == TRUE
  if (sensitivity){
    parameter = as.character(sens.parameter)
    # if the parameter does not represent a probability
    if (par.type == "NP"){
      # if the parameter is a single value
      if (is.null(ncol(get(parameter)))){
        old.param = get(parameter)
        # assigns a value to the character name stored in "parameter" and multiplies this with the sensitivity value
        assign(parameter,get(parameter)*sensitivity.value)
        rel.change = mean(get(parameter)/old.param,na.rm=TRUE)
        old.param = old.param
        new.param = get(parameter)
        # if the parameter is a data.frame, the second column will be changed in the sens.analysis  
      } else {
        old.param = get(parameter)[,2]
        assign(parameter,data.frame(get(parameter)[,1],get(parameter)[,2]*sensitivity.value))
        rel.change = mean(get(parameter)[,2]/old.param,na.rm=TRUE)
        old.param = mean(old.param,na.rm=TRUE)
        new.param = mean(get(parameter)[,2],na.rm=TRUE)
      }
      # if parameter represents a probability    
    } else {
      # if the parameter is a single value
      if (is.null(ncol(get(parameter)))){
        old.param = get(parameter)
        assign(parameter,min(1,get(parameter)*sensitivity.value))
        rel.change = mean(get(parameter)/old.param,na.rm=TRUE)
        old.param = old.param
        new.param = get(parameter)
      } else {
        old.param = get(parameter)[,2]
        assign(parameter,data.frame(get(parameter)[,1],pmin(1,get(parameter)[,2]*sensitivity.value)))
        rel.change = mean(get(parameter)[,2]/old.param,na.rm=TRUE)
        old.param = mean(old.param,na.rm=TRUE)
        new.param = mean(get(parameter)[,2],na.rm=TRUE)
      }
    }
  }
  ######################################################################
  
  alphavect=alpha[,2]; alphavect[is.na(alpha[,2])]=0
  alphavect_1=1-alpha[,2]; alphavect_1[is.na(alpha[,2])]=0
  betavect=beta[,2]; betavect[is.na(beta[,2])]=0
  betavect_1=1-beta[,2]; betavect_1[is.na(beta[,2])]=0
    
  ########################
  # Output variable:
  ########################
  nbcolres = 4+3+6*nsteps
  # in RESULTS probability of transfer at country level is stored in
  RESULT = matrix(rep(NA,nbcolres*nbMS),ncol=nbcolres) 
  RESULT = as.data.frame(RESULT)
  # in RESULTS exposure at country level is stored in
  RESULT_E = matrix(rep(NA,nbcolres*nbMS),ncol=nbcolres) 
  RESULT_E = as.data.frame(RESULT_E)
  
  #  col 1:  MS
  #  col 2:  Inf vol arriving
  #  col 3:  Inf vol not detected
  #  col 4:  proba disseminate ports
  #  col 5:  proba disseminate during MS exchange 1 of RW
  #  col 6:  proba disseminate during MS exchange 1 of SW
  #  col 7:  proba disseminate from RW at F1 (before tranformation)
  #  col 8:  proba disseminate from SW (including treatment of RW) when stored around factories (F2) before 2nd transformation
  #  col 9:  proba disseminate from residues around factories F1 and F2 
  #  col 10:  proba disseminate from FP (including survival to 2nd transformation)
  
  RESULT[,1]= EUcountries
  RESULT_E[,1]= EUcountries
  
  
  #####################
  # Step 0 
  # in CO
  #####################
  IexportRW = ICO*SCO
  IexportSW = ICO*SCO*SF1CO
  
  #####################
  # Step 1 
  # Arrival
  #####################
  # Convert mass via wooddensity to volume 
  RW[,2]=RW[,2]/(wooddensity*10^(-3)) # converts metric tons to m3
  RW = as.data.frame(RW)    
  RW[is.na(RW[,2]),2]=0 # set to 0 for which no data is available
  
  SW[,2]=SW[,2]/(wooddensity*10^(-3)) # converts metric tons to m3
  SW = as.data.frame(SW)    
  SW[is.na(SW[,2]),2]=0 # set to 0 for which no data is available
  SW = SW[match(RW[,1],SW[,1]),] # makes that the SW data.frame matches the size of the RW data.frame
  
  RW1 = RW[,2]*IexportRW*(rhob + (1-rhob)*(1-pb)) # arrival of roundwood
  SW1 = SW[,2]*IexportSW # arrival of sawn wood
  RESULT[,2]=RW1
  RESULT[,3]=SW1
  RESULT_E[,2]= RW1
  RESULT_E[,3]= SW1
  
  #####################
  # Step 2
  # Detection
  #####################
  qRW = Q(Vtot=RW[,2],Vinf=RW1,quarantine=T,type="RW",consignment.size=consignment.size,insp.per.cons=insp.per.cons, sample.size=sample.size,extract.eff = extract.eff,freq.insp.lots=freq.insp.lots) # infested round wood escaping detection
  qRW[is.na(qRW)] = 0
  V2RW = RW1*(1-qRW)
  RESULT[,4]=V2RW
  RESULT_E[,4]=V2RW
  V2RW = data.frame(EUcountries,V2RW)
  # inspection
  q = Q(Vtot=SW[,2],Vinf=SW1,quarantine=T,type="SW",consignment.size=consignment.size,insp.per.cons=insp.per.cons, sample.size=sample.size,extract.eff = extract.eff,freq.insp.lots=freq.insp.lots)
  # if no volume
  q[is.na(q)] = 0
  V2SW = SW1*(1-q)  # infested sawn wood escaping detection
  RESULT[,5]=V2SW
  RESULT_E[,5]=V2SW
  V2SW = data.frame(EUcountries,V2SW)
  
  
  #####################
  # Step 3
  # Transfer around ports
  #####################
  
  # probability to escape from the ports
  #--------------------------------------
  V2RWkj = V2RW[,2][match(portflow$country,V2RW[,1])]*portflow$relativeflow # infested volume roundwood entering in each port after inspection
  #P3kj = (1-(1-E)^(V2kj))*probadispport*probhost
  
  V2SWkj = V2SW[,2][match(portflow$country,V2SW[,1])]*portflow$relativeflow # infested volume sawnwood entering in each port after inspection
  
  Phost = probhost[,2] # host cover around the port as defined by the dispersal capacity of the vector
  
  PescRWkj = 1-((1-ERW)^(V2RWkj*Nind.RW)) # Probability of at least one escape in Rounw Wood
  PescSWkj = 1-((1-ESW)^(V2SWkj*Nind.SW)) # Probability of at least one escape in Sawn Wood
  VescRWkj = V2RWkj*Nind.RW*ERW
  VescSWkj = V2SWkj*Nind.SW*ESW
  P3RWkj = (1-(1-Phost)^(V2RWkj*ERW*probadispport[,2]*Nind.RW))
  P3SWkj = (1-(1-Phost)^(V2SWkj*ESW*probadispport[,2]*Nind.SW))
  V3RWkj = V2RWkj*Nind.RW*ERW*probadispport[,2]*Phost
  V3SWkj = V2SWkj*Nind.SW*ESW*probadispport[,2]*Phost
  
  # probability to escape from at least one port in each country and infested volume that encounters a host
  #-------------------------------------------------------------
  # P for round wood (RW)
  P3RWj = tapply(P3RWkj,portflow$country,probacountry, na.rm=TRUE) #,na.rm=T
  PRWport= P3RWj[match(RW[,1], names(P3RWj))]
  names(PRWport) = RW[,1]
  PRWport[is.na(PRWport)]=0
  # V for round wood (RW)
  V3RWj = tapply(V3RWkj,portflow$country,sum,na.rm=TRUE) #,na.rm=T # sum of propagules
  VRWport= V3RWj[match(RW[,1], names(V3RWj))]
  names(VRWport) = RW[,1]
  VRWport[is.na(VRWport)]=0
  
  RESULT[,6]=PRWport
  RESULT_E[,6]=VRWport
  
  # P for sawn wood (SW)
  P3SWj = tapply(P3SWkj,portflow$country,probacountry,na.rm=TRUE) #,na.rm=T
  PSWport= P3SWj[match(SW[,1], names(P3SWj))]
  names(PSWport) = SW[,1]
  PSWport[is.na(PSWport)]=0
  # V for sawn wood (SW)
  V3SWj = tapply(V3SWkj,portflow$country,sum,na.rm=TRUE) #,na.rm=T  # sum of propagules
  VSWport= V3SWj[match(SW[,1], names(V3SWj))]
  names(VSWport) = SW[,1]
  VSWport[is.na(VSWport)]=0
  
  RESULT[,7]=PSWport
  RESULT_E[,7]=VSWport
  #####################
  # Step 4
  # Wood flow
  #####################
  
  listwoodflow = list()
  RWenterj = V2RW[,2]
  SWenterj = V2SW[,2]
  for (t in 1:nsteps)
  {
  
  # Step 4-1: intra EU trade
  # Wood which is not transformed in the country where it entered
  #----------------------------------------------------------------
  RWexitj = RWenterj*(alphavect_1) 			# volume of RW leaving country j (not transformed)
  SWexitj = RWenterj*alphavect*(1-ar)*(betavect_1) + SWenterj*(betavect_1)	# volume of SW leaving country j (not transformed)
  
  # Dissemination during transportation of RW
  #------------------------------------------
  TRANSITIONRWVolume = TRANSITIONRW * RWexitj *relative.flow.EU.RW[,2] # matrix giving the volume going each MS to any other MSs, corrected for the amount of wood that is exported outside EU
   # melt this matrix into 3 col data frame (col1: departure MS, col2: arrival MS, col3: volume)
  VECTORRWVolume = suppressMessages(melt(data.frame(V2RW[,1],TRANSITIONRWVolume)))
  nameVECTORRWVolume  = paste(VECTORRWVolume[,1],"-",VECTORRWVolume[,2],sep="")
  DURATION = DISTANCE[match(nameVECTORRWVolume,DISTANCE[,1]),2:ncol(DISTANCE)]/(1000*transportationspeedRW*workinghoursperday) # in days
  Phost.transport = countryhost[,2] # proportion host cover at country scale
  # matrices needed for calculation of PtransRW # 1225 rows, 35 colums, each row is a country-country combination
  Phost.mat = matrix(Phost.transport,nrow=nrow(DURATION),ncol=ncol(DURATION),byrow=TRUE) # probability of host for each country in the format of the distance matrix
  PdispMS.mat = matrix(probadispMS[,2],nrow=nrow(DURATION),ncol=ncol(DURATION),byrow=TRUE) # probability of dispersal for each country-country combination
  VECTORRWVolume.mat = matrix(VECTORRWVolume[,3],nrow=nrow(DURATION),ncol=ncol(DURATION),byrow=FALSE)
  # probabilities of escape and transfer
  Eduration = 1-(1-EtransRW)^DURATION  # probability that at least one escape will take place during transport for one unit of infested wood 
  Pescapetmp =1-(1-EtransRW)^(DURATION*VECTORRWVolume.mat*Nind.RW)
  Pescape.trans.RW = apply(Pescapetmp,2,probacountry,na.rm=TRUE)
  PtransRWtmp = 1-(1-Phost.mat)^(EtransRW*VECTORRWVolume.mat[,3]*Nind.RW*DURATION*PdispMS.mat) 
  PtransRW =  apply(PtransRWtmp,2,probacountry,na.rm=TRUE)
  # Number of escapes
  Vescapetmp =EtransRW*DURATION*VECTORRWVolume.mat*Nind.RW
  Vescape.transRW = apply(Vescapetmp,2,sum,na.rm=TRUE) 
  VtransRWtmp = Vescapetmp * Phost.mat * PdispMS.mat
  VtransRW = apply(VtransRWtmp,2,sum,na.rm=TRUE) # number of propagules transfer from RW to host during transport
  
  RESULT[,7+(t-1)*6+1]=PtransRW 
  RESULT_E[,7+(t-1)*6+1]=VtransRW 
  # Dissemination during transportation of SW
  #------------------------------------------
  TRANSITIONSWVolume = TRANSITIONSW * SWexitj * relative.flow.EU.SW[,2] # matrix giving the volume going each MS to any other MSs, corrected for amount of SW that is exported to outside EU
  # melt this matrix into 3 col data frame (col1: departure MS, col2: arrival MS, col3: volume)
  VECTORSWVolume = suppressMessages(melt(data.frame(V2SW[,1],TRANSITIONSWVolume)))
  nameVECTORSWVolume  = paste(VECTORSWVolume[,1],"-",VECTORSWVolume[,2],sep="")
  DURATION = DISTANCE[match(nameVECTORSWVolume,DISTANCE[,1]),2:ncol(DISTANCE)]/(1000*transportationspeedSW*workinghoursperday) # in days
  Phost.transport = countryhost[,2]
  # matrices needed for calculation of PtransRW # 1225 rows, 35 colums, each row is a country-country combination
  Phost.mat = matrix(Phost.transport,nrow=nrow(DURATION),ncol=ncol(DURATION),byrow=TRUE) # probability of host for each country in the format of the distance matrix
  PdispMS.mat = matrix(probadispMS[,2],nrow=nrow(DURATION),ncol=ncol(DURATION),byrow=TRUE) # probability of dispersal for each country-country combination
  VECTORSWVolume.mat = matrix(VECTORSWVolume[,3],nrow=nrow(DURATION),ncol=ncol(DURATION),byrow=FALSE)  
  # probabilities of escape and transfer
  Eduration = 1-(1-EtransSW)^DURATION  # probability that at least one escape will take place during transport for one unit of infested wood 
  Pescapetmp =1-(1-EtransSW)^(DURATION*VECTORSWVolume.mat*Nind.SW)
  Pescape.trans.SW = apply(Pescapetmp,2,probacountry,na.rm=TRUE)
  PtransSWtmp = 1-(1-Phost.mat)^(EtransSW*VECTORSWVolume.mat[,3]*Nind.SW*DURATION*PdispMS.mat) 
  PtransSW =  apply(PtransSWtmp,2,probacountry,na.rm=TRUE)
  #  Number of escapes
  Vescapetmp =EtransSW*DURATION*VECTORSWVolume.mat*Nind.SW
  Vescape.transSW = apply(Vescapetmp,2,sum,na.rm=TRUE) 
  VtransSWtmp = Vescapetmp * Phost.mat * PdispMS.mat
  VtransSW = apply(VtransSWtmp,2,sum,na.rm=TRUE) # number of propagules transfer from RW to host during transport

  RESULT[,7+(t-1)*6+2]=PtransSW
  RESULT_E[,7+(t-1)*6+2]=VtransSW
  # Step 4-2: transformation chain and dissemination along the chain
  #----------------------------------------------------------------
  # proportion going to 1st transformation
  RWF1 = RWenterj*alphavect # volume of RW undergoing 1st transformation
  
  # Dissemination from RW when stored around factories (F1) before  transformation (NUTS2 region level)
  #-----------------------------------------------------------------------------------------------------
  nuts2factories$volimportedj = RWF1[match(nuts2factories$country,RW[,1])]
  nuts2factories$volimportedNUTS2 = nuts2factories$volimportedj * nuts2factories[,3]
 
  area.factory = nuts2factories[,2]*pi*(d^2) # area that is covered by vectors released at factories in nuts2region
  Phost = (nuts2host[,2]) * pmin(area.factory/nuts2area,1) # host cover in nut2region * proportion of area covered by vector
  PfactRWnuts2 = 1-(1-Phost)^(nuts2factories$volimportedNUTS2*Nind.RW*EfactRW[,2]*probadispnuts2[,2]) # probability to disseminate from RW per NUTS2 region
  PfactRW=tapply(PfactRWnuts2,nuts2factories$country,probacountry,na.rm=TRUE) # probability to disseminate from RW per country
  Pescape.fact.RW = 1-(1-EfactRW[,2])^(nuts2factories$volimportedNUTS2*Nind.RW)
  Vescape.fact.RW = EfactRW[,2]*nuts2factories$volimportedNUTS2*Nind.RW
  VfactRWnuts2 = Vescape.fact.RW*probadispnuts2[,2]*Phost # number of propagules that disseminates from RW per NUTS2 region
  VfactRW = tapply(VfactRWnuts2,nuts2factories$country,sum,na.rm=TRUE)
  
  RESULT[,7+(t-1)*6+3]=PfactRW
  RESULT_E[,7+(t-1)*6+3]=VfactRW
  # Dissemination from SW (including treatment of RW) when stored around factories (F2) before transformation (NUTS2 region level)
  #------------------------------------------------------------------------------------------------------------------------------
  SWF2 = SWenterj + (RWF1 *(1-ar)) * ST * SF1 # volume of SW imported + produced after 1st transformation
  
  nuts2factories$volimportedj = (b*(betavect)*SWF2)[match(nuts2factories$country,RW[,1])]
  nuts2factories$volimportedNUTS2 = nuts2factories$volimportedj * nuts2factories[,3]
  
  area.factory = nuts2factories[,2]*pi*(d^2) # area that is covered by vectors released at factories in nuts2region
  Phost = (nuts2host[,2]) * pmin(area.factory/nuts2area,1) # host cover in nut2region * proportion of area covered by vector
  PfactSWnuts2 = 1-(1-Phost)^(nuts2factories$volimportedNUTS2*Nind.SW*EfactSW[,2]*probadispnuts2[,2])
  PfactSW=tapply(PfactSWnuts2,nuts2factories$country,probacountry,na.rm=TRUE) # probability to disseminate from SW per country
  Pescape.fact.SW = 1-(1-EfactSW[,2])^(nuts2factories$volimportedNUTS2*Nind.SW)
  # number of propagules transfering to host
  Vescape.fact.SW = EfactSW[,2]*nuts2factories$volimportedNUTS2*Nind.SW
  VfactSWnuts2 = Vescape.fact.SW*probadispnuts2[,2]*Phost  # volume that disseminates from SW per NUTS2 region
  VfactSW = tapply(VfactSWnuts2,nuts2factories$country,sum,na.rm=TRUE)  
  
  RESULT[,7+(t-1)*6+4]=PfactSW
  RESULT_E[,7+(t-1)*6+4]=VfactSW
  # sawn Wood which is not transformed in the country where it entered
  #--------------------------------------------------------------------
  SWexitj = SWF2 *(betavect_1) # volume of SW leaving country j (not transformed)
  
  SWenterj = t(TRANSITIONSW) %*% SWexitj
  
  # Dissemination from residues around factories F1 and F2 (NUTS2 region level)
  #----------------------------------------------------------------------------
  Residues = RWF1 *ar * Sr + SWF2*b*betavect*arsw*Sr
  Residues = Residues[match(nuts2factories$country, RW[,1])]
  ResiduesNUTS2 = Residues * nuts2factories[,3]
  
  area.factory = nuts2factories[,2]*pi*(dres^2) # area that is covered by vectors released at factories in nuts2region
  Phost = (nuts2host[,2]) * pmin(area.factory/nuts2area,1) # host cover in nut2region * proportion of area covered by vector
  Pres.nuts2 = 1-(1-Phost)^(ResiduesNUTS2*Eres*probadispnuts2[,2]*Nind.res)  # probability to disseminate from residues per NUTS2 region
  Pres=tapply(Pres.nuts2,nuts2factories$country,probacountry,na.rm=TRUE) # probability to disseminate from residues per country
  Pescape.res.nuts2 = 1-(1-Eres)^(ResiduesNUTS2*Nind.res)
  
  Vescape.res.nuts2 = ResiduesNUTS2*Eres*Nind.res 
  Vres.nuts2 = Vescape.res.nuts2*probadispnuts2[,2]*Phost # volume thatdisseminates from residues per NUTS2 region
  Vres=tapply(Vres.nuts2,nuts2factories$country,sum,na.rm=TRUE) # volume thatdisseminates from residues per country

  RESULT[,7+(t-1)*6+5]= Pres
  RESULT_E[,7+(t-1)*6+5]= Vres
  # Dissemination from FP (including treatment of SW) all over the country (country level)
  #-------------------------------------------------------------------------------
  FP = SWF2*betavect*(b*(1-arsw)*ST*SF2 + (1-b)) 
  FP = FP[match(nuts2factories$country, RW[,1])]
  FP.nuts2 = FP * popsize.data[,2]
  Phost = nuts2host[,2]
  PFP.nuts2 = 1-(1-Phost)^(FP.nuts2*EFP*as.data.frame(probadispnuts2[,2])*Nind.FP)
  PFP = tapply(PFP.nuts2,nuts2factories$country,sum,na.rm=TRUE)
  Pescape.FP.nuts2 = 1-(1-EFP)^(FP.nuts2*Nind.FP)
  Vescape.FP.nuts2 = EFP*FP.nuts2*Nind.FP
  VFP.nuts2 = Vescape.FP.nuts2*probadispnuts2[,2]*Phost
  VFP = tapply(VFP.nuts2,nuts2factories$country,sum,na.rm=TRUE) # volume thatdisseminates from final product per country
  

  RESULT[,7+(t-1)*6+6]= PFP 
  RESULT_E[,7+(t-1)*6+6]= VFP
  # output that any transfer takes place...
  
  listoutput =      list(P3RWkj,P3SWkj,V3RWkj,V3SWkj,PRWport,PSWport,VRWport,VSWport,VescRWkj,
                         VescSWkj,Pescape.trans.RW,Pescape.trans.SW,Vescape.transRW,Vescape.transSW,
                         PtransRW,PtransSW,VtransRW,VtransSW,Pescape.fact.RW,Pescape.fact.SW,Vescape.fact.RW,   
                         Vescape.fact.SW,PfactRW,PfactSW,VfactRW,VfactSW,PfactRWnuts2,PfactSWnuts2,
                         VfactRWnuts2,VfactSWnuts2,Pescape.res.nuts2,Vescape.res.nuts2,Pres,Vres,
                         Pres.nuts2,Vres.nuts2,Pescape.FP.nuts2,Vescape.FP.nuts2,PFP,VFP,PFP.nuts2,VFP.nuts2)
  names(listoutput) = c("P3RWkj","P3SWkj","V3RWkj","V3SWkj","PRWport","PSWport","VRWport","VSWport","VescRWkj",
                      "VescSWkj","Pescape.trans.RW","Pescape.trans.SW","Vescape.transRW","Vescape.transSW",
                      "PtransRW","PtransSW","VtransRW","VtransSW","Pescape.fact.RW","Pescape.fact.SW","Vescape.fact.RW",   
                      "Vescape.fact.SW","PfactRW","PfactSW","VfactRW","VfactSW","PfactRWnuts2","PfactSWnuts2",
                      "VfactRWnuts2","VfactSWnuts2","Pescape.res.nuts2","Vescape.res.nuts2","Pres","Vres",
                      "Pres.nuts2","Vres.nuts2","Pescape.FP.nuts2","Vescape.FP.nuts2","PFP","VFP","PFP.nuts2","VFP.nuts2")
    
  listwoodflow[[t]] = listoutput
  
  # We update the volume of RW entering in each MS
  RWenterj = t(TRANSITIONRW) %*% RWexitj
  }
  
  # Probability of dissemination across all the steps (exchanges between countries including transportation + transformation):
  #--------------------------------------------------------------------------------------------------------------------------------
  

  PtransRWallsteps_1=1
  PtransSWallsteps_1=1
  PfactRWallsteps_1=1
  PfactRWnuts2allsteps_1 = 1
  PfactSWallsteps_1=1
  PfactSWnuts2allsteps_1 = 1
  PFPallsteps_1=1
  Presallsteps_1=1
  Presnuts2allsteps_1 = 1
  for (t in 1:nsteps)
  {
  PtransRWallsteps_1 = PtransRWallsteps_1*(1-listwoodflow[[t]]$PtransRW)
  PtransSWallsteps_1 = PtransSWallsteps_1*(1-listwoodflow[[t]]$PtransSW)
  PfactRWallsteps_1 = PfactRWallsteps_1*(1-listwoodflow[[t]]$PfactRW)
  PfactRWnuts2allsteps_1 = PfactRWnuts2allsteps_1*(1-listwoodflow[[t]]$PfactRWnuts2)
  PfactSWallsteps_1 = PfactSWallsteps_1*(1-listwoodflow[[t]]$PfactSW)
  PfactSWnuts2allsteps_1 = PfactSWnuts2allsteps_1*(1-listwoodflow[[t]]$PfactSWnuts2)
  PFPallsteps_1 = PFPallsteps_1*(1-listwoodflow[[t]]$PFP)
  Presallsteps_1 = Presallsteps_1*(1-listwoodflow[[t]]$Pres)
  Presnuts2allsteps_1 = Presnuts2allsteps_1*(1-listwoodflow[[t]]$Presnuts2)
  }
  
  # do for all elements in listout (stored in listwoodflow)
  results.on = list()
  # ports
  results.on[[1]] = data.frame(portname=probadispport[,1],P3RWkj) # probability that at least one transfer takes place at ports
  results.on[[2]] = data.frame(portname=probadispport[,1],P3SWkj) # probability that at least one transfer takes place at ports
  results.on[[3]] = data.frame(portname=probadispport[,1],V3RWkj) # number of propagules that transfer around ports from RW per port
  results.on[[4]] = data.frame(portname=probadispport[,1],V3SWkj) # number of propagules that transfer around ports from SW per port
  results.on[[5]] = PRWport # probability that at least one transfer takes place at ports per country
  results.on[[6]] = PSWport # probability that at least one transfer takes place at ports per country
  results.on[[7]] = VRWport # number of propagules that transfer around ports from RW per country 
  results.on[[8]] = VSWport # number of propagules that transfer around ports from SW per country
  results.on[[9]] = VescRWkj # number of propagules that escape around ports from RW per country
  results.on[[10]] = VescSWkj # number of propagules that escape around ports from SW per country
  # transport
  results.on[[11]] = data.frame(EUcountries,prob.t(nsteps, listwoodflow,"Pescape.trans.RW")) # probability that at least one escape takes place from RW during transport
  results.on[[12]] = data.frame(EUcountries,prob.t(nsteps, listwoodflow,"Pescape.trans.SW")) # probability that at least one escape takes place from SW during transport
  results.on[[13]] = data.frame(EUcountries,prob.t(nsteps, listwoodflow,"Vescape.transRW"))  # number of propagules that that escape from RW during transport
  results.on[[14]] = data.frame(EUcountries,prob.t(nsteps, listwoodflow,"Vescape.transSW")) # number of propagules that that escape from RW during transport
  results.on[[15]] = data.frame(EUcountries,prob.t(nsteps, listwoodflow,"PtransRW")) # probability that at least one transfer takes place from RW 
  results.on[[16]] = data.frame(EUcountries,prob.t(nsteps, listwoodflow,"PtransSW")) # probability that at least one transfer takes place from RW 
  results.on[[17]] = data.frame(EUcountries,prob.t(nsteps, listwoodflow,"VtransRW")) # number of propagules that that transfer from RW during transport
  results.on[[18]] = data.frame(EUcountries,prob.t(nsteps, listwoodflow,"VtransSW")) # number of propagules that that transfer from SW during transport
  # factories RW & SW
  results.on[[19]] = data.frame(nuts2factories[,1],prob.t(nsteps, listwoodflow,"Pescape.fact.RW")) # probability that at least one escape takes place from RW near factories at NUTS2 level
  results.on[[20]] = data.frame(nuts2factories[,1],prob.t(nsteps, listwoodflow,"Pescape.fact.SW")) # probability that at least one escape takes place from SW near factories at NUTS2 level
  results.on[[21]] = data.frame(nuts2factories[,1],prob.t(nsteps, listwoodflow,"Vescape.fact.RW"))  # number of propagules that escape from RW near factories at NUTS2 level  
  results.on[[22]] = data.frame(nuts2factories[,1],prob.t(nsteps, listwoodflow,"Vescape.fact.SW"))  # number of propagules that escape from SW near factories at NUTS2 level    
  results.on[[23]] = data.frame(EUcountries,prob.t(nsteps, listwoodflow,"PfactRW")) # probability that at least on propagule transfers from RW near factories at country level  
  results.on[[24]] = data.frame(EUcountries,prob.t(nsteps, listwoodflow,"PfactSW")) # probability that at least on propagule transfers from SW near factories at country level  
  results.on[[25]] = data.frame(EUcountries,prob.t(nsteps, listwoodflow,"VfactRW")) # number of propagules that transfer from RW near factories at country level  
  results.on[[26]] = data.frame(EUcountries,prob.t(nsteps, listwoodflow,"VfactSW")) # number of propagules that transfer from SW near factories at country level  
  results.on[[27]] = data.frame(nuts2factories[,1],prob.t(nsteps, listwoodflow,"PfactRWnuts2")) # probability that at least one transfer takes place from RW near factories at NUTS2 level
  results.on[[28]] = data.frame(nuts2factories[,1],prob.t(nsteps, listwoodflow,"PfactSWnuts2")) # probability that at least one transfer takes place from SW near factories at NUTS2 level
  results.on[[29]] = data.frame(nuts2factories[,1],prob.t(nsteps, listwoodflow,"VfactRWnuts2")) # number of propagules that transfer from RW near factories at nuts2 level  
  results.on[[30]] = data.frame(nuts2factories[,1],prob.t(nsteps, listwoodflow,"VfactSWnuts2")) # number of propagules that transfer from SW near factories at nuts2 level  
  # residues
  results.on[[31]] = data.frame(nuts2factories[,1],prob.t(nsteps, listwoodflow,"Pescape.res.nuts2")) # probability that at least one propagule escape takes place from residues near factories at NUTS2 level
  results.on[[32]] = data.frame(nuts2factories[,1],prob.t(nsteps, listwoodflow,"Vescape.res.nuts2"))  # number of propagules escapes from residues at nuts2 level
  results.on[[33]] = data.frame(EUcountries,prob.t(nsteps, listwoodflow,"Pres")) # probability that at least one transfer takes place from residues near factories at country level
  results.on[[34]] = data.frame(EUcountries,prob.t(nsteps, listwoodflow,"Vres")) # number of propagules that transfer from residues at country level
  results.on[[35]] = data.frame(nuts2factories[,1],prob.t(nsteps, listwoodflow,"Pres.nuts2")) # probability that at least one propagule transfers from residues near factories at country level
  results.on[[36]] = data.frame(nuts2factories[,1],prob.t(nsteps, listwoodflow,"Vres.nuts2")) # number of propagules that transfer from residues at country level
  # final product
  results.on[[37]] = data.frame(nuts2factories[,1],prob.t(nsteps, listwoodflow,"Pescape.FP.nuts2")) # probability that at least one escape will take place from FP at NUTS2 level
  results.on[[38]] = data.frame(nuts2factories[,1],prob.t(nsteps, listwoodflow,"Vescape.FP.nuts2")) # number of propagules that escapes from FP at NUTS2 level
  results.on[[39]] = data.frame(EUcountries,prob.t(nsteps, listwoodflow,"PFP"))  # probability that at least one propagule transfers from FP at country level
  results.on[[40]] = data.frame(EUcountries,prob.t(nsteps, listwoodflow,"VFP")) # number of propagules that transfer from FP at country level
  results.on[[41]] = data.frame(nuts2factories[,1],prob.t(nsteps, listwoodflow,"PFP.nuts2")) # probability that at least one propagule escapes from FP at NUTS2 level
  results.on[[42]] = data.frame(nuts2factories[,1],prob.t(nsteps, listwoodflow,"VFP.nuts2")) # number of propagules that transfer from FP at NUTS2 level
  # nuts2 size
  results.on[[43]] = nuts2host.orig[,1:2] # nuts code and area of this nuts region (km2)
  names(results.on) = c("P3RWkj","P3SWkj","V3RWkj","V3SWkj","PRWport","PSWport","VRWport","VSWport","VescRWkj",
                        "VescSWkj","Pescape.trans.RW","Pescape.trans.SW","Vescape.transRW","Vescape.transSW",
                        "PtransRW","PtransSW","VtransRW","VtransSW","Pescape.fact.RW","Pescape.fact.SW","Vescape.fact.RW",   
                        "Vescape.fact.SW","PfactRW","PfactSW","VfactRW","VfactSW","PfactRWnuts2","PfactSWnuts2",
                        "VfactRWnuts2","VfactSWnuts2","Pescape.res.nuts2","Vescape.res.nuts2","Pres","Vres",
                        "Pres.nuts2","Vres.nuts2","Pescape.FP.nuts2","Vescape.FP.nuts2","PFP","VFP","PFP.nuts2","VFP.nuts2","nuts2area")
  
  
  
  # Overall probability of dissemination across all the steps, including dissemination around ports:
  #-------------------------------------------------------------------------------------------------
  Pport = 1-(1-PRWport)*(1-PSWport)
  probaallsteps = 1-(1-Pport)*PtransRWallsteps_1*PtransSWallsteps_1*PfactRWallsteps_1*PfactSWallsteps_1*Presallsteps_1*PFPallsteps_1
  
  
  #------------------------------------------------------
  # plotting the results:
  
  result.MS<-data.frame(names(listwoodflow[[1]]$Pres),
                        PRWport,PSWport,
                        (1-PtransRWallsteps_1),
                        (1-PtransSWallsteps_1),
                        (1-PfactRWallsteps_1),
                        (1-PfactSWallsteps_1),
                        (1-PFPallsteps_1),
                        (1-Presallsteps_1),
                        probaallsteps)
  

  names(result.MS) = c("MS","Pport_RW","Pport_SW","PtransRW","PtransSW","PfactRW","PfactSW","PFP","Pres","TOTAL")

results.MS_E$TOTAL = apply(result.MS_E[,-1],1,sum,na.rm=TRUE)
  
  if (sensitivity){out.sens = data.frame(parameter,old.param,new.param,rel.change); names(out.sens)=c("param","old value","new value","rel.change");
                   results = list(); results[[1]] = result.MS; results[[2]]=results.on} else {
                       results = list(); results[[1]] = result.MS; results[[2]]= results.on}
  return(results)
  
}

#### END MODEL

  

# RUN the model (using default settings)
result.NUTS2 = Conif_wood_PM()
# change one of the arguments: for example the number of times wood is recirculated in europe
result.NUTS2 = Conif_wood_PM(nsteps=3)
# the object result.NUTS2 contains the model output
# first element of the list contains probability that at least one propagule will transfer to a host at country level
# the second element of the list contains several measures: see lines 769-822 for what each element in this list represents

# the ports are located in a NUTS2 region. To assign the escapes and transfers of propagules of ports to a nuts2 region on can use the file:
port.NUTS = read.csv(paste(persofolder,"port.NUTS.combination.csv",sep=""))
  
  