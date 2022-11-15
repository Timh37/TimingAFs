# Victor Malagon Santos, March 2022;

# Compute GPD parameters based on a median threshold percentile of 0.988.

# load required packages.
# options(timeout=9999999)
# devtools::install_github("rjaneUCF/MultiHazard")
library(MultiHazard)
library(texmex)
# library(interp)
library(ncdf4)

# Remove loaded data (if any)
rm(list = (ls()))

# Load data
ndir_data <- "/Users/vmalagonsantos/GitHub/extremes_timing/DATA/19_4_2022_detrended_deseason_coastal_minyr30"
nfiles <- list.files( ndir_data,
                      pattern = "daily")

nfiles_out <- gsub("\\..*","",nfiles) # remove string after '.'
nfiles_out <- substring(nfiles_out,11,nchar(nfiles_out)) #ignore first 10 characters (daily_max_)

scale  <- matrix(NA, length(nfiles)) # create matrix to add selected threshold through the loop
shape  <- matrix(NA, length(nfiles))  
thres  <- matrix(NA, length(nfiles)) 
rate  <- matrix(NA, length(nfiles)) 


for (i in 1:length(nfiles)){ # loop thru files
print(i)
DataFile <- paste( ndir_data,'/', nfiles[i], sep="") # define datafile path

date_obs<- read.csv( DataFile, header = FALSE, sep = ",", skip = 1) # load data
colnames(date_obs) <- c("Date","Value")

# if(i == 455){ # manual fix for Tenerife's tide gauge.
# date_obs[,2] <- date_obs[,2] + 8
# }

# Decluster data
Decluster_SW_Thres95<-function(Data,Window_Width,u){
  #N<-1:length(Data[,1])
  N<-which(Data[,2]>quantile(na.omit(Data[,2]),u))
  EVENT.INDEX<-rep(0,length(Data[,1]))
  EVENT<-rep(0,length(Data[,1]))
  ID1   = 1
  for(i in 1:length(N)){                    
    ID  = N[which(N>=max(1,N[i]-(Window_Width-1)/2) & N<=min(max(N),N[i]+(Window_Width-1)/2))]
    AUX = Data[ID,2]
    MAX = max(AUX,na.rm=T)  
    IDMAX = ifelse(MAX==-Inf,1,which(AUX==MAX))
    if(is.na(MAX)==F & ID[IDMAX]==N[i]){
      EVENT.INDEX[ID1]<-N[i]
      EVENT[ID1]<-MAX
      ID1 = ID1+1  
    }  
  }
  EVENT<-EVENT[EVENT.INDEX>0]
  EVENT.INDEX<-EVENT.INDEX[EVENT.INDEX>0]
  Declustered<-rep(NA,nrow(Data))
  Declustered[EVENT.INDEX]<-Data[EVENT.INDEX,2]
  res<-list("N"=N,"Detrend"=Data[,2],"Declustered"=Declustered,"EventID"=EVENT.INDEX)
  return(res)
} # edit function to include desired u
date_obs_dec <- Decluster_SW_Thres95(Data=date_obs,Window_Width=3, u=0.95) # u=95 declusters data only above the 95th percentile


remove(GPD_fit_out)
#tryCatch({

GPD_MLE<-function(Data){
  mod = evm(na.omit(Data), th = quantile(Data, 0))
  xi = mod[[1]][2]
  sigma = exp(mod[[1]][1])
  u = as.numeric(quantile(Data,0))
  res<-list(xi=xi,sigma=sigma,u=u)
  return(res)
} # edit GPD_MLE function to consider only the 98.8th percentile

obs <- na.omit(date_obs[,-1])
Event = date_obs_dec$Declustered
th <- quantile(obs, 0.988, names = FALSE)
GPD_fit_out = GPD_MLE(na.omit(Event[Event > th]))


GPD_MLE_Boot<-function(Data){
  mod        = GPD_MLE(Data)
  xi         = mod$xi
  sigma      = mod$sigma
  u          = mod$u
  res        = list(xi=xi,sigma=sigma,u=u)
  return(res)
} # bootstrapping function to compute uncertainties

N_Sim = 10^4 # number of boostrapped sampled
N_Years = length(obs)/365.25

rate[i] = length(na.omit(Event[Event > th]))/N_Years
scale[i]  <- GPD_fit_out$sigma # create matrix to add selected threshold thru the loop
shape[i]  <- GPD_fit_out$xi
thres[i] <- GPD_fit_out$u

Event = na.omit(Event)
Event = sort(Event)
Exceedence = Event[Event > th]
Estimate = unlist(GPD_MLE_Boot(Exceedence))
BOOT = array(0, dim = c(N_Sim, length(Estimate)))
for (I in 1:N_Sim) {
  BOOT[I, ] = unlist(GPD_MLE_Boot(rgpd(length(Exceedence),
                                       xi = as.numeric(Estimate[1]), sigma = as.numeric(Estimate[2]),
                                       u = as.numeric(Estimate[3]))))
  while (BOOT[I, 1] == 0) {
    try({
      BOOT[I, ] = unlist(GPD_MLE_Boot(sample(Exceedence,
                                             length(Exceedence), replace = T)))
    }, silent = FALSE)
  }
} # get parameter estimates for each boostrapped sample

colnames(BOOT) <- c('Shape','Scale','Threshold') # store boostrapeed estimates


## save output in netCDF file for each tide gauge
ncname <- nfiles_out[i]
ncfname <- paste("/Users/vmalagonsantos/GitHub/extremes_timing/OUTPUT/t988/gpd_988_", ncname, ".nc", sep="")
dname <- "gpd"  

# define dimensions
est_dim <- ncdim_def("est_d","est_dimension_for_selected_threshold", as.numeric(1) )
nboot_dim <- ncdim_def("nboot","number_of_boostrapped_samples", as.double(1:N_Sim))

# define variables
fillvalue <- 1e32

dlname <- "sigma"
sigma.def <- ncvar_def("sigma","-",list(est_dim),fillvalue,dlname,prec="single")

dlname <- "sigma_boot"
sigmab.def <- ncvar_def("sigma_boot","-",list(nboot_dim),fillvalue,dlname,prec="single")

dlname <- "xi"
xi.def <- ncvar_def("xi","-",list(est_dim),fillvalue,dlname,prec="single")

dlname <- "xi_boot"
xib.def <- ncvar_def("xi_boot","-",list(nboot_dim),fillvalue,dlname,prec="single")

dlname <- "thres"
thres.def <- ncvar_def("thres","m",list(est_dim),fillvalue,dlname,prec="single")

# create netCDF file and put arrays
ncout <- nc_create(ncfname,list(sigma.def,
                                sigmab.def,
                                xi.def,
                                xib.def,
                                thres.def
                                ),
                   force_v4=TRUE)

# put variables

ncvar_put(ncout,sigma.def,data.frame(GPD_fit_out$sigma))
ncvar_put(ncout,sigmab.def,BOOT[,2])
ncvar_put(ncout,xi.def,data.frame(GPD_fit_out$xi))
ncvar_put(ncout,xib.def,BOOT[,1])
ncvar_put(ncout,thres.def,data.frame(GPD_fit_out$u))

# add global attributes
ncatt_put(ncout,0,"title", paste("GPD fit (0.988) appplied to" , nfiles_out[i]))

# close the file, writing data to disk
nc_close(ncout)

} 

# save best parameters estimates for all tide gauges in a csv file
data_out <- data.frame(nfiles_out, scale, shape, thres, rate)
colnames(data_out) <- c('Station name','Scale','Shape','Threshold', 'Rate')
write.csv(data_out, file = '/Users/vmalagonsantos/GitHub/extremes_timing/OUTPUT/t988/gpd_par_988.csv',row.names = FALSE)
