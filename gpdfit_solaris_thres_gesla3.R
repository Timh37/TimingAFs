# Victor Malagon Santos, March 2022;

# GPD threshold selection based on Solari et al. (2017).

## load required packages ------------------------------------------------------
# options(timeout=9999999)
# devtools::install_github("rjaneUCF/MultiHazard")
library(MultiHazard)
library(texmex)
# library(interp)
library(ncdf4)

# Remove loaded data (if any)
rm(list = (ls()))

## Prepare data for loading ----------------------------------------------------
ndir_data <- "/Users/vmalagonsantos/GitHub/extremes_timing/DATA/19_4_2022_detrended_deseason_coastal_minyr30"
nfiles <- list.files( ndir_data,
                      pattern = "daily") # get files names

# Output file names
nfiles_out <- gsub("\\..*","",nfiles) # remove string after '.'
nfiles_out <- substring(nfiles_out,11,
                        nchar(nfiles_out)) # ignore first 10 characters (daily_max_)

mquant <- .95 # min quantile from which to decluster and perform threshold selection

## loop through files - Perform analysis ---------------------------------------
for (i in 1:length(nfiles)){
print(i)
  
# define datafile path
DataFile <- paste( ndir_data,'/', nfiles[i], sep="")

# load data
data_hourly <- read.csv( DataFile, header = FALSE, sep = ",", skip = 1)
colnames(data_hourly) <- c("Date","Value")
# data_hourly[,2] <- data_hourly[,2] + 8

# Decluster data
Decluster_SW_Thres<-function(Data,Window_Width,u){
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
} # this function applies declustering to values above mquant
data.declust.SW95 <- Decluster_SW_Thres(Data=data_hourly, Window_Width=3, u= mquant)

# Test threshold using Solari's method - Plot sensitivity analysis 
cidr <- getwd() 
mkfldr <- "/OUTPUT/PLOTS"
dir.create(file.path(cidr, mkfldr), recursive = TRUE) # create output folder in working directory
mypath <- file.path(cidr, mkfldr,paste("sens_", nfiles_out[i], ".jpg", sep = "")) # path to store threshold sensitivity plots

jpeg(file=mypath)
remove(data_Solari95)
# tryCatch({
data_Solari95<-GPD_Threshold_Solari(Event=data.declust.SW95$Declustered,
                                       Data=na.omit(data_hourly[,-1]),N_Sim = 50,
                                       mu= 365.25,  Min_Quantile = mquant)

dev.off()

# }, error=function(e){})

if (exists('S22_data_Solari95')){

# The candidate threshold according to the Solari approach can be extracted as follows.
# GPD_Threshold_Solari_Sel() provides a closer examination of the GPD fit associated with a candidate threshold Thres by comparing the parameter estimates and estimated return levels with those obtained with the other thresholds tested in GPD_Threshold_Solari().
# grey histrograms: MLE of all selected threshold; Black transparent histograms: boostrapped uncertainties of selected threshold

nbootf <- 10^4 # number fo boostrapping samples for uncertainty analysis
mypath <- file.path(cidr, mkfldr,paste("gpdfit_", nfiles_out[i], ".jpg", sep = "")) # path for GPD plots
jpeg(file=mypath)

remove(GPD_fit_out)
#tryCatch({

GPD_fit_out <- GPD_Threshold_Solari_Sel(Event=data.declust.SW95$Declustered,Data=na.omit(data_hourly[,-1]),
                         Solari_Output=data_Solari95,
                         Thres=data_Solari95$Candidate_Thres,N_Sim = nbootf,
                         Alpha = 0.05,
                         mu= 365.25,
                         RP_Min=5)
dev.off()

#}, error=function(e){})


vthres <- GPD_fit_out$Estimate[3] # find treshold
ind95 <- which.min(abs(data_Solari95$Thres_Candidate - GPD_fit_out$Estimate[3])) #find index

## Store data in netCDF files --------------------------------------------------
if (exists('GPD_fit_out')){
  
  
  
ncname <- nfiles_out[i]
ncfname <- paste(cidr,"/OUTPUT/gpd_solary_", ncname, ".nc", sep="")
dname <- "gpd"  

GPDest <- c(10,20,30,40,50,60,70,80,90,100,200,300,400,500,600,700,800,900,1000)
quantgpd <- c(2.5, 97.5)

# define dimensions
dec_dim <- ncdim_def("id_timeseries","number_of_observations_in_original_time_series",as.double(1:length(data.declust.SW95$Declustered))) 
rp_dim <- ncdim_def("return_period","years",as.double(GPDest))
quant_dim <- ncdim_def("quantiles","quantiles", as.double(quantgpd) )
est_dim <- ncdim_def("est_d","est_dimension_for_selected_threshold", as.numeric(1) )
nboot_dim <- ncdim_def("nboot","number_of_boostrapped_samples", as.double(1:nbootf))

# define variables
fillvalue <- 1e32

dlname <- "declustered_exceedances"
dec.def <- ncvar_def("dec_ex","m",list(dec_dim),fillvalue,dlname,prec="single")

dlname <- "return_period_estimates"
RPest.def <- ncvar_def("est","m",list(rp_dim),fillvalue,dlname,prec="single")

dlname <- "return_period_CIlower2_5"
RPlow.def <- ncvar_def("lower","m",list(rp_dim),fillvalue,dlname,prec="single")

dlname <- "return_period_CIupper97_5"
RPup.def <- ncvar_def("upper","m",list(rp_dim),fillvalue,dlname,prec="single")

dlname <- "return_period_boot"
RPboot.def <- ncvar_def("boot","m",list(rp_dim,nboot_dim),fillvalue,dlname,prec="single")

dlname <- "xi"
xi.def <- ncvar_def("xi","-",list(est_dim),fillvalue,dlname,prec="single")

dlname <- "xi_confidence_interval"
xici.def <- ncvar_def("xi_ci","-",list(quant_dim),fillvalue,dlname,prec="single")

dlname <- "xi_boot"
xib.def <- ncvar_def("xi_boot","-",list(nboot_dim),fillvalue,dlname,prec="single")

dlname <- "sigma"
sigma.def <- ncvar_def("sigma","-",list(est_dim),fillvalue,dlname,prec="single")

dlname <- "sigma_confidence_interval"
sigmaci.def <- ncvar_def("sigma_ci","-",list(quant_dim),fillvalue,dlname,prec="single")

dlname <- "sigma_boot"
sigmab.def <- ncvar_def("sigma_boot","-",list(nboot_dim),fillvalue,dlname,prec="single")

dlname <- "thres"
thres.def <- ncvar_def("thres","m",list(est_dim),fillvalue,dlname,prec="single")

dlname <- "rate"
rate.def <- ncvar_def("rate","-",list(est_dim),fillvalue,dlname,prec="single")

dlname <- "percentile"
prc.def <- ncvar_def("prctile","-",list(est_dim),fillvalue,dlname,prec="single")

# create netCDF file and put arrays
ncout <- nc_create(ncfname,list(dec.def,
                                RPest.def,
                                RPlow.def,
                                RPup.def,
                                RPboot.def,
                                xi.def,
                                xici.def,
                                xib.def,
                                sigma.def,
                                sigmaci.def,
                                sigmab.def,
                                thres.def,
                                rate.def,
                                prc.def),
                   force_v4=TRUE)

# put variables
ncvar_put(ncout,dec.def,data.declust.SW95$Declustered)
ncvar_put(ncout,RPest.def,data.frame(GPD_fit_out$Estimate)[5:23,1])
ncvar_put(ncout,RPlow.def,GPD_fit_out$CI_Lower[5:23])
ncvar_put(ncout,RPup.def,GPD_fit_out$CI_Upper[5:23])
ncvar_put(ncout,RPboot.def,GPD_fit_out$BOOT[,5:23])
ncvar_put(ncout,xi.def,data.frame(GPD_fit_out$Estimate[1]))
ncvar_put(ncout,xici.def,data.frame(GPD_fit_out$CI_Lower[1],
                                  GPD_fit_out$CI_Upper[1]))
ncvar_put(ncout,xib.def,GPD_fit_out$BOOT[,1])
ncvar_put(ncout,sigma.def,data.frame(GPD_fit_out$Estimate[2]))
ncvar_put(ncout,sigmaci.def,data.frame(GPD_fit_out$CI_Lower[2],
                                     GPD_fit_out$CI_Upper[2]))
ncvar_put(ncout,sigmab.def,GPD_fit_out$BOOT[,2])
ncvar_put(ncout,thres.def,data.frame(GPD_fit_out$Estimate[3]))
ncvar_put(ncout,rate.def,data.frame(GPD_fit_out$Estimate[4]))
ncvar_put(ncout,prc.def,data.frame(data_Solari95$Thres_Candidate_Quantile[ind95][1]))

ncatt_put(ncout,0,"title", paste("Theshold selection based on Solari et al. (2017) and associated GPD fit appplied to" , nfiles_out[i]))

# close the file, writing data to disk
nc_close(ncout)
}
}

} 

