#The R script of methods used to compute non-stationary gev/gpd for tasmax variable
#The script contains : 
#   - Selection of data
#   - Computation of Non-stationary gev with
#       - Global Mean Surface Temperature (GMST) as covariate 
#       - European Mean Surface Temperature (EMST) as covariate
#   - Computation of return levels 
#Return all the data as the gev as well as the one corresponding to the final one after selection

library(tidyverse)
library(ggplot2)
library(ncdf4)
library(extRemes)
source("caldat.R")
source("gev_gdp.R")
library(dplyr)
library(patchwork)
library(zoo)
library(Dict)
fig_path = '/home/estimr3/lhasbini/ANDRA/R_fig/'
data_path = '/home/estimr3/lhasbini/ANDRA/R_data/'
CMIP6_path = "/home/estimr3/lhasbini/CMIP6/"

#Select and open data 
args=(commandArgs(TRUE))
print(args)
i=1
if(length(args)>0){
  return_period = as.logical(args[i]);i=i+1 ## Return period in year
  variable = args[i];i=i+1 ## Save the name of the variable
  method = args[i];i=i+1  ## R2D2 or CDFt method
  return_horizon = c(as.integer(args[i]));i=i+1 ## Save the horizon
}else{ ## Default option 
    return_period= 100
    variable = 'tasmax'
    return_horizon = c(2100)
    method='R2D2'
}
rl = TRUE

if(method=='R2D2'){
  CMIP6_Ajust_path = "/home/estimr3/yiou/ANDRA2022/CMIP6_R2D2/GrandEst/"
  groups = c('EC-Earth-Consortium', 'MPI-M', 'AS-RCEC', 'BCC', 'CCCma', 'CNRM-CERFACS', 'CSIRO',  'IPSL', 'MRI', 'NOAA-GFDL')  
  models = c('EC-Earth3', 'MPI-ESM1-2-LR', 'TaiESM1', 'BCC-CSM2-MR', 'CanESM5', 'CNRM-ESM2-1', 'ACCESS-ESM1-5', 'IPSL-CM6A-LR', 'MRI-ESM2-0', 'GFDL-ESM4') 
} else if(method=='CDFt'){
  CMIP6_Ajust_path = "/home/estimr3/yiou/ANDRA2022/CMIP6_CDFt-L-1V-0L/GrandEst/"
  groups = c('EC-Earth-Consortium', 'MPI-M', 'AS-RCEC', 'BCC', 'CCCma', 'CNRM-CERFACS', 'CSIRO',  'IPSL', 'MRI')  
  models = c('EC-Earth3', 'MPI-ESM1-2-LR', 'TaiESM1', 'BCC-CSM2-MR', 'CanESM5', 'CNRM-ESM2-1', 'ACCESS-ESM1-5', 'IPSL-CM6A-LR', 'MRI-ESM2-0') 
  
}

## Initializatin of global parameters 
stations = c('Bure', 'St Dizier','Cirfontaines')
scenarios = c('historical', 'ssp126', 'ssp245', 'ssp370', 'ssp585')

tas_global = data_CMIP6_fldmean_yearmean_to_df(models, CMIP6_path, 'tas', scenarios, 'Global')
tas_Europe = data_CMIP6_fldmean_yearmean_to_df(models, CMIP6_path, 'tas', scenarios, 'Europe')

daily_var_CMIP6_GE = data_CMIP6_GE_to_df(models, CMIP6_path, variable, scenarios)

daily_var_adjust_stations = data_adjust_to_df(groups, models, CMIP6_Ajust_path, variable, scenarios, stations, method)
if(variable=='tasmax'){
  yearly_var_CMIP6_GE = data_CMIP_to_varn_max(daily_var_CMIP6_GE)
  yearly_var_adjust_stations = data_to_varn_max(daily_var_adjust_stations)
} else if(variable=='tasmin'){
  yearly_var_CMIP6_GE = data_to_varn_min(daily_var_CMIP6_GE)
  yearly_var_adjust_stations = data_to_varn_min(daily_var_adjust_stations)
}

ns_models = c('M110','M000','M-100', 'M-110', 'M100')
alpha_llh = 0.1

#Compute GEV with GMST and EMST as covariate 

##### Mutate data #####
spar = 0.9
tas_global_spline = data.frame(matrix(ncol = ncol(tas_global)+1, nrow = 0))
colnames(tas_global_spline) = append(colnames(tas_global),c('covariate'))
for(mod in models){
  for(sc in c('ssp126', 'ssp245', 'ssp370', 'ssp585')){
    data_loop = tas_global%>%filter(scenario%in%c('historical', sc), CMIP_model==mod)
    data_loop = data_loop %>% mutate(var = scale(var))
    model_spline_TX_GMST_loop = smooth.spline(x=data_loop$year, y = data_loop$var, spar=spar)
    data_loop$covariate  = predict(model_spline_TX_GMST_loop, data_loop$year)$y
    data_loop$scenario = sc
    tas_global_spline = rbind(tas_global_spline, data_loop)
  }
  sc = 'historical'
  data_loop = tas_global%>%filter(scenario== sc, CMIP_model==mod)
  data_loop_scenarios = tas_global_spline%>%filter(scenario%in%c('ssp126', 'ssp245', 'ssp370', 'ssp585'), CMIP_model==mod, year<2015)
  data_loop = merge(data_loop, aggregate(covariate~year, data_loop_scenarios, mean))
  data_loop$scenario=sc
  tas_global_spline = rbind(tas_global_spline, data_loop)
}

tas_europe_spline = data.frame(matrix(ncol = ncol(tas_Europe)+1, nrow = 0))
colnames(tas_europe_spline) = append(colnames(tas_Europe),c('covariate'))
for(mod in models){
  for(sc in c('ssp126', 'ssp245', 'ssp370', 'ssp585')){
    data_loop = tas_Europe%>%filter(scenario%in%c('historical', sc), CMIP_model==mod)
    data_loop = data_loop %>% mutate(var = scale(var))
    model_spline_TX_EMST_loop = smooth.spline(x=data_loop$year, y = data_loop$var, spar=spar)
    data_loop$covariate  = predict(model_spline_TX_EMST_loop, data_loop$year)$y
    data_loop$scenario = sc
    tas_europe_spline = rbind(tas_europe_spline, data_loop)
  }
  sc = 'historical'
  data_loop = tas_Europe%>%filter(scenario== sc, CMIP_model==mod)
  data_loop_scenarios = tas_europe_spline%>%filter(scenario%in%c('ssp126', 'ssp245', 'ssp370', 'ssp585'), CMIP_model==mod, year<2015)
  data_loop = merge(data_loop, aggregate(covariate~year, data_loop_scenarios, mean))
  data_loop$scenario=sc
  tas_europe_spline = rbind(tas_europe_spline, data_loop)
}

tas_global_spline_mutate = tas_global_spline %>% select(-var)
tas_europe_spline_mutate = tas_europe_spline %>% select(-var)
yearly_var_adjust_stations_GMST = merge(yearly_var_adjust_stations, tas_global_spline_mutate, by=c("year", "scenario", "CMIP_model"))
yearly_var_adjust_stations_EMST = merge(yearly_var_adjust_stations, tas_europe_spline_mutate, by=c("year", "scenario", "CMIP_model"))
yearly_var_CMIP_adjust_stations_GMST = merge(yearly_var_CMIP6_GE, tas_global_spline_mutate, by=c("year", "scenario", "CMIP_model"))
yearly_var_CMIP_adjust_stations_EMST = merge(yearly_var_CMIP6_GE, tas_europe_spline_mutate, by=c("year", "scenario", "CMIP_model"))

#### Computation of GEV & Return levels ####

### Note that here, confidence interval are computed thanks to bootstraping with the parameter 'confidence_internval'
ns_model_GMST.gev = gev_ns_rl_param_fct(yearly_var_adjust_stations_GMST, models, ns_models, confidence_interval='boot', sigma_exp=TRUE, rl=rl, return_period, return_horizon, variable)
if(rl){
  ns_model_GMST.gev.params = ns_model_GMST.gev$params
  ns_model_GMST.gev.rl = ns_model_GMST.gev$return_levels

  saveRDS(ns_model_GMST.gev.params, file=paste(data_path, "gev_ns_",variable,"_adjust_",method,"_params_GMST.Rda", sep=""))
  saveRDS(ns_model_GMST.gev.rl, file=paste(data_path, "gev_ns_",variable,"_adjust_",method,'_',return_period,"y_rl_",return_horizon[1],"_GMST.Rda", sep=""))
} else{
  ns_model_GMST.gev.params = ns_model_GMST.gev
  saveRDS(ns_model_GMST.gev.params, file=paste(data_path, "gev_ns_",variable,"_adjust_",method,"_params_GMST.Rda", sep=""))
}

ns_model_EMST.gev = gev_ns_rl_param_fct(yearly_var_adjust_stations_EMST, models, ns_models, confidence_interval='boot', sigma_exp=TRUE, rl=rl, return_period, return_horizon, variable)

if(rl){
  ns_model_EMST.gev.params = ns_model_EMST.gev$params
  ns_model_EMST.gev.rl = ns_model_EMST.gev$return_levels

  saveRDS(ns_model_EMST.gev.params, file=paste(data_path, "gev_ns_",variable,"_adjust_",method,"_params_EMST.Rda", sep=""))
  saveRDS(ns_model_EMST.gev.rl, file=paste(data_path, "gev_ns_",variable,"_adjust_",method,'_',return_period,"y_rl_",return_horizon[1],"_EMST.Rda", sep=""))

} else {
  ns_model_EMST.gev.params = ns_model_EMST.gev
  saveRDS(ns_model_EMST.gev.params, file=paste(data_path, "gev_ns_",variable,"_adjust_",method,"_params_EMST.Rda", sep=""))
}


# #### Computation without using the historical data

# ns_model_GMST_no_hist.gev = gev_ns_rl_param_fct(yearly_var_adjust_stations_GMST, models, ns_models, confidence_interval='boot', sigma_exp=TRUE, rl=rl, return_period, return_horizon, variable, with_historical=FALSE)
# if(rl){
#   ns_model_GMST_no_hist.gev.params = ns_model_GMST_no_hist.gev$params
#   ns_model_GMST_no_hist.gev.rl = ns_model_GMST_no_hist.gev$return_levels
#   
#   saveRDS(ns_model_GMST_no_hist.gev.params, file=paste(data_path, "gev_ns_no_hist_",variable,"_adjust_",method,"_params_GMST.Rda", sep=""))
#   saveRDS(ns_model_GMST_no_hist.gev.rl, file=paste(data_path, "gev_ns_no_hist_",variable,"_adjust_",method,"_",return_period,"y_rl_",return_horizon[1],"_GMST.Rda", sep=""))
#   
# } else {
#   ns_model_GMST_no_hist.gev.params = ns_model_GMST_no_hist.gev
#   saveRDS(ns_model_GMST_no_hist.gev.params, file=paste(data_path, "gev_ns_no_hist_",variable,"_adjust_",method,"_params_GMST.Rda", sep=""))
# }
# 
# ns_model_EMST_no_hist.gev = gev_ns_rl_param_fct(yearly_var_adjust_stations_EMST, models, ns_models, confidence_interval='boot', sigma_exp=TRUE, rl=rl, return_period, return_horizon, variable,with_historical=FALSE)
# if(rl){
#   ns_model_EMST_no_hist.gev.params = ns_model_EMST_no_hist.gev$params
#   ns_model_EMST_no_hist.gev.rl = ns_model_EMST_no_hist.gev$return_levels
#   
#   saveRDS(ns_model_EMST_no_hist.gev.params, file=paste(data_path, "gev_ns_no_hist_",variable,"_adjust_",method,"_params_EMST.Rda", sep=""))
#   saveRDS(ns_model_EMST_no_hist.gev.rl, file=paste(data_path, "gev_ns_no_hist_",variable,"_adjust_",method,"_",return_period,"y_rl_",return_horizon[1],"_EMST.Rda", sep=""))
# } else {
#   ns_model_EMST_no_hist.gev.params = ns_model_EMST_no_hist.gev
#   saveRDS(ns_model_EMST_no_hist.gev.params, file=paste(data_path, "gev_ns_no_hist_",variable,"_adjust_",method,"_params_EMST.Rda", sep=""))
# }

#### Selection procedure ####

ns_model_GMST_final.gev = selection_models(daily_var_adjust_stations, yearly_var_adjust_stations_GMST, ns_model_GMST.gev.params, ns_rl = ns_model_GMST.gev.rl, variable)
ns_model_GMST_final.gev.params = ns_model_GMST_final.gev$params
ns_model_GMST_final.gev.rl = ns_model_GMST_final.gev$return_levels

ns_model_EMST_final.gev = selection_models(daily_var_adjust_stations, yearly_var_adjust_stations_EMST, ns_model_EMST.gev.params, ns_rl = ns_model_EMST.gev.rl, variable)
ns_model_EMST_final.gev.params = ns_model_EMST_final.gev$params
ns_model_EMST_final.gev.rl = ns_model_EMST_final.gev$return_levels

#Save final params and return levels
saveRDS(ns_model_GMST_final.gev.params, file=paste(data_path, "gev_final_ns_",variable,"_adjust_",method,"_params_GMST.Rda", sep=""))
saveRDS(ns_model_GMST_final.gev.rl, file=paste(data_path, "gev_final_ns_",variable,"_adjust_",method,"_",return_period,"y_rl_",return_horizon[1],"_GMST.Rda", sep=""))

saveRDS(ns_model_EMST_final.gev.params, file=paste(data_path, "gev_final_ns_",variable,"_adjust_",method,"_params_EMST.Rda", sep=""))
saveRDS(ns_model_EMST_final.gev.rl, file=paste(data_path, "gev_final_ns_",variable,"_adjust_",method,"_",return_period,"y_rl_",return_horizon[1],"_EMST.Rda", sep=""))

#### Computation without the most extreme points ####

# Data files without the most extreme
daily_var_adjust_stations_no_max = data.frame(daily_var_adjust_stations)
daily_var_adjust_stations_max = data.frame(matrix(ncol = ncol(daily_var_adjust_stations_no_max), nrow = 0))
colnames(daily_var_adjust_stations_max) = colnames(daily_var_adjust_stations_no_max)
daily_var_adjust_stations_diff_max = data.frame(matrix(ncol=4, nrow = 0))
colnames(daily_var_adjust_stations_diff_max) = c('CMIP_model', 'scenario', 'station', 'var')

yearly_var_adjust_stations_no_max = data.frame(yearly_var_adjust_stations)
yearly_var_adjust_stations_max = data.frame(matrix(ncol = ncol(yearly_var_adjust_stations_no_max), nrow = 0))
colnames(yearly_var_adjust_stations_max) = colnames(yearly_var_adjust_stations_no_max)
yearly_var_adjust_stations_diff_max = data.frame(matrix(ncol = 4, nrow = 0))
colnames(yearly_var_adjust_stations_diff_max) = c('CMIP_model', 'scenario', 'station', 'var')

for(CMIP_mod in models){
  for(sc in c('ssp126', 'ssp245', 'ssp370', 'ssp585')){
    for(st in stations){
      #Remove all the year corresponding to the daily maximum
      if(variable %in% c('tasmax')){
        df_max_loop = which(daily_var_adjust_stations_no_max$var == max(daily_var_adjust_stations %>% filter(station==st, scenario==sc, CMIP_model==CMIP_mod) %>% select(var))
                            & daily_var_adjust_stations_no_max$CMIP_model == CMIP_mod & daily_var_adjust_stations_no_max$station ==st & daily_var_adjust_stations_no_max$scenario==sc)
        yearly_var_max_loop = which(yearly_var_adjust_stations_no_max$var == max(yearly_var_adjust_stations %>% filter(station==st, scenario==sc, CMIP_model==CMIP_mod) %>% select(var))
                                    & yearly_var_adjust_stations_no_max$CMIP_model == CMIP_mod & yearly_var_adjust_stations_no_max$station ==st & yearly_var_adjust_stations_no_max$scenario==sc)
      }
      else if(variable %in% c('tasmin')){
        df_max_loop = which(daily_var_adjust_stations_no_max$var == min(daily_var_adjust_stations %>% filter(station==st, scenario==sc, CMIP_model==CMIP_mod) %>% select(var))
                            & daily_var_adjust_stations_no_max$CMIP_model == CMIP_mod & daily_var_adjust_stations_no_max$station ==st & daily_var_adjust_stations_no_max$scenario==sc)
        yearly_var_max_loop = which(yearly_var_adjust_stations_no_max$var == min(yearly_var_adjust_stations %>% filter(station==st, scenario==sc, CMIP_model==CMIP_mod) %>% select(var))
                                    & yearly_var_adjust_stations_no_max$CMIP_model == CMIP_mod & yearly_var_adjust_stations_no_max$station ==st & yearly_var_adjust_stations_no_max$scenario==sc)
      }

      daily_var_adjust_stations_max = rbind(daily_var_adjust_stations_max, daily_var_adjust_stations_no_max[c(df_max_loop), ])

      year_max_loop = daily_var_adjust_stations_no_max[df_max_loop,]$year
      df_year_max_loop = which(daily_var_adjust_stations_no_max$station == st &
                                 daily_var_adjust_stations_no_max$scenario == sc &
                                 daily_var_adjust_stations_no_max$CMIP_model == CMIP_mod &
                                 daily_var_adjust_stations_no_max$year == year_max_loop)
      #daily_var_adjust_stations_no_max = daily_var_adjust_stations_no_max[-c(df_year_max_loop),]
      daily_var_adjust_stations_no_max[c(df_year_max_loop),'var']=NaN

      #Remove yearly maximum
      yearly_var_adjust_stations_max = rbind(yearly_var_adjust_stations_max, yearly_var_adjust_stations_no_max[c(yearly_var_max_loop), ])
      #yearly_var_adjust_stations_no_max = yearly_var_adjust_stations_no_max[-c(yearly_var_max_loop), ]
      yearly_var_adjust_stations_no_max[c(yearly_var_max_loop), 'var'] = NaN

      if(variable %in% c('tasmax')){
        diff_loop = max(yearly_var_adjust_stations %>% filter(station==st, scenario==sc, CMIP_model==CMIP_mod) %>% select(var)) -
          max(yearly_var_adjust_stations_no_max %>% filter(station==st, scenario==sc, CMIP_model==CMIP_mod) %>% select(var))
      }
      else if(variable %in% c('tasmin')){
        diff_loop = min(yearly_var_adjust_stations %>% filter(station==st, scenario==sc, CMIP_model==CMIP_mod) %>% select(var)) -
          min(yearly_var_adjust_stations_no_max %>% filter(station==st, scenario==sc, CMIP_model==CMIP_mod) %>% select(var))
      }

      yearly_var_adjust_stations_diff_max = rbind(yearly_var_adjust_stations_diff_max, data.frame(CMIP_model=CMIP_mod,
                                                                                  scenario=sc,
                                                                                  station=st,
                                                                                  var = diff_loop))
    }
  }
}

yearly_var_adjust_stations_no_max_GMST = merge(yearly_var_adjust_stations_no_max, tas_global_spline_mutate, by=c("year", "scenario", "CMIP_model"))
yearly_var_adjust_stations_no_max_EMST = merge(yearly_var_adjust_stations_no_max, tas_europe_spline_mutate, by=c("year", "scenario", "CMIP_model"))

# GEV + save params

ns_model_GMST_no_max.gev = gev_ns_rl_param_fct(yearly_var_adjust_stations_no_max_GMST, models, ns_models, confidence_interval='boot', sigma_exp=TRUE, rl=rl, return_period, return_horizon, variable)
ns_model_GMST_no_max.gev.params = ns_model_GMST_no_max.gev$params
ns_model_GMST_no_max.gev.rl = ns_model_GMST_no_max.gev$return_levels

ns_model_EMST_no_max.gev = gev_ns_rl_param_fct(yearly_var_adjust_stations_no_max_EMST, models, ns_models, confidence_interval='boot', sigma_exp=TRUE, rl=rl, return_period, return_horizon, variable)
ns_model_EMST_no_max.gev.params = ns_model_EMST_no_max.gev$params
ns_model_EMST_no_max.gev.rl = ns_model_EMST_no_max.gev$return_levels

saveRDS(ns_model_GMST_no_max.gev.params, file=paste(data_path, "gev_ns_no_max_",variable,"_adjust_",method,"_params_GMST.Rda", sep=""))
saveRDS(ns_model_GMST_no_max.gev.rl, file=paste(data_path, "gev_ns_no_max_",variable,"_adjust_",method,"_",return_period,"y_rl_",return_horizon[1],"_GMST.Rda", sep=""))

saveRDS(ns_model_EMST_no_max.gev.params, file=paste(data_path, "gev_ns_no_max_",variable,"_adjust_",method,"_params_EMST.Rda", sep=""))
saveRDS(ns_model_EMST_no_max.gev.rl, file=paste(data_path, "gev_ns_no_max_",variable,"_adjust_",method,"_",return_period,"y_rl_",return_horizon[1],"_EMST.Rda", sep=""))

# #Selection procedure + save
ns_model_GMST_no_max_final.gev = selection_models(daily_var_adjust_stations, yearly_var_adjust_stations_GMST, ns_model_GMST_no_max.gev.params, ns_rl = ns_model_GMST_no_max.gev.rl, variable)
ns_model_GMST_no_max_final.gev.params = ns_model_GMST_no_max_final.gev$params
ns_model_GMST_no_max_final.gev.rl = ns_model_GMST_no_max_final.gev$return_levels

ns_model_EMST_no_max_final.gev = selection_models(daily_var_adjust_stations, yearly_var_adjust_stations_EMST, ns_model_EMST_no_max.gev.params, ns_rl = ns_model_EMST_no_max.gev.rl, variable)
ns_model_EMST_no_max_final.gev.params = ns_model_EMST_no_max_final.gev$params
ns_model_EMST_no_max_final.gev.rl = ns_model_EMST_no_max_final.gev$return_levels

saveRDS(ns_model_GMST_no_max_final.gev.params, file=paste(data_path, "gev_final_ns_no_max_",variable,"_adjust_",method,"_params_GMST.Rda", sep=""))
saveRDS(ns_model_GMST_no_max_final.gev.rl, file=paste(data_path, "gev_final_ns_no_max_",variable,"_adjust_",method,"_",return_period,"y_rl_",return_horizon[1],"_GMST.Rda", sep=""))

saveRDS(ns_model_EMST_no_max_final.gev.params, file=paste(data_path, "gev_final_ns_no_max_",variable,"_adjust_",method,"_params_EMST.Rda", sep=""))
saveRDS(ns_model_EMST_no_max_final.gev.rl, file=paste(data_path, "gev_final_ns_no_max_",variable,"_adjust_",method,"_",return_period,"y_rl_",return_horizon[1],"_EMST.Rda", sep=""))

#
#GEV + save params in the case with regression without historical data
# ns_model_GMST_no_hist_no_max.gev = gev_ns_rl_param_fct(yearly_var_adjust_stations_no_max_GMST, models, ns_models, confidence_interval='boot', sigma_exp=TRUE, rl=rl, return_period, return_horizon, variable, with_historical=FALSE)
# ns_model_GMST_no_hist_no_max.gev.params = ns_model_GMST_no_hist_no_max.gev$params
# ns_model_GMST_no_hist_no_max.gev.rl = ns_model_GMST_no_hist_no_max.gev$return_levels
# 
# ns_model_EMST_no_hist_no_max.gev = gev_ns_rl_param_fct(yearly_var_adjust_stations_no_max_EMST, models, ns_models, confidence_interval='boot', sigma_exp=TRUE, rl=rl, return_period, return_horizon, variable, with_historical=FALSE)
# ns_model_EMST_no_hist_no_max.gev.params = ns_model_EMST_no_hist_no_max.gev$params
# ns_model_EMST_no_hist_no_max.gev.rl = ns_model_EMST_no_hist_no_max.gev$return_levels
# 
# saveRDS(ns_model_GMST_no_hist_no_max.gev.params, file=paste(data_path, "gev_ns_no_hist_no_max_",variable,"_adjust_",method,"_params_GMST.Rda", sep=""))
# saveRDS(ns_model_GMST_no_hist_no_max.gev.rl, file=paste(data_path, "gev_ns_no_hist_no_max_",variable,"_adjust_",method,"_",return_period,"y_rl_",return_horizon[1],"_GMST.Rda", sep=""))
# 
# saveRDS(ns_model_EMST_no_hist_no_max.gev.params, file=paste(data_path, "gev_ns_no_hist_no_max_",variable,"_adjust_",method,"_params_EMST.Rda", sep=""))
# saveRDS(ns_model_EMST_no_hist_no_max.gev.rl, file=paste(data_path, "gev_ns_no_hist_no_max_",variable,"_adjust_",method,"_",return_period,"y_rl_",return_horizon[1],"_EMST.Rda", sep=""))

# #Selection procedure + save with no hist models
# ns_model_GMST_no_hist_no_max_final.gev = selection_models(daily_var_adjust_stations, yearly_var_adjust_stations_GMST, ns_model_GMST_no_hist_no_max.gev.params, ns_rl = ns_model_GMST_no_hist_no_max.gev.rl, variable)
# ns_model_GMST_no_hist_no_max_final.gev.params = ns_model_GMST_no_hist_no_max_final.gev$params
# ns_model_GMST_no_hist_no_max_final.gev.rl = ns_model_GMST_no_hist_no_max_final.gev$return_levels
#
# ns_model_EMST_no_hist_no_max_final.gev = selection_models(daily_var_adjust_stations, yearly_var_adjust_stations_EMST, ns_model_EMST_no_hist_no_max.gev.params, ns_rl = ns_model_EMST_no_hist_no_max.gev.rl, variable)
# ns_model_EMST_no_hist_no_max_final.gev.params = ns_model_EMST_no_hist_no_max_final.gev$params
# ns_model_EMST_no_hist_no_max_final.gev.rl = ns_model_EMST_no_hist_no_max_final.gev$return_levels
#
# saveRDS(ns_model_GMST_no_hist_no_max_final.gev.params, file=paste(data_path, "gev_final_ns_no_hist_no_max_",variable,"_adjust_",method,"_params_GMST.Rda", sep=""))
# saveRDS(ns_model_GMST_no_hist_no_max_final.gev.rl, file=paste(data_path, "gev_final_ns_no_hist_no_max_",variable,"_adjust_",method,"_",return_period,"y_rl_",return_horizon[1],"_GMST.Rda", sep=""))
#
# saveRDS(ns_model_EMST_no_hist_no_max_final.gev.params, file=paste(data_path, "gev_final_ns_no_hist_no_max_",variable,"_adjust_",method,"_params_EMST.Rda", sep=""))
# saveRDS(ns_model_EMST_no_hist_no_max_final.gev.rl, file=paste(data_path, "gev_final_ns_no_hist_no_max_",variable,"_adjust_",method,"_",return_period,"y_rl_",return_horizon[1],"_EMST.Rda", sep=""))
#
#
# #
