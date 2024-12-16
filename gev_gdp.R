#Part I : Function structuring the data 

#The function extract data from a specific model and \
#return it as a global dataframe containing all stations and scenarios.
#WARNING only works for SAFRAN adjusted variable saved in the form :
#"tasmaxAdjust_France_IPSL_IPSL-CM6A-LR_historical_R2D2-L-NV-2L_SAFRAN_day_18500101-20141231_gridGE"
#Input : 
#   CMIP_groups : List of the groups to consider, ex : 'IPSL'
#   CMIP_models  : List of the models corresponding to the previous groups, ex : 'IPSL-CM6A-LR'
#   data_path : "/home/estimr3/yiou/ANDRA2022/CMIP6_R2D2/GrandEst/"
#   var : time of variable we want to extract ex: 'tas', 'tasmax', 'tasmin'
#   scenarios : List containing the name of the scenarios. Must be the form with which they are saved
#               ex : c('historical', 'ssp126', 'ssp245', 'ssp370', 'ssp585')
#   stations : List of stations name contained in the netcdf
#   season : 'JJA' for summer, 'DJF' for winter, 'all' if not specification
data_adjust_to_df = function(CMIP_groups, CMIP_models, data_path, var, scenarios, stations, method='R2D2', season = 'all'){
  df_all = data.frame(matrix(ncol = 6, nrow=0))
  for(id_mod in 1:length(CMIP_models)){
    CMIP_group_id = CMIP_groups[id_mod]
    CMIP_model_id = CMIP_models[id_mod]
    for(sc in scenarios){
      #Open the netcdf 
      if(sc=='historical'){
        period = '18500101-20141231'
        year_min = 1850
      }
      else{
        if(CMIP_model_id %in% c('TaiESM1','IPSL-CM6A-LR') & method=='R2D2'){
          period = '20150101-20991231' 
        }
        else{
          period = '20150101-21001231' 
        }
        year_min = 2015
      }
      
      if(method=='R2D2'){
        name = paste(data_path, var, "Adjust_France_", CMIP_group_id, '_', CMIP_model_id, '_', sc, '_R2D2-L-NV-2L_SAFRAN_day_', period, '_gridGE.nc', sep="")
      }
      else if(method=='CDFt'){
        name = paste(data_path, var, "Adjust_France_", CMIP_group_id, "_", CMIP_model_id, "_", sc, "_CDFt-L-1V-0L_day_",  period, '_gridGE.nc', sep="")
      }
      
      loop_nc = nc_open(name)
      
      #Extract time 
      time_loop = loop_nc$dim$time$vals
      time_unit = loop_nc$dim$time$unit
      
      if(time_unit %in% c("days since 1850-01-01", "days since 1850-1-1 00:00:00", "days since 1850-01-01 0:0:0.0", "days since 1850-01-01 00:00:00")){
        conv.time_loop = caldat(time_loop+julday(1,1,1850))
      }
      else if(time_unit %in% c("days since 0001-01-01", "days since 0001-01-01 00:00:00")){
        conv.time_loop = caldat(time_loop+julday(1,1,0001))
      }
      else if(time_unit %in% c("days since 2015-01-01", "days since 2015-01-01 00:00:00")){
        conv.time_loop = caldat(time_loop+julday(1,1,2015))
      }
      
      #Extract variable
      var_loop = ncvar_get(loop_nc, paste(var, 'Adjust', sep=''))
      if(var %in% c('tas', 'tasmax', 'tasmin')){
        var_loop = var_loop -273
      }
      
      #Select specific season
      if(season=='JJA'){
        var_loop = var_loop[,conv.time_loop$month %in% c(6:8)]
        years_loop = conv.time_loop$year[conv.time_loop$month %in% c(6:8)]
        months_loop = conv.time_loop$month[conv.time_loop$month %in% c(6:8)]
        days_loop = conv.time_loop$day[conv.time_loop$month %in% c(6:8)]
        time_loop = julday(days_loop, months_loop, years_loop)
      }
      else if(season=='DJF'){
        var_loop = var_loop[,conv.time_loop$month %in% c(1,2,12)]
        years_loop = conv.time_loop$year[conv.time_loop$month %in% c(1,2,12)]
        months_loop = conv.time_loop$month[conv.time_loop$month %in% c(1,2,12)]
        days_loop = conv.time_loop$day[conv.time_loop$month %in% c(1,2,12)]
        time_loop = julday(days_loop, months_loop, years_loop)
      }
      else{
        years_loop = conv.time_loop$year
        months_loop = conv.time_loop$month
        days_loop = conv.time_loop$day
        time_loop = julday(days_loop, months_loop, years_loop)
      }
      
      #Save data to global netcdf
      rownames(var_loop) = stations
      n_stations = length(stations)
      df_loop = as.data.frame(t(var_loop)) %>%
        gather(station, var, stations) %>%
        mutate(scenario = sc, year = rep(years_loop, n_stations), jul_time = rep(time_loop, n_stations), CMIP_model=CMIP_model_id)
      df_all = rbind(df_all, df_loop %>% filter(year >= year_min))
    }
  }
  return(df_all)
}

#Same function as above but for regular CMIP6 data 
#
#Input : 
#   CMIP_models : List of models ex = 'IPSL-CM6A-LR'
#   data_path : "/home/estimr3/lhasbini/CMIP6/"
#   var_id : Variable name, ex : 'tasmax'
#   scenarios : List of scenarios
data_CMIP6_fldmean_yearmean_to_df = function(CMIP_models, data_path, var_id, scenarios, region){
  df_all = data.frame(matrix(ncol = 4, nrow = 0))
  colnames(df_all) = c('var', 'year', 'scenario', 'CMIP_model')
  for(CMIP_model_id in CMIP_models){
    for(sc in scenarios){
      #df_loop = data.frame(matrix(ncol = 3, nrow = 0))
      #colnames(df_loop) = colnames(df_all)
      #Selection the correct reference period
      if(sc=='historical'){
        perref = '18500101-20141231'
        start_jul = julday(1,1,1850)
      }
      else{
        perref = '20150101-21001231'
        start_jul = julday(1,1,2015)
      }
      name = paste(data_path, var_id, '_fldmean_yearmean_', region, '_', CMIP_model_id, '_', sc, '_', perref, '.nc', sep='')
      nc_loop = nc_open(name)
      
      #Fix time
      time_loop = time_unit = nc_loop$dim$time$vals
      time_unit = nc_loop$dim$time$unit
      
      if(time_unit=="days since 1850-01-01" | 
         time_unit=="days since 1850-01-01 0:0:0.0" | 
         time_unit=="days since 1850-01-01 00:00:00" |
         time_unit=="days since 1850-1-1" |
         time_unit=="days since 1850-1-1 00:00:00"){
        conv.time_loop = caldat(time_loop+julday(1,1,1850))
      }
      else if(time_unit=="days since 0001-01-01"| time_unit=="days since 0001-01-01 00:00:00"){
        conv.time_loop = caldat(time_loop+julday(1,1,0001))
      }
      else if(time_unit=="days since 2015-01-01" | time_unit=="days since 2015-01-01 00:00:00"){
        conv.time_loop = caldat(time_loop+julday(1,1,2015))
      }
      
      #Extract variable
      var = ncvar_get(nc_loop, var_id)
      if(var_id %in% c('tas', 'tasmax', 'tasmin')){
        var = var -273
      }
      
      #Extract necessary variables 
      #df_loop$var = var_loop
      df_loop = as.data.frame(var)
      df_loop$year = conv.time_loop$year
      df_loop$scenario = sc
      df_loop$CMIP_model = CMIP_model_id
      
      #Save to global dataframe
      df_all = rbind(df_all, df_loop)
    } 
  }
  return(df_all)
}

#Extract variable over the GrandEst grid point in CMIP6 data
#Input : 
#   CMIP_models : List of models ex = 'IPSL-CM6A-LR'
#   data_path : "/home/estimr3/lhasbini/CMIP6/"
#   var_id : Variable name, ex : 'tasmax'
#   scenarios : List of scenarios
#   season : 'JJA' for summer, 'DJF' for winter, 'all' if not specification
data_CMIP6_GE_to_df = function(CMIP_models, data_path, var_id, scenarios, season = 'all'){
  df_all = data.frame(matrix(ncol = 6, nrow = 0))
  colnames(df_all) = c('var', 'jul_time', 'year', 'scenario', 'CMIP_model', 'station')
  for(CMIP_model_id in CMIP_models){
    for(sc in scenarios){
      #df_loop = data.frame(matrix(ncol = 3, nrow = 0))
      #colnames(df_loop) = colnames(df_all)
      #Selection the correct reference period
      if(sc=='historical'){
        perref = '18500101-20141231'
        start_jul = julday(1,1,1850)
        year_min = 1850
      }
      else{
        perref = '20150101-21001231'
        start_jul = julday(1,1,2015)
        year_min = 2015
      }
      name = paste(data_path, var_id, '_GrandEst_', CMIP_model_id, '_', sc, '_', perref, '.nc', sep='')
      nc_loop = nc_open(name)
      
      #Fix time
      time_loop = nc_loop$dim$time$vals
      time_unit = nc_loop$dim$time$unit
      
      if(time_unit=="days since 1850-01-01" | 
         time_unit=="days since 1850-01-01 0:0:0.0" | 
         time_unit=="days since 1850-01-01 00:00:00" |
         time_unit=="days since 1850-1-1" |
         time_unit=="days since 1850-1-1 00:00:00"){
        conv.time_loop = caldat(time_loop+julday(1,1,1850))
      }
      else if(time_unit=="days since 0001-01-01"| time_unit=="days since 0001-01-01 00:00:00"){
        conv.time_loop = caldat(time_loop+julday(1,1,0001))
      }
      else if(time_unit=="days since 2015-01-01" | time_unit=="days since 2015-01-01 00:00:00"){
        conv.time_loop = caldat(time_loop+julday(1,1,2015))
      }
      
      #Extract variable
      var = ncvar_get(nc_loop, var_id)
      if(var_id %in% c('tas', 'tasmax', 'tasmin')){
        var = var -273
      }
      
      #Select specific season
      if(season=='JJA'){
        var = var[conv.time_loop$month %in% c(6:8)]
        years_loop = conv.time_loop$year[conv.time_loop$month %in% c(6:8)]
        months_loop = conv.time_loop$month[conv.time_loop$month %in% c(6:8)]
        days_loop = conv.time_loop$day[conv.time_loop$month %in% c(6:8)]
        time_loop = julday(days_loop, months_loop, years_loop)
      }
      else if(season=='DJF'){
        var = var[conv.time_loop$month %in% c(1,2,12)]
        years_loop = conv.time_loop$year[conv.time_loop$month %in% c(1,2,12)]
        months_loop = conv.time_loop$month[conv.time_loop$month %in% c(1,2,12)]
        days_loop = conv.time_loop$day[conv.time_loop$month %in% c(1,2,12)]
        time_loop = julday(days_loop, months_loop, years_loop)
      }
      else{
        years_loop = conv.time_loop$year
        months_loop = conv.time_loop$month
        days_loop = conv.time_loop$day
        time_loop = julday(days_loop, months_loop, years_loop)
      }
      
      #Extract necessary variables 
      #df_loop$var = var_loop
      df_loop = as.data.frame(var)
      df_loop$jul_time = time_loop
      df_loop$year = years_loop
      df_loop$scenario = sc
      df_loop$CMIP_model = CMIP_model_id
      df_loop$station = 'GrandEst'
      
      #Save to global dataframe
      df_all = rbind(df_all, df_loop %>% filter(year >= year_min))
    } 
  }
  return(df_all)
}

#The next function will take input the dataframe created with the previous function and return the n year yearly max 
# var_ndays correspond the the yearly maximum of the n days running mean
#Input
#   df : The dataframe from the previous function
#   n : Numer of days on which the running mean must be done 
data_to_varn_max = function(df, n=1){
  TXn = data.frame(matrix(ncol = 5, nrow = 0))
  colnames(TXn) = c('CMIP_model','station', 'var', 'scenario', 'year')
  
  CMIP_models = unique(df$CMIP_model)
  scenarios = unique(df$scenario)
  stations = unique(df$station)
  
  for (mod in CMIP_models){
    for(sc in scenarios){
      for(st in stations){
        df_loop = df %>% filter(scenario == sc, station==st, CMIP_model==mod)
        years_loop = unique(df_loop$year)
        TXn_loop = data.frame(matrix(nrow=length(years_loop),ncol=1))
        
        for(iy in 1:length(years_loop)){
          y = years_loop[iy]
          data_loop_year = df_loop %>% filter(year == y)
          TXn_loop[iy, 1] = max(rollmean(data_loop_year$var, n))
        }
      
        df_loop.max = data.frame(matrix(ncol = 5, nrow = length(years_loop)))
        x <- c("station", "var", "scenario", "year", "CMIP_model")
        
        colnames(df_loop.max) <- x
        df_loop.max$var = TXn_loop[,1]
        df_loop.max$station = st
        df_loop.max$scenario = sc
        df_loop.max$year = years_loop
        df_loop.max$CMIP_model = mod
        TXn = rbind(TXn, df_loop.max)
      }
    } 
  }
  return(TXn)
}

#The next function works similarly to the previous function but return the n year yearly min
# var_ndays correspond the the yearly minimal of the n days running mean
#Input
#   df : The dataframe from the previous function
#   n : Numer of days on which the running mean must be done 
data_to_varn_min = function(df, n=1){
  TXn = data.frame(matrix(ncol = 5, nrow = 0))
  colnames(TXn) = c('CMIP_model','station', 'var', 'scenario', 'year')
  
  CMIP_models = unique(df$CMIP_model)
  scenarios = unique(df$scenario)
  stations = unique(df$station)
  
  for (mod in CMIP_models){
    for(sc in scenarios){
      for(st in stations){
        df_loop = df %>% filter(scenario == sc, station==st, CMIP_model==mod)
        years_loop = unique(df_loop$year)
        TXn_loop = data.frame(matrix(nrow=length(years_loop),ncol=1))
        
        for(iy in 1:length(years_loop)){
          y = years_loop[iy]
          data_loop_year = df_loop %>% filter(year == y)
          TXn_loop[iy, 1] = min(rollmean(data_loop_year$var, n))
        }
        
        df_loop.min = data.frame(matrix(ncol = 5, nrow = length(years_loop)))
        x <- c("station", "var", "scenario", "year", "CMIP_model")
        
        colnames(df_loop.min) <- x
        df_loop.min$var = TXn_loop[,1]
        df_loop.min$station = st
        df_loop.min$scenario = sc
        df_loop.min$year = years_loop
        df_loop.min$CMIP_model = mod
        TXn = rbind(TXn, df_loop.min)
      }
    } 
  }
  return(TXn)
}

#Same as above but for CMIP data without Station Specification
#Input
#   df : The dataframe from the previous function
#   n : Numer of days on which the running mean must be done 
data_CMIP_to_varn_max = function(df, n=1){
  TXn = data.frame(matrix(ncol = 4, nrow = 0))
  colnames(TXn) = c('CMIP_model', 'var', 'scenario', 'year')
  
  CMIP_models = unique(df$CMIP_model)
  scenarios = unique(df$scenario)
  stations = unique(df$station)
  
  for (mod in CMIP_models){
    for(sc in scenarios){
        df_loop = df %>% filter(scenario == sc, CMIP_model==mod)
        years_loop = unique(df_loop$year)
        TXn_loop = data.frame(matrix(nrow=length(years_loop),ncol=1))
        
        for(iy in 1:length(years_loop)){
          y = years_loop[iy]
          data_loop_year = df_loop %>% filter(year == y)
          TXn_loop[iy, 1] = max(rollmean(data_loop_year$var, n))
        }
        
        df_loop.max = data.frame(matrix(ncol = 4, nrow = length(years_loop)))
        x <- c("var", "scenario", "year", "CMIP_model")
        
        colnames(df_loop.max) <- x
        df_loop.max$var = TXn_loop[,1]
        df_loop.max$scenario = sc
        df_loop.max$year = years_loop
        df_loop.max$CMIP_model = mod
        TXn = rbind(TXn, df_loop.max)
    } 
  }
  return(TXn)
}

#PART II : Functions computing parameters

#The function below compute the spline according to a specific threshold and to the variable under 'covariate'
#Input :
#   df : Datafram with all data (station, scenario, var, year, CMIP_model)
#   percentile_threshold : number (ex 0.95)
#   spar = 0.9
spline_threshold_from_year = function(df_input, percentile_threshold, spar=0.9){
  #Initialize parameters
  stations = unique(df_input$station)
  scenarios = unique(df_input$scenario)
  CMIP_models = unique(df_input$CMIP_model)
  
  #Create dataframes
  df_u_t = data.frame(matrix(ncol = ncol(df_input)+1, nrow = 0))
  colnames(df_u_t) = append(colnames(df_input),c('u'))
  
  for(CMIP_mod in CMIP_models){
    for(sc in scenarios){
      for(st in stations){
        data_loop.quantile = df_input %>% filter(station == st, scenario == sc, CMIP_model==CMIP_mod)
        time_loop = data_loop.quantile$year
        
        Y = tapply(data_loop.quantile$var, time_loop, quantile, probs=percentile_threshold)
        CMIP_model_spline_loop = smooth.spline(x=unique(time_loop), y=Y, spar=spar)
        data_loop.quantile$u  = predict(CMIP_model_spline_loop, time_loop)$y
        df_u_t = rbind(df_u_t, data_loop.quantile)
      }
    }
  }
  return(df_u_t)
}

#Same as before but spline done regarding to the variable 'covariate'
spline_threshold_from_covariate = function(df_input, percentile_threshold, spar=0.9){
  #Initialize parameters
  stations = unique(df_input$station)
  scenarios = unique(df_input$scenario)
  CMIP_models = unique(df_input$CMIP_model)
  
  #Create dataframes
  df_u_t = data.frame(matrix(ncol = ncol(df_input)+1, nrow = 0))
  colnames(df_u_t) = append(colnames(df_input),c('u'))
  
  for(CMIP_mod in CMIP_models){
    for(sc in scenarios){
      for(st in stations){
        data_loop.quantile = df_input %>% filter(station == st, scenario == sc, CMIP_model==CMIP_mod)
        
        cov_loop = data_loop.quantile$covariate
        C = tapply(data_loop.quantile$var, cov_loop, quantile, probs=percentile_threshold)
        model_spline_loop = smooth.spline(x=unique(cov_loop), y=C, spar=spar)
        data_loop.quantile$u  = predict(model_spline_loop, cov_loop)$y
        df_u_t = rbind(df_u_t, data_loop.quantile)
      }
    } 
  }
  return(df_u_t)
}

#The following function will compute the normalized year add it to the column 'covariate'
year_norm_to_covarite = function(df_input){
  #Initialize parameters
  stations = unique(df_input$station)
  scenarios = unique(df_input$scenario)
  CMIP_models = unique(df_input$CMIP_model)
  
  #Create dataframes
  df_year_norm_cov = data.frame(matrix(ncol = ncol(df_input)+1, nrow = 0))
  colnames(df_year_norm_cov) <- append(colnames(df_input),c('covariate'))
  
  for(CMIP_mod in CMIP_models){
    for(sc in scenarios){
      for(st in stations){
        data_loop = df_input %>% filter(station == st, scenario == sc, CMIP_model==CMIP_mod)
        time_loop = data_loop$year
        data_loop$covariate = (time_loop-mean(time_loop))/var(time_loop)
        df_year_norm_cov = rbind(df_year_norm_cov, data_loop)
      }
    }
  }
  return(df_year_norm_cov)
}

#PART III : Functions performing gev/gpd

#The function return the value of the return level as well as the confidence interval(is specified)
#Input : 
#   location : Mu, location parameter
#   scale : Sigma, scale parameter
#   shape : Xhi, shape parameter
#   return period : List of return peiod (in year)
gev_theoretical_rl = function(location, scale, shape, return_period){
  lim = abs(shape)
  return_levels = data.frame(matrix(nrow=length(return_period),ncol=1))
  #First case if shape==0
  if(lim>0.0001){
    for(i in 1:length(return_period)){
      m = return_period[i]
      return_levels[i,1] = location - scale/shape*(1-(-log(1-1/m))^(-shape))
    }
  }
  #Second case if shape !=0
  else{
    for(i in 1:length(return_period)){
      m = return_period[i]
      return_levels[i,1] = location - scale*log(-log(1-1/m))
    }
  }
  return(return_levels[,1])
}

## The following function compute return level with gev method for given characteristics
#Input : 
#   df_max : list of maximums (yearly). Can be a dataframe
#   rl : Boolean. If TRUE, compute the return levels 
#   return_period : list
#   theoretical : Boolean. If TRUE, add a column for the theoretical return levels 
#   scearios : list
#   stations : list

gev_rl_fct = function(df_max, rl=TRUE, return_period, theoretical=FALSE, scenarios, stations, CMIP_models){
  if(rl){
    #Dataframes for return levels
    df.gev.rl = data.frame(matrix(ncol = 7, nrow = 0))
    colnames(df.gev.rl) = c("return_period", "low", "return_level", "high", "scenario", "station", "CMIP_model")
  }
  
  #Dataframes for parameters
  params.gev = data.frame(matrix(ncol=6, nrow=0))
  x_param = c("loc", "scale", "shape", "scenario", "station", "CMIP_model")
  colnames(params.gev) = x_param
  
  #If yes, add a columns for the computation of theoretical return levels 
  if(theoretical){
    df_theo = data.frame(matrix(ncol=1, nrow=0))
    colnames(df_theo) = c('return_level_theo')
    df.gev.rl = cbind(df.gev.rl, df_theo)
  }
  
  for(CMIP_mod in CMIP_models){
    for(sc in scenarios){
      for(st in stations){
        #Selection of the data
        data_loop.gev = df_max %>% 
          filter(station==st, scenario==sc, CMIP_model==CMIP_mod)
        
        #Computation of the gev 
        df_var_loop.gev = fevd(data_loop.gev$var)
        
        #Save parameters 
        look.gev = summary(df_var_loop.gev, silent = TRUE)
        loc_loop.gev = look.gev$par[1]
        scale_loop.gev = look.gev$par[2]
        shape_loop.gev = look.gev$par[3]
        params.gev = rbind(params.gev, data.frame(loc=loc_loop.gev, scale = scale_loop.gev, shape = shape_loop.gev, scenario=sc, station=st, CMIP_model=CMIP_mod))
        
        #Computation of return levels
        if(rl){
          df_var_loop.gev.rl = return.level(df_var_loop.gev, return.period=return_period, do.ci=TRUE)
          
          #Save the return levels and confidence interval in a separate df
          df_loop.gev.rl = data.frame(matrix(ncol = ncol(df.gev.rl), nrow = length(return_period)))
          colnames(df_loop.gev.rl) = colnames(df.gev.rl)
          
          df_loop.gev.rl$return_period = return_period
          df_loop.gev.rl$station = st
          df_loop.gev.rl$scenario = sc
          df_loop.gev.rl$CMIP_model = CMIP_mod
          df_loop.gev.rl$low= df_var_loop.gev.rl[,1]
          df_loop.gev.rl$return_level = df_var_loop.gev.rl[,2]
          df_loop.gev.rl$high = df_var_loop.gev.rl[,3]
        }
        
        #Computation of theoretical return levels  
        if(theoretical){
          df_loop.gev.rl$return_level_theo = gev_theoretical_rl(loc_loop.gev, scale_loop.gev, shape_loop.gev, return_period)
        }
        #Save to the previous data
        if(rl){
          df.gev.rl = rbind(df.gev.rl, df_loop.gev.rl)
        }
      }
    }
  }

  if(rl){
    return(list(params = params.gev, return_levels = df.gev.rl))
  }
  else{
    return(params.gev)
  }
}

#The following function compute parameters and return levels with gev method in the non-stationary case, as a function of a given covariate name
#5 models will be created : 
#   M(000) : Stationary
#   M(100) : loc varying lineary as a fct of the covariate
#   M(110) : Loc varying lineary and scale varying exponentially as a fct of the covariate 
#   M(-100) : Linear loc removed before, stationary GEV
#   M(-110) : Linear loc removed before, GEV exponential in scale 
#Input :
#   df_max_with_cov : list of maximums (yearly). Assume that the data frame already contains the covariate in the column 'covariate'
#   rl : Boolean. If TRUE, compute the return levels 
#   return_period : list
#   stations : list
#   
gev_ns_rl_param_fct = function(df_max_with_cov, CMIP_models, ns_model, confidence_interval = 'model', sigma_exp=TRUE, rl=TRUE, return_period=100, return_horizon=c(2100), variable='tasmax', with_historical=TRUE){
  if(rl){
    #Dataframe of return level of ns models
    ns_models.gev.rl = data.frame(matrix(ncol = 8, nrow = 0))
    colnames(ns_models.gev.rl) = c('station', 'scenario', 'year', 'CMIP_model', 'start_cov_value', 'return_period', 'return_level', 'ns_model')
  }
  
  #Dataframe of the nc models
  ns_models.gev.params = data.frame(matrix(ncol = 25, nrow = 0))
  colnames(ns_models.gev.params) = c('station', 'scenario', 'CMIP_model', 'ns_model',  
                                     'loc_0', 'loc_0_low', 'loc_0_high', 'loc_1', 'loc_1_low', 'loc_1_high', 
                                     'scale_0', 'scale_0_low', 'scale_0_high', 'scale_1', 'scale_1_low', 'scale_1_high', 
                                     'shape', 'shape_low', 'shape_high', 'nllh',
                                     'a', 'a_std', 'b', 'b_std', 'r2')
  
  stations = unique(df_max_with_cov$station)
  #stations = c('Cirfontaines', 'Bure', 'St Dizier')
  
  for(st in stations){
    print(st)
    for(CMIP_mod in CMIP_models){
      print(CMIP_mod)
      scenarios = unique(df_max_with_cov %>% filter(station==st,CMIP_model==CMIP_mod) %>% select(scenario))
      scenarios = scenarios[scenarios!= 'historical']
      for(sc in scenarios){
        print(sc)
        if(with_historical){
          data_loop = df_max_with_cov %>% filter(station==st, scenario%in%c('historical', sc), CMIP_model==CMIP_mod)
        }
        else{
          data_loop = df_max_with_cov %>% filter(station==st, scenario%in%c(sc), CMIP_model==CMIP_mod)
        }
        data_loop_gev = data.frame(data_loop)
        data_loop_gev = na.omit(data_loop_gev)
        #data_loop_gev = data_loop[complete.cases(data_loop), ]
        
        if(variable == 'tasmin'){
          data_loop_gev = data_loop_gev %>% mutate(var = -var)
          data_loop = data_loop %>% mutate(var = -var)
        } 
        
        lm_mod = lm(var~covariate, data=data_loop_gev)
        lm_mod_sum = summary(lm_mod, silent=TRUE)
        var_norm = data_loop_gev$var - predict(lm_mod)+rnorm(nrow(data_loop_gev), 0, 0.01)
        r2_loop = lm_mod_sum$r.squared
        
        #Add somme jitter to make sure that there is not 2 data points with the exact same value
        data_loop_gev = data_loop_gev %>% mutate(var = var + rnorm(nrow(data_loop_gev), 0, 0.01))
        
        for(ns_mod in ns_model){
          #print(ns_mod)
          df_ns_model =  data.frame(station = st, 
                                    scenario = sc, 
                                    CMIP_model = CMIP_mod, 
                                    ns_model= ns_mod, 
                                    #Location
                                    loc_0 = 0,
                                    loc_0_low = 0,
                                    loc_0_high = 0,
                                    loc_1 = 0,
                                    loc_1_low = 0, 
                                    loc_1_high = 0, 
                                    #Scale
                                    scale_0 = 0,
                                    scale_0_low = 0, 
                                    scale_0_high = 0, 
                                    scale_1 = 0,
                                    scale_1_low = 0, 
                                    scale_1_high = 0,
                                    #Shape
                                    shape = 0,
                                    shape_low = 0, 
                                    shape_high = 0,
                                    #Log Likelihood
                                    nllh = 0, 
                                    #Coef of the regression
                                    a=0,
                                    a_std = 0,
                                    b=0, 
                                    b_std = 0, 
                                    r2=0)
          #### SAVE PARAMETER
          
          df_ns_model$r2 = as.numeric(r2_loop)
          
          if(ns_mod=='M-100'){
            model = fevd(var_norm)
            
            df_ns_model$loc_0 = as.numeric(model$results$par[1])
            df_ns_model$scale_0 = as.numeric(model$results$par[2])
            df_ns_model$shape = as.numeric(model$results$par[3])
            df_ns_model$nllh = as.numeric(model$results$value)
            df_ns_model$a = lm_mod$coefficients[2]
            df_ns_model$a_std = lm_mod_sum$coefficients[4]
            df_ns_model$b = lm_mod$coefficients[1]
            df_ns_model$b_std = lm_mod_sum$coefficients[3]
          }
          else if(ns_mod=='M-110'){
            cov = data_loop_gev$covariate
            model = fevd(var_norm, scale.fun = ~cov, use.phi=sigma_exp)

            df_ns_model$loc_0 = as.numeric(model$results$par[1])
            df_ns_model$scale_0 = as.numeric(model$results$par[2])
            df_ns_model$scale_1 = as.numeric(model$results$par[3])
            df_ns_model$shape = as.numeric(model$results$par[4])
            df_ns_model$nllh = as.numeric(model$results$value)
            df_ns_model$a = lm_mod$coefficients[2]
            df_ns_model$a_std = lm_mod_sum$coefficients[4]
            df_ns_model$b = lm_mod$coefficients[1]
            df_ns_model$b_std = lm_mod_sum$coefficients[3]
          }
          else if(ns_mod == 'M000'){
            model = fevd(var, data=data_loop_gev)
            df_ns_model$loc_0 = as.numeric(model$results$par[1])
            df_ns_model$scale_0 = as.numeric(model$results$par[2])
            df_ns_model$shape = as.numeric(model$results$par[3])
            df_ns_model$nllh = as.numeric(model$results$value)
          }
          else if(ns_mod=='M100'){
            model = fevd(var, data=data_loop_gev, location.fun = ~covariate)
            df_ns_model$loc_0 = as.numeric(model$results$par[1])
            df_ns_model$loc_1 = as.numeric(model$results$par[2])
            df_ns_model$scale_0 = as.numeric(model$results$par[3])
            df_ns_model$shape = as.numeric(model$results$par[4])
            df_ns_model$nllh = as.numeric(model$results$value)
          }
          else if(ns_mod=='M110'){
            model = fevd(var, data=data_loop_gev, location.fun = ~covariate, scale.fun = ~covariate, use.phi = sigma_exp)
            df_ns_model$loc_0 = as.numeric(model$results$par[1])
            df_ns_model$loc_1 = as.numeric(model$results$par[2])
            df_ns_model$scale_0 = as.numeric(model$results$par[3])
            df_ns_model$scale_1 = as.numeric(model$results$par[4])
            df_ns_model$shape = as.numeric(model$results$par[5])
            df_ns_model$nllh = as.numeric(model$results$value)
          }
          else if(ns_mod=='M-1-10'){
            loess_mean = loess(var~year, span=0.17, data=data_loop_gev)
            mean_X = predict(loess_mean)
            var_X = (data_loop_gev$var - mean_X)^2
            loess_var = loess(var_X ~ data_loop_gev$year, span=0.17)
            s_X = predict(loess_var)
            Y = (data_loop$var - mean_X)/s_X
            data_loop$m_X = predict(loess_mean, newdata = data_loop$year)
            data_loop$s_X = predict(loess_var, newdata = data_loop$year)
            
            model = fevd(Y)
            df_ns_model$loc_0 = as.numeric(model$results$par[1])
            df_ns_model$scale_0 = as.numeric(model$results$par[2])
            df_ns_model$shape = as.numeric(model$results$par[3])
            df_ns_model$nllh = as.numeric(model$results$value)
          }


          ####COMPUTE AND SAVE CONFIDENCE INTERVAL
          #summary(model)
          if(df_ns_model$shape <= 0 & is.null(parcov.fevd(model))==FALSE){
          if(confidence_interval=='model'){
            ci_model = summary(model, silent=TRUE)$se.theta
            if(ns_mod %in% c('M000', 'M-100', 'M-1-10')){
              df_ns_model$loc_0_low = model$results$par[1]-ci_model[1]
              df_ns_model$loc_0_high = model$results$par[1]+ci_model[1]
              df_ns_model$scale_0_low = model$results$par[2]-ci_model[2] 
              df_ns_model$scale_0_high = model$results$par[2]+ci_model[2] 
              df_ns_model$shape_low = model$results$par[3]-ci_model[3] 
              df_ns_model$shape_high = model$results$par[3]+ci_model[3]
            }
            else if(ns_mod=='M-110'){
              df_ns_model$loc_0_low = model$results$par[1]-ci_model[1]
              df_ns_model$loc_0_high = model$results$par[1]+ci_model[1]
              df_ns_model$scale_0_low = model$results$par[2]-ci_model[2]
              df_ns_model$scale_0_high = model$results$par[2]+ci_model[2]
              df_ns_model$scale_1_low = model$results$par[3]-ci_model[3] 
              df_ns_model$scale_1_high = model$results$par[3]+ci_model[3] 
              df_ns_model$shape_low = model$results$par[4]-ci_model[4] 
              df_ns_model$shape_high = model$results$par[4]+ci_model[4]
            }
            else if(ns_mod=='M100'){
              df_ns_model$loc_0_low = model$results$par[1]-ci_model[1]
              df_ns_model$loc_0_high = model$results$par[1]+ci_model[1]
              df_ns_model$loc_1_low = model$results$par[2]-ci_model[2]
              df_ns_model$loc_1_high = model$results$par[2]+ci_model[2]
              df_ns_model$scale_0_low = model$results$par[3]-ci_model[3] 
              df_ns_model$scale_0_high = model$results$par[3]+ci_model[3] 
              df_ns_model$shape_low = model$results$par[4]-ci_model[4] 
              df_ns_model$shape_high = model$results$par[4]+ci_model[4]
            }
            else if(ns_mod=='M110'){
              df_ns_model$loc_0_low = model$results$par[1]-ci_model[1]
              df_ns_model$loc_0_high = model$results$par[1]+ci_model[1]
              df_ns_model$loc_1_low = model$results$par[2]-ci_model[2]
              df_ns_model$loc_1_high = model$results$par[2]+ci_model[2]
              df_ns_model$scale_0_low = model$results$par[3]-ci_model[3] 
              df_ns_model$scale_0_high = model$results$par[3]+ci_model[3] 
              df_ns_model$scale_1_low = model$results$par[4]-ci_model[4] 
              df_ns_model$scale_1_high = model$results$par[4]+ci_model[4] 
              df_ns_model$shape_low = model$results$par[5]-ci_model[5] 
              df_ns_model$shape_high = model$results$par[5]+ci_model[5]
            }
          }
          else if(confidence_interval=='boot'){
            ci_model = ci(model, type = 'parameter', method="boot")
            if(ns_mod %in% c('M000', 'M-100', 'M-1-10')){
              df_ns_model$loc_0_low = ci_model[1,1]
              df_ns_model$loc_0_high = ci_model[1,3]
              df_ns_model$scale_0_low = ci_model[2,1]
              df_ns_model$scale_0_high = ci_model[2,3]
              df_ns_model$shape_low = ci_model[3,1] 
              df_ns_model$shape_high = ci_model[3,3]
            }
            else if(ns_mod=='M-110'){
              df_ns_model$loc_0_low = ci_model[1,1]
              df_ns_model$loc_0_high = ci_model[1,3]
              df_ns_model$scale_0_low = ci_model[2,1]
              df_ns_model$scale_0_high = ci_model[2,3]
              df_ns_model$scale_1_low = ci_model[3,1]
              df_ns_model$scale_1_high = ci_model[3,3]
              df_ns_model$shape_low = ci_model[4,1] 
              df_ns_model$shape_high = ci_model[4,3]
            }
            else if(ns_mod=='M100'){
              df_ns_model$loc_0_low = ci_model[1,1]
              df_ns_model$loc_0_high = ci_model[1,3]
              df_ns_model$loc_1_low = ci_model[2,1]
              df_ns_model$loc_1_high = ci_model[2,3]
              df_ns_model$scale_0_low = ci_model[3,1]
              df_ns_model$scale_0_high = ci_model[3,3]
              df_ns_model$shape_low = ci_model[4,1] 
              df_ns_model$shape_high = ci_model[4,3]
            }
            else if(ns_mod=='M110'){
              df_ns_model$loc_0_low = ci_model[1,1]
              df_ns_model$loc_0_high = ci_model[1,3]
              df_ns_model$loc_1_low = ci_model[2,1]
              df_ns_model$loc_1_high = ci_model[2,3]
              df_ns_model$scale_0_low = ci_model[3,1]
              df_ns_model$scale_0_high = ci_model[3,3]
              df_ns_model$scale_1_low = ci_model[4,1]
              df_ns_model$scale_1_high = ci_model[4,3]
              df_ns_model$shape_low = ci_model[5,1] 
              df_ns_model$shape_high = ci_model[5,3]
            }
          }
          
          # if(variable == 'tasmin'){
          #   df_ns_model = df_ns_model %>% mutate(loc_0 = -loc_0, 
          #                                        loc_0_low = -loc_0_low , 
          #                                        loc_0_high = -loc_0_high , 
          #                                        loc_1 = -loc_1 , 
          #                                        loc_1_low = -loc_1_low, 
          #                                        loc_1_high = -loc_1_high, 
          #                                        a = -a, 
          #                                        b = -b)
          # }
          #print(df_ns_model)
          ns_models.gev.params = rbind(ns_models.gev.params, df_ns_model)
            if(rl){
              #Compute return levels only in the case when convergence is obtained
              if(return_horizon > max(data_loop$year)){
                return_horizon_loop = max(data_loop$year)
              }
              else{
                return_horizon_loop = return_horizon
              }
              cov_horizon = data_loop %>% filter(year %in% return_horizon_loop) %>% select(covariate)
              adjust_rl = 0
              
              if(ns_mod=='M000'){
                v1 = make.qcov(model)
              }
              else if(ns_mod == 'M100'){
                v1 = make.qcov(model, vals = list(mu1 = as.numeric(cov_horizon)))
              }
              else if(ns_mod == 'M110'){
                v1 = make.qcov(model, vals = list(mu1 = as.numeric(cov_horizon), phi1 = as.numeric(cov_horizon)))
              }
              else if(ns_mod == 'M-100'){
                v1 = make.qcov(model)
                adjust_rl = df_ns_model$a * cov_horizon + df_ns_model$b
              }
              else if(ns_mod=='M-110'){
                v1 = make.qcov(model, vals = list(phi1 = as.numeric(cov_horizon)))
                adjust_rl = df_ns_model$a * cov_horizon + df_ns_model$b
              }
              else if(ns_mod=='M-1-10'){
                v1 = make.qcov(model)
                m_loop = data_loop %>% filter(year %in% return_horizon_loop) %>% select(m_X)
                s_loop = data_loop %>% filter(year %in% return_horizon_loop) %>% select(s_X)
              }
              

              rl_loop = return.level(model, return.period = return_period, qcov = v1, do.ci=TRUE)

              #Compute and save return levels 
              col_names_rl = c('CMIP_model', 'scenario', 'station', 'ns_model', 'return_horizon', 'return_period', 'return_level', 'return_level_low', 'return_level_high')
              model_loop.rl = data.frame(matrix(ncol = length(col_names_rl), nrow = length(return_horizon)))
              colnames(model_loop.rl) = col_names_rl
              if(ns_mod=='M-1-10'){
                model_loop.rl$return_level = as.numeric(rl_loop[2]*s_loop + m_loop)
                model_loop.rl$return_level_low = as.numeric(rl_loop[1]*s_loop + m_loop)
                model_loop.rl$return_level_high = as.numeric(rl_loop[3]*s_loop + m_loop)
              }
              else {
                if(variable=='tasmax'){
                  model_loop.rl$return_level = as.numeric(rl_loop[2] +adjust_rl)
                  model_loop.rl$return_level_low = as.numeric(rl_loop[1] +adjust_rl)
                  model_loop.rl$return_level_high = as.numeric(rl_loop[3] +adjust_rl)}
                else if(variable=='tasmin'){
                  model_loop.rl$return_level = as.numeric(-rl_loop[2] -adjust_rl)
                  model_loop.rl$return_level_low = as.numeric(-rl_loop[1] -adjust_rl)
                  model_loop.rl$return_level_high = as.numeric(-rl_loop[3] -adjust_rl)}
              }
              model_loop.rl$return_horizon = return_horizon_loop
              model_loop.rl$return_period = return_period
              model_loop.rl$station = st
              model_loop.rl$scenario = sc
              model_loop.rl$CMIP_model = CMIP_mod
              model_loop.rl$ns_model = ns_mod
              
              ns_models.gev.rl = rbind(ns_models.gev.rl, model_loop.rl)
              #}
            }
          }
        }
      }
    }
  }
  if(rl){
    return(list(params = ns_models.gev.params, return_levels = ns_models.gev.rl))
  }
  else{
    return(params = ns_models.gev.params)
  }
}

#The function below perform the selection procedure based on shape parameter and the upper threshold value
#The criteriums are the following ones : 
#     No more than 10 values (in daily) above the threshold (50 for tasmax)
#     Bounded distribution : xi<0
#     Upper bound below a given threshould (70 for tasmax)
#     Lower bound above -(À for tasmin)
selection_models = function(tas_data, TX_data, ns_params, ns_rl = FALSE, variable){
  CMIP_models = unique(ns_params$CMIP_model)
  stations = unique(ns_params$station)
  
  ns_params_final = data.frame(ns_params)
  if(ns_rl != FALSE){
    ns_rl_final = data.frame(ns_rl)
  }
  
  for(CMIP_mod in CMIP_models){
    for(st in stations){
      scenarios_loop = unique(ns_params_final %>% filter(CMIP_model==CMIP_mod, station==st) %>% select(scenario))
      for(sc in scenarios_loop[[1]]){
        #First test to check if no more than 10 daily data points are above the threshold
        if(variable=='tasmax' & nrow(tas_data %>% filter(CMIP_model==CMIP_mod, station==st, scenario==sc, var>50)) >= 10){
          ns_params_final = subset(ns_params_final, CMIP_model!=CMIP_mod | scenario!=sc| station != st)
          if(ns_rl != FALSE){
            ns_rl_final = subset(ns_rl_final, CMIP_model!=CMIP_mod | scenario!=sc| station != st)
          }
        }
        else if(variable=='tasmin' & nrow(tas_data %>% filter(CMIP_model==CMIP_mod, station==st, scenario==sc, var < -30)) >= 10){
          ns_params_final = subset(ns_params_final, CMIP_model!=CMIP_mod | scenario!=sc| station != st)
          if(ns_rl != FALSE){
            ns_rl_final = subset(ns_rl_final, CMIP_model!=CMIP_mod | scenario!=sc| station != st)
          }
        }
        
        ns_mod_loop = unique(ns_params_final %>% filter(CMIP_model==CMIP_mod, station==st, scenario==sc) %>% select(ns_model))
        for(ns_mod in ns_mod_loop[[1]]){
          years_loop = unique(TX_data %>%
                                filter(CMIP_model==CMIP_mod, scenario==sc, station==st) %>%
                                select(year))
          cov_loop = TX_data %>% filter(CMIP_model==CMIP_mod, scenario==sc, station==st) %>% select(covariate)
          
          ns_params_loop = ns_params_final %>% filter(CMIP_model==CMIP_mod, scenario==sc, station==st, ns_model==ns_mod)
          loc_loop = ns_params_loop$loc_0 + ns_params_loop$loc_1*cov_loop + ns_params_loop$a*cov_loop + ns_params_loop$b
          if(ns_mod %in% c('M000', 'M100', 'M-100')){
            scale_loop = ns_params_loop$scale_0
          }
          else{
            scale_loop = exp(ns_params_loop$scale_0 + ns_params_loop$scale_1*cov_loop)
          }
          shape_loop = ns_params_loop$shape
          #upper_bound_loop = loc_loop - scale_loop/shape_loop
  
          #Test for selection of final models
          if(shape_loop > 0 | is.na(shape_loop)){
            ns_params_final = subset(ns_params_final, CMIP_model!=CMIP_mod | scenario!=sc| station != st | ns_model !=ns_mod)
          }
          else if(variable=='tasmax' & max(loc_loop - scale_loop/shape_loop) >= 70){
            ns_params_final = subset(ns_params_final, CMIP_model!=CMIP_mod | scenario!=sc| station != st | ns_model !=ns_mod)
          }
          else if(variable=='tasmin' & max(loc_loop - scale_loop/shape_loop) >= 50){
            ns_params_final = subset(ns_params_final, CMIP_model!=CMIP_mod | scenario!=sc| station != st | ns_model !=ns_mod)
          }
          
          if(ns_rl != FALSE){
            if(shape_loop > 0 | is.na(shape_loop)){
              ns_rl_final = subset(ns_rl_final, CMIP_model!=CMIP_mod | scenario!=sc| station != st | ns_model !=ns_mod)
            }
            else if(variable=='tasmax' & max(loc_loop - scale_loop/shape_loop) >= 70){
              ns_rl_final = subset(ns_rl_final, CMIP_model!=CMIP_mod | scenario!=sc| station != st | ns_model !=ns_mod)
            }
            else if(variable=='tasmin' & max(loc_loop - scale_loop/shape_loop) >= 50){
              ns_rl_final = subset(ns_rl_final, CMIP_model!=CMIP_mod | scenario!=sc| station != st | ns_model !=ns_mod)
            }
          }
        }
      }
    }
  }
  return(list(params = ns_params_final, return_levels =  ns_rl_final))
}

#The function below compute the GEV fit as a function of time 
#Input : 
#   ns_model.gev.params : It corresponds to the output of the gev_ns_rl_params_fct
#   TX_data_with_cov : Data used as the input of gev_ns_rl_params_fct
gev_regression = function(ns_model.gev.params, TX_data_with_cov){
  ns_model.evol = data.frame(matrix(nrow=0, ncol=17))
  colnames(ns_model.evol) = c('year', 'station', 'scenario', 'CMIP_model', 'ns_model', 'loc', 'loc_low', 'loc_high', 'scale', 'scale_low', 'scale_high','shape','shape_low', 'shape_high', 'upper_bound', 'upper_bound_low', 'upper_bound_high')
  
  stations= unique(ns_model.gev.params$station)
  CMIP_models = unique(ns_model.gev.params$CMIP_model)
  for(CMIP_mod in CMIP_models){
    for(st in stations){
      scenarios_loop = unique(ns_model.gev.params %>% filter(CMIP_model==CMIP_mod, station==st) %>% select(scenario))
      scenarios_loop = scenarios_loop[scenarios_loop!= 'historical']
      for(sc in scenarios_loop){
        ns_models_loop = unique(ns_model.gev.params %>% filter(CMIP_model==CMIP_mod, station==st, scenario==sc) %>% select(ns_model))
        for(ns_mod in ns_models_loop[[1]]){
          #GMST
          years_loop = unique(TX_data_with_cov %>%
                                filter(CMIP_model==CMIP_mod, scenario==sc, station==st) %>%
                                select(year))
          cov_loop = TX_data_with_cov %>% filter(CMIP_model==CMIP_mod, scenario==sc, station==st) %>% select(covariate)
          ns_model_loop = ns_model.gev.params %>% filter(CMIP_model==CMIP_mod, scenario==sc, station==st, ns_model==ns_mod)
          
          if(length(ns_model_loop$shape)==0){}
          else if(is.na(ns_model_loop$shape) | ns_model_loop$shape >0 ){}
          else{
            if(ns_mod %in% c('M000', 'M100', 'M-100')){      
              df_loop = data.frame(year = years_loop,
                                        station=st,
                                        scenario=sc,
                                        CMIP_model = CMIP_mod,
                                        ns_model=ns_mod,
                                        loc=ns_model_loop$loc_0 + ns_model_loop$loc_1*cov_loop +
                                          ns_model_loop$a*cov_loop + ns_model_loop$b,
                                        loc_low = ns_model_loop$loc_0_low + ns_model_loop$loc_1_low*cov_loop +
                                          (ns_model_loop$a - ns_model_loop$a_std)*cov_loop + (ns_model_loop$b - ns_model_loop$b_std),
                                        loc_high = ns_model_loop$loc_0_high + ns_model_loop$loc_1_high*cov_loop +
                                          (ns_model_loop$a + ns_model_loop$a_std)*cov_loop + (ns_model_loop$b + ns_model_loop$b_std),
                                        scale = ns_model_loop$scale_0 + ns_model_loop$scale_1*cov_loop,
                                        scale_low = ns_model_loop$scale_0_low + ns_model_loop$scale_1_high*cov_loop,
                                        scale_high = ns_model_loop$scale_0_low + ns_model_loop$scale_1_high*cov_loop,
                                        shape = ns_model_loop$shape,
                                        shape_low = ns_model_loop$shape_low,
                                        shape_high = ns_model_loop$shape_high,
                                        upper_bound = 0,
                                        upper_bound_low = 0,
                                        upper_bound_high = 0)}
            else{df_loop = data.frame(year = years_loop,
                                           station=st,
                                           scenario=sc,
                                           CMIP_model = CMIP_mod,
                                           ns_model=ns_mod,
                                           loc=ns_model_loop$loc_0 + ns_model_loop$loc_1*cov_loop +
                                             ns_model_loop$a*cov_loop + ns_model_loop$b,
                                           loc_low = ns_model_loop$loc_0_low + ns_model_loop$loc_1_low*cov_loop +
                                             (ns_model_loop$a - ns_model_loop$a_std)*cov_loop + (ns_model_loop$b - ns_model_loop$b_std),
                                           loc_high = ns_model_loop$loc_0_high + ns_model_loop$loc_1_high*cov_loop +
                                             (ns_model_loop$a + ns_model_loop$a_std)*cov_loop + (ns_model_loop$b + ns_model_loop$b_std),
                                           scale = exp(ns_model_loop$scale_0 + ns_model_loop$scale_1*cov_loop),
                                           scale_low = exp(ns_model_loop$scale_0_low + ns_model_loop$scale_1_high*cov_loop),
                                           scale_high = exp(ns_model_loop$scale_0_low + ns_model_loop$scale_1_high*cov_loop),
                                           shape = ns_model_loop$shape,
                                           shape_low = ns_model_loop$shape_low,
                                           shape_high = ns_model_loop$shape_high,
                                           upper_bound = 0,
                                           upper_bound_low = 0,
                                           upper_bound_high = 0)
            }
            colnames(df_loop) = c('year', 'station', 'scenario', 'CMIP_model', 'ns_model', 'loc', 'loc_low', 'loc_high', 'scale', 'scale_low', 'scale_high','shape','shape_low', 'shape_high', 'upper_bound', 'upper_bound_low', 'upper_bound_high')
            
            upper_bound_loop = df_loop$loc - df_loop$scale/df_loop$shape
            upper_bound_low_loop = df_loop$loc_low - df_loop$scale_high/df_loop$shape_low
            
            if(df_loop$shape_high<0){
              upper_bound_high_loop = df_loop$loc_high - df_loop$scale_low/df_loop$shape_high
            }
            else{
              upper_bound_high_loop = max(upper_bound_loop) + 100
            }
            
            df_loop = df_loop%>% mutate(upper_bound = upper_bound_loop,
                                                  upper_bound_low = upper_bound_low_loop,
                                                  upper_bound_high = upper_bound_high_loop)
            
            ns_model.evol = rbind(ns_model.evol, df_loop)
          }
        }
      }
    }
  }
  return(ns_model.evol)
}
