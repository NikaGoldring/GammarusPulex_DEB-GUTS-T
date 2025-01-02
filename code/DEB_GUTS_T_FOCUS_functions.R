## Function to prepare work environment
checkPackages <- function(required.packages = c("rlang","tidyverse","scales", "gridExtra","patchwork","ggridges","gstat",
                                                "colorspace","ggpubr", "labeling", "farver"), library.loc){
  .libPaths(new = c(library.loc,.libPaths()))
  if(length(installed.packages(lib.loc = library.loc))==0){
    packages <- installed.packages(lib.loc = paste0(R.home(),"/library"))
  }else{packages <- rbind(installed.packages(lib.loc = library.loc),installed.packages(lib.loc = paste0(R.home(),"/library")))}
  
  required.packages <- required.packages
  not.installed <- required.packages[!required.packages %in% packages[,"Package"]]
  if(length(not.installed)==0){
    print("Required packages installed.")
    return(TRUE)}else{
      print("Additional packages needed, please install the following packages and re-run function:")
      print(not.installed)
      return(FALSE)
    }
}

loadPackages <- function(required.packages = c("rlang","tidyverse","scales", "gridExtra","patchwork","ggridges","gstat",
                                               "colorspace","ggpubr", "labeling", "farver"), library.loc,required.packages.installed){
  .libPaths(new = c(library.loc,.libPaths()))
  
  if(required.packages.installed){
    lapply(required.packages,function(x){require(x, character.only = T,lib.loc = library.loc)})}else{
      print("Additional packages needed, please install the required packages and re-run function checkPackages.")
    }
}

## Function to read in data that is created from DEB Gammarus with Hans' Smalltalk implementation
## Applying filter to get the desired subset
readData <- function(data.location, 
                     guts.model.version,
                     ignore.version,
                     ignore.chemical = F,
                     filter.concentrations = T, 
                     filter.temp.amplitudes = T,
                     desired.exposure.concentration,
                     desired.temp.amplitude,
                     application.pulse.shift = F){
  
  # List all folders with matching names for the guts.model.version (can be a list of multiple entries)
  # Returns a vector with folders and full paths to folders
  results.list <- lapply(guts.model.version,function(x){ 
    if(application.pulse.shift){
      c(list.dirs(path = data.location,recursive = F)[grepl(x,x = list.dirs(path = data.location,recursive = F,full.names = F))],
        list.dirs(path = paste0(data.location,"pulse_shift"),recursive = F)[grepl(x,x = list.dirs(path = paste0(data.location,"pulse_shift"),recursive = F,full.names = F))])}
    else{
        list.dirs(path = data.location,recursive = F)[grepl(x,x = list.dirs(path = data.location,recursive = F,full.names = F))]
        }}) %>% 
    do.call(rbind,.) %>% as.vector()
  
  # Exclude (therefore it's called ignore) the folders of either constant = "verylow", or pulsed exposure. 
  if(ignore.version=="verylow"){
    results.list <- results.list[!grepl(ignore.version,results.list)]
  }else{
    results.list <- results.list[grepl("verylow",results.list)]
  }
  
  # Exclude (therefore it's called ignore) the folders of FPF simulations, or IMI simulations.
  if(ignore.chemical=="FPF"){
    results.list <- results.list[!grepl(ignore.chemical,results.list)]
    exposure.chemical <- "IMI"
  }else{
    results.list <- results.list[grepl("FPF",results.list)]
    exposure.chemical <- "FPF"
  }
  
  # Create a data list, looping (with lapply) through all results folders and extracting required information on 
  # number of individuals over time (mean and sd of replicates)
  df <- lapply(results.list,function(x){
    
    # extract information on all scenarios from the .modelscript in a result.list folder
    model.script <- list.files(path = x,pattern = "modelscript" ,full.names = T)
    scenarios <- read.delim(model.script,header = F)
    
    # create a vector with numbers of same length as scenarios
    scenario.number <- 1:length(as.matrix(scenarios))
    
    # Is filter.concentrations set to TRUE, then filter scenario.number to keep only desired concentrations
    if(filter.concentrations){
      
      scenario.numbers <- scenario.number[grepl(x = as.matrix(scenarios), pattern = paste0("exposureConc: ",desired.exposure.concentration))]
      
    }
    
    # Is filter.temp.amplitudes set to TRUE, then filter scenario.number to keep only desired temperature amplitudes
    if(filter.temp.amplitudes){
      
      scenario.numbers.2 <- scenario.number[grepl(x = as.matrix(scenarios),pattern = paste0("envTamp: ",desired.temp.amplitude))]
      
    }
    
    # Create a vector "ls" with all folder names, these are the folders where the text output files are stored
    ls <- paste0("x1s",1:length(list.dirs(path = paste0(x,"/x1/"),recursive = F,full.names = F)))
    
    # match scenario.numbers and scenario.numbers.2 to keep only ls entries that match the desired concentrations and temp amplitudes 
    if(filter.concentrations & filter.temp.amplitudes){
      ls <- ls[scenario.numbers.2[scenario.numbers.2 %in% scenario.numbers]]
      scenarios <- as.matrix(scenarios)[scenario.numbers.2[scenario.numbers.2 %in% scenario.numbers]]}
    if(filter.concentrations & !filter.temp.amplitudes){
      ls <- ls[scenario.numbers]
      scenarios <- as.matrix(scenarios)[scenario.numbers]
    }
    if(!filter.concentrations & filter.temp.amplitudes){
      ls <- ls[scenario.numbers.2]
      scenarios <- as.matrix(scenarios)[scenario.numbers.2]
    }
    
    # Extract some run info
    run.info <- read.csv(list.files(path = paste0(x[1],"/x1/x1s1/"),
                                    pattern = "ModelSystem_control.csv",full.names = T),
                         header = F,stringsAsFactors = F)[,1:2]
    temp.info <- pivot_wider(run.info,names_from = V1, values_from = V2) %>% as.data.frame()
    
    #### Loop through ls (i.e., the folders containing DEB simulation outputs)
    output <- lapply(ls, function(y) {  # for each element y in ls 
      # Store the temperature information used in the simulations 
      envT <- read.delim(list.files(path = paste0(x,"/x1/",y), pattern = "environment",full.names = T)[1])
        
      # Only look for folders with "stagestructureEmbJuvAdultsInds"
      # these contain, for each day in the simulation, the number of individuals.
      ls.data.files <- list.files(path = paste0(x,"/x1/",y), pattern = "stageStructureEmbrJuvAdultsInds",full.names = T)
      
      input <- lapply(ls.data.files, function(z){ 
        ind <- read.delim(file = z)[,4] # More specifically: only keep column 4 of the data which is the total number of all stages combined.
        ind
          }
      ) %>% do.call(cbind,.) %>% as.matrix()
        
      # some scenarios don't survive the first period, these should only used if >=3 of the replicates survive the first year
      input <- as.matrix(input[,apply(as.matrix(input[(temp.info$`startApplicationYear:` - temp.info$`startYear:`)*365,]),2,FUN = function(x) x>0)])
        
      # calculate the means and SDs for each of the replicates
      if(application.pulse.shift){
        if(ncol(input)>=3){
          output.scenario <- data.frame(mean = apply(input,1,FUN = mean), sd = apply(input,1,FUN = sd),
                                          envT = envT$temperature,
                                          scenario.id = y,
                                          model.version = substr(strsplit(strsplit(x,split = "/")[[1]][length(strsplit(x,split = "/")[[1]])],split = "fGammarus")[[1]][1],start = 3,
                                                                 stop = nchar(strsplit(strsplit(x,split = "/")[[1]][length(strsplit(x,split = "/")[[1]])],split = "fGammarus")[[1]][1])),
                                          time_shift = as.numeric(strsplit(strsplit(strsplit(x,split = "/")[[1]][length(strsplit(x,split = "/")[[1]])],
                                                                                    split = "_")[[1]][5],"t")[[1]][2]))}
        else{
          output.scenario <- data.frame(mean = NA, sd = NA,
                                          envT = envT$temperature,
                                          scenario.id = y,
                                          model.version = substr(strsplit(strsplit(x,split = "/")[[1]][length(strsplit(x,split = "/")[[1]])],split = "fGammarus")[[1]][1],start = 3,
                                                                 stop = nchar(strsplit(strsplit(x,split = "/")[[1]][length(strsplit(x,split = "/")[[1]])],split = "fGammarus")[[1]][1])),
                                          time_shift = as.numeric(strsplit(strsplit(strsplit(x,split = "/")[[1]][[length(strsplit(x,split = "/")[[1]])]],
                                                                                    split = "_")[[1]][5],"t")[[1]][2]))
        }
      }
      else{
        if(ncol(input)>=3){
          output.scenario <- data.frame(mean = apply(input,1,FUN = mean), sd = apply(input,1,FUN = sd),
                                          envT = envT$temperature,
                                          scenario.id = y,
                                          model.version = substr(strsplit(strsplit(x,split = "/")[[1]][length(strsplit(x,split = "/")[[1]])],split = "fGammarus")[[1]][1],start = 3,
                                                                 stop = nchar(strsplit(strsplit(x,split = "/")[[1]][length(strsplit(x,split = "/")[[1]])],split = "fGammarus")[[1]][1])), # Hard coded file path, needs adjustment  ## Maybe use data.location that is provided as function argument
                                          exposure.chemical = exposure.chemical,
                                          T.scenario = sub(".*_", "", data.location) # get T.sceanrio from data.location string
                                        )} 
        else{
          output.scenario <- data.frame(mean = NA, sd = NA,
                                        envT = envT$temperature,
                                        scenario.id = y,
                                        model.version = substr(strsplit(strsplit(x,split = "/")[[1]][length(strsplit(x,split = "/")[[1]])],split = "fGammarus")[[1]][1],start = 3,
                                                                 stop = nchar(strsplit(strsplit(x,split = "/")[[1]][length(strsplit(x,split = "/")[[1]])],split = "fGammarus")[[1]][1])), # Hard coded file path, needs adjustment  
                                        exposure.chemical = exposure.chemical,
                                        T.scenario = sub(".*_", "", data.location) # get T.sceanrio from data.location string
                                        )}
      }
      output.scenario
    })
    
    # Return the output
    return(list(scenarios = scenarios, data = output, run.info = run.info))
  })
  return(df)
}

## Helper function to extract and filter scenarios
extract_and_filter_scenarios <- function(df.list, 
                                         desired.exposure.concentration){
  scens <- as.matrix(df.list[[1]]$scenarios$V1)
  scens <- apply(scens, 1, function(z) {
    as.numeric(
      substr(
        unlist(strsplit(z, " "))[8], # Adjust index based on input structure
        start = 1,
        stop = nchar(unlist(strsplit(z, " "))[8]) - 1
      )
    )
  })
  scens <- data.frame(exposureConcentration = scens)
  
  # Change scens to applied filter 
  if (!is.null(desired.exposure.concentration)) {
    select.scens <- scens$exposureConcentration %in% desired.exposure.concentration
    scens <- scens[select.scens, , drop = FALSE]
  }
return(list(scens=scens, select.scens = select.scens)) #Returns all or selected scens as df and select.scen as logi for filtering
}

## Helper function to process a list of data frames
process_model_data <- function(df.list, 
                               desired.exposure.concentration, 
                               time.range) {

  ## Get scenario info to define application window
  startYear <- df.list[[1]]$run.info$V2[df.list[[1]]$run.info$V1 == "startYear:"]
  endYear <- df.list[[1]]$run.info$V2[df.list[[1]]$run.info$V1 == "endYear:"]
  startApplicationYear <- df.list[[1]]$run.info$V2[df.list[[1]]$run.info$V1 == "startApplicationYear:"]
  endApplicationYear <- df.list[[1]]$run.info$V2[df.list[[1]]$run.info$V1 == "endApplicationYear:"]
  
  days.application <- seq(as.Date(paste(startApplicationYear, "/1/1", sep = "")),
                          as.Date(paste(endApplicationYear, "/12/31", sep = "")), "days")
  days.full.period <- seq(as.Date(paste(startYear, "/1/1", sep = "")),
                          as.Date(paste(endYear, "/12/31", sep = "")), "days")
  
  ## Prepare df for each model version to apply filters 
  # Get filtered scenario infos
  filtered <- extract_and_filter_scenarios(df.list, desired.exposure.concentration)
  scens <- filtered$scens
  select.scens <-filtered$select.scens
  
  # Individual filtered data frames
  df.SD <- df.list[[1]]$data[select.scens]
  df.SDT <- df.list[[2]]$data[select.scens]
  df.SDTStd <- df.list[[3]]$data[select.scens]
  
  ##Select data of filtered df within application window only
  df.SD <- lapply(1:nrow(scens),function(x){
    out <- df.SD[[x]]
    if(!nrow(out)==1){
      
      out$date <- days.full.period
      out <- out[out$date %in% days.application,]
      out$exposureConc <- scens$exposureConcentration[x]
    }
    else{
      out$date <- NA
      out$exposureConc <- scens$exposureConcentration[x]
    }
    out}) %>% do.call(rbind,.)
  
  df.SDT <- lapply(1:nrow(scens),function(x){
    out <- df.SDT[[x]]
    if(!nrow(out)==1){
      out$date <- days.full.period
      out <- out[out$date %in% days.application,]
      out$exposureConc <- scens$exposureConcentration[x]
    }
    else{
      out$date <- NA
      out$exposureConc <- scens$exposureConcentration[x]
    }
    out}) %>% do.call(rbind,.)
  
  df.SDTStd <- lapply(1:nrow(scens),function(x){
    out <- df.SDTStd[[x]]
    if(!nrow(out)==1){
      out$date <- days.full.period
      out <- out[out$date %in% days.application,]
      out$exposureConc <- scens$exposureConcentration[x]
    }
    else{
      out$date <- NA
      out$exposureConc <- scens$exposureConcentration[x]
    }
    out}) %>% do.call(rbind,.)
  
  ##Finally get data within defined time range only
  h <- as.numeric(time.range)
  # if(!grepl(",",h)){
  #   h <- as.numeric(h)}else{
  #     h <- as.numeric(unlist(strsplit(h,",")))
  #   }
  
  df.SD <- df.SD[format(df.SD$date,"%Y") %in% h,]
  df.SD$year <- format(df.SD$date,"%Y")
  
  df.SDT <- df.SDT[format(df.SDT$date,"%Y") %in% h,]
  df.SDT$year <- format(df.SDT$date,"%Y")
  
  df.SDTStd <- df.SDTStd[format(df.SDTStd$date,"%Y") %in% h,]
  df.SDTStd$year <- format(df.SDTStd$date,"%Y")

  return(list(SD = df.SD, SDT= df.SDT, SDTStd = df.SDTStd)) ##Returns filtered dfs
}

## Function to get population size at the end of the scenarios, as well as quantile information
popSize <- function (simulation.data.list, application.pulse.shift=F){
  #### this part is not adjusted yet.
  if(application.pulse.shift){
    scens <- as.matrix(x$scenarios)
    scens <- apply(scens,1,function(z) as.numeric(apply(as.matrix(unlist(strsplit(z," "))[c(8,11,13)]),1, ### what is 13? ->pulse shift
                                                        function(x) substr(x,start = 1,nchar(x)-1))))
    scens <- t(scens)
    scens <- as.data.frame(scens)
    names(scens) <- c("exposureConcentration","tempAmplitude","pulseShift")
    
    out <- lapply(simulation.data.list,function(x){
      
      # extract information on SD vs SDT (IT vs ITT)
      m.type <- strsplit(x$scenarios[1][1,],split = "Gammarus")[[1]][1]
      m.type <- substr(m.type,start = 4, stop = nchar(m.type))
      
      
      df1 <- lapply(x$data, function(y){
        if(sum(is.na(y$mean))==0){
          # Calculate quantiles ###hardcoded numbers needs adjustment!!! here 3 years warm up and 6 years without application in the end ##### HARDCODE - needs change
          # We can use the ModeslSystem control csv for this actually, this holds start and end year of application 
          q.df <- t(as.matrix(quantile(y$mean[1096:(nrow(y)-2190)],
                                       seq(0.01,1,0.01)))) #Calculating 100 quantiles
          colnames(q.df) <- paste0(m.type,"_Q",gsub(pattern = "%","",colnames(q.df)))
          
          # calculate last day endpoint
          ld <- as.matrix(y$mean[nrow(y)])
          colnames(ld) <- paste0(m.type,"_LastDayN")
          
          as.data.frame(cbind(ld,q.df))}
        else{
          q.df <- t(as.matrix(rep(NA,100)))
          colnames(q.df) <- paste0(m.type,"_Q",1:100)
          ld <- as.matrix(NA)
          colnames(ld) <- paste0(m.type,"_LastDayN")
          
          as.data.frame(cbind(ld,q.df))
        }
      }) %>% do.call(rbind,.)
      df1 <- cbind(as.data.frame(df1),scens)
    })
    return(out)}else{
      scens <- simulation.data.list[[1]]$scenarios
      scens <- apply(scens,1, #the 1 as the second argument in apply function means the function is being applied to each row 
                     function(z) as.numeric(  #save as numeric value
                       apply(as.matrix(unlist(strsplit(z," ")) #splitting scenario string by spaces then unlisting them into a vector
                                       [c(8)]),1,     #selecting element(s) of the vector ##Needs adjustment based on input string
                             #here: 8 = exposure concentration
                             function(x) substr(x, start = 1, nchar(x) - 1))))  # removes the % from each line (don't know why yet..)
      
      scens <- t(scens)  #transpose 
      scens <- as.data.frame(scens)
      names(scens) <- c("exposureConcentration")
      
      out <- lapply(simulation.data.list,function(x){
        
        # extract information on SD vs SDT (IT vs ITT)
        m.type <- strsplit(x$scenarios[1][1,],split = "Gammarus")[[1]][1]
        m.type <- substr(m.type,start = 4, stop = nchar(m.type))
        
        
        df1 <- lapply(x$data, function(y){
          if(sum(is.na(y$mean))==0){
            # Calculate quantiles
            q.df <- t(as.matrix(quantile(y$mean[1096:(nrow(y)-2190)],seq(0.01,1,0.01))))
            colnames(q.df) <- paste0(m.type,"_Q",gsub(pattern = "%","",colnames(q.df)))
            
            # calculate last day endpoint
            ld <- as.matrix(y$mean[nrow(y)])
            colnames(ld) <- paste0(m.type,"_LastDayN")
            
            as.data.frame(cbind(ld,q.df))}
          else{
            q.df <- t(as.matrix(rep(NA,100)))
            colnames(q.df) <- paste0(m.type,"_Q",1:100)
            ld <- as.matrix(NA)
            colnames(ld) <- paste0(m.type,"_LastDayN")
            
            as.data.frame(cbind(ld,q.df))
          }
        }) %>% do.call(rbind,.)
        df1
      }) %>% do.call(cbind,.)
      out <- as.data.frame(out)
      out <- cbind(out,scens)
      rownames(out) <- NULL
      return(out)
    }
}


checkParamArrheniusCorrection <- function(){
  
  # quick function to calculate temperature pattern
  temps <- function(tAmplitude, length.time.series){
    15 - tAmplitude * cos(((c(1:length.time.series) - 31)/365)*2*pi)
  }
  
  hbs <- function(temps){
    hb.sd <- 0.007848 * exp((7450/(20+273.15))-(7450/(temps + 273.15)))
    hb.it <- 0.01157 * exp((11020/(20+273.15))-(11020/(temps + 273.15)))
    kd.sd <- 0.02982 * exp((7450/(20+273.15))-(7450/(temps + 273.15)))
    kd.it <- 1e-6 * exp((11020/(20+273.15))-(11020/(temps + 273.15)))
    data.frame(hb.sd = hb.sd,hb.it = hb.it, kd.sd = kd.sd, kd.it = kd.it)
  }
  
  hb.tamp1 <- hbs(temps = temps(tAmplitude = 1,length.time.series = 365))
  hb.tamp11 <- hbs(temps = temps(10,365))
  
  plot(1:365,hb.tamp11$hb.it,type = "l")
  abline(h = 0.01157,lty = 6)
  lines(1:365,hb.tamp1$hb.it,col = 'red')
  legend("topright",legend=c("hb_Tamp = 10", "hb_Tamp = 1","hb_Tref_20"),
         col=c("black", "red","black"), lty=c(1,1,6), cex=0.8, title = "Model: IT")
  
  plot(1:365,hb.tamp11$hb.sd,type = "l")
  abline(h = 0.007848,lty = 6)
  lines(1:365,hb.tamp1$hb.sd,col = 'red')
  legend("topright",legend=c("hb_Tamp = 10", "hb_Tamp = 1","hb_Tref_20"),
         col=c("black", "red","black"), lty=c(1,1,6), cex=0.8, title = "Model: SD")
  
  kd.tamp1 <- hbs(temps = temps(tAmplitude = 1,length.time.series = 365))
  kd.tamp11 <- hbs(temps = temps(10,365))
  
  plot(1:365,kd.tamp11$kd.it,type = "l")
  abline(h = 1e-6,lty = 6)
  lines(1:365,kd.tamp1$kd.it,col = 'red')
  legend("topright",legend=c("kd_Tamp = 10", "kd_Tamp = 1","kd_Tref_20"),
         col=c("black", "red","black"), lty=c(1,1,6), cex=0.8, title = "Model: IT")
  
  plot(1:365,kd.tamp11$kd.sd,xaxt = "n",type = "l",
       xlab = "Day of the year [d]",ylab = "Dominant rate constant [d-1]")
  axis(1, at = seq(0,360,20))
  abline(h = 0.02982,lty = 6)
  abline(v = 100, lty = 2)
  lines(1:365,kd.tamp1$kd.sd,col = 'red')
  legend("topright",legend=c("Amplitude = 10", "Amplitude = 1","Reference 20C"),
         col=c("black", "red","black"), lty=c(1,1,6), cex=0.8, title = "Model: SD")
}


#### Hier verder: application.pulse.shift toevoegen, df_SD list item 1:4 SD -20 20 40 60, df_SD list item 5:8 SDT.
### Plots voor verschil tussen de twee: dus 1+5, 2+6, 3+7, 4+8
### OOK plots voor verschil binnen de groepen? 1:4 en 5:8?

## Function for plotting two model versions
plotTAmpPopDynamics <- function(df_SD.list,
                                exposure.type,
                                desired.exposure.concentrations,
                                desired.temp.amplitudes,
                                time.range,
                                application.pulse.shift = F,
                                y.trans = F){
  
  if(application.pulse.shift){
    # full data set
    df <- df_SD.list
    
    # select all scenarios matching concentrations and T-amplitudes ->> This should go into a seperate function that is called here and in other functions, to avoid correcting it in multiple places ic code needs changing
    scens <- as.matrix(df[[1]]$scenarios$V1)
    scens <- apply(scens,1, #the 1 as the second argument in apply function means the function is being applied to each row 
                   function(z) as.numeric(  #save as numeric value
                     apply(as.matrix(unlist(strsplit(z," ")) #splitting scenario string by spaces then unlisting them into a vector
                                     [c(8)]),1,     #selecting element(s) of the vector ##Needs adjustment based on input string
                           #here: 8 = exposure concentration
                           function(x) substr(x, start = 1, nchar(x) - 1))))  # removes the % from each line (don't know why yet..)
    scens <- t(scens)  #transpose 
    scens <- as.data.frame(scens)
    names(scens) <- c("exposureConcentration")
    
    # create boolean for scenarios
    select.scens <- scens$exposureConcentration %in% desired.exposure.concentrations & scens$tempAmplitude %in% desired.temp.amplitudes
    scens <- scens[select.scens,]
    
    # select only application window data and add scenario info
    ## Note: this seems trivial, but sadly isn't....
    ## Reason: simulations appear to use actual date information, which means we have to account for leap-years
    ## this can only be done by making a time-series of actual dates, for which we need the run control info
    startYear <- df[[1]]$run.info$V2[df[[1]]$run.info$V1=="startYear:"]
    endYear <- df[[1]]$run.info$V2[df[[1]]$run.info$V1=="endYear:"]
    startApplicationYear <- df[[1]]$run.info$V2[df[[1]]$run.info$V1=="startApplicationYear:"]
    endApplicationYear <- df[[1]]$run.info$V2[df[[1]]$run.info$V1=="endApplicationYear:"]
    
    days.warmup <- seq(as.Date(paste(startYear,"/1/1",sep = "")),as.Date(paste(startApplicationYear-1,"/12/31",sep = "")),"days")
    days.application <- seq(as.Date(paste(startApplicationYear,"/1/1",sep = "")),as.Date(paste(endApplicationYear,"/12/31",sep = "")),"days")
    days.recovery <- seq(as.Date(paste(endApplicationYear + 1,"/1/1",sep = "")),as.Date(paste(endYear,"/12/31",sep = "")),"days")
    days.full.period <- seq(as.Date(paste(startYear,"/1/1",sep = "")),as.Date(paste(endYear,"/12/31",sep = "")),"days")
    
    # select all dataframes matching the desired concentrations and T-amplitudes by matching multiple string patterns
    df <- lapply(1:length(df),function(x) df[[x]]$data[select.scens][[1]])
    types <- lapply(df, function(x) unique(x$model.version)) %>% do.call(rbind,.)
    
    df.SD <- lapply(which(types == "SD"),function(x){
      out <- df[[x]]
      if(!nrow(out)==1){
        out$date <- days.full.period
        out <- out[out$date %in% days.application,]
        out$exposureConc <- scens$exposureConcentration
      }
      else{
        out$date <- NA
        out$exposureConc <- scens$exposureConcentration[x]
      }
      out}) %>% do.call(rbind,.)
    df.SDT <- lapply(which(types == "SDT"),function(x){
      out <- df[[x]]
      if(!nrow(out)==1){
        out$date <- days.full.period
        out <- out[out$date %in% days.application,]
        out$exposureConc <- scens$exposureConcentration
      }
      else{
        out$date <- NA
        out$exposureConc <- scens$exposureConcentration[x]
      }
      out}) %>% do.call(rbind,.)
    
    h <- time.range
    if(!grepl(",",h)){
      h <- as.numeric(h)}else{
        h <- as.numeric(unlist(strsplit(h,",")))
      }
    
    df.SD <- df.SD[format(df.SD$date,"%Y") %in% h,]
    df.SD$year <- format(df.SD$date,"%Y")
    df.SD$time_shift[is.na(df.SD$time_shift)] <- 0
    
    df.SDT <- df.SDT[format(df.SDT$date,"%Y") %in% h,]
    df.SDT$year <- format(df.SDT$date,"%Y")
    df.SDT$time_shift[is.na(df.SDT$time_shift)] <- 0
    
    # # data frame with temperature data
    # Ts <- data.frame(x = 1:365, apply(as.matrix(unique(df.SD$envTamp)),1,function(x){temps(x,365)})) 
    # Ts <- pivot_longer(Ts,cols = -1,names_to = "grp",values_to = "Temperature")
    # if(length(unique(df.SD$envTamp))==1){
    #   Ts$grp <- "X1"
    # }
    # Ts <- left_join(Ts,data.frame(grp = paste0("X",1:length(unique(df.SD$envTamp))),envTamp = unique(df.SD$envTamp)))
    
    # create a colour scale for the exposure concentrations
    cols <- grDevices::colorRampPalette(colors = c("darkblue","red"))
    cols <- cols(length(unique(df.SD$time_shift)))
    
    p1 <- ggplot(df.SD) +
      # geom_rect(aes(xmin = 100,ymin = 0, xmax = 200, ymax = max(df.SD$mean)),fill = "grey75") +
      geom_ribbon(aes(x = as.numeric(format(date,"%j")), ymin = mean - sd,ymax = mean + sd, group = as.factor(time_shift),
                      fill = as.factor(time_shift)),alpha = 0.25, show.legend = F) +
      geom_line(aes(x = as.numeric(format(date,"%j")), y = mean, group = as.factor(time_shift),
                    colour = as.factor(time_shift)),lwd = 1) +
      scale_fill_manual(values = cols) +
      scale_colour_manual(values = cols) +
      guides(colour = "none",linetype = "none") + 
      scale_x_continuous("Day of the year",breaks = seq(0,360,20),guide = guide_axis(n.dodge = 2)) +
      ylab("Mean population size") +
      geom_vline(xintercept = c(100),linetype = 6,lwd = 0.5, show.legend = F) +
      # facet_wrap(.~envTamp,nrow = length(h)) + 
      ggtitle("DEB-GUTS") + 
      coord_cartesian(ylim = c(0,1100),expand = T) +
      theme_pubr() 
    
    p2 <- ggplot(df.SDT) +
      # geom_rect(aes(xmin = 100,ymin = 0, xmax = 200, ymax = max(df.SD$mean)),fill = "grey75") +
      geom_ribbon(aes(x = as.numeric(format(date,"%j")), ymin = mean - sd,ymax = mean + sd, group = as.factor(time_shift),
                      fill = as.factor(time_shift)),alpha = 0.25, show.legend = F) +
      geom_line(aes(x = as.numeric(format(date,"%j")), y = mean, group = as.factor(time_shift),
                    colour = as.factor(time_shift)),lwd = 1) +
      scale_fill_manual(values = cols) +
      scale_y_continuous("Mean population size") +
      scale_x_continuous("Day of the year",breaks = seq(0,360,20),guide = guide_axis(n.dodge = 2)) +
      scale_colour_manual(paste0("Shift application\nstart (days)"), values = cols) +
      #guides(linetype = guide_legend(title = "Temperature amplitude")) + 
      geom_vline(xintercept = c(100),linetype = 6,lwd = 0.5, show.legend = F) +
      # facet_wrap(.~envTamp,nrow = length(h)) + 
      ggtitle("DEB-GUTS-T") + 
      theme_pubr(legend = "right") +
      coord_cartesian(ylim = c(0,1100),expand = T)
    
    # Let's also look at relative differences
    df.SD <- lapply(unique(df.SD$time_shift), function(x){
      out <- df.SD[df.SD$time_shift == x,]
      out$mean <- round((out$mean / df.SD$mean[df.SD$time_shift==0])*100,digits = 1)
      out}) %>% do.call(rbind,.)
    # Let's also look at relative differences
    df.SDT <- lapply(unique(df.SDT$time_shift), function(x){
      out <- df.SDT[df.SDT$time_shift == x,]
      out$mean <- round((out$mean / df.SDT$mean[df.SDT$time_shift==0])*100,digits = 1)
      out}) %>% do.call(rbind,.)
    
    p3 <- ggplot(df.SD) +
      # geom_rect(aes(xmin = 100,ymin = 0, xmax = 200, ymax = max(df.SD$mean)),fill = "grey75") +
      geom_line(aes(x = as.numeric(format(date,"%j")), y = mean, group = as.factor(time_shift),
                    colour = as.factor(time_shift)),lwd = 1) +
      scale_colour_manual(values = cols) +
      scale_x_continuous("Day of the year",breaks = seq(0,360,20),guide = guide_axis(n.dodge = 2)) +
      guides(colour = "none",linetype = "none") + 
      ylab("Relative population size (%)") +
      geom_vline(xintercept = c(100),linetype = 6,lwd = 0.5, show.legend = F) + 
      ggtitle("DEB-GUTS") + 
      coord_cartesian(ylim = c(0,500),expand = T) +
      theme_pubr() 
    
    p4 <- ggplot(df.SDT) +
      # geom_rect(aes(xmin = 100,ymin = 0, xmax = 200, ymax = max(df.SD$mean)),fill = "grey75") +
      geom_line(aes(x = as.numeric(format(date,"%j")), y = mean, group = as.factor(time_shift),
                    colour = as.factor(time_shift)),lwd = 1) +
      scale_y_continuous("Relative population size (%)") +
      scale_x_continuous("Day of the year",breaks = seq(0,360,20),guide = guide_axis(n.dodge = 2)) +
      scale_colour_manual(paste0("Shift application\nstart (days)"), values = cols) +
      #guides(linetype = guide_legend(title = "Temperature amplitude")) + 
      geom_vline(xintercept = c(100),linetype = 6,lwd = 0.5, show.legend = F) +
      ggtitle("DEB-GUTS-T") + 
      theme_pubr(legend = "right") +
      coord_cartesian(ylim = c(0,500),expand = T)
    
    list(SD = p1, SDT = p2, SD_Rel = p3, SDT_Rel = p4)}
  else{
    ################################################################
    ################# From here is the original code ##########################
    # full data set
    df <- df_SD.list
    
    # select all scenarios matching concentrations and T-amplitudes ->> This should go into a seperate function that is called here and in other functions, to avoid correcting it in multiple places ic code needs changing
    scens <- as.matrix(df[[1]]$scenarios$V1)
    scens <- apply(scens,1, #the 1 as the second argument in apply function means the function is being applied to each row 
                   function(z) as.numeric(  #save as numeric value
                     apply(as.matrix(unlist(strsplit(z," ")) #splitting scenario string by spaces then unlisting them into a vector
                                     [c(8)]),1,     #selecting element(s) of the vector ##Needs adjustment based on input string
                           #here: 8 = exposure concentration
                           function(x) substr(x, start = 1, nchar(x) - 1))))  # removes the % from each line (don't know why yet..)
    scens <- data.frame(exposureConcentration = scens)
    
    # create boolean for scenarios
    if(!is.null(desired.exposure.concentrations)){
      select.scens <- scens$exposureConcentration %in% desired.exposure.concentrations
    scens <- scens[select.scens,]}else{select.scens <- rep(T,length(scens$exposureConcentration))}
    
    
    # select all dataframes matching the desired concentrations and T-amplitudes by matching multiple string patterns
    df.SD <- df[[1]]$data[select.scens]
    df.SDT <- df[[2]]$data[select.scens]
    
    # select only application window data and add scenario info
    ## Note: this seems trivial, but sadly isn't....
    ## Reason: simulations appear to use actual date information, which means we have to account for leap-years
    ## this can only be done by making a time-series of actual dates, for which we need the run control info
    ### Here we already have the code that gets the infos from the csv!!!
    startYear <- df[[1]]$run.info$V2[df[[1]]$run.info$V1=="startYear:"]
    endYear <- df[[1]]$run.info$V2[df[[1]]$run.info$V1=="endYear:"]
    startApplicationYear <- df[[1]]$run.info$V2[df[[1]]$run.info$V1=="startApplicationYear:"]
    endApplicationYear <- df[[1]]$run.info$V2[df[[1]]$run.info$V1=="endApplicationYear:"]
    
    days.warmup <- seq(as.Date(paste(startYear,"/1/1",sep = "")),as.Date(paste(startApplicationYear-1,"/12/31",sep = "")),"days")
    days.application <- seq(as.Date(paste(startApplicationYear,"/1/1",sep = "")),as.Date(paste(endApplicationYear,"/12/31",sep = "")),"days")
    days.recovery <- seq(as.Date(paste(endApplicationYear + 1,"/1/1",sep = "")),as.Date(paste(endYear,"/12/31",sep = "")),"days")
    days.full.period <- seq(as.Date(paste(startYear,"/1/1",sep = "")),as.Date(paste(endYear,"/12/31",sep = "")),"days")
    
    df.SD <- lapply(1:length(scens),function(x){
      out <- df.SD[[x]]
      if(!nrow(out)==1){
        out$date <- days.full.period
        out <- out[out$date %in% days.application,]
        out$exposureConc <- scens[x]
      }
      else{
        out$date <- NA
        out$exposureConc <- scens[x]
      }
      out}) %>% do.call(rbind,.)
    df.SDT <- lapply(1:length(scens),function(x){
      out <- df.SDT[[x]]
      if(!nrow(out)==1){
        out$date <- days.full.period
        out <- out[out$date %in% days.application,]
        out$exposureConc <- scens[x]
      }
      else{
        out$date <- NA
        out$exposureConc <- scens[x]
      }
      out}) %>% do.call(rbind,.)
    
    h <- time.range
    
    df.SD <- df.SD[format(df.SD$date,"%Y") %in% h,]
    df.SD$year <- format(df.SD$date,"%Y")
    
    df.SDT <- df.SDT[format(df.SDT$date,"%Y") %in% h,]
    df.SDT$year <- format(df.SDT$date,"%Y")
    
    # # data frame with temperature data
    # Ts <- data.frame(x = 1:365, apply(as.matrix(unique(df.SD$envTamp)),1,function(x){temps(x,365)}))
    # Ts <- pivot_longer(Ts,cols = -1,names_to = "grp",values_to = "Temperature")
    # Ts <- left_join(Ts,data.frame(grp = paste0("X",1:length(unique(df.SD$envTamp))),envTamp = unique(df.SD$envTamp)))
    
    # create a colour scale for the exposure concentrations
    cols <- grDevices::colorRampPalette(colors = c("forestgreen","red"))
    ncols <- ifelse(is.null(desired.exposure.concentrations),length(select.scens),length(desired.exposure.concentrations))
    cols <- cols(ncols)
    
    if(exposure.type == "pulsed"){
      p1 <- ggplot(df.SD) +
        geom_ribbon(aes(x = as.numeric(format(date,"%j")), ymin = mean - sd,ymax = mean + sd, group = as.factor(exposureConc),
                        fill = as.factor(exposureConc)),alpha = 0.25, show.legend = F) +
        # geom_rect(aes(xmin = 100,ymin = 0, xmax = 200, ymax = max(df.SD$mean)),fill = "grey75") +
        geom_line(aes(x = as.numeric(format(date,"%j")), y = mean, group = as.factor(exposureConc),
                      colour = as.factor(exposureConc)),lwd = 1) +
        scale_fill_manual(values = cols) +
        scale_colour_manual(values = cols) +
        guides(colour = "none",fill = "none",linetype = "none") + 
        xlab("Day of the year") + ylab("Mean population size") +
        geom_vline(xintercept = c(100,120,140,160,180,200),linetype = 6,lwd = 0.5, show.legend = F) +
        ggtitle("DEB-GUTS") + 
        coord_cartesian(ylim = c(0,max(df.SD$mean,df.SDT$mean)),expand = T) +
        theme_pubr() 
      
      p2 <- ggplot(df.SDT) +
        geom_ribbon(aes(x = as.numeric(format(date,"%j")), ymin = mean - sd,ymax = mean + sd, group = as.factor(exposureConc),
                        fill = as.factor(exposureConc)),alpha = 0.25, show.legend = F) +
        # geom_rect(aes(xmin = 100,ymin = 0, xmax = 200, ymax = max(df.SD$mean)),fill = "grey75") +
        geom_line(aes(x = as.numeric(format(date,"%j")), y = mean, group = as.factor(exposureConc),
                      colour = as.factor(exposureConc)),lwd = 1) +
        scale_y_continuous("Mean population size") +
        scale_x_continuous("Day of the year") +
        scale_fill_manual(values = cols) +
        scale_colour_manual(paste0("Exposure\nconcentration (","\u00B5","g/L)"), values = cols) +
        ggtitle("DEB-GUTS-T") + 
        theme_pubr(legend = "right") +
        coord_cartesian(ylim = c(0,max(df.SD$mean,df.SDT$mean)*1.1),expand = T) +
        geom_vline(xintercept = c(100,120,140,160,180,200),linetype = 6,lwd = 0.5, show.legend = F)
    }
    if(exposure.type=="constant"){
      if(y.trans == "log10"){
        p1 <- ggplot(df.SD) + #Plotting the population dynamics of the DEB-GUTS version
          geom_ribbon(aes(x = date, ymin = mean - sd,ymax = mean + sd, group = as.factor(exposureConc),
                          fill = as.factor(exposureConc)),alpha = 0.25, show.legend = F) +  #color the area of uncertainty
          geom_line(aes(x = date, y = mean, group = as.factor(exposureConc),
                        colour = as.factor(exposureConc)),lwd = 1) +                        #line with the mean of the simulations
          scale_x_date(date_breaks = "6 months") +
          scale_fill_manual(values = cols) +
          scale_colour_manual(values = cols) +
          scale_y_continuous(trans = "log10") +
          guides(colour = "none",fill = "none",linetype = "none") + 
          xlab("Day of the year") + ylab("Mean population size") +
          ggtitle("DEB-GUTS") + 
          coord_cartesian(ylim = c(min(df.SD$mean,df.SDT$mean), max(df.SD$mean,df.SDT$mean)*1.1),expand = T) +
          theme_pubr(x.text.angle = 45)
        
        p2 <- ggplot(df.SDT) + #Plotting the population dynamics of the DEB-GUTS-T version
          geom_ribbon(aes(x = date, ymin = mean - sd,ymax = mean + sd, group = as.factor(exposureConc),
                          fill = as.factor(exposureConc)),alpha = 0.25, show.legend = F) +  #color the area of uncertainty
          geom_line(aes(x = date, y = mean, group = as.factor(exposureConc),
                        colour = as.factor(exposureConc)),lwd = 1) +                        #line with the mean of the simulations
          scale_x_date(date_breaks = "6 months") +
          scale_fill_manual(values = cols) +
          scale_colour_manual(values = cols) +
          scale_y_continuous(trans = "log10")+
          guides(colour = "none",fill = "none",linetype = "none") + 
          xlab("Day of the year") + ylab("Mean population size") +
          ggtitle("DEB-GUTS-T") + 
          coord_cartesian(ylim = c(min(df.SD$mean,df.SDT$mean), max(df.SD$mean,df.SDT$mean)*1.1),expand = T) +
          theme_pubr(x.text.angle = 45)
      }
      else{
        p1 <- ggplot(df.SD) + #Plotting the population dynamics of the DEB-GUTS version
          geom_ribbon(aes(x = date, ymin = mean - sd,ymax = mean + sd, group = as.factor(exposureConc),
                          fill = as.factor(exposureConc)),alpha = 0.25, show.legend = T) +  #color the area of uncertainty
          geom_line(aes(x = date, y = mean, group = as.factor(exposureConc),
                        colour = as.factor(exposureConc)),lwd = 1) +                        #line with the mean of the simulations
          scale_x_date(date_breaks = "6 months") +
          scale_fill_manual(values = cols) +
          scale_colour_manual(values = cols) +
          # scale_y_continuous(sec.axis = sec_axis(~ . / (max(df.SD$mean,df.SDT$mean) / max(df.SD$envT)), name = "Temperature (Celsius)")) + #second y axis for temperature profile
          guides(colour = "none",fill = "none",linetype = "none") + 
          xlab("Day of the year") + ylab("Mean population size") +
          ggtitle("DEB-GUTS") + 
          coord_cartesian(ylim = c(0, max(df.SD$mean,df.SDT$mean)*1.1),expand = T) +
          theme_pubr(x.text.angle = 45)
        
        p2 <- ggplot(df.SDT) + #Plotting the population dynamics of the DEB-GUTS-T version
          geom_ribbon(aes(x = date, ymin = mean - sd,ymax = mean + sd, group = as.factor(exposureConc),
                          fill = as.factor(exposureConc)),alpha = 0.25, show.legend = T) +  #color the area of uncertainty
          geom_line(aes(x = date, y = mean, group = as.factor(exposureConc),
                        colour = as.factor(exposureConc)),lwd = 1) +                        #line with the mean of the simulations
          scale_x_date(date_breaks = "6 months") +
          scale_fill_manual(values = cols) +
          scale_colour_manual(values = cols) +
          # scale_y_continuous(sec.axis = sec_axis(~ . / (max(df.SD$mean,df.SDT$mean) / max(df.SD$envT)), name = "Temperature (Celsius)")) + #second y axis for temperature profile
          #guides(colour = "none",fill = "none",linetype = "none") + 
          xlab("Day of the year") + ylab("Mean population size") +
          ggtitle("DEB-GUTS-T") + 
          theme(legend.position = "right") +
          coord_cartesian(ylim = c(0, max(df.SD$mean,df.SDT$mean)*1.1),expand = T) +
          theme_pubr(x.text.angle = 45)

      }
    }
    # Plot for environmental temperature
    p3 <- ggplot(df.SD) +
      geom_line(aes(x = date, y = envT), colour = "black", lwd = 1) +
      scale_x_date(date_breaks = "6 months",  date_labels = "") +
      xlab("") +  # Remove x-axis label for the top plot
      ylab("Temperature (Celsius)") +
      theme(panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text.x=element_blank())+
      ggtitle("Environmental Temperature")
    
    
    list(SD = p1, SDT = p2, p_envT= p3)  
  }
}

# New function for figures that compare the model simulations of GUTS,our current GUTS-T and the GUTS-T-Std   
plotModelComparison <- function(df_SD.list,
                                desired.exposure.concentrations,
                                time.range
                                ){
  
  ## Prepare df 
  # full data set
  df <- df_SD.list

  ## Get exposure concentrations 
  # select all scenarios matching concentrations  ->> This should go into a seperate function that is called here and in other functions, to avoid correcting it in multiple places ic code needs changing
  scens <- as.matrix(df[[1]]$scenarios$V1)
  scens <- apply(scens,1, #the 1 as the second argument in apply function means the function is being applied to each row 
                 function(z) as.numeric(  #save as numeric value
                   apply(as.matrix(unlist(strsplit(z," ")) #splitting scenario string by spaces then unlisting them into a vector
                                   [c(8)]),1,     #selecting element(s) of the vector ##Needs adjustment based on input string
                         #here: 8 = exposure concentration
                         function(x) substr(x, start = 1, nchar(x) - 1))))  # removes the % from each line (don't know why yet..)
  scens <- data.frame(exposureConcentration = scens)
  
  
  # create boolean to filter for scenarios
  if(!is.null(desired.exposure.concentrations)){
    select.scens <- scens$exposureConcentration %in% desired.exposure.concentrations
    scens <- scens[select.scens, , drop = FALSE]}else{select.scens <- rep(T,length(scens$exposureConcentration))}
  
  # select all dataframes matching the desired concentrations by matching multiple string patterns
  # Filter for exposure scenarios
  df.SD <- df[[1]]$data[select.scens]
  df.SDT <- df[[2]]$data[select.scens]
  df.SDTStd <- df[[3]]$data[select.scens]
  
  # # Get model types in df.list
  # df_temp <- lapply(1:length(df),function(x) df[[x]]$data[select.scens][[1]])
  # types <- lapply(df_temp, function(x) unique(x$model.version)) %>% do.call(rbind,.)
  
  # select only application window data and add scenario info  
  ## Note: this seems trivial, but sadly isn't....
  ## Reason: simulations appear to use actual date information, which means we have to account for leap-years
  ## this can only be done by making a time-series of actual dates, for which we need the run control info
  ### Here we already have the code that gets the infos from the csv!!!
  startYear <- df[[1]]$run.info$V2[df[[1]]$run.info$V1=="startYear:"]
  endYear <- df[[1]]$run.info$V2[df[[1]]$run.info$V1=="endYear:"]
  startApplicationYear <- df[[1]]$run.info$V2[df[[1]]$run.info$V1=="startApplicationYear:"]
  endApplicationYear <- df[[1]]$run.info$V2[df[[1]]$run.info$V1=="endApplicationYear:"]
  
  days.warmup <- seq(as.Date(paste(startYear,"/1/1",sep = "")),as.Date(paste(startApplicationYear-1,"/12/31",sep = "")),"days")
  days.application <- seq(as.Date(paste(startApplicationYear,"/1/1",sep = "")),as.Date(paste(endApplicationYear,"/12/31",sep = "")),"days")
  days.recovery <- seq(as.Date(paste(endApplicationYear + 1,"/1/1",sep = "")),as.Date(paste(endYear,"/12/31",sep = "")),"days")
  days.full.period <- seq(as.Date(paste(startYear,"/1/1",sep = "")),as.Date(paste(endYear,"/12/31",sep = "")),"days")
  # apply to each model type
  df.SD <- lapply(1:nrow(scens),function(x){
    out <- df.SD[[x]]
    if(!nrow(out)==1){
      out$date <- days.full.period
      out <- out[out$date %in% days.application,]
      out$exposureConc <- scens$exposureConcentration[x]
    }
    else{
      out$date <- NA
      out$exposureConc <- scens$exposureConcentration[x]
    }
    out}) %>% do.call(rbind,.)
  
  df.SDT <- lapply(1:nrow(scens),function(x){
    out <- df.SDT[[x]]
    if(!nrow(out)==1){
      out$date <- days.full.period
      out <- out[out$date %in% days.application,]
      out$exposureConc <- scens$exposureConcentration[x]
    }
    else{
      out$date <- NA
      out$exposureConc <- scens$exposureConcentration[x]
    }
    out}) %>% do.call(rbind,.)
  
  df.SDTStd <- lapply(1:nrow(scens),function(x){
    out <- df.SDTStd[[x]]
    if(!nrow(out)==1){
      out$date <- days.full.period
      out <- out[out$date %in% days.application,]
      out$exposureConc <- scens$exposureConcentration[x]
    }
    else{
      out$date <- NA
      out$exposureConc <- scens$exposureConcentration[x]
    }
    out}) %>% do.call(rbind,.)
  
  h <- time.range
  
  df.SD <- df.SD[format(df.SD$date,"%Y") %in% h,]
  df.SD$year <- format(df.SD$date,"%Y")
  
  df.SDT <- df.SDT[format(df.SDT$date,"%Y") %in% h,]
  df.SDT$year <- format(df.SDT$date,"%Y")
  
  df.SDTStd <- df.SDTStd[format(df.SDTStd$date,"%Y") %in% h,]
  df.SDTStd$year <- format(df.SDTStd$date,"%Y")
  
  # Plotting ############################################################################################ adjust to get plot as in paper draft 
  # Combine all data into a single data frame
  combined_df <- bind_rows(df.SD, df.SDT, df.SDTStd)
  
  # Create a dynamic title
  model_versions <- unique(combined_df$model.version)
  esposure_chemical <- unique(combined_df$exposure.chemical)
  dynamic_title <- paste("Model Comparisons:", paste(model_versions, collapse = ", "), "for", paste(esposure_chemical, collapse = ", "), "exposure")
  
  
  # Plot the data for all model types per exposure concentration
  p1 <- ggplot(combined_df, aes(x = date, y = mean, color = model.version)) +
    geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd, fill = model.version), #SD shade behind line
                alpha = 0.25, color = NA, show.legend = F) +
    geom_line() +
    facet_wrap(~ exposureConc, ncol = 1) + # Facet by exposure concentration
    labs(title = dynamic_title,
         x = "Date",
         y = "Mean",
         color = "Model Version"
        ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "bottom"
        )
  #list(p1)
} 


plotScenarioComparison <- function(df_D3ref.list, 
                                   df_Hn2150.list,
                                   relative_diff,
                                   desired.exposure.concentration,
                                   time.range){ 
  
  ## Prepare combined df 
  df_D3ref <- process_model_data(df_D3ref.list, desired.exposure.concentration, time.range)
  df_Hn2150 <- process_model_data(df_Hn2150.list, desired.exposure.concentration, time.range)
  
  #Also relative differences, easier to do here instead of after the combining
  df_D3ref$SD$mean_rel      <- df_D3ref$SD$mean      / df_D3ref$SD$mean
  df_D3ref$SD$sd_rel        <- df_D3ref$SD$sd        / df_D3ref$SD$sd # not sure if this correct!

  df_D3ref$SDT$mean_rel     <- df_D3ref$SDT$mean     / df_D3ref$SD$mean
  df_D3ref$SDT$sd_rel       <- df_D3ref$SDT$sd       / df_D3ref$SD$sd

  df_D3ref$SDTStd$mean_rel  <- df_D3ref$SDTStd$mean  / df_D3ref$SD$mean
  df_D3ref$SDTStd$sd_rel    <- df_D3ref$SDTStd$sd    / df_D3ref$SD$sd

  df_Hn2150$SD$mean_rel     <- df_Hn2150$SD$mean     / df_Hn2150$SD$mean
  df_Hn2150$SD$sd_rel       <- df_Hn2150$SD$sd       / df_Hn2150$SD$sd # not sure if this correct!

  df_Hn2150$SDT$mean_rel    <- df_Hn2150$SDT$mean    / df_Hn2150$SD$mean
  df_Hn2150$SDT$sd_rel      <- df_Hn2150$SDT$sd      / df_Hn2150$SD$sd

  df_Hn2150$SDTStd$mean_rel <- df_Hn2150$SDTStd$mean / df_Hn2150$SD$mean
  df_Hn2150$SDTStd$sd_rel   <- df_Hn2150$SDTStd$sd   / df_Hn2150$SD$sd
  
  # Combine df for plotting
  combined_df <- bind_rows(df_D3ref$SD,     df_Hn2150$SD,
                           df_D3ref$SDT,    df_Hn2150$SDT,
                           df_D3ref$SDTStd, df_Hn2150$SDTStd)

  
  
  
  # Plotting ############################################################################################ 
  # Create a dynamic title
  model_versions <- unique(combined_df$model.version)
  esposure_chemical <- unique(combined_df$exposure.chemical)
  T_scenarios <- unique(combined_df$T.scenario)
  dynamic_title <- paste("Model Comparisons:", paste(model_versions, collapse = ", "), "for", paste(esposure_chemical, collapse = ", "), "exposure")
  
  
  # Plot the data for all model types per exposure concentration
  if (relative_diff == F) { #Plotting absolute values 
    p1 <- ggplot(combined_df, aes(x = date, y = mean, color = model.version)) +
      geom_line() +
      geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd, fill = model.version), #SD shade behind line
                  alpha = 0.25, color = NA, show.legend = F) +
      facet_grid(exposureConc ~ factor(T.scenario, levels = T_scenarios)) + # Facet by exposure concentration and T.scenario (the factor argument is just to ensure that D3 scenario is left of plot)
      labs(title = dynamic_title,
           x = "Date",
           y = "Mean",
           color = "Model Version"
      ) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "bottom",
            panel.spacing = unit(1, "lines")
      )
  }else{ #Plotting relative to control 
    p1 <- ggplot(combined_df, aes(x = date, y = mean_rel, color = model.version)) +
      geom_line() +
      geom_ribbon(aes(ymin = mean_rel - sd_rel, ymax = mean_rel + sd_rel, fill = model.version), #SD shade behind line
                  alpha = 0.25, color = NA, show.legend = F) +
      facet_grid(exposureConc ~ factor(T.scenario, levels = T_scenarios), scales = "free_y") + # Facet by exposure concentration and T.scenario (the factor argument is just to ensure that D3 scenario is left of plot)
      geom_hline(yintercept = 1, linetype = "dashed", color = "gray") +  # Horizontal line at y = 1 (control baseline)      
      labs(title = dynamic_title,
           x = "Date",
           y = "Relative difference",
           color = "Model Version"
      ) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "bottom",
            panel.spacing = unit(1, "lines")
      )
  }
  
}


plotPopQuantilesTvsNoT <- function(popsize.data.frame, 
                                   desired.exposure.concentrations,
                                   label.txt,
                                   model.type){
  df <- popsize.data.frame[popsize.data.frame$exposureConcentration %in% desired.exposure.concentrations,]
  df <- df[order(df$exposureConcentration,decreasing = T),]
  
  columns <- names(df)[grepl(pattern = paste(paste0("Q",c(10,50,90),"$"),collapse = "|"),x =   names(df))]
  
  # create a colour scale for the exposure concentrations
  cols <- grDevices::colorRampPalette(colors = c("forestgreen","red"))
  cols <- cols(length(desired.exposure.concentrations))
  
  p1 <- ggplot(data = df, aes(x = .data[[columns[4]]], y = .data[[columns[1]]], 
                              colour = as.factor(exposureConcentration), size = tempAmplitude)) + 
    geom_point() + guides(size = "none", colour = "none") +
    scale_size_continuous(breaks = unique(df$tempAmplitude)) +
    scale_colour_manual("Temperature amplitude",values = cols) +
    geom_line(data = data.frame(x = c(0,max(c(df[,columns[1]],df[,columns[4]]),na.rm = T)), 
                                y = c(0,max(c(df[,columns[1]],df[,columns[4]]),na.rm = T))),
              aes(x = x, y = y), colour = "black", inherit.aes = F) +
    xlab(paste0("Q10 mean population\nsize for ",model.type,"-T")) +
    ylab(paste0("Q10 mean population\nsize for ",model.type)) + 
    coord_equal()+ theme_pubr() + 
    geom_text(data = data.frame(x = -Inf, y = Inf, label = label.txt[1]),aes(x = -Inf, y = Inf, label = label),
              hjust = -0.5,
              vjust = 1.5, inherit.aes = F)
  
  p2 <- ggplot(data = df, aes(x = .data[[columns[5]]], y = .data[[columns[2]]], 
                              colour = as.factor(exposureConcentration), size = tempAmplitude)) + 
    geom_point() + guides(size = "none", colour = "none") +
    scale_size_continuous(breaks = unique(df$tempAmplitude)) +
    scale_colour_manual("Exposure concentration",values = cols) +
    geom_line(data = data.frame(x = c(0,max(c(df[,columns[2]],df[,columns[5]]),na.rm = T)), 
                                y = c(0,max(c(df[,columns[2]],df[,columns[5]]),na.rm = T))),
              aes(x = x, y = y), colour = "black", inherit.aes = F) +
    xlab(paste0("Q50 mean population\nsize for ",model.type,"-T")) +
    ylab(paste0("Q50 mean population\nsize for ",model.type)) +
    coord_equal()+ theme_pubr()+
    geom_text(data = data.frame(x = -Inf, y = Inf, label = label.txt[2]),aes(x = -Inf, y = Inf, label = label),
              hjust = -0.5,
              vjust = 1.5, inherit.aes = F)
  
  p3 <- ggplot(data = df, aes(x = .data[[columns[6]]], y = .data[[columns[3]]], 
                              colour = as.factor(exposureConcentration), size = tempAmplitude)) + 
    geom_point() + guides(size = guide_legend(title = "Temparature\namplitude (Celsius)")) + 
    scale_size_continuous(breaks = unique(df$tempAmplitude)) +
    scale_colour_manual(paste0("Exposure\nconcentration (","\u00B5","g/L)"),values= cols) +
    geom_line(data = data.frame(x = c(0,max(c(df[,columns[3]],df[,columns[6]]),na.rm = T)), 
                                y = c(0,max(c(df[,columns[3]],df[,columns[6]]),na.rm = T))),
              aes(x = x, y = y), colour = "black", inherit.aes = F) +
    xlab(paste0("Q90 mean population\nsize for ",model.type,"-T")) +
    ylab(paste0("Q90 mean population\nsize for ",model.type)) +
    coord_equal() + theme_pubr(legend =  "right") + theme(legend.justification = "top") +
    geom_text(data = data.frame(x = -Inf, y = Inf, label = label.txt[3]),aes(x = -Inf, y = Inf, label = label),
              hjust = -0.5,
              vjust = 1.5, inherit.aes = F)
  
  # create cumulative distribution plots
  df <- pivot_longer(data = df[,-c(1,102)],cols = -c(exposureConcentration,tempAmplitude),names_to = "Q",values_to = "AvgPopSize")
  df$m.type <- do.call(rbind,lapply(strsplit(df$Q,"_"),function(x) x[1]))[,1]
  df$x <- as.numeric(do.call(rbind,lapply(strsplit(df$Q,"Q"),function(x) x[2]))[,1])
  df$grp <- paste0(df$exposureConcentration,"_",df$tempAmplitude)
  df$exposureConcentration <- paste0(df$exposureConcentration," (","\u00B5","g/L)")
  
  p4 <- ggplot(data = df[df$m.type == "SD" | df$m.type == "IT" ,], 
               aes(x = AvgPopSize, y = x, colour = tempAmplitude, group = grp)) + 
    geom_line() + 
    scale_y_continuous(breaks = seq(0,100,10)) +
    scale_colour_gradient("Temperature\namplitude (Celsius)",low = "blue",high = "red") + 
    theme_pubr(legend = "right") + guides(colour = "none", x = guide_axis(angle = 90)) + 
    ylab("Cumulative frequency") + xlab("Mean population size") +
    geom_hline(yintercept = c(10,50,90)) +
    facet_wrap(.~exposureConcentration) + ggtitle("DEB-GUTS")
  
  p5 <- ggplot(data = df[df$m.type == "SDT" | df$m.type == "ITT" ,], 
               aes(x = AvgPopSize, y = x, colour = tempAmplitude, group = grp)) + 
    geom_line() +
    scale_y_continuous(breaks = seq(0,100,10)) +
    scale_colour_gradient("Temperature\namplitude (Celsius)",low = "blue",high = "red") + 
    theme_pubr(legend = "right") + guides(x = guide_axis(angle = 90)) +
    ylab("Cumulative frequency") + xlab("Mean population size") +
    geom_hline(yintercept = c(10,50,90)) +
    facet_wrap(.~exposureConcentration) + ggtitle("DEB-GUTS-T")
  
  list(Quantile.plots = p1 + p2 + p3, Cumulative.plots = p4 + p5)
  
}


