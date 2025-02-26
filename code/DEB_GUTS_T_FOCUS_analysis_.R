source("./code/DEB_GUTS_T_FOCUS_functions.R")
pks.installed <- checkPackages(library.loc = "C:/ProgramData/R/R_4.4.1") # Adjust to your settings
loadPackages(library.loc = "C:/ProgramData/R/R_4.4.1", # Adjust to your settings
             required.packages.installed = pks.installed)


## constant exposures ###################

# Choose T profiles #########################################
T_D3ref = "./data/runsOct2024_D3ref"   #Adjust to file location 
# Read data #################################################
df_SD_D3 <- readData(data.location = T_D3ref , #Adjust to  temperature scenario
                  guts.model.version = "SD", 
                  ignore.version = "pulsed", 
                  ignore.chemical = "FPF",
                  filter.concentrations = F,
                  filter.temp.amplitudes = F,
                  desired.exposure.concentration = F,
                  desired.temp.amplitude = F,
                  application.pulse.shift = F)

df_IT_D3 <- readData(data.location = T_D3ref , #Adjust to  temperature scenario
                  guts.model.version = "IT", 
                  ignore.version = "pulsed",
                  ignore.chemical = "FPF",
                  filter.concentrations = F,
                  filter.temp.amplitudes = F,
                  desired.exposure.concentration = F,
                  desired.temp.amplitude = 1, 
                  application.pulse.shift = F)

# Choose T profiles #########################################
T_Hn2150 = "./data/runsOct2024_Hn_2150" #Adjust to file location
# Read data #################################################
df_SD_Hn <- readData(data.location = T_Hn2150 , 
                     guts.model.version = "SD", 
                     ignore.version = "pulsed", 
                     ignore.chemical = "FPF",
                     filter.concentrations = F,
                     filter.temp.amplitudes = F,
                     desired.exposure.concentration = F,
                     desired.temp.amplitude = F,
                     application.pulse.shift = F)

df_IT_Hn <- readData(data.location = T_Hn2150 , 
                     guts.model.version = "IT", 
                     ignore.version = "pulsed",
                     ignore.chemical = "FPF",
                     filter.concentrations = F,
                     filter.temp.amplitudes = F,
                     desired.exposure.concentration = F,
                     desired.temp.amplitude = 1, 
                     application.pulse.shift = F)


## Plot of population dynamics in relation to temp amplitude and selected years
#SD-constant 
p1 <- plotTAmpPopDynamics(df_SD.list = df_SD_D3, 
                          exposure.type = "constant",
                          desired.exposure.concentrations = c("0","0.4","0.8","1.2"),
                          desired.temp.amplitudes = F,
                          time.range = c(1989:2016),
                          y.trans = F   #"log10" or F ,  # switch to log transform y axis
                          )

p_SD <- p1$p_envT/p1$SD + plot_layout(heights = c(1, 2)) 
p_SDT <- p1$p_envT/p1$SDT + plot_layout(heights = c(1, 2)) 

(p_SD | p_SDT ) + plot_annotation(title = "SD-IMI") 

#################### Plotting new
#SD-constant
p1 <- plotModelComparison(df_SD.list = df_SD_D3, 
                          desired.exposure.concentrations =  c("0","0.4","0.8","1.2"), # when including control, don't use 0.0 but only 0, else it doesn't match 
                          time.range = c(1989:2016)
                          )
p1


p1 <- plotScenarioComparison(df_D3ref.list = df_SD_D3, 
                             df_Hn2150.list = df_SD_Hn,
                             relative_diff = F, # F for absolute mean values, T for relative differences to control
                             desired.exposure.concentration = c("0","0.4","0.8","1.2"), # when including control, don't use 0.0 but only 0, else it doesn't match 
                             compare_models = F, # F for comparison within model version, i.e., all concentrations in one plot 
                             time.range = c(1989:2016))
p1

############# Plotting parameter values IMI
## Prepare combined df 
df_D3ref_SD <- process_model_data(df.list = df_SD_D3,
                                     desired.exposure.concentration =  c("0","0.4","0.8","1.2"),
                                     time.range = c(1989:2016))


# Plot the parameters SDT
p_SDT <- plot_temp_corrected_parameters(
  df = df_D3ref_SD,
  T_A = 27240,
  kd = 0.01,
  mw = 1.045e-16 ,
  bw = 0.098 ,
  method = "SDT"  
)

ggarrange(plotlist = p_SDT, ncol = 2, nrow = 2) 

# Plot the parameters SDTstd
p_SDTStd <- plot_temp_corrected_parameters(
  df = df_D3ref_SD,
  T_A = 100,
  kd = 0.030, 
  mw = 1.483e-13 ,
  bw = 0.004 , 
  method = "SDTstd" 
)

ggarrange(plotlist = p_SDTStd, ncol = 2, nrow = 2) 

####### Prepare combined df for IT version
df_D3ref_IT <- process_model_data(df.list = df_IT_D3,
                                     desired.exposure.concentration =  c("0","0.4","0.8","1.2"),
                                     time.range = c(1989:2016))
# Plot the parameters ITT
p_ITT <- plot_temp_corrected_parameters(
  df = df_D3ref_IT,
  T_A = 1919,
  kd = 0.01, 
  mw = 0.620,
  bw = NA,
  method = "ITT"  
)

ggarrange(plotlist = p_ITT, ncol = 2, nrow = 2) 

# Plot the parameters ITTstd
p_ITTStd <- plot_temp_corrected_parameters(
  df = df_D3ref_IT,
  T_A = 2529,
  kd = 0.01, 
  mw = 0.596,
  bw = NA, 
  method = "ITTStd" 
)

ggarrange(plotlist = p_ITTStd, ncol = 2, nrow = 2) 


#############################################################################################


#IT-constant (NOTE: code elements still named SD but uses IT data now!, maybe this should be changed at somepoint to avoid confusion)
p2 <- plotTAmpPopDynamics(df_SD.list = df_IT,
                          exposure.type = "constant",
                          desired.exposure.concentrations = NULL,
                          desired.temp.amplitudes = F,
                          time.range = c(1989:2016),
                          y.trans = F   #"log10" or F ,  # switch to log transform y axis
                          )

p_SD <- p2$p_envT/p2$SD + plot_layout(heights = c(1, 2)) 
p_SDT <- p2$p_envT/p2$SDT + plot_layout(heights = c(1, 2)) #+ theme(legend.position="right") 

(p_SD | p_SDT) + plot_annotation(title = "IT-IMI") 



##############################################################
## constant exposures for FPF ###################


# Choose T profiles #########################################
T_D3ref = "C:/Users/magol001/OneDrive - Wageningen University & Research/Git_AMD/GammarusPulex_DEB-GUTS-T/data/runsOct2024_D3ref"   #Adjust to file location 
T_Hn2150 = "C:/Users/magol001/OneDrive - Wageningen University & Research/Git_AMD/GammarusPulex_DEB-GUTS-T/data/runsOct2024_Hn_2150" #Adjust to file location

# Read data #################################################
df_SD <- readData(data.location = T_D3ref , #Adjust to  temperature scenario
                  guts.model.version = "SD", 
                  ignore.version = "pulsed", 
                  ignore.chemical = "IMI",
                  filter.concentrations = F,
                  filter.temp.amplitudes = F,
                  desired.exposure.concentration = F,
                  desired.temp.amplitude = F,
                  application.pulse.shift = F)

df_IT <- readData(data.location = T_D3ref , #Adjust to  temperature scenario
                  guts.model.version = "IT", 
                  ignore.version = "pulsed",
                  ignore.chemical = "IMI",
                  filter.concentrations = F,
                  filter.temp.amplitudes = F,
                  desired.exposure.concentration = F,
                  desired.temp.amplitude = 1, 
                  application.pulse.shift = F)


## Plot of population dynamics in relation to temp amplitude and selected years
#SD-constant 
p3 <- plotTAmpPopDynamics(df_SD.list = df_SD, 
                          exposure.type = "constant",
                          desired.exposure.concentrations = NULL,
                          desired.temp.amplitudes = F,
                          time.range = c(1989:2016),
                          y.trans = F   #"log10" or F ,  # switch to log transform y axis
                          )

p_SD <- p3$p_envT/p3$SD + plot_layout(heights = c(1, 2)) 
p_SDT <- p3$p_envT/p3$SDT + plot_layout(heights = c(1, 2)) 

(p_SD | p_SDT) + plot_annotation(title = "SD-FPF") 


#IT-constant (NOTE: code elements still named SD but uses IT data now!, maybe this should be changed at somepoint to avoid confusion)
p4 <- plotTAmpPopDynamics(df_SD.list = df_IT,
                          exposure.type = "constant",
                          desired.exposure.concentrations = NULL,
                          desired.temp.amplitudes = F,
                          time.range = c(1989:2016),
                          y.trans = F   #"log10" or F ,  # switch to log transform y axis
                          )


p_SD <- p4$p_envT/p4$SD + plot_layout(heights = c(1, 2)) 
p_SDT <- p4$p_envT/p4$SDT + plot_layout(heights = c(1, 2)) #+ theme(legend.position="right") 

(p_SD | p_SDT) + plot_annotation(title = "IT-FPF") 


p3$SD + p4$SDT + plot_annotation(title = "SD vs ITT for FPF")
p3$SDT + p4$SD + plot_annotation(title = "SDT vs IT for FPF")


ref_plot <- plot(x = p1[["SD"]][["data"]][["date"]], y = p1[["SD"]][["data"]][["envT"]],
                    col = "black")

Hn2150_plot <- plot(x = p1[["SD"]][["data"]][["date"]], y = p1[["SD"]][["data"]][["envT"]],
              col = "black")


## Population size at end of the simulation ##################################################################################################################
# NEEDS ADJUSTMENTS TOO
df_PopsizeSD <- popSize(simulation.data.list = df_SD,application.pulse.shift = F)
df_PopsizeIT <- popSize(simulation.data.list = df_IT,application.pulse.shift = F)

# plot of quantiles of different scenarios
p3 <- plotPopQuantilesTvsNoT(popsize.data.frame = df_PopsizeSD,
                             desired.exposure.concentrations = c(0,0.1,0.2,0.4),
                             label.txt = c("A","B","C"),model.type = "SD")

p4 <- plotPopQuantilesTvsNoT(popsize.data.frame = df_PopsizeIT,desired.exposure.concentrations = c(0,0.1,0.2,0.4),
                                   label.txt = c("D","E","F"),model.type = "IT")

p3$Quantile.plots / p4$Quantile.plots + plot_layout(guides = "collect")

p3$Cumulative.plots
p4$Cumulative.plots




##############################################################
############### pulsed exposures #############################
##############################################################
## Read data #################################################
df_SD <- readData(data.location = "C:/...",
                  guts.model.version = "SD",ignore.version = "verylow", 
                  filter.concentrations = F,
                  filter.temp.amplitudes = F,
                  desired.exposure.concentration = F,
                  desired.temp.amplitude = 1,application.pulse.shift = F)

df_IT <- readData(data.location = "C:/...",
                  guts.model.version = "IT",ignore.version = "verylow", 
                  filter.concentrations = F,
                  filter.temp.amplitudes = F,
                  desired.exposure.concentration = F,
                  desired.temp.amplitude = 1,application.pulse.shift = F)


## Plot of population dynamics in relation to temp amplitude and selected years
#SD-constant first year
p1 <- plotTAmpPopDynamics(df_SD.list = df_SD,exposure.type = "pulsed",
                          desired.exposure.concentrations = c(0,10,20,30,40),desired.temp.amplitudes = c(1,5,10), time.range = 1992)

p1$SD + p1$SDT

#IT-constant first year (NOTE: code elements still named SD but uses IT data now!, maybe this should be changed at somepoint to avoid confusion)
p2 <- plotTAmpPopDynamics(df_SD.list = df_IT,exposure.type = "pulsed",
                          desired.exposure.concentrations = c(0,10,20,30,40),desired.temp.amplitudes = c(1,5,10), time.range = 1992)

p2$SD + p2$SDT


## Population size at end of the simulation ##################################################################################################################
df_PopsizeSD <- popSize(simulation.data.list = df_SD,application.pulse.shift = F)
df_PopsizeIT <- popSize(simulation.data.list = df_IT,application.pulse.shift = F)

# plot of quantiles of different scenarios
p3 <- plotPopQuantilesTvsNoT(popsize.data.frame = df_PopsizeSD,
                             desired.exposure.concentrations = c(0,10,20,30,40),
                             label.txt = c("A","B","C"),model.type = "SD")

p4 <- plotPopQuantilesTvsNoT(popsize.data.frame = df_PopsizeIT,
                                desired.exposure.concentrations = c(0,10,20,30,40),
                                label.txt = c("D","E","F"),model.type = "IT")

p3$Quantile.plots / p4$Quantile.plots + plot_layout(guides = "collect")

p3$Cumulative.plots
p4$Cumulative.plots

