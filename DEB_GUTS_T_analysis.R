pks.installed <- checkPackages(library.loc = "C:/ProgramData/R/...") # Adjust to your settings
loadPackages(library.loc = "C:/ProgramData/R/...", # Adjust to your settings
             required.packages.installed = pks.installed)

##############################################################
############### constant exposures ###########################
##############################################################
## Read data #################################################
df_SD <- readData(data.location = "C:/...", #Adjust to file location
                  guts.model.version = "SD", ignore.version = "pulsed", 
                  filter.concentrations = F,
                  filter.temp.amplitudes = F,
                  desired.exposure.concentration = F,
                  desired.temp.amplitude = 1,application.pulse.shift = F)

df_IT <- readData(data.location = "C:/...",
                  guts.model.version = "IT", ignore.version = "pulsed", 
                  filter.concentrations = F,
                  filter.temp.amplitudes = F,
                  desired.exposure.concentration = F,
                  desired.temp.amplitude = 1, application.pulse.shift = F)


## Plot of population dynamics in relation to temp amplitude and selected years
#SD-constant first year
p1 <- plotTAmpPopDynamics(df_SD.list = df_SD,exposure.type = "constant",
                          desired.exposure.concentrations = c(0,0.1,0.2,0.4),desired.temp.amplitudes = c(1,5,10),time.range = 1992)

p1$SD + p1$SDT

#IT-constant first year (NOTE: code elements still named SD but uses IT data now!, maybe this should be changed at somepoint to avoid confusion)
p2 <- plotTAmpPopDynamics(df_SD.list = df_IT,exposure.type = "constant",
                          desired.exposure.concentrations = c(0,0.1,0.2,0.4),desired.temp.amplitudes = c(1,5,10),time.range = 1992)

p2$SD + p2$SDT


## Population size at end of the simulation ##################################################################################################################
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


##############################################################
############### SD model Mod Kd constant exposure ############
##############################################################
## Read data #################################################
df_SD <- readData(data.location = "C:/.../modKD/",
                  guts.model.version = "SD",ignore.version = "pulsed", 
                  filter.concentrations = F,
                  filter.temp.amplitudes = F,
                  desired.exposure.concentration = 30,
                  desired.temp.amplitude = 1)

## Population size at end of the simulation ##################################################################################################################
df_Popsize <- popSize(simulation.data.list = df_SD)

# plot of population dynamics in relation to temp amplitude and selected years
p1 <- plotTAmpPopDynamics(df_SD.list = df_SD,desired.exposure.concentrations = c(0,0.1,0.2,0.4),
                          desired.temp.amplitudes = c(1,5,10),exposure.type = "constant",time.range = 1992)

p1$SD + p1$SDT

##############################################################
############### SD model Mod Kd pulsed exposure ##############
##############################################################
## Read data ##################################################
df_SD <- readData(data.location = "C:/.../modKD/",
                  guts.model.version = "SD",ignore.version = "verylow", 
                  filter.concentrations = F,
                  filter.temp.amplitudes = F,
                  desired.exposure.concentration = 30,
                  desired.temp.amplitude = 1)

## Population size at end of the simulation ##################################################################################################################
df_Popsize <- popSize(simulation.data.list = df_SD)

# plot of population dynamics in relation to temp amplitude and selected years
p1 <- plotTAmpPopDynamics(df_SD.list = df_SD,desired.exposure.concentrations = c(0,10,20,30,40),
                          desired.temp.amplitudes = c(1,5,10),exposure.type = "pulsed",time.range = 1992)

p1$SD + p1$SDT


################################################################
############### Pulse shifted data sD model ####################
################################################################
## Read data ###################################################
df_SD <- readData(data.location = "C:/.../pulse_shift",
                  guts.model.version = "SD",
                  ignore.version = "verylow",
                  filter.concentrations = F,
                  filter.temp.amplitudes = F,
                  desired.exposure.concentration = 20,
                  desired.temp.amplitude = c(1,10),
                  application.pulse.shift = T)

# plot of population dynamics in relation to temp amplitude and selected years
p1 <- plotTAmpPopDynamics(df_SD.list = df_SD,desired.exposure.concentrations = c(20),
                          desired.temp.amplitudes = c(10),exposure.type = "pulsed",time.range = 1992,application.pulse.shift = T)

p1$SD + p1$SDT

p1$SD_Rel + p1$SDT_Rel

################################################################
############### Arrhenius correction - figure ##################
################################################################
################################################################
checkParamArrheniusCorrection()

