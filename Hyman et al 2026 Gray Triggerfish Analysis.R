# ==============================================================================
## Title: Depth and water temperature drive elevated post-release mortality of gray triggerfish (Balistes capriscus)
## Authors: A. Challen Hyman, Chloe Ramsay, Sean Wilms, and Thomas K. Frazer
## Year: 2025
##
## Description:
##   This script fits a series of Cormack–Jolly–Seber (CJS) tag-recapture models
##   using Bayesian inference (Stan) to estimate post-release mortality of
##   gray triggerfish as a function of depth, temperature, and size
##   selectivity. The script reproduces all tables and figures referenced in
##   Hyman et al. (2026) except conceptual tables/figures and maps.
##
## Key outputs:
##   - Model fits saved as .rds files
##   - Tables 2–4 (model comparison and parameter summaries)
##   - Figures 3–6 and Supplemental Figures S1–S4
##
## Data requirements:
##   All input CSV files must be present in the working directory:
##     - GTF_Data_for_CJS.csv
##     - FHO_GTF.csv
##     - Rec_GTF.csv
##     - GT Selectivity.csv
##     - Longterm GTF.csv
##
## Reproducibility notes:
##   - Random seed is fixed (set.seed(1234))
##   - Stan models must be compiled locally
##   - Results depend on Stan version and C++ toolchain
##
## Last updated: December 19th, 2025
#------------------~-------------------~------------------~--------------------#
################################### Libraries ##################################
#------------------~-------------------~------------------~--------------------#
## Syntax packages
suppressMessages(library(dplyr))          ## Data wrangling
suppressMessages(library(tidyr))          ## Data wrangling
suppressMessages(library(tidyverse))      ## Data wrangling
suppressMessages(library(lubridate))      ## Date handling
suppressMessages(library(readxl))         ## Excel import (used for selectivity inputs)
suppressMessages(library(tidyselect))     ## Core data manipulation and I/O

## Spatial packages (used for long-term observer data)
suppressMessages(library(sp))
suppressMessages(library(sf))

## Graphics packages
suppressMessages(library(cowplot))        ## Visualization of multi-panel figures
suppressMessages(library(bayesplot))      ## Visualization of Bayesian model outputs
suppressMessages(library(xtable))         ## Table export
suppressMessages(library(ggpubr))         ## Visualization of multi-panel figures
suppressMessages(library(ggnewscale))     ## Visualization of multiple legends in ggplot
suppressMessages(library(ggdist))         ## Visualization
suppressMessages(library(scales))         ## Visualization

## Modeling packages
suppressMessages(library(rstan))          ## Bayesian model fitting
suppressMessages(library(loo))            ## Model selection with LOO

#------------------~-------------------~------------------~--------------------#
################################## Housekeeping ################################
#------------------~-------------------~------------------~--------------------#
## Clear the workspace to avoid accidental dependence on cached objects
rm(list=ls(all=TRUE)) 

## Fix RNG state for full reproducibility of simulations and posterior summaries
set.seed(1234)

#------------------~-------------------~------------------~--------------------#
############################## User-defined functions ##########################
#------------------~-------------------~------------------~--------------------#

## Negation of %in% for cleaner subsetting
"%nin%" <- Negate("%in%")

## Custom ggplot theme used throughout all figures
My_theme <- function(){
  theme_bw(base_family = 'serif')%+replace%
    theme(axis.text = element_text(size = 25, color = 1),
          strip.text.x = element_text(size = 25, margin = margin(0.25,0,0.25,0, "cm")),
          strip.text.y = element_text(size = 25, margin = margin(0, 0.25, 0, 0.25, "cm"), angle = 270),
          axis.title = element_text(size = 30),
          legend.text = element_text(size = 25),
          legend.title = element_text(size = 30))}
theme_set(My_theme())  

## Wrapper function to compute Maximum a posteriori (MAP) and central credible interval
## x  : vector of posterior samples
## CI : width of Bayesian credible interval (default = 80%)
## Returns: named vector (lower CI, MAP, upper CI)
MAP <- function(x, CI = 0.8) {
  ## x is a vector of posterior samples
  ## CI is the desired Bayesian credible interval
  
  # Calculate quantile bounds
  LCL <- (1 - CI) / 2
  UCL <- 1 - LCL
  CI_bounds <- quantile(x, probs = c(LCL, UCL))
  
  # MAP estimate from kernel density
  dens <- density(x)
  map_est <- dens$x[which.max(dens$y)]
  
  # Return vector: lower CI, MAP, upper CI
  return(c(CI_lower = CI_bounds[1], MAP = map_est, CI_upper = CI_bounds[2]))
}

## Modified softplus transformation
## Used to impose threshold-like nonlinear effects for depth and temperature
## while maintaining differentiability for Stan sampling. Included here for 
## conceptualization of model outputs.
softplus_modified <- function(x, theta){
  log1p(1+exp(x - theta))-log1p(1+exp(-theta))
}

#------------------~-------------------~------------------~--------------------#
################################## Load data ###################################
#------------------~-------------------~------------------~--------------------#
## NOTE: setwd() is intentionally user-specific and should be changed externally
setwd("C:/Users/ichal/OneDrive/Documents/Hyman gag grouper models/Hyman et al 2025 Discard mortality")


## Main tag-recapture dataset
## One row per released fish
## Required columns include: Recap, Start_Index, End_Index, Depth, Temp, FLM, layer
GTF_Data <- read.csv("GTF_Data_for_CJS.csv")

## For-hire effort proxy (observer data)
FHO  <- read.csv("FHO_GTF.csv")  

## Private-recreational effort proxy (SRFS)
Rec <- read.csv("Rec_GTF.csv") 

## Size selectivity curve
Selectivity <- read.csv("GT Selectivity.csv")

## Long-term observer data for projections
GTF_LT <- read.csv("Longterm GTF.csv")%>%filter(VesselType=="C")
#------------------~-------------------~------------------~--------------------#
################################## Load models #################################
#------------------~-------------------~------------------~--------------------#
## Toggle whether models are re-fit or loaded from disk
run <- T
if (run) {
  g0  <- stan_model("CJS_GTF_0.stan")
  g1  <- stan_model("CJS_GTF_1.stan")
  g2  <- stan_model("CJS_GTF_2.stan")
}

#------------------~-------------------~------------------~--------------------#
################################ Format data ###################################
#------------------~-------------------~------------------~--------------------#
## Standardize effort covariates to unit scale
FHO <- FHO[1:12,-1]
FHO <- FHO/max(FHO)

Rec <- Rec[1:12,-1]
Rec <- Rec/max(Rec)

## Reformat location data for spatial indexing
Location_index <- data.frame(layer = sort(unique(GTF_Data$layer)), ID = 1:16)
GTF_Data$loc <- Location_index$ID[match(GTF_Data$layer, Location_index$layer)]

## Recapture indexing (0 = censored)
GTF_Data$Recap_index <- ifelse(GTF_Data$Recap==1, GTF_Data$End_Index, 0)
GTF_Data$Time_at_large <- ifelse(GTF_Data$Recap==0, max(GTF_Data$Time_at_large, na.rm = T), GTF_Data$Time_at_large)

## Fractional exposure within release interval (for scaling, used internally in models)
GTF_Data$Frac <- as.numeric((90-(as.Date(GTF_Data$SERIESDATE)-as.Date(GTF_Data$StartDate)))/90)

## Add Selectivity based on Garner et. al 2017
GTF_Data$Selectivity <- Selectivity$P[match(GTF_Data$FLM*1000, Selectivity$Lengths)]

## Long-term covariates
GTF_LT$Temp <- GTF_LT$X0_m_Temp
GTF_LT$Year <- year(as.Date(GTF_LT$SERIESDATE))
GTF_LT$Year[which(GTF_LT$Year==2021)] <- 2020
GTF_LT$Month <- month(as.Date(GTF_LT$SERIESDATE), label = T)

#------------------~-------------------~------------------~--------------------#
############################### Summary data ###################################
#------------------~-------------------~------------------~--------------------#
## Summary statistics reported in Results
Summary <- GTF_Data%>%reframe(Depth = c(min(Depth), median(Depth),max(Depth)),
                              SST = c(min(Temp), median(Temp), max(Temp)))%>%t()
colnames(Summary) <- c("Min", "Med", "Max")

## Creates summary statistics referenced in results section
print(xtable(Summary),only.contents=TRUE, include.rownames=FALSE, 
      include.colnames=T, floating=F, sanitize.rownames.function = identity,
      sanitize.text.function = identity)

#------------------~-------------------~------------------~--------------------#
################################## Figure 3 ####################################
#------------------~-------------------~------------------~--------------------#
## This creates Figure 3A: recapture proportion vs depth
Sum_depth <- GTF_Data%>%mutate(Depth = round(Depth/3)*3)%>%
  filter(Depth >19)%>%
  group_by(Depth)%>%
  summarise(Recap = mean(Recap),
            N = length(Depth))%>%
  ggplot()+
  geom_point(aes(Depth, Recap, size = N), alpha = 0.9, col = 'dodgerblue4')+
  ylab("")+labs(size = "Sample size")+theme(legend.position = "top")+xlab("Depth (m)")+
  scale_size_continuous(range = c(2, 6))

## This creates Figure 3B: recapture proportion vs SST
Sum_SST <- GTF_Data%>%
  mutate(SST = round(Temp))%>%
  filter(SST >15)%>%
  group_by(SST)%>%
  summarise(Recap = mean(Recap),
            N = length(Depth))%>%
  ggplot()+
  geom_point(aes(SST, Recap, size = N), alpha = 0.9, show.legend = F, , col = 'dodgerblue4')+labs(size = "Sample size", x = "SST (\u00B0C)")+
  ylab("")+ylim(0.025,0.15)+scale_size_continuous(range = c(2, 6))

## Assemble Figure 3
annotate_figure(ggarrange(Sum_depth, Sum_SST, ncol = 1, labels = c("A", "B"), 
                          font.label = list(size = 30),
                          label.x = 0.02,
                          label.y = c(0.87, 0.99)),
                left = text_grob("Proportion recaptured", 
                                 rot = 90, vjust = 1, size = 35, family = 'serif'))
ggsave("Hyman et al 2025 GTF Summary.png", dpi = 500, device = "png")

#------------------~-------------------~------------------~--------------------#
################################### Modeling ###################################
#------------------~-------------------~------------------~--------------------#
## Prepare data
stan_data <- list(
  N = nrow(GTF_Data),                     ## Number of rows in dataset
  R = GTF_Data$Recap_index,               ## Recapture status for each fish row (recapture time or 0)             
  frac = GTF_Data$Frac,                   ## Time at large
  T = 12,                                 ## Total time intervals
  L = 16,                                 ## Total number of locations
  t0 = GTF_Data$Start_Index,              ## Release times
  loc = GTF_Data$loc,                     ## Release locations
  Rec = Rec,                              ## Private recreational effort
  FHO = FHO,                              ## For-hire fishing pressure (using observer data as a proxy)
  Depth = GTF_Data$Depth,                 ## Depth
  Selectivity = GTF_Data$Selectivity,     ## Size selectivity
  Temp = GTF_Data$Temp                    ## Temperature (SST)
)

## Run MCMC
if(run ==T){
  ## Use all cores
  options(mc.cores = parallel::detectCores()) 
  
  ## Run and save each model to hard drive
  g_0 <- sampling(g0, data = stan_data, chains = 4, iter = 10000, refresh = 1000)
  g_0@stanmodel@dso <- new("cxxdso")
  saveRDS(g_0, file = "GTF_0.rds") 
  
  g_1 <- sampling(g1, data = stan_data, chains = 4, iter = 10000, refresh = 1000)
  g_1@stanmodel@dso <- new("cxxdso")
  saveRDS(g_1, file = "GTF_1.rds") 
  
  g_2 <- sampling(g2, data = stan_data,cores = 4,  chains = 4, iter = 10000, refresh = 1000)
  g_2@stanmodel@dso <- new("cxxdso")
  saveRDS(g_2, file = "GTF_2.rds") 
} else {
  g_0 <- readRDS("GTF_0.rds")  
  g_1 <- readRDS("GTF_1.rds")  
  g_2 <- readRDS("GTF_2.rds")  
}

#------------------~-------------------~------------------~--------------------#
################################### Table 2 ####################################
#------------------~-------------------~------------------~--------------------#
## Extract log_lik from the generated quantities block
log_lik_matrix_0 <- rstan::extract(g_0, pars = "log_lik")$log_lik
log_lik_matrix_1 <- rstan::extract(g_1, pars = "log_lik")$log_lik
log_lik_matrix_2 <- rstan::extract(g_2, pars = "log_lik")$log_lik

## Compute approximate leave-one-out cross-validation
loo_0 <- loo(log_lik_matrix_0)
loo_1 <- loo(log_lik_matrix_1)
loo_2 <- loo(log_lik_matrix_2)

## Run ELPD comparisons based on PSIS-LOO-IC
Comparisons <- loo_compare(loo_0,
                           loo_1,
                           loo_2
)%>%as.data.frame()%>%round(.,2)

## Assemble Table 2
Comparisons$Model <- paste0("$g_",c(2,1,0),"$")
Comparisons <- Comparisons[,c(9,7,3,1,2)]
colnames(Comparisons) <- c("$Model$", 
                           "$LOO$", 
                           "$ELPD_{LOO}$", 
                           "$\\Delta_{ELPD}$",
                           "$SE_{\\Delta_{ELPD}}$")
row.names(Comparisons) <- NULL

## Print Table 2 in Latex code
print(xtable(Comparisons),only.contents=TRUE, include.rownames=FALSE, 
      include.colnames=T, floating=F, sanitize.rownames.function = identity,
      sanitize.text.function = identity)

#------------------~-------------------~------------------~--------------------#
################################### Table 3 ####################################
#------------------~-------------------~------------------~--------------------#
## Extract best-fitting model estimates
Model_values <- extract(g_2)

## Cggregate values to data frame
Parameter_chains <- data.frame(
  alpha_0 = Model_values$alpha[,1],
  alpha_1 = Model_values$alpha[,2],
  alpha_2 = Model_values$alpha[,3],
  gamma = Model_values$gamma,
  beta_D = Model_values$beta_D,
  beta_T = Model_values$beta_T,
  theta_D = Model_values$Depth_thresh,
  theta_T = Model_values$Temp_thresh
)

## Calculate summary statistics
Sum_table <- Parameter_chains%>%apply(., 2, function(x){
  quantile(x, c(0.1, 0.5, 0.9))
})%>%t()%>%as.data.frame()

## Prepare table for use in Latex
rownames(Sum_table) <- paste0("$\\", rownames(Sum_table), "$")

## Print Table 3 in Latex code
print(xtable(Sum_table),only.contents=TRUE, include.rownames=T, 
      include.colnames=T, floating=F, sanitize.rownames.function = identity,
      sanitize.text.function = identity)

#------------------~-------------------~------------------~--------------------#
################################## Figure 4 ####################################
#------------------~-------------------~------------------~--------------------#
## Generate prediction grid for Depth and Temperature
Sim_data <- expand.grid(Depth = seq(20, 70, length = 100),
                        Temp = c(15, 22.5, 30))

## Define number of posterior draws to simulate
N <- 1000 

## Set baseline survival
S <- 1

## Simulate post-release mortality across posterior draws
Recapture <- NULL
for (i in 1:N){
  ## Compute additive effects on the log-scale using the softplus transformation
  Temp_pred  <- Parameter_chains$beta_T[i]*softplus_modified(Sim_data$Temp, theta = Parameter_chains$theta_T[i])
  Depth_pred <- Parameter_chains$beta_D[i]*softplus_modified(Sim_data$Depth, theta = Parameter_chains$theta_D[i])
  
  ## Compute the baseline recapture probability, scaled by survival S
  RR <- exp(Depth_pred + Temp_pred)
  Recap_sim <- RR*S
  
  ## Convert to post-release mortality
  Recapture <- rbind(Recapture, 1-Recap_sim)
  
  ## Optional: print progress every 100 iterations to reduce console clutter
  if (i %% 100 == 0) print(paste("Simulation draw:", i))
}

## Summarize posterior predictions for each Depth-Temp combination
Sim_data$Lwr_DM <- apply(Recapture, 2, quantile, probs = 0.1, na.rm = TRUE)  # 10th percentile
Sim_data$Est_DM <- apply(Recapture, 2, quantile, probs = 0.5, na.rm = TRUE)  # Median
Sim_data$Upr_DM <- apply(Recapture, 2, quantile, probs = 0.9, na.rm = TRUE)  # 90th percentile


## Plot conditional effects using ggplot2
ggplot(Sim_data)+
  geom_line(aes(Depth, Est_DM, col = as.factor(Temp)), lwd = 1)+
  geom_ribbon(aes(Depth, ymin = Lwr_DM, ymax = Upr_DM, fill = as.factor(Temp)), alpha = 0.5)+
  scale_color_manual(values = c("#1f77b4",
                                "orange",
                                "darkred"),
                     name = "SST (\u00B0C)") +
  scale_fill_manual(values = c("#1f77b4",
                               "orange",
                               #'red',
                               "darkred"),name = "SST (\u00B0C)")+
  ylab("Post-release mortality")+xlab("Depth (m)")+
  theme(legend.position = "top")+
  facet_wrap(~Temp, labeller = label_bquote(.(Temp)~"\u00B0C"))

## Save figure
ggsave("Hyman et al 2025 GTF CF.png", dpi = 500, device = "png")

#------------------~-------------------~------------------~--------------------#
################################## Figure 5 ####################################
#------------------~-------------------~------------------~--------------------#
## Projected monthly and annual values

## Simulate post-release mortality for each observation in GTF_LT
Recapture <- NULL
for (i in 1:N){
  ## Log-scale additive effects for Temp and Depth
  Temp_pred  <- Parameter_chains$beta_T[i]*softplus_modified(GTF_LT$Temp, theta = Parameter_chains$theta_T[i])
  Depth_pred <- Parameter_chains$beta_D[i]*softplus_modified(GTF_LT$Depth, theta = Parameter_chains$theta_D[i])
  
  ## Baseline recapture probability scaled by survival S
  RR <- exp(Depth_pred + Temp_pred)
  Recap_sim <- RR*S
  
  ## Convert to post-release mortality
  Recapture <- rbind(Recapture, 1-Recap_sim)
}

## Combine posterior simulations with original data
OP_data <- cbind(GTF_LT[,c("Month","Year", "Temp", "Depth")], as.data.frame(t(Recapture)))%>%
  mutate(Year = as.factor(Year))

## Summarize annual post-release mortality (Years > 2020)
annual_data <- OP_data %>%
  filter(as.numeric(as.character(Year)) > 2020) %>%
  select(Year, 5:ncol(.)) %>%   # Select posterior draws
  group_by(Year) %>%
  summarise(across(everything(), ~ mean(.x, na.rm = TRUE))) %>%  # mean across observations
  pivot_longer(cols = 2:ncol(.), names_to = "draw", values_to = "Value")

## Compute credible intervals per year
annual_summary <- annual_data %>%
  group_by(Year) %>%
  summarise(
    min = quantile(Value, 0.1),
    med = median(Value),
    max = quantile(Value, 0.9)
  )

## Create annual figure (Figure 5a)
annual_plot <- ggplot(annual_data) +
  geom_violin(aes(x = Year, y = Value), fill = 'lightblue', col = 'dodgerblue4', lwd = 1) +
  geom_linerange(data = annual_summary, aes(x = Year, ymin = min, ymax = max), lwd = 3, col = 'dodgerblue4') +
  geom_point(data = annual_summary, aes(x = Year, y = med), size = 7, col = 'dodgerblue4', shape = 21, fill = 'lightblue', stroke = 1.5) +
  ylab("") +
  xlab("Year") +
  ylim(0, 0.8) +
  theme(axis.ticks.y = element_blank(),
        axis.text.y = element_blank())


## Summarize monthly post-release mortality (Months > 2020)
monthly_data <- OP_data %>%
  filter(as.numeric(as.character(Year)) > 2020) %>%
  select(Month, 5:ncol(.)) %>% 
  group_by(Month) %>%
  summarise(across(everything(), ~ mean(.x, na.rm = TRUE))) %>%
  rowwise() %>%
  mutate(
    row_Q10 = quantile(c_across(2:ncol(.)), 0.1, na.rm = TRUE),
    row_med  = median(c_across(2:ncol(.)), na.rm = TRUE),
    row_Q90 = quantile(c_across(2:ncol(.)), 0.9, na.rm = TRUE),
    Date = as.Date(paste0("2024-", Month, "-01"), format = '%Y-%m-%d')
  ) %>%
  ungroup()

## Create monthly figure (Figure 5b)
monthly_plot <- ggplot(monthly_data) +
  geom_linerange(aes(x = Date, ymin = row_Q10, ymax = row_Q90), lwd = 3, col = 'dodgerblue4') +
  geom_point(aes(x = Date, y = row_med), size = 5, col = 'dodgerblue4', shape = 21, fill = 'lightblue', stroke = 1.5) +
  ylab("Post-release mortality") +
  scale_x_date(date_breaks = "1 month", date_labels = "%b") +
  xlab("Month") +
  ylim(0, 0.8)

## Combine annual and monthly plots
combined_plot <- ggarrange(monthly_plot, annual_plot, labels = c("A", "B"), font.label = list(color = "black", size = 25))

# Save combined figure
ggsave("GTF Projection.png", dpi = 900, device = "png", width = 20, height = 7)

#------------------~-------------------~------------------~--------------------#
################################## Figure 6 ####################################
#------------------~-------------------~------------------~--------------------#
## Summarize OP_data across posterior draws by Year
long_term_data  <- OP_data%>%.[,c(2,5:1004)]%>%group_by(Year)%>%st_drop_geometry()%>%
  summarise(across(everything(), ~ mean(.x, na.rm = TRUE)))%>%
  pivot_longer(2:1001, names_to = "draw", values_to = "Value")

## Compute annual credible intervals (10th, 50th, 90th percentile)
long_term_summary <- long_term_data %>%group_by(Year)%>%
  summarize(min = quantile(Value, 0.1),
            med = median(Value),
            max = quantile(Value, 0.9))

## Adjust year labels for plotting
long_term_summary <- long_term_summary %>%
  mutate(
    Year = as.character(Year),
    Year = ifelse(Year == "2020", "2020 \u2013 2021", Year)  # Combine 2020-2021 due to COVID-19 pandemic
  )


## Create long-term trend plot
long_term_plot <- ggplot(long_term_summary) +
  geom_linerange(aes(x = Year, ymin = min, ymax = max), lwd = 3, col = 'dodgerblue4') +
  geom_point(aes(x = Year, y = med), size = 7, col = 'dodgerblue4', shape = 21, fill = 'lightblue', stroke = 1.5) +
  xlab("Year") +
  ylab("Post-release mortality") +
  theme_minimal(base_size = 16)

## Save figure
ggsave("GTF LT Projection.png", dpi = 900, device = "png", width = 20, height = 7)
#------------------~-------------------~------------------~--------------------#
################################### Table 4 ####################################
#------------------~-------------------~------------------~--------------------#
## Initialize empty table to store results
Table_4 <- NULL

## Loop over baseline survival scenarios
for (S in c(1, 0.925, 0.85)){
  
  ## Filter long-term dataset for years 2022-2024
  GTF_2022_2024 <- GTF_LT%>%filter(Year >2021)

  ## Initialize matrix to store recapture simulations
  Recapture <- NULL
  for (i in 1:N){
    ## Compute temperature and depth effects using softplus transformation
    Temp_pred  <- Parameter_chains$beta_T[i]*softplus_modified(GTF_2022_2024$Temp, theta = Parameter_chains$theta_T[i])
    Depth_pred <- Parameter_chains$beta_D[i]*softplus_modified(GTF_2022_2024$Depth, theta = Parameter_chains$theta_D[i])
    
    ## Combine baseline survival (S) with depth and temperature effects
    RR <- exp(log(S) + Depth_pred + Temp_pred)
    
    ## Convert to mortality
    Recap_sim <- RR
    Recapture <- rbind(Recapture, 1-Recap_sim)
  }
  
  ## Combine simulation results with covariates
  OP_data <- cbind(GTF_2022_2024[,c("Month","Year", "Temp", "Depth")], as.data.frame(t(Recapture)))%>%
    mutate(Year = as.factor(Year))
  
  ## Summarize annual post-release mortality (10%, 50%, 90% quantiles)
  annual_DM <- OP_data%>%.[,c(2,5:1004)]%>%group_by(Year)%>%st_drop_geometry()%>%
    summarise(across(everything(), ~ mean(.x, na.rm = TRUE)))%>%
    pivot_longer(2:1001, names_to = "inter", values_to = "Value")%>%
    group_by(Year)%>%
    summarize(min = quantile(Value, 0.1, na.rm = T),
              med = median(Value, na.rm = T),
              max = quantile(Value, 0.9, na.rm = T))%>%as.data.frame()
  colnames(annual_DM) <- c("Year","10%", "50%","90%")
  
  ## Add baseline survival scenario for tracking
  annual_DM$baseline <- S
  
  ## Append to master table
  Table_4 <- rbind(Table_4, annual_DM)
}

## Reshape table for publication format
Table_4 <- Table_4%>%
  pivot_wider(names_from = baseline, values_from = c(`10%`, `50%`, `90%`))%>%
  .[,c(1,2,5,8,3,6,9,4,7,10)]

## Print as Latex code
print(xtable(Table_4),only.contents=TRUE, include.rownames=FALSE, 
      include.colnames=T, floating=F, sanitize.rownames.function = identity,
      sanitize.text.function = identity)

#------------------~-------------------~------------------~--------------------#
########################## Supplemental figures ################################
#------------------~-------------------~------------------~--------------------#
## Figure S1
## Initialize selective data frame
Selectivity_grid <- expand.grid(Hooktype = c("C"),
                                Lengths = c(1:1000),
                                Selectivity = NA)

## Set parameter values from Garner et al. (2017)
alpha = 0.045
beta = 0.254
theta = 332.8

## Calculate expected selectivity using exponential-logistic selectivity curve
## σ = exp(βα(θ−l))) / [1-β(1-exp(α(θ-l)))]
Selectivity_grid$Selectivity <- exp(beta*alpha*(theta-Selectivity_grid$Lengths))/(1-beta*(1-exp(alpha*(theta-Selectivity_grid$Lengths))))

## Plot selectivity curve for conceptualization
ggplot(Selectivity_grid)+
  geom_line(aes(Lengths, Selectivity), lwd = 1)+
  ylab(expression(sigma))+xlab("Fork length (mm)")

## Save figure
ggsave("Selectivity_GTF.png", dpi = 900, device = "png", width = 20, height = 7)

#------------------~-------------------#
## Figure S2

## Set reference depths for projections
Depth <- 10:70

## Set varying beta_D values for conceptualization
beta_D <- c(-0.01, -0.05, -0.1)

## Set varying theta_D values for conceptualization
theta_D <- 25

## Intialize object to store values
Conceptual <- NULL
for (i in 1:3){
  ## For each beta_D value, calculate relative post-release survival
  psi <- exp(beta_D[i]*softplus_modified(x = Depth, theta = theta_D))
  
  ## Save as data frame with reference parameter values for tracking
  dat <- data.frame(Depth = Depth, psi = psi, beta = beta_D[i])
  
  ## Save to master data frame
  Conceptual <- rbind(Conceptual, dat)
}

## Plot results
ggplot(Conceptual)+
  geom_point(aes(Depth, psi, col = as.factor(beta)))+
  ylab("Post-release survival")+labs(col = expression(beta[D]))+
  geom_vline(xintercept = theta_D)+
  geom_text(
    data = data.frame(x = 20, y = 0.9),
    aes(x = 20, y = 0.0), label = expression(theta[D]),
    size = 10, family = 'serif'
  )+
  scale_color_viridis_d()

## Save plot
ggsave("Conceptual_GTF.png", dpi = 900, device = "png", width = 20, height = 7)

#------------------~-------------------#
## Figure S3

## Extract phi, eta, psi for posterior predictive checks
phi_post <- Model_values$phi              ## dims: iter x T x L
h_recap_post <- Model_values$h_recap      ## dims: iter x T x L
psi_post <- Model_values$psi              ## dims: iter x N

## Extract values from Stan data list derivedf from from GTF dataset (above)
t0 <- stan_data$t0                        ## Tag-release times 
frac <- GTF_Data$Frac                     ## Fractional exposure in first interval
loc = stan_data$loc                       ## Spatial strata (locations)
Tmax = stan_data$T                        ## Censor time
Selectivity = stan_data$Selectivity       ## Selectivity based on Garner et al. (2017)

## Helper: compute per-draw simulated recapture times (0 = censored i.e., never recaptured)
simulate_recapture <- function(phi_post, h_recap_post, psi_post, t0, frac, loc, Tmax, Selectivity, n_iter) {
  N <- length(t0)
  recaps_sim <- matrix(NA, nrow = n_iter, ncol = N)
  
  for (iter in 1:n_iter) {
    for (i in 1:N) {
      t_start <- t0[i]
      l <- loc[i]
      eta_post <- 1-exp(-h_recap_post[iter,,l]*Selectivity[i])
      ## Survive post-release?
      if(rbinom(1,1, psi_post[iter,i])==0){
        recaps_sim[iter,i] <- 0  # censored
        next
      }
      
      ## Survive time step?
      surv <- phi_post[iter,t_start,l]^(max(frac[i],0))
      if(rbinom(1,1, surv)==0){
        recaps_sim[iter,i] <- 0  # censored
        next
      }
      
      ## Recaptured at time step?
      rec <- 1 - (1 - eta_post[t_start])^(max(frac[i],0))
      if(rbinom(1,1, rec)){
        recaps_sim[iter,i] <- t_start  # recaptured
        recaptured <- TRUE
        next
      }
      
      recaptured <- FALSE
      T <- min(t_start+4, Tmax)
      for (t in (t_start+1):T){
        ## Survive t?
        surv <- phi_post[iter,t,l]
        if(rbinom(1,1, surv)==0){
          recaps_sim[iter,i] <- 0  # censored
          break
        }
        ## Recaptured at time step?
        rec <- 1 - (1 - eta_post[t])
        if(rbinom(1,1, rec)==1){
          recaps_sim[iter,i] <- t  # recaptured
          recaptured <- TRUE
          break
        }
      }
      if (!recaptured) recaps_sim[iter,i] <- 0  # censored
    }
  }
  return(recaps_sim)
}


## Simulate recaptures based on covariates in GTF data frame to compare to true values
recap_sim <- simulate_recapture(phi_post, h_recap_post, psi_post, stan_data$t0, stan_data$frac, stan_data$loc, stan_data$T, stan_data$Selectivity, n_iter = 4000)

## Summarize simulated recaptures
sim_means <- apply(ceiling(recap_sim/13), 1, mean)  # average predicted recapture time per fish

## Figure S3a: Histogram of observed vs predicted
Hist_plot <- ggplot(data.frame(
  predicted = sim_means
)) +
  geom_histogram(aes(x = predicted), fill = "#B3CDE0", alpha = 0.5, col = "#000A39") +
  geom_vline(xintercept = mean(GTF_Data$Recap), lwd = 1, col = "#000A39")+
  labs(y = "Frequency", x = "Proportion recaptured")


## Figure S3b: scatter plot of observed vs mean predicted proportion of fish recaptured
df_ppc <- tibble(
  predicted = apply(ceiling(recap_sim/13), 2, mean),
  observed = GTF_Data$Recap
)

# Bin predicted values
df_ppc2 <- df_ppc %>%
  mutate(pred_bin = cut(predicted, breaks = seq(0, 1, by = 0.001))) %>%
  group_by(pred_bin) %>%
  summarise(
    n = n(),
    obs_mean = mean(observed),
    pred_mean = mean(predicted)
  )%>%filter(n > 10)

# Plot calibration curve
scatter_plot <- ggplot(df_ppc2, aes(x = pred_mean, y = obs_mean, weight = n)) +
  geom_point(aes(), fill = "#B3CDE0", alpha = 0.9, col = "#000A39", shape = 21, size = 5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "#000A39", lwd = 1) +
  labs(x = "Predicted Recapture Rate", y = "Observed Recapture Rate")

## Figure S3c: Comparison of the empirical distribution of observed recapture 
## times to the distributions of 1,000 scans from the posterior predictive distribution. 
y = stan_data$R            ## True recaptures
yrep = recap_sim[1:4000,]  ## First 4000 simulations based on posterior draws

Density_plot <- ppc_dens_overlay(y = y, yrep = yrep) +
  ggplot2::labs(
    x = "Recaptures time step (0 = censored)",
    y = "Kernel density"
  )+theme(legend.position = 'top')+scale_x_continuous(breaks = c(0,3,6,9))

## Figure S3d: Empirical cumulative distribution function (ECDF) of observed recapture times
## compared with ECDFs simulated from posterior predictive draws.
ECDF_plot <- ppc_ecdf_overlay(y, yrep)+
  ggplot2::labs(
    x = "Recaptures time step (0 = censored)",
    y = "ECDF of recaptures"
  )+theme(legend.position = 'top')+scale_x_continuous(breaks = c(0,3,6,9))


## Arrange to create Figure S3
ggarrange(Hist_plot, scatter_plot,
          Density_plot, ECDF_plot, common.legend = T, labels = c("A", "B", "C", "D"),
          font.label = list(size = 20, face = "bold", family = "serif") )
ggsave("Hyman et al 2025 GTF diagnostics.png", dpi = 900, device = "png")


#------------------~-------------------#
## Figure S4

## Aggregate simulated and observed data to data frame
df_ppc <- tibble(
  predicted = apply(ceiling(recap_sim/13), 2, mean),
  observed = GTF_Data$Recap,
  Temp = GTF_Data$Temp,
  Depth = GTF_Data$Depth,
  Size = GTF_Data$FLM*1000
)

## Figure S4a: Observed vs predicted proportion of fish recaptured as a function of depth
Depth_sim <- df_ppc %>%
  filter(Depth >=19)%>%
  mutate(Var = round(Depth/3)*3)%>%
  group_by(Var) %>%
  summarise(
    n = n(),
    obs_mean = mean(observed),
    pred_mean = mean(predicted)
  )%>%filter(n > 10)%>%
  ggplot() +
  geom_point(aes(x = Var, y = obs_mean, fill = "Observed"), alpha = 0.9, col = "#000A39", shape = 21, size = 5) +
  geom_point(aes(x = Var, y = pred_mean, fill = "Predicted"), alpha = 0.9, col = "#000A39", shape = 21, size = 5) +
  scale_fill_manual(
    name = "", values = c("Observed" = "#000A39", "Predicted" = "#B3CDE0")
  )+ylim(0,0.15)+
  labs(x = "Depth (m)", y = "Recapture Rate")+theme(legend.position = "top")

## Figure S4b: Observed vs predicted proportion of fish recaptured as a function of SST
Temp_sim <- df_ppc %>%
  mutate(Var = round(Temp))%>%
  group_by(Var) %>%
  summarise(
    n = n(),
    obs_mean = mean(observed),
    pred_mean = mean(predicted)
  )%>%filter(n > 10)%>%
  ggplot() +
  geom_point(aes(x = Var, y = obs_mean, fill = "Observed"), alpha = 0.9, col = "#000A39", shape = 21, size = 5) +
  geom_point(aes(x = Var, y = pred_mean, fill = "Predicted"), alpha = 0.9, col = "#000A39", shape = 21, size = 5) +
  scale_fill_manual(
    name = "", values = c("Observed" = "#000A39", "Predicted" = "#B3CDE0")
  )+ylim(0,0.15)+
  labs(x = "SST (\u00B0C)", y = "")+theme(legend.position = "top",
                                          axis.title.y = element_blank(),
                                          axis.text.y = element_blank(),
                                          axis.ticks.y = element_blank(),
                                          plot.margin = margin(t = 5, r = 5, b = 5, l = 0))

## Figure S4c: Observed vs predicted proportion of fish recaptured as a function of size
Size_sim <- df_ppc %>%
  mutate(Var = round(Size/1)*1)%>%
  group_by(Var) %>%
  summarise(
    n = n(),
    obs_mean = mean(observed),
    pred_mean = mean(predicted)
  )%>%filter(n > 10)%>%
  ggplot() +
  geom_point(aes(x = Var, y = obs_mean, fill = "Observed"), alpha = 0.9, col = "#000A39", shape = 21, size = 5) +
  geom_point(aes(x = Var, y = pred_mean, fill = "Predicted"), alpha = 0.9, col = "#000A39", shape = 21, size = 5) +
  scale_fill_manual(
    name = "", values = c("Observed" = "#000A39", "Predicted" = "#B3CDE0")
  )+scale_y_continuous(labels = label_number(accuracy = 0.01))+
  labs(x = "Size (mm)", y = "Recapture Rate")+theme(legend.position = "top")


## Extract legend from first plot
legend <- get_legend(Depth_sim)


## Remove legends from the other plots
p1_noleg <- Depth_sim + theme(legend.position = "none")
p2_noleg <- Temp_sim + theme(legend.position = "none")
p3_noleg <- Size_sim + theme(legend.position = "none")


## Combine plots (2x2 grid with legend)
combined1 <- plot_grid(
  p1_noleg, p2_noleg, labels = c("A", "B"), label_size = 25,
  label_x = c(0.1,0),
  ncol = 2, rel_widths = c(1.1, 1)
)

combined <- plot_grid(
  combined1, p3_noleg, labels = c("", "C"), label_size = 25,
  label_x = c(0,0.052),
  ncol = 1, align = "!hv"
)


## Add legend below or to the right
final_plot <- plot_grid(
  legend,combined, 
  ncol = 1,              # legend below
  rel_heights = rev(c(1, 0.1))  # adjust spacing
)

## Print plot
final_plot

## Save plot
ggsave("Hyman et al 2025 GTF diagnostics 2.png", dpi = 900, device = "png")
