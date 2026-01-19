###read in data excluding communities with contamination 
vel_dat<-velvet_week1_6ForANAlysis_ng
str(vel_dat)

library(dplyr)
library (ggplot2)
# create new factor that contains treatment and week information 
vel_dat <- mutate(vel_dat, tr = paste(treatment, week)) 

#overall abundances of each isolate for average networks per treatment 
by_tr <- vel_dat %>% group_by(tr)
by_tr %>% summarise(
  b = mean(b),
  o = mean(o),
  p = mean(p),
)
#estimate isolate abundance without plasmids (low abundance in the network)
by_np <- vel_dat %>% group_by(tr)
by_np%>% summarise(
  b0 = mean(b0),
  o0 = mean(o0),
  p0 = mean(p0),
)

by_r <- vel_dat %>% 
  na.omit() %>%
  group_by(tr)

#estimate mean infection rate (based on the rate already in the table)
by_r %>% summarise(
  bpb12 = mean(bEr),
  opb12 = mean(oEr),
  ppb12 = mean(pEr),
  bpb5 = mean(bGr),
  opb5 = mean(oGr),
  ppb5 = mean(pGr),
  bpkjk5 = mean(bPr),
  opkjk5 = mean(oPr),
  ppkjk5 = mean(pPr),
)
#read in networks (need individual csv files for each network)
setwd("/Users/ds306/Documents/3h3p_ms/")

files <- list.files(pattern = ".csv") 
files

n1 <- read.csv("~/Documents/3h3p_ms/n1.csv", header=FALSE)
n6 <- read.csv("~/Documents/3h3p_ms/n6.csv", header=FALSE)

w1 <- read.csv("~/Documents/3h3p_ms/w1.csv", header=FALSE)
w6 <- read.csv("~/Documents/3h3p_ms/w6.csv", header=FALSE)

s1 <- read.csv("~/Documents/3h3p_ms/s1.csv", header=FALSE)
s6 <- read.csv("~/Documents/3h3p_ms/s6.csv", header=FALSE)

# add names WS= white species to get relative abundance right across the networks
colnames(n1) <- c("pB12","pB5","pkjk5")
rownames(n1) <- c("Bor","Och","Pse","WS")

colnames(n6) <- c("pB12","pB5","pkjk5")
rownames(n6) <- c("Bor","Och","Pse","WS")

colnames(w1) <- c("pB12","pB5","pkjk5")
rownames(w1) <- c("Bor","Och","Pse","WS")

colnames(w6) <- c("pB12","pB5","pkjk5")
rownames(w6) <- c("Bor","Och","Pse","WS")

colnames(s1) <- c("pB12","pB5","pkjk5")
rownames(s1) <- c("Bor","Och","Pse","WS")

colnames(s6) <- c("pB12","pB5","pkjk5")
rownames(s6) <- c("Bor","Och","Pse","WS")

#add white species different to most abundant network (the "white species" accounts for differences in network sizes and allows visual comparison)
#n 1 and 6:	60.68928	110.054582
#s 1 and 6:	121.266	47.246
#w 1 and 6:	0	78.20605

#lower abundance (these are the average isolate numbers without plasmid)
#treatment n
low.n1a<-c(0.6,92,22.9,60.6)
low.n1<-as.numeric(low.n1a)
names(low.n1) <- rownames(n1)

low.n6a<-c(6.78,86.2,0.111,110.1)
low.n6<-as.numeric(low.n6a)
names(low.n6) <- rownames(n6)

#treatmnet w
low.w1a<-c(1,0.875,2, 0)
low.w1<-as.numeric(low.w1a)
names(low.w1) <- rownames(w1)

low.w6a<-c(7.29,7.14,0.143, 78.2)
low.w6<-as.numeric(low.w6a)
names(low.w6) <- rownames(w6)

#treatment s
low.s1a<-c(0.6,1.4,1.4, 121.3)
low.s1<-as.numeric(low.s1a)
names(low.s1) <- rownames(s1)

low.s6a<-c(1.2,1.4,0,47.3)
low.s6<-as.numeric(low.s6a)
names(low.s6) <- rownames(s6)

#plot networks
require(bipartite)
par(mfrow=c(3,2))
plotweb(n1, method = "cca",low.abun=low.n1,col.interaction ="#ABD9E9",low.spacing = 0.1,labsize = 1.6,empty=F,
        bor.col.interaction="#ABD9E9",high.abun.col="green",low.abun.col="gray64", arrow="down", high.lablength=6, low.lablength=3)
#mtext("n1",3, line=0,adj = -0.1 )
plotweb(n6, method = "cca",low.abun=low.n6,col.interaction ="#ABD9E9",low.spacing = 0.1,labsize = 1.6,empty=F,
        bor.col.interaction="#ABD9E9",high.abun.col="green",low.abun.col="gray64", arrow="down", high.lablength=6, low.lablength=3)
#mtext("n6",3, line=0,adj = -0.1 )
plotweb(w1, method = "cca",low.abun=low.w1, col.interaction ="#FDAE61",low.spacing = 0.1,labsize = 1.6,empty=F,
        bor.col.interaction="#FDAE61",high.abun.col="green",low.abun.col="gray64", arrow="down", high.lablength=6, low.lablength=3)
#mtext("w1",3, line=0,adj = -0.1 )
plotweb(w6, method = "cca",low.abun=low.w6, col.interaction ="#FDAE61",low.spacing = 0.1,labsize = 1.6,empty=F,
        bor.col.interaction="#FDAE61",high.abun.col="green",low.abun.col="gray64", arrow="down", high.lablength=6, low.lablength=3)
#mtext("w6",3, line=0,adj = -0.1 )
plotweb(s1, method = "cca",low.abun=low.s1, col.interaction ="#D73027",low.spacing = 0.1,labsize = 1.6,empty=F,
        bor.col.interaction="#D73027",high.abun.col="green",low.abun.col="gray64", arrow="down", high.lablength=6, low.lablength=3)
#mtext("s1",3, line=0,adj = -0.1 )
plotweb(s6, method = "cca",low.abun=low.s6, col.interaction ="#D73027",low.spacing = 0.1,labsize = 1.6,empty=F,
        bor.col.interaction="#D73027",high.abun.col="green",low.abun.col="gray64", arrow="down", high.lablength=6, low.lablength=3)
#mtext("s6",3, line=0,adj = -0.1 )

## extract network metrics, we will use connectance (conn) and average plasmid generality (GqM) 
require(vegan)
MetricData<-vel_dat%>%
  transform(shannon = diversity(select(vel_dat,b,o,p)),
            simpson = diversity(select(vel_dat,b,o,p), index ="simpson"),
            hostp=rowSums(vel_dat[,c("b",	"o",	"p")]!= 0)*3,
            conn =rowSums(vel_dat[,c("bE",	"oE",	"pE",	"bG",	"oG",	"pG",	"bP",	"oP",	"pP")]!= 0)/(h=rowSums(vel_dat[,c("b",	"o",	"p")]!= 0)*3),
            GqE=2^(diversity(select(vel_dat, bE,oE,pE))),
            GqG=2^(diversity(select(vel_dat,bG,oG,pG))),
            GqP=2^(diversity(select(vel_dat,bP,oP,pP))),
            GqM= ((GqE=2^(diversity(select(vel_dat, bE,oE,pE)))) + (GqG=2^(diversity(select(vel_dat,bG,oG,pG)))) + (GqP=2^(diversity(select(vel_dat,bP,oP,pP)))))/3,
            GqBo=2^(diversity(select(vel_dat, bE,bG,bP))),
            GqOc=2^(diversity(select(vel_dat,oE,oG,oP))),
            GqPs=2^(diversity(select(vel_dat,pE,pG,pP))))


#save the data file 
write.csv(MetricData,"/Users/DS306/Documents/metrics.csv", row.names = FALSE)
#

#### order treatment levels 
MetricData$treatment = factor(MetricData$treatment , levels=c("n","w","s"))
MetricData$tr = factor(MetricData$tr , levels=c("n 1","n 6","w 1","w 6","s 1","s 6"))

#############################################
#####bmrs model##############################done 3 March 25

# without specified priors 
lme1.G <- brm(
  GqM ~ week * treatment + (1 | repl),
  data = MetricData,
  family = gaussian(),  # Use an appropriate family for your data
  chains = 4,           # Number of Markov chains (increase for better estimates)
  iter = 2000,          # Number of iterations per chain
  warmup = 500,         # Burn-in period
  cores = parallel::detectCores()) # Uses multiple cores for faster sampling
  
summary(lme1.G)

plot(lme1.G)  # Diagnostic plots
pp_check(lme1.G)  # Posterior predictive checks locking good 

loo(lme1.G)  # Leave-One-Out Cross-Validation
waic(lme1.G) # Widely Applicable Information Criterion

# add weak priors 
prior <- prior(normal(0, 10), class = "b") +
  prior(cauchy(0, 2), class = "sd")
lme1.Ga <- brm(GqM ~ week * treatment + (1 | repl), 
              data = MetricData, 
              family = gaussian(),
              prior = prior)

plot(lme1.Ga)  # Diagnostic plots
pp_check(lme1.Ga)  # Posterior predictive checks locking good 

#looks like weak priors improve the model but not by much
loo(lme1.Ga)
loo(lme1.G)

#test specific hypotheses
hypothesis(lme1.Ga, "week < 0")
hypothesis(lme1.Ga, "treatments < treatmentw")
hypothesis(lme1.Ga, "week:treatmentw < 0")  # Test interaction effect

#model results 
summary(lme1.Ga)
#plots
plot(conditional_effects(lme1.Ga), points = TRUE)
#as table
posterior_summary(lme1.Ga)

###############plot posterior distribution for treatment levels 

# Ensure required libraries are loaded
library(tidyverse)
library(brms)

# Create new data for weeks 1 & 6
new_data <- expand.grid(
  week = c(1, 6), 
  treatment = unique(MetricData$treatment)
)

# Get posterior predictions (ignore random effects)
posterior_preds <- posterior_epred(lme1.Ga, newdata = new_data, allow_new_levels = TRUE, re_formula = NA)

# Convert to tidy format
posterior_df <- posterior_preds %>%
  as_tibble() %>%
  mutate(iteration = row_number()) %>%   # Add iteration number
  pivot_longer(cols = -iteration, names_to = "newdata_index", values_to = "predicted") %>%
  mutate(
    week = rep(new_data$week, times = nrow(posterior_preds)),  # Repeat properly
    treatment = rep(new_data$treatment, times = nrow(posterior_preds))
  )


######
# Reorder the 'week' variable to ensure that week 6 is below week 1
posterior_df$week <- factor(posterior_df$week, levels = c(1, 6))

####including observed values 

#Ensure you have a dataset (observed_data) that contains the observed values with corresponding week and treatment labels
library(dplyr)
observed_dataA<-MetricData  %>%
  select(week, treatment, GqM)

# Define treatment-specific colors
treatment_colors <- c("#ABD9E9", "#FDAE61", "#D73027")

# Assuming 'observed_data' contains observed values with 'week' and 'treatment' columns
ggplot(posterior_df, aes(x = predicted, fill = treatment)) +
  geom_density(alpha = 0.5) +  # Density plot with transparency
  
  # Adjust y-position dynamically for each treatment
  geom_point(data = observed_dataA, aes(x = GqM, y = as.numeric(treatment) * -0.5, color = treatment), 
             position = position_jitter(width = 0.02, height = 0), size = 2) +  # Spread points horizontally
  
  scale_fill_manual(values = treatment_colors) +  # Custom treatment-specific colors
  scale_color_manual(values = treatment_colors) +  # Match observed points to treatment colors
  
  facet_wrap(~week, scales = "fixed", ncol = 1) +  # Arrange week panels in a single column
  
  labs(
    title = "",
    x = "Predicted Value for Generality",
    y = "Density",
    fill = "Treatment",
    color = "Observed"
  ) +
  theme_minimal()

##############plot with posterior distribution and boxplots
# Define treatment-specific colors
treatment_colors <- c("#ABD9E9", "#FDAE61", "#D73027")

plot2_bmrs<-ggplot(posterior_df, aes(x = predicted, fill = treatment)) +
  geom_density(alpha = 0.5) +  # Density plot with transparency
  
  # Add boxplots for observed data
    geom_boxplot(data = observed_dataA, aes(
    x = GqM, 
    y = as.numeric(factor(treatment)) * -0.5,  # Ensure correct spacing per treatment
    fill = treatment, 
    group = interaction(treatment, week)
  ), 
  width = 0.3, alpha = 0.5, outlier.shape = NA) +  # Set width and remove outliers
  
  # Add observed points (jittered)
  geom_point(data = observed_dataA, aes(
    x = GqM, 
    y = as.numeric(factor(treatment)) * -0.5, 
    color = treatment
  ), 
  position = position_jitter(width = 0.1, height = 0), size = 2) + 
  
  scale_fill_manual(values = treatment_colors) +  # Custom treatment-specific colors
  scale_color_manual(values = treatment_colors) +  # Match observed points to treatment colors
  
  # Rename facet labels
  facet_wrap(~week, scales = "fixed", ncol = 1, 
             labeller = labeller(week = c("1" = "Week 1", "6" = "Week 6"))) +  
  
  labs(
    title = "",
    x = "Predicted Value for Generality",
    y = "Density",
    fill = "",   # Remove "Treatment" label in legend
    color = ""   # Remove "Observed" label in legend
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 14),  
    strip.text = element_text(hjust = 0, size = 13, face = "italic", color = "gray30"),  
    legend.position = "none"  # 🔹 Removes legend (treatment label)
  )

######conn
#################
#############################################
#####bmrs model##############################done 3 March 25

# without specified priors 
lme1.C <- brm(
  conn ~ week * treatment + (1 | repl),
  data = MetricData,
  family = gaussian(),  # Use an appropriate family for your data
  chains = 4,           # Number of Markov chains (increase for better estimates)
  iter = 2000,          # Number of iterations per chain
  warmup = 500,         # Burn-in period
  cores = parallel::detectCores()) # Uses multiple cores for faster sampling

summary(lme1.C)

plot(lme1.C)  # Diagnostic plots
pp_check(lme1.C)  # Posterior predictive checks locking good 

loo(lme1.C)  # Leave-One-Out Cross-Validation
waic(lme1.C) # Widely Applicable Information Criterion

# add weak priors 
prior <- prior(normal(0, 10), class = "b") +
  prior(cauchy(0, 2), class = "sd")
lme1.Ca <- brm(conn ~ week * treatment + (1 | repl), 
               data = MetricData, 
               family = gaussian(),
               prior = prior)

plot(lme1.Ca)  # Diagnostic plots
pp_check(lme1.Ca)  # Posterior predictive checks locking good 

#looks like weak priors improve the model but not by much
loo(lme1.Ca)
loo(lme1.C)

#test specific hypotheses
hypothesis(lme1.Ca, "week < 0")
hypothesis(lme1.Ca, "treatments < treatmentw")
hypothesis(lme1.Ca, "week:treatmentw < 0")  # Test interaction effect

#model results 
summary(lme1.Ca)
#plots
plot(conditional_effects(lme1.Ca), points = TRUE)
#as table
posterior_summary(lme1.Ca)

##### plot posterior distribution for single parameters 

library(posterior)

################plot posterior distribution for treatment levels 
# Ensure required libraries are loaded
library(tidyverse)
library(brms)

# Create new data for weeks 1 & 6
new_dataC <- expand.grid(
  week = c(1, 6), 
  treatment = unique(MetricData$treatment)
)

# Get posterior predictions (ignore random effects)
posterior_predsC <- posterior_epred(lme1.Ca, newdata = new_dataC, allow_new_levels = TRUE, re_formula = NA)

# Convert to tidy format
posterior_dfC <- posterior_predsC %>%
  as_tibble() %>%
  mutate(iteration = row_number()) %>%   # Add iteration number
  pivot_longer(cols = -iteration, names_to = "newdata_index", values_to = "predicted") %>%
  mutate(
    week = rep(new_dataC$week, times = nrow(posterior_predsC)),  # Repeat properly
    treatment = rep(new_dataC$treatment, times = nrow(posterior_predsC))
  )

# Reorder the 'week' variable to ensure that week 6 is below week 1
posterior_df$week <- factor(posterior_df$week, levels = c(1, 6))

#Ensure you have a dataset (observed_data) that contains the observed values with corresponding week and treatment labels
library(dplyr)
observed_data<-MetricData  %>%
  select(week, treatment, conn)

############## Figure 2 #######################################################
###############################################################################
##############plot with posterior distribution and boxplots
# Define treatment-specific colors
treatment_colors <- c("#ABD9E9", "#FDAE61", "#D73027")

plot1_bmrs<-ggplot(posterior_dfC, aes(x = predicted, fill = treatment)) +
  geom_density(alpha = 0.5) +  # Density plot with transparency
  
  # Add boxplots for observed data
  geom_boxplot(data = observed_data, aes(
    x = conn, 
    y = as.numeric(factor(treatment)) * -1.5,  # Ensure correct spacing per treatment
    fill = treatment, 
    group = interaction(treatment, week)
  ), 
  width = 1.2, alpha = 0.5, outlier.shape = NA) +  # Set width and remove outliers
  
  # Add observed points (jittered)
  geom_point(data = observed_data, aes(
    x = conn, 
    y = as.numeric(factor(treatment)) * -1.5, 
    color = treatment
  ), 
  position = position_jitter(width = 0.1, height = 0), size = 2) + 
  
  scale_fill_manual(values = treatment_colors) +  # Custom treatment-specific colors
  scale_color_manual(values = treatment_colors) +  # Match observed points to treatment colors
  
  # Rename facet labels
  facet_wrap(~week, scales = "fixed", ncol = 1, 
             labeller = labeller(week = c("1" = "Week 1", "6" = "Week 6"))) +  
  
  labs(
    title = "",
    x = "Predicted Value for Connectance",
    y = "Density",
    fill = "",   # Remove "Treatment" label in legend
    color = ""   # Remove "Observed" label in legend
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 14),  
    strip.text = element_text(hjust = 0, size = 13, face = "italic", color = "gray30"),  
    legend.position = "none"  # 🔹 Removes legend (treatment label)
  )


plot2_bmrs<-ggplot(posterior_df, aes(x = predicted, fill = treatment)) +
  geom_density(alpha = 0.5) +  # Density plot with transparency
  
  # Add boxplots for observed data
  geom_boxplot(data = observed_dataA, aes(
    x = GqM, 
    y = as.numeric(factor(treatment)) * -1.5,  # Ensure correct spacing per treatment
    fill = treatment, 
    group = interaction(treatment, week)
  ), 
  width = 1.2, alpha = 0.5, outlier.shape = NA) +  # Set width and remove outliers
  
  # Add observed points (jittered)
  geom_point(data = observed_dataA, aes(
    x = GqM, 
    y = as.numeric(factor(treatment)) * -1.5, 
    color = treatment
  ), 
  position = position_jitter(width = 0.1, height = 0), size = 2) + 
  
  scale_fill_manual(values = treatment_colors) +  # Custom treatment-specific colors
  scale_color_manual(values = treatment_colors) +  # Match observed points to treatment colors
  
  # Rename facet labels
  facet_wrap(~week, scales = "fixed", ncol = 1, 
             labeller = labeller(week = c("1" = "Week 1", "6" = "Week 6"))) +  
  
  labs(
    title = "",
    x = "Predicted Value for Generality",
    y = "Density",
    fill = "",   # Remove "Treatment" label in legend
    color = ""   # Remove "Observed" label in legend
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 14),  
    strip.text = element_text(hjust = 0, size = 13, face = "italic", color = "gray30"),  
    legend.position = "none"  # 🔹 Removes legend (treatment label)
  )

### make combined plot for Figure 2 in the ms 
library("cowplot")
plot_grid(plot1_bmrs,plot2_bmrs,
          labels = c("",""),
          label_size = 18, label_fontfamily = NULL, label_fontface = "plain",
          ncol = 1, nrow = 2, align = "v")

#save as 5x10 pdf
################################################################################
################################################################################

####test for differences in individual plasmidgeneralism depending on treatment x time x plasmidID

#read in csv table that has only 1 column for Generalist and included plasmid ID (pkjk5,pb5,pb12)
Pl_m<-Plas_metrics
str(Pl_m)

#create a new random factor that represents each of the microbial communities over time
Pl_m <- mutate(Pl_m, repl = paste(treatment, replicate)) 
# create new factor that contains treatment, week, plasID  information 
Pl_m <- mutate(Pl_m, tr = paste(treatment, week, plasID)) 

#### analysing networks metrics 
Pl_m$treatment = factor(Pl_m$treatment , levels=c("n","w","s"))

######################################################################
# use brms on the same model 
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

# without specified priors 
library(brms)

# Ensure categorical variables
Pl_m$plasID <- as.factor(Pl_m$plasID)
Pl_m$repl <- as.factor(Pl_m$repl)

# Fit the Bayesian model
lme1.Plas <- brm(
  G ~ week * treatment * plasID + (1 | repl),
  data = Pl_m,
  family = gaussian(),
  chains = 4,
  iter = 4000,  # More iterations for better estimates
  warmup = 1000,  # More burn-in to improve convergence
  cores = parallel::detectCores(),
)

summary(lme1.Plas)

plot(lme1.Plas)  # Diagnostic plots
pp_check(lme1.Plas)  # Posterior predictive checks locking good 

loo(lme1.Plas)  # Leave-One-Out Cross-Validation

# add weak priors 

# Ensure categorical variables
Pl_m$plasID <- as.factor(Pl_m$plasID)
Pl_m$repl <- as.factor(Pl_m$repl)

# Define weak priors
priors <- c(
  set_prior("normal(0, 10)", class = "Intercept"),  # Prior for the intercept
  set_prior("normal(0, 5)", class = "b"),          # Priors for fixed effects
  set_prior("student_t(3, 0, 10)", class = "sd"),  # Prior for random effect SD
  set_prior("student_t(3, 0, 10)", class = "sigma") # Prior for residual SD
)

# Fit the Bayesian model with weakly informative priors
lme1.Plas <- brm(
  G ~ week * treatment * plasID + (1 | repl),
  data = Pl_m,
  family = gaussian(),
  prior = priors,  # Adding priors 
  chains = 4,
  iter = 6000,  
  warmup = 1500,
  cores = parallel::detectCores(),
)

#check priors 
prior_summary(lme1.Plas)
#posterior predictive model check 
pp_check(lme1.Plas)

#looks like weak priors improve the model but not by much
loo(lme1.Plas) #good after increasing inter and warmup

#model results 
summary(lme1.Plas)
#plots
plot(conditional_effects(lme1.Plas), points = TRUE)
#as table
posterior_summary(lme1.Plas)

##### plot posterior distribution for single parameters 

library(posterior)

################plot posterior distribution for treatment levels 
# Ensure required libraries are loaded
library(tidyverse)
library(brms)

##observed data 
# dataset (observed_data) that contains the observed values with corresponding week, treatment and PlasID labels
library(dplyr)
observed_dataPlas<-Pl_m  %>%
  select(week, treatment, plasID,G)

######################### predicted values
# Assuming 'new_data_Plas' is created with the relevant combinations of 'PlasID', 'treatment', and 'week'.
new_data_Plas <- expand.grid(
  week = c(1, 6),  # Example for week 1 and 6
  treatment = unique(Pl_m$treatment),  # List of treatments in your data
  plasID = unique(Pl_m$plasID)  # List of PlasID values
)

# Generate posterior predictions (allowing for new levels and including the random effects)
posterior_preds_Plas <- posterior_epred(lme1.Plas, newdata = new_data_Plas, allow_new_levels = TRUE)

# Convert posterior predictions to a tidy format
posterior_df_Plas <- posterior_preds_Plas %>%
  as_tibble() %>%
  mutate(iteration = row_number()) %>%
  pivot_longer(cols = -iteration, names_to = "newdata_index", values_to = "predicted") %>%
  mutate(
    week = rep(new_data_Plas$week, times = nrow(posterior_preds_Plas)),  # Properly repeat for correct matching
    treatment = rep(new_data_Plas$treatment, times = nrow(posterior_preds_Plas)),
    plasID = rep(new_data_Plas$plasID, times = nrow(posterior_preds_Plas))
  )

#####################
###############plot with observed and predicted values for plasID on generality 
# Define treatment-specific colors
treatment_colors <- c("#ABD9E9", "#FDAE61", "#D73027")

ggplot(posterior_df_Plas, aes(x = predicted, fill = treatment, group = interaction(treatment, plasID))) +
  geom_density(alpha = 0.5) +  # Density plot with transparency
  
  # Add boxplots for observed data
  geom_boxplot(data = observed_dataPlas, aes(
    x = G, 
    y = as.numeric(factor(treatment)) * -1.5,  # Ensure correct spacing per treatment
    fill = treatment, 
    group = interaction(treatment, plasID, week)  # Group by plasID too
  ), 
  width = 1.2, alpha = 0.5, outlier.shape = NA) +  # Slightly wider boxplots
  
  # Add observed points (jittered)
  geom_point(data = observed_dataPlas, aes(
    x = G, 
    y = as.numeric(factor(treatment)) * -1.5, 
    color = treatment
  ), 
  position = position_jitter(width = 0.15, height = 0), size = 2) +  # More jitter for clarity
  
  scale_fill_manual(values = treatment_colors) +  # Custom treatment-specific colors
  scale_color_manual(values = treatment_colors) +  # Match observed points to treatment colors
  
  # Facet by week AND PlasID
  facet_grid(week ~ plasID, scales = "fixed", 
             labeller = labeller(week = c("1" = "Week 1", "6" = "Week 6"))) +  
  
  labs(
    title = "",
    x = "Predicted Value for Generality",
    y = "Density",
    fill = "",   # Remove "Treatment" label in legend
    color = ""   # Remove "Observed" label in legend
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 14),  
    strip.text = element_text(hjust = 0, size = 13, face = "italic", color = "gray30"),  
    legend.position = "none"  # 🔹 Removes legend (treatment label)
  )

#################################################################################
#################################################################################
#create table for the model PlasID model outputs
# Get summary of the model (including fixed effects and random effects)
model_summary <- summary(lme1.Plas)

# Extract fixed effects and random effects
fixed_effects <- model_summary$fixed
random_effects <- model_summary$random

# Print the fixed effects
fixed_effects

library(broom.mixed)

# Tidy the model output for easy use in tables
tidy_results <- tidy(lme1.Plas)
print(tidy_results)

library(knitr)

# Create a summary table with kable
kable(tidy_results, format = "markdown", caption = "Summary of Fixed Effects and Random Effects")
# This will print a clean, tabular format in the console
print(kable(tidy_results, caption = "Summary of Fixed Effects and Random Effects"))

# Load required libraries
library(officer)
library(flextable)

# Assuming 'tidy_results' is your model output (from brms or lme4 model)
# Create a Word document
doc <- read_docx()

# Add a title to the document
doc <- doc %>% 
  body_add_par("Model Summary Table", style = "heading 1") %>% 
  body_add_par("Summary of Fixed Effects and Random Effects", style = "heading 2")

# Convert the tidy_results to a flextable (formatted table)
flextable_results <- flextable(tidy_results)

# Add the table to the Word document
doc <- doc %>% body_add_flextable(flextable_results)

# Save the document to a Word file
print(doc, target = "model_summary_table.docx")



#####code for plot
# Define treatment-specific colors
treatment_colors <- c("#ABD9E9", "#FDAE61", "#D73027")

ggplot(posterior_df_Plas, aes(x = predicted, fill = treatment)) +
  geom_density(alpha = 0.5) +  # Density plot with transparency
  
  # Add boxplots for observed data
  geom_boxplot(data = observed_dataPlas, aes(
    x = G, 
    y = as.numeric(factor(treatment)) * -1.5,  # Ensure correct spacing per treatment
    fill = treatment, 
    group = interaction(treatment, week, plasID)  # Account for plasID
  ), 
  width = 1.2, alpha = 0.5, outlier.shape = NA) +  # Slightly wider boxplots
  
  # Add observed points (jittered)
  geom_point(data = observed_dataPlas, aes(
    x = G, 
    y = as.numeric(factor(treatment)) * -1.5, 
    color = treatment
  ), 
  position = position_jitter(width = 0.15, height = 0), size = 2) +  # More jitter for clarity
  
  scale_fill_manual(values = treatment_colors) +  # Custom treatment-specific colors
  scale_color_manual(values = treatment_colors) +  # Match observed points to treatment colors
  
  # Facet by week
  facet_wrap(~week, scales = "fixed", ncol = 1, 
             labeller = labeller(week = c("1" = "Week 1", "6" = "Week 6"))) +  
  
  labs(
    title = "",
    x = "Predicted Value for Generality",
    y = "Density",
    fill = "",   # Remove "Treatment" label in legend
    color = ""   # Remove "Observed" label in legend
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 14),  
    strip.text = element_text(hjust = 0, size = 13, face = "italic", color = "gray30"),  
    legend.position = "none"  # 🔹 Removes legend (treatment label)
  )

