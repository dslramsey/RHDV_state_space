library(tidyverse)
library(lubridate)
library(gridExtra)
library(rstan)
library(MCMCvis)
library(RColorBrewer)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# Read Data ---------------

sera_age<- read_csv("./Data/serology_age_data.csv",na=c("","NA"))


time_periods = sera_age %>% group_by(Site_code) %>% 
  summarise(start_time=min(Time), end_time=max(Time), n=end_time-start_time+1)

#---------------------------------------------------------
av<- function(x){as.vector(x)}

y_rcva_j<- sera_age$RCVA_juvenile
n_rcva_j<- sera_age$n_RCVA_juvenile

y_rcva_a<- sera_age$RCVA_adult
n_rcva_a<- sera_age$n_RCVA_adult

y_rhd1_j<- sera_age$RHD1_juvenile
n_rhd1_j<- sera_age$n_RHD1_juvenile

y_rhd1_a<- sera_age$RHD1_adult
n_rhd1_a<- sera_age$n_RHD1_adult

y_rhd2_j<- sera_age$RHD2_juvenile
n_rhd2_j<- sera_age$n_RHD2_juvenile

y_rhd2_a<- sera_age$RHD2_adult
n_rhd2_a<- sera_age$n_RHD2_adult

nsites=max(sera_age$Site_code)

site=av(sera_age$Site_code)
Time<- av(sera_age$Time)

start_ind<- c(1,cumsum(time_periods$n)[-nsites] + 1)
end_ind<- cumsum(time_periods$n)

N<- nrow(sera_age)
NT<- max(Time)


#======================================================================================
# Mean and site age-specific trends
#======================================================================================
data<- list(y_rcva_j=y_rcva_j,n_rcva_j=n_rcva_j,
            y_rcva_a=y_rcva_a,n_rcva_a=n_rcva_a,
            y_rhd1_j=y_rhd1_j,n_rhd1_j=n_rhd1_j,
            y_rhd1_a=y_rhd1_a,n_rhd1_a=n_rhd1_a,
            y_rhd2_j=y_rhd2_j,n_rhd2_j=n_rhd2_j,
            y_rhd2_a=y_rhd2_a,n_rhd2_a=n_rhd2_a,
            N=N,NT=NT,time_ind=Time,nsites=nsites,site_ind=site)            


# MCMC settings
ni <- 3000
nt <- 1
nb <- 1000
nc <- 5


params <- c("prev_rcva_j","prev_rcva_a","prev_rhd1_j","prev_rhd1_a","prev_rhd2_j",
             "prev_rhd2_a","A","sigma_proc","sigma_site","mean_rcva_j","mean_rcva_a",
             "mean_rhd1_j","mean_rhd1_a","mean_rhd2_j","mean_rhd2_a")


## Call Stan from R
out <- stan("./stan/serology_model_trends.stan",
            data = data,
            init = "0", pars = params,
            chains = nc, iter = ni, warmup = nb, thin = nt, refresh=50,
            control=list(adapt_delta=0.95))


MCMCsummary(out, params=c("A","sigma_proc","sigma_site"), n.eff=TRUE)

#-----------------------------------------------------------------
# Average age-specific trends
#-----------------------------------------------------------------
pd <- position_dodge(0.2) # move them .05 to the left and right
pal<- rev(brewer.pal(9, "Oranges"))


preds_rcva_j<- as_tibble(MCMCsummary(out,"mean_rcva_j",HPD=TRUE, hpd_prob=0.89)) 
preds_rcva_j<- preds_rcva_j %>% select(c(1,3,4)) %>% 
  mutate(age_class="juvenile",strain="RCVA", Time=1:14)    

preds_rcva_a<- as_tibble(MCMCsummary(out,"mean_rcva_a",HPD=TRUE, hpd_prob=0.89)) 
preds_rcva_a<- preds_rcva_a %>% select(c(1,3,4)) %>% 
  mutate(age_class="adult",strain="RCVA",Time=1:14)       

preds_rhd1_j<- as_tibble(MCMCsummary(out,"mean_rhd1_j",HPD=TRUE, hpd_prob=0.89)) 
preds_rhd1_j<- preds_rhd1_j %>% select(c(1,3,4)) %>% 
  mutate(age_class="juvenile",strain="RHDV", Time=1:14)

preds_rhd1_a<- as_tibble(MCMCsummary(out,"mean_rhd1_a",HPD=TRUE, hpd_prob=0.89)) 
preds_rhd1_a<- preds_rhd1_a %>% select(c(1,3,4)) %>% 
  mutate(age_class="adult",strain="RHDV",Time=1:14)    

preds_rhd2_j<- as_tibble(MCMCsummary(out,"mean_rhd2_j",HPD=TRUE, hpd_prob=0.89)) 
preds_rhd2_j<- preds_rhd2_j %>% select(c(1,3,4)) %>% 
  mutate(age_class="juvenile",strain="RHDV2",Time=1:14)     

preds_rhd2_a<- as_tibble(MCMCsummary(out,"mean_rhd2_a",HPD=TRUE, hpd_prob=0.89)) 
preds_rhd2_a<- preds_rhd2_a %>% select(c(1,3,4)) %>% 
  mutate(age_class="adult",strain="RHDV2",Time=1:14)     

preds<- bind_rows(preds_rcva_j, preds_rcva_a, preds_rhd1_j, preds_rhd1_a, preds_rhd2_j,preds_rhd2_a)
preds<- preds %>% mutate(age_class=factor(age_class, levels=c("juvenile","adult")))
names(preds)[2:3]<- c("lcl","ucl")


win.graph(10,5)
preds %>% ggplot(aes(Time, mean, color=fct_rev(age_class))) +
  geom_errorbar(aes(ymin=lcl,ymax=ucl), width=0.21, position=pd) +  
  geom_line(size=1.5,position=pd) +
  geom_point(size=3, position=pd, shape=21, fill="white") +
  facet_wrap(~ strain) +
  scale_color_manual(values = pal[c(2,6)]) +
  ylim(0, 1) +
  scale_x_continuous(breaks=1:14,labels=as.character(0:(14-1))) +
  annotate("rect", xmin=1.5, xmax=5.5, ymin=0, ymax=1, alpha=0.2) +
  annotate("rect", xmin=9.5, xmax=13.5, ymin=0, ymax=1, alpha=0.2) +
  labs(x="Time since RHDV2 arrival (quarters)",y="Seroprevalence",color="Age Class") +
  theme_bw() +
  theme(axis.text.x = element_text(size=10),
        axis.title = element_text(face="bold",size=12),
        strip.background = element_rect(fill="white"))

#-----------------------------------------------------------------
# Site and age-specific trends
#-----------------------------------------------------------------
time.mat<- select(sera_age, Time, Site_code, Site)

preds_rcva_j<- as_tibble(MCMCsummary(out,"prev_rcva_j",HPD=TRUE, hpd_prob=0.89)) 
preds_rcva_j<- preds_rcva_j %>% select(c(1,3,4)) %>% 
  mutate(age_class="juvenile",strain="RCVA") %>% bind_cols(time.mat)    

preds_rcva_a<- as_tibble(MCMCsummary(out,"prev_rcva_a",HPD=TRUE, hpd_prob=0.89)) 
preds_rcva_a<- preds_rcva_a %>% select(c(1,3,4)) %>% 
  mutate(age_class="adult",strain="RCVA") %>% bind_cols(time.mat)      

preds_rhd1_j<- as_tibble(MCMCsummary(out,"prev_rhd1_j",HPD=TRUE, hpd_prob=0.89)) 
preds_rhd1_j<- preds_rhd1_j %>% select(c(1,3,4)) %>% 
  mutate(age_class="juvenile",strain="RHDV") %>% bind_cols(time.mat)    
preds_rhd1_a<- as_tibble(MCMCsummary(out,"prev_rhd1_a",HPD=TRUE, hpd_prob=0.89)) 
preds_rhd1_a<- preds_rhd1_a %>% select(c(1,3,4)) %>% 
  mutate(age_class="adult",strain="RHDV") %>% bind_cols(time.mat)    

preds_rhd2_j<- as_tibble(MCMCsummary(out,"prev_rhd2_j",HPD=TRUE, hpd_prob=0.89)) 
preds_rhd2_j<- preds_rhd2_j %>% select(c(1,3,4)) %>% 
  mutate(age_class="juvenile",strain="RHDV2") %>% bind_cols(time.mat)     
preds_rhd2_a<- as_tibble(MCMCsummary(out,"prev_rhd2_a",HPD=TRUE, hpd_prob=0.89)) 
preds_rhd2_a<- preds_rhd2_a %>% select(c(1,3,4)) %>% 
  mutate(age_class="adult",strain="RHDV2") %>% bind_cols(time.mat)     

preds<- bind_rows(preds_rcva_j, preds_rcva_a, preds_rhd1_j, preds_rhd1_a, preds_rhd2_j,preds_rhd2_a)
preds<- preds %>% mutate(age_class=factor(age_class, levels=c("juvenile","adult")))
names(preds)[2:3]<- c("lcl","ucl")


win.graph(10,10)
preds %>% filter(strain %in% c("RHDV")) %>% 
  ggplot(aes(Time,mean, color=fct_rev(age_class))) +
  geom_errorbar(aes(ymin=lcl,ymax=ucl), width=0.21, position=pd) +  
  geom_line(size=1,position=pd) +
  geom_point(size=1.5, position=pd, shape=21, fill="white") +
  facet_wrap(~ Site,  nrow=6) +
  scale_color_manual(values = pal[c(2,6)]) +
  ylim(0, 1) +
  scale_x_continuous(breaks=1:N,labels=as.character(0:(N-1))) +
  annotate("rect", xmin=1.5, xmax=5.5, ymin=0, ymax=1, alpha=0.2) +
  annotate("rect", xmin=9.5, xmax=13.5, ymin=0, ymax=1, alpha=0.2) +
  labs(x="Time since RHDV2 arrival (quarters)",y="Seroprevalence",color="Age Class") +
  theme_bw() +
  theme(strip.background = element_blank(), 
        strip.text.x=element_text(hjust=0.05),
        panel.border = element_rect(colour = "black"),
        panel.grid.major=element_line(colour="grey90", size=0.15),
        panel.grid.minor=element_line(colour="grey90", size=0.05),
        axis.text.x=element_text(size=9),
        axis.title = element_text(face="bold",size=12))
        
