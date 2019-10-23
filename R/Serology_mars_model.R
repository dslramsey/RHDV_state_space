library(tidyverse)
library(lubridate)
library(gridExtra)
library(rstan)
library(MCMCvis)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# Read Data ---------------

sera<- read_csv("./Data/serology_data.csv",na=c("","NA"))

time_periods = sera %>% group_by(Site_code) %>% 
  summarise(start_time=min(time_ind), end_time=max(time_ind), n=n())

#--------------------------------
av<- function(x){as.vector(x)}

i_rcva<- which(sera$n_RCVA > 0)
y_rcva<- sera$RCVA[i_rcva]
n_rcva<- sera$n_RCVA[i_rcva]

i_rhd1<- which(sera$n_RHD1 > 0)
y_rhd1<- sera$RHD1[i_rhd1]
n_rhd1<- sera$n_RHD1[i_rhd1]

i_rhd2<- which(sera$n_RHD2 > 0)
y_rhd2<- sera$RHD2[i_rhd2]
n_rhd2<- sera$n_RHD2[i_rhd2]

nsites=max(sera$Site_code)
season<- av(sera$Season_code)

site=av(sera$Site_code)

start_ind<- c(1,cumsum(time_periods$n)[-nsites] + 1)
end_ind<- cumsum(time_periods$n)

N=nrow(sera)
N_rcva<- length(i_rcva)
N_rhd1<- length(i_rhd1)
N_rhd2<- length(i_rhd2)

Data<- list(y_rcva=y_rcva,n_rcva=n_rcva,N_rcva=N_rcva,i_rcva=i_rcva,
            y_rhd1=y_rhd1,n_rhd1=n_rhd1,N_rhd1=N_rhd1,i_rhd1=i_rhd1,
            y_rhd2=y_rhd2,n_rhd2=n_rhd2,N_rhd2=N_rhd2,i_rhd2=i_rhd2,
            N=N,nsites=nsites,season=season,start=start_ind,end=end_ind)
            
            
# MCMC settings
ni <- 4000
nt <- 1
nb <- 2000
nc <- 5


params <- c("prev_rcva","prev_rhd1","prev_rhd2","site",
            "beta","eta","sigma_proc","Sigma")
            
set.seed(123) # for reproducability

out <- stan("./stan/serology_model_mars.stan",
            data = Data,
            init = "0", pars = params,
            chains = nc, iter = ni, warmup = nb, thin = nt, refresh=50)


MCMCsummary(out, params=c("beta","eta","sigma_proc"), n.eff=TRUE)

MCMCsummary(out, params=c("Sigma"))

#-----------------------------------------------------------------

preds_rcva<- as_tibble(MCMCsummary(out,"prev_rcva")) 
preds_rcva<- preds_rcva %>% select("mean", q5="2.5%", q95="97.5%")
preds_rhd1<- as_tibble(MCMCsummary(out,"prev_rhd1")) 
preds_rhd1<- preds_rhd1 %>% select("mean", q5="2.5%", q95="97.5%")
preds_rhd2<- as_tibble(MCMCsummary(out,"prev_rhd2")) 
preds_rhd2<- preds_rhd2 %>% select("mean", q5="2.5%", q95="97.5%")

dat<- sera %>% mutate(Prevalence = RCVA/n_RCVA) %>% select(Site, time_ind, Prevalence)
preds_rcva<- bind_cols(dat, preds_rcva) 
dat<- sera %>% mutate(Prevalence = RHD1/n_RHD1) %>% select(Site, time_ind, Prevalence)
preds_rhd1<- bind_cols(dat, preds_rhd1) 
dat<- sera %>% mutate(Prevalence = RHD2/n_RHD2) %>% select(Site, time_ind, Prevalence)
preds_rhd2<- bind_cols(dat, preds_rhd2) 

preds_rcva<- mutate(preds_rcva, Strain="RCVA") 
preds_rhd1<- mutate(preds_rhd1, Strain="RHDV") 
preds_rhd2<- mutate(preds_rhd2, Strain="RHDV2")

preds<- bind_rows(preds_rcva, preds_rhd1, preds_rhd2)

pp<- sort(unique(sera$period))
plabs<- round(pp[seq(4,length(pp),4)])
mint<- min(preds$time_ind)
maxt<- max(preds$time_ind)

#--------------
# All strains
win.graph(12,12)
preds %>% 
  ggplot(aes(x=time_ind, y=mean, color=Strain))+
  geom_line(size=0.75)+
  scale_colour_manual(values=c("grey60","black","orange")) +
  ylab("Prevalence")+
  xlab("Time")+
  scale_y_continuous(breaks=c(0, 0.25, 0.5, 0.75, 1), lim=c(0, 1)) +
  scale_x_continuous(breaks=seq(4,28,4), minor_breaks=seq(mint, maxt, 1), 
                     labels=plabs, lim=c(mint, maxt) ) +
  facet_wrap(~ Site, nrow=6) +
  theme_bw()+
  theme(strip.background = element_blank(), 
        strip.text.x=element_text(hjust=0.05),
        panel.border = element_rect(colour = "black"),
        panel.grid.major=element_line(colour="grey50", size=0.15),
        panel.grid.minor=element_line(colour="grey90", size=0.05),
        axis.text.x=element_text(size=8))

#--------------
win.graph(10,12)
preds %>% filter(Strain %in% c("RHDV2")) %>% 
  ggplot(aes(x=time_ind, y=mean))+
  geom_line(size=0.7)+
  geom_ribbon(aes(ymin=q5,ymax=q95),alpha=0.1)+
  ylab("Prevalence")+
  xlab("Time")+
  scale_y_continuous(breaks=c(0, 0.25, 0.5, 0.75, 1), lim=c(0, 1)) +
  scale_x_continuous(breaks=seq(4,28,4), minor_breaks=seq(mint, maxt, 1), 
                     labels=plabs, lim=c(mint, maxt) ) +
  geom_point(aes(x=time_ind, y=Prevalence), cex=1.5) +
  facet_wrap(~ Site,  nrow=6) +
  theme_bw()+
  theme(strip.background = element_blank(), 
        strip.text.x=element_text(hjust=0.05),
        panel.border = element_rect(colour = "black"),
        panel.grid.major=element_line(colour="grey50", size=0.15),
        panel.grid.minor=element_line(colour="grey90", size=0.05),
        axis.text.x=element_text(size=8),
        legend.position="none")

