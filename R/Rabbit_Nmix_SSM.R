library(tidyverse)
library(lubridate)
library(gridExtra)
library(rstan)
library(MCMCvis)
library(gridExtra)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

#---------------------------------
spot<- read_csv("./Data/spotlight_data.csv",na=c("","NA"))


time_periods = spot %>% group_by(Site,Site_code) %>% summarise(start_time=min(time_ind),
                                                          end_time=max(time_ind),
                                                          n=n())

K5<- spot %>% group_by(Site, Site_code) %>% summarise(K5_site =  unique(K5_site))
#----------------------------------

row.sort<- function(x, y) {
  # function to sort rows of matrix 'x' 
  # according to rank 'y'  
  n<- dim(x)[1]
  res<- list()
  for(i in 1:n) {
    res[[i]]<- x[i,y[i,]] 
  }
  do.call('rbind',res)
}

#------------------------------------

y = as.matrix(select(spot,Night1,Night2,Night3))
rank<- t(apply(y,1,order,decreasing=T))
y<- row.sort(y,rank)            
ncounts<- apply(y,1, function(x) length(x[!is.na(x)]))
y[is.na(y)]<- -1

nsites=max(spot$Site_code)
season<- spot$Season_code

rhdv2<- spot$rhdv2
K5_rel<- spot$K5_rel
K5_site<- K5$K5_site

start_ind<- c(1,cumsum(time_periods$n)[-nsites] + 1)
end_ind<- cumsum(time_periods$n)

trans_length<- spot$Distance
TL<- log(trans_length)

site<- as.vector(spot$Site_code)

N=nrow(spot)
max_y<- apply(y, 1, max, na.rm=T)


Data<- list(y=y,N=N,ncounts=ncounts,nsites=nsites,site=site,season=season,rhdv2=rhdv2,
            Rel=K5_rel,K5=K5_site,start=start_ind,end=end_ind,TL=TL,max_y=max_y)


# MCMC settings
# ** short run:  Longer warmup and iters needed for inference **
ni <- 1000  
nt <- 1
nb <- 500
nc <- 1

params <- c("mu_gam","sigma_gam","mu_b","sigma_b","mu_k5","sigma_k5","beta","eta",
            "gam","delta","sigma_site","sigma_proc","r_mean","site_r","mu_rabbits","p",
            "EQB","EQA","DE","mu_p","sigma_p")

## Call Stan from R
out <- stan("./stan/rabbit_Nmix_model.stan",
            data = Data,
            init = "0", pars = params,
            chains = nc, iter = ni, warmup = nb, thin = nt, refresh=50)


MCMCsummary(out, params=c("gam","eta","delta", "sigma_site","sigma_proc"))
            
MCMCsummary(out, params=c("r_mean","site_r","mu_b"))

MCMCsummary(out, params=c("mu_p","sigma_p","p"))

#=======================================================
# Rabbit trends at each site
#=======================================================

tindex<- spot %>% mutate(rn=row_number()) %>% select(rn,Site_code,time_ind,Season_code)

preds<- MCMCsummary(out,"mu_rabbits")
preds<- rename(preds, q5='2.5%',q50='50%', q95='97.5%')

preds<- preds %>% mutate(mean=exp(mean), q5=exp(q5), q50=exp(q50), q95=exp(q95),
                         Site=tindex$Site_code, Time=tindex$time_ind, 
                         Season=factor(tindex$Season_code),
                                       Site=factor(Site,labels=time_periods$Site))
                    

summ<- spot %>% group_by(Site,Time=time_ind) %>% 
  summarise(Index=mean(c(Night1,Night2,Night3),na.rm=T)/Distance) %>% ungroup() 


rhd_arr<- spot %>% group_by(Site) %>% summarise(RHD2=min(time_ind[rhdv2==1]))
                                                
K5_arr<- K5 %>% mutate(K5_rel = if_else(K5_site==2,25,-1))                                                


preds<- left_join(preds, summ)

pp<- sort(unique(spot$period))
pp<- round(pp[seq(4,28,4)])
xxx<- data.frame(start=seq(2, 28, 4), end=seq(3, 29, 4),ymn=0.1,ymx=1500)

win.graph(12,12)
preds %>% 
  ggplot(aes(x=Time, y=mean))+
  geom_line(colour="black", size=1)+
  geom_ribbon(aes(ymin=q5,ymax=q95, colour=NULL),alpha=0.25, fill="grey40")+
  ylab(expression(paste(log[10],"(Rabbits k",m^{-1},")")))+
  xlab("Time")+
  scale_y_log10(breaks=c(1, 10, 100, 1000), lim=c(0.1, 1500.1))+
  scale_x_continuous(breaks=seq(4, 28, 4), minor_breaks=seq(1, 28, 1), 
                     labels=pp, lim=c(1, 28) )+
  geom_point(aes(x=Time, y=Index), cex=1, shape=1) +
  facet_wrap(~ Site,  nrow=6) +
  geom_vline(aes(xintercept = RHD2), rhd_arr, colour="red", linetype=2, size=1) +
  geom_vline(aes(xintercept = K5_rel), K5_arr, colour="blue", linetype=2, size=1) +
  theme_bw()+
  theme(strip.background = element_blank(),
        strip.text.x=element_text(hjust=0.05),
        panel.border = element_rect(colour = "black"),
        panel.grid.major=element_line(colour="grey70", size=0.15),
        panel.grid.minor=element_line(colour="grey90", size=0.05),
        axis.text.x=element_text(size=8))

#=======================================================
# RHDV2 impact
#=======================================================
preds<- MCMCsummary(out,c("gam","mu_gam"))
preds<- rename(preds, q5='2.5%',q50='50%', q95='97.5%')

pred0<- preds["mu_gam",]
predn<- bind_cols(preds[1:nsites,], rhd_arr)

predn<- predn %>% mutate(Sig=as.factor(if_else(q95 < 0 ,1, 0)),
                         Site=lvls_reorder(Site, c(17,18,1,9,14,16,2,4,7,10,11,13,3,8,5,15,6,12)))

win.graph(10,5)
p1<- predn  %>% 
  ggplot(aes(x = Site, y = mean, colour=Sig)) +
  geom_crossbar(aes(ymin=q5, ymax=q95), width=0.5) +
  coord_cartesian(ylim = c(-1.5, 1.5)) +
  guides(colour=FALSE) +
  scale_colour_manual(values=c("black","red")) +
  scale_x_discrete() +
  geom_hline(yintercept = 0, linetype=2) +
  geom_vline(xintercept=c(2.5,6.5,12.5,14.5,16.5),linetype=1,colour="blue") +
  labs(x="Site",y="Change in growth rate") +
  annotate("text", x=c(1.5,4.5,9.5,13.5,15.5,17.5),y=rep(1.5,6),
           label=c("Qld","Vic","NSW","ACT","SA","WA")) +
  annotate("text", x=c(3,6,11,12,15,16,17,18),y=rep(-1.5,8),
           label=c("*")) +
  theme_bw()+
  theme(strip.background = element_blank(), 
        panel.border = element_rect(colour = "black"),
        panel.grid.major=element_line(colour="grey90", size=0.15),
        panel.grid.minor=element_line(colour="grey90", size=0.05),
        axis.text.x=element_text(size=10, angle=45, hjust=1),
        axis.title = element_text(face="bold",size=14))


p2<- pred0 %>% 
  ggplot(aes(x="Overall",y=mean)) +
  geom_crossbar(aes(ymin=q5, ymax=q95),width=0.25) +
  geom_hline(yintercept = 0, linetype=2) +
  coord_cartesian(ylim = c(-1, 1)) +
  labs(x="Overall",y="Change in growth rate") +
  theme_bw()+
  theme(strip.background = element_blank(), 
        panel.border = element_rect(colour = "black"),
        panel.grid.major=element_line(colour="grey90", size=0.15),
        panel.grid.minor=element_line(colour="grey90", size=0.05),
        axis.text.x=element_blank(),
        axis.title = element_text(face="bold",size=12))

grid.arrange(p1, p2, ncol=2, widths=c(0.8,0.2))

#=======================================================
# K5 impact
#=======================================================
preds<- MCMCsummary(out,"delta")
preds<- rename(preds, q5='2.5%',q50='50%', q95='97.5%')

predn<- bind_cols(preds, K5_arr)

predn<- predn %>% mutate(K5=if_else(K5_site == 2,"Release","Non-release"),
                         Site=lvls_reorder(Site, c(18,17,16,1,9,14,4,7,10,2,13,11,8,3,5,15,6,12)))

preds<- MCMCsummary(out,"mu_k5")
preds<- rename(preds, q5='2.5%',q50='50%', q95='97.5%')
pred0<- mutate(preds, Release=c("Non-release","Release"))


win.graph(10,5)
p1<- predn  %>% 
  ggplot(aes(x = Site, y = mean, colour=K5)) +
  geom_crossbar(aes(ymin=q5, ymax=q95), width=0.5) +
  coord_cartesian(ylim = c(-1, 1)) +
  guides(colour=FALSE) +
  scale_colour_manual(values=c("black","red")) +
  scale_x_discrete() +
  geom_hline(yintercept = 0, linetype=2) +
  geom_vline(xintercept=c(2.5,6.5,12.5,14.5,16.5),linetype=1,colour="blue") +
  labs(x="Site",y="Change in growth rate") +
  annotate("text", x=c(1.5,4.5,9.5,13.5,15.5,17.5),y=rep(1,6),
           label=c("Qld","Vic","NSW","ACT","SA","WA")) +
  theme_bw()+
  theme(strip.background = element_blank(), 
        panel.border = element_rect(colour = "black"),
        panel.grid.major=element_line(colour="grey90", size=0.15),
        panel.grid.minor=element_line(colour="grey90", size=0.05),
        axis.text.x=element_text(size=10, angle=45, hjust=1),
        axis.text.y=element_text(size=10),
        axis.title = element_text(face="bold",size=14))


p2<- pred0 %>% 
  ggplot(aes(x=Release,y=mean,colour=Release)) +
  geom_crossbar(aes(ymin=q5, ymax=q95),width=0.5) +
  geom_hline(yintercept = 0, linetype=2) +
  guides(colour=FALSE) +
  scale_colour_manual(values=c("black","red")) +
  coord_cartesian(ylim = c(-1, 1)) +
  labs(x="Overall",y="Change in growth rate") +
  theme_bw()+
  theme(strip.background = element_blank(), 
        panel.border = element_rect(colour = "black"),
        panel.grid.major=element_line(colour="grey90", size=0.15),
        panel.grid.minor=element_line(colour="grey90", size=0.05),
        axis.text.x=element_text(size=10),
        axis.text.y=element_text(size=10),
        axis.title = element_text(face="bold",size=12))

grid.arrange(p1, p2, ncol=2, widths=c(0.7,0.3))

