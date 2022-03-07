## Signaling climate resilience to municipal bond markets: 
## Does membership in adaptation-focused voluntary clubs affect bond rating?
## Accepted in Climatic Change
## Final check: Mar 4, 2022

## Inhwan Ko and Aseem Prakash
## Department of Political Science, University of Washington

#rm(list=ls())
# always run this first when reproducing (including set.seed)
library(Amelia);library(simcf);library(tile);library(MASS);library(plm);library(tidyverse);library(ggpubr);library(ggplot2)

set.seed(2020)

data <- read_csv("data.csv") # load dataset
unique(data$city[data$rc100city==1]) # inspect the number of 100rc cities
colnames(data)
data_ame <- data %>% # selecting variables to be used for estimating models
  dplyr::select(city, year, bond1_score, bond2_score, bond3_score,
                population, pcincome, homeown, unemp,
                fema_state, fema_state_fire, fema_county, fema_county_fire) %>% 
                # adding wildfire only as more variables will not make imputation converge
  pdata.frame(index=c("city","year")) %>% mutate(year=as.integer(year)) 

m <- 100

colnames(data_ame)
# this takes a long time!
# we would recommend loading the Rds file directly to replicate
# last checked with November 22, 2021

# a.out <- amelia( # multiple imputation
#  data_ame, cs=1, ts=2, polytime=3, intercs = F,
#  ords=c(3:5), logs=c(6:7), lags=c(3:5), leads=c(3:5), m=m)

a.df <- list(NULL)
for (i in 1:m) { # calculating standardized average scores 
  a.df[[i]] <- a.out$imputations[[i]]
  
  bar1 <- mean(a.df[[i]]$bond1_score)
  sd1 <- sd(a.df[[i]]$bond1_score)
  bar2 <- mean(a.df[[i]]$bond2_score)
  sd2 <- sd(a.df[[i]]$bond2_score)
  bar3 <- mean(a.df[[i]]$bond3_score)
  sd3 <- sd(a.df[[i]]$bond3_score)
  
  a.df[[i]]$bond1_zscore <- (a.df[[i]]$bond1_score-bar1)/sd1
  a.df[[i]]$bond2_zscore <- (a.df[[i]]$bond2_score-bar2)/sd2
  a.df[[i]]$bond3_zscore <- (a.df[[i]]$bond3_score-bar3)/sd3
  a.df[[i]]$bond_zscore <- (a.df[[i]]$bond1_zscore 
                            + a.df[[i]]$bond2_zscore 
                            + a.df[[i]]$bond3_zscore) /3
}

for (i in 1:m) { # adding variables from original dataset to 100 imputed datasets
  a.df[[i]] <- a.df[[i]] %>%
    mutate(iclei=data$iclei, rc100=data$rc100, c40=data$c40, statecode=data$statecode) %>%
    group_by(city) %>% 
    mutate(year=c(1995:2018)) %>% 
    dplyr::select(city, year, statecode, bond1_score, bond2_score, bond3_score,
                  bond_zscore, iclei, c40, rc100,
                  population, pcincome, homeown, unemp, 
                  fema_state, fema_state_fire, fema_county, fema_county_fire) 
  a.df[[i]]$rc100city <- data$rc100city # 1 if a city has been in RC100, 0 otherwise
}

for (i in 1:m) { # adding variables from original dataset to 100 imputed datasets
  a.df[[i]]$fema_county_hur <- data$fema_county_hur
  a.df[[i]]$fema_county_snow <- data$fema_county_snow
  a.df[[i]]$fema_county_flood <- data$fema_county_flood
}

all_a.df <- NULL
for (i in 1:m) {
 all_a.df <- rbind(all_a.df, a.df[[i]])
}
 
max(data$bond1_score, na.rm=T)
max(data$bond2_score, na.rm=T)
max(data$bond3_score, na.rm=T)

####################### Original Model (Figure 1) ###########################
### Note: this model is replaced by the model in line 166 as per reviewers' suggestion
### which includes all control variables for disaster
### We leave this model result for transparency
#############################################################################

formula1a <- bond_zscore ~ lag(bond_zscore) + rc100 + c40 + iclei +
  log(population) + log(pcincome) + homeown + unemp + fema_county_fire

model1a <- simbetas1a <- simbeta1a <- NULL
for (i in 1:m) { # create a matrix that has 10000 simulated coefficients 
  # for all explanatory variables
  model1a <- plm(formula1a, effect="twoway", model="within", data=a.df[[i]])
  simbeta1a <- mvrnorm(10000/m, coefficients(model1a), vcovSCC(model1a, group="statecode"))
  simbetas1a <- rbind(simbetas1a, simbeta1a)
}

simbetas1a <- as.data.frame(simbetas1a)

##################### sensitivity check #################################

########### S&P measure as DV ###############
formula2a <- bond1_score ~ lag(bond1_score) + rc100 + c40 + iclei +
  log(population) + log(pcincome) + homeown + unemp + fema_county_fire
model2a <- simbetas2a <- simbeta2a <- NULL
for (i in 1:m) { # create a matrix that has 10000 simulated coefficients 
  # for all explanatory variables
  model2a <- plm(formula2a, effect="twoway", model="within", data=a.df[[i]])
  simbeta2a <- mvrnorm(10000/m, coefficients(model2a), vcovSCC(model2a, group="statecode"))
  simbetas2a <- rbind(simbetas2a, simbeta2a)
}
simbetas2a <- as.data.frame(simbetas2a)
quantile(simbetas2a[,2], c(0.025, 0.975))
mean(simbetas2a[,2])
mean(simbetas2a[,2])/sd(simbetas2a[,2])
mean(simbetas2a[,3])/sd(simbetas2a[,3])
mean(simbetas2a[,4])/sd(simbetas2a[,4])

########### Moody's measure as DV ###############
formula3a <- bond2_score ~ lag(bond2_score) + rc100 + c40 + iclei +
  log(population) + log(pcincome) + homeown + unemp + fema_county_fire

model3a <- simbetas3a <- simbeta3a <- NULL
for (i in 1:m) { # create a matrix that has 10000 simulated coefficients 
  # for all explanatory variables
  model3a <- plm(formula3a, effect="twoway", model="within", data=a.df[[i]])
  simbeta3a <- mvrnorm(10000/m, coefficients(model3a), vcovSCC(model3a, group="statecode"))
  simbetas3a <- rbind(simbetas3a, simbeta3a)
}
simbetas3a <- as.data.frame(simbetas3a)
quantile(simbetas3a[,2], c(0.025, 0.975))
mean(simbetas3a[,2])
mean(simbetas3a[,2])/sd(simbetas3a[,2])
mean(simbetas3a[,3])/sd(simbetas3a[,3])
mean(simbetas3a[,4])/sd(simbetas3a[,4])

########### Fitch measure as DV ###############
formula4a <- bond3_score ~ lag(bond3_score) + rc100 + c40 + iclei +
  log(population) + log(pcincome) + homeown + unemp + fema_county_fire

model4a <- simbetas4a <- simbeta4a <- NULL
for (i in 1:m) { # create a matrix that has 10000 simulated coefficients 
  # for all explanatory variables
  model4a <- plm(formula4a, effect="twoway", model="within", data=a.df[[i]])
  simbeta4a <- mvrnorm(10000/m, coefficients(model4a), vcovSCC(model4a, group="statecode"))
  simbetas4a <- rbind(simbetas4a, simbeta4a)
}
simbetas4a <- as.data.frame(simbetas4a)
quantile(simbetas4a[,2], c(0.025, 0.975))
mean(simbetas4a[,2])
mean(simbetas4a[,2])/sd(simbetas4a[,2])
mean(simbetas4a[,3])/sd(simbetas4a[,3])
mean(simbetas4a[,4])/sd(simbetas4a[,4])

########### County-level disaster: all measures ###############
formula1b <- bond_zscore ~ lag(bond_zscore) + rc100 + c40 + iclei +
  log(population) + log(pcincome) + homeown + unemp + 
  fema_county_fire + fema_county_hur + fema_county_snow + fema_county_flood

##############################################################################
## the reviewers suggested we use this model as a main model (Feb 7, 2022): 
## therefore, we used the figure 1 to incorporate this result instead:
##############################################################################

model1b <- simbetas1b <- simbeta1b <- NULL
for (i in 1:m) { # create a matrix that has 10000 simulated coefficients 
  # for all explanatory variables
  model1b <- plm(formula1b, effect="twoway", model="within", data=a.df[[i]])
  simbeta1b <- mvrnorm(10000/m, coefficients(model1b), vcovSCC(model1b, group="statecode"))
  simbetas1b <- rbind(simbetas1b, simbeta1b)
}
simbetas1b <- as.data.frame(simbetas1b)
quantile(simbetas1b[,3], c(0.025, 0.975))
mean(simbetas1b[,3]) / sd(simbetas1b[,3])
mean(simbetas1b[,2])


############ Interaction effects #################
formula1c <- bond_zscore ~ lag(bond_zscore) + 
  rc100*c40 + rc100*iclei + c40*iclei + 
  log(population) + log(pcincome) + homeown + unemp + 
  fema_county_fire + fema_county_hur + fema_county_snow

model1c <- simbetas1c <- simbeta1c <- NULL
for (i in 1:m) { # create a matrix that has 10000 simulated coefficients 
  # for all explanatory variables
  model1c <- plm(formula1c, effect="twoway", model="within", data=a.df[[i]])
  simbeta1c <- mvrnorm(10000/m, coefficients(model1c), vcovSCC(model1c, group="statecode"))
  simbetas1c <- rbind(simbetas1c, simbeta1c)
}
simbetas1c <- as.data.frame(simbetas1c)
colnames(simbetas1c)
mean(simbetas1c[,2]); sd(simbetas1c[,2]); mean(simbetas1c[,2])/sd(simbetas1c[,2])  # 100RC 
mean(simbetas1c[,3]); sd(simbetas1c[,3]) # C40
mean(simbetas1c[,4]); sd(simbetas1c[,4]) # ICLEI
mean(simbetas1c[,12]); sd(simbetas1c[,12]) # 100RC*C40
mean(simbetas1c[,13]); sd(simbetas1c[,13]) # 100RC*ICLEI
mean(simbetas1c[,14]); sd(simbetas1c[,14]) # C40*ICLEI

######## three-way interaction effects ###########
formula1d <- bond_zscore ~ lag(bond_zscore) + 
  rc100*c40 + rc100*iclei + c40*iclei +  rc100*c40*iclei +
  log(population) + log(pcincome) + homeown + unemp + 
  fema_county_fire + fema_county_hur + fema_county_snow

model1d <- simbetas1d <- simbeta1d <- NULL
for (i in 1:m) { # create a matrix that has 10000 simulated coefficients 
  # for all explanatory variables
  model1d <- plm(formula1d, effect="twoway", model="within", data=a.df[[i]])
  simbeta1d <- mvrnorm(10000/m, coefficients(model1d), vcovSCC(model1d, group="statecode"))
  simbetas1d <- rbind(simbetas1d, simbeta1d)
}
simbetas1d <- as.data.frame(simbetas1d)
colnames(simbetas1d)
mean(simbetas1d[,2]); sd(simbetas1d[,2]); mean(simbetas1d[,2])/sd(simbetas1d[,2]) # 100RC
mean(simbetas1d[,3]); sd(simbetas1d[,3]) # C40
mean(simbetas1d[,4]); sd(simbetas1d[,4]) # ICLEI
mean(simbetas1d[,12]); sd(simbetas1d[,12]) # 100RC*C40
mean(simbetas1d[,13]); sd(simbetas1d[,13]) # 100RC*ICLEI
mean(simbetas1d[,14]); sd(simbetas1d[,14]) # C40*ICLEI
mean(simbetas1d[,15]); sd(simbetas1d[,15]) # 100RC*C40*ICLEI

############ without lagged DV #################
formula1e <- bond_zscore ~  rc100 + c40 + iclei +
  log(population) + log(pcincome) + homeown + unemp + fema_county_fire

model1e <- simbetas1e <- simbeta1e <- NULL
for (i in 1:m) { # create a matrix that has 10000 simulated coefficients 
  # for all explanatory variables
  model1e <- plm(formula1e, effect="twoway", model="within", data=a.df[[i]])
  simbeta1e <- mvrnorm(10000/m, coefficients(model1e), vcovSCC(model1e, group="statecode"))
  simbetas1e <- rbind(simbetas1e, simbeta1e)
}
simbetas1e <- as.data.frame(simbetas1e)
colnames(simbetas1e)
mean(simbetas1e[,1]); sd(simbetas1e[,1]); mean(simbetas1e[,1])/sd(simbetas1e[,1])
mean(simbetas1e[,2]); sd(simbetas1e[,2]); mean(simbetas1e[,2])/sd(simbetas1e[,2])
mean(simbetas1e[,3]); sd(simbetas1e[,3]); mean(simbetas1e[,3])/sd(simbetas1e[,3])


####################### Figure 1 ###########################

simmat_fixed <- simbetas1b
coefs <- c(mean(simmat_fixed[,2]), 
           mean(simmat_fixed[,3]),
           mean(simmat_fixed[,4]), 
           mean(simmat_fixed[,1]),
           mean(simmat_fixed[,5]), 
           mean(simmat_fixed[,6]), 
           mean(simmat_fixed[,7]),
           mean(simmat_fixed[,8]), 
           mean(simmat_fixed[,9]),
           mean(simmat_fixed[,10]),
           mean(simmat_fixed[,11]),
           mean(simmat_fixed[,12])
           
)

lower <- c(quantile(simmat_fixed[,2], 0.025), 
           quantile(simmat_fixed[,3], 0.025),
           quantile(simmat_fixed[,4], 0.025),
           quantile(simmat_fixed[,1], 0.025),
           quantile(simmat_fixed[,5], 0.025),
           quantile(simmat_fixed[,6], 0.025),
           quantile(simmat_fixed[,7], 0.025),
           quantile(simmat_fixed[,8], 0.025),
           quantile(simmat_fixed[,9], 0.025),
           quantile(simmat_fixed[,10], 0.025),
           quantile(simmat_fixed[,11], 0.025),
           quantile(simmat_fixed[,12], 0.025)
)

upper <- c(quantile(simmat_fixed[,2], 0.975), 
           quantile(simmat_fixed[,3], 0.975),
           quantile(simmat_fixed[,4], 0.975),
           quantile(simmat_fixed[,1], 0.975),
           quantile(simmat_fixed[,5], 0.975),
           quantile(simmat_fixed[,6], 0.975),
           quantile(simmat_fixed[,7], 0.975),
           quantile(simmat_fixed[,8], 0.975),
           quantile(simmat_fixed[,9], 0.975),
           quantile(simmat_fixed[,10], 0.975),
           quantile(simmat_fixed[,11], 0.975),
           quantile(simmat_fixed[,12], 0.975)
)

se <- c(sd(simmat_fixed[,2]), 
        sd(simmat_fixed[,3]),
        sd(simmat_fixed[,4]),
        sd(simmat_fixed[,1]),
        sd(simmat_fixed[,5]),
        sd(simmat_fixed[,6]),
        sd(simmat_fixed[,7]),
        sd(simmat_fixed[,8]),
        sd(simmat_fixed[,9]),
        sd(simmat_fixed[,10]),
        sd(simmat_fixed[,11]),
        sd(simmat_fixed[,12])
)

coefs <- round(coefs, 3)
lwr <- round(lower, 3)
upr <- round(upper, 3)
se <- round(se, 3)

result <- data.frame(pe=coefs,
                     lower=lwr,
                     upper=upr,
                     se=se)
rownames(result) <- colnames(simbetas1b)[c(2,3,4,1,5,6,7,8,9,10,11,12)]
result

labels <- c("100RC membership \n [.015, .278]", 
            "C40 membership \n [-.003, .228]", 
            "ICELI membership \n [-.030, .192]",
            "Bond score, lagged by one year \n [.241, .530]",
            "Population, logged \n [-.337, .115]", 
            "Per capita income, logged \n [-.504, .586]", 
            "Homeownership (%) \n [-.004, .013]",
            "Unemployment (%) \n [.000, .000]", 
            "County-level fire emergencies, count \n [-.061, .027]",
            "County-level hurricane emergencies, count \n [-.058, .070]",
            "County-level snowstorm emergencies, count \n [-.159, .090]",
            "County-level flood emergencies, count \n [-.110, .174]"
)

col2 <- RColorBrewer::brewer.pal(9,"Paired")
trace1 <- ropeladder(x=coefs,
                     lower=lower,
                     upper=upper,
                     labels=labels, fontsize=17,
                     plot=1, col=col2, size=0.65, mark="shaded")
tile(trace1,
     limits=c(-0.55,0.55),
     gridlines=list(type="xt", add=T, lty="solid"),
     topaxis=list(add=T, at=seq(-0.4, 0.4, by=0.2),
                  labels=c("-0.4","-0.2","0","0.2","0.4")),
     maintitle=list(labels="Figure 1. Margial effects of each covariate on average bond score (95% CI)"),
     width=list(spacer=0.1, plot=1.2, leftborder=0),
     height=list(plot=5,spacer=6,plottitle=3,xaxistitle=3.5,topaxistitle=3.5),
     output=list(outfile="figure1", width=12, height=18, pointsize=15)
)

######################## Figure 2 ################################

pd <- list(NULL) # plotdata => pd
for (i in 1:m) {
  pd <- rbind(pd, a.df[[i]])
}

allpd <- pd %>% 
  group_by(rc100city, year) %>%   
  summarise(score=mean(bond_zscore), 
            snp=mean(bond1_score), 
            moody=mean(bond2_score),
            fitch=mean(bond3_score))

allplot <- allpd %>% 
  ggplot(aes(year, score, color=as.factor(rc100city), 
             linetype=as.factor(rc100city))) + geom_line(size=1) +
  geom_vline(xintercept=2013) +  
  annotate("text", x=2009.5, y=-0.3, col="black", 
           label="100RC join", fontface=2, size=5) +
  annotate("text", x=2001, y=0.23, col="black", 
           label="Average of three rating scores \n standardized", fontface=2, size=5) +
  xlab("Year") +
  ylab("Standardized bond rating score \n (all three rating companies)") +
  scale_color_manual(name="100RC", values=col[1:2],
                     labels=c("Non-member (N=69)", "Member (N=21)")) +
  scale_linetype_manual(name="100RC", values=c("dotted","solid"),
                        labels=c("Non-member (N=69)", "Member (N=21)")) +
  theme_classic(base_size=15)

snpplot <- allpd %>% 
  ggplot(aes(year, snp, color=as.factor(rc100city),
             linetype=as.factor(rc100city))) + geom_line(size=1) +
  geom_vline(xintercept=2013) +  
  annotate("text", x=2009.5, y=18.5, col="black", 
           label="100RC join", fontface=2, size=5) +
  annotate("text", x=2000, y=20, col="black", 
           label="S&P bond rating score", fontface=2, size=5) +
  xlab("Year") +
  ylab("S&P bond rating score") +
  scale_color_manual(name="100RC", values=col[1:2],
                     labels=c("Non-member (N=69)", "Member (N=21)")) +
  scale_linetype_manual(name="100RC", values=c("dotted","solid"),
                        labels=c("Non-member (N=69)", "Member (N=21)")) +
  theme_classic(base_size=15)

moodyplot <- allpd %>% 
  ggplot(aes(year, moody, color=as.factor(rc100city),
             linetype=as.factor(rc100city))) + geom_line(size=1) +
  geom_vline(xintercept=2013) +  
  annotate("text", x=2009.5, y=17.25, col="black", 
           label="100RC join", fontface=2, size=5) +
  annotate("text", x=2001, y=19, col="black", 
           label="Moody's bond rating score", fontface=2, size=5) +
  xlab("Year") +
  ylab("Moody's bond rating score") +
  scale_color_manual(name="100RC", values=col[1:2],
                     labels=c("Non-member (N=69)", "Member (N=21)")) +
  scale_linetype_manual(name="100RC", values=c("dotted","solid"),
                        labels=c("Non-member (N=69)", "Member (N=21)")) +
  theme_classic(base_size=15)

fitchplot <- allpd %>% 
  ggplot(aes(year, fitch, color=as.factor(rc100city),
             linetype=as.factor(rc100city))) + geom_line(size=1) +
  geom_vline(xintercept=2013) +  
  annotate("text", x=2009.8, y=15.8, col="black", 
           label="100RC join", fontface=2, size=5) +
  annotate("text", x=2000, y=17.6, col="black", 
           label="Fitch bond rating score", fontface=2, size=5) +
  xlab("Year") +
  ylab("Fitch bond rating score") +
  scale_color_manual(name="100RC", values=col[1:2],
                     labels=c("Non-member (N=69)", "Member (N=21)")) +
  scale_linetype_manual(name="100RC", values=c("dotted","solid"),
                        labels=c("Non-member (N=69)", "Member (N=21)")) +
  theme_classic(base_size=15)

figure2 <- ggarrange(allplot, snpplot, moodyplot, fitchplot,
                     common.legend=T, legend="bottom") 

figure2 <- annotate_figure(figure2,
                           top=text_grob("Figure 2. Small multiples of parallel trends in bond ratings between 100RC members group non-members group", face="bold", size=17)) +
  theme(plot.margin=unit(c(1,1,1,1),"cm"))

width<-16
ggsave("figure2.pdf", width=width, height=width/1.618)

################ DiD estimates ###################

did_formula_all <- bond_zscore ~ post + rc100city + post*rc100city +
  log(population) + log(pcincome) + homeown + unemp + fema_county_fire
did_formula_1 <- bond1_score ~ post + rc100city + post*rc100city +
  log(population) + log(pcincome) + homeown + unemp + fema_county_fire
did_formula_2 <- bond2_score ~ post + rc100city + post*rc100city +
  log(population) + log(pcincome) + homeown + unemp + fema_county_fire
did_formula_3 <- bond3_score ~ post + rc100city + post*rc100city +
  log(population) + log(pcincome) + homeown + unemp + fema_county_fire

did_model_all <- did_model_1 <- did_model_2 <- did_model_3 <- NULL
did_beta_all <- did_beta_1 <- did_beta_2  <- did_beta_3  <- NULL
did_betas_all <- did_betas_1 <- did_betas_2 <- did_betas_3 <- NULL

for (i in 1:m) {
  did_data <- a.df[[i]]
  did_data$post <- ifelse(did_data$year>=2013, 1, 0)
  did_model_all <- lm(did_formula_all,did_data)
  did_model_1 <- lm(did_formula_1,did_data)
  did_model_2 <- lm(did_formula_2,did_data)
  did_model_3 <- lm(did_formula_3,did_data)
  
  did_beta_all <- mvrnorm(10000/m, coefficients(did_model_all), vcovHC(did_model_all, type="HC4", cluster="statecode"))
  did_betas_all <- rbind(did_betas_all, did_beta_all)
  
  did_beta_1 <- mvrnorm(10000/m, coefficients(did_model_1), vcovHC(did_model_all, type="HC4", cluster="statecode"))
  did_betas_1 <- rbind(did_betas_1, did_beta_1)
  
  did_beta_2 <- mvrnorm(10000/m, coefficients(did_model_2), vcovHC(did_model_all, type="HC4", cluster="statecode"))
  did_betas_2 <- rbind(did_betas_2, did_beta_2)
  
  did_beta_3 <- mvrnorm(10000/m, coefficients(did_model_3), vcovHC(did_model_all, type="HC4", cluster="statecode"))
  did_betas_3 <- rbind(did_betas_3, did_beta_3)
  
}
colnames(did_betas_all)


## Table 2 estimates
quantile(did_betas_all[,9], c(0.025, 0.975))
quantile(did_betas_1[,9], c(0.025, 0.975))
quantile(did_betas_2[,9], c(0.025, 0.975))
quantile(did_betas_3[,9], c(0.025, 0.975))


################# interpretation ######################################
sample(1:100) # 36 in the first random draw
datayear <- a.df[[36]] %>% 
  select(city,year,bond_zscore) %>% 
  filter(year==2013) %>% arrange(desc(bond_zscore))

## using the datayear object to interpret as shown in 4. Results section

### missing observations information
missingcount <- data %>% 
  select(city, year, bond1, bond2, bond3) %>% 
  group_by(city) %>% 
  summarize(bond1NA = sum(is.na(bond1)),
            bond2NA = sum(is.na(bond2)),
            bond3NA = sum(is.na(bond3)))

1- (sum(is.na(missingcount$bond1)) / nrow(missingcount))
1- (sum(is.na(missingcount$bond2)) / nrow(missingcount))
1- (sum(is.na(missingcount$bond3)) / nrow(missingcount))
