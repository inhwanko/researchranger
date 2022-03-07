# Using simcf and ggplot2 & tile to explore an estimated logistic regression
# Voting example using 2000 NES data after King, Tomz, and Wittenberg
# Inhwan Ko (Based on Chris Adolph, Brian Leung, and Kenya Amano)


# Clear memory
rm(list=ls())

# Load libraries
library(MASS)
library(simcf)
library(tile)
library(RColorBrewer)
library(tidyverse)  # Should come after MASS
source("theme_caviz.R") #ggplot theme 
#library(WhatIf)
# Load data
#file <- "nes00a.csv"
#data <- read.csv(file, header=TRUE)
data <- read.csv("nes00a.csv", header = TRUE)

# Estimate logit model using optim()
# Construct variables and model objects
y <- data$vote00
x <- cbind(data$age,data$age^2,data$hsdeg,data$coldeg)

# Likelihood function for logit
llkLogit <- function(param,y,x) {
  os <- rep(1,length(x[,1]))
  x <- cbind(os,x)
  b <- param[ 1 : ncol(x) ]
  xb <- x%*%b
  sum( y*log(1+exp(-xb)) + (1-y)*log(1+exp(xb)));
               # optim is a minimizer, so min -ln L(param|y)
}

# Fit logit model using optim
lsResult <- lm(y~x)  # use ls estimates as starting values
stval <- lsResult$coefficients  # initial guesses
logitResultOpt <- optim(stval,llkLogit,method="BFGS",hessian=TRUE,y=y,x=x)
                   # call minimizer procedure
peOpt <- logitResultOpt$par   # point estimates
vcOpt <- solve(logitResultOpt$hessian)  # var-cov matrix
seOpt <- sqrt(diag(vcOpt))    # standard errors
llOpt <- -logitResultOpt$value  # likelihood at maximum

# Estimate logit model using glm()
# Set up model formula and model specific data frame
model <- vote00 ~ age + I(age^2) + hsdeg + coldeg
mdata <- extractdata(model, data, na.rm=TRUE)

# Run logit & extract results
logitResult <- glm(model, family=binomial, data=mdata)
peGlm <- logitResult$coefficients  # point estimates
vcGlm <- vcov(logitResult)         # var-cov matrix


## Simulate quantities of interest using simcf
##
## We could do this from the optim or glm results;
## here, we do from glm

# Simulate parameter distributions
sims <- 10000
simbetas <- mvrnorm(sims, peGlm, vcGlm)

# Set up counterfactuals:  all ages, each of three educations
xhyp <- seq(18,97,1)
nscen <- length(xhyp)
nohsScen <- hsScen <- collScen <- cfMake(model, mdata, nscen)
for (i in 1:nscen) {
  # No High school scenarios (loop over each age)
  nohsScen <- cfChange(nohsScen, "age", x = xhyp[i], scen = i)
  nohsScen <- cfChange(nohsScen, "hsdeg", x = 0, scen = i)
  nohsScen <- cfChange(nohsScen, "coldeg", x = 0, scen = i)

  # HS grad scenarios (loop over each age)
  hsScen <- cfChange(hsScen, "age", x = xhyp[i], scen = i)
  hsScen <- cfChange(hsScen, "hsdeg", x = 1, scen = i)
  hsScen <- cfChange(hsScen, "coldeg", x = 0, scen = i)

  # College grad scenarios (loop over each age)
  collScen <- cfChange(collScen, "age", x = xhyp[i], scen = i)
  collScen <- cfChange(collScen, "hsdeg", x = 1, scen = i)
  collScen <- cfChange(collScen, "coldeg", x = 1, scen = i)
}

# Simulate expected probabilities for all scenarios
nohsSims <- logitsimev(nohsScen, simbetas, ci=0.95)
hsSims <- logitsimev(hsScen, simbetas, ci=0.95)
collSims <- logitsimev(collScen, simbetas, ci=0.95)

######################################################################
# Get 3 nice colors for traces
col <- brewer.pal(3,"Dark2")
######################################################################


######################################################################
# ggplot
######################################################################

# Data wrangling in tidyverse
nohsSims_tb <- nohsSims %>%
  bind_rows() %>%
  mutate(
    xhyp = xhyp,  # add "xhyp" as a covariate
    edu = "nohs"  # add a column to identify which edu. scenario
  )

hsSims_tb <- hsSims %>%
  bind_rows() %>%
  mutate(xhyp = xhyp, edu = "hs")

collSims_tb <- collSims %>%
  bind_rows() %>%
  mutate(xhyp = xhyp, edu = "coll")

allSims_tb <- bind_rows(nohsSims_tb, hsSims_tb, collSims_tb)

# Visualize expected probabilities using ggplot2
# Basic: 
allSims_tb %>%
  ggplot(aes(x = xhyp, y = pe, ymax = upper, ymin = lower, colour = edu, fill = edu)) +
  geom_line() +
  geom_ribbon(alpha = 0.2, linetype = 0) +
  labs(y = "Probability of Voting", x = "Age of Respondent")

# More elaborate way if you are interested in approximating Chris's graph...
col <- brewer.pal(3, "Dark2")

allSims_tb %>%
  ggplot(aes(x = xhyp, y = pe, ymax = upper, ymin = lower, colour = edu, fill = edu)) +
  geom_line() +
  geom_ribbon(alpha = 0.2, linetype = 0) +
  annotate(geom = "text", x = 55, y = 0.26, label = "Less than HS", col = col[1]) +
  annotate(geom = "text", x = 49, y = 0.56, label = "High School", col = col[2]) +
  annotate(geom = "text", x = 30, y = 0.87, label = "College", col = col[3]) +
  scale_color_manual(values = rev(col)) +
  scale_fill_manual(values = rev(col)) +
  scale_x_continuous(breaks = seq(20, 90, 10)) +
  scale_y_continuous(breaks = seq(0, 1, 0.2), limits = c(0, 1), expand = c(0, 0)) +
  labs(y = "Probability of Voting", x = "Age of Respondent") +
  
theme_caviz_hgrid +    # Brian's theme
  
theme(
  panel.background = element_rect(fill = NA),
  axis.line.x.bottom = element_line(size = 0.5),
  axis.ticks.length = unit(0.5, "char"),
  legend.position = "none"
)

######################################################################
# tile
######################################################################
# Set up lineplot traces of expected probabilities
#
# When recycling this code, omit the extrapolate input
# if you are unsure how to use it correctly
nohsTrace <- lineplot(x=xhyp,
                      y=nohsSims$pe,
                      lower=nohsSims$lower,
                      upper=nohsSims$upper,
                      col=col[1],
                      #extrapolate=list(data=mdata[,2:ncol(mdata)],
                      #                  cfact=nohsScen$x[,2:ncol(hsScen$x)],
                      #                 omit.extrapolated=TRUE),
                      plot=1)

hsTrace <- lineplot(x=xhyp,
                    y=hsSims$pe,
                    lower=hsSims$lower,
                    upper=hsSims$upper,
                    col=col[2],
                    #extrapolate=list(data=mdata[,2:ncol(mdata)],
                    #                   cfact=hsScen$x[,2:ncol(hsScen$x)],
                    #                   omit.extrapolated=TRUE),
                    plot=1)

collTrace <- lineplot(x=xhyp,
                      y=collSims$pe,
                      lower=collSims$lower,
                      upper=collSims$upper,
                      col=col[3],
                      #extrapolate=list(data=mdata[,2:ncol(mdata)],
                      #                  cfact=collScen$x[,2:ncol(hsScen$x)],
                      #                 omit.extrapolated=TRUE),
                      plot=1)

# Set up traces with labels and legend
labelTrace <- textTile(labels=c("Less than HS", "High School", "College"),
                       x=c( 55,    49,     30),
                       y=c( 0.26,  0.56,   0.87),
                       col=col,
                       plot=1)

legendTrace <- textTile(labels=c("Logit estimates:", "95% confidence", "interval is shaded"),
                        x=c(82, 82, 82),
                        y=c(0.2, 0.15, 0.10),
                        cex=0.9,
                        plot=1)

# Plot traces using tile
tile(nohsTrace,
     hsTrace,
     collTrace,
     labelTrace,
     legendTrace,
     limits=c(18,94,0,1),
     xaxis=list(at=c(20,30,40,50,60,70,80,90)),
     yaxis=list(label.loc=-0.5, major=FALSE),
     xaxistitle=list(labels="Age of Respondent"),
     yaxistitle=list(labels="Probability of Voting"),
     width=list(null=5,yaxistitle=4,yaxis.labelspace=-0.5)#,
     #output=list(file="educationEV",width=5.5)
     )



################################################################
#
# Now consider a new specification adding the variable
# "ever married", or marriedo
#
# We will estimate this new model with glm(), then
# simulate new scenarios for marrieds and non-marrieds
#
# We could also rerun the age x education scenarios if we wanted


# Estimate logit model using glm()
# Set up model formula and model specific data frame
model2 <- vote00 ~ age + I(age^2) + hsdeg + coldeg + marriedo
mdata2 <- extractdata(model2, data, na.rm=TRUE)

# Run logit & extract results
logitM2 <- glm(model2, family=binomial, data=mdata2)
peM2 <- logitM2$coefficients  # point estimates
vcM2 <- vcov(logitM2)         # var-cov matrix


# Simulate parameter distributions
sims <- 10000
simbetasM2 <- mvrnorm(sims, peM2, vcM2)


# Set up counterfactuals:  all ages, each of three educations
xhyp <- seq(18,97,1)
nscen <- length(xhyp)
marriedScen <- notmarrScen <- cfMake(model2, mdata2, nscen)
for (i in 1:nscen) {
  
  # Married (loop over each age)
  # Note below the careful use of before scenarios (xpre) and after scenarios (x)
  #  - we will use the marriedScen counterfactuals in FDs and RRs as well as EVs
  marriedScen <- cfChange(marriedScen, "age", x = xhyp[i], xpre= xhyp[i], scen = i)
  marriedScen <- cfChange(marriedScen, "marriedo", x = 1, xpre= 0, scen = i)

  # Married (loop over each age)
  notmarrScen <- cfChange(notmarrScen, "age", x = xhyp[i], scen = i)
  notmarrScen <- cfChange(notmarrScen, "marriedo", x = 0, scen = i)
}

# Simulate expected probabilities for all scenarios
marriedSims <- logitsimev(marriedScen, simbetasM2, ci=0.95)
notmarrSims <- logitsimev(notmarrScen, simbetasM2, ci=0.95)

# Simulate first difference of voting wrt marriage
marriedFD <- logitsimfd(marriedScen, simbetasM2, ci=0.95)

# Simulate relative risk of voting wrt marriage
marriedRR <- logitsimrr(marriedScen, simbetasM2, ci=0.95)



######################################################################
# ggplot
######################################################################
# Data wrangling in tidyverse
AllMarried <- 
  bind_rows(  
  
  marriedSims %>%
    bind_rows() %>%
    mutate(notmarried = 0,
           xhyp = xhyp,
           type = "EV"),
  
  notmarrSims %>%
    bind_rows() %>%
    mutate(notmarried = 1,
           xhyp = xhyp,
           type = "EV"),
  
  marriedFD %>%
    bind_rows() %>%
    mutate(notmarried = NA,
           xhyp = xhyp,
           type = "FD"),
  
  marriedRR %>%
    bind_rows() %>%
    mutate(notmarried = NA,
           xhyp = xhyp,
           type = "RR"),
  
  
) %>% 
  mutate(notmarried = factor(notmarried))
  
  
  

AllMarried %>% 
  filter(type == "EV") %>% 
ggplot(aes(x = xhyp, y = pe, 
           color = notmarried, fill = notmarried)) +
  # Point estimates lines
  geom_line(show.legend = FALSE) +
  scale_color_manual(values = c(col[1], col[2]), 
                     labels = c("Currently Married", "Not Married")) +
  # CIs for white voters
  geom_ribbon(aes(ymin = lower, ymax = upper),
              alpha = 0.5, linetype = 0, show.legend = FALSE) +
  scale_fill_manual(values = c(col[1], NA)) +
  # CIs for non-white voters
  geom_line(aes(y = upper, linetype = notmarried), show.legend = FALSE) +
  geom_line(aes(y = lower, linetype = notmarried), show.legend = FALSE) +
  scale_linetype_manual(values = c(0, 2)) +
  # Label
  annotate(geom = "text", x = 40, y = 0.86, label = "Currently Married", col = col[1]) +
  annotate(geom = "text", x = 49, y = 0.46, label = "Not Married", col = col[2]) +
  # Other adjustments
  scale_x_continuous(breaks = seq(20, 90, 10)) +
  scale_y_continuous(breaks = seq(0, 1, 0.2), limits = c(0, 1), 
                     expand = c(0, 0)) +
  theme_caviz_hgrid +
  theme(legend.position = c(0.2, 0.13), 
        legend.key.size = unit(0.2, "cm"),
        panel.grid.major.y = element_blank())+
  labs(y = "Probability of Voting", x = "Age of Respondent")



# First Difference
AllMarried %>% 
  filter(type == "FD") %>% 
  ggplot(aes(x = xhyp, y = pe)) +

  # Point estimates lines
  geom_line(color = col[1]) +
  # CIs for white voters
  geom_ribbon(aes(ymin = lower, ymax = upper),
              alpha = 0.5, linetype = 0, 
              show.legend = FALSE, fill = col[1]) +
  # Other adjustments
  scale_x_continuous(breaks = seq(20, 90, 10)) +
  scale_y_continuous(breaks = seq(-0.1, 0.5, 0.1), limits = c(-0.1, 0.5),
                    )+
  
  geom_hline(yintercept = 0)+
                     
  theme_caviz_hgrid +
  theme(legend.position = c(0.2, 0.13), 
        legend.key.size = unit(0.2, "cm"),
        panel.grid.major.y = element_blank())+
  labs(y = "Difference in Probability of Voting", x = "Age of Respondent")
  
  #+geom_text(aes(60, 0.2), label="First difference bewteen married vs unmarried",
  #          size=5)


# Relative Risk
AllMarried %>% 
  filter(type == "RR") %>% 
  ggplot(aes(x = xhyp, y = pe)) +
  
  # Point estimates lines
  geom_line(color = col[1]) +
  # CIs for white voters
  geom_ribbon(aes(ymin = lower, ymax = upper),
              alpha = 0.5, linetype = 0, 
              show.legend = FALSE, fill = col[1]) +
  # Other adjustments
  scale_x_continuous(breaks = seq(20, 90, 10)) +
  scale_y_continuous(breaks = seq(0.9, 1.5, 0.1), limits = c(0.9, 1.5),
  )+
  
  geom_hline(yintercept = 1)+
  
  theme_caviz_hgrid +
  theme(legend.position = c(0.2, 0.13), 
        legend.key.size = unit(0.2, "cm"),
        panel.grid.major.y = element_blank())+
  labs(y = "Relative Risk of Voting", x = "Age of Respondent")


######################################################################
# tile
######################################################################

## Make plots using tile

# Get 3 nice colors for traces
col <- brewer.pal(3,"Dark2")

# Set up lineplot traces of expected probabilities
marriedTrace <- lineplot(x=xhyp,
                         y=marriedSims$pe,
                         lower=marriedSims$lower,
                         upper=marriedSims$upper,
                         col=col[1],
                #         extrapolate=list(data=mdata2[,2:ncol(mdata2)],
                #           cfact=marriedScen$x[,2:ncol(marriedScen$x)],
                #          omit.extrapolated=TRUE),
                         plot=1)

notmarrTrace <- lineplot(x=xhyp,
                         y=notmarrSims$pe,
                         lower=notmarrSims$lower,
                         upper=notmarrSims$upper,
                         col=col[2],
                         ci = list(mark="dashed"),
                #         extrapolate=list(data=mdata2[,2:ncol(mdata2)],
                #           cfact=notmarrScen$x[,2:ncol(notmarrScen$x)],
                #           omit.extrapolated=TRUE),
                         plot=1)


# Set up traces with labels and legend
labelTrace <- textTile(labels=c("Currently Married", "Not Married"),
                       x=c( 35,    53),
                       y=c( 0.8,  0.56),
                       col=col,
                       plot=1)

legendTrace <- textTile(labels=c("Logit estimates:", "95% confidence", "interval is shaded"),
                        x=c(80, 80, 80),
                        y=c(0.2, 0.15, 0.10),
                        cex=0.9,
                        plot=1)

# Plot traces using tile
tile(marriedTrace,
     notmarrTrace,
     labelTrace,
     legendTrace,
     limits=c(18,94,0,1),
     xaxis=list(at=c(20,30,40,50,60,70,80,90)),
     yaxis=list(label.loc=-0.5, major=FALSE),
     xaxistitle=list(labels="Age of Respondent"),
     yaxistitle=list(labels="Probability of Voting"),
     width=list(null=5,yaxistitle=4,yaxis.labelspace=-0.5)#,
     #output=list(file="marriedEV",width=5.5)
     )



# Plot First Difference

# Set up lineplot trace of relative risk
marriedFDTrace <- lineplot(x=xhyp,
                         y=marriedFD$pe,
                         lower=marriedFD$lower,
                         upper=marriedFD$upper,
                         col=col[1],
                  #       extrapolate=list(data=mdata2[,2:ncol(mdata2)],
                  #         cfact=marriedScen$x[,2:ncol(marriedScen$x)],
                  #         omit.extrapolated=TRUE),
                         plot=1)


# Set up baseline: for first difference, this is 0
baseline <- linesTile(x=c(18,94),
                      y=c(0,0),
                      plot=1)

# Set up traces with labels and legend
labelFDTrace <- textTile(labels=c("Married compared \n to Not Married"),
                       x=c( 40),
                       y=c( 0.20),
                       col=col[1],
                       plot=1)

legendFDTrace <- textTile(labels=c("Logit estimates:", "95% confidence", "interval is shaded"),
                        x=c(80, 80, 80),
                        y=c(-0.02, -0.05, -0.08),
                        cex=0.9,
                        plot=1)

# Plot traces using tile
tile(marriedFDTrace,
     labelFDTrace,
     legendFDTrace,
     baseline,
     limits=c(18,94,-0.1,0.5),
     xaxis=list(at=c(20,30,40,50,60,70,80,90)),
     yaxis=list(label.loc=-0.5, major=FALSE),
     xaxistitle=list(labels="Age of Respondent"),
     yaxistitle=list(labels="Difference in Probability of Voting"),
     width=list(null=5,yaxistitle=4,yaxis.labelspace=-0.5)#,
     #output=list(file="marriedFD",width=5.5)
     )


# Plot Relative Risk

# Set up lineplot trace of relative risk
marriedRRTrace <- lineplot(x=xhyp,
                         y=marriedRR$pe,
                         lower=marriedRR$lower,
                         upper=marriedRR$upper,
                         col=col[1],
                 #        extrapolate=list(data=mdata2[,2:ncol(mdata2)],
                 #           cfact=marriedScen$x[,2:ncol(marriedScen$x)],
                 #           omit.extrapolated=TRUE),
                         plot=1)


# Set up baseline: for relative risk, this is 1
baseline <- linesTile(x=c(18,94),
                      y=c(1,1),
                      plot=1)

# Set up traces with labels and legend
labelRRTrace <- textTile(labels=c("Married compared \n to Not Married"),
                       x=c( 55),
                       y=c( 1.25),
                       col=col[1],
                       plot=1)

legendRRTrace <- textTile(labels=c("Logit estimates:", "95% confidence", "interval is shaded"),
                          x=c(80, 80, 80),
                          y=c(0.98, 0.95, 0.92),
                          cex=0.9,
                          plot=1)

# Plot traces using tile
tile(marriedRRTrace,
     labelRRTrace,
     legendRRTrace,
     baseline,
     limits=c(18,94,0.9,1.5),
     xaxis=list(at=c(20,30,40,50,60,70,80,90)),
     yaxis=list(label.loc=-0.5, major=FALSE),
     xaxistitle=list(labels="Age of Respondent"),
     yaxistitle=list(labels="Relative Risk of Voting"),
     width=list(null=5,yaxistitle=4,yaxis.labelspace=-0.5)#,
     #output=list(file="marriedRR",width=5.5)
     )

 
