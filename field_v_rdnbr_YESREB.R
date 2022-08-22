## TRIM DATA TO 95TH PERCENTILE OF RDNBR ##quantile(data$RDNBR, c(0.025,0.975))
##95 PERCENT OF THE DATA LIES BETWEEN  -269 and 1747.44 ##data1 = data[data$RDNBR >= -269.14 & data$RDNBR <= 1474.44,] ##write.csv(data1,"C:/Users/sabaj/OneDrive/Desktop/manuscript_2/data_ms2_trim.csv")
library(gamlss)
respdata= read.csv("C:/Users/sabaj/OneDrive/Desktop/manuscript_2/data_ms2_trim.csv")
respdata= na.omit(respdata)
respdata$prop_deepchar = respdata$meandeepchar/3
respdata$prop_needle = respdata$needle/7
respdata$propCBI = respdata$PlotCBI_mean/3
#----------------------------------------------------------------
CBI= gamlss(propCBI ~ RDNBR*reburn,nu.formula=~RDNBR*reburn, tau.formula=~RDNBR*reburn, data=respdata,family=BEINF)
CCloss= gamlss(prop_lccloss ~ RDNBR*reburn, nu.formula=~RDNBR*reburn, tau.formula=~RDNBR*reburn, data=respdata,family=BEINF)
needle= gamlss(prop_needle ~ RDNBR*reburn, nu.formula=~RDNBR*reburn, tau.formula=~RDNBR*reburn, data=respdata,family=BEINF)
BAkill=  gamlss(prop_Kba ~ RDNBR*reburn,nu.formula=~RDNBR*reburn, tau.formula=~RDNBR*reburn, data=respdata,family=BEINF)
STEMkill= gamlss(prop_Kstem ~ RDNBR*reburn,nu.formula=~RDNBR*reburn, tau.formula=~RDNBR*reburn, data=respdata,family=BEINF)
bolesc= gamlss(prop_meanbolescorch ~ RDNBR*reburn,nu.formula=~RDNBR*reburn, tau.formula=~RDNBR*reburn, data=respdata,family=BEINF)
char.ht = gamlss(prop_meancharht ~ RDNBR*reburn,nu.formula=~RDNBR*reburn, tau.formula=~RDNBR*reburn, data=respdata,family=BEINF)
deepchar =gamlss(prop_deepchar ~ RDNBR*reburn,nu.formula=~RDNBR*reburn, tau.formula=~RDNBR*reburn, data=respdata,family=BEINF)
surfchar = gamlss(prop_surfcharmean ~ RDNBR*reburn,nu.formula=~RDNBR*reburn, tau.formula=~RDNBR*reburn, data=respdata,family=BEINF)
####-------------------------------------------------------------------------------------------------------------
pred <- predictAll(surfchar, data=respdata, type="response")
# create duplicate/dummy dataframe to hold bootstrapping simulations
new_data <- respdata
# create dataset for plotting
plot_data <- data.frame("RDNBR" = rep(seq(from=-98, to=1252, length.out=300),3), # range of covariate you want to plot
"reburn"=c(rep("non reburn",300),rep("non SR",300),rep("SR",300)))

# set up bootstrapping
ndata <- nrow(respdata) # number of observations in dataset
nplot <- nrow(plot_data) # number of predictions for plotting;nplot
boots <- 1000 # number of bootstrap iterations
# matrix to hold results for plotting
yest <- matrix(NA, nrow=nplot, ncol=boots);
# bootstrap predictions
for(i in 1:boots){
  # simulate from the model
  new_data$y <- rBEINF(n=ndata, mu=pred$mu, sigma=pred$sigma, nu=pred$nu, tau=pred$tau); 
  # re-fit the model using simulated data
  ymod <- gamlss(y~RDNBR*reburn, nu.formula= ~RDNBR*reburn, tau.formula= ~RDNBR*reburn, 
                 data=new_data, family=BEINF())
  # Extract new predictions for plot data from re-fitted model
  new_pred <- predictAll(ymod, newdata=plot_data, data=new_data, type="response") 
  new_mu <- new_pred$mu
  new_nu <- new_pred$nu
  new_tau <- new_pred$tau
  new_p0 <- new_nu / (1+new_nu+new_tau) # probability of zeros
  new_p1 <- new_tau / (1+new_nu+new_tau) # probability of ones
  # calculate expected value for new predictions and store in matrix
  yest[,i] <- (1-new_p0)*(new_p1+(1-new_p1)*new_mu)
}
# summarize bootstrap results - extract mean and upper/lower quantiles;
summ_CI <- matrix(NA, nrow=nplot, ncol=3)
for(i in 1:nplot){
  summ_CI[i,1] <- mean(yest[i,], na.rm=TRUE) 
  summ_CI[i,2] <- quantile(yest[i,],probs=0.025, na.rm=TRUE) # lower bound
  summ_CI[i,3] <- quantile(yest[i,],probs=0.975, na.rm=TRUE) # upper bound
}
summ_CI <- data.frame("mean_surf" = summ_CI[,1],
                      "lower_surf" = summ_CI[,2],
                      "upper_surf" = summ_CI[,3])

summ_CI <- cbind(summ_CI,plot_data)

head(summ_CI)

write.csv(summ_CI,
"C:/Users/sabaj/OneDrive/Desktop/manuscript_2/reburns/GEE_composite/SurfCharvRDNBR.csv")



###------------------- how to get AUCS-------------------------------------------------------------------------------
library(pROC)
###dichotomize the response variable at thresholds CANOPY COVER LOSS-------------------------------------------
respdata$dic.C1<-ifelse(respdata$propCBI> 0.05, 1, 0)
respdata$dic.C2<-ifelse(respdata$propCBI > 0.275, 1, 0)
respdata$dic.C3<-ifelse(respdata$propCBI > 0.5, 1, 0) 
respdata$dic.C4<-ifelse(respdata$propCBI> 0.725, 1, 0)
respdata$dic.C5<-ifelse(respdata$propCBI > 0.95, 1, 0)

prob.CBI <- predict(CBI, data=respdata, type="response")

cbi.th05<-roc(respdata$dic.C1, prob.CBI); cbi.th05
cbi.th3<-roc(respdata$dic.C2, prob.CBI); cbi.th3
cbi.th5<- roc(respdata$dic.C3, prob.CBI);cbi.th5
cbi.th7<-roc(respdata$dic.C4, prob.CBI);cbi.th7
cbi.th95<-roc(respdata$dic.C5, prob.CBI);cbi.th95

###dichotomize the response variable at thresholds CANOPY COVER LOSS-------------------------------------------
respdata$dic.C1<-ifelse(respdata$propCBI> 0.05, 1, 0)
respdata$dic.C2<-ifelse(respdata$propCBI > 0.275, 1, 0)
respdata$dic.C3<-ifelse(respdata$propCBI > 0.5, 1, 0) 
respdata$dic.C4<-ifelse(respdata$propCBI> 0.725, 1, 0)
respdata$dic.C5<-ifelse(respdata$propCBI > 0.95, 1, 0)

prob.CC <- predict(CCloss, data=respdata, type="response")

closs.th05<-roc(respdata$dic.C1, prob.CC); closs.th05
closs.th3<-roc(respdata$dic.C2, prob.CC); closs.th3
closs.th5<- roc(respdata$dic.C3, prob.CC);closs.th5
closs.th7<-roc(respdata$dic.C4, prob.CC);closs.th7
closs.th95<-roc(respdata$dic.C5, prob.CC);closs.th95

#needle------------------------------------------------------------------
respdata$dic.N1<-ifelse(respdata$prop_needle > 0.05, 1, 0)
respdata$dic.N2<-ifelse(respdata$prop_needle > 0.275, 1, 0)
respdata$dic.N3<-ifelse(respdata$prop_needle > 0.5, 1, 0) 
respdata$dic.N4<-ifelse(respdata$prop_needle > 0.725, 1, 0)
respdata$dic.N5<-ifelse(respdata$prop_needle > 0.95, 1, 0)

prob.NE <- predict(needle, data=respdata, type="response")

NE.th05<-roc(respdata$dic.N1, prob.NE); NE.th05
NE.th3<-roc(respdata$dic.N2, prob.NE); NE.th3
NE.th5<- roc(respdata$dic.N3, prob.NE); NE.th5
NE.th7<-roc(respdata$dic.N4, prob.NE);NE.th7
NE.th95<-roc(respdata$dic.N5, prob.NE);NE.th95

#Basal area killed -----------------------------------------------------------
respdata$dic.BA1<-ifelse(respdata$prop_Kba > 0.05, 1, 0)
respdata$dic.BA2<-ifelse(respdata$prop_Kba > 0.275, 1, 0)
respdata$dic.BA3<-ifelse(respdata$prop_Kba > 0.5, 1, 0) 
respdata$dic.BA4<-ifelse(respdata$prop_Kba > 0.725, 1, 0)
respdata$dic.BA5<-ifelse(respdata$prop_Kba > 0.95, 1, 0)

prob.BA <- predict(BAkill, data=respdata, type="response")

BA.th05<-roc(respdata$dic.BA1, prob.BA); BA.th05
BA.th3<-roc(respdata$dic.BA2, prob.BA); BA.th3
BA.th5<- roc(respdata$dic.BA3, prob.BA); BA.th5
BA.th7<-roc(respdata$dic.BA4, prob.BA);BA.th7
BA.th95<-roc(respdata$dic.BA5, prob.BA);BA.th95

#stems killed -----------------------------------------------------------
respdata$dic.ST1<-ifelse(respdata$prop_Kstem > 0.05, 1, 0)
respdata$dic.ST2<-ifelse(respdata$prop_Kstem > 0.275, 1, 0)
respdata$dic.ST3<-ifelse(respdata$prop_Kstem > 0.5, 1, 0) 
respdata$dic.ST4<-ifelse(respdata$prop_Kstem > 0.725, 1, 0)
respdata$dic.ST5<-ifelse(respdata$prop_Kstem > 0.95, 1, 0)

prob.ST <- predict(STEMkill, data=respdata, type="response")

ST.th05<-roc(respdata$dic.ST1, prob.ST); ST.th05
ST.th3<-roc(respdata$dic.ST2, prob.ST); ST.th3
ST.th5<- roc(respdata$dic.ST3, prob.ST); ST.th5
ST.th7<-roc(respdata$dic.ST4, prob.ST);ST.th7
ST.th95<-roc(respdata$dic.ST5, prob.ST);ST.th95

#BOLE SCORCH -----------------------------------------------------------
respdata$dic.BO1<-ifelse(respdata$prop_meanbolescorch > 0.05, 1, 0)
respdata$dic.BO2<-ifelse(respdata$prop_meanbolescorch > 0.275, 1, 0)
respdata$dic.BO3<-ifelse(respdata$prop_meanbolescorch > 0.5, 1, 0) 
respdata$dic.BO4<-ifelse(respdata$prop_meanbolescorch > 0.725, 1, 0)
respdata$dic.BO5<-ifelse(respdata$prop_meanbolescorch > 0.95, 1, 0)

prob.BO <- predict(bolesc, data=respdata, type="response")

BO.th05<-roc(respdata$dic.BO1, prob.BO); BO.th05
BO.th3<-roc(respdata$dic.BO2, prob.BO); BO.th3
BO.th5<- roc(respdata$dic.BO3, prob.BO); BO.th5
BO.th7<-roc(respdata$dic.BO4, prob.BO);BO.th7
BO.th95<-roc(respdata$dic.BO5, prob.BO);BO.th95

#CHAR HEIGHT -----------------------------------------------------------
respdata$dic.CH1<-ifelse(respdata$prop_meancharht > 0.05, 1, 0)
respdata$dic.CH2<-ifelse(respdata$prop_meancharht > 0.275, 1, 0)
respdata$dic.CH3<-ifelse(respdata$prop_meancharht > 0.5, 1, 0) 
respdata$dic.CH4<-ifelse(respdata$prop_meancharht > 0.725, 1, 0)
respdata$dic.CH5<-ifelse(respdata$prop_meancharht > 0.95, 1, 0)

prob.CH <- predict(char.ht, data=respdata, type="response")

CH.th05<-roc(respdata$dic.CH1, prob.CH); CH.th05
CH.th3<-roc(respdata$dic.CH2, prob.CH); CH.th3
CH.th5<- roc(respdata$dic.CH3, prob.CH); CH.th5
CH.th7<-roc(respdata$dic.CH4, prob.CH);CH.th7
CH.th95<-roc(respdata$dic.CH5, prob.CH);CH.th95
#DEEP CHAR -----------------------------------------------------------
respdata$dic.DE1<-ifelse(respdata$prop_deepchar > 0.05, 1, 0)
respdata$dic.DE2<-ifelse(respdata$prop_deepchar > 0.275, 1, 0)
respdata$dic.DE3<-ifelse(respdata$prop_deepchar > 0.5, 1, 0) 
respdata$dic.DE4<-ifelse(respdata$prop_deepchar > 0.725, 1, 0)
respdata$dic.DE5<-ifelse(respdata$prop_deepchar > 0.95, 1, 0)

prob.DE <- predict(deepchar, data=respdata, type="response")

DE.th05<-roc(respdata$dic.DE1, prob.DE);DE.th05
DE.th3<-roc(respdata$dic.DE2, prob.DE);DE.th3
DE.th5<-roc(respdata$dic.DE3, prob.DE);DE.th5
DE.th7<-roc(respdata$dic.DE4, prob.DE);DE.th7
DE.th95<-roc(respdata$dic.DE5, prob.DE);DE.th95

#surface CHAR -----------------------------------------------------------
respdata$dic.SU1<-ifelse(respdata$prop_surfcharmean > 0.05, 1, 0)
respdata$dic.SU2<-ifelse(respdata$prop_surfcharmean > 0.275, 1, 0)
respdata$dic.SU3<-ifelse(respdata$prop_surfcharmean > 0.5, 1, 0) 
respdata$dic.SU4<-ifelse(respdata$prop_surfcharmean > 0.725, 1, 0)
respdata$dic.SU5<-ifelse(respdata$prop_surfcharmean > 0.95, 1, 0)

prob.su <- predict(surfchar, data=respdata, type="response")

SU.th05<-roc(respdata$dic.SU1, prob.su);SU.th05
SU.th3<-roc(respdata$dic.SU2, prob.su);SU.th3
SU.th5<-roc(respdata$dic.SU3, prob.su);SU.th5
SU.th7<-roc(respdata$dic.SU4, prob.su);SU.th7
SU.th95<-roc(respdata$dic.SU5, prob.su);SU.th95
