
#Script to apply designmatch to lactate study
library(MIMICbook)
library(dplyr)
library(DMwR)
library(Hmisc)
library(tableone)
library(designmatch)
library(gurobi)
library(corrplot)

rm(list = ls())

setwd("C://Users//Kristen//Documents//HST")

dat = read.csv("basic_set_LevoTime.csv")

#fill in missing ethnicity data
dat$ethnicity[is.na(dat$ethnicity)] <- "OTHER"

#lump ethnicity data
levels(dat$ethnicity) <- list("ASIAN"=c("ASIAN", "ASIAN - CHINESE", "ASIAN - KOREAN", "ASIAN - ASIAN INDIAN", "ASIAN - FILIPINO"),                              
                              "BLACK/AFRICAN"=c("BLACK/AFRICAN","BLACK/AFRICAN AMERICAN","BLACK/CAPE VERDEAN", "BLACK/HAITIAN"),
                              "HISPANIC/LATINO"=c("HISPANIC OR LATINO", "HISPANIC/LATINO - PUERTO RICAN", "HISPANIC/LATINO - DOMINICAN", "HISPANIC/LATINO - CENTRAL AMERICAN (OTHER)"),                   
                              "OTHER"=c("OTHER","PATIENT DECLINED TO ANSWER","UNABLE TO OBTAIN","UNKNOWN/NOT SPECIFIED","NA","AMERICAN INDIAN/ALASKA NATIVE FEDERALLY RECOGNIZED TRIBE","MULTI RACE ETHNICITY","MIDDLE EASTERN"),
                              "WHITE"=c("PORTUGUESE","WHITE","WHITE - RUSSIAN","WHITE - BRAZILIAN","WHITE - OTHER EUROPEAN","WHITE - EASTERN EUROPEAN"))


#Remove unnecessary columns
data_ana = dat 
data_ana$clear_type <- NULL
data_ana$first_lact_time <- NULL
data_ana$vaso_time <- NULL
data_ana$angus <- NULL
data_ana$LR <- NULL
data_ana$Albumin <- NULL
data_ana$Levophed_time <- NULL 
data_ana$X <- NULL
data_ana$lact_change <- data_ana$last_lact - data_ana$first_lact
data_ana$lact_per_change <- data_ana$lact_change/data_ana$first_lact

data_ana$clear_int <- as.numeric(data_ana$lact_change > 0)

#fill in missing fio2 values
data_ana$fio2[data_ana$fio2 == 0] <- 21
data_ana$fio2[data_ana$fio2 == 20] <- 21
#unsure what a number greater than zero but less than 20 would mean, replace with NA
data_ana$fio2[data_ana$fio2 < 21] <- NA 

#calculate the ratio of spO2 to fiO2
data_ana$stof <- data_ana$spo2_mean/data_ana$fio2

#replace ages greater than 100
data_ana$age_on_admit[data_ana$age_on_admit > 100 ] <- 91.4

#use KNN to impute missing creatinine values
data_imp <- data_ana[,c("apsiii","age_on_admit","renal_failure","liver_disease","congestive_heart_failure","creatinine_mean","heartrate_mean","meanbp_mean","spo2_mean","vent","NS","fio2")]
data_imp$gender[data_ana$gender == "M"] <- 0
data_imp$gender[data_ana$gender == 'F'] <- 1
creatine <- knnImputation(data_imp)
data_ana$creatinine_mean <- creatine$creatinine_mean

#keep only complete cases, i.e. remove cases where heart rate, MAP or spO2 is missing
data_ana <- data_ana[complete.cases(data_ana),]
skew_vars <- c("age_on_admit", "apsiii", "first_lact", "creatinine_mean", "heartrate_mean", "meanbp_mean", "spo2_mean","NS", "stof","fio2")
data_ana <- data_ana[order(data_ana$vaso, decreasing = TRUE),]

data_table <- data_ana
data_table$explicit_sepsis <- as.factor(data_table$explicit_sepsis)
data_table$renal_failure <- as.factor(data_table$renal_failure)
data_table$liver_disease <- as.factor(data_table$liver_disease)
data_table$congestive_heart_failure <- as.factor(data_table$congestive_heart_failure)
data_table$vent <- as.factor(data_table$vent)
summary_tab <- CreateTableOne(data=data_table, strata = c("vaso"))
print(summary_tab, nonnormal = skew_vars)

data_ana[,c("apsiii", "renal_failure", "liver_disease", "congestive_heart_failure", "vent", "vaso", "gender", "ethnicity","admission_type","admission_location")] = 
  lapply(data_ana[,c("apsiii", "renal_failure", "liver_disease", "congestive_heart_failure", "vent", "vaso", "gender","ethnicity","admission_type","admission_location")], as.numeric)

attach(data_ana)

#data visualization to check pre-process
par(mfrow = c(3,3))
hist(apsiii, main = "APS III")
hist(age_on_admit, main = "Age")
hist(NS, main = "Normal Saline")
hist(creatinine_mean, main = "Creatinine")
hist(heartrate_mean, main = "Heart Rate")
hist(meanbp_mean, main = "MAP")
hist(fio2, main = "FiO2")
hist(spo2_mean, main = "SpO2")
hist(stof,main="S to F ratio")
par(mfrow = c(1,1))
hist(first_lact,main = "First Lactate")

#treatment indicator
t_ind = vaso

cor_data <- data_ana[,c("apsiii","gender","ethnicity", "admission_type", "admission_location" ,"age_on_admit"
                        ,"renal_failure","liver_disease","congestive_heart_failure","creatinine_mean","heartrate_mean"
                        ,"meanbp_mean","spo2_mean","vent","NS","fio2","first_lact")]
cor_test_vaso <- cor(cor_data[t_ind == 1,])
corrplot(cor_test_vaso, method="circle")

cor_test_cont <- cor(cor_data[t_ind == 0,])
corrplot(cor_test_cont, method="circle")

corrplot(cor_test_vaso - cor_test_cont, method = "circle")

## Perform matching

#distance matrix 
dist_mat = NULL

#subset matching weight
subset_weight = 1

#Moment balance: constrain the differences in means to be at most 0.05 standard deviations apart
mom_covs = cbind(age_on_admit, ethnicity, gender, first_careunit, admission_type, admission_location, apsiii, renal_failure,
                 liver_disease, congestive_heart_failure, creatinine_mean, 
                 heartrate_mean, meanbp_mean, spo2_mean, vent, 
                 NS, fio2, stof, first_lact, explicit_sepsis)
mom_tols = round(absstddif(mom_covs, t_ind, 0.05), 2)
mom = list(covs = mom_covs, tols = mom_tols)

# Fine balance
fine_covs = cbind(gender, liver_disease, first_careunit, explicit_sepsis)
fine = list(covs = fine_covs)

# Near-exact matching
near_exact_covs = cbind(apsiii, first_lact)
near_exact_devs = cbind(106,11)
near_exact = list(covs = near_exact_covs, devs = near_exact_devs)

# Exact matching
exact_covs = cbind(gender, first_careunit)
exact = list(covs = exact_covs)

#Solver options
t_max = 60*5
solver = "gurobi"
approximate = 1
solver = list(name = solver, t_max = t_max, approximate = approximate, round_cplex = 0, trace = 0)

out = bmatch(t_ind = t_ind, dist_mat = dist_mat, subset_weight = subset_weight, mom = mom, fine = fine, near_exact = near_exact, exact = exact, solver = solver)

t_id = out$t_id
c_id = out$c_id


max(abs(data_ana[c_id,"first_lact"] - data_ana[t_id,"first_lact"]))
max(abs(data_ana[c_id,"apsiii"] - data_ana[t_id,"apsiii"]))

loveplot(mom_covs, t_id, c_id, v_line = 0.1, legend_position = "topright")

meantab(mom_covs, t_ind, t_id, c_id)

for (i in 1:ncol(fine_covs)){
  print(finetab(fine_covs[,i],t_id,c_id))
}

table(exact_covs[t_id,]==exact_covs[c_id,])

all_ind = cbind(t_id,c_id)
match_summary_tab <- CreateTableOne(data=data_table[all_ind,], strata = c("vaso"))
print(match_summary_tab, nonnormal = skew_vars, exact = c("admission_type", "admission_location"))

## Perform analysis

c_id_before = which(t_ind == 0)
t_id_before = which(t_ind == 1)

#data visualization to check pre-process
par(mfrow = c(3,3))
ecdfplot(apsiii, t_id_before, c_id_before, main = "APS III Before Matching")
ecdfplot(age_on_admit, t_id_before, c_id_before, main = "Age Before Matching")
ecdfplot(NS, t_id_before, c_id_before, main = "Normal Saline Before Matching")
ecdfplot(creatinine_mean, t_id_before, c_id_before, main = "Creatinine Before Matching")
ecdfplot(heartrate_mean, t_id_before, c_id_before, main = "Heart Rate Before Matching")
ecdfplot(meanbp_mean, t_id_before, c_id_before, main = "MAP Before Matching")
ecdfplot(fio2, t_id_before, c_id_before, main = "FiO2 Before Matching")
ecdfplot(spo2_mean, t_id_before, c_id_before, main = "SpO2 Before Matching")
ecdfplot(stof,t_id_before, c_id_before, main="S to F ratio Before Matching")

par(mfrow = c(3,3))
ecdfplot(apsiii, t_id, c_id, main = "APS III After Matching")
ecdfplot(age_on_admit, t_id, c_id, main = "Age After Matching")
ecdfplot(NS, t_id, c_id, main = "Normal Saline After Matching")
ecdfplot(creatinine_mean, t_id, c_id, main = "Creatinine After Matching")
ecdfplot(heartrate_mean, t_id, c_id, main = "Heart Rate After Matching")
ecdfplot(meanbp_mean, t_id, c_id, main = "MAP After Matching")
ecdfplot(fio2, t_id, c_id, main = "FiO2 After Matching")
ecdfplot(spo2_mean, t_id, c_id, main = "SpO2 After Matching")
ecdfplot(stof,t_id, c_id, main="S to F ratio After Matching")
par(mfrow = c(1,1))


#look at differences in APS III in matched pairs
aps_diff <- data_ana[t_id,"apsiii"] - data_ana[c_id,"apsiii"]
hist(aps_diff)


plot(data_ana[t_id,"first_lact"],data_ana[c_id,"first_lact"],xlab = "Treated, first lactate", ylab = "Control, first lactate")
abline(a = 0, b = 1)
abline(a = -1, b = 1, lty = 2)
abline(a = 1, b= 1, lty = 2)


plot(data_ana[t_id,"apsiii"],data_ana[c_id,"apsiii"])
abline(a = 0, b = 1)
abline(a = -15, b = 1, lty = 2)
abline(a = 15, b= 1, lty = 2)

match_data = data.frame(Received_vaso = data_ana[t_id,"clear_int"], No_vaso = data_ana[c_id,"clear_int"])

tab.match1 <- table(match_data$Received_vaso,match_data$No_vaso,dnn = c("Vasopressin","Matched Control"))
tab.match1
tab.match1[2,1]/tab.match1[1,2]
paste("95% Confint", round(exp(c(log(tab.match1[2,1]/tab.match1[1,2]) - qnorm(0.975)*sqrt(1/tab.match1[1,2] +1/tab.match1[2,1]),log(tab.match1[2,1]/tab.match1[1,2]) + qnorm(0.975)*sqrt(1/tab.match1[1,2] +1/tab.match1[2,1])) ),2))
mcnemar.test(tab.match1)

lact_vals_data = data.frame(Received_vaso_frst = data_ana[t_id,"first_lact"], Received_vaso_last = data_ana[t_id,"last_lact"],
                            No_vaso_first = data_ana[c_id,"first_lact"], No_vaso_last = data_ana[c_id,"last_lact"])

match_lact_change <- data_ana[t_id,"lact_change"] - data_ana[c_id,"lact_change"]
d <- density(match_lact_change)
plot(d,main = "Estimated Density of Change in Lactate Pair Differences",xlim = c(-15,15))
abline(v = 0:0.1, col = "black", lty=2)

match_los <- data_ana[t_id,"los_hospital"] - data_ana[c_id,"los_hospital"]
plot(density(match_los), main = "Estimated Density of LOS Pair Differences",xlim = c(-65,65))
abline(v = 0:0.1, col = "black", lty=2)

boxplot(match_lact_change, horizontal = TRUE, main = "Boxplot of Pair Differences",ylim = c(-15,15))
abline(v = 0:0.1, col = "black", lty=2)


hist(match_lact_change, breaks = seq(from = -10, to = 20, by = 1),xlim = c(-10,20))
abline(v = 0:0.1, col = "black", lty=2)

plot(data_ana[t_id,"lact_change"],data_ana[c_id,"lact_change"],pch = 16, 
     xlab = "Change in Lactate for Treated Group", ylab = "Change in Lactate for Control Group")
abline(v = 0, col = "black")
abline(h = 0, col = "black")
abline(a = 0, b = 1)
abline(a = 0.5, b = 1, col = "black", lty = 2)
abline(a = -0.5, b = 1, col = "black", lty = 2)
abline(a = 1, b = 1, col = "blue", lty = 2)
abline(a = -1, b = 1, col = "blue", lty = 2)
abline(a = 2, b = 1, col = "red", lty = 2)
abline(a = -2, b = 1, col = "red", lty = 2)
abline(a = 3, b = 1, col = "green", lty = 2)
abline(a = -3, b = 1, col = "green", lty = 2)

mean(match_lact_change)

par(mfrow = c(1,1))
ms <- array(c(data_ana[t_id,"first_lact"],data_ana[c_id,"first_lact"]), c(41,2))
plot(apply(ms,c(1),mean), match_lact_change, xlab = "Average first lactate for mached pair", ylab = expression(paste("Lactate", Delta, "treated - ", Delta, "control")))
abline(h = 0, col = "black", lty=2)


par(mfrow = c(1,2))
plot(data_ana[t_id,"first_lact"], data_ana[c_id,"first_lact"],type = "n",xlim = c(0,13),ylim = c(0,12))
text(data_ana[t_id,"first_lact"], data_ana[c_id,"first_lact"], labels = out$group_id)
abline(a = 0, b = 1)
plot(data_ana[t_id,"last_lact"], data_ana[c_id,"last_lact"],type = "n",xlim = c(0,13),ylim = c(0,12))
text(data_ana[t_id,"last_lact"], data_ana[c_id,"last_lact"], labels = out$group_id)
abline(a = 0, b = 1)
par(mfrow = c(1,1))

match_data2 = data.frame(Received_vaso = data_ana[t_id,"hospital_expire_flag"], No_vaso = data_ana[c_id,"hospital_expire_flag"])

tab.match2 <- table(match_data2$Received_vaso,match_data2$No_vaso,dnn = c("Vasopressin","Matched Control"))
tab.match2
tab.match2[2,1]/tab.match2[1,2]
paste("95% Confint", round(exp(c(log(tab.match2[2,1]/tab.match2[1,2]) - qnorm(0.975)*sqrt(1/tab.match2[1,2] +1/tab.match2[2,1]),log(tab.match2[2,1]/tab.match2[1,2]) + qnorm(0.975)*sqrt(1/tab.match2[1,2] +1/tab.match2[2,1])) ),2))
mcnemar.test(tab.match2)

match_data3 = data.frame(Received_vaso = data_ana[t_id,"los_hospital"], No_vaso = data_ana[c_id,"los_hospital"])

mean(data_ana[t_id,"los_hospital"])
mean(data_ana[c_id,"los_hospital"])


wilcox.test(data_ana[c_id,"lact_change"],data_ana[t_id,"lact_change"], paired = TRUE)
wilcox.test(data_ana[c_id,"los_hospital"],data_ana[t_id,"los_hospital"], paired = TRUE)

hist(data_ana$first_lact[data_ana$hospital_expire_flag == 0])#, breaks = c(0,2,4,6,8,10,12,14,16,18,20,22), col=rgb(0,0,1,0.5), main = "First Lactate Values")
hist(data_ana$first_lact[data_ana$hospital_expire_flag == 1])#, breaks = c(0,2,4,6,8,10,12,14,16,18,20,22), col=rgb(1,0,0,0.5), add = T)
mean(data_ana$first_lact[data_ana$hospital_expire_flag == 0])
mean(data_ana$first_lact[data_ana$hospital_expire_flag == 1])
median(data_ana$first_lact[data_ana$hospital_expire_flag == 0])
median(data_ana$first_lact[data_ana$hospital_expire_flag == 1])
t.test(data_ana$first_lact[data_ana$hospital_expire_flag == 0], data_ana$first_lact[data_ana$hospital_expire_flag == 1])
t.test(log(data_ana$first_lact[data_ana$hospital_expire_flag == 0]), log(data_ana$first_lact[data_ana$hospital_expire_flag == 1]))

match_all <- data_ana[all_ind,]
hist(match_all$first_lact[match_all$hospital_expire_flag == 1])#, breaks = c(0,1,2,3,4,5,6,7,8),col=rgb(1,0,0,0.5))
hist(match_all$first_lact[match_all$hospital_expire_flag == 0])#, breaks = c(0,1,2,3,4,5,6,7,8), add=T,col=rgb(0,0,1,0.5))

par(mfrow = c(1,1))
hist(data_ana[c_id,"lact_change"],breaks = seq(-10,20),col=rgb(1,0,0,0.5),
     xlab = "Change in lactate (mmol/L)", main = "Change in lactate for the control and treated groups")
hist(data_ana[t_id,"lact_change"],breaks = seq(-10,20), add=T,col=rgb(0,0,1,0.5),xlab = "Change in lactate (mmol/L)")

hist(data_ana[t_id,"lact_change"],breaks = seq(-10,20), col=rgb(0,0,1,0.5))
hist(data_ana[c_id,"lact_change"],breaks = seq(-10,20), add=T, col=rgb(1,0,0,0.5))

median(data_ana[c_id,"lact_change"])
median(data_ana[t_id,"lact_change"])

png("ScatterLac_Sensitivity.png",width = 6.5, height = 5, units = "in", res = 300)
pdf <- data.frame(ydata = data_ana[c_id,"lact_change"], xdata = data_ana[t_id,"lact_change"])
pdf$shape[pdf$ydata <= 0 & pdf$xdata <= 0] <- "dd"
pdf$shape[pdf$ydata > 0 & pdf$xdata <= 0] <- "ud"
pdf$shape[pdf$ydata <= 0 & pdf$xdata > 0] <- "du"
pdf$shape[pdf$ydata > 0 & pdf$xdata > 0] <- "uu"
ggplot(pdf, aes(x = xdata, y = ydata, shape = shape)) + geom_point() + 
  geom_abline(intercept = 0) + 
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) + 
  scale_x_continuous(name = "Treated Patient Change in Serum Lactate (mmol/L)", limits = c(-7,7)) +
  scale_y_continuous(name = "Control Patient Change in Serum Lactate (mmol/L)", limits = c(-7,7)) +
  theme_bw() + annotate("text", x = 2.5, y = 5, label = "n = 6 concordant\n positive pairs") +
  annotate("text", x = 4, y = -5, label = "n = 14 discordant\n pairs") +
  annotate("text", x = -4, y = 5, label = "n = 3 discordant\n pairs") +
  annotate("text", x = -2.5, y = -5, label = "n = 18 concordant\n negative pairs")
dev.off()
