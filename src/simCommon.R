library(JMbayes)
library(splines)
library(survival)
library(MASS)
library(doParallel)

load("Rdata/data.RData")

MAX_FAIL_TIME = 14
max_cores = 4

surv_formula = ~ age2 + gender + BElength_cat + esophagitis
long_formula = ~ time + age2 + gender + BElength_cat + esophagitis

gammas = jmfit_barrett_value_expit$statistics$postMeans$gammas
beta_sox2 = jmfit_barrett_value_expit$statistics$postMeans$betas1
beta_p53 = jmfit_barrett_value_expit$statistics$postMeans$betas2
beta_LGD = jmfit_barrett_value_expit$statistics$postMeans$betas3

generateBaselineData = function(n_sub){
  #Generate gender
  gender = relevel(factor(sample(c("male", "female"), size = n_sub,replace = T, c(0.7337559, 1-0.7337559))), 
                   ref = "male")
  age2 = rnorm(n = n_sub, mean = -1.0338 + 0.43662 * ifelse(gender=="female", 1, 0), sd = 0.936395)
  id = 1:n_sub
  
  BElength_cat_prob = plogis(3.2072 + rnorm(n_sub, mean=0, sd=3.992))
  esophagitis_log_odds_notime = -3.69772 + rnorm(n_sub, mean=0, sd = 1.829)
  
  #every week
  exo_times = seq(0, MAX_FAIL_TIME, by = 7/365)
  patient_data = data.frame(id=rep(id, each=length(exo_times)), 
                            age2 = rep(age2, each=length(exo_times)),
                            gender = rep(gender, each=length(exo_times)),
                            time = rep(exo_times, n_sub))
  patient_data$BElength_cat = sapply(rep(BElength_cat_prob, each=length(exo_times)), rbinom, n=1, size=1)
  patient_data$esophagitis = sapply(plogis(rep(esophagitis_log_odds_notime, each=length(exo_times)) - 0.09237 * patient_data$time), 
                                    rbinom, n=1, size=1)
  ###########################
  # Generate the longitudinal dataset
  ###########################
  D = jmfit_barrett_value_expit$statistics$postMeans$D
  b <- mvrnorm(n_sub, mu = rep(0, nrow(D)), D)
  
  gammas = jmfit_barrett_value_expit$statistics$postMeans$gammas
  beta_sox2 = jmfit_barrett_value_expit$statistics$postMeans$betas1
  beta_p53 = jmfit_barrett_value_expit$statistics$postMeans$betas2
  beta_LGD = jmfit_barrett_value_expit$statistics$postMeans$betas3
  
  hazardFunc = function (time, patient_id) {
    b_subject = b[patient_id,]
    patient_data = patient_data[patient_data$id %in% patient_id,]
    nearest_time_indices = sapply(time, FUN = function(t){which.min(abs(patient_data$time - t))})
    patient_data = patient_data[nearest_time_indices,]
    patient_data$time = time
    
    wGamma = model.matrix(surv_formula, patient_data)[,-1] %*% gammas
    baselinehazard = exp(splineDesign(jmfit_barrett_value_expit$control$knots, time, 
                                      ord = jmfit_barrett_value_expit$control$ordSpline, outer.ok = T) %*% jmfit_barrett_value_expit$statistics$postMeans$Bs_gammas)
    
    long_X = model.matrix(long_formula, patient_data)
    true_sox2_prob = plogis(long_X %*% beta_sox2 + b_subject[1])
    true_p53_prob = plogis(long_X %*% beta_p53 + b_subject[2])
    true_LGD_prob = plogis(long_X %*% beta_LGD + b_subject[3])
    
    y_Alpha = cbind(true_sox2_prob, true_p53_prob, true_LGD_prob) %*% jmfit_barrett_value_expit$statistics$postMeans$alphas
    
    baselinehazard * exp(wGamma + y_Alpha)
  }
  
  invSurvival <- function (t, u, patient_id) {
    log(u) + integrate(hazardFunc, lower = 0.001, upper = t, patient_id)$value
  }
  
  pSurvTime = function(survProb, patient_id){
    Low = 1e-05
    Up <- MAX_FAIL_TIME
    Root <- try(uniroot(invSurvival, interval = c(Low, Up), 
                        u = survProb, patient_id=patient_id, tol = 1e-3)$root, TRUE)
    if(inherits(Root, "try-error")){
      return(NA)
    }else{
      return(Root)
    }
  }
  
  survProbs = runif(n=n_sub, 0, 1)
  patient_data$prog_time = rep(sapply(1:n_sub, function(i){pSurvTime(survProbs[i], i)}),
                               each = length(exo_times))
  
  
  longX_matrix = model.matrix(long_formula, patient_data)
  
  #Now we add longitudinal data to this patient
  patient_data$sox2 = sapply(plogis(longX_matrix %*% beta_sox2 + rep(b[,1], each=length(exo_times))), 
                             rbinom, n=1, size=1)
  patient_data$p53 = sapply(plogis(longX_matrix %*% beta_p53 + rep(b[,2], each=length(exo_times))), 
                            rbinom, n=1, size=1)
  patient_data$LGD = sapply(plogis(longX_matrix %*% beta_LGD + rep(b[,3], each=length(exo_times))), 
                            rbinom, n=1, size=1)
  
  retList = list(patient_data = patient_data, b = b)
  return(retList)
}

fitJointModel = function(patient_data, n_training, 
                         cens_start_time = 0, cens_end_time = 30){
  training_data = patient_data[patient_data$id %in% 1:n_training,]
  test_data = patient_data[patient_data$id > n_training,]
  
  ct = makeCluster(max_cores)
  registerDoParallel(ct)
  training_data = foreach(i=1:n_training, .combine='rbind', 
                          .export=c("simulateProtocol", "MAX_FAIL_TIME"),
                          .packages = c("splines", "JMbayes")) %do%{
                            simulateProtocol(training_data[training_data$id == i,])
                          }
  stopCluster(ct)
  
  
  
  fitted_jm = NULL
  
  return(list(training_data = training_data, test_data=test_data,
              fitted_jm = fitted_jm))
}

#Send 1 row of data
simulateProtocol = function(patient_data){
  available_time_list = patient_data$time
  
  prog_time = patient_data$prog_time[1]
  
  #Choose 6 months of data
  ret_data = patient_data[which.min(abs(available_time_list-0)),]
  
  latest_gap = 3
  next_time = 0.5
  repeat{
    cur_time = next_time
    new_data = patient_data[which.min(abs(available_time_list-cur_time)),]
    new_data$time = cur_time
    ret_data = rbind(ret_data, new_data)
    
    if(new_data$LGD == 1){
      if(latest_gap==3){
        next_time = cur_time + 0.5
      }else{
        next_time = cur_time + 1
      }
    }else{
      if(latest_gap==0.5){
        next_time = cur_time + 1
      }else{
        next_time = cur_time + 3
      }
    }
    
    latest_gap = next_time - cur_time
    
    if(cur_time > min(prog_time, MAX_FAIL_TIME, na.rm = T)){
      if(!is.na(prog_time)){
        ret_data$LGD[nrow(ret_data)] = NA
      }else{
        ret_data = ret_data[-nrow(ret_data),]
      }
      break
    }
  }
  
  return(ret_data)
}

simulateRiskBased = function(patient_data, fitted_jm, threshold){
  pDynSurv = function(survProb, patientDs){
    invDynSurvLastTime <- function (time, u, patientDs, maxRiskDt) {
      u - round(survfitJM(fitted_jm, patientDs, survTimes = time)$summaries[[1]][1, "Mean"])
    }
    
    Low = max(patientDs$time) + 1e-05
    Up <- MAX_FAIL_TIME
    
    tries = tries + 1
    Root <- try(uniroot(invDynSurvLastTime, interval = c(Low, Up), 
                        u = survProb, patientDs = patientDs)$root, TRUE)
    
    if(inherits(Root, "try-error")){
      return(MAX_FAIL_TIME)
    }else{
      return(Root)
    }
  }
  
  ############## main algo starts here ############
  available_time_list = patient_data$time
  prog_time = patient_data$prog_time[1]
  
  #Choose 6 months of data
  ret_data = patient_data[sapply(c(0,0.5), FUN = function(x){
    which.min(abs(available_time_list - x))
  })]
  
  ret_data$time = c(0,0.5)
  
  repeat{
    new_time = pDynSurv(survProb = threshold, patientDs = ret_data)
    new_data = patient_data[which.min(abs(available_time_list-new_time)),]
    new_data$time = new_time
    ret_data = rbind(ret_data, new_data)
    
    if(tail(ret_data$time,1) > min(prog_time, MAX_FAIL_TIME, na.rm = T)){
      if(!is.na(prog_time)){
        ret_data$LGD[nrow(ret_data)] = NA
      }else{
        ret_data = ret_data[-nrow(ret_data),]
      }
      break
    }
  }
  
  return(ret_data)
}