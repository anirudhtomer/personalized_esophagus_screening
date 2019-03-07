source("src/simCommon.R")

dataSetNums = 4:20

n_sub = 1500
n_training = 1000
last_seed = 104

getNextSeed = function(last_seed){
  return(last_seed + 1)
}

FIXED = "Fixed"
thresholds = 1 - seq(0, 0.05, by= 0.005)[-1]
for(i in dataSetNums){
  print(paste("******** Started working on Data Set: ", i, "*******"))
  
  last_seed = getNextSeed(last_seed)
  repeat{
    print(paste("Using seed:", last_seed))
    set.seed(last_seed)
    new_data = try(generateBaselineData(n_sub = n_sub),T)
    if(inherits(new_data, "try-error")){
      print(new_data)
      last_seed = getNextSeed(last_seed)
      print("Error generating simulation data: trying again")
    }else{
      joint_model_data = try(fitJointModel(new_data$patient_data, n_training, 
                                           cens_start_time = 0, cens_end_time = 30),T)
      if(inherits(joint_model_data, "try-error")){
        print(joint_model_data)
        last_seed = getNextSeed(last_seed)
        print("Error fitting JM: trying again")
      }else{
        joint_model_data$seed = last_seed
        joint_model_data$b = new_data$b
        rm(new_data)
        break
      }
    }
  }
  
  saveName = paste0("joint_model_data_seed_",joint_model_data$seed,"_simNr_",i, ".Rdata")
  save(joint_model_data, file = paste0("Rdata/simulation/", saveName))
  
  #plus 1 because of FIXED
  schedule_results = do.call(rbind, replicate(length(thresholds)+1, 
                                              joint_model_data$test_data[!duplicated(joint_model_data$test_data$id),], 
                                              simplify = F))
  schedule_results = schedule_results[order(schedule_results$id, decreasing = F),]
  #automatically repeats itself correctly
  schedule_results$methodName = c(FIXED,thresholds)
  schedule_results$nb = schedule_results$offset = NA
  
  # Then we do the FIXED schedule
  print("Running Fixed schedule")
  schedule_results[schedule_results$methodName == FIXED, c("nb", "offset")] = runFixedSchedule(joint_model_data$test_data)
  print("Done running Fixed schedule")
   
  #Then we do the schedule with Dyn. Risk of GR
  for(threshold in thresholds){
     riskMethodName = as.character(threshold)
     
     print(paste("Running",riskMethodName,"schedule"))
     schedule_results[schedule_results$methodName == riskMethodName, c("nb", "offset")] = runRiskBasedSchedule(joint_model_data$fitted_jm, joint_model_data$test_data, threshold)
     print(paste("Done running",riskMethodName,"schedule"))
  }
   
  # print(paste("********* Saving the results ******"))
  # 
  joint_model_data$schedule_results = schedule_results
  saveName = paste0("joint_model_data_seed_",joint_model_data$seed,"_simNr_",i, ".Rdata")
  save(joint_model_data, file = paste0("Rdata/simulation/", saveName))
  rm(schedule_results)
  rm(joint_model_data)
}
