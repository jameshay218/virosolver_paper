expected_files <- paste0("Testing_VIRO_",1:900,".Rout")
real_files <- list.files("~/Downloads/TESTING_CAPACITY_VIRO/")
setdiff(expected_files,real_files)
expected_files[!(expected_files %in% real_files)]
testing_viro_reruns <- which(!(expected_files %in% real_files))

expected_files <- paste0("changing_testing_",1:1200,".Rout")
real_files <- list.files("~/Downloads/CHANGING_TESTING/")
setdiff(expected_files,real_files)
expected_files[!(expected_files %in% real_files)]
changing_testing_reruns <- which(!(expected_files %in% real_files))

expected_files <- paste0("Testing_EN2_",1:1400,".Rout")
real_files <- list.files("~/Downloads/TESTING_CAPACITY_EN2/")
setdiff(expected_files,real_files)
expected_files[!(expected_files %in% real_files)]
testing_en2_reruns <- which(!(expected_files %in% real_files))

expected_files <- paste0("OutSim_",1:450,".Rout")
real_files <- list.files("~/Downloads/FIT_MA_SIM/")
setdiff(expected_files,real_files)
expected_files[!(expected_files %in% real_files)]
fit_ma_sim_reruns <- which(!(expected_files %in% real_files))


## Testing_viro reruns
## 6.symptomatic_testing_virosolver.R
## bash_scripts/testing_simulations_virosolver.sh
dput(testing_viro_reruns)
length(testing_viro_reruns)
# DONE

## Changing testing reruns
## 6.changing_test_rates.R
## bash_scripts/changing_testing_simulations.sh
dput(changing_testing_reruns)
length(changing_testing_reruns)
## DONE

## Testing_EN2 reruns
## 6.symptomatic_testing_epinow2.R
## bash_scripts/testing_simulations_epinow2.sh
dput(testing_en2_reruns)
length(testing_en2_reruns)
## DONE


## FIT_MA_GP reruns
## scripts/cluster_JH/3.full_sim_recovery_cluster.R
dput(fit_ma_sim_reruns)
length(fit_ma_sim_reruns)
## LEAVE IT - NOT WORTH IT