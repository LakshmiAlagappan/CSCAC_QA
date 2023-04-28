# Set the working directory to be the location where the project folder is saved
wdir <- paste0("C:/Users/alagl/OneDrive - Wilmar International Limited/Desktop/IPP Research/Projects/CSCAC_QA")
setwd(wdir)
# Setting the directory for the plots generated to be saved
plot_counter = 1
plot_results_dir = paste0(wdir,"/Plot_Results")
# Directory where data for the project is stored
rds_data_dir = paste0(wdir,"/Data")
# directory where the model results will be stored
model_results_dir = paste0(wdir, "/Models_Results")

dir.create(plot_results_dir)
dir.create(model_results_dir)

source("C:/Users/alagl/OneDrive - Wilmar International Limited/Desktop/IPP Research/Projects/helpers.R")
#source("./preprocessing.R")

Batch_1 = readRDS(paste0(rds_data_dir,"/Batch_1.rds"))
Batch_2 = readRDS(paste0(rds_data_dir,"/Batch_2.rds"))

Batch_1_pp = readRDS(paste0(rds_data_dir,"/Batch_1_pp.rds"))
Batch_2_pp = readRDS(paste0(rds_data_dir,"/Batch_2_pp.rds"))


ref = "10_100"

Train_all = readRDS(paste0(rds_data_dir,"/Train_all.rds"))
Train = readRDS(paste0(rds_data_dir,"/Train.rds"))
Transfer = readRDS(paste0(rds_data_dir,"/Transfer_5.rds")) #Transfer spec for 5 batches
Test = readRDS(paste0(rds_data_dir,"/Test_5.rds")) #Test spec for 5 batches

# Example for 1 Test Batch #
bcm = build_bc_model(Train, Transfer[[3]], "PDS")
thresh = calc_threshold(unlist(bcm$Transfer_res$res$depth_score),1)
tr_pds = test_bc_model(bcm, Test[[3]], thresh, "PDS")
dt = calc_deviation_tolerance(tr, thresh)
tr$ans_df
tr$dt

bcm = build_bc_model(Train, Transfer[[3]], "CCA")
thresh = calc_threshold(unlist(bcm$Transfer_res$res$depth_score),1)
tr_pcpds = test_bc_model(bcm, Test[[3]], thresh, "CCA")
dt = calc_deviation_tolerance(tr, thresh)
tr$ans_df
tr$dt


#source("./Project_2/binary_classification.R")
#source("./Project_2/comp_classification.R")


