####################
####################

model_results_dir = paste0(wdir, "/Models_Results")
type = "CCA"
model_results_dir = paste0(model_results_dir, "/", type)
dir.create(model_results_dir)

# Example for 1 Test Batch #
bcm = build_bc_model(Train, Transfer[[3]], type)
thresh = calc_threshold(unlist(bcm$Transfer_res$res$depth_score),1)
tr_pds = test_bc_model(bcm, Test[[3]], thresh, type)
dt = calc_deviation_tolerance(tr, thresh)
tr$ans_df
tr$dt

n = as.numeric(length(Test))
bc_model_list = lapply(1:n, function(i){build_bc_model(Train, Transfer[[i]], type)})
ans_list = list()
for (i in 1:n){
  bcm = bc_model_list[[i]]
  print(unlist(bcm$Transfer_res$res$depth_score))
  thresh = calc_threshold(unlist(bcm$Transfer_res$res$depth_score),1)
  a = test_bc_model(bcm, Test[[i]], 1, type)$res$res
  print(thresh)
  a = a[order(as.numeric(a$depth_score)),][1:10,]
  print(a)
  dp = thresh #readline()
  tr = test_bc_model(bcm, Test[[i]], as.numeric(dp), type)
  dt = calc_deviation_tolerance(tr, dp)
  print(dt)
  ans = tr$ans_df
  ans$thresh = dp
  ans$dev = dt
  ans_list[[i]] = ans
}


##############################################
# Sensitivity and Specificity Analysis
fname = paste0(model_results_dir, "/basic_bc_model_list_", type, ".rds")
if (file.exists(fname)){
  bc_model_list = readRDS(fname)
}else{
  bc_model_list = lapply(1:length(Transfer), function(i){build_bc_model(Train, Transfer[[i]], type)})
  saveRDS(bc_model_list, fname)
}

sr = seq(0,3,0.1)
ans = find_ideal_step(bc_model_list, Test, sr, type)
fname = paste0(model_results_dir, "/find_step_0_3_0p1_", type, ".rds")
saveRDS(ans, fname)
ans = readRDS(fname)

step = ans$step
fname = paste0(model_results_dir, "/threshold_step_", step, ".rds")
fname_1 = paste0(model_results_dir, "/Test_res_list_step_", step, ".rds")
threshold_list = readRDS(fname)
test_res = readRDS(fname_1)

# Calculation of Deviation tolerance
fname = paste0(model_results_dir, "/deviation_tolerance", "_", type, ".rds")
deviation_tolerance_list = lapply(1:length(test_res), function(i){
  calc_deviation_tolerance(test_res[[i]], threshold_list[[i]])})
saveRDS(deviation_tolerance_list, fname)
