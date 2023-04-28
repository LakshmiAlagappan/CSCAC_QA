####################
match.numeric <- function(x, table) {
are.equal <- function(x, y) isTRUE(all.equal(x, y))
match.one <- function(x, table)
  match(TRUE, vapply(table, are.equal, logical(1L), x = x))
vapply(x, match.one, integer(1L), table)
}
#################### Loading the 3 versions of CSCAC Results ######
model_results_dir = paste0(wdir, "/Models_Results")
k_tune = readRDS(paste0(model_results_dir, "/PDS/find_step_0_3_0p1_PDS.rds"))
ii = as.numeric(match.numeric(as.numeric(k_tune$step), as.numeric(seq(0,3,0.1))))
cscac_pds = k_tune$res_summary_list[[ii]]

k_tune = readRDS(paste0(model_results_dir, "/CCA/find_step_0_3_0p1_CCA.rds"))
ii = as.numeric(match.numeric(as.numeric(k_tune$step), as.numeric(seq(0,3,0.1))))
cscac_cca = k_tune$res_summary_list[[ii]]

k_tune = readRDS(paste0(model_results_dir, "/PCPDS/find_step_0_3_0p1_PCPDS.rds"))
ii = as.numeric(match.numeric(as.numeric(k_tune$step), as.numeric(seq(0,3,0.1))))
cscac_pcpds = k_tune$res_summary_list[[ii]]

n = as.numeric(length(Test))
###### Creation of folder for storing result ######
model_results_dir = paste0(model_results_dir, "/Comparison")
dir.create(model_results_dir)

# Processing the data for classification for comparison
Tr = Train_all
Tr$all_label[Tr$all_label == ref] = "Pass"
Tr$all_label[Tr$all_label != ref & Tr$all_label!="Pass"] = "Fail"

Tr_Pass = Tr[Tr$all_label == "Pass",]
Tr_Fail = Tr[Tr$all_label == "Fail",]
Tr_Pass_rep = Tr_Pass[rep(c(1,2,3,4,5,6), 40),]
Tr = rbind(Tr_Pass_rep, Tr_Fail)

Te = Test
for (i in 1:n){
  x = Te[[i]]$all_label
  x[Te[[i]]$all_label == ref] = "Pass"
  x[Te[[i]]$all_label != ref & Te[[i]]$all_label !="Pass"] = "Fail"
  Te[[i]]$all_label = x
}

Tr_C = Train
Tr_C$all_label[Tr_C$all_label == ref] = "Pass"
Tr_C$all_label[Tr_C$all_label != ref & Tr_C$all_label!="Pass"] = "Fail"

Trs = Transfer
for (i in 1:n){
  x = Trs[[i]]$all_label
  x[Trs[[i]]$all_label == ref] = "Pass"
  x[Trs[[i]]$all_label != ref & Trs[[i]]$all_label !="Pass"] = "Fail"
  Trs[[i]]$all_label = x
}


###################3 Method 1 - PCLDA 4 Comp ######################3
pclda_res = list()
for (i in 1:n){
  pclda_res[[i]] = build_test_pclda(Tr, Te[[i]], "all_label", 4, FALSE)
}
saveRDS(pclda_res, paste0(model_results_dir, "/comp_pclda_res.rds"))

#########################################################

###################3 Method 1 - PCLDA 50 Comp ######################3
pclda_res_50 = list()
for (i in 1:n){
  pclda_res_50[[i]] = build_test_pclda(Tr, Te[[i]], "all_label", 50, FALSE)
}
saveRDS(pclda_res_50, paste0(model_results_dir, "/comp_pclda_res_50.rds"))
#########################################################

###################3 Method 2 - PDS-PCLDA 4 Comp ######################3


pds_pclda_res = list()
for (i in 1:n){
  sim_pds_mod = build_pds_model(Tr_C, Transfer[[i]], 3, 2, "all_label")
  sim_pds_adj = test_pds_model(Te[[i]], sim_pds_mod$model, FALSE)
  pds_pclda_res[[i]] = build_test_pclda(Tr, sim_pds_adj, "all_label", 4,FALSE)
}
saveRDS(pds_pclda_res, paste0(model_results_dir, "/comp_pds_pclda_res.rds"))
#########################################################


###################3 Method 3 - PLS 4 Comp ######################3
pls_all_res = list()
for (i in 1:n){
  pls_all_res[[i]] = build_test_pls_classification(Tr, Te[[i]], 10, "pno_purity", "all_label")
}
saveRDS(pls_all_res, paste0(model_results_dir, "/comp_pls_res.rds"))

#########################################################

pds_pls_res = list()
for (i in 1:n){
  sim_pds_mod = build_pds_model(Tr_C, Transfer[[i]], 3, 2, "all_label")
  sim_pds_adj = test_pds_model(Te[[i]], sim_pds_mod$model, FALSE)
  pds_pls_res[[i]] = build_test_pls_classification(Tr, sim_pds_adj, 10, "pno_purity", "all_label")
}
saveRDS(pds_pls_res, paste0(model_results_dir, "/comp_pds_pls_res.rds"))

#########################################################


######## Compiling the results ############

cscac_pds_sens = cscac_pds$sens_num/cscac_pds$sens_den
cscac_pds_spec = cscac_pds$spec_num/cscac_cca$spec_den

cscac_cca_sens = cscac_cca$sens_num/cscac_cca$sens_den
cscac_cca_spec = cscac_cca$spec_num/cscac_cca$spec_den

cscac_pcpds_sens = cscac_pcpds$sens_num/cscac_pcpds$sens_den
cscac_pcpds_spec = cscac_pcpds$spec_num/cscac_pcpds$spec_den

pclda_4_sens = unlist(lapply(pclda_res, function(x){calc_sens_spec(x$cm, "Pass")}))
pclda_4_spec = unlist(lapply(pclda_res, function(x){calc_sens_spec(x$cm, "Fail")}))

pclda_50_sens = unlist(lapply(pclda_res_50, function(x){calc_sens_spec(x$cm, "Pass")}))
pclda_50_spec = unlist(lapply(pclda_res_50, function(x){calc_sens_spec(x$cm, "Fail")}))

pds_pclda_res_sens = unlist(lapply(pds_pclda_res, function(x){calc_sens_spec(x$cm, "Pass")}))
pds_pclda_res_spec = unlist(lapply(pds_pclda_res, function(x){calc_sens_spec(x$cm, "Fail")}))

pls_res_sens = unlist(lapply(pls_all_res, function(x){calc_sens_spec(x$cm, "Pass")}))
pls_res_spec = unlist(lapply(pls_all_res, function(x){calc_sens_spec(x$cm, "Fail")}))

remove_numeric_0 = function(x){
  if (length(x)==0){
    return (0)
  }
  else{
    return (x)
  }
}
pds_pls_res_sens = unlist(lapply(pds_pls_res, function(x){calc_sens_spec(x$cm, "Pass")}))
pds_pls_res_spec = unlist(lapply(lapply(pds_pls_res, function(x){calc_sens_spec(x$cm, "Fail")}), function(x){remove_numeric_0(x)}))


comp_df = data.frame(batch = 1:n, cscac_pds_sens, cscac_pds_spec,
                     cscac_cca_sens, cscac_cca_spec,
                     cscac_pcpds_sens, cscac_pcpds_spec,
                     pclda_4_sens, pclda_4_spec,
                     pds_pclda_res_sens, pds_pclda_res_spec,
                     pls_res_sens, pls_res_spec,
                     pds_pls_res_sens, pds_pls_res_spec)

saveRDS(comp_df, paste0(model_results_dir, "/comp_all_methods.rds"))
