Batch_1 = readRDS(paste0(rds_data_dir,"/Batch_1.rds"))
Batch_2 = readRDS(paste0(rds_data_dir,"/Batch_2.rds"))

Batch_1_pp = snv(correct_bsc(trim(Batch_1, 1, 1751)))
Batch_2_pp = snv(correct_bsc(trim(Batch_2, 1, 1751)))
saveRDS(Batch_1_pp, paste0(rds_data_dir, "/Batch_1_pp.rds"))
saveRDS(Batch_2_pp, paste0(rds_data_dir, "/Batch_2_pp.rds"))
Batches = rbind(Batch_1_pp, Batch_2_pp)

####################
ref = "10_100"
Train_all = Batch_1_pp[Batch_1_pp$spec_batch == "20200831",]
Train = Train_all[Train_all$all_label == ref,]
Train = Train[1:4,]

batches = sort(unique(append(Batch_1_pp$spec_batch, Batch_2_pp$spec_batch)))
batches = batches[-1]

each_batch = lapply(batches, function(x){Batches[Batches$spec_batch == x,]})
transfer = lapply(each_batch, function(x){x[x$all_label == ref, ]})
test_adult = lapply(each_batch, function(x){x[x$all_label != ref, ]})

transfer_4 = lapply(1:length(transfer), function(i){
  transfer[[i]][ks_sample(transfer[[i]], "all_label", 4),]})

test_pure = lapply(1:length(transfer), function(i){
  transfer[[i]][!ks_sample(transfer[[i]], "all_label", 4),]})
test_all = lapply(1:length(test_pure), function(i){
  rbind(test_pure[[i]], test_adult[[i]])})

saveRDS(Train_all, paste0(rds_data_dir, "/Train_all.RDS"))
saveRDS(Train, paste0(rds_data_dir, "/Train.RDS"))
saveRDS(transfer_4, paste0(rds_data_dir, "/Transfer.RDS"))
saveRDS(test_all, paste0(rds_data_dir, "/Test.RDS"))



