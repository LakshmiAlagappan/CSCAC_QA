###### Preprocessing Functions ######

split_transfer_test = function(df, lst="class", nov="14", n=4){
  #Modifying the label to say "Novel
  df[[lst]][df[[lst]]==nov] = "Novel"
  #Splitting df into those which are novel and those that are not
  nov = df[ df[[lst]]=="Novel",]
  rest = df[ df[[lst]]!="Novel",]
  uni_list = unique(rest[[lst]]) #unique class labels
  
  transfer_all = list()
  test_all = list()
  test_all[["Novel"]] = nov #putting novel class totally into test
  #stratified random sampling
  for (i in uni_list){
    each_df = rest[rest[[lst]] == i,]
    cc = sample(nrow(each_df),n,FALSE)
    transfer_all[[i]] = each_df[cc,]
    test_all[[i]] = each_df[-(cc),]
  }
  transfer = do.call(rbind, transfer_all)
  test = do.call(rbind, test_all)
  return(list(transfer, test))
}

trim = function(df, start, end){
  df_trim = df
  #triming the df in the columns
  df_trim$NIR = df$NIR[,start:end]
  return(df_trim)
}

correct_bsc = function(df, plot_flag = FALSE, ...){
  #using the package baseline
  df_bsc = df
  l = length(list(...))
  if (l!=0){
    bl = baseline::baseline(df$NIR, ...)
  }
  else{
    #default
    bl = baseline::baseline(df$NIR, method = 'als',lambda = 6, p = 0.01, maxit = 100)
    #bl = baseline::baseline(df$NIR, method = "modpolyfit", degree = 4, tol = 0.001, rep = 1000)
  }
  if (plot_flag){
    plot_save(baseline::plot(bl), "bsc") #plot will be saved in the plot_results_dir
  }
  df_bsc$NIR = baseline::getCorrected(bl)
  colnames(df_bsc$NIR) = colnames(df$NIR)
  rownames(df_bsc$NIR) = rownames(df$NIR)
  #full.bl = baseline(df$NIR, method, ..., lambda = 6, p = 0.01, maxit = 100, "als")
  #full.bl = baseline(as.matrix(dlist[[2]]), degree = 4, tol = 0.001, rep = 1000, "modpolyfit")
  return (df_bsc)
}

ip_each_spec  =  function(wn_s,abs_s,wn_p, fill = NA) {
  # using package stats
  if (is.null(fill) && (wn_p[1] < wn_s[1] ||
                        wn_p[length(wn_p)] > wn_s[length(wn_s)])) {
    stop("Extrapolation attempted with fill == NULL")
  }
  selector = wn_p <= wn_s[1] & wn_p >= wn_s[length(wn_s)]
  abs_ip = numeric(length(wn_p))
  if (!is.null(fill)){
    abs_ip[!selector] = fill
  }
  if (sum(selector) < 1) {
    NULL
  } else if (sum(selector) <= 25) {
    abs_ip[selector] = stats::spline(x = wn_s,y = abs_s,xout = wn_p[selector],...)$y
  } else {
    abs_ip[selector] = stats::approx(x = wn_s,y = abs_s,xout = wn_p[selector])$y
  }
  return(abs_ip)
}

interpolate  =  function(sec_df, pri_df){
  wn_p = as.numeric(colnames(pri_df$NIR))
  wn_s = as.numeric(colnames(sec_df$NIR))
  abs_set_s = sec_df$NIR
  abs_set_ip = matrix(nrow = dim(abs_set_s)[1], ncol = dim(pri_df$NIR)[2])
  for (i in 1:dim(abs_set_s)[1]){
    abs_set_ip[i,] = ip_each_spec(wn_s, abs_set_s[i,], wn_p)
  }
  colnames(abs_set_ip) = colnames(pri_df$NIR)
  rownames(abs_set_ip) = rownames(abs_set_s)
  sec_adj_df = sec_df
  sec_adj_df$NIR = abs_set_ip
  return(sec_adj_df)
} 

snv = function(df){
  #one by one pre-processing
  df_snv = df
  df_snv$NIR = mdatools::prep.snv(df$NIR)
  colnames(df_snv$NIR) = colnames(df$NIR)
  rownames(df_snv$NIR) = rownames(df$NIR)
  return(df_snv)
} 

KS = function(df, n){
  num_samples = nrow(df)
  pos = logical(num_samples)
  if (num_samples == n){
    npos = seq(1,num_samples)
    pos[npos] = TRUE
    return (pos)
  } else{
    pc12 = (prcomp(df$NIR, center = TRUE))$x[,1:2]
    ken = prospectr::kenStone(as.matrix(pc12),n)
    npos = ken$model
    pos[npos] = TRUE
    return(pos)
  }
}

ks_sample = function(df,lst, n){
  #lst_num is the number based on which msc is grouped
  print (dim(df))
  print(df[[lst]])
  uni_list = unique(df[[lst]])
  repr_pos = unlist(lapply(uni_list, function(i) {KS(df[df[[lst]] == i,], n)}))
  return(repr_pos)
}

###### Index Functions ######

calc_accuracy = function(GT,Pred){
  acc = (sum(GT==Pred)/length(GT))*100
  return (acc)
}

calc_var = function(c,group){
  var = sum(unlist(lapply(seq(1:nrow(group)), function(i){(group[i,1]-c)^2})))
  return(var)
}

calc_wcss = function (pc_i, lst){
  pc_df = data.frame(pc =  pc_i, target = factor(lst))
  point_spl = lapply(unique(pc_df$target), function(x){pc_df[pc_df$target == x,]})
  centroid = unlist(lapply(point_spl, function(x){mean(x[,1])}))
  wcv = sum(unlist(lapply(seq(1:length(centroid)), 
                          function(i) {calc_var(centroid[i], point_spl[[i]])})))
  return(wcv)
}

calc_bcss = function(pc_i, lst){
  pc_df = data.frame(pc =  pc_i, target = factor(lst))
  point_spl = lapply(unique(pc_df$target), function(x){pc_df[pc_df$target == x,]})
  centroid = unlist(lapply(point_spl, function(x){mean(x[,1])}))
  var = sum(dist(centroid)^2)
  return (var)
}

calc_factor_contrib = function(pc_i, contrib_factor, target_factor){
  pc_df = data.frame(pc = pc_i, 
                     contrib = factor(contrib_factor),
                     target = factor(target_factor))
  
  uni_target_list = unique(pc_df$target)
  varr = 0
  for (ti in uni_target_list){
    current_target = pc_df[pc_df$target==ti,]
    contrib_list = unique(current_target$contrib)
    contrib_split = lapply(contrib_list, function(x){current_target[current_target$contrib == x,]})
    med_list = lapply(contrib_split, function(x){median(x[,1])})
    dist_mat = dist(unlist(med_list))^2
    varr = varr + sum(dist_mat)
  }
  return (varr)
}

calc_indexes = function(df,name, contrib, target){
  print("....Getting Indexes....")
  tot = data.frame(matrix(nrow = 4,ncol = 7))
  pc_df = prcomp(df$NIR,TRUE)$x
  for (i in seq(4)){
    wcss_ = calc_wcss(pc_df[,i], target)
    dunn_ = clusterCrit::intCriteria(as.matrix(pc_df[,i]), as.integer(target), "Dunn")$dunn
    bcss_ = calc_bcss(pc_df[,i], target)
    batch = calc_factor_contrib(pc_df[,i], contrib,target)
    d = data.frame("name" = name, "pc" = i, "wcss" = round(wcss_,4), "bcss" = round(bcss_,4),
                   "div" = round(bcss_/wcss_,4), "batch cont" = round(batch, 4), "dunn" = round(dunn_, 4))
    tot[i,] = d
  }
  colnames(tot) = c("Name", "PC", "WCSS", "BCSS", "B/W", "Batch Contrib", "Dunn")
  return (tot)
}

calc_canberra_score <- function(x, y){
  a1 <-  abs(x - y)
  a2 <- abs(x) + abs(y)
  a3 = a1 / a2
  sc = sum(a1 / a2, na.rm = T)
  return(sc)
}

calc_cms = function(median, test_sample, type = "canberra"){
  return(calc_canberra_score(median, test_sample))
  }

calc_depth_score = function(cm_now, pop){
  Nt = length(pop) + 1
  pop_score = median(sqrt(pop)/Nt)
  my_score = sqrt(cm_now)/Nt
  depth_score = my_score/pop_score
  return(depth_score)
}

calc_threshold = function(x_list, step){
  a = median(x_list) + (step*sd(x_list))
  return(a)
}

calc_deviation_tolerance = function(res, threshold_list){
  dev_tol_list = list()
  for (i in 1:length(res)){
    res_h = res[[i]]
    threshold_here = threshold_list[[i]]
    ordered = res_h$res$res
    ordered = ordered[order(as.numeric(ordered$depth_score)),]
    labels = unlist(ordered[ordered$depth_score<=threshold_here,]$pgt)
    if (length(labels) == 0){
      dev_tol_list[[i]] = 0
    }
    else{
      labels_splt = strsplit(labels, "_")
      labels_p = unlist(lapply(labels_splt, function(x){as.numeric(x[2])}))
      dev_tol = 100 - min(labels_p)
      dev_tol_list[[i]] = dev_tol
    }
  }
  return (dev_tol_list)
}

calc_deviation_tolerance = function(res, threshold){
  ordered = res$res$res
  ordered = ordered[order(as.numeric(ordered$depth_score)),]
  labels = unlist(ordered[ordered$depth_score<=threshold,]$pgt)
  if (length(labels) == 0){
    dev_tol = 0
  }
  else{
    labels_splt = strsplit(labels, "_")
    labels_p = unlist(lapply(labels_splt, function(x){as.numeric(x[2])}))
    dev_tol = 100 - min(labels_p)
  }
  return (dev_tol)
}

calc_reg_accuracy = function(GT,Pred, x){
  df_abs = abs(as.numeric(GT)-as.numeric(Pred))
  df_val = df_abs<=x
  ans = round(sum(df_val)/length(df_val),2)
  ans_fmt = paste0(sum(df_val), "/", length(df_val))
  return (list(correct = sum(df_val), acc = ans, str_fmt = ans_fmt))
}

calc_sens_spec = function(cm, attr){ 
  dt = as.data.frame(cm)
  tot = sum(dt[dt$GT == attr,3])
  pred = dt[dt$GT == attr & dt$Pred == attr, 3]
  return (pred/tot)
}

###### Models/Downstream Functions ######

build_test_pclda = function(Tr_df, Te_df, lst, ncomp, plot_flag = TRUE){
  #using package MASS
  Tr_n = length(Tr_df[[lst]])
  Te_n = length(Te_df[[lst]])
  type = c(rep("Train", Tr_n), rep("Test", Te_n))
  
  dat = rbind(Tr_df, Te_df)
  pca.t = prcomp(dat$NIR, center = TRUE)
  pcs = data.frame(pca.t$x[,1:ncomp], class = factor(dat[[lst]]))
  
  tr = pcs[1:Tr_n,]
  te = pcs[(1+Tr_n):length(type),]
  tr.lda = MASS::lda(class ~ . , data = tr)
  te.predict = predict(tr.lda, newdata = te)
  pred = as.character(te.predict$class)
  acc = calc_accuracy(Te_df[[lst]], pred)
  acc_str = sprintf("%d/%d", sum(pred == Te_df[[lst]]), length(pred))
  cm = table("GT"=Te_df[[lst]], "Pred"=te.predict$class)
  
  n_ind = Te_df[[lst]]=="Novel"
  mac = calc_accuracy(Te_df[[lst]][!n_ind], pred[!n_ind])
  nac = calc_accuracy(Te_df[[lst]][n_ind], pred[n_ind])
  
  if (plot_flag == TRUE){
    p1 = plot_pca(pca.t, color_group = type, shape_group = dat[[lst]])
    plot_save(p1, "PCLDA_PCA")
  }
  return (list(predicted = pred, accuracy = acc, formatted = acc_str, cm = cm, mac = mac, nac = nac))
}

build_test_simca = function(Tr_df, Te_df, lst, ncomp){
  uni_list = unique(Tr_df[[lst]])
  l = length(uni_list)
  sim_list = lapply(1:l, function(i)
  {mdatools::simca(dplyr::filter(Tr_df, Tr_df[[lst]] == uni_list[i])$NIR, 
                   paste0(uni_list[i]), ncomp = ncomp)})
  mm = mdatools::simcam(sim_list)
  pred = predict(mm, Te_df$NIR, as.character(Te_df[[lst]]))
  #summary(pred)
  
  pred_all = list()
  or = rownames(as.data.frame(pred$c.pred[1,,]))
  for (i in 1:length(Te_df[[lst]])){
    p = match(1, pred$c.pred[i,,])
    if (is.na(p)){
      pp = "Novel"
    }
    else{
      pp = or[p]
    }
    pred_all = append(pred_all, pp)
  }
  pred_all_list = unlist(pred_all)
  acc = calc_accuracy(pred_all_list, Te_df[[lst]])
  acc_str = sprintf("%d/%d", sum(pred_all_list == Te_df[[lst]]), length(pred_all_list))
  cm = table("GT"=Te_df[[lst]], "Pred"=pred_all_list)
  
  
  n_ind = Te_df[[lst]]=="Novel"
  mac = calc_accuracy(Te_df[[lst]][!n_ind], pred_all_list[!n_ind])
  nac = calc_accuracy(Te_df[[lst]][n_ind], pred_all_list[n_ind])
  
  return (list(predicted = pred_all_list,accuracy = acc, formatted=acc_str, cm=cm, mac = mac, nac = nac))
}

build_test_pls_classification = function(Tr, Te, ncp, qu_label, cl_label){
  pls_mod = build_pls_model(Tr$NIR, Tr[[qu_label]], ncp, FALSE)
  pls_res = test_pls_model(pls_mod, Te$NIR, Te[[qu_label]], FALSE)
  pred_y = pls_res$res$pred_y
  pred_y = ceiling(pred_y / 0.5) * 0.5
  pred_y[pred_y<0] = 0
  pred_y[pred_y>100] = 100
  class_y = pred_y
  class_y[pred_y==100] = "Pass"
  class_y[pred_y!=100] = "Fail"
  pls_res$class_res = class_y
  cm = table("GT"=Te[[cl_label]], "Pred"=class_y)
  return(list(pls_mod = pls_mod, pls_res = pls_res, cm = cm))
}

build_pop_classification_model = function(ref_df,lst){
  uni_list = unique(ref_df[[lst]])
  ref_pool_list = lapply(uni_list, function(x){ref_df[ref_df[[lst]] == x,]})
  med_list = lapply(ref_pool_list, function(x){apply(x$NIR,2,median)})
  pop_list = lapply(1:length(ref_pool_list), 
                    function(rp_i){unlist(
                      lapply(1:nrow(ref_pool_list[[rp_i]]),
                             function(i){calc_cms(med_list[[rp_i]],
                                                   ref_pool_list[[rp_i]]$NIR[i,])}
                      ))})
  
  ret_list = list(ref_df = ref_df, med_list = med_list, pop_list = pop_list,
                  lst = lst, target_order = uni_list)
  return(ret_list)
}

test_pop_classification_model = function(test_df, pop_model, plot_flag = TRUE){
  wn_l = ncol(test_df$NIR)
  num = length(pop_model$target_order)
  cms_mat = matrix(nrow = nrow(test_df), ncol = num)
  colnames(cms_mat) = pop_model$target_order
  
  res = list()
  for (i in 1:nrow(test_df)){
    cms_list = unlist(lapply(1:num, function(x){
      calc_cms(pop_model$med_list[[x]], test_df$NIR[i,])}))
    cms_mat[i,] = cms_list
    min_cms = min(cms_list)
    min_ind = which(cms_list == min_cms)
    
    initial_pred = pop_model$target_order[min_ind]
    depth_sc = calc_depth_score(min_cms, pop_model$pop_list[[min_ind]])
    res[[i]] = list(sno = i, min_cms = min_cms, min_ind = min_ind, 
                    initial_pred = initial_pred, gt = test_df[i,][[pop_model$lst]],depth_sc = depth_sc)
    cat("\r",paste(round(i/nrow(test_df)*100)," %",sep = ""))
  }
  
  res_df = as.data.frame(Reduce(rbind, res))
  
  if (plot_flag){
    type = append(rep("Ref DF", nrow(pop_model$ref_df)), 
                  rep("Test", nrow(test_df)))
    original = rbind(pop_model$ref_df,test_df)
    p1 = plot_pca(prcomp(original$NIR, TRUE), 
                  color_group = type, 
                  shape_group = original[[pop_model$lst]])
    
    plot_save(p1, "pop_model_testing")
    
    p2 = plot_heatmap(as.data.frame(cms_mat))
    plot_save(p2, "heat_map")
  }
  return(list(test_df = test_df, cms_mat = cms_mat, res = res_df))
}

test_misclassification_detection = function(testing_res, pop_model, depth_cutoff, plot_flag = FALSE){
  res_df = testing_res$res
  pos = round(as.numeric(res_df$depth_sc),2) > depth_cutoff
  res_df$new_pred = res_df$initial_pred
  if (sum(pos) != 0){
    res_df[pos,]$new_pred = "Novel"
  }
  pred = unlist(res_df$new_pred)
  gt = unlist(res_df$gt)
  acc = calc_accuracy(gt,pred)
  acc_str = sprintf("%d/%d", sum(pred == gt), length(pred))
  cm = table("GT" = gt, "Pred" = pred)
  testing_res$res = res_df
  
  n_ind = gt =="Novel"
  mac = calc_accuracy(gt[!n_ind], pred[!n_ind])
  nac = calc_accuracy(gt[n_ind], pred[n_ind])
  
  if (plot_flag){
    yy = unlist(res_df$depth_sc)
    p4 = plot_line(data.frame(1:length(yy), sort(yy)), c("Simple Plot", "Test sample #", "Depth Score", "Line"))
    plot_save(p4, "depth_sc_values sorted")
    
    type = append(rep("Primary", nrow(pop_model$ref_df)), 
                  rep("Test", nrow(testing_res$test_df)))
    original = rbind(pop_model$ref_df,testing_res$test_df)
    p1 = plot_pca(prcomp(original$NIR, TRUE), 
                  color_group = type, 
                  shape_group = original[[pop_model$lst]])
    
      
    type = append(rep("Primary", nrow(pop_model$ref_df)), 
                  rep("Test_Adj_ab4_cutoff", nrow(old_adj_test_df)))
    adj = rbind(pop_model$ref_df,old_adj_test_df)
    p2 = plot_pca(prcomp(adj$NIR, TRUE), 
                  color_group = type, 
                  shape_group = adj[[pop_model$lst]])
    
    type = append(rep("Primary", nrow(pop_model$ref_df)), 
                  rep("Test_Adj_aft_cutoff", nrow(testing_res$adj_test_df)))
    adj = rbind(pop_model$ref_df,testing_res$adj_test_df)
    p3 = plot_pca(prcomp(adj$NIR, TRUE), 
                  color_group = type, 
                  shape_group = adj[[pop_model$lst]])
    
    plot_save(cowplot::plot_grid(p1,p2,p3,ncol = 3), "CSPDS_Test", bh=4, bw = 6)
    plot_save(p1, "CSPDS_Test_SB", 8, 4)
    
  }
  testing_res$predicted = pred
  testing_res$accuracy = acc
  testing_res$formatted = acc_str
  testing_res$cm = cm
  testing_res$mac = mac
  testing_res$nac = nac
  return(testing_res)
}


build_bc_model = function(train_now, transfer_now, type="PDS"){
  Tr_model = build_pop_classification_model(train_now, "all_label")
  
  if (type == "PDS"){
    Transfer_model = build_cspds_model(train_now, transfer_now, 3, 2, "all_label", FALSE)
    Transfer_res = test_cspds_model(transfer_now,Transfer_model$model,Tr_model,FALSE)
  }
  else if (type == "CCA"){
    Transfer_model = build_cscca_model(train_now, transfer_now, 0, "all_label", FALSE)
    Transfer_res = test_cscca_model(transfer_now,Transfer_model$model,Tr_model,FALSE)
  }
  else{
    Transfer_model = build_cspcpds_model(train_now, transfer_now, 1,1,8,"all_label", FALSE)
    Transfer_res = test_cspcpds_model(transfer_now, Transfer_model$model, Tr_model, FALSE)
  }
  
  return(list(Tr_model = Tr_model, Transfer_model = Transfer_model, Transfer_res = Transfer_res))
}

test_bc_model = function(bc_model, test_now, threshold, type="PDS"){
  if (type == "PDS"){
    res = test_cspds_model(test_now,bc_model$Transfer_model$model,bc_model$Tr_model,FALSE)
  }
  else if (type == "CCA"){
    res = test_cscca_model(test_now,bc_model$Transfer_model$model,bc_model$Tr_model,FALSE) 
  }
  else if (type == "PCPDS"){
    res = test_cspcpds_model(test_now,bc_model$Transfer_model$model,bc_model$Tr_model,FALSE)
  }
  else{
    print("None")
    res = test_pop_classification_model(test_now, bc_model$Tr_model, plot_flag = FALSE)
  }
  res$res$pgt = res$res$gt
  res$res$gt[res$res$gt != "10_100"] = "Novel"
  
  pds_res = test_misclassification_detection(res, bc_model$Tr_model, threshold, FALSE)
  pds_res_accuracy = pds_res$accuracy
  numerator = as.numeric(strsplit(pds_res$formatted,"/")[[1]][1])
  denominator = as.numeric(strsplit(pds_res$formatted,"/")[[1]][2])
  
  c1 = pds_res$res[pds_res$res$gt == "Novel",]
  c2 = pds_res$res[pds_res$res$gt != "Novel",]
  spec = sum(c1$new_pred=="Novel")/sum(c1$gt == "Novel")
  sens = sum(c2$new_pred=="10_100")/sum(c2$gt == "10_100")
  
  df = data.frame(acuracy = pds_res_accuracy, correct_pred = numerator, total = denominator, 
                  sens = sens, sens_num = sum(c2$new_pred == "10_100"), sens_den = sum(c2$gt == "10_100"),
                  spec = spec, spec_num = sum(c1$new_pred == "Novel"), spec_den = sum(c1$gt == "Novel"))
  return(list(res = pds_res, ans_df = df))
}

find_ideal_step = function(bc_model_list, Test_list, seq_range, type){
  df = data.frame()
  res_summary_list = list()
  for (ii in 1:length(seq_range)){
    fname = paste0(model_results_dir, "/threshold_step_", seq_range[[ii]], ".rds")
    fname_1 = paste0(model_results_dir, "/Test_res_list_step_", seq_range[[ii]], ".rds")
    
    if (file.exists(fname) & file.exists(fname_1)){
      print(paste0(seq_range[[ii]], "_Exists"))
      threshold_list = readRDS(fname)
      Test_res_list = readRDS(fname_1)
    }
    else{
      print(paste0("Running_", seq_range[[ii]]))
      threshold_list = lapply(1:length(bc_model_list), function (i){
        calc_threshold(unlist(bc_model_list[[i]]$Transfer_res$res$depth_score),seq_range[[ii]])})
      Test_res_list = lapply(1:length(bc_model_list), function(i){
        test_bc_model(bc_model_list[[i]], Test_list[[i]], threshold_list[[i]], type)})
      saveRDS(threshold_list, fname)
      saveRDS(Test_res_list, fname_1)
    }
    res_summary = do.call(rbind, lapply(1:length(bc_model_list), function(i){Test_res_list[[i]]$ans_df}))
    print(res_summary)
    res_summary_list[[ii]] = res_summary
    tot_sen = sum(res_summary$sens_num)/sum(res_summary$sens_den)
    tot_spec = sum(res_summary$spec_num)/sum(res_summary$spec_den)
    temp = data.frame(ii, seq_range[[ii]],tot_sen, tot_spec)
    print(temp)
    df = rbind(df, temp)
  }
  
  print("Here are the results: ")
  print(df)
  print("Look through the results and based on your business strategy, Determine a step for the cut-off")
  step = readline(prompt="Enter the step for the cut_off: " )
  step = as.numeric(step)
  return(list(ans_df = df, step = step, res_summary_list = res_summary_list))
}

build_pls_model = function(Tr_NIR,Tr_label,ncp,cv){
  library(pls)
  df = data.frame(datas = I(as.matrix(Tr_NIR)), y = I(as.matrix(Tr_label)))
  if (cv){
    model = pls::plsr(y ~ datas, ncomp = ncp, data = df, method = "oscorespls", center = TRUE, validation = "CV")
    RMSEP_list = as.numeric(pls::RMSEP(model)$val[2,1,])
  }
  else{
    model = pls::plsr(y ~ datas, ncomp = ncp, data = df, method = "oscorespls", center = TRUE)
    RMSEP_list = as.numeric(pls::RMSEP(model)$val[,1,])
  }
  return(model)
}

test_pls_model = function(model, Te_NIR,Te_label, plot_flag=FALSE){
  df_test = data.frame(datas = I(as.matrix(Te_NIR)), y = I(as.matrix(Te_label)))
  RMSEP_list = as.numeric(RMSEP(model, newdata = df_test)$val[,1,])
  print(RMSEP_list)
  n_min_ = which(RMSEP_list == min(RMSEP_list))
  if (n_min_ <= 1){
    n_min = 1
  }
  else{
    n_min = n_min_-1
  }
  
  res_scores_y = list(scores = predict(model, newdata = df_test, comps =  1:n_min, 
                                       type = "scores", center = TRUE),
                      pred_y = predict(model, newdata = df_test, comps =  1:n_min, center = TRUE))
  
  R2 = as.numeric(R2(model, newdata = df_test, comps = 1:n_min)$val[,1,])
  RMSEP_a = as.numeric(RMSEP_list[(n_min+1)])
  
  a_0 = calc_reg_accuracy(as.numeric(Te_label), as.numeric(res_scores_y$pred_y), 0)
  a_0p5 = calc_reg_accuracy(as.numeric(Te_label), as.numeric(res_scores_y$pred_y), 0.5)
  a_1 = calc_reg_accuracy(as.numeric(Te_label), as.numeric(res_scores_y$pred_y), 1)
  a_2.5 = calc_reg_accuracy(as.numeric(Te_label), as.numeric(res_scores_y$pred_y), 2.5)
  a_5 = calc_reg_accuracy(as.numeric(Te_label), as.numeric(res_scores_y$pred_y), 5)
  a_10 = calc_reg_accuracy(as.numeric(Te_label), as.numeric(res_scores_y$pred_y), 10)
  
  reg_acc_df = data.frame(x = c(0,0.5,1,2.5,5,10), 
                          acc = c(a_0$correct, a_0p5$correct, a_1$correct,
                                  a_2.5$correct, a_5$correct, a_10$correct))
  
  if (plot_flag){
    
    tit = paste0("AtN:", n_min, " RMSEP:",
                 round(RMSEP_a,3), " R2:", round(R2,3))
    
    p2 = plot_dot(data.frame(Te_label, res_scores_y[[2]]),
                     c(tit,"TRUE Values", "Predicted Values")) +
      ggplot2::geom_abline(size = 1, colour = "green") +
      ggplot2::geom_abline(intercept = c(1,2,-1,-2,5,-5), slope = 1, size = 1, color = "red", linetype = "dashed")
    
    p3 = plot_line(data.frame(1:length(RMSEP_list), RMSEP_list),
                      c(paste0("Best RMSEP ",  round(RMSEP_a,3) , ""),"# Components", "RMSEP")) +
      ggplot2::geom_vline(xintercept = n_min+1, colour = "red", size = 1.5)
    
    p4_summary = cowplot::plot_grid(p2, p3, ncol = 2, nrow = 1)
    plot_save(p4_summary, "pls_summary")
    colnames(reg_acc_df) = c("Abs Diff", "Num Samples Pred within")
    p6 = plot_line(reg_acc_df,c("Within", "Abs Diff", "Num", " "))
    plot_save(p6, "within")
  }
  
  reg_acc_df_for = data.frame(x = c(0,0.5,1,2.5,5,10), acc = c(a_0$str_fmt, a_0p5$str_fmt, 
                                                               a_1$str_fmt,a_2.5$str_fmt, a_5$str_fmt, a_10$str_fmt))
  
  return(list(res = res_scores_y, rmsep = RMSEP_a, rec_acc = reg_acc_df_for))
}

###### Batch Correction Functions ######
build_pds_tm = function(pri_subset, sec_subset, mean_win_size, ncomp){
  #using package pls
  i = mean_win_size
  k = i - 1
  Fmat = matrix(0,nrow = ncol(pri_subset),ncol = ncol(pri_subset) - (2*i) + 2)
  intercept = c()
  while (i <= (ncol(pri_subset) - k)){
    window_slave = sec_subset[,(i - k):(i + k)]
    fit = pls::plsr(pri_subset[,i] ~ window_slave,ncomp = ncomp, scale = F, method = "oscorespls") 
    coef_reg = as.numeric(coef(fit, ncomp = ncomp, intercept = TRUE))
    intercept = c(intercept,coef_reg[1])
    coefs = coef_reg[2:length(coef_reg)]
    Fmat[(i - k):(i + k),i - k] = t(coefs)
    rm(coef_reg,fit,coefs)
    i = i + 1
  }
  Fmat = as.matrix(data.frame(matrix(0,nrow = ncol(pri_subset),ncol = k), Fmat,
                 matrix(0,nrow = ncol(pri_subset),ncol = k))) #Padding to include 0s
  Bvec = c(rep(0,k),intercept,rep(0,k)) 
  tm = list(Fmat = Fmat , Bvec = Bvec)
  return(tm)
}

build_pds_model = function(pri_calib_df,sec_calib_df,pds_mean_win,ncomp,lst="class", plot_flag = FALSE){
  s_wn_l = ncol(sec_calib_df$NIR)
  transfer_model  =  build_pds_tm(pri_calib_df$NIR, sec_calib_df$NIR, pds_mean_win, ncomp)
  
  adj_sec = sec_calib_df$NIR %*% transfer_model$Fmat
  adj_sec = adj_sec
  adj_sec = sweep(adj_sec, 2, transfer_model$Bvec, "+")
  if (pds_mean_win != 1){
    adj_sec = cbind(sec_calib_df$NIR[,1:(pds_mean_win - 1)],
                    adj_sec[,pds_mean_win:(s_wn_l - pds_mean_win)],
                    sec_calib_df$NIR[,(s_wn_l - pds_mean_win + 1):s_wn_l])
  }
  colnames(adj_sec) = colnames(sec_calib_df$NIR)
  adj_sec_calib_df = sec_calib_df
  adj_sec_calib_df$NIR = adj_sec
  pri_adj_sec_calib_df = rbind(pri_calib_df, adj_sec_calib_df)
  simple_pds_model = list(tm = transfer_model, pds_mean_win = pds_mean_win, 
                          ncomp = ncomp, pri_adj_sec_calib_df = pri_adj_sec_calib_df, lst = lst)
  
  if (plot_flag){
    p1 = plot_line(data.frame(as.numeric(colnames(pri_calib_df$NIR)),
                                      primary_spec = pri_calib_df$NIR[1,],
                                      secondary_spec = sec_calib_df$NIR[1,],
                                      adj_secondary_spec = adj_sec_calib_df$NIR[1,]))
    plot_save(p1, "PDS_spec", 8,4)
    
    type = append(rep("Primary Calib", nrow(pri_calib_df)), 
                  rep("Secondary Calib", nrow(sec_calib_df)))
    original = rbind(pri_calib_df,sec_calib_df)
    p2 = plot_pca(prcomp(original$NIR, TRUE), 
                  color_group = type, shape_group = original[[lst]])
    
    type = append(rep("Primary Calib", nrow(pri_calib_df)), 
                  rep("Secondary_Adj Calib", nrow(adj_sec_calib_df)))
    p3 = plot_pca(prcomp(pri_adj_sec_calib_df$NIR, TRUE), 
                  color_group = type, shape_group = pri_adj_sec_calib_df[[lst]])
    plot_save(cowplot::plot_grid(p2,p3,ncol = 2), "PDS_PCA", 8,4)
    
  }
  return (list(model = simple_pds_model, adj_sec_calib_df = adj_sec_calib_df))
}

test_pds_model = function(test_df, simple_pds_model, plot_flag = FALSE){
  wn_l = ncol(test_df$NIR)
  adj_test = test_df$NIR %*% simple_pds_model$tm$Fmat
  adj_test = sweep(adj_test, 2, simple_pds_model$tm$Bvec, "+")
  
  if (simple_pds_model$pds_mean_win != 1 ){
    adj_test = cbind(test_df$NIR[,1:(simple_pds_model$pds_mean_win - 1),drop = FALSE],
                     adj_test[,simple_pds_model$pds_mean_win:(wn_l - simple_pds_model$pds_mean_win),drop = FALSE],
                     test_df$NIR[,(wn_l - simple_pds_model$pds_mean_win + 1):wn_l,drop = FALSE])
  }
  colnames(adj_test) = colnames(test_df$NIR)
  adj_test_df = test_df
  adj_test_df$NIR = adj_test
  
  if (plot_flag){
    type = append(rep("Primary-AdjSecondary", nrow(simple_pds_model$pri_adj_sec_calib_df)), 
                  rep("Test", nrow(test_df)))
    original = rbind(simple_pds_model$pri_adj_sec_calib_df,test_df)
    p1 = plot_pca(prcomp(original$NIR, TRUE), 
                  color_group = type, 
                  shape_group = original[[simple_pds_model$lst]])
    
    type = append(rep("Primary-AdjSecondary", nrow(simple_pds_model$pri_adj_sec_calib_df)), 
                  rep("Test_Adj", nrow(adj_test_df)))
    adj = rbind(simple_pds_model$pri_adj_sec_calib_df,adj_test_df)
    p2 = plot_pca(prcomp(adj$NIR, TRUE), 
                  color_group = type, 
                  shape_group = adj[[simple_pds_model$lst]])
    plot_save(cowplot::plot_grid(p1,p2,ncol = 2), "PDS_Test", 8,18)
  }
  return(adj_test_df)
}

build_cca_tm = function(pri_subset, sec_subset){
  #using package mixOmics
  transfer_model = mixOmics::rcc(X = pri_subset,Y = sec_subset,method = "shrinkage",ncomp = ncol(pri_subset))
  return(transfer_model)
}

build_cca_model = function(pri_calib_df,sec_calib_df,k,lst="class", plot_flag=FALSE){
  transfer_model = build_cca_tm(pri_calib_df$NIR, sec_calib_df$NIR)
  ncomp = length(transfer_model$cor[transfer_model$cor >= k]); #k = 0, so ncomp = 1751
  
  sec_W = transfer_model$loadings$Y; #1751x1751
  sec_W = sec_W[,1:ncomp]; 
  
  pri_W = transfer_model$loadings$X; #1751x1751
  pri_W = pri_W[,1:ncomp];
  
  pri_L = as.matrix(pri_calib_df$NIR %*% pri_W)
  sec_L = as.matrix(sec_calib_df$NIR %*% sec_W)
  
  f1 = (MASS::ginv(sec_L)) %*% as.matrix(pri_L);
  f2 = (MASS::ginv(pri_L)) %*% as.matrix(pri_calib_df$NIR);
  
  sec_K = sec_calib_df$NIR %*% sec_W;
  adj_sec = sec_K %*% f1 %*% f2;
  adj_sec_calib_df = sec_calib_df
  adj_sec_calib_df$NIR = adj_sec
  
  pri_adj_sec_calib_df = rbind(pri_calib_df, adj_sec_calib_df)
  simple_ctcca_model = list(tm = transfer_model, f1 = f1, f2 = f2,
                            k = k, pri_adj_sec_calib_df = pri_adj_sec_calib_df, lst = lst)
  
  if (plot_flag){
    p1 = plot_line(data.frame(as.numeric(colnames(pri_calib_df$NIR)),
                                      primary_spec = pri_calib_df$NIR[1,],
                                      secondary_spec = sec_calib_df$NIR[1,],
                                      adj_secondary_spec = adj_sec_calib_df$NIR[1,]))
    plot_save(p1, "Simple_Model_CCA_spec", 8,4)
    type = append(rep("Primary Calib", nrow(pri_calib_df)), 
                  rep("Secondary Calib", nrow(sec_calib_df)))
    original = rbind(pri_calib_df,sec_calib_df)
    p2 = plot_pca(prcomp(original$NIR, TRUE), 
                  color_group = type, shape_group = original[[lst]])
    
    type = append(append(rep("Primary Calib", nrow(pri_calib_df)), rep("Secondary Calib", nrow(sec_calib_df))),
                  rep("Secondary_Adj Calib", nrow(adj_sec_calib_df)))
    all = rbind(pri_calib_df, sec_calib_df, adj_sec_calib_df)
    p3 = plot_pca(prcomp(all$NIR, TRUE), 
                  color_group = type, shape_group = all[[lst]])
    plot_save(cowplot::plot_grid(p2,p3,ncol = 2), "Simple_Model_CCA_PCA", 8, 4)
  }
  return(list(model = simple_ctcca_model, adj_sec_calib_df = adj_sec_calib_df))
}

test_cca_model = function(test_df, simple_ctcca_model, plot_flag=FALSE){
  transfer_model = simple_ctcca_model$tm
  ncomp = length(transfer_model$cor[transfer_model$cor >= simple_ctcca_model$k]);
  sec_W = transfer_model$loadings$Y;
  sec_W = sec_W[,1:ncomp];
  
  test_K = test_df$NIR %*% sec_W;
  adj_test = test_K %*% simple_ctcca_model$f1 %*% simple_ctcca_model$f2;
  adj_test_df = test_df
  adj_test_df$NIR = adj_test
  
  if (plot_flag){
    type = append(rep("Primary-AdjSecondary", nrow(simple_ctcca_model$pri_adj_sec_calib_df)), 
                  rep("Test", nrow(test_df)))
    
    original = rbind(simple_ctcca_model$pri_adj_sec_calib_df,test_df)
    p1 = plot_pca(prcomp(original$NIR, TRUE),
                  color_group = type,
                  shape_group = original[[simple_ctcca_model$lst]])
    
    type = append(rep("Primary-AdjSecondary", nrow(simple_ctcca_model$pri_adj_sec_calib_df)), 
                  rep("Test_Adj", nrow(adj_test_df)))
    adj = rbind(simple_ctcca_model$pri_adj_sec_calib_df,adj_test_df)
    p2 = plot_pca(prcomp(adj$NIR, TRUE), 
                  color_group = type, 
                  shape_group = adj[[simple_ctcca_model$lst]])
    plot_save(cowplot::plot_grid(p1,p2,ncol = 2), "Simple_Transfer_CCA_Test", 8,18)
  }
  return(adj_test_df)
}

build_pcpds_model = function(pri_calib_df,sec_calib_df,pds_mean_win,ncomp,npc,lst="class", plot_flag = TRUE){
  all = rbind(pri_calib_df, sec_calib_df)
  pri_pc = pri_calib_df
  sec_pc = sec_calib_df
  pc_obj = prcomp(all$NIR, TRUE)
  num_pcs_here = ncol(pc_obj$x)
  npc_final = min(num_pcs_here, npc)
  pc_spcae = pc_obj$x[,1:npc_final]
  num_pri = length(pri_pc$NIR[,1])
  num_sec = length(sec_pc$NIR[,1])
  pri_pc$NIR = pc_spcae[1:num_pri,]
  sec_pc$NIR = pc_spcae[(num_pri+1):(num_pri + num_sec),]
  colnames(pri_pc$NIR) = 1:npc_final
  colnames(sec_pc$NIR) = 1:npc_final
  
  res = build_pds_model(pri_pc, sec_pc, pds_mean_win, ncomp,"class", plot_flag)
  adj_full = res$adj_sec_calib_df
  adj_full$NIR = adj_full$NIR %*% t(pc_obj$rotation[,1:npc_final])
  adj_full$NIR = scale(adj_full$NIR, center = -pc_obj$center, scale = FALSE)
  
  colnames(adj_full$NIR) = 1:length(adj_full$NIR[1,])
  res$adj_sec_calib_df = adj_full
  
  res$model[["pc_obj"]] = pc_obj
  res$model[["npc"]] = npc_final
  
  return(res)
}


test_pcpds_model = function(test_df, simple_pcpds_model, plot_flag = TRUE){
  test_pc = test_df
  test_pc$NIR = predict(simple_pcpds_model[["pc_obj"]], test_df$NIR)[,1:simple_pcpds_model[["npc"]],drop = FALSE]
  colnames(test_pc$NIR) = 1:simple_pcpds_model[["npc"]]
  res = test_pds_model(test_pc, simple_pcpds_model, plot_flag)
  res_full = res
  res_full$NIR = res$NIR %*% t(simple_pcpds_model[["pc_obj"]]$rotation[,1:simple_pcpds_model[["npc"]]])
  res_full$NIR = scale(res_full$NIR, center = -simple_pcpds_model[["pc_obj"]]$center, scale = FALSE)
  colnames(res_full$NIR) = 1:length(res_full$NIR[1,])
  return(res_full)
}


build_cspds_model = function(pri_calib_df,sec_calib_df,pds_mean_win,ncomp,lst="class", plot_flag = TRUE){
  s_wn_l = ncol(sec_calib_df$NIR)
  target_order = unique(pri_calib_df[[lst]])
  list_cs_ct = lapply(target_order, function(x){
    build_pds_model(pri_calib_df[pri_calib_df[[lst]] == x,],
                    sec_calib_df[sec_calib_df[[lst]] == x,],
                    pds_mean_win, ncomp, plot_flag = FALSE)})
  
  list_model = lapply(list_cs_ct, function(x){x[[1]]})
  adj_sec_calib_df = Reduce(rbind,lapply(list_cs_ct, function(x) {x[[2]]}))
  pri_adj_sec_calib_df = rbind(pri_calib_df, adj_sec_calib_df)
  cspds_model = list(list_model = list_model, 
                      pri_adj_sec_calib_df = pri_adj_sec_calib_df,
                      lst = lst,
                      target_order = target_order)
  
  if (plot_flag){
    p1 = plot_line(data.frame(as.numeric(colnames(pri_calib_df$NIR)),
                                      primary_spec = pri_calib_df$NIR[1,],
                                      secondary_spec = sec_calib_df$NIR[1,],
                                      adj_secondary_spec = adj_sec_calib_df$NIR[1,]))
    plot_save(p1, "CSPDS_Spec",8,4)
    
    type = append(rep("Primary Calib", nrow(pri_calib_df)), 
                  rep("Secondary Calib", nrow(sec_calib_df)))
    original = rbind(pri_calib_df,sec_calib_df)
    p2 = plot_pca(prcomp(original$NIR, TRUE), 
                  color_group = type, 
                  shape_group = original[[lst]])
    
    
    type = append(rep("Primary Calib", nrow(pri_calib_df)), 
                  rep("Secondary_Adj Calib", nrow(sec_calib_df)))
    p3 = plot_pca(prcomp(pri_adj_sec_calib_df$NIR, TRUE), 
                  color_group = type, shape_group = pri_adj_sec_calib_df[[lst]])
    plot_save(cowplot::plot_grid(p2,p3,ncol = 2), "CSPDS_PCA", 8,4)
    
  }
  return(list(model = cspds_model, adj_sec_calib_df = adj_sec_calib_df))
}


test_cspds_model = function(test_df, cspds_model, pop_model, plot_flag = TRUE){
  wn_l = ncol(test_df$NIR)
  num_model = length(cspds_model$list_model)
  
  adj_test = matrix(nrow = nrow(test_df), ncol = wn_l)
  colnames(adj_test) = colnames(test_df$NIR)
  
  cms_mat = matrix(nrow = nrow(test_df), ncol = num_model)
  colnames(cms_mat) = cspds_model$target_order
  
  res = list()
  for (i in 1:nrow(test_df)){
    adj_list = lapply(cspds_model$list_model,
                      function(x){test_pds_model(test_df[i,], x, plot_flag = plot_flag)})
    cms_list = unlist(lapply(1:num_model, function(i){
      calc_cms(pop_model$med_list[[i]], adj_list[[i]]$NIR[1,])}))
    cms_mat[i,] = cms_list
    min_cms = min(cms_list)
    min_ind = which(cms_list==min_cms)
    adj_test[i,] = adj_list[[min_ind]]$NIR
    
    initial_pred = pop_model$target_order[min_ind]
    depth_score = calc_depth_score(min_cms, pop_model$pop_list[[min_ind]])
    res[[i]] = list(sno = i, min_cms = min_cms, min_ind = min_ind, 
                    initial_pred = initial_pred, gt= test_df[i,][[cspds_model$lst]],depth_score = depth_score)
    cat("\r",paste(round(i/nrow(test_df)*100)," %",sep = ""))
  }
  
  res_df = as.data.frame(Reduce(rbind, res))
  
  adj_test_df = test_df
  adj_test_df$NIR = adj_test
  
  if (plot_flag){
    
    p3 = plot_heatmap(as.data.frame(cms_mat))
    plot_save(p3, "heat_map")
    
  }
  return(list(test_df = test_df, adj_test_df = adj_test_df, cms_mat = cms_mat, res = res_df))
}

build_cscca_model = function(pri_calib_df,sec_calib_df,k,lst="class", plot_flag){
  #For every class, build model (Split Class before model)
  target_order = unique(pri_calib_df[[lst]])
  list_cs_ct = lapply(target_order, function(x){
    build_cca_model(pri_calib_df[pri_calib_df[[lst]] == x,],
                       sec_calib_df[sec_calib_df[[lst]] == x,],
                       k, plot_flag = FALSE)})
  
  list_model = lapply(list_cs_ct, function(x){x[[1]]})
  adj_sec_calib_df = Reduce(rbind,lapply(list_cs_ct, function(x) {x[[2]]}))
  pri_adj_sec_calib_df = rbind(pri_calib_df, adj_sec_calib_df)
  cscca_model = list(list_model = list_model, 
                      pri_adj_sec_calib_df = pri_adj_sec_calib_df,
                      lst = lst,
                      target_order = target_order)
  
  if (plot_flag){
    p1 = plot_line(data.frame(as.numeric(colnames(pri_calib_df$NIR)),
                                      primary_spec = pri_calib_df$NIR[1,],
                                      secondary_spec = sec_calib_df$NIR[1,],
                                      adj_secondary_spec = adj_sec_calib_df$NIR[1,]))
    plot_save(p1, "CSCCA_Spec_1",8,4)
    
    type = append(rep("Primary Calib", nrow(pri_calib_df)), 
                  rep("Secondary Calib", nrow(sec_calib_df)))
    original = rbind(pri_calib_df,sec_calib_df)
    p2 = plot_pca(prcomp(original$NIR, TRUE), 
                  color_group = type, 
                  shape_group = original[[lst]])
    
    
    type = append(rep("Primary Calib", nrow(pri_calib_df)), 
                  rep("Secondary_Adj Calib", nrow(sec_calib_df)))
    p3 = plot_pca(prcomp(pri_adj_sec_calib_df$NIR, TRUE), 
                  color_group = type, shape_group = pri_adj_sec_calib_df[[lst]])
    plot_save(cowplot::plot_grid(p2,p3,ncol = 2), "CSCCA_1_PCA", 4,8)
    
  }
  return(list(model = cscca_model, adj_sec_calib_df = adj_sec_calib_df))
}

test_cscca_model = function(test_df, cscca_model, pop_model, plot_flag = TRUE){
  wn_l = length(test_df$NIR[1,])
  num_model = length(cscca_model$list_model)
  
  adj_test = matrix(nrow = nrow(test_df), ncol = wn_l)
  colnames(adj_test) = colnames(test_df$NIR)
  
  cms_mat = matrix(nrow = nrow(test_df), ncol = num_model)
  colnames(cms_mat) = cscca_model$target_order
  
  res = list()
  for (i in 1:nrow(test_df)){
    adj_list = lapply(cscca_model$list_model,
                      function(x){test_cca_model(test_df[i,], x, plot_flag = FALSE)})
    cms_list = unlist(lapply(1:num_model, function(i){
      calc_cms(pop_model$med_list[[i]], adj_list[[i]]$NIR[1,])}))
    
    cms_mat[i,] = cms_list
    min_cms = min(cms_list)
    min_ind = which(cms_list==min_cms)
    adj_test[i,] = adj_list[[min_ind]]$NIR
    
    initial_pred = pop_model$target_order[min_ind]
    depth_score = calc_depth_score(min_cms, pop_model$pop_list[[min_ind]])
    res[[i]] = list(sno = i, min_cms = min_cms, min_ind = min_ind, 
                    initial_pred = initial_pred, gt= test_df[i,][[cscca_model$lst]],depth_score = depth_score)
    cat("\r",paste(round(i/nrow(test_df)*100)," %",sep=""))
  }
  
  res_df = as.data.frame(Reduce(rbind, res))
  
  adj_test_df = test_df
  adj_test_df$NIR = adj_test
  
  if (plot_flag){
    
    p3 = plot_heatmap(as.data.frame(cms_mat))
    plot_save(p3, "Transfer_CSCCA_heat_map")
  }
  
  return(list(test_df = test_df, adj_test_df = adj_test_df, cms_mat = cms_mat, res = res_df))
}

build_cspcpds_model = function(pri_calib_df,sec_calib_df,pds_mean_win,ncomp,npc,lst="class", plot_flag = TRUE){
  target_order = unique(pri_calib_df[[lst]])
  list_cs_ct = lapply(target_order, function(x){
    build_pcpds_model(pri_calib_df[pri_calib_df[[lst]] == x,],
                       sec_calib_df[sec_calib_df[[lst]] == x,],
                       pds_mean_win, ncomp,npc, plot_flag = FALSE)})
  
  list_model = lapply(list_cs_ct, function(x){x[[1]]})
  adj_sec_calib_df = Reduce(rbind,lapply(list_cs_ct, function(x) {x[[2]]}))
  pri_adj_sec_calib_df = rbind(pri_calib_df, adj_sec_calib_df)
  cspcpds_model = list(list_model = list_model, 
                        pri_adj_sec_calib_df = pri_adj_sec_calib_df,
                        lst = lst,
                        target_order = target_order)
  
  if (plot_flag){
    p1 = plot_line(data.frame(as.numeric(colnames(pri_calib_df$NIR)),
                                      primary_spec = pri_calib_df$NIR[1,],
                                      secondary_spec = sec_calib_df$NIR[1,],
                                      adj_secondary_spec = adj_sec_calib_df$NIR[1,]))
    plot_save(p1, "cspcpds_Spec",4,8)
    
    type = append(rep("Primary Calib", nrow(pri_calib_df)), 
                  rep("Secondary Calib", nrow(sec_calib_df)))
    original = rbind(pri_calib_df,sec_calib_df)
    p2 = plot_pca(prcomp(original$NIR, TRUE), 
                  color_group = type, 
                  shape_group = original[[lst]])
    
    
    type = append(rep("Primary Calib", nrow(pri_calib_df)), 
                  rep("Secondary_Adj Calib", nrow(sec_calib_df)))
    p3 = plot_pca(prcomp(pri_adj_sec_calib_df$NIR, TRUE), 
                  color_group = type, shape_group = pri_adj_sec_calib_df[[lst]])
    plot_save(cowplot::plot_grid(p2,p3,ncol = 2), "cspcpds_PCA", 4,8)
    
  }
  return(list(model = cspcpds_model, adj_sec_calib_df = adj_sec_calib_df))
}

test_cspcpds_model = function(test_df, cspcpds_model, pop_model, plot_flag = TRUE){
  wn_l = ncol(test_df$NIR)
  num_model = length(cspcpds_model$list_model)
  
  adj_test = matrix(nrow = nrow(test_df), ncol = wn_l)
  colnames(adj_test) = colnames(test_df$NIR)
  
  cms_mat = matrix(nrow = nrow(test_df), ncol = num_model)
  colnames(cms_mat) = cspcpds_model$target_order
  
  res = list()
  for (i in 1:nrow(test_df)){
    adj_list = lapply(cspcpds_model$list_model,
                      function(x){test_pcpds_model(test_df[i,], x, plot_flag = FALSE)})
    cms_list = unlist(lapply(1:num_model, function(i){
      calc_cms(pop_model$med_list[[i]], adj_list[[i]]$NIR[1,])}))
    cms_mat[i,] = cms_list
    min_cms = min(cms_list)
    min_ind = which(cms_list==min_cms)
    adj_test[i,] = adj_list[[min_ind]]$NIR
    
    initial_pred = pop_model$target_order[min_ind]
    depth_score = calc_depth_score(min_cms, pop_model$pop_list[[min_ind]])
    res[[i]] = list(sno = i, min_cms = min_cms, min_ind = min_ind, 
                    initial_pred = initial_pred, gt= test_df[i,][[cspcpds_model$lst]],depth_score = depth_score)
    cat("\r",paste(round(i/nrow(test_df)*100)," %",sep=""))
  }
  
  res_df = as.data.frame(Reduce(rbind, res))
  
  adj_test_df = test_df
  adj_test_df$NIR = adj_test
  
  if (plot_flag){
    
    p3 = plot_heatmap(as.data.frame(cms_mat))
    plot_save(p3, "heat_map")
    
  }
  return(list(test_df = test_df, adj_test_df = adj_test_df, cms_mat = cms_mat, res = res_df))
}

###### Plotting Functions ######

theme_Publication  =  function(base_size=16, base_family="serif") {
  ggthemes::theme_foundation(base_size=base_size, base_family=base_family) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold",
                                                      size = ggplot2::rel(1.2), hjust = 0.5),
                   text = ggplot2::element_text(),
                   panel.background = ggplot2::element_rect(colour = NA),
                   plot.background = ggplot2::element_rect(colour = NA),
                   panel.border = ggplot2::element_rect(colour = NA),
                   axis.title = ggplot2::element_text(size = ggplot2::rel(1)),
                   axis.title.y = ggplot2::element_text(angle=90,vjust =2),
                   axis.title.x = ggplot2::element_text(vjust = -0.2),
                   axis.text = ggplot2::element_text(), 
                   axis.line = ggplot2::element_line(colour="black"),
                   axis.ticks = ggplot2::element_line(),
                   panel.grid.major = ggplot2::element_blank(), #element_line(colour="#f0f0f0"),
                   panel.grid.minor = ggplot2::element_blank(),
                   legend.key = ggplot2::element_rect(colour = NA),
                   legend.text = ggplot2::element_text(size = ggplot2::rel(1), color = "black"),      
                   legend.position = "right",
                   legend.direction = NULL,         
                   legend.justification = "center",
                   legend.key.size= ggplot2::unit(1, "cm"),
                   legend.margin = ggplot2::margin(0, 0, 0, 0, "cm"),
                   legend.title = ggplot2::element_text(),
                   plot.margin=ggplot2::unit(c(10,5,5,5),"mm"),
                   strip.background=ggplot2::element_rect(colour="#f0f0f0",fill="#f0f0f0"),
                   strip.text = ggplot2::element_text(face="bold")
    )
}

plot_save = function(p, name, bh = 3, bw = 4){
  fname = paste0(plot_results_dir, "/", plot_counter,"_",name, ".tiff")
  cowplot::save_plot(fname, p,  base_height = bh, base_width = bw)
  assign("plot_counter", (plot_counter+1), envir = .GlobalEnv)
}

plot_save = function(p, name, bh = 3, bw = 4){
  fname = paste0(plot_results_dir, "/", name, ".png")
  cowplot::save_plot(fname, p,  base_height = bh, base_width = bw)
  assign("plot_counter", (plot_counter+1), envir = .GlobalEnv)
}

plot_line = function(df, pal = palette("default"), 
                     ll=c("Simple Plot", expression(paste("Wavenumber (", cm^{-1},")")),
                                      "Absorbance", "Specs")) {
  #using package reshape2, ggplot2
  col_df = colnames(df)
  col_df[1] = "x"
  colnames(df) = col_df
  df_m = reshape2::melt(df, id.vars = "x")
  df_m$value = as.numeric(df_m$value)
  p = ggplot2::ggplot(df_m, ggplot2::aes(x, value)) +
    ggplot2::geom_line(ggplot2::aes(colour = factor(variable)), size = 1) +
    ggplot2::scale_color_manual(values = pal)+
    theme_Publication() +
    ggplot2::labs(title = ll[1],
         x = ll[2],
         y = ll[3],
         color = ll[4])
  return(p)
}

plot_dot = function(df, pal = palette("default"), 
                    ll=c("Simple Plot", expression(paste("Wavenumber (", cm^{-1},")")),
                              "Absorbance", "Specs")) {
  #using package reshape2, ggplot2
  col_df = colnames(df)
  col_df[1] = "x"
  colnames(df) = col_df
  df_m = reshape2::melt(df, id.vars = "x")
  df_m$value = as.numeric(df_m$value)
  p = ggplot2::ggplot(df_m, ggplot2::aes(x, value)) +
    ggplot2::geom_point(ggplot2::aes(colour = factor(variable)), size = 1.25) +
    ggplot2::scale_color_manual(values = pal)+
    theme_Publication() +
    ggplot2::labs(title = ll[1],
                  x = ll[2],
                  y = ll[3],
                  color = ll[4])
  return(p)
}


plot_bar = function(df, pal = palette("default"), 
                    ll=c("Simple Plot", "Wavenumber", "Absorbance", "Specs")) {
  col_df = colnames(df)
  col_df[1] = "x"
  colnames(df) = col_df
  df_m = reshape2::melt(df, id.vars = "x")
  df_m$value = as.numeric(df_m$value)
  p = ggplot2::ggplot(df_m, ggplot2::aes(x = x, y = value, fill = factor(variable))) + 
    ggplot2::scale_fill_manual(values = pal)+
    ggplot2::geom_bar(stat = "identity", position = "dodge", color = "black") +
    theme_Publication() +
    ggplot2::labs(title = ll[1],
         x = ll[2],
         y = ll[3],
         fill = ll[4])
  return(p)
}


plot_scatter = function(df,quant=FALSE, pal = palette("default"),
                        ll = c("Simple Scatter Plot","Comp1","Comp2",
                                                "Color Factor", "Shape Factor")){
  colnames(df) = c("Comp1", "Comp2", "Col_gr", "Sh_gr")
  if (quant) {
    color_group = df[, 3]
  }
  else {
    color_group = factor(df[, 3])
  }
  shape_group = factor(df[, 4])
  sh_chars = c(0, seq(1:25), 33,35,36,37,38,42,43,63,64)
  p = ggplot2::ggplot(df, ggplot2::aes(x = Comp1, y = Comp2)) +
    ggplot2::geom_point(ggplot2::aes(shape = shape_group, color = color_group),
               size = 2.5, alpha = 0.8, stroke = 1.5) +
    ggplot2::geom_hline(yintercept = 0, lty = 2) +
    ggplot2::geom_vline(xintercept = 0, lty = 2) +
    ggplot2::scale_shape_manual(values = sh_chars[1:length(unique(shape_group))])
  if (quant) {
    #p1 = p + ggplot2::scale_color_gradient(low = "#0A823F", high = "#A9230E")
    p1 = p + ggplot2::scale_color_gradient(low = pal[1], high = pal[2])
  }
  else {
    #p1 = p + ggplot2::scale_color_manual(values = seq(1, length(unique(color_group))))
    p1 = p + ggplot2::scale_color_manual(values = pal)
  }
  p2 = p1 + theme_Publication() +
    ggplot2::labs(title = ll[1],x = ll[2], y = ll[3],color = ll[4], shape = ll[5])
  return(p2)
}


plot_pca = function(pca, c1=1, c2=2, quant=FALSE, pal = palette("default"),
                    color_group, shape_group,
                    ll = c("PCA Plot", "Colour Factor", "Shape Factor")) {
  eig_value = (pca$sdev)^2
  var_pc = eig_value[1:5] / sum(eig_value)
  pc_labels = sprintf("PC%d (%.2f%%)", 1:5, var_pc * 100)
  df = data.frame("PC1" = pca$x[, c1],
                  "PC2" = pca$x[, c2],
                  "Col_gr" = color_group,
                  "Sh_gr" = shape_group)
  p = plot_scatter(df, quant,pal,
                   c(ll[1], pc_labels[c1], pc_labels[c2], ll[2], ll[3])) +
    ggplot2::theme(axis.line = ggplot2::element_blank()) #+
  #guides(shape = guide_legend(nrow=2), color = guide_legend(nrow = 2))
  #plot_save(p, "PCA")
  return(p)
}

plot_heatmap = function(df, pal = palette("default"), ll =  c("Heatmap", "Samples", "Population")){
  df$row = rownames(df)
  melted_df = reshape2::melt(df, "row")
  p = ggplot2::ggplot(data = melted_df, ggplot2::aes(x = row, y = variable, fill = value)) +
    ggplot2::geom_tile(ggplot2::aes(fill = value), color = "white",lwd = 1.5,linetype = 1) +
    ggplot2::geom_text(ggplot2::aes(label = round(value, 4)))+
    ggplot2::scale_fill_gradient(low = pal[1], high = pal[2])+
    ggplot2::guides(fill = ggplot2::guide_colourbar(title = ll[1]))+
    ggplot2::labs(title = ll[[1]], x = ll[[2]], y = ll[[3]]) +
    theme_Publication() +
    ggplot2::theme(axis.line = ggplot2::element_blank())
  return(p)
}

plot_tile = function(df, pal = palette("default"), ll =  c("Tile", "True Class", "Predicted Class")){
  df$row = rownames(df)
  melted_df = reshape2::melt(df, "row")
  p = ggplot2::ggplot(data = melted_df, ggplot2::aes(x = row, y = variable)) +
    ggplot2::geom_tile(fill = pal, color = "white",lwd = 1.5,linetype = 1) +
    ggplot2::geom_text(ggplot2::aes(label = round(value, 4)), size = 5)+
    ggplot2::labs(title = ll[[1]], x = ll[[2]], y = ll[[3]]) + 
    theme_Publication(base_size = 16) +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 14)) +
    ggplot2::theme(axis.line = ggplot2::element_blank())+
    ggplot2::theme(axis.ticks =  ggplot2::element_blank())+
    ggplot2::theme(axis.title.x = ggplot2::element_blank(),
                   axis.title.y = ggplot2::element_blank())+
    ggplot2::theme(plot.margin = ggplot2::unit(c(0,0,0.5,0.5), "cm"))
  return(p)
}

plot_crossval = function(df, ll =  c("Type", "Test Batches", "Train Batches")){
  df$row = rownames(df)
  melted_df = reshape2::melt(df, "row")
  p = ggplot2::ggplot(data = melted_df, ggplot2::aes(x=row, y=variable)) +
    #ggplot2::geom_tile(ggplot2::aes(fill = value)) +
    ggplot2::geom_tile(fill = "white", color = "black", lwd = 1, linetype = 1) +
    ggplot2::geom_text(ggplot2::aes(label = round(value, 2)))+
    #scale_fill_gradient(low = "#FF0D0D", high = "#69B34C")+
    #scale_fill_gradient(low = "#ED2938", high = "#3BCA6D")+
    #ggplot2::scale_fill_gradient(low = "#EA4335", high = "#34A853")+
    ggplot2::labs(title=ll[[1]], x=ll[[2]], y = ll[[3]])+
    theme_Publication()+
    ggplot2::theme(axis.line = ggplot2::element_blank())
  return(p)
}

plot_summary = function(df, ll =  c("Type", "Standardization Models", 
                                    "Classification Models")){
  df$row = factor(rownames(df), levels = rownames(df))
  melted_df = reshape2::melt(df, "row")
  p = ggplot2::ggplot(data = melted_df, ggplot2::aes(x=row, y=variable)) +
    ggplot2::geom_tile(fill = "white", color = "black", lwd = 1, linetype = 1) +
    ggplot2::geom_text(ggplot2::aes(label = round(value, 2)), size = 6)+
    ggplot2::labs(title=ll[[1]], x=ll[[2]], y = ll[[3]])+
    theme_Publication(base_size = 20)+
    ggplot2::theme(axis.line = ggplot2::element_blank())
  return(p)
}
