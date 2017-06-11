#!/usr/bin/env Rscript

library(optparse);
library(data.table);
# library(ggplot2);
# library(grDevices);
# library(scales);
# library(stats);
# library(ggsignif);

get_args <- function() 
{
  option_list <- list(
    make_option('--table', dest = 'table', action = 'store'),
    make_option('--output_dir', dest = 'output_dir', action = 'store'),
    make_option('--log', dest = 'log', action = 'store')
  );
  
  return(parse_args(OptionParser(option_list = option_list)));
}


load_test_args <- function()
{
  args$table <- '/home/charles/clones/IT_Miner/Out/Bacillus_subtilis/ASM904v1/Output_of_each_step/Step06-With_fake_complements_list.csv'
  #args$table <- '/home/charles/clones/IT_Miner/Out/Escherichia_coli/U00096.2/Output_of_each_step/Step06-With_fake_complements_list.csv';
  args$output_dir <- '/home/charles/clones/IT_Miner/Out/Escherichia_coli/U00096.2/Figures';
  args$log <- '/home/charles/clones/IT_Miner/Out/Escherichia_coli/U00096.2/Figures/test_R.log';
  
  return(args);
}


print_keywords <- function()
{
  print('Some key words will be used for sake of convenience in this analysis:')
  print(' * single: real IT which has no predicted reverse complement');
  print(' * fake  : faked reverse complement of a single IT');
  print(' * least : real IT which is closer from its upstream gene than its complement is');
  print(' * most  : real IT which is farther from its upstream gene than its complement is');
}


get_index <- function(tab)
{
  fake <- tab$`RNIE Mode` == 'Fake'; 
  real <- !fake;
  
  least <- tab$`Compl. Dist. Class` == 'Least Distant';
  most <- !least;
  
  # True couples must have an id which begins with a number
  first_char_id <- substr(tab$`ID complementary couple`, 1, 1);
  real_cmpl <- grepl('[[:digit:]]', first_char_id);
  real_single <- real != real_cmpl;
  
  real_least <- real_cmpl == TRUE & least == real_cmpl; 
  real_most <- real_cmpl == TRUE & most == real_cmpl;
  
  return(data.table(real, fake, least, most,
                    real_cmpl, real_single, real_least, real_most));
}


compute_cutoff <- function(dist_best_grp, dist_worst_grp) 
{
  asc_perc <- quantile(dist_best_grp, seq(from = 0.01, to = 0.99, by = 0.01));
  des_perc <- quantile(dist_worst_grp, seq(from = 0.99, to = 0.01, by = -0.01));
  
  diff <- abs(asc_perc - des_perc);
  cutoff_idx <- which.min(diff);
  cutoff_distance <- mean(c(asc_perc[cutoff_idx], des_perc[cutoff_idx]));
  
  return(round(cutoff_distance));
}


get_cutoff <- function(distance, idxs)
{
  pval_tab <- data.table(i = character(5), j = character(5), h1 = character(5), pval = double(5));
  
  # Check 1: single & least must have both lower upDistance than fake & most
  list_i <- list(single = idxs$real_single, least = idxs$real_least);
  list_j <- list(fake = idxs$fake, most = idxs$real_most);
  k = 1;
  for(i in 1:2) {
    for(j in 1:2) {
      pval <- wilcox.test(distance[list_i[[i]]], distance[list_j[[j]]], alternative = 'less')$p.value;
      pval_tab[k] <- data.table(names(list_i[i]), names(list_j[j]), 'less', pval);
      k <- k + 1;
    }
  }
  if(any(pval_tab$pval > 0.05)) {
    return(list(exit_code = -1, cutoff = NULL, pval_tab = pval_tab));
  }
  
  # Check 2: fake & most must not follow significantly different distributions
  pval <- wilcox.test(distance[list_j[[1]]], distance[list_j[[2]]])$p.value;
  pval_tab[k] <- data.table(names(list_j[1]), names(list_j[2]), 'not_equal', pval);
  if(pval < 0.05) {
    return(list(exit_code = -2, cutoff = NULL, pval_tab = pval_tab));
  }

  # Compute Cutoff.
  small_distances <- c(distance[list_i[[1]]], distance[list_i[[2]]]); 
  long_distances <- c(distance[list_j[[1]]], distance[list_j[[2]]]);
  cutoff <- compute_cutoff(small_distances, long_distances)
  return(list(exit_code = 0, cutoff = cutoff, pval_tab = pval_tab));
}


print_cutoff_results <- function(cutoff_out)
{
  if(cutoff_out$exit_code == 0) {
    print(sprintf('Cutoff distance from upstream gene: %dnt', cutoff_out$cutoff));
  } else if(cutoff_out$exit_code == -1) {
    print('The first condition of application of a cutoff distance is not satisfied:');
    print('single & least do not have both lower distance from upstream gene than fake & most');
    print(cutoff_out$pval_tab[1:4,]);
  } else {
    print('The second condition of application of a cutoff distance is not satisfied:');
    print('fake & most follow significantly different distributions');
    print(cutoff_out$pval_tab[5,])
  }
}

make_boxplot_each_distrib <- function(tab, idxs, pval_tab) 
{ 
  # TO DO 
  foo <- NULL;
}
 
#############################################################################
# MAIN
#############################################################################
args <- get_args(); if(is.null(args$table)) { args <- load_test_args(); }

print_keywords();

tab <- fread(args$table, sep = '\t', header = TRUE); 
tab <- tab[`termStart-UpGeneEnd` >= 0]; 

idxs <- get_index(tab);

cutoff_out <- get_cutoff(as.numeric(tab$`termStart-UpGeneEnd`), idxs);
print_cutoff_results(cutoff_out);

