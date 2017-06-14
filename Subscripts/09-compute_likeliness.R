#!/usr/bin/env Rscript

library(scales);
library(stats);
library(optparse);
library(data.table);
library(ggplot2);

get_args <- function() 
{
  option_list <- list(
    make_option('--table', dest = 'table', action = 'store'),
    make_option('--fig_dir', dest = 'fig_dir', action = 'store'),
    make_option('--out_file', dest = 'out_file', action = 'store')
  );
  
  return(parse_args(OptionParser(option_list = option_list)));
}


load_test_args <- function()
{
  # args$table <- '/home/charles/clones/IT_Miner/Out/Bacillus_subtilis/ASM904v1/Output_of_each_step/Step06-With_fake_complements_list.csv';
  # args$fig_dir <- '/home/charles/clones/IT_Miner/Out/Bacillus_subtilis/ASM904v1/Figures';
  
  # args$table <- '/home/charles/clones/IT_Miner/Out/Escherichia_coli/U00096.2/Output_of_each_step/Step06-With_fake_complements_list.csv';
  # args$fig_dir <- '/home/charles/clones/IT_Miner/Out/Escherichia_coli/U00096.2/Figures';
  
  args$table <- '/home/charles/clones/IT_Miner/Out/Salmonella_enterica/Output_of_each_step/Step08-Filtered_list.csv';
  args$fig_dir <- '/home/charles/clones/IT_Miner/Out/Salmonella_enterica/Figures';
  args$out_file <- '/home/charles/clones/IT_Miner/Out/Salmonella_enterica/Output_of_each_step/final_list.csv';
  
  return(args);
}

compute_bit_coeff <- function(bit_score)
{
  max_bit <- round(max(bit_score));
  max_bit_coeff <- 1;
  
  bit_score_levels <- levels(factor(round(bit_score)));
  bit_score_coeff <- log2(seq(from = 2^(0), to = 2^(max_bit_coeff), 
                               length.out = max_bit + 1));
  bit_score_coeff <- bit_score_coeff[as.integer(bit_score_levels) + 1];
  
  bit_tab <- data.table(bit_score_levels, bit_score_coeff);
  return(bit_tab)
}

add_class_column <- function()
{
  class <- vector(mode = 'character', length = nrow(idxs));
  class[idxs$pred_single] <- 'Single';
  class[idxs$fake] <- 'Fake';
  class[idxs$least] <- 'Least';
  class[idxs$most] <- 'Most';
  class <- factor(class, levels = c('Single', 'Fake', 'Least', 'Most'), ordered = TRUE);
  
  tab <- cbind(tab, class);
  return(tab);
}

#############################################################################
# MAIN
#############################################################################
args <- get_args(); if(is.null(args$table)) { args <- load_test_args(); }
out_file <<- args$out_file;

tab <<- fread(args$table, sep = '\t', header = TRUE);
tab$`Bit Score` <- round(tab$`Bit Score`);

bit_tab <- compute_bit_coeff(tab$`Bit Score`);
vector_bit_coeff <- bit_tab[levels(factor(tab$`Bit Score`))]$bit_score_coeff

N <- tab[, .N, by = `termStart-UpGeneEnd`]$N
sum_bit_coeff <- tab[, function (x) sum(x), by = `termStart-UpGeneEnd`]

