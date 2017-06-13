#!/usr/bin/env Rscript

library(scales);
library(stats);
library(optparse);
library(data.table);
library(ggplot2);
# library(grDevices);

# library(ggsignif);

get_args <- function() 
{
  option_list <- list(
    make_option('--table', dest = 'table', action = 'store'),
    make_option('--output_dir', dest = 'output_dir', action = 'store'),
    make_option('--cutoff', dest = 'cutoff_method', action = 'store'),
    make_option('--log', dest = 'log', action = 'store')
  );
  
  return(parse_args(OptionParser(option_list = option_list)));
}


load_test_args <- function()
{
  # args$table <- '/home/charles/clones/IT_Miner/Out/Bacillus_subtilis/ASM904v1/Output_of_each_step/Step06-With_fake_complements_list.csv';
  # args$output_dir <- '/home/charles/clones/IT_Miner/Out/Bacillus_subtilis/ASM904v1/Figures';
  
  # args$table <- '/home/charles/clones/IT_Miner/Out/Escherichia_coli/U00096.2/Output_of_each_step/Step06-With_fake_complements_list.csv';
  # args$output_dir <- '/home/charles/clones/IT_Miner/Out/Escherichia_coli/U00096.2/Figures';
  
  args$table <- '/home/charles/clones/IT_Miner/Out/Salmonella_enterica/Output_of_each_step/Step06-With_fake_complements_list.csv';
  args$output_dir <- '/home/charles/clones/IT_Miner/Out/Salmonella_enterica/Figures';
  
  args$cutoff <- 'conservative';
  
  return(args);
}


print_keywords <- function()
{
  cat('\nSome key words will be used for sake of convenience in this analysis:\n')
  cat(' - single: predicted IT which has no predicted reverse complement\n');
  cat(' - fake  : faked reverse complement of a single IT\n');
  cat(' - least : predicted IT which is closer from its upstream gene than its predicted complement is\n');
  cat(' - most  : predicted IT which is farther from its upstream gene than its predicted complement is\n');
}


get_index <- function(tab)
{
  fake <- tab$`RNIE Mode` == 'Fake'; 
  pred <- !fake;
  
  all_least <- tab$`Compl. Dist. Class` == 'Least Distant';
  all_most <- !all_least;
  
  # True couples must have an id which begins with a number
  first_char_id <- substr(tab$`ID complementary couple`, 1, 1);
  pred_cmpl <- grepl('[[:digit:]]', first_char_id);
  pred_single <- pred != pred_cmpl;
  
  least <- pred_cmpl == TRUE & all_least == pred_cmpl; 
  most <- pred_cmpl == TRUE & all_most == pred_cmpl;
  
  return(data.table(pred, fake, all_least, all_most,
                    pred_cmpl, pred_single, least, most));
}


compare_distrib <- function(x, y, alternative = 'inequal', allowed_fluctuation = 0.25)
{
  # NB: To compare distributions, the KS test is not chosen because too sensitive
  # to size effect. A rather very simple test is applied, relying on the differences
  # between first quartile, median and third quartile... Kernel tests are too 
  # heavy in term of implementation for what we really want here. 
  
  key_perc_x <- log10(quantile(x, seq(from = 0.25, to = 0.75, by = 0.25)));
  key_perc_y <- log10(quantile(y, seq(from = 0.25, to = 0.75, by = 0.25)));
  diff <- key_perc_x - key_perc_y;
  
  if(alternative == 'inequal') {
    if(any(abs(diff) > allowed_fluctuation)) {
      return(TRUE);
    } else {
      return(FALSE);
    }
  } else if(alternative == 'less') {
    if(any(diff > 0) && any(abs(diff) < allowed_fluctuation)) {
      return(FALSE);
    } else {
      return(TRUE);
    }
  }
}

compute_cutoff <- function(likely_distances, unlikely_distances) 
{
  # the cutoff is defined as the distance which includes a max of likely distances
  # and exclude a max of unlikely distances
  
  asc_perc <- quantile(likely_distances, seq(from = 0.01, to = 0.99, by = 0.01));
  des_perc <- quantile(unlikely_distances, seq(from = 0.99, to = 0.01, by = -0.01));
  
  diff <- abs(asc_perc - des_perc);
  cutoff_idx <- which.min(diff);
  cutoff <- mean(c(asc_perc[cutoff_idx], des_perc[cutoff_idx]));
  
  return(round(cutoff));
}


get_cutoff <- function(tab, idxs, cutoff_method)
{
  # Initialization
  distance = as.numeric(tab$`termStart-UpGeneEnd`);
  outcome = vector(mode = "logical", length = 8); outcome <- TRUE;
  likely_idxs <- unlikely_idxs <- vector(mode = 'logical', length = nrow(tab));
  
  list_likely_idxs <- list(single = idxs$pred_single, 
                           least = idxs$least,
                           least_convergence = (idxs$least & tab$`Compl. Genomic Class` == 'Convergence'),
                           least_codirectionality = (idxs$least & tab$`Compl. Genomic Class` == 'Co-directionality'),
                           most_convergence = (idxs$most & tab$`Compl. Genomic Class` == 'Convergence'));
  list_unlikely_idxs <- list(fake = idxs$fake, 
                             most = idxs$most,
                             least_divergence = (idxs$least & tab$`Compl. Genomic Class` == 'Divergence'),
                             most_codirectionality = (idxs$most & tab$`Compl. Genomic Class` == 'Co-directionality'),
                             most_divergence = (idxs$most & tab$`Compl. Genomic Class` == 'Divergence'));
  
  # Requirement A (for every cutoff methods): single & least must have both lower upDistance than fake & most
  k = 1;
  for(i in 1:2) {
    for(j in 1:2) {
      current_outcome <- compare_distrib(distance[list_likely_idxs[[i]]], distance[list_unlikely_idxs[[j]]], alternative = 'less');
      outcome[k] <- current_outcome;
      k <- k + 1;
    }
  }
  if(any(outcome[1:4] == FALSE)) {
    return(list(exit_code = 'A', cutoff = NULL));
  }

  if(cutoff_method != 'inclusive') {
    for(i in 1:2) { likely_idxs[list_likely_idxs[[i]]] <- TRUE; }
    for(j in 1:2) { unlikely_idxs[list_unlikely_idxs[[j]]] <- TRUE; }
    conservative_cutoff <- compute_cutoff(distance[likely_idxs], distance[unlikely_idxs]);
    
    if(cutoff_method == 'conservative') {
      return(list(exit_code = 0, cutoff = conservative_cutoff));
    }
  }
  
  likely_idxs[] <- unlikely_idxs[] <- FALSE;
  if(cutoff_method != 'conservative') {
    # Requirement B: least in convergence & co-directionality must have both lower upDistance than least in divergence
    for(i in 3:4) {
      if(!is.null(list_unlikely_idxs[[3]])) {
        current_outcome <- compare_distrib(distance[list_likely_idxs[[i]]], distance[list_unlikely_idxs[[3]]], alternative = 'less');
        outcome[k] <- current_outcome;
      }
      k <- k + 1;
    }
    if(any(outcome[5:6] == FALSE)) {
      return(list(exit_code = 'B', cutoff = conservative_cutoff));
    }
    
    # Requirement C: most in convergence must have lower upDistance than both most in co-directionality & divergence
    for(j in 4:5) {
      if(!is.null(list_unlikely_idxs[[j]])) {
        current_outcome <- compare_distrib(distance[list_likely_idxs[[5]]], distance[list_unlikely_idxs[[j]]], alternative = 'less');
        outcome[k] <- current_outcome;
      }
      k <- k + 1;
    }
    if(any(outcome[7:8] == FALSE)) {
      return(list(exit_code = 'C', cutoff = conservative_cutoff));
    }
    
    for(i in 1:5) { if(i != 2) { likely_idxs[list_likely_idxs[[i]]] <- TRUE; }}
    for(j in 1:5) { if(j != 2) { unlikely_idxs[list_unlikely_idxs[[j]]] <- TRUE; }}
    inclusive_cutoff <- compute_cutoff(distance[likely_idxs], distance[unlikely_idxs]);
    
    if(cutoff_method == 'average') {
      return(list(exit_code = 0, cutoff = round(mean(c(conservative_cutoff, inclusive_cutoff)))));
    } else {
      return(list(exit_code = 0, cutoff = inclusive_cutoff));
    }
  }
}


print_cutoff_results <- function(cutoff_out)
{
  if(cutoff_out$exit_code == 0) {
    cat(sprintf('\tCutoff distance from upstream gene: %d nt\n', cutoff_out$cutoff));
  } else if(cutoff_out$exit_code == 'A') {
    print_keywords();
    cat('\nThe first condition of application of a cutoff distance is not satisfied:\n');
    cat(' - \"single\" & \"least\" do not seem to have both lower distance\nfrom upstream gene than \"fake\" & \"most\"\n');
  } else if(cutoff_out$exit_code == 'B') {
    print_keywords();
    cat('\nThe second condition of application of an inclusive cutoff distance is not satisfied:\n');
    cat(' - \"least\" in convergence & co-directionality do not seem to have both\n   lower distance from upstream gene than \"least\" in divergence\n');
    cat(sprintf('A conservative cutoff will be therefore applied: %d nt\n', cutoff_out$cutoff));
  } else {
    print_keywords();
    cat('\nThe third condition of application of an inclusive cutoff distance is not satisfied:\n');
    cat(' - \"most\" in convergence does not seem to have lower distance from\n   upstream gene than both \"most\" in co-directionality & divergence\n');
    cat(sprintf('A conservative cutoff will be therefore applied: %d nt\n', cutoff_out$cutoff));
  }
}


add_class_column <- function(tab, idxs)
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

tab <- fread(args$table, sep = '\t', header = TRUE); 
tab <- tab[`termStart-UpGeneEnd` > 0]; 

idxs <- get_index(tab);

cutoff_out <- get_cutoff(tab, idxs, args$cutoff);
print_cutoff_results(cutoff_out);

tab <- add_class_column(tab, idxs);

##############################################################################
# Boxplot 1
##############################################################################
filename <- file.path(args$output_dir, "1-Boxplot_distance_from_upgene_by_complement.svg")
svg(filename, width = 16, height = 10);

bp_tab <- data.table(class = tab$class, 
                     distance = tab$`termStart-UpGeneEnd`, 
                     group = tab$`Compl. Genomic Class`);
ymax <- max(bp_tab$`distance`);
cutoff_line <- data.frame(x = c(0,5), y = cutoff_out$cutoff, Cutoff = factor(paste('\n', cutoff_out$cutoff, 'nt\n')));

ggplot(bp_tab, aes(x = class, y = distance)) +
  ggtitle('Distance of IT from upstream gene') +
  geom_boxplot(notch = TRUE) +
  geom_jitter(shape = 16, position = position_jitter(0.15), aes(colour = group)) +
  xlab('') +
  ylab('Distance from upstream STOP codon (nt)') +
  stat_summary(fun.y = mean, geom = "point", colour = 'red', size = 3) +
  stat_summary(fun.data = mean_cl_normal, geom = 'errorbar', width = 0.5, colour = 'red') +
  scale_x_discrete(labels = c('SINGLE:\nIT with no predicted\ncomplement',
                              'FAKE:\nFake reverse complement\n to single IT',
                              'LEAST:\nIT closer from\nits upstream gene\nthan its complement is',
                              'MOST:\nIT farther from\nits upstream gene\nthan its complement is')) +
  scale_y_log10(breaks = trans_breaks('log10', function(x) 10^x), labels = trans_format('log10', math_format(10^.x))) +
  annotation_logticks(sides = 'l') +
  scale_colour_discrete(name = "Genomic Context\n",
                        labels = c('\nCo-directionality\n--->_||_--->\n',
                                   '\nConvergence\n--->_||_<---\n',
                                   '\nDivergence\n<---_||_--->\n')) +
  annotate('text', x = 1, y = as.numeric(ymax * 1.4), label = paste('N=', nrow(tab[idxs$pred_single]))) +
  annotate('text', x = 2, y = ymax * 1.4, label = paste('N=', nrow(tab[idxs$fake]), '\n', nrow(tab[idxs$pred_single])-nrow(tab[idxs$fake]), 'discarded (because intraCDS)')) +
  annotate('text', x = 3, y = ymax * 1.4, label = paste('N=', nrow(tab[idxs$least]))) +
  annotate('text', x = 4, y = ymax * 1.4, label = paste('N=', nrow(tab[idxs$most]))) +
  geom_vline(xintercept = 2.5, size = 1, linetype = 'dashed') +
  geom_line(aes(x, y, linetype = Cutoff), cutoff_line, size = 2, colour = 'red') +
  theme(legend.title = element_text(size = rel(1.5)), legend.text = element_text(size = rel(1.2))) +
  theme(axis.text.x = element_text(size = 15, margin = unit(c(5, 0, 0, 0), 'mm'))) +
  theme(axis.text.y = element_text(size = 13)) +
  theme(axis.title.y = element_text(size = rel(1.8), margin = unit(c(0, 5, 0, 0), 'mm'))) +
  theme(plot.title = element_text(hjust = 0.55, size = rel(2), margin = unit(c(0,0,7,0), 'mm'))) +
  theme(plot.margin = unit(c(1, 1, 1, 1), 'cm'));

junk <- dev.off();
