#!/usr/bin/env Rscript

#############################################################################
# LOAD LIBRARIES
#############################################################################
library(optparse);
library(data.table);
# library(ggplot2);
# library(grDevices);
# library(scales);
# library(stats);
# library(ggsignif);


#############################################################################
# READ & EXTRACT THE OPTIONS
#############################################################################
option_list <- list(
  make_option('--table', dest = 'table', action = 'store'),
  make_option('--output_dir', dest = 'output_dir', action = 'store'),
  make_option('--log', dest = 'log', action = 'store')
);
opt <- parse_args(OptionParser(option_list = option_list));

# For testing
if(1 == 1 && is.null(opt$list_term)) {
  opt$table <- "/home/charles/clones/IT_Miner/Out/Escherichia_coli/U00096.2/Output_of_each_step/Step06-With_fake_complements_list.csv";
  opt$output_dir <- "/home/charles/clones/IT_Miner/Out/Escherichia_coli/U00096.2/Figures";
  opt$log <- "/home/charles/clones/IT_Miner/Out/Escherichia_coli/U00096.2/Figures/test_R.log";
}

##############################################################################
# LOAD THE LIST OF TERMINATORS 
##############################################################################
tab <- fread(opt$table, sep = '\t', header = TRUE);

##############################################################################
# Get idx of ITs based on:
#   1) if they are real or faked 
#   2) if the real ones are predicted to have a reverse complement or not
#   3) if partner in real couples are either Least or Most Distant from upGene
##############################################################################
fake_idx <- tab$`RNIE Mode` == 'Fake';
true_idx <- tab$`RNIE Mode` != 'Fake';

# True couples must have an id which begins with a number
first_char_id <- substr(tab$`ID complementary couple`, 1, 1);
cmpl_idx <- grepl('[[:digit:]]', first_char_id);
single_idx <- true_idx != cmpl_idx;

cmpl_least_idx
# least_idx <- which(tab$`Compl. Dist. Class` == "Least Distant");
# most_idx <- which(tab$`Compl. Dist. Class` == "Most Distant");
# compl_idx <- c(least_idx, most_idx);

##############################################################################
# Boxplot of the complements / groups based on distance / color on geno_class
##############################################################################
boxplot_tab <- data.table(class = tab[compl_idx]$`RNIE Mode`,
                          distance = tab[compl_idx]$`termStart-UpGeneEnd`,
                          color_grp = tab[compl_idx]$`Compl. Genomic Class`)

file <- file.path(opt$output_dir, "08-Boxplot_fake_based_on_geno_context");
svg(file, width = 12, height = 10);

give.n <- function(x){
  return(c(y = mean(x), label = length(x)));
}

distance = boxplot_tab$distance;
ymax = max(distance);

ggplot(boxplot_tab, aes(x = class, y = distance)) +
  ggtitle("Distance of terminator from upstream gene") +
  geom_boxplot(notch = TRUE) +
  geom_jitter(shape = 16, position = position_jitter(0.15), aes(colour = color_grp)) +
  #scale_color_gradient2(name = "Couple ID", midpoint = 10, low="darkgreen", mid = "yellow", high = "red") +
  #scale_color_gradient(low="darkgreen", high = "red") +
  stat_summary(fun.y = mean, geom = "point", colour = "red", size = 3) +
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", width = 0.5, colour = "red") +
  annotate("text", x = 1, y = ymax*1.4, label = paste("N=", length(distance)/2)) +
  annotate("text", x = 2, y = ymax*1.4, label = paste("N=", length(distance)/2)) +
  scale_x_discrete(labels = c("Fake\nComplement", "Gene\nComplement")) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) +
  annotation_logticks(sides = "l") +
  xlab("") +
  ylab("Distance from upstream gene (nt)") +
  theme(axis.text.x = element_text(size = 15, margin = unit(c(5,0,0,0), "mm"))) +
  theme(axis.text.y = element_text(size = 13)) +
  theme(axis.title.y = element_text(size = rel(1.8), margin = unit(c(0,5,0,0), "mm"))) +
  theme(plot.title = element_text(hjust = 0.3, size = rel(2), margin = unit(c(0,0,7,0), "mm"))) +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))
junk <- dev.off();

# ##############################################################################
# # Retrieve Index of complemented terminators
# ##############################################################################
# non_compl_idx <- which(tab$`Reverse Complement Attribute` == "");
# fake_idx <- non_compl_idx[which(tab[non_compl_idx]$`Genome` == "")];
# true_idx <- non_compl_idx[which(tab[non_compl_idx]$`Genome` != "")];
# compl_idx <- which(tab$`Reverse Complement Attribute` != "");
# 
# ##############################################################################
# # Create a column "class"
# #############################################################################
# class_col <- vector(mode = "character", length = nrow(tab))
# class_col[compl_idx] <- "Complements";
# class_col[true_idx] <- "No-Complement";
# class_col[fake_idx] <- "Simulated Complement";
# 
# ##############################################################################
# # Gradient of color depending on distance
# ##############################################################################
# max_distance <- max(tab$`Term. start dist. from stop codon`);
# n_colors <- 20;
# scale <- 10^(seq(from = 0, to = log10(max_distance+1), length.out = n_colors));
# color_grp <- unlist(lapply(tab$`Term. start dist. from stop codon`, 
#                              function(x) min(which(x < scale)) - 1));
# color_grp[fake_idx] <- color_grp[true_idx];
# 
# ##############################################################################
# # Boxplot of the non_compl term (& their simul. compl)
# ##############################################################################
# boxplot_tab <- data.table(class = class_col[non_compl_idx],
#                           distance = tab[non_compl_idx]$`Term. start dist. from stop codon`,
#                           color_grp = color_grp[non_compl_idx],
#                           complement_class = tab[non_compl_idx]$`Complement Class`);
# file <- file.path(opt$output_dir, "01-Boxplot_non_compl_distances_from_upstream_gene.svg");
# svg(file, width = 8, height = 10);
# give.n <- function(x){
#   return(c(y = mean(x), label = length(x)));
# }
# distance = boxplot_tab$distance;
# ymax = max(distance);
# ggplot(boxplot_tab, aes(x = class, y = distance)) +
#   ggtitle("Distance of terminator from upstream gene") +
#   geom_boxplot(notch = TRUE) +
#   geom_jitter(shape = 16, position = position_jitter(0.15), aes(colour = color_grp)) +
#   scale_color_gradient2(name = "Couple ID", midpoint = 10, low="darkgreen", mid = "yellow", high = "red") +
#   #scale_color_gradient(low="darkgreen", high = "red") +
#   stat_summary(fun.y = mean, geom = "point", colour = "red", size = 3) +
#   stat_summary(fun.data = mean_cl_normal, geom = "errorbar", width = 0.5, colour = "red") +
#   annotate("text", x = 1, y = ymax*1.4, label = paste("N=", length(distance)/2)) +
#   annotate("text", x = 2, y = ymax*1.4, label = paste("N=", length(distance)/2)) + 
#   scale_x_discrete(labels = c("Terminators with no\npredicted complement", "Simulated Complement")) +
#   scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) +
#   annotation_logticks(sides = "l") +
#   xlab("") +
#   ylab("Distance from upstream gene (nt)") +
#   theme(axis.text.x = element_text(size = 15, margin = unit(c(5,0,0,0), "mm"))) +
#   theme(axis.text.y = element_text(size = 13)) +
#   theme(axis.title.y = element_text(size = rel(1.8), margin = unit(c(0,5,0,0), "mm"))) +
#   theme(plot.title = element_text(hjust = 0.3, size = rel(2), margin = unit(c(0,0,7,0), "mm"))) +
#   theme(plot.margin = unit(c(1,1,1,1), "cm"))
# junk <- dev.off();
# 
# ##############################################################################
# # Boxplot non_compl vs compl
# ##############################################################################
# real_compl_class <- tab$`Complement Class`;
# real_compl_class[non_compl_idx] <- NA;
# boxplot_tab <- data.table(class = class_col[c(true_idx, compl_idx)],
#                           distance = tab[c(true_idx, compl_idx)]$`Term. start dist. from stop codon`,
#                           legend_grp = real_compl_class[c(true_idx, compl_idx)]);
# boxplot_tab$class <- factor(boxplot_tab$class, levels = c("No-Complement", "Complements"), ordered = T);
# 
# file <- file.path(opt$output_dir, "02-Boxplot_non_compl_vs_compl_distances_from_upstream_gene.svg");
# svg(file, width = 8, height = 10);
# give.n <- function(x){
#   return(c(y = mean(x), label = length(x)));
# }
# ymax = max(distance);
# ggplot(boxplot_tab, aes(x = class, y = distance)) +
#   ggtitle("Distance of terminator from upstream gene depending on\nwhether or not predicted to have a complementary terminator") +
#   geom_boxplot(notch = TRUE) +
#   geom_jitter(shape = 16, position = position_jitter(0.17), aes(color = legend_grp)) +
#   scale_colour_discrete(name = "Within each couple:", breaks = c("Least Distant", "Most Distant"),
#                         labels = c("Least distant", "Most distant")) +
#   theme(legend.title = element_text(size = rel(1)), legend.text = element_text(size = rel(1))) +
#   stat_summary(fun.y = mean, geom = "point", colour = "red", size = 3) +
#   stat_summary(fun.data = mean_cl_normal, geom = "errorbar", width = 0.5, colour = "red") +
#   annotate("text", x = 1, y = ymax*1.4, label = paste("N=", length(boxplot_tab$distance[true_idx]))) +
#   annotate("text", x = 2, y = ymax*1.4, label = paste("N=", length(boxplot_tab$distance[compl_idx])))  +
#   scale_x_discrete(labels = c("Terminators with no\npredicted complement", "Complementary\nterminators")) +
#   scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) +
#   annotation_logticks(sides = "l") +
#   xlab("") +
#   ylab("Distance from upstream gene (nt)") +
#   theme(axis.text.x = element_text(size = 15, margin = unit(c(5,0,0,0), "mm"))) +
#   theme(axis.text.y = element_text(size = 13)) +
#   theme(axis.title.y = element_text(size = rel(1.8), margin = unit(c(0,5,0,0), "mm"))) +
#   theme(plot.title = element_text(hjust = 0.3, size = rel(1.5), margin = unit(c(0,0,7,0), "mm"))) +
#   theme(plot.margin = unit(c(1,1,1,1), "cm"))
# junk <- dev.off();
# 
# ##################################################################
# # Comparisons of distribution
# ##################################################################
# least_idx <- which(real_compl_class == "Least Distant");
# most_idx <- which(real_compl_class == "Most Distant");
# new_class_col <- class_col;
# new_class_col[least_idx] <- "Least Distant Complement";
# new_class_col[most_idx] <- "Most Distant Complement";
# true_vs_fake <- wilcox.test(tab[true_idx]$`Term. start dist. from stop codon`, 
#                              tab[fake_idx]$`Term. start dist. from stop codon`)$p.value;
# true_vs_least <- wilcox.test(tab[true_idx]$`Term. start dist. from stop codon`, 
#                              tab[least_idx]$`Term. start dist. from stop codon`)$p.value;
# true_vs_most <- wilcox.test(tab[true_idx]$`Term. start dist. from stop codon`, 
#                              tab[most_idx]$`Term. start dist. from stop codon`)$p.value;
# fake_vs_least <- wilcox.test(tab[fake_idx]$`Term. start dist. from stop codon`, 
#                              tab[least_idx]$`Term. start dist. from stop codon`)$p.value;
# fake_vs_most <- wilcox.test(tab[fake_idx]$`Term. start dist. from stop codon`, 
#                             tab[most_idx]$`Term. start dist. from stop codon`)$p.value;
# 
# ##################################################################
# # All boxplots
# ##################################################################
# boxplot_tab <- data.table(class = new_class_col,
#                           distance = tab$`Term. start dist. from stop codon`);
# boxplot_tab$class <- factor(boxplot_tab$class, 
#                             levels = c("No-Complement", "Simulated Complement", "Least Distant Complement", 
#                                        "Most Distant Complement"), ordered = T);
# 
# file <- file.path(opt$output_dir, "03-Boxplot_all_distances_from_upstream_gene.svg");
# svg(file, width = 16, height = 10);
# give.n <- function(x){
#   return(c(y = mean(x), label = length(x)));
# }
# ymax = max(distance);
# bp <- ggplot(boxplot_tab, aes(x = class, y = distance)) +
#   ggtitle("Comparisons of distributions of distances\nbetween terminators and upstream genes") +
#   geom_boxplot(notch = TRUE, aes(colour = class)) +
#   geom_jitter(shape = 16, position = position_jitter(0.17)) +
#   # scale_colour_discrete(name = "Within each couple:", breaks = c("Least Distant", "Most Distant"),
#   #                       labels = c("Least distant", "Most distant")) +
#   # theme(legend.title = element_text(size = rel(1)), legend.text = element_text(size = rel(1))) +
#   stat_summary(fun.y = mean, geom = "point", colour = "red", size = 3) +
#   stat_summary(fun.data = mean_cl_normal, geom = "errorbar", width = 0.5, colour = "red") +
#   annotate("text", x = 1, y = ymax*1.4, label = paste("N=", length(boxplot_tab$distance[true_idx]))) +
#   annotate("text", x = 2, y = ymax*1.4, label = paste("N=", length(boxplot_tab$distance[fake_idx])))  +
#   annotate("text", x = 3, y = ymax*1.4, label = paste("N=", length(boxplot_tab$distance[least_idx]))) +
#   annotate("text", x = 4, y = ymax*1.4, label = paste("N=", length(boxplot_tab$distance[most_idx])))  +
#   scale_x_discrete(labels = c("Terminator with no\npredicted complement", 
#                               "Simulated\nComplementary Terminator",
#                               "Least Distant\nComplementary Terminator",
#                               "Most Distant\nComplementary Terminator")) +
#   scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) +
#   annotation_logticks(sides = "l") +
#   xlab("") +
#   ylab("Distance from upstream gene (nt)") +
#   theme(axis.text.x = element_text(size = 15, margin = unit(c(5,0,0,0), "mm"))) +
#   theme(axis.text.y = element_text(size = 13)) +
#   theme(axis.title.y = element_text(size = rel(1.8), margin = unit(c(0,5,0,0), "mm"))) +
#   theme(plot.title = element_text(hjust = 0.5, size = rel(2), margin = unit(c(0,0,7,0), "mm"))) +
#   theme(plot.margin = unit(c(1,1,1,1), "cm"));
# 
# bar_true_least <- data.frame(a = c(1, 1, 3, 3), 
#                             b = c(100000,120000,120000,100000));
# bar_fake_most <- data.frame(a = c(2, 2, 4, 4), 
#                             b = c(200000,240000,240000,200000));
# bp + geom_line(data = bar_true_least, aes(x = a, y = b)) + 
#   annotate("text", x = 2, y = 150000, label = "s. (pval= 1.04e-11)") + 
#   geom_line(data = bar_fake_most, aes(x = a, y = b)) + 
#   annotate("text", x = 3, y = 300000, label = "n.s. (pval= 0.64)")
# junk <- dev.off();
# 
# 
# ##############################################################################
# # Gradient of color depending on distance
# ##############################################################################
# max_distance <- max(tab[least_idx]$`Term. start dist. from stop codon`);
# n_colors <- 20;
# scale <- 10^(seq(from = 0, to = log10(max_distance+1), length.out = n_colors));
# color_grp <- vector(mode = "integer", length = length(compl_idx))
# color_grp_least <- unlist(lapply(tab[least_idx]$`Term. start dist. from stop codon`, 
#                            function(x) min(which(x < scale)) - 1));
# boxplot_tab <- data.table(distance = tab[compl_idx]$`Term. start dist. from stop codon`,
#                           class = new_class_col[compl_idx],
#                           id = tab[compl_idx]$`ID couple`)
# new_least_idx <- which(boxplot_tab$`class` == "Least Distant Complement");
# new_most_idx <- which(boxplot_tab$`class` == "Most Distant Complement");
# color_grp[new_least_idx] = color_grp_least;
# for(i in 1:length(color_grp)) {
#   if(color_grp[i] == 0) {
#     color_grp[i] = color_grp[i-1]
#   }
# }
# boxplot_tab = cbind(boxplot_tab, color_grp);
# 
# 
# ##############################################################################
# # Boxplot of the non_compl term (& their simul. compl)
# ##############################################################################
# file <- file.path(opt$output_dir, "04-Boxplot_compl_distances_from_upstream_gene.svg");
# svg(file, width = 8, height = 10);
# give.n <- function(x){
#   return(c(y = mean(x), label = length(x)));
# }
# distance = boxplot_tab$distance;
# ymax = max(distance);
# ggplot(boxplot_tab, aes(x = class, y = distance)) +
#   ggtitle("Distance of terminator from upstream gene") +
#   geom_boxplot(notch = TRUE) +
#   geom_jitter(shape = 16, position = position_jitter(0.15), aes(colour = color_grp)) +
#   scale_color_gradient2(name = "Couple ID", midpoint = 10, low="darkgreen", mid = "yellow", high = "red") +
#   #scale_color_gradient(low="darkgreen", high = "red") +
#   stat_summary(fun.y = mean, geom = "point", colour = "red", size = 3) +
#   stat_summary(fun.data = mean_cl_normal, geom = "errorbar", width = 0.5, colour = "red") +
#   annotate("text", x = 1, y = ymax*1.4, label = paste("N=", length(distance)/2)) +
#   annotate("text", x = 2, y = ymax*1.4, label = paste("N=", length(distance)/2)) + 
#   scale_x_discrete(labels = c("Least\nDistant Complement", "Most\nDistant Complement")) +
#   scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) +
#   annotation_logticks(sides = "l") +
#   xlab("") +
#   ylab("Distance from upstream gene (nt)") +
#   theme(axis.text.x = element_text(size = 15, margin = unit(c(5,0,0,0), "mm"))) +
#   theme(axis.text.y = element_text(size = 13)) +
#   theme(axis.title.y = element_text(size = rel(1.8), margin = unit(c(0,5,0,0), "mm"))) +
#   theme(plot.title = element_text(hjust = 0.3, size = rel(2), margin = unit(c(0,0,7,0), "mm"))) +
#   theme(plot.margin = unit(c(1,1,1,1), "cm"))
# junk <- dev.off();
# 
# ################################################################################
# # Find threshold (include the most best terminator, exclude the most worst terminator)
# ################################################################################
# best_idx <- c(true_idx, least_idx);
# worst_idx <- c(fake_idx, most_idx);
# best_perc <- quantile(tab[best_idx]$`Term. start dist. from stop codon`, 
#          seq(from = 0.01, to = 0.99, by = 0.01));
# worst_perc <- quantile(tab[worst_idx]$`Term. start dist. from stop codon`,
#          seq(from = 0.99, to = 0.01, by = -0.01));
# diff <- abs(best_perc - worst_perc);
# threshold <- which.min(diff);
# cutoff_distance <- round(mean(c(best_perc[threshold], c(worst_perc[threshold]))));
# 
# ################################################################################
# # Boxplot with threshold
# ################################################################################
# new_class_col[best_idx] <- "Non-complemented & Least Distant Complement";
# new_class_col[worst_idx] <- "Simulated & Most Distant Complement";
# boxplot_tab <- data.table(class = new_class_col,
#                           distance = tab$`Term. start dist. from stop codon`);
# file <- file.path(opt$output_dir, "05-Boxplot_cutoff_distance.svg");
# svg(file, width = 16, height = 10);
# give.n <- function(x){
#   return(c(y = mean(x), label = length(x)));
# }
# distance = boxplot_tab$distance
# ymax = max(distance);
# ggplot(boxplot_tab, aes(x = class, y = distance)) +
#   ggtitle("Cutoff Distance which includes a maximum of good predictions\nand exclude a maximum of bad predictions") +
#   geom_boxplot(notch = TRUE, aes(colour = class)) +
#   geom_jitter(shape = 16, position = position_jitter(0.17)) +
#   # scale_colour_discrete(name = "Within each couple:", breaks = c("Least Distant", "Most Distant"),
#   #                       labels = c("Least distant", "Most distant")) +
#   # theme(legend.title = element_text(size = rel(1)), legend.text = element_text(size = rel(1))) +
#   stat_summary(fun.y = mean, geom = "point", colour = "red", size = 3) +
#   stat_summary(fun.data = mean_cl_normal, geom = "errorbar", width = 0.5, colour = "red") +
#   annotate("text", x = 1, y = ymax*1.4, label = paste("N=", length(boxplot_tab$distance[best_idx]))) +
#   annotate("text", x = 2, y = ymax*1.4, label = paste("N=", length(boxplot_tab$distance[worst_idx])))  +
#   scale_x_discrete(labels = c("Non-complemented &\nLeast Distant Complementary\nTerminators", 
#                               "Simulated &\n Most Distant Complementary\nTerminators")) +
#   scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) +
#   annotation_logticks(sides = "l") +
#   xlab("") +
#   ylab("Distance from upstream gene (nt)") +
#   theme(axis.text.x = element_text(size = 15, margin = unit(c(5,0,0,0), "mm"))) +
#   theme(axis.text.y = element_text(size = 13)) +
#   theme(axis.title.y = element_text(size = rel(1.8), margin = unit(c(0,5,0,0), "mm"))) +
#   theme(plot.title = element_text(hjust = 0.5, size = rel(2), margin = unit(c(0,0,7,0), "mm"))) +
#   theme(plot.margin = unit(c(1,1,1,1), "cm")) +
#   geom_hline(yintercept = cutoff_distance, colour = "red", size = 2)
#   # theme(panel.background = element_rect(fill = 'green', colour = 'red')) +
#   # theme(axis.line = element_blank());
# junk <- dev.off()

