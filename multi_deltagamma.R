library(ape);
# You need to use the file read.nexus.R which is a modified version of the 
# one distributed with the package APE
source("read.nexus.R");



#########################################
# MULTI DELTAGAMMA PROCESSING 
#########################################

# parameters
# We will use as input a nexus file containing 1000 trees from 
# the posterior distibution of Bayesian analyses
# These trees were sampled randomly after discarding burnin
# loop through the 1000 trees

# 1. tree_files should be a vector with the names of the tree files to
#    process:   c("Lep_tree01.nex", "Lep_tree02.nex" and so on)
tree_files = c("Lep_10_01.trees", "Lep_10_02.trees", "Lep_10_03.trees");

# 2. number of trees contained in each tree file
number_trees = c(1:9);

# 3. vector of interval.width values to proccess on each file
interval.widths = c(1,2,3,4);


for( i in 1:length(tree_files) ) {
	cat("\n\n## Working on tree_file: ", tree_files[i]);
	pdf(file=paste(tree_files[i], ".pdf", sep=""));
	# plot cool histogram
	par(mfcol=c(2,2));
	for( j in 1:length(interval.widths) ) {
		cat("\n\t# Working on interval.width: ", interval.widths[j]);
		all_break_points <- run_delta.gamma(tree_files[i], number_trees,
										interval.widths[j]);
		# plot cool histogram
		hist(all_break_points, breaks=80, 
                    main=paste("interval.width", interval.widths[j]), xlab="Mya");
	}
	dev.off();
}
