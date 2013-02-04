#Code to run GAM method from McInnes et al. 2011

library(ape);
library("RSvgDevice");
# You need to use the file read.nexus.R which is a modified version of the 
# one distributed with the package APE
source("read.nexus.R");


#David Orme actually wrote this code based on my initial work as a need way to plot the results of the gamma method, so the annotations are his. I can help explain if anything is unclear.

# ------------------------------------------------------------------------------
# function to do deltaGamma analysis
# read the function delta.gamma in memory
# @input: 
# @params: one fully resolved phylogenetic tree
# @params: interval.width value (integer) 
# @output: deltaGamma object
#
# if you want to change the x axis of your plot. Change the line 188 accordingly
delta.gamma <- function(phy, n.interval=10, windows.width=5, interval.width=NULL){
		
	# branching times - aren't always in sequence so need to make them
	bt <- sort(branching.times(phy), decreasing=TRUE)
	
	# order here is slightly confusing
	# - the bt values are time since present
	# - ape plots trees with zero at the root node	
	# quick check - nodes exactly at the present
	# - not clear what should happen here (do we add a zero length internode?)
	#   and it is only likely to happen on simulated trees so... abort
	#if(any(bt <= 0)) stop('Branching times found exactly at present - simulated tree?')
		
	# get the breaks between the intervals
	# - this is guaranteed to line up with the deepest node
	# - with n.intervals - also lines up at shallowest
	if(! is.null(interval.width)){
		# is it sensible?
		if(max(bt) < (2 * windows.width - 1) * interval.width) {
			stop('Interval width is too wide to permit even a single delta gamma, \n',
			     'given the number of intervals in the calculation window')
		} 
		# extend in this case to make sure the last interval runs
		# over the present and then trim back to the present day
		breaks <- seq(max(bt), 0 - interval.width, by=-interval.width)
		breaks <- ifelse(breaks > 0, breaks, 0)
	} 
	else {
		breaks     <- seq(max(bt), 0, length=n.interval + 1)
		interval.width <- breaks[1] - breaks[2]
	}

	nBreaks <- length(breaks)
	
	# get the time ranges in each window 
	# - also add a complete span for whole thing to get overall gamma
	win.start <- c(breaks[1], breaks[1:(nBreaks-windows.width)])
	win.end   <- c(breaks[nBreaks], breaks[(windows.width + 1):nBreaks])
	nWin <- length(win.start) - 1
	
	# get a logical matrix of which nodes fall in each window
	windows <- outer(win.start ,bt, '>=') & ! outer(win.end, bt, '>')
	
	# turn that into a list of nodes _within_ the window
	btByWin <- apply(windows, 1, function(X) bt[X])
	
	# now need to clip those into the window 
	# - what to do with nodes that fall _exactly_ at the start or end
	# - this is only inevitable at the windows on the first node (whole tree and first window)
	#   and the code simply stops this happening for those rows
	for(win in seq_along(btByWin)){
		btByWin[[win]] <- c(btByWin[[win]], win.end[win])
		if(win > 2) btByWin[[win]] <- c(win.start[win], btByWin[[win]])
	}	
	
	##  ## NB with this new structure, now quite easy to implement other (probably wronger) modes
	##  ## as it stands, using mode 4 (aka Lynsey's mode) - currently thought to be right
	##  ## calculation modes schematic:
	##  bt <- c(0.1672, 0.2341, 0.2926, 0.4318, 0.474, 0.5885, 0.6187, 0.6374, 0.6531, 0.7648, 1.3424)
	##  plotRug <- function(x, y, len=0.1, ...) arrows(x, y-len, x, y+len, code=0, ...)
  ##  
	##  plot(0,0, type='n', ylim=c(0.5,4.5), xlim=c(0,1.5), xlab='Time', ylab='Options')
  ##  
	##  arrows(c(bt[1], bt[1], bt[1],0), 1:4, c(bt[11],1,bt[10],1), 1:4, code=0, col='grey', lwd=2)
	##  arrows(c(bt[1], bt[1], bt[1],0), 1:4, c(bt[10],bt[10],bt[9],bt[10]), 1:4, code=0, lwd=2)
	##  abline(v=c(0,1), col='grey')
  ##  
	##  for(opt in 1:4){
	##  plotRug(bt, rep(opt,11))
	##  }
  ##  
	##  points(x=c(0,1,1), y=c(4,4,2), col='red')
	
	## check the number of nodes in each window
	nNode <- sapply(btByWin, length)
	if(any(nNode <= 2)){
		warning('The window width and node distribution give rise to windows that do not contain any real nodes.')
	}

	# gamma calculations

	# - internode intervals
	gByWin <- lapply(btByWin, function(X) abs(diff(X)))
	
	# number of lineages at start (m) and end (n) of intervals
	# - look at which internode interval each break falls in
	nLin <- findInterval(-breaks, -bt) +1
	m <- c(nLin[1], nLin[1:(nBreaks-windows.width)])
	n   <- c(nLin[nBreaks], nLin[(windows.width + 1):nBreaks])
	
	# sequences of numbers
	nByWin <- mapply(seq, m, n)
	
	# products of lineage number and internode differences
	ng <- mapply('*', nByWin, gByWin)
	
	# T - sum of n*g values within the window
	T <- sapply(ng, sum)
	
	# sum of all but the last of the cumulative sums of the ng within the window
	sum.cum.ng <- sapply(ng, function(X) sum(cumsum(X[-length(X)])))
	
	# iv) gamma!
	gamma <- ((1/(n-m) * sum.cum.ng) - (T/2)) / (T * sqrt(1/(12*(n-m))))
	
	# window summary dropping the row of the complete gamma
	winSum <- data.frame(win.start=win.start, win.end = win.end, m=m, n=n, n.nodes=nNode, gamma=gamma)
	overall <- winSum[1,]
	winSum  <- winSum[-1,]
	
	# focal interval summary
	focal.start <- breaks[c(-(1:(windows.width -1)), -((nBreaks-windows.width + 1):nBreaks))]
	# get the delta gammas for the focal intervals where it can be calculated
	delta.gamma <- with(winSum, gamma [windows.width:nWin] - gamma[1:(nWin - windows.width +1)])
	focalSum <- data.frame(focal.start=focal.start, focal.end = focal.start+interval.width, delta.gamma=delta.gamma)

	RET <- list(breaks=breaks, branching.times=bt, window.summary=winSum, focal.interval.summary=focalSum, overall.gamma=overall)
	attr(RET, 'windows.width') <- windows.width
	attr(RET, 'interval.width') <- interval.width
	
	return(RET)
}




# ------------------------------------------------------------------------------
# function to run the deltaGamma analysis over 1000 trees
# @input: 
# @params: NEXUS file containing 1000 trees
# @params: number_trees vector c(1:1000);
# @params: interval.width value (integer)
# @output: significant_bursts variable containing all splits estimated on every tree
#
run_delta.gamma <- function(trees_file, number_trees, interval.width){
    for( i in 1:length(number_trees) ) {
		cat("\nDoing tree", i)
		tree<-read.nexus(trees_file, count=i)#Here is a test tree that I used in my simulations built using b = 1, d = 0.1, M1 = 500, M2 = 1000, T1 = 55, T2 = 55
		
		# remove outgroups
		#tree <- drop.tip(tree, tip);
	
		tree<-multi2di(tree) #make the tree binary if needed
	
		tree$edge.length<-round(tree$edge.length, 1)
	
		lynsTreeDG  <- delta.gamma(tree, interval.width=interval.width) #running the function, to get the delta gammas you need to call lynsTreeDG$focal.interval.summary,
														#specifically lynsTreeDG$focal.interval.summary$delta.gamma
														# interval.width=1 by default
		## ---------------------------------------
		## add deltagama to my all_delta.gama list
		all_delta.gamma[[i]]  <- lynsTreeDG;
	
		#The following code is good for making nice plots
		brnch      <- lynsTreeDG$branching.times
		treeHeight <- max(brnch)
		overlap <- attr(lynsTreeDG, 'windows.width')
		focal <- attr(lynsTreeDG, 'interval.width')
		breaks <- lynsTreeDG$breaks
		nWin <- length(breaks) - 1
		rem <- 30   #rem corresponds to the number of My to remove at the beginning and end of the simulation
	
		ss<-which(lynsTreeDG$focal.interval.summary$focal.start < (treeHeight - rem))
		ee<-which(lynsTreeDG$focal.interval.summary$focal.start > rem)
		DGs<-lynsTreeDG$focal.interval.summary[min(ss):max(ee),]
		testing<-pnorm(DGs$delta.gamma, mean=0,sd=1.96)
	
		#which(testing < 0.05) + rem  #prints out the significant intervals
#		plot(delta.gamma ~ I(treeHeight-focal.start ), data=lynsTreeDG$focal.interval.summary, bty="l",las=1,cex.lab=1.4,xlab="",ylab="",type="n",xlim=c(200,0))  
		#points(delta.gamma ~ I(treeHeight-focal.start ), data=lynsTreeDG$focal.interval.summary,pch=20,col="darkgrey")   #NB the intervals outside the cutoffs should probably not be analysed...
	
		#abline(h=0)  
#		abline(v=c(rem, lynsTreeDG$overall.gamma$win.start-rem))  
		#highlighting the significant points
		#points((treeHeight-(DGs$focal.start[which(testing < 0.05)])),DGs$delta.gamma[which(testing < 0.05)],col="white",pch=20)
#		points((treeHeight-(DGs$focal.start[which(testing < 0.05)])),DGs$delta.gamma[which(testing < 0.05)],lwd=0.5,pch=25)
	
		## ---------------------------------------
		significant_bursts[[i]] <- treeHeight - DGs$focal.start[which(testing < 0.05)];
	
		abline(v=c(rem, lynsTreeDG$overall.gamma$win.start-rem))  
	}
	return(significant_bursts);
}
	
## function multi_deltagamma to run deltagamma on a group of tree_files
## using several interval.widths and plotting the histograms nicely
run_multi_deltagamma <- function(tree_files, number_trees, interval.widths) {
    for( i in 1:length(tree_files) ) {
        devSVG(file=paste(tree_files[i], ".svg", sep=""));
        plot.new();
        all_significant_bursts <- list();
        cat("\n\n##Working on treefile: ", tree_files[i])
        
        for( j in 1:length(interval.widths) ) {
            cat("\n\n## Working on interval.width: ", interval.widths[j])
            significant_bursts <- run_delta.gamma(tree_files[i], number_trees, 
                                      interval.width=interval.widths[j])
            all_break_poninst <- c();
            for( l in 1:length(significant_bursts)) {
                all_break_points <- c(all_break_points, significant_bursts[[l]])
            }
            
            # plot cool histogram
            hist(all_break_points, 
                    breaks=80, 
                    main=paste("interval.width", interval.widths[j]), 
                    xlab="Mya");
        }
        dev.off();
    }
}

# save all deltaGamma results in disk for later use
# to open these files in R use load("filename.txt");
save(all_delta.gamma, file="all_delta_gama.txt", ascii=TRUE);
save(significant_bursts, file="significant_bursts.txt", ascii=TRUE);
	
print("All the delta.gamma values per tree are in the variable: all_delta.gamma");
print("All significant bursts per tree are in the variable: significant_bursts");

all_break_points <- c();


for (i in 1:length(significant_bursts) ) {
    all_break_points <- c(all_break_points, significant_bursts[[i]])
}

unique_break_points <- unique(all_break_points);



#########################################
# PROCESSING 
#########################################

# parameters
# We will use as input a nexus file containing 1000 trees from 
# the posterior distibution of Bayesian analyses
# These trees were sampled randomly after discarding burnin
# loop through the 1000 trees

# 1. tree_files should be a vector with the names of the tree files to
#    process:   c("Lep_tree01.nex", "Lep_tree02.nex" and so on)
tree_files = c("Lep_10_01.trees", "Lep_10_02.trees", "Lep_10_03.trees",
                "Lep_10_04.trees", "Lep_10_05.trees");

# 2. number of trees contained in each tree file
number_trees = c(1:9);

# 3. vector of interval.width values to proccess on each file
interval.widths = c(1,2,3,4);
all_delta.gamma <- list();
significant_bursts <- list();

# plot cool histogram
par(mfcol=c(2,3))
run_multi_deltagamma(tree_files, number_trees, interval.widths);
