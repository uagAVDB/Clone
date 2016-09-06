#####################################################################
# 																	#
#		LTR5'-LTR3' detection: Merging into single positions		#
#																	#	
#																	#
#####################################################################


## Step 1. : Load the data
## Step 2. : Prepare three GRanges objects [IS.3LTR, IS.LTR5 and a window.gr]
## Step 3. : Separate the IS.LTR into four categories [LTR3 alone, LTR5 alone, LTR3LTR5, LTR3LTRXLTRX]
## Step 4. : Report a file containing the new positions, a stat.txt file and an error file containing LTR3LTRXLTRX positions

## Variables 
# IStable = a table containing the IS information
# name = Name of the output
# window.up = Search for IS opposite to LTR3 in a window of ... (upstream) (default = 20)
# window.down = Search for IS opposite to LTR3 in a window of ... (downstream) (default = 20)
# output = Ex: "/path/to/output/files"

# NB: Assume that IRanges package is installed

mergeISbyLTR <- function(IStable, name, window.up=10, window.down=10, output){
	
	require(IRanges)
	require(GenomicRanges)

	print(paste(" ## 1. Prepare GRanges objects: ", name, " ##", sep=""))
		
	viral.pos <- GRanges()	# Output
	error <- GRanges()

	if(nrow(IStable) >= 1){ 
	
		IS <- makeGRangesFromDataFrame(data.frame(
		seqnames=IStable$chr, 
		start=c(IStable$pos),
		end=c(IStable$pos)+1,
		LTR=IStable$LTR, 
		orientation=IStable$orientation, 
		readCount=IStable$readCount, 
		Upstream_eid=IStable$Upstream_eid,
		Upstream_gene=IStable$Upstream_gene,
		distance_up= IStable$distance_up,
		Downstream_eid=IStable$Downstream_eid,
		Downstream_gene=IStable$Downstream_gene,
		distance_down=IStable$distance_down,
		Genomic_Context=IStable$Genomic_Context,
		context=IStable$context),
		keep.extra.columns=T, ignore.strand=T)
		
		IS.3LTR <- IS[which(IS$LTR=="3LTR")]
		IS.5LTR <- IS[which(IS$LTR=="5LTR")]
	
		# Add to the LTR5 position the window up-down in order to subset on the range of LTR5
		start(IS.5LTR) <- start(IS.5LTR)-abs(window.up)
		end(IS.5LTR) <- end(IS.5LTR)+abs(window.down)
	
		print(" ## 2. Find overlaps between GRanges LTR3/LTR5 ## ")
			
		if(length(IS) > 1 & length(IS.3LTR) >= 1){
	
			pb <- txtProgressBar(0, max=length(IS.3LTR), style=3)
	
			for(k in 1:length(IS.3LTR)){			

				# Core function: find if the LTR3 is overlapping some LTR5 and confirm that they have the same orientation 
				LTR5 <- subsetByOverlaps(IS.5LTR, IS.3LTR[k])[IS.3LTR[k]$orientation == subsetByOverlaps(IS.3LTR[k], IS.5LTR)$orientation]
		
				## No LTR5 found? 
				# Keep the LTR3 position as it is.
				if(length(LTR5)==0){
					viral.pos <- c(viral.pos, IS.3LTR[k])
		
			
				## Found one LTR3 & one corresponding LTR5? 
				# Output the readCount of the two positions which have max(readCountLTR5, readCountLTR3)
				# The position is not modified ==> ALWAYS LTR3
				} else if(length(LTR5)==1) {
					
					IS.3LTR[k]$readCount <- IS.3LTR[k]$readCount <- max(LTR5$readCount, IS.3LTR[k]$readCount)
					viral.pos <- c(viral.pos, IS.3LTR[k])

					levels(viral.pos$LTR) <- c(levels(viral.pos$LTR), "LTR3LTR5")
					viral.pos[k]$LTR <- as.factor("LTR3LTR5")

			
				## Found one LTR3 and more than one LTR5 around?
				# Output the overlapping LTR with maximum readCount
				# Create an error file reporting such bizarre cases
				} else if(length(LTR5) >= 2) {
					print(paste("### Check potential confounding error in file: line:", k, "###", sep=" "))
			
						IS.3LTR[k]$readCount <- IS.3LTR[k]$readCount <- max(LTR5$readCount, IS.3LTR[k]$readCount)
						viral.pos <- c(viral.pos, IS.3LTR[k])

 						levels(viral.pos$LTR) <- c(levels(viral.pos$LTR), "LTR3LTR5")
						viral.pos[k]$LTR <- as.factor("LTR3LTR5")

					error <- c(error, append(IS.3LTR[k], LTR5)) 
				}
		
				setTxtProgressBar(pb, k)
			}

			## Add the information of the 5LTRs that are alone
			viral.pos <- c(viral.pos, IS.5LTR[which(!overlapsAny(IS.5LTR, viral.pos))])

			print(" ## 3. Re-Arrange the dataframe ## ")
	
			## Transform the GRange object into a dataframe, remove extra columns and calculate mean position

			viral.pos <- as.data.frame(viral.pos)[,-c(3:5)]
			#viral.pos$start <- round(apply(as.data.frame(viral.pos)[,2:3], 1, mean))
			#viral.pos <- viral.pos[,-3]
			colnames(viral.pos) <- c("chr", "Approx.pos", "LTR", "Orientation", "Clones", "Upstream_ID", "Upstream_Gene", "DistanceUp", "Downstream_ID", "Downstream_Gene", "DistanceDown", "Genomic_Context", "Context")
	
		# If only 0/1 IS detected or Only LTR5
		} else if(length(IS) < 1 | length(IS.3LTR) < 1){
		
			print(" Only LTR5 or LTR3 detected ")
			viral.pos <- as.data.frame(IS)[,c(1,2,4,7,3,8,9,10,11,12,13,14,15)]
			colnames(viral.pos) <- c("chr", "Approx.pos", "LTR", "Orientation", "Clones", "Upstream_ID", "Upstream_Gene", "DistanceUp",
			"Downstream_ID", "Downstream_Gene", "DistanceDown", "Genomic_Context", "Context")
			} else if(length(IS)==0){
				viral.pos <- data.frame(chr=as.numeric(0), Approx.pos=as.numeric(0), LTR=as.numeric(0), Orientation=as.numeric(0), 
				Clones=as.numeric(0), Upstream_ID=as.numeric(0), Upstream_Gene=as.numeric(0), DistanceUp=as.numeric(0), 
				Downstream_ID=as.numeric(0), Downstream_Gene=as.numeric(0), DistanceDown=as.numeric(0), Genomic_Context=as.numeric(0), 
				Context=as.numeric(0))
			}
		
	# If the input file is empty
	} else if(nrow(IStable)==0){
		viral.pos <- data.frame(chr=as.numeric(0), Approx.pos=as.numeric(0), LTR=as.numeric(0), Orientation=as.numeric(0), 
		Clones=as.numeric(0), Upstream_ID=as.numeric(0), Upstream_Gene=as.numeric(0), DistanceUp=as.numeric(0), 
		Downstream_ID=as.numeric(0), Downstream_Gene=as.numeric(0), DistanceDown=as.numeric(0), Genomic_Context=as.numeric(0), 
		Context=as.numeric(0))
	}
	
	print(" ## 4. Compute statistics ## ")

	## Compute various statistics and output them
	LTR3LTR5.perc <- 0
	LTR3LTR5readAbondance.perc  <- 0
	NumberIS <- 0
	LTR3alone.perc <- 0
	LTR3readAbondance.perc  <- 0
	LTR5alone.perc <- 0
	LTR5readAbondance.perc  <- 0
			
	LTR3LTR5.perc <- nrow(viral.pos[which(viral.pos$LTR=="LTR3LTR5"),])/nrow(viral.pos)
	LTR3LTR5readAbondance.perc  <- sum(viral.pos$Clones[which(viral.pos$LTR=="LTR3LTR5")])/sum(viral.pos$Clones)
	NumberIS <- nrow(viral.pos)
	LTR3alone.perc <- nrow(viral.pos[which(viral.pos$LTR=="3LTR"),])/nrow(viral.pos)
	LTR3readAbondance.perc  <- sum(viral.pos$Clones[which(viral.pos$LTR=="3LTR")])/sum(viral.pos$Clones)
	LTR5alone.perc <- nrow(viral.pos[which(viral.pos$LTR=="5LTR"),])/nrow(viral.pos)
	LTR5readAbondance.perc  <- sum(viral.pos$Clones[which(viral.pos$LTR=="5LTR")])/sum(viral.pos$Clones)
	
	stat <- data.frame(NumberIS_Unique=NumberIS, LTR3alone.perc=LTR3alone.perc, LTR5alone.perc=LTR5alone.perc, LTR3LTR5.perc=LTR3LTR5.perc,
	LTR3readAbondance.perc=LTR3readAbondance.perc, LTR5readAbondance.perc=LTR5readAbondance.perc, LTR3LTR5readAbondance.perc=LTR3LTR5readAbondance.perc)
	
	print(" ## 5. Output the results ## ")
			
	# Reorder based on number of clones
	viral.pos <- viral.pos[order(viral.pos$Clones, decreasing=T),]
			
	dir.create(paste(output,"/", "Stats", sep=""), showWarnings = FALSE)
	dir.create(paste(output,"/", "Errors", sep=""), showWarnings = FALSE)
	dir.create(paste(output,"/", "MergedIS", sep=""), showWarnings = FALSE)

	write.table(stat, paste(output, "/", "Stats/", name, "_Abudance.stats.txt", sep=""), sep="\t", row.names=F, quote=F)
	write.table(viral.pos, paste(output, "/", "MergedIS/", name, "_Insertions_U_orientation.MergedPositions.txt", sep=""), sep="\t", row.names=F, quote=F)
			
	if(length(error)>0){
		write.table(error, paste(output, "/", "Errors/", name, "_SpuriousPositions.txt", sep=""), sep="\t", row.names=F)
	}	
}


#files <- list.files("./tsv", full.names=T)
#a <- sapply(strsplit(list.files("./tsv"), split="="), function(x) paste(x[1:2], collapse="="))
#b <- sapply(strsplit(sapply(strsplit(list.files("./tsv"), split="="), function(x) x[3]), "_"), function(x) x[1])
#name <- paste(a,b, sep="=")

#for(i in 1:length(files)){
#	mergeISbyLTR(files[[i]], name[i], 10,10, "./")
#	}