# A. G. Sangster

# This program takes the downloaded MAF and reduces the data to relevant features and cases with 
# confirmed stage and grade in clinical data. Then subsets the data according to the number of 
# calls a mutation has received and implements a FDR. 

# This program requires 'DDR_CM_geneList.csv' in the same directory as this program

# There is the option to generate figures and save the resulting subsetted data.




#------------------------------------------------------------------------------------------------ Set Preferences

# If this program was not run directly after downloading the data, the data will need to be 
# loaded to the environment. 
# Note: to load the data using this preference, it must be saved by TCGA-BLCA-download.R
load_data = FALSE

# show the number of cases, genes and mutations as the data is reduced
show_data_reduction = TRUE

# show the venn diagram figures that compare the mutation variants found by the different MAF
# and show the co-efficient of variation histogram that compares the MAF
show_MAF_figures = TRUE

# save the least, mid and most conservative data subsets in their own directory
save_data_subsets = TRUE

# save the supplementary tables of the data 
# supplementary table 1 has all mutations, all cases and shows which pipelines found them
# supplementary table 2 has a summary table for each data subset
save_supp = TRUE


#------------------------------------------------------------------------------------------------ Load Packages and Data
library(VennDiagram)


if (load_data){
  thisDir <- getwd()
  setwd(paste(thisDir,"TCGA-Data",sep="/"))
  mutect2 <- read.csv("TCGA-BLCA-mutect2.csv", stringsAsFactors = F)
  varscan2 <- read.csv("TCGA-BLCA-varscan2.csv", stringsAsFactors = F)
  somaticsniper <- read.csv("TCGA-BLCA-somaticsniper.csv", stringsAsFactors = F)
  muse <- read.csv("TCGA-BLCA-muse.csv", stringsAsFactors = F)
  setwd(thisDir)
}

#------------------------------------------------------------------------------------------------ Helper Functions

# reduce the number of columns in the data
# keep hugo gene name (1), patient-sample identifier(16), 
# mutation variant(35) and mutation impact(94)
# add a shortened bar code to the identify the patient
column_reduction <- function(dat){
  
  #shorten each barcode to identify each patient and not the entire sample
  sBar <- c()
  for (ii in 1:(dim(dat)[1])){
    bar <- dat[ii,16]
    sBar <- c(sBar, substr(bar,1,12))
  }
  
  out <- data.frame("Short_Bars" = sBar, dat[,c(1,16,35,94)])
  return(out)
}



# remove cases: TCGA-FJ-A3Z9, TCGA-ZF-AA4N, TCGA-CF-A9FM, TCGA-BL-A0C8, TCGA-HQ-A2OE
# using shortened bar codes, this was done because of missing values from the clinical data
clinical_reduction <- function(data){
  
  cases <- c("TCGA-FJ-A3Z9", "TCGA-ZF-AA4N", "TCGA-CF-A9FM", "TCGA-BL-A0C8", "TCGA-HQ-A2OE")
  out <- data[!(data$Short_Bar %in% cases),]
}



# Use DDR_CM_geneList.csv to identify genes involved in DDR and/or CM processed
# remove genes that are not included in this list of DDR and CM genes
DDR_CM_reduction <- function(dat){
  
  DDR_CM <- read.csv("DDR_CM_geneList.csv")
  #remove citations and class of genes
  DDR_CM <- DDR_CM$Gene
  
  #Select genes that are involved in DDR and/or CM processes
  out <- dat[dat$Hugo_Symbol %in% DDR_CM,]
  
  return(out)
}



# False Discovery Rate (FDR)
# identify genes that are mutated in more than 7% of cases
FDR <- function(dat){
  
  # make a list of genes and bar codes
  mut <- cbind(dat$Short_Bars, dat$Hugo_Symbol)
  
  # remove multiple mutations in same gene in a case
  mut <- unique(mut)
  
  # find gene mutation frequency 
  freq <- table(mut[,2])/(length(unique(dat[,1])))
  # find which genes are mutated in more than 7% of cases
  freq <- freq[freq>0.065]
  
  # return those gene names
  out <- row.names(freq)
  
  return(out)
}




# make supplementary table 2
# given the data post-processing 
supp2 <- function(data, DDR_CM){
  
  #get gene list and initialize output
  genes <- unique(data$Hugo_Symbol)
  out <- data.frame()
  
  for (ii in 1:length(genes)){
    #set gene
    gene <- genes[ii]
    #get gene class
    class <- DDR_CM[DDR_CM$Gene==gene,]$Class
    #determine number of mutations
    n <- sum(data$Hugo_Symbol==gene)
    #determine number of cases with the mutation
    bars <- length(unique(data[data$Hugo_Symbol==gene,]$Short_Bars))
    #determine percentage of cases with the mutation
    per <- round(100*bars/length(unique(data$Short_Bars)),2)
    
    out <- rbind(out, cbind(gene, class, n, bars, per))
  }
  
  #order the table by percent as decreasing
  final <- out[order(as.numeric(out$per), decreasing =TRUE),]
  
  #change column names to more meaningful
  colnames(final) <- c("Gene", "DDR CM Class", "Number of Mutations", "Number of cases with Mutation", "% Cases with Mutation")
  
  return(final)
}







#------------------------------------------------------------------------------------------------ Data Reduction


# colum reduction
muse <- column_reduction(muse)
mutect2 <- column_reduction(mutect2)
varscan2 <- column_reduction(varscan2)
somaticsniper <- column_reduction(somaticsniper)


# clinical reduction
muse <- clinical_reduction(muse)
mutect2 <- clinical_reduction(mutect2)
varscan2 <- clinical_reduction(varscan2)
somaticsniper <- clinical_reduction(somaticsniper)


# make data subgroups (least conservative, mid conservative and most conservaitve)
# add the name of the mutation calling pipeline to the data entry
muse <- cbind(muse,"muse")
mutect2 <- cbind(mutect2, "mutect2")
varscan2 <- cbind(varscan2, "varscan2")
somaticsniper <- cbind(somaticsniper, "somaticsniper")

# merge the data into one data frame that includes all mutations 
muse_mutect2 <- merge(muse, mutect2, by=c("Short_Bars", "Hugo_Symbol", "HGVSc", "Tumor_Sample_Barcode", "IMPACT"), all = TRUE)
varscan2_somaticsniper <- merge(varscan2, somaticsniper, by=c("Short_Bars", "Hugo_Symbol", "HGVSc", "Tumor_Sample_Barcode", "IMPACT"), all=TRUE)
all_mutations <- merge(muse_mutect2, varscan2_somaticsniper, by=c("Short_Bars", "Hugo_Symbol", "HGVSc", "Tumor_Sample_Barcode", "IMPACT"), all=TRUE)

# determine how many mutation calls each row of all_mutations has
num_calls <- rowSums(!is.na(all_mutations[,6:9]))

# least conservative:  keep all mutations and remove the columns that indicate mutation calling pipeline
least <- all_mutations[,1:5]
# mid conservaitve:  only keep mutations found by at least 2 different mutation calling pipelines
mid <- all_mutations[num_calls>1,1:5]
# most conservative:  only keep mutations found by all mutation calling pipelines 
most <- all_mutations[num_calls==4,1:5]

#save the numbers for 'show_data_reduction'
least_nums <- c(length(unique(least$Short_Bars)),length(unique(least$Hugo_Symbol)),dim(least)[1])
mid_nums <- c(length(unique(mid$Short_Bars)),length(unique(mid$Hugo_Symbol)),dim(mid)[1])
most_nums <- c(length(unique(most$Short_Bars)),length(unique(most$Hugo_Symbol)),dim(most)[1])
subgroups <- data.frame(least_nums,mid_nums,most_nums)
colnames(subgroups) = c("Least Conservative", "Mid Conservative", "Most Conservative")
rownames(subgroups) = c("Cases","Genes", "Mutations")


# filter mutation impact
leastHM <- least[least$IMPACT=="HIGH" | least$IMPACT=="MODERATE",]
leastH <- least[least$IMPACT=="HIGH",]
midHM <- mid[mid$IMPACT=="HIGH" | mid$IMPACT=="MODERATE",]
mostHM <- most[most$IMPACT=="HIGH" | most$IMPACT=="MODERATE",]

#save the numbers for 'show_data_reduction'
leastHM_nums <- c(length(unique(leastHM$Short_Bars)),length(unique(leastHM$Hugo_Symbol)),dim(leastHM)[1])
leastH_nums <- c(length(unique(leastH$Short_Bars)),length(unique(leastH$Hugo_Symbol)),dim(leastH)[1])
midHM_nums <- c(length(unique(midHM$Short_Bars)),length(unique(midHM$Hugo_Symbol)),dim(midHM)[1])
mostHM_nums <- c(length(unique(mostHM$Short_Bars)),length(unique(mostHM$Hugo_Symbol)),dim(mostHM)[1])
impact <- data.frame(leastHM_nums,leastH_nums,midHM_nums,mostHM_nums)
colnames(impact) = c("Least-High and Moderate", "Least-High", "Mid-High and Moderate", "Most-High and Moderate")
rownames(impact) = c("Cases","Genes", "Mutations")


# filter DDR and CM genes
leastHM <- DDR_CM_reduction(leastHM)
leastH <- DDR_CM_reduction(leastH)
midHM <- DDR_CM_reduction(midHM)
mostHM <- DDR_CM_reduction(mostHM)

#save the numbers for 'show_data_reduction'
leastHM_nums <- c(length(unique(leastHM$Short_Bars)),length(unique(leastHM$Hugo_Symbol)),dim(leastHM)[1])
leastH_nums <- c(length(unique(leastH$Short_Bars)),length(unique(leastH$Hugo_Symbol)),dim(leastH)[1])
midHM_nums <- c(length(unique(midHM$Short_Bars)),length(unique(midHM$Hugo_Symbol)),dim(midHM)[1])
mostHM_nums <- c(length(unique(mostHM$Short_Bars)),length(unique(mostHM$Hugo_Symbol)),dim(mostHM)[1])
DDR_CM_nums <- data.frame(leastHM_nums,leastH_nums,midHM_nums,mostHM_nums)
colnames(DDR_CM_nums) = c("Least-High and Moderate", "Least-High", "Mid-High and Moderate", "Most-High and Moderate")
rownames(DDR_CM_nums) = c("Cases","Genes", "Mutations")



# FDR
leastHM <- leastHM[leastHM$Hugo_Symbol%in%FDR(leastHM),]
leastH <- leastH[leastH$Hugo_Symbol%in%FDR(leastH),]
midHM <- midHM[midHM$Hugo_Symbol%in%FDR(midHM),]
mostHM <- mostHM[mostHM$Hugo_Symbol%in%FDR(mostHM),]

#save the numbers for 'show_data_reduction'
leastHM_nums <- c(length(unique(leastHM$Short_Bars)),length(unique(leastHM$Hugo_Symbol)),dim(leastHM)[1])
leastH_nums <- c(length(unique(leastH$Short_Bars)),length(unique(leastH$Hugo_Symbol)),dim(leastH)[1])
midHM_nums <- c(length(unique(midHM$Short_Bars)),length(unique(midHM$Hugo_Symbol)),dim(midHM)[1])
mostHM_nums <- c(length(unique(mostHM$Short_Bars)),length(unique(mostHM$Hugo_Symbol)),dim(mostHM)[1])
FDR_nums <- data.frame(leastHM_nums,leastH_nums,midHM_nums,mostHM_nums)
colnames(FDR_nums) = c("Least-High and Moderate", "Least-High", "Mid-High and Moderate", "Most-High and Moderate")
rownames(FDR_nums) = c("Cases","Genes", "Mutations")



if (show_data_reduction){
  cat("\nAll mutations sorted into subgroups\n")
  print(subgroups)
  cat("\n\nSubgroups filtered for mutation impact\n")
  print(impact)
  cat("\n\nSubgroups filtered for mutation impact and DDR-CM genes\n")
  print(DDR_CM_nums)
  cat("\n\nSubgroups filtered (as above) with FDR\n")
  print(FDR_nums)
}


if (save_data_subsets){
  
  # make directory for data subsets
  dir.create("TCGA-Data-Subsets")
  thisDir <- getwd()
  setwd(paste(thisDir,"TCGA-Data-Subsets",sep="/"))
  
  write.csv(leastHM, "leastHM.csv", row.names = F)
  write.csv(leastH, "leastH.csv", row.names = F)
  write.csv(midHM, "midHM.csv", row.names = F)
  write.csv(mostHM, "mostHM.csv", row.names = F)
  
  # back to initial working directory
  setwd(thisDir)
}


if (save_supp){
  
  DDR_CM <- read.csv("DDR_CM_geneList.csv")
  
  # make directory for data subsets
  dir.create("TCGA-Supplementary")
  thisDir <- getwd()
  setwd(paste(thisDir,"TCGA-Supplementary",sep="/"))
  
  write.csv(all_mutations, "Supplementary2-collectively.csv", row.names = F)
  
  leastHM_sum <- supp2(leastHM, DDR_CM)
  leastH_sum <- supp2(leastH, DDR_CM)
  midHM_sum <- supp2(midHM, DDR_CM)
  mostHM_sum <- supp2(mostHM, DDR_CM)
  
  write.csv(leastHM_sum, "Supplementary1-1.csv", row.names = F)
  write.csv(leastH_sum, "Supplementary1-2.csv", row.names = F)
  write.csv(midHM_sum, "Supplementary1-3.csv", row.names = F)
  write.csv(mostHM_sum, "Supplementary1-4.csv", row.names = F)
  
  # back to initial working directory
  setwd(thisDir)
  
}





#------------------------------------------------------------------------------------------------ MAF figures

#finds the coefficient of variance with all NA values counted as 0
coeffOfVarianceNAas0 <- function(data){
  data <- as.numeric(data)
  
  # mean
  s <- 0
  n <- dim(array(data))
  for(i in 1:n){
    if(!is.na(data[i])){
      s <- s + data[i]
    }
  }
  m <- (s/n)
  
  #standard deviation
  s <- 0
  for(i in 1:n){
    if(!is.na(data[i])){
      s <- s + ((data[i] - m)**2)
    }
  }
  sd <- sqrt((s/(n-1)))
  
  #coeff of variation
  var <- sd/m
  return (var)
}



makeVenn <- function(bar){
  #paste the hugo symbol and mutation variant together to make unique mutation identifiers
  #do this for each mutation calling method 
  #m2m = mutec2 mutations, mm = muse mutations, ssm = somatic sniper mutations, v2m = varscan2 mutations
  m2m <- paste(mutect2[mutect2$Tumor_Sample_Barcode==bar,]$Hugo_Symbol, mutect2[mutect2$Tumor_Sample_Barcode==bar,]$HGVSc)
  mm <- paste(muse[muse$Tumor_Sample_Barcode==bar,]$Hugo_Symbol, muse[muse$Tumor_Sample_Barcode==bar,]$HGVSc)
  ssm <- paste(somaticsniper[somaticsniper$Tumor_Sample_Barcode==bar,]$Hugo_Symbol, somaticsniper[somaticsniper$Tumor_Sample_Barcode==bar,]$HGVSc)
  v2m <- paste(varscan2[varscan2$Tumor_Sample_Barcode==bar,]$Hugo_Symbol, varscan2[varscan2$Tumor_Sample_Barcode==bar,]$HGVSc)
  
  #make a sorted list of unique mutation identifiers (dupicates removed)
  ml <- sort(unique(c(m2m, mm, v2m, ssm)))
  
  #make an empty list with a spot for each unique mutation identifier
  Mutect2 <- rep(0,length(ml))
  Muse <- rep(0,length(ml))
  SomaticSniper <- rep(0,length(ml))
  Varscan2 <- rep(0,length(ml))
  
  #put those lists in a data frame
  mtable <- data.frame(Mutect2,Muse,SomaticSniper,Varscan2)
  rownames(mtable) <- ml
  
  #fill the data frame to generate a matrix that indicates which mutation variants are present 
  #in each mutation calling pipeline
  for(i in 1:length(ml)){
    gene <- ml[i]
    if(gene%in%m2m){
      mtable[i,1] <- 1
    }
    if(gene%in%mm){
      mtable[i,2] <- 1
    }
    if(gene%in%ssm){
      mtable[i,3] <- 1
    }
    if(gene%in%v2m){
      mtable[i,4] <- 1
    }
  }
  
  
  
  #a1 = mutect2
  a1 <- sum(mtable$Mutect2)
  #a2 = muse
  a2 <- sum(mtable$Muse)
  #a3 = somatic sniper
  a3 <- sum(mtable$SomaticSniper)
  #a4 = varscan2
  a4 <- sum(mtable$Varscan2)
  
  #intersect a1 and a2
  n12 <- sum(mtable$Mutect2==1 & mtable$Muse==1)
  #intersect a1 and a3
  n13 <- sum(mtable$Mutect2==1 & mtable$SomaticSniper==1)
  #intersect a1 and a4
  n14 <- sum(mtable$Mutect2==1 & mtable$Varscan2==1)
  
  #intersect a2 and a3
  n23 <- sum(mtable$Muse==1 & mtable$SomaticSniper==1)
  #intersect a2 and a4
  n24 <- sum(mtable$Muse==1 & mtable$Varscan2==1)
  #intersect a3 and a4
  n34 <- sum(mtable$SomaticSniper==1 & mtable$Varscan2==1)
  
  #intersect a1,a2,a3
  n123 <- sum(mtable$Mutect2==1 & mtable$Muse==1 & mtable$SomaticSniper==1)
  #intersect a1,a2,a4
  n124 <- sum(mtable$Mutect2==1 & mtable$Muse==1 & mtable$Varscan2==1)
  #intersect a1,a3,a4
  n134 <- sum(mtable$Mutect2==1 & mtable$SomaticSniper==1 & mtable$Varscan2==1)
  #intersect a2,a3,a4
  n234 <- sum(mtable$Muse==1 & mtable$SomaticSniper==1 & mtable$Varscan2==1)
  
  #interest a1,a2,a3,a4 (all)
  n1234 <- sum(mtable$Mutect2==1 & mtable$Muse==1 & mtable$SomaticSniper==1 & mtable$Varscan2==1)
  
  #make the venn diagram
  grid.newpage()
  draw.quad.venn(a1,a2,a3,a4,n12,n13,n14,n23,n24,n34,n123,n124,n134,n234,n1234,
                 category = c("Mutect2", "Muse", "SomaticSniper", "Varscan2"), 
                 fill = c("white", "white", "white", "white"),
                 cex = 1.5, cat.cex = 1.1)
  
  return(mtable)
}



if (show_MAF_figures){
  
  #generate a table that contains the number of mutation found by each method for each case
  mutectF <- data.frame(sort(table(mutect2$Tumor_Sample_Barcode), decreasing = T))
  museF <- data.frame(sort(table(muse$Tumor_Sample_Barcode), decreasing = T))
  somF <- data.frame(sort(table(somaticsniper$Tumor_Sample_Barcode), decreasing = T))
  varF <- data.frame(sort(table(varscan2$Tumor_Sample_Barcode), decreasing = T))
  tab <- merge(x = mutectF, y = museF, by = "Var1", all = TRUE)
  colnames(tab) <- c("Var1","mutect2","muse")
  tab <- merge(x = tab, y = somF, by = "Var1", all = TRUE)
  tab <- merge(x = tab, y = varF, by = "Var1", all = TRUE)
  colnames(tab) <- c("Var1","mutect2","muse","somaticSniper","varscan2")
  
  
  #generate a list of cov values
  varCol <- NULL
  for(i in 1:dim(tab)[1]){
    varCol <- rbind(varCol, coeffOfVarianceNAas0(tab[i,2:5]))
  }
  #make coefficient of variation histogram
  hist(varCol, breaks=30, main="Difference in Number of Mutation Calls \n Between Calling Pipelines per Sample", xlab="Coefficient of Variation",cex.main=2, cex.lab=1.5, cex.axis=1.5)
  
  
  #make venn diagrams
  #these are the cases selected for figures
  #the function (makeVenn) can be called with any valid barcode
  bars <- c("TCGA-XF-A9SL-01A-11D-A391-08", "TCGA-K4-A83P-01A-11D-A34U-08", "TCGA-XF-A9SI-01A-11D-A391-08", "TCGA-DK-A1AF-01A-11D-A13W-08", "TCGA-LT-A8JT-01A-11D-A364-08", "TCGA-DK-A6AW-01A-11D-A30E-08")
  
  for (i in 1:6){
    makeVenn(bars[i])
  }
  
}



