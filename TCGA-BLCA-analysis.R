# A. G. Sangster

# This program takes the prepared data subsets and performs a statistical analysis to determine 
# mutually exclusive mutation patterns between gene pairs. This program generates summary tables that 
# can be saved using 'saveTables' and also generates plots of the random sampling distributions, its mean 
# and the observed conditional probability that can be saved using the 'savePlots' in Preferences section.







#------------------------------------------------------------------------------------------------ Set Preferences

# If this program was not run directly after preparing the data, the data will need to be 
# loaded to the environment. 
# Note: to load the data using this preference, it must be saved by TCGA-BLCA-dataPrep.R
load_data = FALSE

# save summary tables for mutually exclusive relatinoships for each mutation subset
saveTables = TRUE

# make and save plots for relationships that show some mutual exclusivity
# shows RGCP distribution, its mean and the value observed in the data
savePlots = TRUE

# make and save plots for all gene pairings
saveAllPlots = FALSE


#------------------------------------------------------------------------------------------------ Load Packages and Data



if (load_data){

  thisDir <- getwd()
  setwd(paste(thisDir,"TCGA-Data-Subsets",sep="/"))
  
  leastH <- read.csv("leastH.csv", stringsAsFactors = F)
  leastHM <- read.csv("leastHM.csv", stringsAsFactors = F)
  midHM <- read.csv("midHM.csv", stringsAsFactors = F)
  mostHM <- read.csv("mostHM.csv", stringsAsFactors = F)
  
  setwd(thisDir)
  
}





#------------------------------------------------------------------------------------------------ Analysis Functions



#P(notA|B)
#notAgivenB(mutationData, geneOrder) >> conditional probability matrix
# This function takes the mutation data after preparations
# Finds the conditional probability of not finding a mutation in a particular gene
# given that there is a mutation in another gene for all gene pairs
notAgivenB <- function(data, geneOrder){
  
  allBars <- unique(data$Short_Bars)  
  nGenes <- dim(array(geneOrder))
  
  #generate a list for the 
  genesNot <- NULL
  for(a in 1:nGenes){
    genesNot <- c(genesNot, paste("not", geneOrder[a], sep = ""))
  }
  
  #initialize the conditional probability matrix
  m <- array(0, dim = c(nGenes, nGenes),
             dimnames = list(genesNot, geneOrder))
  
  #loop through each gene pairing and determine the conditional probability
  for(i in 1:nGenes){
    gene1 <- geneOrder[i]
    for(j in 1:nGenes){
      if(i!=j){
        gene2 <- geneOrder[j]
        bars1 <- unique(data[data$Hugo_Symbol==gene1,]$Short_Bars)
        bars2 <- unique(data[data$Hugo_Symbol==gene2,]$Short_Bars)
        bars <- intersect(bars1, setdiff(allBars, bars2))
        m[j,i] <- round((dim(array(bars))/(dim(array(bars1))))*100, 2)
      }
    }
  }
  
  return(m)
}





#makeRandom(mutationData) >> mutationData
#This function takes maf and returns a new maf with the mutations randomly distributed
#amongst the samples such that there are the same mutations and the same number of samples
#and the same number of mutations per sample.
makeRandom <- function(maf){
  bars <- as.character(maf$Short_Bars)
  muts <- as.character(maf$Hugo_Symbol)
  newBars <- sample(bars)
  out <- data.frame("Hugo_Symbols"=muts, "Short_Bars"=newBars)
  return(out)
}




#RGCPmatrices(mutationData, N) >> RGCP[,,N]
#This function takes the preprocessed mutation data (MAF)
#The MAF is randomized N times and a CP matrix for each randomized MAF is generated
#The random CP matrices are collected into a 3D data structure and returned
RGCPmatrices <- function(maf, N){
  
  nGenes <- length(unique(maf$Hugo_Symbol))
  #set gene order to be the names of the genes by decreasing frequency
  geneOrder <-  names(sort(table(maf$Hugo_Symbol), decreasing=T))
  
  #generate list of 'not'gene for labels
  genesNot <- NULL
  for(a in 1:nGenes){
    genesNot <- c(genesNot, paste("not", geneOrder[a], sep = ""))
  }
  
  #initalize the 3D data frame to contain N conditional probability matrices
  data <- array(0, dim = c(nGenes, nGenes, N),
                dimnames = list(genesNot,geneOrder))
  data[,,1] <- notAgivenB(maf, geneOrder)
  
  #make random matrices
  for(i in 1:N){
    data[,,i] <- notAgivenB(makeRandom(maf), geneOrder)
  }
  
  return(data)
  
}




# sigmaPaired(mutationData, RGCP[,,N]) >> table
# This function takes a data subset and its correpsonding RGCP matrix
# For both relationships in each gene pairing, the distribution in analyzed by 
# establishing normality (shapiro wilks) and determining the sigma value.
# If both relationships have normal distrbutions and sigma values over 1.5 
# they are gathered into a table and returned as the result
sigmaPaired <- function(dat,RGCP){
  
  # set relevant variables and initialize output
  n <- dim(RGCP)[1]
  geneOrder <-  names(sort(table(dat$Hugo_Symbol), decreasing=T))
  p <- notAgivenB(dat, geneOrder)
  out <- data.frame(matrix(ncol=5, nrow=0))
  
  # for each gene pairing that is not paired with itself
  for (ii in 1:n){
    for (jj in ii:n){
      if (ii!=jj){
        
        # find the distribution for that gene pairing P(notA|B)
        r1 <- RGCP[ii,jj,]
        # determine the sigma value
        sigma1 <- round(abs(mean(r1)-p[ii,jj])/sd(r1),2)
        
        # find the distribution for the gene pairing in the other direction P(notB|A)
        r2 <- RGCP[jj,ii,]
        # determine the sigma value
        sigma2 <- round(abs(mean(r2)-p[jj,ii])/sd(r2),2)
        
        
        #include in output if both sigma values are greater than 1.5:  (sigma1>1.5)&(sigma2>1.5)
        #and only if the relationship is a negative correlation:  (mean(r1)-p[ii,jj]<0)&(mean(r2)-p[jj,ii]<0)
        if ((sigma1>=1.5)&(sigma2>=1.5)&(mean(r1)-p[ii,jj]<0)&(mean(r2)-p[jj,ii]<0)){
          
          #check distribution with shapiro wilks
          st1 <- shapiro.test(RGCP[ii,jj,])
          st2 <- shapiro.test(RGCP[jj,ii,])
          p1 <- st1$p.value
          p2 <- st2$p.value
          
          # if both distrbutions are normal
          if (p1>0.05&p2>0.05){
            
            #include category (double line, single line, or dotted line)
            if ((sigma1>2.75)&(sigma2>2.75)){
              out <- rbind(out, cbind(colnames(RGCP)[ii],colnames(RGCP)[jj],sigma1,sigma2,"double line"))
              
            }else if((sigma1>2)&(sigma2>2)){
              out <- rbind(out, cbind(colnames(RGCP)[ii],colnames(RGCP)[jj],sigma1,sigma2,"single line"))
              
            }else{
              out <- rbind(out, cbind(colnames(RGCP)[ii],colnames(RGCP)[jj],sigma1,sigma2,"dotted line"))
              
            }
          }
        }
      }
    }
  }
  colnames(out) <- c("Gene A", "Gene B", "Sigma for P(notA|B)" , "Sigma for P(notB|A)", "Shown As")
  
  return(out)
}







#------------------------------------------------------------------------------------------------ Figure Functions

# takes mutuation data (dat), the RGCP for that data (RGCP), a name for the PDF that will contain the plots
# and a bool that determine whether all plots will be generated or just the ones that show exclusivity
# This function generates plots that show the RGCP distribution, its mean and the observed conditional probability
# for a given gene pairing. These plots are saved in a PDF in a directory called Results that can be found in the
# current working directory.
makeDistPlots <- function(dat, RGCP, name, all_plots){
  
  # make location for plots to go
  thisDir <- getwd()
  path <- paste(thisDir, "/Results/", name, ".PDF", sep="")
  pdf(file=path) 
  
  # set variables for the loop 
  n <- dim(RGCP)[1] #number of random samples
  p <- notAgivenB(dat, colnames(RGCP)) #observed conditional probability
  
  # look at each mutuation pairing to determine which should be plotted
  for (ii in 1:n){
    for (jj in 1:n){
      if (ii!=jj){
        r <- RGCP[ii,jj,]
        
        sigma <- round(abs(mean(r)-p[ii,jj])/sd(r),2)
        
        #check distribution with shapiro wilks
        st1 <- shapiro.test(RGCP[ii,jj,])
        st2 <- shapiro.test(RGCP[jj,ii,])
        p1 <- st1$p.value
        p2 <- st2$p.value
        
        # if both distrbutions are normal
        if (p1>0.05&p2>0.05){
        
          # if the pairing shows significant exclusivity or if should show all plots
          # significant exclusivity > negative correlation and sigma>1.5
          if (((mean(r)<p[ii,jj])&(sigma>=1.5)) | all_plots) {
              title <- paste(rownames(p)[ii],'-',colnames(p)[jj],'\n',sigma,' sigma',sep="")
              
              minr = min(r)
              maxr = max(r)
              hist(r, xlim=c(min(minr-10,p[ii,jj]),100), col='white', main=title, xlab="Probability")
              #breaks=?
              abline(v=mean(r), col="blue", lwd=2)
              abline(v=p[ii,jj], col="red", lwd=2)
          }
        }
      }
    }
  }
  
  dev.off()
}









#------------------------------------------------------------------------------------------------ Analysis Calls

set.seed(0)
N = 2000

# generate RGCP for N random samples for each data subset
#RGCP_LHM <- RGCPmatrices(leastHM, N)
#RGCP_MiHM <- RGCPmatrices(midHM, N)
#RGCP_MoHM <- RGCPmatrices(mostHM, N)
#RGCP_LH <- RGCPmatrices(leastH, N)

# generate summary tables for mutually exclusive relationships
stLHM <- sigmaPaired(leastHM, RGCP_LHM)
stMiHM <- sigmaPaired(midHM, RGCP_MiHM)
stMoHM <- sigmaPaired(mostHM, RGCP_MoHM)
stLH <- sigmaPaired(leastH, RGCP_LH)


if (saveTables | savePlots) {
  # make directory for results
  dir.create("Results")
}


if (saveTables) {
  
  thisDir <- getwd()
  setwd(paste(thisDir,"Results",sep="/"))
  
  write.csv(stLHM, "leastHM-results.csv", row.names = FALSE)
  write.csv(stMiHM, "midHM-results.csv", row.names = FALSE)
  write.csv(stMoHM, "mostHM-results.csv", row.names = FALSE)
  write.csv(stLH, "leastH-results.csv", row.names = FALSE)
  
  # back to initial working directory
  setwd(thisDir)
  
}



if (savePlots){
  
  makeDistPlots(leastHM, RGCP_LHM, "LeastHM", saveAllPlots)
  makeDistPlots(leastH, RGCP_LH, "LeastH", saveAllPlots)
  makeDistPlots(midHM, RGCP_MiHM, "MidHM", saveAllPlots)
  makeDistPlots(mostHM, RGCP_MoHM, "MostHM", saveAllPlots)
  
}



