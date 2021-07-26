#Gene Analysis Script
#Joshua Anderson
#Data Mining
#HW2


#Get the working directory from the user.
path <- readline("Please enter the filepath for the location where the program and data files have been extracted to.")
path <- gsub("\\\\", "/",path)
setwd(path)

#Import gene relation data set.  Assumes data files are in the working directory.
Genes_relation <- read.csv("Genes_relation.data")
Genes_relation_test <- read.csv("Genes_relation.test")
Keys <- read.csv("Keys.txt")

#Method to return number of tuples with no missing data
fullRows <- function(frame)
{
count <- 0
dims <- dim(frame)
dim1 <- dims[1]
dim2 <- dims[2]
for(i in 1:dim1)
{
  full <- TRUE
  for(j in 1:dim2)
  {
    if(frame[i,j] == "?")
    {full <- FALSE}
  }
  
  if(full == TRUE)
  {count <- count + 1}
}

return(count)
}


#Method to make a subset dataframe with only tuples with no missing data
fullFrame <- function(inFrame)
{
outFrame <- data.frame()
inDims <- dim(inFrame)
rowCount <- inDims[1]
colCount <- inDims[2]
cnames <- names(inFrame)

for(i in 1:rowCount)
{
  full <- TRUE
  tempDF <- inFrame[i,]
  for(j in 1:colCount)
  {
    if(inFrame[i,j] == "?")
      {full <- FALSE}
  }
  
  if(full == TRUE)
    {outFrame <- rbind(outFrame, tempDF)}
  
  
}
rownames(outFrame) <- 1:nrow(outFrame)
colnames(outFrame) <- cnames
return(outFrame)
}


#Method to calc the simple probability of a single attribute having the given value within a dataframe.
#   Returns a vector of (numerator, denominator, prob) that will be used to determine if
#   Laplacian correciton is needed during conditional prob calcs.
simpleProb <- function(inFrame, attribute, value)
{
  dims <- dim(inFrame)
  rowCount <- dims[1]
  cnames <- names(inFrame)
  attNum <- cIndex(cnames, attribute)
  count <- 0
  total <- 0
  
  
  for(i in 1:rowCount)
  {
    if(inFrame[i, attNum] == value)
    {
      count <- count + 1
    }
    
    total <- total + 1
  }
  res <- c(count, total, count/total)
  return(res)
}


#Method to calc the conditional probability of an attribute having a given value assuming a
#   conditional attribute has another given value. i.e. P(att = val | cAtt = cVal)
#   Returns a vector of (numerator, denominator, prob) that will be used to determine if
#   Laplacian correciton is needed.
condProb <- function(inFrame, cAtt, cVal, att, val)
{
  condFrame <- subFrame(inFrame, cAtt, cVal)
  return(simpleProb(condFrame, att, val))
}


#Method to calculate the conditional probability for a sample vector given some conditional attribute value.
#   i.e. X=(x1, x2, ..., xn), return the P(X=X | cAtt = cVal)
#   Laplaction correction applied as needed. Assumes conditional attribute independence naively.
vecCondProb <- function(inFrame, Xatts, Xvals, cAtt, cVal)
{
  IPprobs <- c()
  finProbs <- c()
  
  
  for(i in 1:length(Xvals))
  {
    IPprobs <- append(IPprobs, condProb(inFrame, cAtt, cVal, Xatts[i], Xvals[i]), after=length(IPprobs))
  }
  #print(IPprobs)
  
  #check for and apply Laplacian correction if needed.
  noZero <- TRUE
  vec <- seq(3, length(IPprobs), by=3)
  for(i in vec)
  {
    if(IPprobs[i] == 0)
    {noZero <- FALSE}
  }
  
  if(noZero == FALSE)
  {
    for(i in vec)
    {
      IPprobs[i] <- (IPprobs[i-2] + 1) / (IPprobs[i-1] + 15) #+15 since 15 classes for Localization
    }
  }
  
  #Calculate the final prob assuming conditional attribute independence
  finalProb <- 1
  for(i in vec)
  {finalProb <- finalProb * IPprobs[i]}
  
  
  return(finalProb)
}


#Method to determine a sample vectors class attribute using a naive Bayesian approach and a given full reference frame of data.
classify <- function(inFrame, sample, classAtt)
{
  cnames <- names(inFrame)
  classAttNum <- cIndex(cnames, classAtt)
  cleansample <- cleanSample(sample, classAttNum)
  cleanIndices <- filledIndices(sample, classAttNum)
  subframe <- selectSubFrame(inFrame, cleanIndices)
  subframe <- cbind(subframe, inFrame[classAttNum])
  subcNames <- names(subframe)
  classValList <- attList(inFrame, classAtt)
  classProbs <- c()
  condProbs <- c()
  
  #Gather simple probs for class attribute values based on subframe
  for(i in 1:length(classValList))
  {classProbs <- append(classProbs, simpleProb(subframe, classAtt, classValList[i])[3], after=length(classProbs))}

  
  #Gather associated conditional probs for cleaned sample and subframe
  for(i in 1:length(classValList))
  {condProbs <- append(condProbs, vecCondProb(subframe, subcNames[1:length(subcNames)-1], cleansample, classAtt, classValList[i]), after=length(condProbs))}

  
  #Multiply above vector elements for the final values to check
  checkVals <- c()
  for(i in 1:length(classProbs))
  {checkVals <- append(checkVals, classProbs[i] * condProbs[i], after=length(checkVals))}

  
  #Pull the index of the max value in checkVals to determine classificaiton. Return classValList[max]
  max <- which.max(checkVals)
  
  return(classValList[max])
  #testList <- list(classProbs, condProbs, checkVals, classValList, which.max(checkVals))
  #return(testList)
}


#Method to return a subframe with a selection of columns based on an input vector of indices.
selectSubFrame <- function(inFrame, cols)
{
  outFrame <- as.data.frame(matrix(nrow=dim(inFrame)[1])) #ensures we get matching row numbers when calling cbind.  Must pop off first column before returning.
  cnames <- names(inFrame)
  
  for(i in 1:dim(inFrame)[2])
  {
    check <- FALSE
    for(j in 1:length(cols))
    {
      if(i == cols[j])
      {check <- TRUE}
    }
    if(check == TRUE)
    {outFrame <- cbind(outFrame, inFrame[,i])}
  }
  outFrame <- outFrame[2:dim(outFrame)[2]]
  colnames(outFrame) <- cnames[cols]
  rownames(outFrame) <- 1:nrow(outFrame)

  return(outFrame)
}


#Method to return a subframe including only a certain attribute value.  Helper method for condProb
subFrame <- function(inFrame, att, val)
{
  dims <- dim(inFrame)
  outFrame <- data.frame()
  cnames <- names(inFrame)
  attNum <- cIndex(cnames, att)
  

  for(i in 1:dims[1])
  {
    tempDf <- inFrame[i,]
    if(inFrame[i, attNum] == val)
    {
      outFrame <- rbind(outFrame, tempDf)
    }
  }
  rownames(outFrame) <- 1:nrow(outFrame)
  colnames(outFrame) <- cnames
  
  return(outFrame)
}


#Method to return vector of possible values for an attribute.  Excludes "?"
attList <- function(inFrame, att)
{
  vals <- c()
  cnames <- names(inFrame)
  attNum <- cIndex(cnames, att)
  
  
  for(i in 1:dim(inFrame)[1])
  {
    if(inFrame[i,attNum] != "?")
    {
      if(length(vals) == 0)
      {vals <- c(inFrame[i,attNum])}
    
      else
      {
        contains <- FALSE
        for(j in 1:length(vals))
        {
          if(inFrame[i,attNum] == vals[j])
          {contains <- TRUE}
        }
        if(contains != TRUE)
        {vals <- append(vals, inFrame[i,attNum], after=length(vals))}
      }
    }  
  }
  
  return(vals)
}


#Small helper method to return the column index associated with an attribute name
cIndex <- function(cnames, att)
{
  for(i in 1:length(cnames))
  {
    if(att == cnames[i])
    {
      attNum <- i
    }
  }
  return(attNum)
}


#Helper method to clean a sample vector of "?".  Excludes the class attribute as well.
cleanSample <- function(sample, classAttNum)
{
  cleanSample <- c()
  
  for(i in 1:length(sample))
  {
    if(sample[i] != "?" & i != classAttNum)
    {cleanSample <- append(cleanSample, sample[i], after=length(cleanSample))}
  }
  return(cleanSample)
}


#Helper method to return a vector of indices in a sample vector that are non empty or "?" and not the class attribute
filledIndices <- function(sample, classAttNum)
{
  indices <- c()
  for(i in 1:length(sample))
  {
    if(sample[i] != "?" & i != classAttNum)
    {indices <- append(indices, i, after=length(indices))}
  }
  return(indices)
}


#Method to fill in missing values in originally given reference dataframe based on a full subframe using naive Bayesian approach
fillInFrameBayes <- function(FrameToFill, fullSubFrame)
{
  outFrame <- data.frame()
  dims <- dim(FrameToFill)
  cnames <- names(FrameToFill)
  
  for(i in 1:dims[1]) #Process and fill each row
  {
    inSample <- FrameToFill[i,]
    classAttList <- c()
    classAttIndices <- c()
    full <- TRUE
    
    for(j in 1:length(inSample)) #build list of col names and indices with empty values - these will be the 'class attributes' to be predicted for the row.
    {
      if(inSample[j] == "?")
      {
        classAttList <- append(classAttList, cnames[j], after=length(classAttList))
        classAttIndices <- append(classAttIndices, j, after=length(classAttIndices))
        full <- FALSE
      }
    }
    
    if(full == FALSE) #check if the row was already full, if so skip this step.  Else, fill in missing values with predictions
    {
      tempClassList <- c() #This will be used to hold predicted results so they can be added to the row all at once at the end.  This prevents predictions being used to make predictions.
      for(j in 1:length(classAttList)) #Predict value of each missing attribute and add it to list of predictions
      {
        tempClassList <- append(tempClassList, classify(fullSubFrame, inSample, classAttList[j]), after=length(tempClassList))
      }
      
      for(j in 1:length(tempClassList)) #Replace missing values in inSample
      {
        inSample[classAttIndices[j]] <- tempClassList[j]
      }
    }
    outFrame <- rbind(outFrame, inSample) #add filled row to output frame.  Output frame to keep things separate from original frame just in case.
    
  }# End of row processing and attachment to outFrame
  colnames(outFrame) <- cnames
  rownames(outFrame) <- 1:nrow(outFrame)
  
  return(outFrame)
}


#Method to predict values in test data based on filled in gene_relation data.  Slightly modified version of the previous method.
finalPredictionsBayes <- function(FrameToFill, fullSubFrame)
{
  outFrame <- data.frame(FrameToFill[,1:2]) #Copy all the gene IDs to predict.  Values in second column will be replaced by predictions
  dims <- dim(FrameToFill)
  classAtt <- "Localization"
  
  for(i in 1:dims[1]) #Process and predict localization for each row
  {
    inSample <- FrameToFill[i,]

    outFrame[i,2] <- classify(fullSubFrame, inSample, classAtt) #add result to output frame.

  }# End of row processing and attachment to outFrame
  colnames(outFrame) <- c("Gene ID", "Localization")
  
  return(outFrame)
}


#Method to remove "Function" column from input data sets.  Then checks for duplicate tuples and removes as needed.
removeFunction <- function(inFrame)
{
  dims <- dim(inFrame)
  outFrame <- as.data.frame(inFrame[,1]) #attach all geneIDs to output frame
  cnames <- names(inFrame)
  attNum <- cIndex(cnames, "Function")

  
  for(i in 2:dims[2]) #attach all cols except function col to outFrame
  {
    if(i != attNum)
    {outFrame <- cbind(outFrame, inFrame[,i])}
  }
  
  #Now process outFrame to remove duplicates.
  outFrame <- unique(outFrame)
  colnames(outFrame) <- c("GeneID", "Essential", "Class", "Complex", "Phenotype", "Motif", "Chromosome", "Localization")
  rownames(outFrame) <- 1:nrow(outFrame)
  
  return(outFrame)
}


#Method to return mode of an attribute, excluding "?"
getMode <- function(inFrame, att)
{
  cnames <- names(inFrame)
  dims <- dim(inFrame)
  classAttNum <- cIndex(cnames, att)
  listAtts <- attList(inFrame, att)
  counts <- c()
  
  for(i in 1:length(listAtts))
  {
    curCount <- 0
    for(j in 1:dims[1])
    {
      if(inFrame[j,classAttNum] == listAtts[i])
      {curCount <- curCount + 1}
    }
    counts <- append(counts, curCount, after=length(counts))
  }
  result <- listAtts[which.max(counts)]
  return(result)
}


#Method to return a vector of the modes for each attribute except localization. Function will have been removed already
getModes <- function(inFrame)
{
  dims <- dim(inFrame)
  cnames <- names(inFrame)
  modes <- c()
  
  for(i in 1:dims[2])
  {
    if(cnames[i] != "Localization")
    {modes <- append(modes, getMode(inFrame, cnames[i]), after=length(modes))}
  }
  return(modes)
}


#Method to fill a frame based on attribute modes.  Excludes localization
fillInFrameModes <- function(frameToFill, modes)
{
  cnames <- names(frameToFill)
  dims <- dim(frameToFill)
  outFrame <- data.frame()
  
  for(i in 1:dims[1]) #Fill missing values for each row
  {
    inRow <- frameToFill[i,1:dims[2]-1] #Exclude last column, Localization as this is the to be predicted att
    for(j in 1:length(inRow)) #Fill in each missing value based on mode
    {
      if(inRow[j] == "?")
      {inRow[j] <- modes[j]}
    }
    outFrame <- rbind(outFrame, inRow) #add filled rows to outframe
  }
  outFrame <- cbind(outFrame, frameToFill[,dims[2]]) #Tack back on the localization column
  colnames(outFrame) <- cnames
  rownames(outFrame) <- 1:nrow(outFrame)
  
  return(outFrame)  
}


#Method to find accuracy of classification process
accuracy <- function(keys, results)
{
  correct <- 0
  count <- dim(results)[1]
  
  for(i in 1:dim(results)[1]) #For every result row
  {
    check <- results[i,]
    
    for(j in 1:dim(keys)[1]) #Locate row in keys with matching geneID and check for match
      if(check[1] == keys[j,1])
      {
        if(check[2] == keys[j,2]) #add to count if classification correct
        {
          correct <- correct + 1
        }
      }
  }
  return(correct/count)
}





#Main method to test approach replacing missing vals with mode
main <- function(GeneRelation, GeneRelationTest)
{
  trimmedGeneRelation <- removeFunction(GeneRelation)
  modes <- getModes(trimmedGeneRelation)
  trimmedGeneRelationTest <- removeFunction(GeneRelationTest)
  filledTrimmedGeneRelation <- fillInFrameModes(trimmedGeneRelation, modes)
  filledTrimmedGeneRelationTest <- fillInFrameModes(trimmedGeneRelationTest, modes)
  results <- finalPredictionsBayes(filledTrimmedGeneRelationTest, filledTrimmedGeneRelation)
  Accuracy <- accuracy(Keys, results)
  
  View(Keys)
  View(results)
  View(Accuracy)
  resList <- list(modes, GeneRelation, GeneRelationTest, trimmedGeneRelation, trimmedGeneRelationTest, filledTrimmedGeneRelation, filledTrimmedGeneRelationTest, results, Keys, Accuracy)
  
  return(resList)
}


Results <- main(Genes_relation, Genes_relation_test)