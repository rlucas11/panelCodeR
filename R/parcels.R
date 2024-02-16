################################
#####CREATE PARCEL FUNCTION#####
################################

#Need to structure data by Wave, only including items to parcel
#Relies on 'gtools' package for evaluation of odd and even lines
#pattern should be in the form of c(n1,n2,n3), where n1, n2, and n3 are the number of items per parcel


parcel <- function(data,items,waves,pattern) {
require(gtools)

#Create Temporary Matrix
loadingSum <- matrix(0,items,1)

#Get Average Loading
for (i in 1:waves) {
  start=1+(i-1)*items
  finish=i*items
  loadingSum <- loadingSum+factanal(na.omit(data[,start:finish]), 1)$loadings[,1]
}
loadingAvg <- loadingSum/waves

#Rank Items by Loading
itemList <- data.frame(cbind(loadingAvg,matrix(1:items,items,1)))
itemList <- itemList[order(itemList$X1),]

#Create Pattern Matrix for Creating Parcels
temp <- matrix(nrow=max(pattern),ncol=length(pattern))
for (i in 1:max(pattern)) {
  if (odd(i)) {
    temp[i,] <- ((1+((i-1)*(length(pattern)))):(3+((i-1)*(length(pattern)))))
  }
  else {
    temp[i,] <- ((length(pattern)*i):(length(pattern)*i-2))
  }
}

#Create Parcels
for (i in 1:waves){
  if (odd(max(pattern))) newPattern <- pattern else {
    newPattern <- pattern
    newPattern[1] <- pattern[3]
    newPattern[2] <- pattern[2]
    newPattern[3] <- pattern[1]
  }
  for (j in 1:length(newPattern)) {
    parcel=0
    maxrows <- newPattern[j]  
    for (k in 1:maxrows) {
      parcel <- parcel+data[,((i-1)*items+itemList[temp[k,j],2])]
    }
    indexw <- paste("W",i,sep="")
    indexp <- paste("P",j,sep="")
    index <- paste(indexw,indexp,sep="")
    data[,paste("parcel",index,sep="")] <- parcel/newPattern[j]
  }
}

return(data[,(items*waves+1):(items*waves+length(pattern)*waves)])

}



################################
#####END PARCELS FUNCTION#######
################################
