################################
#####CREATE PARCEL FUNCTION#####
################################

#' Creates parcels from items
#'
#' `parcel()` takes a set of items from a single scale and creates parcels based
#' on the user specification. The function conducts a factor analysis of items
#' at each wave and then calculates the average factor loading for each item.
#' The item with the strongest loading is assigned to the first parcel, the
#' item with the second strongest loading is assigned to the second parcel,
#' and so on until all items are assigned. 
#'
#'
#' @param data Dataframe with just the items to be parcelled.
#' @param items Number of items for the scale.
#' @param waves Number of waves.
#' @param pattern Vector specifying the number of parcels and the number of
#'   items per parcel. For example, to get three parcels with 3 items for the
#'   first parcel, 3 items for the second, and 4 items for the third, specify
#'   `pattern = c(3, 3, 4)`
#'
#' @export
parcel <- function(data,items,waves,pattern) {

#Create Temporary Matrix
loadingSum <- matrix(0,items,1)

#Get Average Loading
for (i in 1:waves) {
  start=1+(i-1)*items
  finish=i*items
  loadingSum <- loadingSum+stats::factanal(na.omit(data[,start:finish]),
                                           1)$loadings[,1]
}
loadingAvg <- loadingSum/waves

#Rank Items by Loading
itemList <- data.frame(cbind(loadingAvg,matrix(1:items,items,1)))
itemList <- itemList[order(itemList$X1),]

#Create Pattern Matrix for Creating Parcels
temp <- matrix(nrow=max(pattern),ncol=length(pattern))
for (i in 1:max(pattern)) {
  if (gtools::odd(i)) {
    temp[i,] <- ((1+((i-1)*(length(pattern)))):(3+((i-1)*(length(pattern)))))
  }
  else {
    temp[i,] <- ((length(pattern)*i):(length(pattern)*i-2))
  }
}

#Create Parcels
for (i in 1:waves){
  if (gtools::odd(max(pattern))) newPattern <- pattern else {
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
