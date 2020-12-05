#' @title Visualize the distribution of the classes
#'
#' @description This function takes the data frame and the mutation column to visualize the frequency of the
#' classes
#'
#' @param data A data frame
#' @param mutations Specify the column name that contains the mutations
#'
#' @examples
#' class_frequency(data = df, mutations = "Mutation")
#'
#' @export

class_frequency <- function(data, mutations){

  # extract the wt in a vector
  wt <- vector("character", length = nrow(data))

  for (i in 1:nrow(data)){
    wt[i] <- stri_sub(data[, "Mutation"][i], from=-stri_length(data[, "Mutation"][i]),
                      to=-stri_length(data[, "Mutation"][i]))
  }

  # extract the mt in a vector
  mt <- vector("character", length = nrow(data))

  for (j in 1:nrow(data)){
    mt[j] <- stri_sub(data[, "Mutation"][j], from=stri_length(data[, "Mutation"][j]),
                      to=stri_length(data[, "Mutation"][j]))
  }

  #create a dataframe from wt, mt
  q <- data.frame(cbind(wt, mt))
  colnames(q)=c("WT", "MT")

  q[, "WT"] <- as.character(q[, "WT"])
  q[, "MT"] <- as.character(q[, "MT"])

  # replace each aminoacid with its class label
  ali = c("G", "A", "P", "V", "L", "I", "M")
  aro = c("F", "W", "Y")
  p = c("S", "T", "C", "N", "Q")
  c = c("R", "H", "K", "D", "E")

  # replace each aminoacid with their classes
  q$WT[q$WT %in% ali] <- "Aliphatic"
  q$WT[q$WT %in% aro] <- "Aromatic"
  q$WT[q$WT %in% p] <- "Polar"
  q$WT[q$WT %in% c] <- "Charged"

  q$MT[q$MT %in% ali] <- "Aliphatic"
  q$MT[q$MT %in% aro] <- "Aromatic"
  q$MT[q$MT %in% p] <- "Polar"
  q$MT[q$MT %in% c] <- "Charged"

  # list of class labels
  CC = c("Aliphatic", "Aromatic", "Polar", "Charged")

  # initialize matrix with 0
  udx <- matrix(0, ncol = 4, nrow = 4)
  colnames(udx) <- CC
  rownames(udx) <- CC

  # create a mapping table
  df <- as.data.frame(matrix(1:4, ncol=4, nrow=1))
  colnames(df) <- CC

  for (m in 1:nrow(q)){

    l = as.character(q[,"WT"][m])
    k = as.character(q[,"MT"][m])

    lv1 = df[1, l]
    kv1 = df[1, k]

    udx[lv1, kv1] = 1 + udx[lv1, kv1]

}
  return(corrplot(udx, method = "color", type = "full", is.corr = FALSE, addgrid.col = TRUE, addCoef.col = "black",
         tl.pos = "lt", tl.cex = 1.2, tl.col="black", tl.srt=0, cl.pos = "n"))

}



