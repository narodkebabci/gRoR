#' @title Visualize the distribution of the mutations
#'
#' @description This function takes the data frame and the mutation column to
#' visualize the frequency of the mutations.
#'
#' @import stringi corrplot
#'
#' @param data A data frame
#' @param mutations Specify the column name that contains the mutations
#'
#' @examples
#' wt2mt(data = df, mutations = "Mutation")
#'
#' @export

wt2mt <- function(data, mutations){

  # extract the wt in a vector
  wt <- vector("character", length = nrow(data))

  for (i in 1:nrow(data)){
    wt[i] <- stri_sub(data[, mutations][i], from=-stri_length(data[, mutations][i]),
                      to=-stri_length(data[, mutations][i]))
  }

  # extract the mt in a vector
  mt <- vector("character", length = nrow(data))

  for (j in 1:nrow(data)){
    mt[j] <- stri_sub(data[, mutations][j], from=stri_length(data[, mutations][j]),
                      to=stri_length(data[, mutations][j]))
  }

  q <- data.frame(cbind(wt, mt))
  colnames(q)=c("WT", "MT")

  # list of AminoAcids
  AA = c("A", "G", "P", "V", "L", "I", "M", "F", "W", "Y", "S", "T", "C", "N", "Q", "R", "H", "K", "D", "E")

  # initialize matrix with 0
  udx <- matrix(0, ncol = 20, nrow = 20)
  colnames(udx) <- AA
  rownames(udx) <- AA

  # create a mapping table
  df <- as.data.frame(matrix(1:20, ncol=20, nrow=1))
  colnames(df) <- AA

  for (m in 1:nrow(q)){

    l = as.character(q[,"WT"][m])
    k = as.character(q[,"MT"][m])

    lv1 = df[1, l]
    kv1 = df[1, k]

    udx[lv1, kv1] = 1 + udx[lv1, kv1]

  }

  return(corrplot(udx, method = "color", is.corr = FALSE, addgrid.col = TRUE, addCoef.col = "black",
                  tl.cex = 1, tl.col="black", tl.srt=45, cl.pos = "n"))

}


