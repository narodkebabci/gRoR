#' @title Build duplets based on mutations
#'
#' @description This function takes the data frame and mutation column to form WT-MT duplets
#'
#' @param data A data frame containing mutations that will be classified
#' @param mutations Specify the column name that contains the mutations
#'
#' @examples
#' duplet_formation(data = df, mutations = "Mutation")
#'
#' @export

duplet_formation <- function(data, mutations){

  # extract the wt in a vector
  WT <- vector("character", length = nrow(data))

  for (i in 1:nrow(data)){
    WT[i] <- stri_sub(data[, mutations][i], from=-stri_length(data[, mutations][i]),
                           to=-stri_length(data[, mutations][i]))
  }

  # extract the mt in a vector
  MT <- vector("character", length = nrow(data))

  for (j in 1:nrow(data)){
    MT[j] <- stri_sub(data[, mutations][j], from=stri_length(data[, mutations][j]),
                           to=stri_length(data[, mutations][j]))
  }

  # create a vector of duplets
  Duplets <- gsub(" ","",paste(WT, MT))

  # add Duplets vector to the original dataframe
  data[, "Duplets"] <- Duplets

  # relocate the position of Duplets vector
  data <- data %>% relocate(Duplets, .after = mutations)

  return(data)

}
