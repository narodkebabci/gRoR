
class_formation <- function(data, mutations){

  # define aminoacid groups
  ali = c("G", "A", "P", "V", "L", "I", "M")
  aro = c("F", "W", "Y")
  p = c("S", "T", "C", "N", "Q")
  c = c("R", "H", "K", "D", "E")

  # define wt to mt classes
  # aliphatic to others
  ali_ali <- list()
  ali_aro <- list()
  ali_p <- list()
  ali_c <- list()

  # aromatic to others
  aro_ali <- list()
  aro_aro <- list()
  aro_p <- list()
  aro_c <- list()

  # polar to others
  p_ali <- list()
  p_aro <- list()
  p_p <- list()
  p_c <- list()

  # charged to others
  c_ali <- list()
  c_aro <- list()
  c_p <- list()
  c_c <- list()

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

  # add WT and MT vectors to the original dataframe
  q <- data.frame(cbind(data, wt, mt))

  # relocate the position of WT & MT vectors
  q <- q %>% relocate(wt, .after = mutations)
  q <- q %>% relocate(mt, .after = wt)

  WT <- "wt"
  MT <- "mt"

  for (i in 1:nrow(q)){

    # aliphatic to others
    if (q[, WT][i] %in% ali & q[, MT][i] %in% ali){
      ali_ali = rbind(ali_ali, q[i,])
    }else if (q[, WT][i] %in% ali & q[, MT][i] %in% aro){
      ali_aro = rbind(ali_aro, q[i,])
    }else if (q[, WT][i] %in% ali & q[, MT][i] %in% p){
      ali_p = rbind(ali_p, q[i,])
    }else if (q[, WT][i] %in% ali & q[, MT][i] %in% c){
      ali_c = rbind(ali_c, q[i,])
    }

    # aromatic to others
    else if (q[, WT][i] %in% aro & q[, MT][i] %in% ali){
      aro_ali = rbind(aro_ali, q[i,])
    }else if (q[, WT][i] %in% aro & q[, MT][i] %in% aro){
      aro_aro = rbind(aro_aro, q[i,])
    }else if (q[, WT][i] %in% aro & q[, MT][i] %in% p){
      aro_p = rbind(aro_p, q[i,])
    }else if (q[, WT][i] %in% aro & q[, MT][i] %in% c){
      aro_c = rbind(aro_c, q[i,])
    }

    # polar to others
    else if (q[, WT][i] %in% p & q[, MT][i] %in% ali){
      p_ali = rbind(p_ali, q[i,])
    }else if (q[, WT][i] %in% p & q[, MT][i] %in% aro){
      p_aro = rbind(p_aro, q[i,])
    }else if (q[, WT][i] %in% p & q[, MT][i] %in% p){
      p_p = rbind(p_p, q[i,])
    }else if (q[, WT][i] %in% p & q[, MT][i] %in% c){
      p_c = rbind(p_c, q[i,])
    }

    # charged to others
    else if (q[, WT][i] %in% c & q[, MT][i] %in% ali){
      c_ali = rbind(c_ali, q[i,])
    }else if (q[, WT][i] %in% c & q[, MT][i] %in% aro){
      c_aro = rbind(c_aro, q[i,])
    }else if (q[, WT][i] %in% c & q[, MT][i] %in% p){
      c_p = rbind(c_p, q[i,])
    }else if (q[, WT][i] %in% c & q[, MT][i] %in% c){
      c_c = rbind(c_c, q[i,])
    }
  }

  # label subsets
  ali_ali$Classes <- "Ali-Ali"
  ali_aro$Classes <- "Ali-Aro"
  ali_p$Classes <- "Ali-P"
  ali_c$Classes <- "Ali-C"
  aro_ali$Classes <- "Aro-Ali"
  aro_aro$Classes <- "Aro-Aro"
  aro_p$Classes <- "Aro-P"
  aro_c$Classes <- "Aro-C"
  p_ali$Classes <- "P-Ali"
  p_aro$Classes <- "P-Aro"
  p_p$Classes <- "P-P"
  p_c$Classes <- "P-C"
  c_ali$Classes <- "C-Ali"
  c_aro$Classes <- "C-Aro"
  c_p$Classes <- "C-P"
  c_c$Classes <- "C-C"

  # create a new data frame with class labels
  df <- data.frame(rbind(ali_ali, ali_aro, ali_p, ali_c,
                         aro_ali, aro_aro, aro_p, aro_c,
                         p_ali, p_aro, p_p, p_c,
                         c_ali, c_aro, c_p, c_c))

  df <- df %>% relocate(Classes, .after = mutations)

  return(df)
}
