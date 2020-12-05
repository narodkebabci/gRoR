#' @title Compilation of a non-redundant dataset based on classes
#'
#' @description This function takes the data frame, class labels, experimental DDG values, and the mutations'
#' biochemical features. With this function, the dataset can be filtered based on class, secondary structure,
#' and relative ASA of the mutations.
#'
#' @param data A redundant dataset of mutations
#' @param classes Column name of labels
#' @param DDG Column name of Experimental DDG values
#' @param SS Default is Null. Specify the column name that contains secondary structure information
#' @param ASA Default is Null. Specify the column name that contains relative ASA information
#'
#' @examples
#' class_reduction(data = df, classes = "Classes", DDG = "Experimental")
#' class_reduction(data = df, classes = "Classes", DDG = "Experimental", SS = "SStructure")
#'
#' @return
#' A list of reduced dataset and a box plot
#'
#' @export

class_reduction <- function(data, classes, DDG, SS = NULL, ASA = NULL){

  # get the unique class names from the data (~16)
  label <- unique(data[, classes])

  idx <- matrix(ncol = ncol(data), nrow = 0)

  if (is.null(ASA) & !is.null(SS)){

    for (l in label){

      df <- data[which(data[, classes] %in% l), ]

      n <- 3

      sstructure <- unique(data[, SS])

      for (st in sstructure){

        df2 <- df[which(df[, SS] %in% st), ]

        # pick the median of df
        if (nrow(df2) %% 2 != 0){
          df2_median <- df2[df2[, DDG] == median(df2[, DDG]),]
        } else {
          df2_median <- df2[(nrow(df2)/2),]
        }

        # append the median to the list
        if (nrow(df2_median) > 1){
          df2_median <- sample_n(df2_median, 1, replace = TRUE)
          idx <- rbind(idx, df2_median)
        } else {
          idx <- rbind(idx, df2_median)
        }

        # pick the max & min values of m1
        if (nrow(df2) >= n){
          idx <- rbind(idx, df2[which.max(df2[, DDG]),], df2[which.min(df2[, DDG]),])
        }else {
          idx <- idx
        }
      }
    }

  }else if (is.null(SS) & !is.null(ASA)){

    for (l in label){

      df <- data[which(data[, classes] %in% l), ]

      n <- 3

      asa_1 <- df %>% filter(ASA < 0.3)

      # pick the median of asa_1
      if (nrow(asa_1) %% 2 != 0){
        asa_1_median <- asa_1[asa_1[, DDG] == median(asa_1[, DDG]),]
      } else {
        asa_1_median <- asa_1[(nrow(asa_1)/2),]
      }

      # append the median to the list
      if (nrow(asa_1_median) > 1){
        asa_1_median <- sample_n(asa_1_median, 1, replace = TRUE)
        idx <- rbind(idx, asa_1_median)
      } else {
        idx <- rbind(idx, asa_1_median)
      }

      # pick the max & min values of asa_1
      if (nrow(asa_1) >= n){
        idx <- rbind(idx, asa_1[which.max(asa_1[, DDG]),], asa_1[which.min(asa_1[, DDG]),])
      }else {
        idx <- idx
      }

      asa_2 <- df %>% filter(ASA >= 0.3 & ASA < 0.7)

      # pick the median of asa_2
      if (nrow(asa_2) %% 2 != 0){
        asa_2_median <- asa_2[asa_2[, DDG] == median(asa_2[, DDG]),]
      } else {
        asa_2_median <- asa_2[(nrow(asa_2)/2),]
      }

      # append the median to the list
      if (nrow(asa_2_median) > 1){
        asa_2_median <- sample_n(asa_2_median, 1, replace = TRUE)
        idx <- rbind(idx, asa_2_median)
      } else {
        idx <- rbind(idx, asa_2_median)
      }

      # pick the max & min values of asa_2
      if (nrow(asa_2) >= n){
        idx <- rbind(idx, asa_2[which.max(asa_2[, DDG]),], asa_2[which.min(asa_2[, DDG]),])
      }else {
        idx <- idx
      }

      asa_3 <- df %>% filter(ASA >= 0.7)

      # pick the median of asa_3
      if (nrow(asa_3) %% 2 != 0){
        asa_3_median <- asa_3[asa_3[, DDG] == median(asa_3[, DDG]),]
      } else {
        asa_3_median <- asa_3[(nrow(asa_3)/2),]
      }

      # append the median to the list
      if (nrow(asa_3_median) > 1){
        asa_3_median <- sample_n(asa_3_median, 1, replace = TRUE)
        idx <- rbind(idx, asa_3_median)
      } else {
        idx <- rbind(idx, asa_3_median)
      }

      # pick the max & min values of asa_3
      if (nrow(asa_3) >= 3){
        idx <- rbind(idx, asa_3[which.max(asa_3[, DDG]),], asa_3[which.min(asa_3[, DDG]),])
      }else {
        idx <- idx
      }
    }

  }else if (!is.null(SS) & !is.null(ASA)){

    sstructure <- unique(data[, SS])

    data[, "Labels"] <- as.character("")

    for (l in label){

      n <- 3

      df <- data[which(data[, classes] %in% l), ]

      for (st in sstructure){

        df2 <- df[which(df[, SS] %in% st), ]

        asa_1 <- df2 %>% filter(ASA < 0.3)

        if(nrow(asa_1) != 0){
          asa_1$Labels <- gsub(" ","-", paste(l, st, "asa_1"))
        } else {
          idx = idx
        }

        asa_2 <- df2 %>% filter(ASA >= 0.3 & ASA < 0.7)

        if(nrow(asa_2) != 0){
          asa_2$Labels <- gsub(" ","-", paste(l, st, "asa_2"))
        } else {
          idx = idx
        }

        asa_3 <- df2 %>% filter(ASA >= 0.7)

        if(nrow(asa_3) != 0){
          asa_3$Labels <- gsub(" ","-", paste(l, st, "asa_3"))
        } else {
          idx = idx
        }

        # pick the median of asa_1
        if (nrow(asa_1) %% 2 != 0){
          asa_1_median <- asa_1[asa_1[, DDG] == median(asa_1[, DDG]),]
        } else {
          asa_1_median <- asa_1[(nrow(asa_1)/2),]
        }

        # append the median to the list
        if (nrow(asa_1_median) > 1){
          asa_1_median <- sample_n(asa_1_median, 1, replace = TRUE)
          idx <- rbind(idx, asa_1_median)
        } else {
          idx <- rbind(idx, asa_1_median)
        }

        # pick the max & min values of asa_1
        if (nrow(asa_1) >= n){
          idx <- rbind(idx, asa_1[which.max(asa_1[, DDG]),], asa_1[which.min(asa_1[, DDG]),])
        }else {
          idx <- idx
        }

        # pick the median of asa_2
        if (nrow(asa_2) %% 2 != 0){
          asa_2_median <- asa_2[asa_2[, DDG] == median(asa_2[, DDG]),]
        } else {
          asa_2_median <- asa_2[(nrow(asa_2)/2),]
        }

        # append the median to the list
        if (nrow(asa_2_median) > 1){
          asa_2_median <- sample_n(asa_2_median, 1, replace = TRUE)
          idx <- rbind(idx, asa_2_median)
        } else {
          idx <- rbind(idx, asa_2_median)
        }

        # pick the max & min values of asa_2
        if (nrow(asa_2) >= n){
          idx <- rbind(idx, asa_2[which.max(asa_2[, DDG]),], asa_2[which.min(asa_2[, DDG]),])
        }else {
          idx <- idx
        }

        # pick the median of asa_3
        if (nrow(asa_3) %% 2 != 0){
          asa_3_median <- asa_3[asa_3[, DDG] == median(asa_3[, DDG]),]
        } else {
          asa_3_median <- asa_3[(nrow(asa_3)/2),]
        }

        # append the median to the list
        if (nrow(asa_3_median) > 1){
          asa_3_median <- sample_n(asa_3_median, 1, replace = TRUE)
          idx <- rbind(idx, asa_3_median)
        } else {
          idx <- rbind(idx, asa_3_median)
        }

        # pick the max & min values of asa_3
        if (nrow(asa_3) >= n){
          idx <- rbind(idx, asa_3[which.max(asa_3[, DDG]),], asa_3[which.min(asa_3[, DDG]),])
        }else {
          idx <- idx
        }

      }
    }

  }else {

    for (l in label){

      df <- data[which(data[, classes] %in% l), ]

      n <- 3

      # pick the median of df
      if (nrow(df) %% 2 != 0){
        df_median <- df[df[, DDG] == median(df[, DDG]),]
      } else {
        df_median <- df[(nrow(df)/2),]
      }

      # append the median to the list
      if (nrow(df_median) > 1){
        df_median <- sample_n(df_median, 1, replace = TRUE)
        idx <- rbind(idx, df_median)
      } else {
        idx <- rbind(idx, df_median)
      }

      # pick the max & min values of m1
      if (nrow(df) > n){
        idx <- rbind(idx, df[which.max(df[, DDG]),], df[which.min(df[, DDG]),])
      }else {
        idx <- idx
      }
    }
  }

  #return(idx)

  # label the datasets to create box plots

  if (nrow(idx) > nrow(data)){
    t <- nrow(data)
    r <- nrow(idx)

    data[, "Labels"] <- paste(paste("Total", "(n=", sep=" "), t, ")", sep="")
    names(idx)[names(idx) == "Labels"] <- paste(paste("Reduced", "(n=", sep=" "), r, ")", sep="")

    all <- rbind(data, idx)
    names(all)[names(all) == DDG] <- "Values"

  }else {

  t <- nrow(data)
  r <- nrow(idx)

  data[, "Labels"] <- paste(paste("Total", "(n=", sep=" "), t, ")", sep="")
  idx[, "Labels"] <- paste(paste("Reduced", "(n=", sep=" "), r, ")", sep="")

  all <- rbind(data, idx)
  names(all)[names(all) == DDG] <- "Values"

  }

  # create Box plot
  bp <- ggplot(all, aes(x=Labels, y=Values)) +
    geom_boxplot(aes(fill = Labels)) +
    labs(x="Labels", y="DDG Values") +
    theme_light() +
    theme(legend.position = "right")

  # return(bp)

  return(list(idx[,-c(ncol(idx))],bp))

}
