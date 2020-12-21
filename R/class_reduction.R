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

      df <- df[order(df[, DDG]),][1:length(df[, DDG]),]

      n <- 3

      sstructure <- unique(data[, SS])

      for (st in sstructure){

        df2 <- df[which(df[, SS] %in% st), ]

        if (nrow(df2) > n){

          # pick the median of df2
          if (nrow(df2) %% 2 != 0){
            df2_median <- df2[df2[, DDG] == median(df2[, DDG]),]
          } else {
            df2_median <- df2[(nrow(df2)/2),]
          }

          # append the median to the list
          if (nrow(df2_median) > 1){
            df2_median <- df2_median[sample(nrow(df2_median), 1), ]
            idx <- rbind(idx, df2_median)
          } else {
            idx <- rbind(idx, df2_median)
          }

          # pick the min of df2
          df_min <- df2[df2[, DDG] == df2[, DDG][1], ]

          if (nrow(df_min) > 1){
            abs_min <- anti_join(df_min, df2_median)
            df_min <- abs_min[sample(nrow(abs_min), 1), ]
            idx <- rbind(idx, df_min)
          }else {
            idx <- rbind(idx, df_min)
          }

          # pick the max of df2
          df_max <- df2[df2[, DDG] == df2[, DDG][nrow(df2)], ]

          if (nrow(df_max) > 1){
            abs_max <- anti_join(df_max, df2_median)
            df_max <- abs_max[sample(nrow(abs_max), 1), ]
            idx <- rbind(idx, df_max)
          }else {
            idx <- rbind(idx, df_max)
          }

        }else{
          idx <- rbind(idx, df2)
        }
      }
    }

  }else if (is.null(SS) & !is.null(ASA)){

    for (l in label){

      df <- data[which(data[, classes] %in% l), ]

      df <- df[order(df[, DDG]),][1:length(df[, DDG]),]

      n <- 3

      df1 <- df[df[, ASA] < 0.30, ]
      df2 <- df[df[, ASA] >= 0.30 & df[, ASA] < 0.70, ]
      df3 <- df[df[, ASA] >= 0.70, ]

      if (nrow(df1) > n){

        # pick the median of df1
        if (nrow(df1) %% 2 != 0){
          df1_median <- df1[df1[, DDG] == median(df1[, DDG]),]
        } else {
          df1_median <- df1[(nrow(df1)/2),]
        }

        # append the median to the list
        if (nrow(df1_median) > 1){
          df1_median <- df1_median[sample(nrow(df1_median), 1), ]
          idx <- rbind(idx, df1_median)
        } else {
          idx <- rbind(idx, df1_median)
        }

        # pick the min of df1
        df1_min <- df1[df1[, DDG] == df1[, DDG][1], ]

        if (nrow(df1_min) > 1){
          abs_min <- anti_join(df1_min, df1_median)
          df1_min <- abs_min[sample(nrow(abs_min), 1), ]
          idx <- rbind(idx, df1_min)
        }else {
          idx <- rbind(idx, df1_min)
        }

        # pick the max of df1
        df1_max <- df1[df1[, DDG] == df1[, DDG][nrow(df1)], ]

        if (nrow(df1_max) > 1){
          abs_max <- anti_join(df1_max, df1_median)
          df1_max <- abs_max[sample(nrow(abs_max), 1), ]
          idx <- rbind(idx, df1_max)
        }else {
          idx <- rbind(idx, df1_max)
        }

      }else{
        idx <- rbind(idx, df1)
      }

      if (nrow(df2) > n){

        # pick the median of df2
        if (nrow(df2) %% 2 != 0){
          df2_median <- df2[df2[, DDG] == median(df2[, DDG]),]
        } else {
          df2_median <- df2[(nrow(df2)/2),]
        }

        # append the median to the list
        if (nrow(df2_median) > 1){
          df2_median <- df2_median[sample(nrow(df2_median), 1), ]
          idx <- rbind(idx, df2_median)
        } else {
          idx <- rbind(idx, df2_median)
        }

        # pick the min of df2
        df2_min <- df2[df2[, DDG] == df2[, DDG][1], ]

        if (nrow(df2_min) > 1){
          abs_min <- anti_join(df2_min, df2_median)
          df2_min <- abs_min[sample(nrow(abs_min), 1), ]
          idx <- rbind(idx, df2_min)
        }else {
          idx <- rbind(idx, df2_min)
        }

        # pick the max of df2
        df2_max <- df2[df2[, DDG] == df2[, DDG][nrow(df2)], ]

        if (nrow(df2_max) > 1){
          abs_max <- anti_join(df2_max, df2_median)
          df2_max <- abs_max[sample(nrow(abs_max), 1), ]
          idx <- rbind(idx, df2_max)
        }else {
          idx <- rbind(idx, df2_max)
        }

      }else{
        idx <- rbind(idx, df2)
      }

      if (nrow(df3) > n){

        # pick the median of df3
        if (nrow(df3) %% 2 != 0){
          df3_median <- df3[df3[, DDG] == median(df3[, DDG]),]
        } else {
          df3_median <- df3[(nrow(df3)/2),]
        }

        # append the median to the list
        if (nrow(df3_median) > 1){
          df3_median <- df3_median[sample(nrow(df3_median), 1), ]
          idx <- rbind(idx, df3_median)
        } else {
          idx <- rbind(idx, df3_median)
        }

        # pick the min of df3
        df3_min <- df3[df3[, DDG] == df3[, DDG][1], ]

        if (nrow(df3_min) > 1){
          abs_min <- anti_join(df3_min, df3_median)
          df3_min <- abs_min[sample(nrow(abs_min), 1), ]
          idx <- rbind(idx, df3_min)
        }else {
          idx <- rbind(idx, df3_min)
        }

        # pick the max of df3
        df3_max <- df3[df3[, DDG] == df3[, DDG][nrow(df3)], ]

        if (nrow(df3_max) > 1){
          abs_max <- anti_join(df3_max, df3_median)
          df3_max <- abs_max[sample(nrow(abs_max), 1), ]
          idx <- rbind(idx, df3_max)
        }else {
          idx <- rbind(idx, df3_max)
        }

      }else{
        idx <- rbind(idx, df3)
      }

    }

  }else if (!is.null(SS) & !is.null(ASA)){

    sstructure <- unique(data[, SS])

    data[, "Labels"] <- as.character("")

    for (l in label){

      n <- 3

      df <- data[which(data[, classes] %in% l), ]

      df <- df[order(df[, DDG]),][1:length(df[, DDG]),]

      for (st in sstructure){

        ddf <- df[which(df[, SS] %in% st), ]

        df1 <- ddf[ddf[, ASA] < 0.30, ]
        df2 <- ddf[ddf[, ASA] >= 0.30 & ddf[, ASA] < 0.70, ]
        df3 <- ddf[ddf[, ASA] >= 0.70, ]

        # create label for each line
        if(nrow(df1) != 0){
          df1$Labels <- gsub(" ","-", paste(l, st, "asa_1"))
        } else {
          idx <- idx
        }

        if(nrow(df2) != 0){
          df2$Labels <- gsub(" ","-", paste(l, st, "asa_2"))
        } else {
          idx <- idx
        }

        if(nrow(df3) != 0){
          df3$Labels <- gsub(" ","-", paste(l, st, "asa_3"))
        } else {
          idx <- idx
        }

        if (nrow(df1) > n){

          # pick the median of df2
          if (nrow(df1) %% 2 != 0){
            df1_median <- df1[df1[, DDG] == median(df1[, DDG]),]
          } else {
            df1_median <- df1[(nrow(df1)/2),]
          }

          # append the median to the list
          if (nrow(df1_median) > 1){
            df1_median <- df1_median[sample(nrow(df1_median), 1), ]
            idx <- rbind(idx, df1_median)
          } else {
            idx <- rbind(idx, df1_median)
          }

          # pick the min of df1
          df1_min <- df1[df1[, DDG] == df1[, DDG][1], ]

          if (nrow(df1_min) > 1){
            abs_min <- anti_join(df1_min, df1_median)
            df1_min <- abs_min[sample(nrow(abs_min), 1), ]
            idx <- rbind(idx, df1_min)
          }else {
            idx <- rbind(idx, df1_min)
          }

          # pick the max of df1
          df1_max <- df1[df1[, DDG] == df1[, DDG][nrow(df1)], ]

          if (nrow(df1_max) > 1){
            abs_max <- anti_join(df1_max, df1_median)
            df1_max <- abs_max[sample(nrow(abs_max), 1), ]
            idx <- rbind(idx, df1_max)
          }else {
            idx <- rbind(idx, df1_max)
          }

        }else{
          idx <- rbind(idx, df1)
        }

        if (nrow(df2) > n){

          # pick the median of df2
          if (nrow(df2) %% 2 != 0){
            df2_median <- df2[df2[, DDG] == median(df2[, DDG]),]
          } else {
            df2_median <- df2[(nrow(df2)/2),]
          }

          # append the median to the list
          if (nrow(df2_median) > 1){
            df2_median <- df2_median[sample(nrow(df2_median), 1), ]
            idx <- rbind(idx, df2_median)
          } else {
            idx <- rbind(idx, df2_median)
          }

          # pick the min of df2
          df2_min <- df2[df2[, DDG] == df2[, DDG][1], ]

          if (nrow(df2_min) > 1){
            abs_min <- anti_join(df2_min, df2_median)
            df2_min <- abs_min[sample(nrow(abs_min), 1), ]
            idx <- rbind(idx, df2_min)
          }else {
            idx <- rbind(idx, df2_min)
          }

          # pick the max of df2
          df2_max <- df2[df2[, DDG] == df2[, DDG][nrow(df2)], ]

          if (nrow(df2_max) > 1){
            abs_max <- anti_join(df2_max, df2_median)
            df2_max <- abs_max[sample(nrow(abs_max), 1), ]
            idx <- rbind(idx, df2_max)
          }else {
            idx <- rbind(idx, df2_max)
          }

        }else{
          idx <- rbind(idx, df2)
        }

        if (nrow(df3) > n){

          # pick the median of df3
          if (nrow(df3) %% 2 != 0){
            df3_median <- df3[df3[, DDG] == median(df3[, DDG]),]
          } else {
            df3_median <- df3[(nrow(df3)/2),]
          }

          # append the median to the list
          if (nrow(df3_median) > 1){
            df3_median <- df3_median[sample(nrow(df3_median), 1), ]
            idx <- rbind(idx, df3_median)
          } else {
            idx <- rbind(idx, df3_median)
          }

          # pick the min of df3
          df3_min <- df3[df3[, DDG] == df3[, DDG][1], ]

          if (nrow(df3_min) > 1){
            abs_min <- anti_join(df3_min, df3_median)
            df3_min <- abs_min[sample(nrow(abs_min), 1), ]
            idx <- rbind(idx, df3_min)
          }else {
            idx <- rbind(idx, df3_min)
          }

          # pick the max of df3
          df3_max <- df3[df3[, DDG] == df3[, DDG][nrow(df3)], ]

          if (nrow(df3_max) > 1){
            abs_max <- anti_join(df3_max, df3_median)
            df3_max <- abs_max[sample(nrow(abs_max), 1), ]
            idx <- rbind(idx, df3_max)
          }else {
            idx <- rbind(idx, df3_max)
          }

        }else{
          idx <- rbind(idx, df3)
        }
      }
    }

  }else {

    for (l in label){

      df <- data[which(data[, classes] %in% l), ]

      df <- df[order(df[, DDG]),][1:length(df[, DDG]),]

      n <- 3

      if (nrow(df) > n){

        # pick the median of df
        if (nrow(df) %% 2 != 0){
          df_median <- df[df[, DDG] == median(df[, DDG]),]
        } else {
          df_median <- df[(nrow(df)/2),]
        }

        # append the median to the list
        if (nrow(df_median) > 1){
          df_median <- df_median[sample(nrow(df_median), 1), ]
          idx <- rbind(idx, df_median)
        } else {
          idx <- rbind(idx, df_median)
        }

        # pick the min of df
        df_min <- df[df[, DDG] == df[, DDG][1], ]

        if (nrow(df_min) > 1){
          abs_min <- anti_join(df_min, df_median)
          df_min <- abs_min[sample(nrow(abs_min), 1), ]
          idx <- rbind(idx, df_min)
        }else {
          idx <- rbind(idx, df_min)
        }

        # pick the max of df
        df_max <- df[df[, DDG] == df[, DDG][nrow(df)], ]

        if (nrow(df_max) > 1){
          abs_max <- anti_join(df_max, df_median)
          df_max <- abs_max[sample(nrow(abs_max), 1), ]
          idx <- rbind(idx, df_max)
        }else {
          idx <- rbind(idx, df_max)
        }

      }else {
        idx <- rbind(idx, df)
      }
    }
  }

  idx <- na.omit(idx)

  # label the datasets to create box plots
  t <- nrow(data)
  r <- nrow(idx)

  data[, "Plot_Labels"] <- paste(paste("Total", "(n=", sep=" "), t, ")", sep="")
  idx[, "Plot_Labels"] <- paste(paste("Reduced", "(n=", sep=" "), r, ")", sep="")

  all <- rbind(data, idx)
  names(all)[names(all) == DDG] <- "Values"

  # create Box plot
  bp <- ggplot(all, aes(x=Plot_Labels, y=Values)) +
    geom_boxplot(aes(fill = Plot_Labels)) +
    labs(x="Labels", y="DDG Values") +
    theme_light() +
    theme(legend.position = "right")

  # return(bp)

  return(list(idx[,-c(ncol(idx))],bp))

}
