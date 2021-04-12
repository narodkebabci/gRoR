#' @title Compilation of a non-redundant dataset based on classes
#'
#' @description This function filters the dataset based on class labels,
#' secondary structure, and relative ASA of the mutations. This function takes
#' four inputs: the data frame to reduce, class labels, experimental DDG values,
#' and the specified difference between intervals in kcal/mol. Secondary structure
#' and relative ASA of the mutations are also applicable for reduction.
#'
#' @importFrom stats median na.omit
#' @import dplyr ggplot2
#'
#' @param data A redundant dataset of mutations
#' @param classes Column name of labels
#' @param DDG Column name of Experimental DDG values
#' @param r difference between each interval
#' @param SS Default is Null. Specify the column name that contains secondary structure information
#' @param ASA Default is Null. Specify the column name that contains relative ASA information
#'
#' @examples
#' ranged_class_reduction(data = df, classes = "Classes",
#'                        DDG = "Experimental", r = 2)
#' ranged_class_reduction(data = df, classes = "Classes",
#'                        DDG = "Experimental", r = 2, SS = "SStructure")
#' ranged_class_reduction(data = df, classes = "Classes",
#'                        DDG = "Experimental", r = 2, ASA = "ASA")
#'
#' @return
#' A list of reduced dataset and a box plot
#'
#' @export

ranged_class_reduction <- function(data, classes, DDG, r, SS = NULL, ASA = NULL){

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

        if (nrow(df2) != 0){

          # create intervals for each sub-classes

          dfmin <- min(df2[, DDG])
          dfmax <- max(df2[, DDG])

          # set x1 to the min value of the class and xn to the max value of the class
          x1 <- dfmin # set min
          xn <- dfmax # set max

          # create an empty, numeric vector to save the intervals
          interval <- vector()

          a <- 1

          while (x1 < xn){
            interval[a] <- x1
            x1 <- x1 + r
            a <- a + 1
          }

          interval[a] <- xn

          for (i in 1:(length(interval)-1)){

            if (i < (length(interval)-1)){

              df3 <- df2[df2[, DDG] >= interval[i] & df2[, DDG] < interval[i+1], ]

            }else {

              df3 <- df2[df2[, DDG] >= interval[i] & df2[, DDG] <= interval[i+1], ]

            }

            if (nrow(df3) > n){

              # pick the median of df
              if (nrow(df3) %% 2 != 0){
                df_median <- df3[df3[, DDG] == median(df3[, DDG]),]
              } else {
                df_median <- df3[(nrow(df3)/2),]
              }

              # append the median to the list
              if (nrow(df_median) > 1){
                df_median <- df_median[sample(nrow(df_median), 1), ]
                idx <- rbind(idx, df_median)
              } else {
                idx <- rbind(idx, df_median)
              }

              # pick the min of df
              df_min <- df3[df3[, DDG] == df3[, DDG][1], ]

              if (nrow(df_min) > 1){
                abs_min <- anti_join(df_min, df_median)
                df_min <- abs_min[sample(nrow(abs_min), 1), ]
                idx <- rbind(idx, df_min)
              }else {
                idx <- rbind(idx, df_min)
              }

              # pick the max of df
              df_max <- df3[df3[, DDG] == df3[, DDG][nrow(df3)], ]

              if (nrow(df_max) > 1){
                abs_max <- anti_join(df_max, df_median)
                df_max <- abs_max[sample(nrow(abs_max), 1), ]
                idx <- rbind(idx, df_max)
              }else {
                idx <- rbind(idx, df_max)
              }

            }else {
              idx <- rbind(idx, df3)
            }
          }

        }

        idx <- na.omit(idx)

      }
    }
  }else if (is.null(SS) & !is.null(ASA)){

    for (l in label){

      df <- data[which(data[, classes] %in% l), ]

      df <- df[order(df[, DDG]),][1:length(df[, DDG]),]

      n <- 3

      # create df1, ASA lower than 0.30
      df1 <- df[df[, ASA] <= 0.10, ]

      if (nrow(df1) != 0){

        # create intervals for each sub-classes

        dfmin <- min(df1[, DDG])
        dfmax <- max(df1[, DDG])

        # set x1 to the min value of the class and xn to the max value of the class
        x1 <- dfmin # set min
        xn <- dfmax # set max

        # create an empty, numeric vector to save the intervals
        interval <- vector()

        a <- 1

        while (x1 < xn){
          interval[a] <- x1
          x1 <- x1 + r
          a <- a + 1
        }

        interval[a] <- xn

        for (i in 1:(length(interval)-1)){

          if (i < (length(interval)-1)){

            df1_1 <- df1[df1[, DDG] >= interval[i] & df1[, DDG] < interval[i+1], ]

          }else {

            df1_1 <- df1[df1[, DDG] >= interval[i] & df1[, DDG] <= interval[i+1], ]

          }

          if (nrow(df1_1) > n){

            # pick the median of df
            if (nrow(df1_1) %% 2 != 0){
              df1_median <- df1_1[df1_1[, DDG] == median(df1_1[, DDG]),]
            } else {
              df1_median <- df1_1[(nrow(df1_1)/2),]
            }

            # append the median to the list
            if (nrow(df1_median) > 1){
              df1_median <- df1_median[sample(nrow(df1_median), 1), ]
              idx <- rbind(idx, df1_median)
            } else {
              idx <- rbind(idx, df1_median)
            }

            # pick the min of df
            df1_min <- df1_1[df1_1[, DDG] == df1_1[, DDG][1], ]

            if (nrow(df1_min) > 1){
              abs_min <- anti_join(df1_min, df1_median)
              df1_min <- abs_min[sample(nrow(abs_min), 1), ]
              idx <- rbind(idx, df1_min)
            }else {
              idx <- rbind(idx, df1_min)
            }

            # pick the max of df
            df1_max <- df1_1[df1_1[, DDG] == df1_1[, DDG][nrow(df1_1)], ]

            if (nrow(df1_max) > 1){
              abs_max <- anti_join(df1_max, df1_median)
              df1_max <- abs_max[sample(nrow(abs_max), 1), ]
              idx <- rbind(idx, df1_max)
            }else {
              idx <- rbind(idx, df1_max)
            }

          }else {
            idx <- rbind(idx, df1_1)
          }
        }
      }

      # create df2, ASA between 0.30 and 0.70
      df2 <- df[df[, ASA] > 0.10 & df[, ASA] < 0.50, ]

      if (nrow(df2) != 0){

        # create intervals for each sub-classes

        dfmin <- min(df2[, DDG])
        dfmax <- max(df2[, DDG])

        # set x1 to the min value of the class and xn to the max value of the class
        x1 <- dfmin # set min
        xn <- dfmax # set max

        # create an empty, numeric vector to save the intervals
        interval <- vector()

        a <- 1

        while (x1 < xn){
          interval[a] <- x1
          x1 <- x1 + r
          a <- a + 1
        }

        interval[a] <- xn

        for (j in 1:(length(interval)-1)){

          if (j < (length(interval)-1)){

            # df3 <- df2 %>% dplyr::filter(Experimental >= interval[i] & Experimental < interval[i+1])
            df2_1 <- df2[df2[, DDG] >= interval[j] & df2[, DDG] < interval[j+1], ]

          }else {

            df2_1 <- df2[df2[, DDG] >= interval[j] & df2[, DDG] <= interval[j+1], ]

          }

          if (nrow(df2_1) > n){

            # pick the median of df
            if (nrow(df2_1) %% 2 != 0){
              df2_median <- df2_1[df2_1[, DDG] == median(df2_1[, DDG]),]
            } else {
              df2_median <- df2_1[(nrow(df2_1)/2),]
            }

            # append the median to the list
            if (nrow(df2_median) > 1){
              df2_median <- df2_median[sample(nrow(df2_median), 1), ]
              idx <- rbind(idx, df2_median)
            } else {
              idx <- rbind(idx, df2_median)
            }

            # pick the min of df
            df2_min <- df2_1[df2_1[, DDG] == df2_1[, DDG][1], ]

            if (nrow(df2_min) > 1){
              abs_min <- anti_join(df2_min, df2_median)
              df2_min <- abs_min[sample(nrow(abs_min), 1), ]
              idx <- rbind(idx, df2_min)
            }else {
              idx <- rbind(idx, df2_min)
            }

            # pick the max of df
            df2_max <- df2_1[df2_1[, DDG] == df2_1[, DDG][nrow(df2_1)], ]

            if (nrow(df2_max) > 1){
              abs_max <- anti_join(df2_max, df2_median)
              df2_max <- abs_max[sample(nrow(abs_max), 1), ]
              idx <- rbind(idx, df2_max)
            }else {
              idx <- rbind(idx, df2_max)
            }

          }else {
            idx <- rbind(idx, df2_1)
          }
        }
      }

      # create df3, ASA greater than 0.70
      df3 <- df[df[, ASA] >= 0.50, ]

      if (nrow(df3) != 0){

        # create intervals for each sub-classes

        dfmin <- min(df3[, DDG])
        dfmax <- max(df3[, DDG])

        # set x1 to the min value of the class and xn to the max value of the class
        x1 <- dfmin # set min
        xn <- dfmax # set max

        # create an empty, numeric vector to save the intervals
        interval <- vector()

        a <- 1

        while (x1 < xn){
          interval[a] <- x1
          x1 <- x1 + r
          a <- a + 1
        }

        interval[a] <- xn

        for (k in 1:(length(interval)-1)){

          if (k < (length(interval)-1)){

            # df3 <- df2 %>% dplyr::filter(Experimental >= interval[i] & Experimental < interval[i+1])
            df3_1 <- df3[df3[, DDG] >= interval[k] & df3[, DDG] < interval[k+1], ]

          }else {

            df3_1 <- df3[df3[, DDG] >= interval[k] & df3[, DDG] <= interval[k+1], ]

          }

          if (nrow(df3_1) > n){

            # pick the median of df
            if (nrow(df3_1) %% 2 != 0){
              df3_median <- df3_1[df3_1[, DDG] == median(df3_1[, DDG]),]
            } else {
              df3_median <- df3_1[(nrow(df3_1)/2),]
            }

            # append the median to the list
            if (nrow(df3_median) > 1){
              df3_median <- df3_median[sample(nrow(df3_median), 1), ]
              idx <- rbind(idx, df3_median)
            } else {
              idx <- rbind(idx, df3_median)
            }

            # pick the min of df
            df3_min <- df3_1[df3_1[, DDG] == df3_1[, DDG][1], ]

            if (nrow(df3_min) > 1){
              abs_min <- anti_join(df3_min, df3_median)
              df3_min <- abs_min[sample(nrow(abs_min), 1), ]
              idx <- rbind(idx, df3_min)
            }else {
              idx <- rbind(idx, df3_min)
            }

            # pick the max of df
            df3_max <- df3_1[df3_1[, DDG] == df3_1[, DDG][nrow(df3_1)], ]

            if (nrow(df3_max) > 1){
              abs_max <- anti_join(df3_max, df3_median)
              df3_max <- abs_max[sample(nrow(abs_max), 1), ]
              idx <- rbind(idx, df3_max)
            }else {
              idx <- rbind(idx, df3_max)
            }

          }else {
            idx <- rbind(idx, df3_1)
          }
        }
      }

      idx <- na.omit(idx)

    }
  }else if (!is.null(SS) & !is.null(ASA)){

    sstructure <- unique(data[, SS])

    for (l in label){

      n <- 3

      df <- data[which(data[, classes] %in% l), ]

      df <- df[order(df[, DDG]),][1:length(df[, DDG]),]

      for (st in sstructure){

        ddf <- df[which(df[, SS] %in% st), ]

        df1 <- ddf[ddf[, ASA] <= 0.10, ]

        if (nrow(df1) != 0){

          # create intervals for each sub-classes

          dfmin <- min(df1[, DDG])
          dfmax <- max(df1[, DDG])

          # set x1 to the min value of the class and xn to the max value of the class
          x1 <- dfmin # set min
          xn <- dfmax # set max

          # create an empty, numeric vector to save the intervals
          interval <- vector()

          a <- 1

          while (x1 < xn){
            interval[a] <- x1
            x1 <- x1 + r
            a <- a + 1
          }

          interval[a] <- xn

          for (i in 1:(length(interval)-1)){

            if (i < (length(interval)-1)){

              df1_1 <- df1[df1[, DDG] >= interval[i] & df1[, DDG] < interval[i+1], ]

            }else {

              df1_1 <- df1[df1[, DDG] >= interval[i] & df1[, DDG] <= interval[i+1], ]

            }

            if (nrow(df1_1) > n){

              # pick the median of df
              if (nrow(df1_1) %% 2 != 0){
                df1_median <- df1_1[df1_1[, DDG] == median(df1_1[, DDG]),]
              } else {
                df1_median <- df1_1[(nrow(df1_1)/2),]
              }

              # append the median to the list
              if (nrow(df1_median) > 1){
                df1_median <- df1_median[sample(nrow(df1_median), 1), ]
                idx <- rbind(idx, df1_median)
              } else {
                idx <- rbind(idx, df1_median)
              }

              # pick the min of df
              df1_min <- df1_1[df1_1[, DDG] == df1_1[, DDG][1], ]

              if (nrow(df1_min) > 1){
                abs_min <- anti_join(df1_min, df1_median)
                df1_min <- abs_min[sample(nrow(abs_min), 1), ]
                idx <- rbind(idx, df1_min)
              }else {
                idx <- rbind(idx, df1_min)
              }

              # pick the max of df
              df1_max <- df1_1[df1_1[, DDG] == df1_1[, DDG][nrow(df1_1)], ]

              if (nrow(df1_max) > 1){
                abs_max <- anti_join(df1_max, df1_median)
                df1_max <- abs_max[sample(nrow(abs_max), 1), ]
                idx <- rbind(idx, df1_max)
              }else {
                idx <- rbind(idx, df1_max)
              }

            }else {
              idx <- rbind(idx, df1_1)
            }
          }
        }

        df2 <- ddf[ddf[, ASA] > 0.10 & ddf[, ASA] < 0.50, ]

        if (nrow(df2) != 0){

          # create intervals for each sub-classes

          dfmin <- min(df2[, DDG])
          dfmax <- max(df2[, DDG])

          # set x1 to the min value of the class and xn to the max value of the class
          x1 <- dfmin # set min
          xn <- dfmax # set max

          # create an empty, numeric vector to save the intervals
          interval <- vector()

          a <- 1

          while (x1 < xn){
            interval[a] <- x1
            x1 <- x1 + r
            a <- a + 1
          }

          interval[a] <- xn

          for (j in 1:(length(interval)-1)){

            if (j < (length(interval)-1)){

              # df3 <- df2 %>% dplyr::filter(Experimental >= interval[i] & Experimental < interval[i+1])
              df2_1 <- df2[df2[, DDG] >= interval[j] & df2[, DDG] < interval[j+1], ]

            }else {

              df2_1 <- df2[df2[, DDG] >= interval[j] & df2[, DDG] <= interval[j+1], ]

            }

            if (nrow(df2_1) > 3){

              # pick the median of df
              if (nrow(df2_1) %% 2 != 0){
                df2_median <- df2_1[df2_1[, DDG] == median(df2_1[, DDG]),]
              } else {
                df2_median <- df2_1[(nrow(df2_1)/2),]
              }

              # append the median to the list
              if (nrow(df2_median) > 1){
                df2_median <- df2_median[sample(nrow(df2_median), 1), ]
                idx <- rbind(idx, df2_median)
              } else {
                idx <- rbind(idx, df2_median)
              }

              # pick the min of df
              df2_min <- df2_1[df2_1[, DDG] == df2_1[, DDG][1], ]

              if (nrow(df2_min) > 1){
                abs_min <- anti_join(df2_min, df2_median)
                df2_min <- abs_min[sample(nrow(abs_min), 1), ]
                idx <- rbind(idx, df2_min)
              }else {
                idx <- rbind(idx, df2_min)
              }

              # pick the max of df
              df2_max <- df2_1[df2_1[, DDG] == df2_1[, DDG][nrow(df2_1)], ]

              if (nrow(df2_max) > 1){
                abs_max <- anti_join(df2_max, df2_median)
                df2_max <- abs_max[sample(nrow(abs_max), 1), ]
                idx <- rbind(idx, df2_max)
              }else {
                idx <- rbind(idx, df2_max)
              }

            }else {
              idx <- rbind(idx, df2_1)
            }
          }
        }

        df3 <- ddf[ddf[, ASA] >= 0.50, ]

        if (nrow(df3) != 0){

          # create intervals for each sub-classes

          dfmin <- min(df3[, DDG])
          dfmax <- max(df3[, DDG])

          # set x1 to the min value of the class and xn to the max value of the class
          x1 <- dfmin # set min
          xn <- dfmax # set max

          # create an empty, numeric vector to save the intervals
          interval <- vector()

          a <- 1

          while (x1 < xn){
            interval[a] <- x1
            x1 <- x1 + r
            a <- a + 1
          }

          interval[n] <- xn

          for (k in 1:(length(interval)-1)){

            if (k < (length(interval)-1)){

              df3_1 <- df3[df3[, DDG] >= interval[k] & df3[, DDG] < interval[k+1], ]

            }else {

              df3_1 <- df3[df3[, DDG] >= interval[k] & df3[, DDG] <= interval[k+1], ]

            }

            if (nrow(df3_1) > 3){

              # pick the median of df
              if (nrow(df3_1) %% 2 != 0){
                df3_median <- df3_1[df3_1[, DDG] == median(df3_1[, DDG]),]
              } else {
                df3_median <- df3_1[(nrow(df3_1)/2),]
              }

              # append the median to the list
              if (nrow(df3_median) > 1){
                df3_median <- df3_median[sample(nrow(df3_median), 1), ]
                idx <- rbind(idx, df3_median)
              } else {
                idx <- rbind(idx, df3_median)
              }

              # pick the min of df
              df3_min <- df3_1[df3_1[, DDG] == df3_1[, DDG][1], ]

              if (nrow(df3_min) > 1){
                abs_min <- anti_join(df3_min, df3_median)
                df3_min <- abs_min[sample(nrow(abs_min), 1), ]
                idx <- rbind(idx, df3_min)
              }else {
                idx <- rbind(idx, df3_min)
              }

              # pick the max of df
              df3_max <- df3_1[df3_1[, DDG] == df3_1[, DDG][nrow(df3_1)], ]

              if (nrow(df3_max) > 1){
                abs_max <- anti_join(df3_max, df3_median)
                df3_max <- abs_max[sample(nrow(abs_max), 1), ]
                idx <- rbind(idx, df3_max)
              }else {
                idx <- rbind(idx, df3_max)
              }

            }else {
              idx <- rbind(idx, df3_1)
            }
          }
        }

        idx <- na.omit(idx)

      }
    }

  }else{

    for (l in label){

      df <- data[which(data[, classes] %in% l), ]

      df <- df[order(df[, DDG]),][1:length(df[, DDG]),]

      n <- 3

      # create intervals for each sub-classes

      dfmin <- min(df[, DDG])
      dfmax <- max(df[, DDG])

      # set x1 to the min value of the class and xn to the max value of the class
      x1 <- dfmin # set min
      xn <- dfmax # set max

      # create an empty, numeric vector to save the intervals
      interval <- vector()

      a <- 1

      while (x1 < xn){
        interval[a] <- x1
        x1 <- x1 + r
        a <- a + 1
      }

      interval[a] <- xn

      for (i in 1:(length(interval)-1)){

        df2 <- df %>% dplyr::filter(Experimental >= interval[i] & Experimental < interval[i+1])

        if (nrow(df2) > n){

          # pick the median of df
          if (nrow(df2) %% 2 != 0){
            df_median <- df2[df2[, DDG] == median(df2[, DDG]),]
          } else {
            df_median <- df2[(nrow(df2)/2),]
          }

          # append the median to the list
          if (nrow(df_median) > 1){
            df_median <- df_median[sample(nrow(df_median), 1), ]
            idx <- rbind(idx, df_median)
          } else {
            idx <- rbind(idx, df_median)
          }

          # pick the min of df
          df_min <- df2[df2[, DDG] == df2[, DDG][1], ]

          if (nrow(df_min) > 1){
            abs_min <- anti_join(df_min, df_median)
            df_min <- abs_min[sample(nrow(abs_min), 1), ]
            idx <- rbind(idx, df_min)
          }else {
            idx <- rbind(idx, df_min)
          }

          # pick the max of df
          df_max <- df2[df2[, DDG] == df2[, DDG][nrow(df2)], ]

          if (nrow(df_max) > 1){
            abs_max <- anti_join(df_max, df_median)
            df_max <- abs_max[sample(nrow(abs_max), 1), ]
            idx <- rbind(idx, df_max)
          }else {
            idx <- rbind(idx, df_max)
          }

        }else {
          idx <- rbind(idx, df2)
        }
      }

      idx <- na.omit(idx)

      # replace xn with the max of the interval
      idx[nrow(idx),] <- df[which.max(df[, DDG]),]

    }

  }

  # label the datasets to create box plots
  q <- nrow(data)
  w <- nrow(idx)

  data[, "Plot_Labels"] <- paste(paste("Total", "(n=", sep=" "), q, ")", sep="")
  idx[, "Plot_Labels"] <- paste(paste("Reduced", "(n=", sep=" "), w, ")", sep="")

  all <- rbind(data, idx)
  names(all)[names(all) == DDG] <- "Values"

  # create Box plot
  bp <- ggplot(all, aes(x=Plot_Labels, y=Values)) +
    geom_boxplot(aes(fill = Plot_Labels)) +
    labs(x="Labels", y="DDG Values") +
    theme_light() +
    theme(legend.position = "right", legend.title = element_blank(),
          axis.title.x = element_blank())

  # return(bp)

  return(list(idx[,-c(ncol(idx))], bp))

}
