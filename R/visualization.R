
visualization <- function(data, DDG, ...){

  # label the datasets to create box plots
  data[, "Labels"] <- "data"

  # Box plot
  bp <- ggplot(data, aes(Labels, DDG)) +
    geom_boxplot(aes(fill = Labels)) +
    theme_minimal() +
    theme(legend.position = "top")

  return(bp)
}
