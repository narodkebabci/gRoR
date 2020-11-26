\name{class_formation}
\alias{class_formation}
\title{Create classes based on mutations}
\description{
This function takes the data frame and mutation column to form amino acid classes
}
\usage{
class_formation(data, mutations)
}
\arguments{
  \item{data}{A data frame containing mutations that will be classified}
  \item{mutations}{Specify the column name that contains the mutations}
}
\examples{
class_formation(data = df, mutations = "Mutation")
}
