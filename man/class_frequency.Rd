% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/class_frequency.R
\name{class_frequency}
\alias{class_frequency}
\title{Visualize the distribution of the classes}
\usage{
class_frequency(data, mutations)
}
\arguments{
\item{data}{A data frame}

\item{mutations}{Specify the column name that contains the mutations}
}
\description{
This function takes the data frame and the mutation column
to visualize the frequency of the classes
}
\examples{
class_frequency(data = df, mutations = "Mutation")

}
