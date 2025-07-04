% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/source_varsel.R
\name{all_subsets}
\alias{all_subsets}
\title{Compute all subsets of a set}
\usage{
all_subsets(set)
}
\arguments{
\item{set}{the set from which to compute all subsets (e.g., \code{1:p})}
}
\value{
a data frame where the rows indicate the \code{2^p} different subsets
and the columns indicate inclusion (logical) for each element in that subset
}
\description{
Given a set of variables, compute the inclusion indicators for
all possible subsets.
}
\references{
Code adapted from <https://www.r-bloggers.com/2012/04/generating-all-subsets-of-a-set/>
}
