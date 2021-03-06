% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FuzzyInferenceSystem.R
\name{fis.builder}
\alias{fis.builder}
\title{TSK FIS builder}
\usage{
fis.builder(
  x.range,
  input.num,
  input.mf.num,
  input.mf.type,
  rule.num = prod(input.mf.num),
  rule.which = NULL,
  defuzzMethod = "default",
  params.ante,
  params.conse
)
}
\arguments{
\item{x.range}{a vector/matrix as the range of input(s)}

\item{input.num}{the number of inputs}

\item{input.mf.num}{a list of the number of membership functions for all inputs}

\item{input.mf.type}{designed for different membershp function types, however, currently, 'T1' for gbellmf, else 'it2gbellmf'}

\item{rule.num}{the number of rules}

\item{rule.which}{selected rules to be used in the full rule list, for example, c(1,2,3) specify the first three rules}

\item{defuzzMethod}{"default"}

\item{params.ante}{parameter settings for initialising antecedent membership functions}

\item{params.conse}{parameter settings for initialising consequent membership functions}
}
\description{
To build a one-output TSK FIS by automatically generating the input membership functions and the fuzzy rules
}
\author{
Chao Chen
}
