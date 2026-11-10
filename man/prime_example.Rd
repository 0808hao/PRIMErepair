\name{prime_example}
\alias{prime_example}
\title{Generate a small synthetic engineering example}
\description{Generate a small synthetic engineering example.}
\usage{
prime_example(n = 240L, seed = 42L)
}
\arguments{
\item{n}{Integer number of synthetic units, from four to one million.}
\item{seed}{Integer seed in [0, .Machine$integer.max]. The caller's RNG kind and state are restored, including after errors.}
}
\value{A list with pool, y and an explicit synthetic-data description.}
\details{No real patient or sleep-subject data are included.}
\examples{
ex <- prime_example(40)
print(ex$pool)
}
