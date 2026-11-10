\name{prime_methods}
\alias{prime_methods}
\title{List supported audit acquisition methods}
\description{List supported audit acquisition methods.}
\usage{
prime_methods()
}
\value{A data frame describing 11 deterministic-support acquisition methods and four statistical confidence-sequence methods.}
\details{Deterministic describes certificate validity; some acquisition orders are randomized. Published methods are explicitly labelled task adaptations. The decision-cell statistical method chooses a public cell, then consumes its uniform random permutation without ranking individual units within the cell.}
\examples{
prime_methods()
}
