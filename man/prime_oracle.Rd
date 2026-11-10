\name{prime_oracle}
\alias{prime_oracle}
\title{Create a retrospective single-use label callback}
\description{Create a retrospective single-use label callback.}
\usage{
prime_oracle(y)
}
\arguments{
\item{y}{Nonempty numeric/logical binary vector, in pool row order. If named, names must be unique unit IDs matching the requested rows.}
}
\value{A stateful function(indices) that rejects repeat or invalid queries.}
\details{This callback is for retrospective experiments. It keeps outcomes outside the public pool; it is not a security sandbox. Real labeling services can implement the same callback contract. In case of service errors, reconcile any external charges before retrying. The returned function is not a reusable callback for independent runs.}
\examples{
oracle <- prime_oracle(c(0, 1, 0))
oracle(c(1, 3))
}
