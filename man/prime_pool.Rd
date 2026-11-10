\name{prime_pool}
\alias{prime_pool}
\title{Construct a locked finite prediction pool}
\description{Construct a locked finite prediction pool.}
\usage{
prime_pool(predictions, cost, reference, groups = NULL, strata = NULL,
  proxy = NULL, ids = NULL, cost_unit = "input channels")
}
\arguments{
\item{predictions}{Numeric matrix of finite probabilities, units by candidates, with unique nonempty column names. At least two candidates. Row order must match all other unit-level inputs.}
\item{cost}{Finite nonnegative cost vector, one per candidate. Named vectors are matched by candidate name; unnamed vectors follow input column order.}
\item{reference}{Name of one reference column. Its cost must be maximal. Alternatives more expensive than this reference are outside this version's contract.}
\item{groups}{Binary matrix of public, prespecified group membership, with named columns and no empty groups. Groups may overlap. Default: one All group. No minimum group-size filter is silently applied.}
\item{strata}{Optional unit-level vector for the deterministic age_stratified comparator. The historical name means supplied strata, not automatic age extraction. Named values match IDs. Statistical cell methods use group membership signatures instead.}
\item{proxy}{Optional locked probabilities used by the decision-focused heuristic and residual CS transformation. Default: reference probabilities. Literature comparators use their explicitly specified reference/ensemble surrogates.}
\item{ids}{Unique nonempty unit IDs. Default: prediction row names, otherwise sequential character IDs. Existing prediction/group row names must match IDs exactly in order.}
\item{cost_unit}{A nonempty label describing the cost scale. This does not validate economic or clinical meaning.}
}
\value{A \code{prime_pool} containing public prediction and bound arrays, candidate metadata and a content fingerprint; it contains no outcomes.}
\details{Candidates are stably sorted by cost, then input order, with the reference last among maximal-cost ties. Returned indices refer to this order. Input-column mapping is retained. The fingerprint detects accidental edits, not malicious tampering. Locking requires an external design in which audit labels have not trained the predictions or defined the groups. All units have equal label acquisition cost. Supports binary Brier loss only.}
\examples{
ex <- prime_example(40)
print(ex$pool)
}
