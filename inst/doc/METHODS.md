# PRIMErepair 0.1.0: method and scope contract

## Target and design

Fix a finite pool of N units, K candidate probability vectors, a reference r,
public binary group membership, nonnegative input costs, and binary outcomes.
Predictions, candidate library, costs, proxy, groups and margins must be locked
before audit labels are queried. Any fitting, calibration or hyperparameter
selection must use separate development information. The software cannot verify
feature availability or enforce this separation outside its own query API.

For candidate k, unit i and unknown label y in {0,1}, define

    d_i,k(y) = (p_i,k - y)^2 - (p_i,r - y)^2
    lo_i,k = min(d_i,k(0), d_i,k(1))
    hi_i,k = max(d_i,k(0), d_i,k(1)).

The target Delta_k,g is the mean d_i,k(y_i) over the known finite set of units
in group g. Positive means worse than reference. The optimization problem is
minimum declared input cost subject to Delta_k,g <= epsilon for every included
group. The reference self-gap is zero. Version 0.1.0 requires the reference to
have maximal cost and uses one scalar margin shared by groups.

Unit label acquisition costs are equal; a queried unit reveals its paired loss
for every candidate and every group containing that unit. Input cost and label
cost are distinct objectives. Costs in the development studies are channel or
feature-family counts, not actual clinical or financial costs.

## Deterministic certificates

For queried set Q, the lower and upper bounds replace each queried unit's support
endpoint with its observed paired loss. Divide summed endpoints by known group
size. These bounds hold for every possible completion of unqueried binary labels,
independently of the acquisition order, proxy accuracy and stopping time.

A candidate is feasible when all upper bounds are <= epsilon + 1e-12. It is
excluded when any lower bound is > epsilon + 1e-12. A feasible candidate has
certified minimum cost once all strictly cheaper candidates are excluded. This
does not identify a uniquely best model or rank every equal-cost candidate.
Among currently certified feasible ties, stable input order is used, with the
reference last. Fixed-preference certification is reported separately.

The decision-focused policy retains Step131's heuristic. A fixed proxy estimates
expected paired loss. If a candidate appears to violate a constraint, query scores
favor raising the lower bound of the most promising exclusion constraint;
otherwise they favor reducing upper bounds of unresolved groups. The heuristic
may be wrong and need not beat width-based or random allocation. It affects
efficiency, not deterministic validity. The proxy-completed mean is not an
unbiased estimate or confidence interval.

Shared-width sums unresolved support-width contributions. Portfolio replay
alternates independent virtual decision-focused and shared-width batches while
paying for a shared unit only once. Each lane consumes only its own history;
cross-lane observations are not injected into its acquisition heuristic early.
The two-lane bound is min(N, 2 B ceil(min(L_decision,L_width)/B)) when L values
refer to complete independent B-batch anchor paths on the same pool/protocol.
Intermediate caps can truncate a pending batch and are part of the protocol.
No universal efficiency improvement is asserted.

## Separate without-replacement statistical mode

This mode implements a task adaptation of the empirical-Bernstein confidence
sequence in Waudby-Smith and Ramdas (2020), Theorem 3.2 and Eqs. 3.2, 3.10, 3.13.
For each region, sample a uniform permutation without replacement. With region
size N_h, known interval [l,u], c=u-l, and local observation number j, set

    A_j = X_j + sum_{i<j} X_i / (N_h-j+1)
    D_j = N_h / (N_h-j+1)
    muhat_WoR,j-1 = sum_{i<j} A_i / sum_{i<j} D_i
    vhat_j-1 = (c^2/4 + sum_{i<j}(X_i - mean(X_1:i))^2) / j
    ell = log(2 J / alpha)
    lambda_j = min(sqrt(2 ell / (vhat_j-1 j log(j+1))), 0.5/c)
    psi_j = [-log(1-c lambda_j) - c lambda_j] / c^2
    V_j = sum_{i<=j} psi_i (X_i - muhat_WoR,i-1)^2.

At j=1, use midpoint (l+u)/2 as the predictable center. Take the running
intersection, at every local sample, of

    [(sum lambda_i A_i - V_j - ell) / sum lambda_i D_i,
     (sum lambda_i A_i + V_j + ell) / sum lambda_i D_i].

J is the number of nonconstant candidate-region constraints. Bonferroni covers
the two tails and these J constraints. Constant regions are exact. The
implementation uses at least 1e-14 for c in numerical divisions, which only
enlarges the supplied range when needed. At exhaustion, intersect with the exact
mean; do not erase earlier failed intervals. The frozen scalar recurrence and
batch recurrence are independently compared in the regression suite.

Residual modes use X_i = d_i(y_i) - pr_i, where pr_i is a public fixed-proxy
expected loss gap. The finite-region mean of pr_i is known and added back.
Residual supports are computed from public endpoint arrays, not from sampled
outcomes.

SRS modes use the uniform global order; its restriction to each prespecified
group is uniform without replacement. Cell modes partition by complete public
group-membership signatures, maintain candidate-by-cell sequences, and combine
cell means using known group fractions. The decision-cell rule may adapt which
cell is queried based on observed labels, but consumes a uniform permutation
within each cell. Simultaneous coverage of every cell prefix supports this
adaptive interleaving; individual within-cell score ranking is not allowed.

Statistical bounds are intersected with deterministic support bounds. If a
conflict larger than 1e-10 occurs, it is recorded and that alpha uses deterministic
bounds subsequently. Small rounding inversions are reconciled by their midpoint.
The procedure does not observe full truth and cannot detect all coverage errors.
Certificates issued before a conflict remain in the record. Benchmark coverage
observers use full truth only outside acquisition.

For each prespecified alpha separately, coverage is simultaneous across the
locked candidate/group constraints and time in one fixed pool. Epsilon values
can share the coverage event. Multiple alpha values are distinct error regimes,
not a license to choose a favorable level after seeing results. Multiple
independent or overlapping audits do not automatically have joint coverage.
The same coverage statement does not certify future populations, hospital-level
generalization, causal validity or clinical safety.

## Literature acquisition and fixed-budget LURE

Hara ensemble uses Eq.11's surrogate-expected absolute pairwise Brier gaps,
summed over candidate pairs, with the mean candidate prediction as surrogate.
The reference-surrogate variant is a controlled adaptation. The Active Testing
reference comparator uses expected reference Brier loss. Static positive weights
produce a weighted-without-replacement order through exponential races.
The query probability q_m is conditional on remaining weights, not the initial
weight or the final inclusion probability. A fixed floor of 1e-10 times mean
scaled weight ensures positivity; changed-unit counts are retained.

These acquisition policies can use the SAME deterministic support stopping rule
as PRIME. Separately, prime_lure computes fixed-budget estimates with

    v_m = 1 + (N-M)/(N-m) * [1 / ((N-m+1) q_m) - 1],

and v_m=1 at M=N. Mean weighted loss estimates overall risk. Group gaps use
the known group size: N/(M n_g) times the weighted sum of group-member gaps.
No clipping or weight normalization is applied. Constrained plugin selection
and minimum overall estimated-risk selection are reported separately. Neither
is a confidence certificate; M must not be chosen by outcome-dependent stopping
if interpreting the estimator under the fixed-budget unbiasedness argument.
Monte Carlo MAE means the mean of per-run absolute errors, not the absolute
error of the mean across repetitions.

## Evidence retained with the release

Step133 aggregate CSVs are an immutable historical evidence record. They do not
contain patient records. Its original readiness table deliberately still shows
the then-pending package checks; the release bundle has a separate current
engineering report. No historical table is rewritten to imply it had already
passed packaging at the time.

The eICU primary deterministic decision-focused mean is 3,494.4 labels versus
5,056 shared-width and 6,854.4 Hara ensemble under the same support criterion.
Five reused hospital splits are not five independent external cohorts. The
Sleep task used 80 subjects across development/audit/test, with 16,000 audit
epochs. Deterministic decision-focused and shared-width both used 384 labels;
the decision-cell alpha=0.05 mode averaged 272. Both retained two feature
families. This is development-stage cross-domain replication. The prior strong
Hara+LURE point-selection result at 2.5% budget remains valid and is not a
certificate comparison. This package does not claim uniform dominance.

Six calibration scenarios times four methods times three alpha values had zero
observed coverage failures in 200 orders per combination. With zero failures,
the one-sided 95% binomial upper bound is 0.014867, so these simulations alone
cannot establish nominal error 0.01 or 0.001. Derivation, implementation tests,
empirical performance and clinical validation are separate forms of evidence.

## References

- Waudby-Smith I, Ramdas A. Confidence sequences for sampling without replacement.
  NeurIPS 2020. https://proceedings.neurips.cc/paper/2020/hash/e96c7de8f6390b1e6c71556e4e0a4959-Abstract.html
- Hara S, Matsuura M, Honda J, Ito S. Active model selection: A variance
  minimization approach. Machine Learning 113:8327-8345 (2024).
  https://doi.org/10.1007/s10994-024-06603-1
- Kossen J, Farquhar S, Gal Y, Rainforth T. Active Testing: Sample-Efficient Model
  Evaluation. PMLR 139:5753-5763 (2021).
  https://proceedings.mlr.press/v139/kossen21a.html

These references supply comparator/CS foundations. PRIME's candidate-cost
objective, public-cell adaptation and implementation details must be described
as such, not attributed as original theorems of the cited authors.
