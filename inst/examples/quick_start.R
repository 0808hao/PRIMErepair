library(PRIMErepair)
ex <- prime_example(240, seed=42)
det <- prime_audit(ex$pool, prime_oracle(ex$y), epsilon=.002,
  method="decision_focused", budgets=c(.1,.5,1), batch=64)
print(det)
cs <- prime_audit(ex$pool, prime_oracle(ex$y), epsilon=.002,
  method="decision_cell_residual_wor_EB", alpha=.05, budgets=c(.1,.5,1), batch=64)
print(cs)
point <- prime_lure(ex$pool, prime_oracle(ex$y), budgets=c(.1,.5), seed=42)
print(point)
oracle_truth <- prime_truth(ex$pool,ex$y)
print(oracle_truth$candidates)
stopifnot(all(det$summary$certified_minimum_cost[det$summary$cap_fraction==1]))
