library(PRIMErepair)
n_checks <- 0L
ok <- function(x,label) {
  if(!isTRUE(x)) stop("FAIL: ",label,call.=FALSE)
  n_checks <<- n_checks+1L
}
equal <- function(x,y) isTRUE(all.equal(x,y,tolerance=1e-12,check.attributes=TRUE))
reject <- function(expr,label) ok(inherits(tryCatch(force(expr),error=identity),"error"),label)
internal <- function(nm) get(nm,asNamespace("PRIMErepair"),inherits=FALSE)

# Fixture comes from the original Step133 implementation, not from this API.
f <- readRDS(system.file("extdata","frozen_step133_fixture.rds",package="PRIMErepair"))
p <- f$public
pool <- prime_pool(p$P,p$costs,colnames(p$P)[p$K],p$G,p$strata,ids=p$ids)
ok(equal(pool$public,p),"public arrays match frozen Step133")
for(m in names(f$det)) {
  a <- prime_audit(pool,prime_oracle(f$y),epsilon=.002,method=m,seed=19,
    budgets=c(.1,.5,1),batch=7)
  s <- a$summary; names(s)[names(s)=="selected_cost"] <- "selected_input_channels"
  ok(equal(s[,names(f$det[[m]]$snapshots)],f$det[[m]]$snapshots),paste(m,"frozen decisions"))
  ok(equal(a$bounds,f$det[[m]]$bounds),paste(m,"frozen bounds"))
  ok(identical(a$queries$index,f$det[[m]]$queried_indices),paste(m,"frozen query path"))
  ok(!anyDuplicated(a$queries$index),paste(m,"unique charging"))
}
for(m in names(f$stat)) {
  a <- prime_audit(pool,prime_oracle(f$y),epsilon=c(.001,.002,.005),method=m,
    seed=19,alpha=c(.05,.01,.001),budgets=c(.1,.5,1),batch=7)
  s <- a$summary; names(s)[names(s)=="selected_cost"] <- "selected_input_channels"
  ok(equal(s[,names(f$stat[[m]]$runs)],f$stat[[m]]$runs),paste(m,"frozen decisions"))
  ok(equal(a$bounds,f$stat[[m]]$final_bounds),paste(m,"frozen bounds"))
  ok(identical(a$queries$index,f$stat[[m]]$queried_indices),paste(m,"frozen query path"))
  ok(!anyDuplicated(a$queries$index),paste(m,"unique charging"))
  if(m=="decision_cell_residual_wor_EB") {
    ss <- internal("cs_setup")(pool$public,m,c(.05,.01,.001),19)
    for(h in seq_len(ncol(ss$R))) {
      actual<-a$queries$index[ss$R[a$queries$index,h]==1]
      expected<-ss$order[ss$R[ss$order,h]==1]
      ok(identical(as.integer(actual),as.integer(head(expected,length(actual)))),"uniform within-cell prefixes")
    }
  }
}
for(m in names(f$lure)) {
  a <- prime_lure(pool,prime_oracle(f$y),method=m,seed=19,
    epsilon=c(.001,.002,.005),budgets=c(.1,.5,1))
  s<-a$summary;names(s)[names(s)=="selected_cost"]<-"selected_input_channels"
  ok(equal(s[,names(f$lure[[m]]$decisions)],f$lure[[m]]$decisions),paste(m,"frozen LURE choices"))
  ok(equal(a$estimates,f$lure[[m]]$estimates),paste(m,"frozen LURE estimates"))
  ok(all(a$summary$certificate=="NONE_POINT_ESTIMATE_ONLY"),"LURE not labelled as a certificate")
}

# Exhaust all possible binary labels for a small pool, with direct full-label truth.
P <- cbind(low=c(.1,.8,.2,.9,.4),mid=c(.2,.7,.3,.8,.5),ref=c(.4,.5,.4,.6,.5))
G <- cbind(All=1,small=c(1,1,0,0,0),overlap=c(0,1,1,1,0))
tiny <- prime_pool(P,c(1,2,3),"ref",G)
Ys <- as.matrix(expand.grid(rep(list(0:1),5)))
for(j in seq_len(nrow(Ys))) {
  y <- as.numeric(Ys[j,]); tr<-prime_truth(tiny,y,.02)
  d <- sweep((P-y)^2,1,(P[,"ref"]-y)^2,"-")
  exact <- sweep(crossprod(d,G),2,colSums(G),"/")
  ok(equal(unname(tr$group_gaps$mean_brier_gap),as.vector(exact)),"independent Brier truth")
  for(m in c("decision_focused","shared_width","portfolio_replay")) {
    z<-prime_audit(tiny,prime_oracle(y),epsilon=.02,method=m,batch=2,budgets=c(.4,1))
    end<-tail(z$summary,1); b<-z$bounds[z$bounds$cap_labels==5,]
    ok(end$certified_minimum_cost && end$selected_cost==tr$minimum_cost &&
      tr$candidates$feasible[match(end$selected_candidate,tr$candidates$candidate)],"exhaustive correct cost certificate")
    ok(all(b$lower<=as.vector(exact)+1e-10 & b$upper>=as.vector(exact)-1e-10),"deterministic bounds contain all completions")
  }
}

# Independent scalar recurrence, followed by stream-vs-batch checks.
manual <- function(x,lo,hi,alpha,J) {
  N<-length(x);c0<-hi-lo; sx<-sse<-wn<-wd<-sn<-sd<-pen<-0;L<-lo;U<-hi;ans<-matrix(NA_real_,N,2)
  for(t in seq_len(N)) {
    center<-if(t==1)(lo+hi)/2 else wn/wd
    lambda<-min(sqrt(2*log(2*J/alpha)/(((c0^2/4+sse)/t)*t*log(t+1))),.5/c0)
    corr<-x[t]+sx/(N-t+1);den<-N/(N-t+1)
    sn<-sn+lambda*corr;sd<-sd+lambda*den
    pen<-pen+(x[t]-center)^2*(-log1p(-c0*lambda)-c0*lambda)/c0^2
    L<-max(L,(sn-pen-log(2*J/alpha))/sd);U<-min(U,(sn+pen+log(2*J/alpha))/sd)
    sx<-sx+x[t];sse<-sse+(x[t]-sx/t)^2;wn<-wn+corr;wd<-wd+den
    if(t==N){L<-max(L,mean(x));U<-min(U,mean(x))};ans[t,]<-c(L,U)
  };ans
}
x<-c(.1,.9,.2,.7,.8,.3);expected<-manual(x,0,1,.05,3)
state<-internal("eb_state")(length(x),0,1,0,.05,3)
for(i in seq_along(x)) {
  state<-internal("eb_update")(state,matrix(x[i],1,1))
  ok(max(abs(c(state$L,state$U)-expected[i,]))<1e-12,"scalar EB recurrence")
}
single<-internal("eb_update")(internal("eb_state")(length(x),0,1,0,.05,3),matrix(x,ncol=1))
ok(isTRUE(all.equal(state,single,tolerance=1e-12,check.attributes=FALSE)),"stream-batch EB numerical identity")
state<-internal("eb_state")(2,0,1,0,.05,1)
state<-internal("eb_update")(state,matrix(0,1,1));state$L[,]<-.5
state<-internal("eb_update")(state,matrix(0,1,1))
ok(state$L[1,1]>state$U[1,1],"exhaustion preserves historical coverage failure")

# Enumerate randomization distribution: two populations, all 120 permutations.
perms<-function(x) if(length(x)==1L) list(x) else unlist(lapply(seq_along(x),
  function(j) lapply(perms(x[-j]),function(y)c(x[j],y))),recursive=FALSE)
orders<-perms(1:5)
for(x in list(c(0,0,0,0,1),c(0,.1,.3,.6,1))) for(a in c(.05,.2,.5)) {
  bad<-vapply(orders,function(o) {
    s<-internal("eb_state")(5,0,1,0,a,1);fail<-FALSE
    for(i in o) {s<-internal("eb_update")(s,matrix(x[i],1,1));
      fail<-fail||s$L[1,1]>mean(x)+1e-12||s$U[1,1]<mean(x)-1e-12}
    fail
  },TRUE)
  ok(mean(bad)<=a+1e-12,"exact finite randomization coverage check")
}

# Boundary and public API behavior.
same<-prime_pool(cbind(a=c(0,1,.5),ref=c(0,1,.5)),c(0,1),"ref")
for(m in prime_methods()$method) {
  a<-prime_audit(same,function(i)stop("Unneeded label queried"),method=m,batch=100)
  ok(nrow(a$queries)==0L && all(a$summary$certified_minimum_cost),paste(m,"zero labels for identical predictions"))
}
one<-prime_pool(matrix(c(.2,.6),1,dimnames=list(NULL,c("a","ref"))),c(0,1),"ref")
for(m in prime_methods()$method) {
  a<-prime_audit(one,prime_oracle(1),method=m,batch=100)
  ok(all(a$summary$certified_minimum_cost),paste(m,"one-unit population"))
}
ex<-prime_example(40)
set.seed(123);old<-.Random.seed;kinds<-RNGkind()
invisible(prime_audit(ex$pool,prime_oracle(ex$y),method="srs"))
ok(identical(old,.Random.seed)&&identical(kinds,RNGkind()),"RNG restored after success")
err<-tryCatch(prime_audit(ex$pool,function(i)stop("service failed"),method="srs"),error=identity)
ok(inherits(err,"prime_query_error")&&length(err$pending_indices)>0&&nrow(err$completed_queries)==0,"service error retains accounting")
ok(identical(old,.Random.seed)&&identical(kinds,RNGkind()),"RNG restored after failure")
rm(".Random.seed",envir=.GlobalEnv)
invisible(prime_audit(ex$pool,prime_oracle(ex$y),method="srs"))
ok(!exists(".Random.seed",envir=.GlobalEnv,inherits=FALSE),"absence of RNG state restored")
ok(equal(prime_truth(ex$pool,ex$y),prime_truth(ex$pool,rev(ex$y))),"named full labels align by ID")
reordered<-prime_pool(P[,c(3,1,2)],c(ref=3,mid=2,low=1),"ref",G)
ok(equal(reordered$public$P,tiny$public$P),"candidate names align with costs")
changed<-ex$pool;changed$public$P[1,1]<-0
reject(prime_audit(changed,prime_oracle(ex$y)),"edited pool fingerprint rejected")
reject(prime_pool(P,c(1,2,-1),"ref"),"negative cost rejected")
reject(prime_pool(P,c(1,4,3),"ref"),"reference not maximal rejected")
reject(prime_pool(P,c(low=1,wrong=2,ref=3),"ref"),"misnamed cost rejected")
reject(prime_pool(P,c(1,2,3),"ref",groups=cbind(All=1,empty=rep(0,5))),"empty group rejected")
reject(prime_pool(P,c(1,2,3),"ref",ids=rep("same",5)),"duplicate IDs rejected")
bad<-P;bad[1,1]<-NA
reject(prime_pool(bad,c(1,2,3),"ref"),"missing prediction rejected")
bad<-P;bad[1,1]<-1.01
reject(prime_pool(bad,c(1,2,3),"ref"),"invalid probability rejected")
bad<-P;colnames(bad)<-c("a","a","ref")
reject(prime_pool(bad,c(1,2,3),"ref"),"duplicate candidates rejected")
for(b in list(0,-1,NA,Inf,1.1,c(.1,.1)))
  reject(prime_audit(ex$pool,prime_oracle(ex$y),budgets=b),"invalid budgets rejected")
for(b in list(0,-1,NA,Inf,1.1))
  reject(prime_audit(ex$pool,prime_oracle(ex$y),batch=b),"invalid batch rejected")
for(a in list(0,1,NA,Inf,c(.05,.05)))
  reject(prime_audit(ex$pool,prime_oracle(ex$y),method="srs_wor_EB",alpha=a),"invalid alpha rejected")
reject(prime_audit(ex$pool,prime_oracle(ex$y),alpha=.05),"alpha cannot silently alter deterministic mode")
reject(prime_audit(ex$pool,prime_oracle(ex$y),epsilon=-.001),"negative margin rejected")
reject(prime_audit(ex$pool,prime_oracle(ex$y),method="typo"),"unknown method rejected")
reject(prime_audit(ex$pool,function(i)rep(2,length(i))),"nonbinary returned labels rejected")
reject(prime_audit(ex$pool,function(i)setNames(rep(0,length(i)),rep("wrong",length(i)))),"misaligned returned labels rejected")
oracle<-prime_oracle(c(0,1,0));invisible(oracle(1))
reject(oracle(1),"oracle prevents duplicate charging")
reject(prime_lure(ex$pool,prime_oracle(ex$y),method="decision_focused"),"LURE method restrictions")

# Complete-pool LURE is exact; SRS weights are exactly one at any fixed budget.
for(M in c(1,4,10)) ok(max(abs(internal("lure_weights")(10,M,1/(10-seq_len(M)+1))-1))<1e-12,"SRS LURE unit weights")
l<-prime_lure(tiny,prime_oracle(as.numeric(Ys[1,])),budgets=1)
t<-prime_truth(tiny,as.numeric(Ys[1,]))
ok(equal(l$estimates$estimate,t$group_gaps$mean_brier_gap),"full-pool LURE equals truth")

# Export contract, no overwriting, and round-trip persistence.
a<-prime_audit(ex$pool,prime_oracle(ex$y),budgets=c(.1,1),batch=7)
out<-tempfile("prime_export_test_");prime_export(a,out)
ok(identical(readRDS(file.path(out,"audit_result.rds")),a),"export round trip")
ok(all(dir.exists(file.path(out,c("01_methodology_tables","02_result_tables")))),"separate table folders")
manifest<-read.csv(file.path(out,"file_manifest.csv"))
ok(identical(unname(tools::md5sum(file.path(out,manifest$file))),manifest$md5),"export file hashes")
reject(prime_export(a,out),"existing export never overwritten")
unlink(out,recursive=TRUE)
out<-tempfile("prime_lure_export_test_");prime_export(l,out)
ok(file.exists(file.path(out,"02_result_tables","lure_estimates.csv")),"LURE export")
unlink(out,recursive=TRUE)
cat("PASS:",n_checks,"checks, including frozen-path regression and independent finite-pool checks.\n")
