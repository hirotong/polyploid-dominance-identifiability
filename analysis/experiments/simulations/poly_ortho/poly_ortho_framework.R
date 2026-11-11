#!/usr/bin/env Rscript
# =====================================================================
# Polyploid orthogonal dominance framework + Figure-1 sweep (idea:010).
# (A) METHOD: per-locus orthonormal dosage contrasts (NOIA generalised to
#     autotetraploid, degrees 1..4 = additive/digenic/trigenic/quadrigenic),
#     assembled into orthogonal relationship kernels A^,D2^,D3^,D4^.
# (B) SWEEP: simulate pure-degree-k truths on a structured autotetraploid
#     population (R reps); compare
#       STANDARD model  A_VanRaden + D_digenic   -> Vd share
#       ORTHOGONAL model A^+D2^+D3^+D4^           -> variance fraction per degree
#     Figure 1A: standard Vd-share by true degree (flat ~0.5 => non-identifiable).
#     Figure 1B: orthogonal attribution matrix (near-diagonal => recovers truth).
# (C) ALTERNATIVE KERNELS (supplement, table tab:robust): on the SAME replicate
#     traits, the AGHmatrix Slater (full-autopolyploid) and Endelman (digenic)
#     kernels are fitted in place of the D kernel, alongside the same additive
#     kernel. Sharing the replicates keeps the own-kernel column of that table
#     identical to Figure 1A. (AGHmatrix's VanRaden additive kernel equals Amat()
#     to machine precision after normalization, so one additive kernel is used.)
# =====================================================================
suppressMessages({library(AlphaSimR); library(BGLR); library(AGHmatrix); library(parallel)})
set.seed(123)
OUT <- if(length(commandArgs(TRUE))>=1) commandArgs(TRUE)[1] else "out"; dir.create(OUT,showWarnings=FALSE,recursive=TRUE)
PL<-4L; H2<-0.6; REPS<-as.integer(Sys.getenv("REPS","50")); NITER<-4000L; BURN<-1000L
NC<-as.integer(Sys.getenv("NC", Sys.getenv("SLURM_CPUS_PER_TASK","8")))
TRUTHS<-c("additive","digenic","trigenic","quadrigenic","hetero")
STRUCTURES<-c("unrelated","halfsib","fullsib")  # n~400 each, from ONE shared founder pool

## ---- (A) orthonormal dosage contrasts under binomial genotype freqs ----
ortho <- function(p, ploidy=PL){
  d<-0:ploidy; w<-choose(ploidy,d)*p^d*(1-p)^(ploidy-d); w<-w/sum(w)
  inner<-function(a,b) sum(w*a*b); Q<-matrix(0,ploidy+1,ploidy+1)
  for(k in 1:(ploidy+1)){ v<-d^(k-1)
    if(k>1) for(j in 1:(k-1)) v<-v-inner(v,Q[,j])/inner(Q[,j],Q[,j])*Q[,j]
    nrm<-sqrt(inner(v,v)); Q[,k]<-if(nrm>0) v/nrm else v }
  Q }  # cols 1..5 = degree 0..4 ; rows d=0..4

## per-marker contrast covariate matrices (degrees 1..4) -> orthogonal kernels
build_ortho_kernels <- function(M, ploidy=PL){
  n<-nrow(M); poly<-which(apply(M,2,function(g){p<-mean(round(g))/ploidy; p>0 && p<1}))
  X<-list(matrix(0,n,length(poly)),matrix(0,n,length(poly)),matrix(0,n,length(poly)),matrix(0,n,length(poly)))
  for(idx in seq_along(poly)){ j<-poly[idx]; g<-round(M[,j]); p<-mean(g)/ploidy; Q<-ortho(p,ploidy)
    for(deg in 1:4) X[[deg]][,idx]<-Q[g+1, deg+1] }
  mkK<-function(Z){ K<-tcrossprod(Z)/ncol(Z); K/mean(diag(K))+diag(1e-6,nrow(Z)) }
  list(A=mkK(X[[1]]), D2=mkK(X[[2]]), D3=mkK(X[[3]]), D4=mkK(X[[4]])) }

## standard kernels
Amat<-function(M,pl=PL){ p<-colMeans(M)/pl; Z<-sweep(M,2,pl*p); d<-sum(pl*p*(1-p)); if(d<=0)d<-1
  K<-tcrossprod(Z)/d; K/mean(diag(K))+diag(1e-6,nrow(M)) }
Dmat<-function(M,pl=PL){ Hh<-(M*(pl-M))/choose(pl,2); Hh<-sweep(Hh,2,colMeans(Hh))
  K<-tcrossprod(Hh); K/mean(diag(K))+diag(1e-6,nrow(M)) }

## ---- ONE shared coalescent founder pool -> three population structures ----
## Single runMacs + single SnpChip => all three structures share identical marker
## loci, so any difference between structures is RELATEDNESS, not allele freq/LD.
set.seed(123)
N0<-450L
fp<-runMacs(nInd=N0,nChr=8,segSites=600,ploidy=PL,species="GENERIC")
SP<-SimParam$new(fp); SP$addSnpChip(nSnpPerChr=300L)
founders<-newPop(fp,simParam=SP)
pop_unrel<-founders[1:400]                                              # unrelated: relatedness ~0
pop_fs   <-randCross(founders[1:40], nCrosses=50, nProgeny=8, simParam=SP)  # full-sib: 50 fam x8, within ~0.5
## paternal half-sib: sire f x 8 distinct dams, 1 progeny each (within-fam ~0.25, between ~0)
sires<-1:50; damsM<-matrix(51:450, nrow=50, byrow=TRUE)
plan<-do.call(rbind, lapply(1:50, function(f) cbind(damsM[f,], rep(sires[f],8))))  # col1=dam, col2=sire(shared)
pop_hs<-makeCross(founders, crossPlan=plan, nProgeny=1L, simParam=SP)
POPS<-list(unrelated=pullSnpGeno(pop_unrel,simParam=SP),
           halfsib  =pullSnpGeno(pop_hs,   simParam=SP),
           fullsib  =pullSnpGeno(pop_fs,   simParam=SP))[STRUCTURES]
## SAME causal loci for every structure (drawn from markers polymorphic in ALL three)
polyOf<-function(M) which(apply(M,2,function(g){p<-mean(round(g))/PL; p>0 && p<1}))
common<-Reduce(intersect, lapply(POPS, polyOf))
set.seed(7); causal<-sort(sample(common,200))
cat(sprintf("founder pool n=%d; %d markers; %d polymorphic in all structures; 200 causal\n",
    N0, ncol(POPS[[1]]), length(common)))

build_gv<-function(M,causal,truth){
  gv<-numeric(nrow(M))
  if(truth=="hetero"){ Hh<-(M[,causal,drop=FALSE]*(PL-M[,causal,drop=FALSE]))/choose(PL,2); return(scale(rowMeans(Hh))[,1]) }
  deg<-switch(truth,additive=1L,digenic=2L,trigenic=3L,quadrigenic=4L)
  for(j in causal){ g<-round(M[,j]); p<-mean(g)/PL; if(p<=0||p>=1) next; Q<-ortho(p,PL); gv<-gv+rnorm(1)*Q[g+1,deg+1] }
  gv }

## prebuild kernels per structure (depend only on M); cache for the parallel sweep
vc<-function(fm,nm) fm$ETA[[nm]]$varU
cat("building kernels per structure...\n")
normK<-function(K){ K/mean(diag(K))+diag(1e-6,nrow(K)) }
aghK<-function(SNP,method,s) tryCatch(
  normK(AGHmatrix::Gmatrix(SNP,method=method,ploidy=PL,missingValue=-9,integer=TRUE)),
  error=function(e){ cat(method,"failed [",s,"]:",conditionMessage(e),"\n"); NULL })
PREP<-lapply(names(POPS), function(s){ M<-POPS[[s]]; SNP<-M; mode(SNP)<-"integer"
  list(M=M, n=nrow(M), nmark=length(polyOf(M)),
       KA=Amat(M), KD=Dmat(M), KO=build_ortho_kernels(M),
       KDsl=aghK(SNP,"Slater",s), KDen=aghK(SNP,"Endelman",s)) })
names(PREP)<-names(POPS)
## relatedness sanity: within- vs between-family off-diagonal of the additive kernel
## (expect fullsib within > halfsib within > between ~ unrelated)
fam_of<-function(s,n) if(s=="unrelated") seq_len(n) else ((seq_len(n)-1L)%/%8L + 1L)
for(s in names(PREP)){ K<-PREP[[s]]$KA; n<-PREP[[s]]$n; fam<-fam_of(s,n)
  od<-row(K)!=col(K); same<-outer(fam,fam,`==`)&od
  wi<-if(any(same)) mean(K[same]) else NA_real_; bw<-mean(K[od & !same])
  cat(sprintf("  [%-9s] n=%d valid_markers=%d  KA within-fam=%.3f between-fam=%.3f\n",
      s, n, PREP[[s]]$nmark, wi, bw)) }

## ---- structure x truth x rep sweep, parallel over jobs (mclapply) ----
## Each job self-seeds from (rep,truth) exactly as before, so results are
## independent of execution order; structure just selects the cached kernels.
jobs<-expand.grid(structure=names(POPS), truth=TRUTHS, rep=1:REPS, stringsAsFactors=FALSE)
run_one<-function(i){
  st<-jobs$structure[i]; truth<-jobs$truth[i]; r<-jobs$rep[i]; P<-PREP[[st]]; n<-P$n
  set.seed(1000*r+which(TRUTHS==truth)); gv<-scale(build_gv(P$M,causal,truth))[,1]
  ve<-var(gv)*(1-H2)/H2; y<-gv+rnorm(n,0,sqrt(ve))
  pfx<-file.path(tempdir(),sprintf("j_%s_%s_%d_",st,truth,r))  # unique per fork: no BGLR temp collisions
  # STANDARD A + digenic D
  f1<-BGLR(y=y,ETA=list(A=list(K=P$KA,model="RKHS"),D=list(K=P$KD,model="RKHS")),nIter=NITER,burnIn=BURN,verbose=FALSE,saveAt=paste0(pfx,"s_"))
  vA<-vc(f1,"A"); vD<-vc(f1,"D")
  # ORTHOGONAL A^+D2^+D3^+D4^
  f2<-BGLR(y=y,ETA=list(A=list(K=P$KO$A,model="RKHS"),D2=list(K=P$KO$D2,model="RKHS"),D3=list(K=P$KO$D3,model="RKHS"),D4=list(K=P$KO$D4,model="RKHS")),nIter=NITER,burnIn=BURN,verbose=FALSE,saveAt=paste0(pfx,"o_"))
  v<-c(vc(f2,"A"),vc(f2,"D2"),vc(f2,"D3"),vc(f2,"D4")); v<-v/sum(v)
  # ALTERNATIVE dominance kernels on the same trait, same additive kernel (fitted
  # after f1 and f2, so those two fits see the same random-number stream as before)
  alt_share<-function(KD,tag){ if(is.null(KD)) return(NA_real_)
    fm<-BGLR(y=y,ETA=list(A=list(K=P$KA,model="RKHS"),D=list(K=KD,model="RKHS")),nIter=NITER,burnIn=BURN,verbose=FALSE,saveAt=paste0(pfx,tag))
    vc(fm,"D")/(vc(fm,"A")+vc(fm,"D")) }
  sl<-alt_share(P$KDsl,"sl_"); en<-alt_share(P$KDen,"en_")
  list(std =data.frame(structure=st,truth=truth,rep=r,Vd_share=vD/(vA+vD)),
       orth=data.frame(structure=st,truth=truth,rep=r,fA=v[1],fD2=v[2],fD3=v[3],fD4=v[4]),
       alt =data.frame(structure=st,truth=truth,rep=r,Vd_own=vD/(vA+vD),Vd_slater=sl,Vd_endelman=en))
}
cat(sprintf("sweep: %d jobs (%d structures x %d truths x %d reps), NC=%d\n",
    nrow(jobs), length(STRUCTURES), length(TRUTHS), REPS, NC))
resl<-mclapply(seq_len(nrow(jobs)), run_one, mc.cores=NC)
std <-do.call(rbind, lapply(resl, `[[`, "std"))
orth<-do.call(rbind, lapply(resl, `[[`, "orth"))
write.csv(std ,file.path(OUT,"std_vd_share.csv")    ,row.names=FALSE)  # cols: structure,truth,rep,Vd_share
write.csv(orth,file.path(OUT,"ortho_attribution.csv"),row.names=FALSE) # cols: structure,truth,rep,fA,fD2,fD3,fD4
alt <-do.call(rbind, lapply(resl, `[[`, "alt"))
write.csv(alt ,file.path(OUT,"alt_kernel_share.csv") ,row.names=FALSE) # cols: structure,truth,rep,Vd_own,Vd_slater,Vd_endelman

## per-structure x truth summary to the log
for(st in names(POPS)){ for(truth in TRUTHS){
  s<-std$Vd_share[std$structure==st & std$truth==truth]
  o<-orth[orth$structure==st & orth$truth==truth,]
  cat(sprintf("  [%-9s] truth=%-11s standard Vd_share=%.2f | ortho A/D2/D3/D4 = %.2f/%.2f/%.2f/%.2f\n",
      st, truth, mean(s), mean(o$fA), mean(o$fD2), mean(o$fD3), mean(o$fD4))) } }
for(st in names(POPS)){ for(truth in TRUTHS){ a<-alt[alt$structure==st & alt$truth==truth,]
  cat(sprintf("  [%-9s] truth=%-11s D-kernel share: own=%.2f  Slater=%.2f  Endelman=%.2f\n",
      st, truth, mean(a$Vd_own), mean(a$Vd_slater,na.rm=TRUE), mean(a$Vd_endelman,na.rm=TRUE))) } }

## Figure 1 (paper/figures/fig_simulation.png) is built by make_fig_simulation.R from these CSVs.
cat("DONE\n")
