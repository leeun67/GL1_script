         seed =  -1
      seqfile =/home/a2431559261/OLDISK/evo/genefami_12_28/pep2/mcmctree/cds.phylip
     treefile =/home/a2431559261/OLDISK/evo/genefami_12_28/pep2/mcmctree/pretime
      outfile =mcmc.out
        ndata =  1
      usedata =  1    * 0: no data; 1:seq like; 2:use in.BV; 3: out.BV
        clock =  3    * 1: global clock; 2: independent rates; 3: correlated rates
      RootAge =  <100.0  * safe constraint on root age, used if no fossil for root.
        model =  0    * 0:JC69, 1:K80, 2:F81, 3:F84, 4:HKY85
        alpha =  0.5    * alpha for gamma rates at sites
        ncatG =  4    * No. categories in discrete gamma
    cleandata =  0    * remove sites with ambiguity data (1:yes, 0:no)?
      BDparas =  1 1 0    * birth, death, sampling
  kappa_gamma =  6 2      * gamma prior for kappa
  alpha_gamma =  1 1      * gamma prior for alpha
  rgene_gamma =  2 2   * gamma prior for overall rates for genes
 sigma2_gamma =  1 10    * gamma prior for sigma^2     (for clock=2 or 3)
     finetune =  0.05 0.1 0.12 0.1 0.3  * times, rates, mixing, paras, RateParas
        print =  1
       burnin =  4000
     sampfreq =  10
      nsample =  20000

