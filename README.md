A repo to reproduce figures in  the snRNA-seq manuscript
-----------------------
## The codes are organized into two part:

 1. Reproduce figures of [Habib et
    al. 2017](https://www.ncbi.nlm.nih.gov/pubmed/28846088),
    aka. droncseq. Including a random forest classifier to classify GABA
    subtypes based on feature genes' expression. 
 2. Correlate our data with these two datasets


## Details 

* Folder [Habib2017](./Habib2017): contains codes to reproduce figures
  of [Habib et al 2017](https://www.ncbi.nlm.nih.gov/pubmed/28846088)
  *  [druncseq_gaba_dataprepare.R](./Habib2017/druncseq_gaba_dataprepare.R)
  and [habib2016.R](./Habib2017/habib2016.R) are codes to generate
  `druncseq_gaba` dataset and `snucseq.gaba.Rdata` to construct Random
  Forest classifiers ([RF_classifier.R](./Habib2017/RF_classifier.R)) based on variable genes 
  * [reproduce_mouse_clusering.R](./Habib2017/reproduce_mouse_clusering.R):
  code to take raw UMI matrix and perform clustering and feature
  selection. 
  * [reproduce_mouse_gaba_subclusering.R](./Habib2017/reproduce_mouse_gaba.R):code
  to subcluster GABA neurons. 
  * [reproduce_mouse.R](./Habib2017/reproduce_mouse.R) and
    [reproduce_human.R](./Habib2017/reproduce_human.R] are just codes to
    reproduce the figures. 


* Folder [figcode](./figcodes): contains two R scripts to produce correlation
  heatmaps between our data and dataset from
  [Habib2017](https://www.ncbi.nlm.nih.gov/pubmed/28846088). 
  * [fig1d.R](./figcodes/fig1d.R): Calculate and plot correlation between our mouse data with
      Habib2017's mouse data. Result is shown in [fig1d.pdf](./figs/fig1d.pdf)
  * [fig2e.R](./figcodes/fig2e.R): correlates our human data with
    Habib2017's human dat. Result is shown in [fig2e.pdf](./figs/fig2e.pdf)
	
  
