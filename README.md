# Spas

We are often interested in the association between the expression level of a set of genes in a cell or tissue and a measured physical trait or phenotype. We may have prior knowledge of which genes likely cause a certain phenotype, and therefore we may want to predict what phenotype a sample might have given its gene expression. Alternatively, we may want to determine whether a set of genes are good predictors of a known phenotype ([1](#ref-1)). In either case, regression is a useful tool to quantify this association ([2](#ref-2)), either linear regression, for a continuous phenotype, or logistic regression, for a binary phenotype. Here we implement our own algorithms for linear and logistic regression optimized for work with sparse data (or other big scale daat) as input. We calculate beta coefficient estimates as well as their associated F statistics (linear regression), z-scores (logistic regression) and significance values. Our algorithm is based on iteratively reweighted least squares (IRWLS)([3](#ref-3)).

## References

<div id="refs" class="references">

<div id="ref-1">

M Asyali, D Colak, O Demirkaya, and M Inan. Gene expression profile classification: A review. Current Bioinformatics, 1(1):53–73, 2008.

</div>
  
<div id="ref-2">

W N van Wieringen, D Kun, R Hampel, and A L Boulesteix. Survival prediction using gene expression data: A review and comparison. Computational Statistics and Data Analysis, 53(5):1590–1603, 2009.

</div>

<div id="ref-3">

Paul W Holland and Roy E Welsch. Robust regression using iteratively reweighted least squares. Communications in Statistics theory and Methods, 6(9):813–827, 1977.
  
</div>
