# Bag of DAGs: Flexible Nonstationary Modeling of Spatiotemporal Dependence 

This is an R package under development for the [manuscript](https://doi.org/10.48550/arXiv.2112.11870) availale in arXiv. 
You can install this package using the following code 
```
devtools::install_github("jinbora0720/bags")
```

Once it is installed, you can reproduce all analyses in the manuscript with the code available in a separate GitHub repository called 
[GBAGs](https://github.com/jinbora0720/GBAGs). 
Refer to `README` in GBAGs repository for further instructions to reproduce the analyses by sections of the manuscript. 


For detailed guidelines of how to use functions of the package `bags`, refer to either help page of each function or 
[this vignette](https://jinbora0720.github.io/media/BAGs/example.html). 
The vignette demonstrates a full analysis of a univariate Bayesian regression model with Gaussian errors whose spatiotemporal random effects 
are modeled with a G-BAG prior. It shows how to 1) generate random data from prior sampling of G-BAGs, 
2) compute a nonstationary covariance matrix derived from G-BAG, 3) run posterior sampling from the model, and 
4) make inferences using posterior samples. 

**Note:** The current version is tested on R version 4.1.0 under Ubuntu 20.04.2.

**FAQ:** In case installation fails, confirm if all mandatory tools are available to compile the package. For more help, check [R for macOS](https://mac.r-project.org/tools/) or [R for Windows](https://cran.r-project.org/bin/windows/Rtools/).
