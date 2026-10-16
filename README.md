# iOmicsPASS+
## Integrative -Omics Predictive Analysis of Subnetwork Signatures (Version II - *An R-package*)

### Installation

Clone or download this repository:

```bash
git clone https://github.com/Hiromikwl/iOmicsPASSplusRShiny.git
```

Install the package from R using `devtools` (run once):

```r
# Install devtools if you don't have it
install.packages("devtools")

# Install iOmicsPASS+ from the cloned directory
devtools::install_local("path/to/iOmicsPASSplusRShiny")
```

Replace `"path/to/iOmicsPASSplusRShiny"` with the actual path where you cloned the repository.

**Dependencies** installed automatically: `Rcpp`, `BH`, `corpcor`, `gplots`, `grid`, `gridExtra`, `gridGraphics`, `huge`, `nnet`, `pheatmap`, `pracma`, `RColorBrewer`, `shiny`, `bslib`, `DT`, `visNetwork`, `igraph`, `htmlwidgets`.

---

### Launching the Shiny GUI

After installation, launch the interactive GUI from the R console:

```r
library(iOmicsPASSplus)
run_iOmicsPASS()
```

This opens the iOmicsPASS+ web interface in your default browser. The GUI guides you through:

1. **Data upload** — upload up to three omics data files and a network file
2. **Network inference** *(optional)* — estimate a partial correlation network via `NetDeconvolute()`
3. **Subnetwork signature analysis** — run iOmicsPASS+ to identify subnetwork signatures discriminating phenotypic groups
4. **Interactive network visualisation** — explore, filter, and export the resulting network
5. **Prediction** — classify new samples using the identified signatures

To start the server without opening a browser (e.g. on a remote/headless server):

```r
run_iOmicsPASS(launch.browser = FALSE)
```

---

### Introduction

**iOmicsPASS+** is a R-package incorporating **iOmicsPASS** (<a href="https://www.nature.com/articles/s41540-019-0099-y" target="_blank">Koh et al., 2019</a>), extended to other types of -omics data allowing for flexibility and increasing usability. It includes several module including a network inference module `NetDeconvolute()` using graphical LASSO (glasso) to estimate a sparse inverse covariance matrix, creating a confounding-free partial correlation network among features from up to three types of -omics datasets.

**iOmicsPASS** has been improved to **iOmicsPASS+** allowing for higher flexibility and enabling applications to different types of omics data. Improvements include:

* **Specification of direction of association**\
  Users may now specify the direction for every pair of interacting or co-varying molecule by adding an additional column in the network file. However, only molecules that show a concordance in the sign of correlation in the empirical data as the user-specified direction of association will be considered.
  
* **Allows for a single network and input data**\
  Previously, at least two data and two networks were required as input. Now, users can input only one single data and create co-expressions among the variables in the data with a single network file.
  
* **Addition of a Network estimation module** `NetDeconvolute()`\
  Estimates a correlation network, linking the different features from up to three different data, using graphical LASSO (glasso) to estimate a sparse inverse covariance matrix, creating a confounding-free partial correlation network
  
* **New functions to help users compile and run iOmicsPASS using R**\
  Functions included in the R package facilitate users to build `INSTALL.iOmicsPASS()`, create input parameter file `createInputParam()`, create prior probabilities `createPrior()` and run the software `run.iOmicsPASS()` in the R-console.
  
* **Addition of a Prediction module** `Predict.iOmicsPASS()`\
  Uses the network signatures identified in the subentwork discovery module `run.iOmicsPASS()`to assign new samples to the phenotypic groups.
  
* **Adjustment for clinical information**\
  Users can incorporate clinical information such as age, gender and BMI, to modify the prior class probabilities used for assigning samples to the different groups.
  
  The figure below illustrates the overview of **iOmicsPASS+**
![Overview_iOmicsPASSplus](https://user-images.githubusercontent.com/37172948/142598004-93b4ed3a-42b9-428d-8620-92267601e840.png)
