# iOmicsPASS+
## Integrative -Omics Predictive Analysis of Subnetwork Signatures (R Shiny GUI)

### System Requirements

| Requirement | Details |
|---|---|
| R version | 4.1 or later (latest recommended) |
| RStudio | Latest Desktop version (recommended IDE) |
| C++ compiler | Rtools 4.x (Windows) · Xcode CLT (macOS) |
| RAM | 4 GB minimum · 8 GB recommended for large datasets |

---

### Installation

#### Windows

**Step 1 — Install R**  
Download and run the Windows installer from CRAN: https://cran.r-project.org/bin/windows/base/

**Step 2 — Install RStudio**  
Download the latest RStudio Desktop: https://posit.co/download/rstudio-desktop/

**Step 3 — Install Rtools (C++ compiler)**  
iOmicsPASS+ contains compiled C++ code. Download Rtools from https://cran.r-project.org/bin/windows/Rtools/ and match the version to your R version (e.g. Rtools4.5 for R ≥ 4.5.0). Run the installer with default options.

**Step 4 — Verify the compiler is detected**  
In the RStudio Console:
```r
if (!requireNamespace("pkgbuild", quietly = TRUE)) install.packages("pkgbuild")
pkgbuild::has_build_tools(debug = TRUE)
```
The output should print `TRUE`. If you see `FALSE`, restart RStudio and try again.

**Step 5 — Install iOmicsPASS+ from GitHub**  
In the RStudio Console:
```r
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
remotes::install_github("Hiromikwl/iOmicsPASSplusRShiny")
```
All R package dependencies install automatically. The first install takes a few minutes.

---

#### macOS

**Step 1 — Install R**  
Download the macOS `.pkg` installer from https://cran.r-project.org/bin/macosx/  
- Apple Silicon (M1/M2/M3/M4): choose the **arm64** build  
- Intel Macs: choose the **x86_64** build

**Step 2 — Install RStudio**  
Download the latest RStudio Desktop: https://posit.co/download/rstudio-desktop/

**Step 3 — Install Xcode Command Line Tools (C++ compiler)**  
Open Terminal and run:
```bash
xcode-select --install
```
A dialog box will appear — click Install. This takes approximately 5–10 minutes and requires an internet connection.

**Step 4 — Verify the compiler is detected**  
In the RStudio Console:
```r
if (!requireNamespace("pkgbuild", quietly = TRUE)) install.packages("pkgbuild")
pkgbuild::has_build_tools(debug = TRUE)
```
Should print `TRUE`. If `FALSE`, restart RStudio after the Xcode CLT installation completes.

**Step 5 — Install iOmicsPASS+ from GitHub**  
In the RStudio Console:
```r
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
remotes::install_github("Hiromikwl/iOmicsPASSplusRShiny")
```
All R dependencies install automatically. First install takes a few minutes.

---

### Launching the Application

Each time you want to use iOmicsPASS+ R Shiny, open RStudio and run:

```r
library(iOmicsPASSplus)
run_iOmicsPASS()
```

> **Note:** The app opens in your default web browser. Keep the RStudio Console window open — closing it stops the app server. If the browser does not open automatically, copy the URL printed in the Console (e.g. `http://127.0.0.1:XXXX`) and paste it into your browser.

---

### Step-by-Step Workflow

**1. Try the example dataset**  
The quickest way to explore all features is to use the bundled Tulip dataset (obese insulin-resistant vs lean insulin-sensitive):
- On the **Home** tab, click **Load Example Data**
- This loads proteomics + miRNA data, PPI and miRNA-target networks, and Human pathway annotations
- Navigate to **Network Setup**, then **iOmicsPASS Analysis** to run the full workflow

**2. Upload your own data**  
Go to the **Data Upload** tab. Files must be plain text (`.txt` or `.csv`), tab- or comma-separated:
- **Data X**: required — rows = features, columns = samples
- **Data Y / Z**: optional second/third omics type
- **Phenotype file**: sample IDs + group/class labels
- **Pathway annotation**: upload your own, or click *Use example (Human pathways)*

> Sample IDs in the phenotype file must exactly match the column headers in your data files.

**3. Configure network files**  
In the **Network Files** tab, upload edge list files with columns: `nodeA`, `nodeB`, `sign` (1 or −1):
- *Within-network*: edges within Data X (e.g. protein–protein interactions)
- *Between-network*: edges between Data X and Data Y (e.g. miRNA–gene targets)

If no network is available, enable the glasso network inference option in the analysis settings.

**4. Run the analysis**  
- Go to **iOmicsPASS Analysis** and configure parameters (lambda range, CV folds, etc.)
- Click **Run iOmicsPASS+**
- Results appear in the **Results** tab once the analysis completes

**5. Explore results**  
- **Network**: interactive visualisation — colour edges by group direction, toggle curved/straight edges
- **Network Enrichment**: bubble plot of pathway enrichment for the subnetwork
- **Pathway Enrichment**: gene-level pathway enrichment with FDR cutoff
- **Download All Results (ZIP)**: exports all tables, plots (PNG/JPEG), and enrichment results

---

### Input File Format Reference

| File | Format requirements |
|---|---|
| Data X / Y / Z | Tab `.txt` or comma-separated `.csv`. Features in rows, samples in columns. Data values start from the second column. |
| Phenotype | Tab `.txt` or comma-separated `.csv`. First column: sample IDs; second column: group labels. Additional columns for covariates (optional). |
| Network (within/between) | 3 columns: Node A, Node B, sign (1 = positive, −1 = negative). One edge per row. |
| Pathway annotation | 3 columns: Feature name, Pathway ID, Pathway annotation. |

---

### Troubleshooting

**Compiler not found**
- Run `pkgbuild::has_build_tools(debug = TRUE)` in RStudio
- Windows: reinstall Rtools and restart RStudio
- macOS: run `xcode-select --install` in Terminal

**Package install fails**
- Ensure you have R 4.1 or later
- Try updating existing packages first: `update.packages(ask = FALSE)`
- Check internet connectivity — GitHub must be reachable

**App does not open in browser**
- Check the RStudio Console for the URL (e.g. `http://127.0.0.1:XXXX`)
- Copy and paste that URL directly into your browser

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
