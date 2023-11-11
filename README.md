# NMRphasing

NMRphasing is an R package designed for the correction of phase errors
in NMR data. While its primary focus is on 1D NMR data, it can also be
employed for 2D and 3D NMR data by processing one 1D NMR data file at a
time.

## I. NMRphasing installation

First of all, please ensure that you have the R package MassSpecWavelet,
which is a dependency for NMRphasing. You can install MassSpecWavelet
using the following example code:

if (!requireNamespace(“BiocManager”, quietly = TRUE))
install.packages(“BiocManager”) BiocManager::install(“MassSpecWavelet”)

NMRphasing is available on CRAN, and we can directly install it in R:

    install.packages("NMRphasing")

Alternatively, let’s proceed to install NMRphasing from GitHub.

    devtools :: install_github(repo = "ajiangsfu/NMRphasing",force = TRUE)  
    ## if you do not have old versions of NMRphasing, please remove force = TRUE

## II. Data format

The input data for NMRphasing can be in one of four formats:

1.  A vector of absorption spectrum.
2.  A complex vector.
3.  A data matrix with two columns of spectrum data, where the first
    column represents the absorption spectrum, and the second column
    represents the dispersion spectrum.
4.  A data frame with two columns of spectrum data, where the first
    column is for the absorption spectrum, and the second column is for
    the dispersion spectrum.

After installing NMRphasing, you can load a subset of example data from
our multiple metabolite spike-in experiment using the following code
snippet:

    library(NMRphasing)
    data("fdat")
    str(fdat)

    ## 'data.frame':    5891 obs. of  2 variables:
    ##  $ frequency_domain: cplx  -252643-294983i -221414-311592i -189411-330984i ...
    ##  $ ppm             : num  4 4 4 4 4 ...

In the code above, we load the dataset ‘fdat’ from the package.
‘fdat$frequency\_domain’ represents the complex vector used in the
following illustration, while ‘freq’ is utilized for plotting purposes.
If your NMR data uses ppm values, you can substitute ‘ppm’ accordingly.

    ## in order to make comparison, absorption part can be extracted
    fdat$Observed_Absorption = Re(fdat$frequency_domain)

    library(ggpubr)

    ## Loading required package: ggplot2

    p1 = ggplot(fdat, aes(x = ppm, y = Observed_Absorption)) +
          geom_line() + theme_bw() + labs(y = "Observed Absorption") +
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            axis.title.x = element_text(size = 6), 
            axis.title.y = element_text(size = 6))
    p1

![](README_files/figure-markdown_strict/unnamed-chunk-4-1.png)

The observed absorption spectrum without phase error correction looks
very bad.

## III. Example Code

In this package, there are currently nine NMR phase error correction
methods implemented. They are as follows: “NLS”,
“SPC\_DANM”,“MPC\_DANM”, “SPC\_EMP”, “MPC\_EMP”,“SPC\_AAM”,“MPC\_AAM”,
“SPC\_DSM” and “MPC\_DSM”.

### 1. NLS

This is our new shrinkage method, which is the fastest phase error
correction method since it does not involve any optimization step.

    fdat$Phased_Absoprtion = NMRphasing(specDatIn = fdat$frequency_domain, method = "NLS") 

    p2 = ggplot(fdat, aes(x = ppm, y = Phased_Absoprtion)) +
          geom_line() + theme_bw() + labs(y = "Phased Absoprtion") +
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            axis.title.x = element_text(size = 6), 
            axis.title.y = element_text(size = 6))
    ggarrange(plotlist = list(p1,p2),labels = c("Before","After"),font.label = list(size = 6),nrow = 2, ncol=1)

![](README_files/figure-markdown_strict/unnamed-chunk-5-1.png)

### 2. SPC\_DANM

The “SPC\_DANM” method employs the traditional single linear model
approach for phase error correction but incorporates a new optimization
function designed to minimize the difference between absolute area and
net area.

    fdat$Phased_Absoprtion = NMRphasing(specDatIn = fdat$frequency_domain, method = "SPC_DANM") 

    p2 = ggplot(fdat, aes(x = ppm, y = Phased_Absoprtion)) +
          geom_line() + theme_bw() + labs(y = "Phased Absoprtion") +
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            axis.title.x = element_text(size = 6), 
            axis.title.y = element_text(size = 6))
    ggarrange(plotlist = list(p1,p2),labels = c("Before","After"),font.label = list(size = 6),nrow = 2, ncol=1)

![](README_files/figure-markdown_strict/unnamed-chunk-6-1.png)

### 3. MPC\_DANM

This approach utilizes our new multiple linear model for phase error
correction along with our optimization function designed to minimize the
difference between absolute area and net area.

    fdat$Phased_Absoprtion = NMRphasing(specDatIn = fdat$frequency_domain, method = "MPC_DANM") 

    p2 = ggplot(fdat, aes(x = ppm, y = Phased_Absoprtion)) +
          geom_line() + theme_bw() + labs(y = "Phased Absoprtion") +
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            axis.title.x = element_text(size = 6), 
            axis.title.y = element_text(size = 6))
    ggarrange(plotlist = list(p1,p2),labels = c("Before","After"),font.label = list(size = 6),nrow = 2, ncol=1)

![](README_files/figure-markdown_strict/unnamed-chunk-7-1.png)

### 4. SPC\_EMP

This is an old phase error correction method based on the traditional
single model approach targeting on an existing optimization function:
entropy minimization with negative peak penalty.

    fdat$Phased_Absoprtion = NMRphasing(specDatIn = fdat$frequency_domain, method = "SPC_EMP") 

    p2 = ggplot(fdat, aes(x = ppm, y = Phased_Absoprtion)) +
          geom_line() + theme_bw() + labs(y = "Phased Absoprtion") +
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            axis.title.x = element_text(size = 6), 
            axis.title.y = element_text(size = 6))
    ggarrange(plotlist = list(p1,p2),labels = c("Before","After"),font.label = list(size = 6),nrow = 2, ncol=1)

![](README_files/figure-markdown_strict/unnamed-chunk-8-1.png)

### 5. MPC\_EMP

This is to apply our new multiple linear model approach for phase error
correction based on an existing optimization function: entropy
minimization with negative peak penalty.

    fdat$Phased_Absoprtion = NMRphasing(specDatIn = fdat$frequency_domain, method = "MPC_EMP") 

    p2 = ggplot(fdat, aes(x = ppm, y = Phased_Absoprtion)) +
          geom_line() + theme_bw() + labs(y = "Phased Absoprtion") +
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            axis.title.x = element_text(size = 6), 
            axis.title.y = element_text(size = 6))
    ggarrange(plotlist = list(p1,p2),labels = c("Before","After"),font.label = list(size = 6),nrow = 2, ncol=1)

![](README_files/figure-markdown_strict/unnamed-chunk-9-1.png)

### 6. SPC\_AAM

This is an old phase error correction method based on the traditional
single model approach with an existing optimization function: absolute
area minimization.

    fdat$Phased_Absoprtion = NMRphasing(specDatIn = fdat$frequency_domain, method = "SPC_AAM") 

    p2 = ggplot(fdat, aes(x = ppm, y = Phased_Absoprtion)) +
          geom_line() + theme_bw() + labs(y = "Phased Absoprtion") +
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            axis.title.x = element_text(size = 6), 
            axis.title.y = element_text(size = 6))
    ggarrange(plotlist = list(p1,p2),labels = c("Before","After"),font.label = list(size = 6),nrow = 2, ncol=1)

![](README_files/figure-markdown_strict/unnamed-chunk-10-1.png)

### 7. MPC\_AAM

This is to apply our new multiple linear model approach for phase error
correction based on an existing optimization function:absolute area
minimization.

    fdat$Phased_Absoprtion = NMRphasing(specDatIn = fdat$frequency_domain, method = "MPC_AAM") 

    p2 = ggplot(fdat, aes(x = ppm, y = Phased_Absoprtion)) +
          geom_line() + theme_bw() + labs(y = "Phased Absoprtion") +
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            axis.title.x = element_text(size = 6), 
            axis.title.y = element_text(size = 6))
    ggarrange(plotlist = list(p1,p2),labels = c("Before","After"),font.label = list(size = 6),nrow = 2, ncol=1)

![](README_files/figure-markdown_strict/unnamed-chunk-11-1.png)

### 8. SPC\_DSM

This is an old phase error correction method based on the traditional
single model approach with an existing optimization function: dispersion
summation minimization.

    fdat$Phased_Absoprtion = NMRphasing(specDatIn = fdat$frequency_domain, method = "SPC_DSM") 

    p2 = ggplot(fdat, aes(x = ppm, y = Phased_Absoprtion)) +
          geom_line() + theme_bw() + labs(y = "Phased Absoprtion") +
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            axis.title.x = element_text(size = 6), 
            axis.title.y = element_text(size = 6))
    ggarrange(plotlist = list(p1,p2),labels = c("Before","After"),font.label = list(size = 6),nrow = 2, ncol=1)

![](README_files/figure-markdown_strict/unnamed-chunk-12-1.png)

### 9. MPC\_DSM

This is to apply our new multiple linear model approach for phase error
correction based on an existing optimization function:dispersion
summation minimization..

    fdat$Phased_Absoprtion = NMRphasing(specDatIn = fdat$frequency_domain, method = "MPC_DSM") 

    p2 = ggplot(fdat, aes(x = ppm, y = Phased_Absoprtion)) +
          geom_line() + theme_bw() + labs(y = "Phased Absoprtion") +
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            axis.title.x = element_text(size = 6), 
            axis.title.y = element_text(size = 6))
    ggarrange(plotlist = list(p1,p2),labels = c("Before","After"),font.label = list(size = 6),nrow = 2, ncol=1)

![](README_files/figure-markdown_strict/unnamed-chunk-13-1.png)

## IV. NMRphasing R package general information

Version: 1.0.3 <Authors@R>: c( person(“Aixiang Jiang”, role = c(“aut”,
“cre”, “cph”), email = “<aijiang@bccrc.ca>”, comment = c(ORCID =
“0000-0002-6153-7595”)) ) Maintainer: Aixiang Jiang <aijiang@bccrc.ca>
Depends: R (&gt;= 4.3.0),stats Suggests: knitr, rmarkdown, ggpubr
VignetteBuilder: knitr Imports: baseline, splines, MassSpecWavelet
Description: There are three distinct approaches for phase error
correction, they are: a single linear model with a choice of
optimization functions, multiple linear models with optimization
function choices and a shrinkage-based method. The methodology is based
on our new algorithms and various references: Binczyk et al. (2015)
<doi:10.1186/1475-925X-14-S2-S5> Chen et al. (2002)
<doi:10.1016/S1090-7807(02)00069-1> de Brouwer (2009)
<doi:10.1016/j.jmr.2009.09.017> Džakula (2000)
<doi:10.1006/jmre.2000.2123> Ernst (1969)
<doi:10.1016/0022-2364(69)90003-1> Liland et al. (2010)
<doi:10.1366/000370210792434350> License: MIT + file LICENSE Encoding:
UTF-8 LazyData: true RoxygenNote: 7.2.3
