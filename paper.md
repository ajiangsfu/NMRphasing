docker run --rm \
    --volume $PWD/paper:/data \
    --user $(id -u):$(id -g) \
    --env JOURNAL=joss \
    openjournals/inara

---
title: "NMRphasing: An R Package for 1D NMR Phase Error Correction"
authors:
  - name: Aixiang Jiang
    affiliation: [1, 2, 3]
  - name: Robert Nadon
    affiliation: [4]

affiliations:
  - index: 1
    name: "Department of Epidemiology, Biostatistics, and Occupational Health, McGill University, Montreal, Quebec, Canada"
  - index: 2
    name: "Department of Pathology and Laboratory Medicine, University of British Columbia, Vancouver, British Columbia, Canada"
  - index: 3
    name: "British Columbia Cancer Centre for Lymphoid Cancer, Vancouver, British Columbia, Canada"
  - index: 4
    name: "Department of Human Genetics, McGill University, Montreal, Quebec, Canada"
---


**Abstract**

Nuclear magnetic resonance (NMR) spectroscopy is widely employed to
identify and quantify chemical compounds. However, the presence of phase
errors in NMR data, where signal timing deviates from the true phase,
requires correction to prevent significant distortions. Current methods,
relying on a single simple linear regression model, face challenges in
correcting nonlinear phase errors, resulting in residual phase errors.
This limitation leads to residual phase distortions that can impact the
accuracy of chemical compound identification and quantification. To
address this challenge, we introduce the NMRphase R package, which
incorporates innovative phase error correction methods capable of
handling all types of phase errors while preserving existing correction
techniques. While designed primarily for 1D NMR data, NMRphase can also
be applied to 2D and 3D NMR data by processing one 1D NMR data file at a
time.
 
 
**Keywords:** NMR, phase error, phase error correction              

  -----------------------------------------------------------------------

**Introduction**

NMR, or Nuclear Magnetic Resonance, is a versatile scientific tool that
leverages signals emitted by nuclei when placed within magnetic fields.
Its significance spans various disciplines, including chemistry,
physics, biology, biochemistry, materials science, and medicine. Since
its inception in the 1920s, NMR has given rise to diverse variants, such
as Magnetic Resonance Imaging (MRI), MR Spectroscopic Imaging (MRSI),
Magnetoencephalography (MEG), pharmacological MRI (phMRI), and
functional MRI (fMRI).

These applications, while intricate, share a common foundation: the
detection and analysis of electromagnetic signals emitted by atomic
nuclei, including protons. When a radiofrequency (RF) field is applied,
protons within the sample become excited, absorbing RF radiation. This
absorbed energy briefly elevates the protons to a higher energy state
and increases their rotation speed. Subsequently, as the RF field
deactivates, the protons relax, returning to a lower energy state, and
emit electromagnetic radiation as a resonance signal. This resonance
signal is essentially a time series of decaying waveforms recorded in
NMR raw data, often referred to as NMR time domain data.

NMR signal analysis considers frequency, strength, and phase. Frequency
is a measure of how often a waveform repeats its pattern within a given
time frame. Strength indicates the maximum deviation of a waveform from
its baseline or zero level. Phase describes the timing or position of a
waveform relative to a reference point, indicating the shift in the
waveform's position along the time axis.

In order to separate signals from the same samples based on their
frequencies, NMR time domain data are transformed into NMR frequency
data usually through Fourier transformation(Fourier, 1822), where
signals should be displayed as sharpened narrow symmetric peaks. Despite
these advantages, when phase measurement errors exist, signal peaks are
distorted, causing alterations in signal shape, peak height, and area,
leading to inaccurate molecular identification and quantification.
Therefore, phase errors must be corrected before any data analysis.

Current existing phase error correction methods rely on a linear model,
almost always a simple linear model. The phase correction value
(opposite in sign to the phase error) serves as an unobservable
independent variable *y*, while a scaled frequency index acts as the
dependent variable *x* (Binczyk et al., 2015; Craig & Marshall, 1988;
Wang et al., 2009). An optimization process identifies the optimal
intercept and slope for the phase error correction model. This model is
then applied to the entire spectrum, making phase error correction at
each point (Binczyk et al., 2015; Chen et al., 2002; Ernst, 1969).

However, a single simple linear model applied to an entire sample\'s NMR
data is insufficient for handling non-linear phase errors. This
limitation leads to remaining phase errors and necessitates
labor-intensive manual correction. Additionally, the effectiveness of
existing methods depends on the optimization function. Various functions
prioritize specific signal aspects, influencing phase estimates and
downstream analysis. While entropy-based functions (Chen et al., 2002)
demonstrate superior performance compared to other methods (Brown et
al., 1989; Craig & Marshall, 1988), no single current function can
optimally address all phase errors, leaving room for improvement.

In this NMRphasing R package
(https://cran.r-project.org/web/packages/NMRphasing/index.html), we
propose novel approaches to correct phase errors in NMR spectroscopy,
addressing non-linear phase errors, and create a new optimization
function.

![Figure 1: Three phase error correction approaches in NMRphasing R package](./Figure1_phaseCorrectionApproaches.pdf)

Figure 1: Three phase error correction approaches in NMRphasing R package

Methods

NMRphasing offers nine phase error correction methods, comprising three
existing methods and six new ones. These methods fall into three
categories: phase error correction with a single simple model, phase
error correction with multiple models, and non-linear shrinkage. The
workflows for these three approaches are illustrated in Figure 1.

**Phase error correction with a single simple model**

When processing phase error correction with a single simple model for a
spectrum, the model can be represented as follows:

$$y = a + bx$$

Here, *y* signifies a phase correction value for an individual data
point, reflecting the phase error in the opposite direction. *a* serves
as an intercept for correcting the global phase error, also known as the
zero-order phase error, while *b* acts as a slope for addressing the
frequency-dependent linear phase error, commonly referred to as the
first-order phase error. The variable *x* represents the frequency
value, ppm value, or the index for a data point.

Since *y* is not observable, an optimization search process is required
to estimate *a* and *b*. To accomplish this, as illustrated in the left
panel of Figure 1, we initiate the process with two initial guesses for
the slope and intercept necessitating an optimization function. In the
current implementation of NMRphasing, five optimization functions are
available, as shown in the following:

• **Absolute area minimization (AAM)**

This method is based on the concept that a phase error free absorption
spectrum, corresponding to the real part of NMR frequency domain data,
should have the minimum integral of absolute intensity values among all
possible absorption spectra with different phase errors (de Brouwer,
2009; Džakula, 2000). The optimization objective function is given by:

$$\left( \widehat{a},\widehat{b} \right) = \underset{(a,b)}{argmin}{\sum_{k = 1}^{N}\left| A_{k}^{'} \right|}$$

Here, $\widehat{a}\ $and $\widehat{b}$ represent the values of *a* and
*b* that minimize the sum, and N is the number of terms in the
summation. *A′~k~* is the absorption value at the k-th data point with
the adjusted phase in the summation.

• **Entropy minimization with penalty (EMP)**

This method is employed to determine optimal parameters that not only
achieve minimal entropy but also penalize negative signals that should
not be present in a phase error free absorption spectrum. Entropy, in
this context, is defined as the negative summation of the absolute
intensity multiplied by the logarithm of the absolute intensity (Binczyk
et al., 2015). The penalty term, introduced to penalize negative values,
involves the square of negative values (de Brouwer, 2009). The
optimization objective function is given by:

$\left( \widehat{a},\widehat{b} \right) = \underset{(a,b)}{argmin}{( - \sum_{k = 1}^{N}{\left| A_{k}^{'} \right|*\ln\left| A_{k}^{'} \right|}}$
$+ \sum_{k = 1}^{N}{(A_{k}^{'}*I(A_{k}^{'} < 0))}^{2})$

Here, the notations are the same as for AAM. Additionally,
$I(A_{k}^{'} < 0)$ is an indicator function that equals 1 when
*A′~k\ ~*\< 0 and 0 otherwise. The formula combines two summations with
logarithmic and quadratic terms, respectively.

• **Dispersion summation minimization (DSM)**

This method is to minimize the integral of the dispersion spectrum
(Binczyk et al., 2015). The optimization objective function is given by:

$$\left( \widehat{a},\widehat{b} \right) = \underset{(a,b)}{argmin}{\sum_{k = 1}^{N}D_{k}^{'}}$$

*D′~k~* denotes the dispersion value at the k-th data point with the
adjusted phase in the summation.

• **Absolute dispersion summation minimization (ADSM)**

This method is founded on the premise that the absolute integral of the
dispersion spectrum, representing the imaginary part of NMR frequency
domain data, should approximate zero in the absence of phase errors
(Chen et al., 2002; Ernst, 1969; Jiang, 2024). The optimization
objective function is given by:

$$\left( \widehat{a},\widehat{b} \right) = \underset{(a,b)}{argmin}{|\sum_{k = 1}^{N}{D_{k}^{'}|}}$$

*D′~k~* denotes the dispersion value at the k-th data point with the
adjusted phase in the summation

• **Delta absolute net minimization (DANM)**

This represents our new method that takes into account both absolute and
net areas. While minimizing the absolute area, the objective is also to
maximize the net area for an absorption spectrum that effectively
minimizes negative values. This approach aims to reduce the difference
between the absolute area under a curve and the net area under a curve.
The optimization objective function is given by:

$$\left( \widehat{a},\widehat{b} \right) = \underset{(a,b)}{argmin}{(\sum_{k = 1}^{N}\left| A_{k}^{'} \right|} - \sum_{k = 1}^{N}A_{k}^{'})$$

Here, $A_{k}^{'}$ represents the absorption value at the k-th data
point.

After introducing the optimization functions provided by the NMRphasing
R package, let\'s delve back into the application of a single model for
phase error correction, as depicted in the left panel of Figure 1.

Once the first initial intercept and slope are determined, the phase
correction value for each data point, *y~k~*, can be calculated using
the simple linear model introduced earlier. Additionally, the absorption
and dispersion values at the k-th data point of the new phased NMR data,
${A^{'}}_{k}$ and ${D^{'}}_{k}$, can be obtained with the following
equations:

${A^{'}}_{k} =$ $A_{k}\cos{{(y}_{k\ })} - D_{k}\sin{{(y}_{k})}$

${D^{'}}_{k} =$ $A_{k_{l}}\sin{{(y}_{k\ })} + D_{k}\cos{{(y}_{k})}$

Here, *A~k~* represents the original or previous absorption value, and
*D~k~* represents the original or previous dispersion value at the k-th
data point.

Subsequently, the initial optimization function value for the entire
spectrum is calculated using a chosen optimization function listed
above. The process is then repeated with the selection of the second
initial intercept and slope. By comparing the two optimization values,
the one deemed either smaller (for a minimization function) or larger
(for a maximization function) is chosen as the current final
optimization value.

Following this, predefined convergence criteria, such as if the
difference between two versions of the function value is smaller than a
given cutoff or after running more than a specified number of
iterations, are used to assess the current optimization function value.
If the criteria are satisfied, the last version of *a* and *b* values
are considered the final phase correction parameters, which are then
employed to generate the final phase-corrected NMR data. If the criteria
are not met, a new iteration begins with a new guess for *a* and *b*,
repeating the process until the predefined criteria are satisfied.

The optimization process is conducted through the "optim" function
(Bélisle, 1992) in the "stats" R package.

**Phase error correction with multiple models**

The approach for phase error correction with multiple models involves
segmenting a spectrum into peak ranges and applying a simple regression
model within each range, as illustrated in the middle part of Figure 1.

As both observed absorption and dispersion spectra may be influenced by
phase errors, we utilize the corresponding magnitude spectrum as a guide
for segmentation. To achieve this, we first calculate the magnitude (M)
spectrum, defined as the square root of the sum of the squared
absorption and dispersion values for each data point:

$$M = \ \sqrt{A^{2} + \ D^{2}}$$

If the dispersion spectrum is unavailable, the Hilbert transformation
can be employed to derive the theoretical dispersion spectrum from the
absorption spectrum. However, it is important to note that this will be
less accurate than observed dispersion, making it preferable to use both
observed absorption and dispersion data.

Once the magnitude spectrum is obtained, major peak locations are
identified using continuous wavelet transform-based pattern matching
through the "peakDetectionCWT" function in the "MassSpecWavelet" R
package (Du et al., 2006). Subsequently, the dividing lines between peak
ranges, referred to as valleys, are defined as the locations of the
global minima between two adjacent identified major peaks. Each peak
range might contain one major peak with or without other peaks.

Following spectrum segmentation, phase error correction is applied
within each peak range using a single model, as described in the
previous section, \"Phase Error Correction with a Single Simple Model.\"
The phased peak ranges are then merged into a complete phase spectrum.

**Phase error correction with non-linear shrinkage**

Non-linear shrinkage also requires spectrum segmentation, as discussed
in the section \"Phase Error Correction with Multiple Models.\" However,
after spectrum segmentation, no regression model or optimization is
necessary. Instead, shrinkage is applied within each peak range as
illustrated in the right part of Figure 1.

While absorption and dispersion intensities might be distorted due to
phase errors, magnitude is not affected by phase, let alone phase
errors. A simple proof is as follows:

$$A = Mcos\theta$$

$$D = Msin\theta$$

$$M = \ \sqrt{A^{2} + \ D^{2}\ }$$

$$= \ \sqrt{{(Mcos\theta)}^{2} + \ {(Msin\theta)}^{2}\ }$$

$$= \ \sqrt{M^{2}\ (\cos^{2}\theta + \ \sin^{2}\theta)}$$

Since $\cos^{2}\theta + \ \sin^{2}\theta \equiv 1$, the magnitude
intensity *M* is unrelated to phase θ, and free of phase errors.
Therefore, the power of magnitude, $P = \ M^{2}$, is also unrelated to
phase θ, and free of phase errors.

Furthermore, the full width at half maximum (FWHM) is the same between a
power peak and its corresponding phase error free absorption peak
(Marshall & Verdun, 1990), and peak height is the same between a
magnitude peak and its corresponding phase error free absorption peak
since $M\  \geq \ A.$

Putting all of this together, our new approach to achieving a phased
absorption spectrum is to use the following formula to shrink power
intensity to be phase error free absorption intensity at the k-th data
point within a peak range:

$${A'}_{k} = P_{k}*\frac{max(\overrightarrow{M})}{max(\overrightarrow{P})}$$

The derived phase error free absorption peak ranges are then merged into
a complete phase spectrum.

**Overview of NMRphasing**

The NMRphasing R package includes 26 functions, divided into two
categories: 12 functions designed for user accessibility and 14
functions intended for developers. Of the 12 accessible functions, 11
are grouped into three distinct categories, each implementing a unique
phase error correction strategy. The 12th function acts as a wrapper,
enabling users to call any of the preceding nine functions by specifying
the appropriate method parameter.

The 12 functions are as follows:

**Phase correction with a single phase correction model (SPC)**

• SPC_AAM: A single model minimizing absolute area.

• SPC_EMP: A single model minimizing entropy with a negative peak
penalty.

• SPC_DSM: A single model minimizing dispersion summation/area.

• SPC_ADSM: A single model minimizing absolute dispersion
summation/area.

• SPC_DANM: A single model minimizing absolute and net area difference.

**Phase correction with multiple phase correction models (MPC)**

• MPC_AAM: Multiple linear phase correction models with absolute area
minimization.

• MPC_EMP: Multiple linear phase correction models that minimize entropy
with a negative peak penalty.

• MPC_DSM: Multiple linear phase correction models with dispersion
summation/area minimization.

• MPC_ADSM: Multiple linear phase correction models with absolute
dispersion summation/area minimization.

• MPC_DANM: Multiple linear phase correction models with delta net area
minimization.

**Phase error free spectrum developed with non-linear shrinkage (NLS)**

• NLS: Non-linear intensity shrinkage.

**Wrap-up function**

• NMRphasing: A wrap-up function capable of calling each of the above
nine functions.

**Examples of R code**

NMRphasing is available on CRAN
(https://cran.r-project.org/web/packages/NMRphasing/) and GitHub
(https://github.com/ajiangsfu/NMRphasing). Install it using one of the
following methods:

Install from CRAN:

**install**.packages(\"NMRphasing\")

or install from GitHub:

devtools :: install_github(repo = \"ajiangsfu/NMRphasing\",force = TRUE)

Note: If you don\'t have old versions of NMRphasing, remove force =
TRUE.

Then you can load the library and example data embedded in NMRphasing
with the following code:

**library**(NMRphasing)\
**data**(\"fdat\", package = \"NMRphasing\")

Use the following to investigate example data structure:

**str**(fdat)

\'data.frame\': 5891 obs. of 2 variables:\
\$ frequency_domain: cplx -252643-294983i -221414-311592i
-189411-330984i \...

\$ ppm : num 4 4 4 4 4 \...

This is a partial of a true NMR spectrum with two columns. The 1st
column is a complex vector of observed NMR data (real: absorption,
imaginary: dispersion). The 2nd column is ppm (parts per million),
indicating the relative frequency position within the magnetic field.

The following are examples of R code to run all nine phase error
correction methods in two different ways: either by directly calling a
method function or by obtaining the same results with the wrap-up
function.

> **NLS**

nlsres1 = NLS(specdat = fdat\$frequency_domain)

> or:

nlsres2 = NMRphasing(specDatIn = fdat\$frequency_domain, **method** =
\"**NLS**\")

The default setting for 'withBC' is TRUE, which tests for baseline bias
based on spline regression on the lowess line. If the maximum of
adjusted squared r is greater than 0.2, baseline correction is performed
with modified polynomial fitting.

If you set 'withBC' as FALSE, then no baseline bias will be tested and
corrected. The example code is:

Nlsres3 = NMRphasing(specDatIn = fdat\$frequency_domain, method =
\"NLS\", withBC = FALSE)

> **SPC_AAM**

saam1 = SPC_AAM(specdat = fdat\$frequency_domain)

> or:

saam2 = NMRphasing(specDatIn = fdat\$frequency_domain, **method** =
\"**SPC_AAM**\")

> **MPC_AAM**

maam1 = MPC_AAM(specdat = fdat\$frequency_domain)

> or:

maam2 = NMRphasing(specDatIn = fdat\$frequency_domain, **method** =
\"**MPC_AAM**\")

> **SPC_EMP**

semp1 = SPC_EMP(specdat = fdat\$frequency_domain)

> or:

semp2 = NMRphasing(specDatIn = fdat\$frequency_domain, **method** =
\"**SPC_EMP**\")

> **MPC_EMP**

memp1 = MPC_EMP(specdat = fdat\$frequency_domain)

> or:

memp2 = NMRphasing(specDatIn = fdat\$frequency_domain, **method** =
\"**MPC_EMP**\")

> **SPC_DANM**

sdanm1 = SPC_DANM(specdat = fdat\$frequency_domain)

or:

sdanm2 = NMRphasing(specDatIn = fdat\$frequency_domain, **method** =
\"**SPC_DANM**\")

**MPC_DANM**

mdanm1 = MPC_DANM(specdat = fdat\$frequency_domain)

or:

mdanm2 = NMRphasing(specDatIn = fdat\$frequency_domain, **method** =
\"**MPC_DANM**\")

**SPC_DSM**

sdsm1 = SPC_DSM(specdat = fdat\$frequency_domain)

or:

sdsm2 = NMRphasing(specDatIn =fdat\$frequency_domain, **method** =
\"**SPC_DSM**\")

**MPC_DSM**

mdsm1 = MPC_DSM(specdat = fdat\$frequency_domain)

or:

mdsm2 = NMRphasing(specDatIn = fdat\$frequency_domain, **method** =
\"**MPC_DSM**\")

**SPC_ADSM**

sadsm1 = SPC_ADSM(specdat = fdat\$frequency_domain)

or:

sadsm2 = NMRphasing(specDatIn =fdat\$frequency_domain, **method** =
\"**SPC_ADSM**\")

**MPC_ADSM**

madsm1 = MPC_ADSM(specdat = fdat\$frequency_domain)

or:

madsm2 = NMRphasing(specDatIn = fdat\$frequency_domain, **method** =
\"**MPC_ADSM**\")

Although the NMR input data in all the above examples are in complex
format, the "NMRphasing" wrap-up function can handle other NMR data
formats as well:

1\) A data matrix or a data frame with two columns of spectrum data,
where the first column is for the absorption spectrum, and the second
column is for the dispersion spectrum.

2\) A vector of the absorption spectrum; in this case, the parameter
\"absorptionOnly\" should be set as TRUE, and it\'s important to note
that the dispersion spectrum is transformed with the Hilbert transform,
which is not as accurate as true observed data. Therefore, it is
recommended to use both observed absorption and dispersion data for the
NMRphasing R package.

The phased absorption spectra from the above examples can be combined
and displayed with the original non-phased spectrum using the following
R code:

fdat = cbind(**fdat**,Re(**fdat**\$frequency_domain),nlsres1, saam1,
maam1, semp1, memp1,\
sdanm1, mdanm1, sdsm1, mdsm1, sadsm1, madsm1)\
library(**ggpubr**)\
plots = lapply(3:12, function(**i**){\
ggplot(**fdat**, aes(**x** = ppm, y = fdat\[,i\])) +\
geom_line() + theme_bw() +\
labs(**y** = \"Absorption\") +\
theme(**panel**.grid.major = element_blank(),\
panel.grid.minor = element_blank(),\
panel.background = element_blank(),

axis.line = element_line(\
**colour** = \"black\"),\
axis.title.x = element_text(**size** = 7),\
axis.title.y = element_text(**size** = 7),\
axis.text = element_text(**size** = 6))\
})\
ggarrange(**plotlist** = plots, nrow = 5, ncol=2,\
labels = c(\"A\", \"B\", \"C\", \"D\", \"E\", \"F\", \"G\", \"H\",
\"I\", \"J\", \"K\", \"L\"),\
vjust = 1)\
ggsave(\"PhaseCorrection.jpeg\", width = 8, height = 9, dpi = 600)

ggsave(\"PhaseCorrection.pdf\", width = 8, height = 9)

The output of the example code is displayed in Figure 2. In this small
NMR data example, the original data appears problematic (Figure 2A),
indicating the need for phase error correction. NLS performs very well
in phase error correction (Figure 2B). SPC_AAM does a poor job on phase
correction (Figure 2C), but MPC_AAM performs better, although some phase
errors remain in the small peak ranges (Figure 2D). SPC_EMP performs
well, although there is a slight distortion in the biggest peak close to
4 ppm (Figure 2E). MPC_EMP seems to perform adequately, but the middle
three significant peaks show a jump in the bottom part (Figure 2F). Both
SPC_DANM and MPC_DANM perform quite well (Figures 2G-H). SPC_DSM
performs poorly (Figure 2I), similar to SPC_AAM (Figure 2C). However,
even with multiple models, DSM, specifically MPC_DSM, cannot correct
phase errors for any peaks (Figure 2J), while MPC_AAM can at least
correct phase errors for the prominent peaks (Figure 2D). Nevertheless,
MPC_DSM is still slightly better than SPC_DSM. When we switch from DSM
to ADSM, both SPC_ADSM and MPC_ADSM perform significantly better than
SPC_DSM and MPC_DSM, respectively.

![Figure 2: Comparison of different phase error correction methods.](./Figure2_PhaseCorrection.pdf)

Figure 2: Comparison of different phase
error correction methods. A: original example data; B: NLS; C: SPC_AAM;
D: MPC_AAM; E: SPC_EMP; F: MPC_EMP; G: SPC_DANM; H: MPC_DANM; I:
SPC_DSM; J: MPC_DSM, K: SPC_ADSM; L: MPC_DSM.

**Conclusions**

In this application note, we introduce our new R package,
\"NMRphasing,\" specifically designed for phase error correction for 1D
NMR data. This package comprises nine distinct phase error correction
methods, categorized into three groups: non-linear shrinkage, a single
linear model, and multiple linear models.

Our innovative non-linear shrinkage approach, NLS, directly derives a
phase error free absorption spectrum from phase error free power and
magnitude spectra, employing different shrinkage coefficients for
distinct peak ranges.

The traditional single phase correction model, SPC, utilizes a
pre-defined optimization function to seek estimates of linear model
parameters (intercept and slope). This model is then employed to correct
phase errors across the entire spectrum. Within this package, we
implement five optimization functions (AAM, ADSM, DANM, DSM, EMP).

Our approach with multiple phase correction models, MPC, segments a
spectrum into peak ranges and applies a linear model to correct phase
errors within each range. The same five optimization functions are
implemented for multiple models as for a single model.

We provide example code for each of these nine methods in two ways:
either by directly calling the method function or invoking it from a
wrap-up function, offering users alternative choices.

Using a partial real-world NMR spectrum, we present the original
spectrum alongside the phased spectra from the nine different methods.
The results indicate that our three new methods---NLS, SPC_DANM, and
MPC_DANM---perform the best, followed by SPC_EMP, SPC_ADSM, MPC_EMP, and
MPC_ADSM, with MPC_AAM slightly behind. MPC_DSM, SPC_AAM, and SPC_DSM
rank as the least effective.

From these observations, it is evident that NLS is the optimal choice.
Overall, MPC outperforms SPC with the same optimization function.
However, if an optimization function is sufficiently robust, SPC with
this function can match the performance of MPC with the same function
(e.g., DANM).

It\'s important to note that these observations are based on a single
example dataset, and their general applicability to other NMR data may
vary. Users are encouraged to apply all methods to their datasets for
testing and decide which one suits their needs. If testing all methods
is challenging, we recommend using the NLS method for phase error
correction.

While designed primarily for 1D NMR data, NMRphasing can also be applied
to 2D and 3D NMR data by processing one 1D NMR data file at a time.

**References**

Bélisle, C. J. P. (1992). Convergence theorems for a class of simulated
annealing algorithms on ℝd. *Journal of Applied Probability*, *29*(4),
885--895. Cambridge Core. https://doi.org/10.2307/3214721

Binczyk, F., Tarnawski, R., & Polanska, J. (2015). Strategies for
optimizing the phase correction algorithms in Nuclear Magnetic Resonance
spectroscopy. *Biomedical Engineering Online*, *14 Suppl 2*(Suppl 2),
S5. https://doi.org/10.1186/1475-925X-14-S2-S5

Brown, D. E., Campbell, T. W., & Moore, R. N. (1989). Automated phase
correction of FT NMR spectra by baseline optimization. *Journal of
Magnetic Resonance (1969)*, *85*(1), 15--23.
https://doi.org/10.1016/0022-2364(89)90315-6

Chen, L., Weng, Z., Goh, L., & Garland, M. (2002). An efficient
algorithm for automatic phase correction of NMR spectra based on entropy
minimization. *Journal of Magnetic Resonance*, *158*(1), 164--168.
https://doi.org/10.1016/S1090-7807(02)00069-1

Craig, E. C., & Marshall, A. G. (1988). Automated Phase Correction of Ft
Nmr-Spectra by Means of Phase Measurement Based on Dispersion Versus
Absorption Relation (Dispa). *Journal of Magnetic Resonance*, *76*(3),
458--475. https://doi.org/10.1016/0022-2364(88)90350-2

de Brouwer, H. (2009). Evaluation of algorithms for automated phase
correction of NMR spectra. *Journal of Magnetic Resonance*, *201*(2),
230--238. https://doi.org/10.1016/j.jmr.2009.09.017

Du, P., Kibbe, W. A., & Lin, S. M. (2006). Improved peak detection in
mass spectrum by incorporating continuous wavelet transform-based
pattern matching. *Bioinformatics (Oxford, England)*, *22*(17),
2059--2065. https://doi.org/10.1093/bioinformatics/btl355

Džakula, Ž. (2000). Phase Angle Measurement from Peak Areas (PAMPAS).
*Journal of Magnetic Resonance*, *146*(1), 20--32.
https://doi.org/10.1006/jmre.2000.2123

Ernst, R. R. (1969). Numerical Hilbert Transform and Automatic Phase
Correction in Magnetic Resonance Spectroscopyt. *Journal of Magnetic
Resonance*, *1*, 7--26.

Fourier, J.-B.-J. (1822). *Théorie analytique de la chaleur*. F. Didot.

Jiang, A. (2024). Phase Error Correction in Magnetic Resonance: A Review
of Models, Optimization Functions, and Optimizers in Traditional
Statistics and Neural Networks. *Preprints*.
https://doi.org/10.20944/preprints202409.2252.v1

Marshall, A. G., & Verdun, F. R. (1990). *Fourier transforms in NMR,
optical, and mass spectrometry: A user's handbook*. Elsevier.

Wang, T., Shao, K., Chu, Q., Ren, Y., Mu, Y., Qu, L., He, J., Jin, C., &
Xia, B. (2009). Automics: An integrated platform for NMR-based
metabonomics spectral processing and data analysis. *BMC
Bioinformatics*, *10*, 83. https://doi.org/10.1186/1471-2105-10-83
