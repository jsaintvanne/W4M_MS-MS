# W4M_MS-MS and Galaxy

<details><summary>
msPurity package
</summary>

## msPurity package
msPurity is a R package developped by Birmingham University searchers. It is associated with Galaxy wrappers to be able to run its functions on the Galaxy platform. We can find on Github each repository : for [R package]("https://github.com/computational-metabolomics/msPurity#mspurity-package-to-assess-precursor-ion-purity-process-fragmentation-spectra-and-perform-spectral-matching"), and for [Galaxy wrappers associated]("https://github.com/computational-metabolomics/mspurity-galaxy"). 
msPurity R package and associated Galaxy tools were developed to : 

 1.  assess the spectral quality of fragmentation spectra by evaluating the "precursor ion purity". 
 2. process fragmentation spectra.
 3. perform spectral matching.

### Functionalities
---------------------------
* Assess the contribution of the targeted precursor of acquired fragmentation spectra by checking isolation windows using a metric called "precursor ion purity" (Works for both LC-MS(/MS) and DI-MS(/MS) data)
* Assess the anticipated “precursor ion purity” (see below) of XCMS LC-MS features and DIMS features where no fragmentation has been acquired
* Map fragmentation spectra to XCMS LC-MS features
* Filter and average MS/MS spectra from an LC-MS/MS dataset
* Create databases of LC-MS(/MS) spectra and associated annotations
* Perform spectral matching of query MS/MS spectra against library MS/MS spectra
* Export fragmentation spectra to MSP format
* Basic processing of DIMS data. Note that these functionalities are not actively developed anymore - see DIMSpy (https://github.com/computational-metabolomics/dimspy) for recommended alternative for DIMS data processing

**What is precursor ion purity?**

What we call "Precursor ion purity" is a measure of the contribution of a selected precursor peak in an isolation window used for fragmentation. The simple calculation involves dividing the intensity of the selected precursor peak by the total intensity of the isolation window. When assessing MS/MS spectra this calculation is done before and after the MS/MS scan of interest and the purity is interpolated at the recorded time of the MS/MS acquisition. Additionally, isotopic peaks can be removed, low abundance peaks are removed that are thought to have limited contribution to the resulting MS/MS spectra and the isolation efficiency of the mass spectrometer can be used to normalise the intensities used for the calculation.


Associated paper  [msPurity: Automated Evaluation of Precursor Ion Purity for Mass Spectrometry Based Fragmentation in Metabolomics. Analytical Chemistry]("http://pubs.acs.org/doi/abs/10.1021/acs.analchem.6b04358")

Use the following links for more details:

* Bioconductor : http://bioconductor.org/packages/msPurity/
* Vignette : https://bioconductor.org/packages/devel/bioc/vignettes/msPurity/inst/doc/msPurity-vignette.html
* Manual : http://bioconductor.org/packages/devel/bioc/manuals/msPurity/man/msPurity.pdf
* Galaxy implementation : https://github.com/computational-metabolomics/mspurity-galaxy
* Bioconda (stable) : https://anaconda.org/bioconda/bioconductor-mspurity
* Conda (dev and testing) : https://anaconda.org/tomnl/bioconductor-mspurity


### Install
---------------
#### Bioconductor

```r
  if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
  BiocManager::install("msPurity")
```

#### Github

```r
  library(devtools)
  install_github('computational-metabolomics/msPurity')
```

### Reference
--------------------
Lawson, T.N., Weber, R.J., Jones, M.R., Chetwynd, A.J., Rodriguez Blanco, G.A., Di Guida, R., Viant, M.R. and Dunn, W.B., 2017. msPurity: Automated Evaluation of Precursor Ion Purity for Mass Spectrometry Based Fragmentation in Metabolomics. Analytical Chemistry.[^1]

[^1]: [Lawson, T.N., Weber, R.J., Jones, M.R., Chetwynd, A.J., Rodriguez Blanco, G.A., Di Guida, R., Viant, M.R. and Dunn, W.B., 2017. msPurity: Automated Evaluation of Precursor Ion Purity for Mass Spectrometry Based Fragmentation in Metabolomics. Analytical Chemistry.]("http://pubs.acs.org/doi/abs/10.1021/acs.analchem.6b04358")
</details>
=======================

<details><summary>
W4M project
</summary>

## W4M project
The [Workflow4Metabolomics]("http://workflow4metabolomics.org"), W4M in short, is a French infrastructure offering software tool processing, analysing and annotating metabolomics data. It is based on the Galaxy platform.
Having no tool able to process MS-MS datas and wishing to integrate one, we are interested in the **msPurity package**. This package is very usefull, especially with its purity score which is very interesting for us. 
As described in this publication[^2], we have a basic protocol for injections describes in the following table :

![A basic injection of samples](https://github.com/jsaintvanne/W4M_MS-MS/blob/dev/images/basic_table.png)

During this protocol, we make pooled QC samples in which we mix **all samples** from **all classes**. That can causes the dilution of our searching compound and we can miss it during the fragmentation selection due to a very low intensity... As explain in the next picture with the A compound :

![A basic manipulation](https://github.com/jsaintvanne/W4M_MS-MS/blob/dev/images/basic_manip.png)
![Results from this basic manipulation](https://github.com/jsaintvanne/W4M_MS-MS/blob/dev/images/basic_results.png)

To avoid this kind of loss, we can propose an easy quite new protocol where we make **pooled QC samples only if they are from the same class**. With this kind of manipulation we should miss less compound because they will be in a good intensity in each pooled QC sample and the fragmentation selection should find them easily : 

![New W4M injection of samples](https://github.com/jsaintvanne/W4M_MS-MS/blob/dev/images/new_table.png)

During these injections, we respect the same protocol as previous one. We just change the end with injection for MS-MS runs. These injections should only concern **samples of the same class** that have been pooled. With this thing, we can obtain the following results :

![New W4M manipulation](https://github.com/jsaintvanne/W4M_MS-MS/blob/dev/images/new_manip.png)

As we can see in the previous picture, pooled QC have now **enough class specific compounds** to be detected then fragmented. We now have more MS-MS fragmented ions to match with all features peakpicked : 
![New W4M results](https://github.com/jsaintvanne/W4M_MS-MS/blob/dev/images/new_results.png)


[^2]: [Broadhurst, D., Goodacre, R., Reinke, S.N. et al. Guidelines and considerations for the use of system suitability and quality control samples in mass spectrometry assays applied in untargeted clinical metabolomic studies.  Metabolomics (2018) 14: 72.]("https://link.springer.com/article/10.1007%2Fs11306-018-1367-3")

</details>



`{note:title=Be Careful}I run it with a little modified msPurity package!!!{note}`
