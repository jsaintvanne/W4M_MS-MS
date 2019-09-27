W4M_MS-MS and Galaxy
====================

msPurity package
----------------

msPurity is a R package developped by Birmingham University searchers. It is associated with Galaxy wrappers to be able to run its functions on the Galaxy platform. We can find on Github each repository : for [R package]("https://github.com/computational-metabolomics/msPurity#mspurity-package-to-assess-precursor-ion-purity-process-fragmentation-spectra-and-perform-spectral-matching"), and for [Galaxy wrappers associated]("https://github.com/computational-metabolomics/mspurity-galaxy"). 
msPurity R package and associated Galaxy tools were developed to: 1) assess the spectral quality of fragmentation spectra by evaluating the "precursor ion purity". 2) process fragmentation spectra. And 3) perform spectral matching.

###Functionalities:

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


Associated paper  `msPurity: Automated Evaluation of Precursor Ion Purity for Mass Spectrometry Based Fragmentation in Metabolomics. Analytical Chemistry <http://pubs.acs.org/doi/abs/10.1021/acs.analchem.6b04358>`_ [1]

Use the following links for more details:

* Bioconductor: http://bioconductor.org/packages/msPurity/
* Vignette: https://bioconductor.org/packages/devel/bioc/vignettes/msPurity/inst/doc/msPurity-vignette.html
* Manual: http://bioconductor.org/packages/devel/bioc/manuals/msPurity/man/msPurity.pdf
* Galaxy implementation: https://github.com/computational-metabolomics/mspurity-galaxy
* Bioconda (stable): https://anaconda.org/bioconda/bioconductor-mspurity
* Conda (dev and testing): https://anaconda.org/tomnl/bioconductor-mspurity


###Install
--------------
####Bioconductor

.. code-block:: r

  if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
  BiocManager::install("msPurity")


####Github


.. code-block:: r

  library(devtools)
  install_github('computational-metabolomics/msPurity')


###Ref

[1] Lawson, T.N., Weber, R.J., Jones, M.R., Chetwynd, A.J., Rodriguez Blanco, G.A., Di Guida, R., Viant, M.R. and Dunn, W.B., 2017. msPurity: Automated Evaluation of Precursor Ion Purity for Mass Spectrometry Based Fragmentation in Metabolomics. Analytical Chemistry.


-----------
Our project
-----------

The [Workflow4Metabolomics]("http://workflow4metabolomics.org"), W4M in short, is a French infrastructure offering software tool processing, analysiong and annotating metabolomics data. It is based on the Galaxy platform.
Having no tool able to process MS-MS datas and wishing to integrate one, we are interested in the msPurity package. This package is very usefull, especially with its purity score which is very interesting for us. 

