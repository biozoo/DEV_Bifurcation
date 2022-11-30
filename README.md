Datasets and documentation of all analytical procedures from a scientific paper, Anticipating the occurrence and type of critical transitions, that will be formally published in Science Advances.

Authors: Florian Grziwotz+, Chun-Wei Chang+, Vasilis Dakos, Egbert van Nes, Markus Schwarzländer, Oliver Kamp, Martin Heßlers, Isao Tokuda, Arndt Telschow*, Chih-hao Hsieh*

+Equal contribution; *Correspondence to: Arndt Telschow (telschow@uni-muenster.de) and Chih-hao Hsieh (chsieh@ntu.edu.tw)

###############################################
# Terms of Use for empirical dataset
 
 (I)   You must acknowledge the use of content.
 (II)  Monitoring data is made available for use in activities of a non-profit nature only.
 (III) Users must contact the database administrator listed in Table S2 before using the dataset for any publications, 
      including conference presentations as well as handouts and presentation materials for meetings such as committees and councils. 
      Co-authorship may be required for some publications depending on the way the time series data is to be used.

###############################################
# Note for analyzing empirical dataset:

   This R file demonstrates the DEV analysis applied in the following five datasets (Fig. 2A-O)
   (1) Microcosm experiment data provided by Dr. Egbert H. van Nes (Fig. 2A-C; microcosm.txt). 
       The complete dataset analyzed in the original paper (Veraart et al. 2012 Nature) has been published on online DANS repository, doi:10.17026/dans-ztg-93aw 
   (2) Voice experiment data are provided by Dr. Isao Tokuda (Fig. 2D-F; voice.txt) 
   (3) Cellular ATP experiment data are offered by Dr. Markus Schwarzländer (Fig. 2G-I; mitochondria.txt)
   (4) Greenhouse earth climate data is publicly available and were downloaded from the World Data Center for Paleoclimatology, National Geophysical Data Center, Boulder, Colorado (http:// www.ncdc.noaa.gov/paleo/data.html). This dataset has also been used for computing the other EWS in a previous study (Dakos et al. PNAS 2008)(Fig. 2J-L; climate.txt)
   (5) Frequency data provided by Bonneville Power Administration (Fig. 2M-O; blackout.txt).
   

###############################################
# Note for analyzing synthetic dataset:
To significantly reduce the computation time running the analysis including 100 replicates, the analyses presented in the R code for producing Fig. 1, 3, and S14 are based on one time series replicate per mathematical model. Due to process noise (as well as observation error), the DEV results produced from analyzing one replicate look more variable and change less smoothly. However, the findings based on analyzing one replicate still reveal the same dynamical behaviours qualitatively (e.g., increasing trend of |DEV|) and quantitatively (e.g., |DEV|->1 at the bifurcation point).

To obtain the results presented in the main text, the R-file can be easily modified in four steps:
1. Calculate all results of this R-file 100 times.
2. Save each set of the results.
3. Average across all results associated with the same time step.
4. Plot the averaged results instead of the results based on one replicate as shown in this file.

The EDM analysis in all R files are based on rEDM package (ver 1.2.3) that can be downloaded via the following link: https://github.com/cran/rEDM/releases/tag/1.2.3
