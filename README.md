#The supplementary materials R script for the scientific research entitled, " Anticipating the occurrence and type of critical transitions"
Authors: Florian Grziwotz+, Chun-Wei Chang+, Vasilis Dakos, Egbert van Nes, Markus Schwarzländer, Martin Hess, Oliver Kamps, Isao Tokuda, Arndt Telschow*, Chih-hao Hsieh*
+Equal contribution; *Correspondence to: Arndt Telschow (telschow@uni-muenster.de) and Chih-hao Hsieh (chsieh@ntu.edu.tw)

###############################################
# Note for analysing synthetic dataset:
To significantly reduce the computation time running the analysis including 100 replicates, the analyses presented in the R code for producing Fig. 1, 3, and S14 are based on one time series replicate per mathematical model. Due to process noise (as well as observation error), the DEV results produced from analyzing one replicate look more variable and change less smoothly. However, the findings based on analyzing one replicate still reveal the same dynamical behaviours qualitatively (e.g., increasing trend of |DEV|) and quantitatively (e.g., |DEV|->1 at the bifurcation point).

# To obtain the results presented in the main text, the R-file can be easily modified in four steps:
1. Calculate all results of this R-file 100 times.
2. Save each set of the results.
3. Average across all results associated with the same time step.
4. Plot the averaged results instead of the results based on one replicate as shown in this file.

###############################################
# Note for analysing synthetic dataset:
This R file demonstrates the DEV analysis applied in the following three datasets (Fig. 2D-L)
 (1) Voice experiment data are provided by Dr. Isao Tokuda (Fig. 2D-F) 
 (2) Cellular ATP experiment data are offered by Dr. Markus Schwarzlander (Fig. 2G-I)
 (3) Greenhouse earth climate data is publicly available and were downloaded from the World Data Center for Paleoclimatology, National Geophysical Data Center, Boulder, Colorado (http:// www.ncdc.noaa.gov/paleo/data.html). This dataset has also been used for computing the other EWS in a previous study (Dakos et al. PNAS 2008)(Fig. 2J-L)

###############################################
# Terms of Use for the dataset (1) & (2)
  I.   You must acknowledge the use of content.
  II.  Monitoring data is made available for use in activities of a non-profit nature only.
  III. Users must contact the database administrator: 
      (1) Dr. Isao Tokuda (voice data; isao@fc.ritsumei.ac.jp)\t
      (2) Dr. Markus Schwarzlander (cellular ATP data; markus.schwarzlander@uni-muenster.de)\t
      before using the dataset for any publications, including conference presentations as well as handouts and presentation materials for meetings such as committees and councils. Co-authorship may be required for some publications depending on the way the time series data is to be used.
###############################################
