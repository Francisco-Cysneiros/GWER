Authors: Yuri A. Araujo, Francisco José A. Cysneiros and Audrey H. M. A. Cysneiros

Year: 2020

Description:
This scripts set, developed in the R language, is used to evaluate 
the performance of the class of geographically weighted elliptical 
regression models (GWER). For its reproduction, it is necessary to
previous installation of some packages using the following commands.

install.packages('spgwr') 
install.packages('gwer')
install.packages('maptools')
install.packages('e1071')
install.packages('spdep') 
install.packages('xtable')
install.packages('robustbase') 


The set of scripts is described below.

- Application_Georgia file contains the routine used for the application of 
GWER models in dataset "GEORGIA" about the census of 1990 in Georgia, 
USA.

- Simulation_AICc file contains the routine used to the simulation study 
where the comparative analysis of the GWER model in relation to the 
GWR model in the context of data with outliers is carried out. 
AICc criterion is used to the selection of bandwidth.
- Simulation_CV file is similar to the simulation study above but CV criterion is considered for  bandwidth selection 

- Envelope contains the  functions to obtain the simulated envelopes 
for standardized residuals in the Application_Georgia file. The Envelope script 
is automatically load by the function source("envelope.R") into the script 
Application_Georgia.  

Note: 
For the function "source" to work normally the both files, Application_Georgia 
and Envelope, must be in the current working directory of the R process.
