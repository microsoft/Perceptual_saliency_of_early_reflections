# Polar Threshold for Source Separation
The scripts in this folder were used to calculate the angular threshold for source separation of concurrent sources in polar direction. They are not needed for the main model to run, the results are stored in `Code/Helper/polarThreshold_full_all_VPs_absoluteLatAngle_lat90added_a396ae7cf8afe4c00387625fa804ebe0.mat`


Required Resources:

* Parallel Computing Toolbox
	* A few parfor loops are used in the scripts. You can get around this by changing all parfor loops to for loops.
* Image Processing Toolbox

* [The HUTUBS HRTF Database](https://depositonce.tu-berlin.de/handle/11303/9429) 
* [The SOFiA sound field analysis toolbox](https://audiogroup.web.th-koeln.de/SOFiA_wiki/WELCOME.html)
* [AKtools](https://www.ak.tu-berlin.de/menue/publications/open_research_tools/aktools/)
* [Auditory Modeling Toolbox](https://amtoolbox.org) (tested with version 0.10.0)

Other Third Party functions are included in the `ThirdParty` folder

`a_startup.m`
Run this script first to temporarily add the required third party toolboxes to your search path

`b_calculateBaumgartner.m`
Run this script to model two concurrent sources and evaluate the perceived location of the summed signal.

`c_optimizeFitRange.m`
Run this script to find the optimum fitting interval for `d_calculatePolarThreshold.m` (minimized total fit error), according to section A.3.

`d_calculatePolarThreshold.m`
Run this script to calculate the polar threshold for source separation by thresholding the Wasserstein Distance between localization patterns of two concurrent sources compared to one source.

`e_averagePolarThresholdLateralAnglesAndAddLat90.m`
Run this script to reformat the results of `d_calculatePolarThreshold.m` to absolute lateral angles as used in `detectReflections.m`