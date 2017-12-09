------------------------------------------
SCRIPTS-I-WROTE-IN-MATLAB READ ME FILE
Alexandra Pulwicki Fall 2017
------------------------------------------

Hello! Here is a quick overview of the winter balance data anyalsis script I wrote. They have varying levels of commenting within the scripts.

------------------------------------------
FUNCTIONS

appendextra.m			Appends non-transect data (initial data work)
commentsearch.m			Pulls data based on comments
DataSubset.m			Subsets the input data 
EuclideanDistance.m		Calculates Euclidean distance
GridUpsizing.m			Increases grid cell size
KrigingR_G.m			Kriging of data for one glacier
KrigingR.m			Kriging of data for all glaciers
LinearRegression_Basic.m	Basic linear regression for all glaciers	
LinearRegression.m		Linear regression for all glaciers with cross validation and 
				model averaging
linspaceNDim.m			
MLRcalval.m			Executes linear regression with cross validation
ObsInCell.m			Average all values within one DEM cell
PlotTopoParameter_DD.m		3 glacier plotting - format 
PlotTopoParameter_IGS.m		3 glacier plotting - format for IGS paper
PlotTopoParameter.m		3 glacier plotting
pulldata.m			Extract data that corresponds to some category
pulldataSWE.m			Extract SWE data that corresponds to some category
RegressionKriging.m		Complete regression kirging on input data
RMSE.m				Calculate RMSE
SampledCell.m			Get values of matrix that correspond to sampling locations
saveFIG_HP.m			Save figure - Hydrological Processes format
saveFIG_IGS.m			Save figure - IGS format
saveFIG.m			Save figure	
SortNSelect.m			Order data and select a subset
SubsetSampleSize.m		
variofitAlex.m			Variogram fitting
variogramAlex.m			Variogram caluclation and fitting
WBnoise.m			Add noise to winter balance data based on std input


------------------------------------------
SCRIPTS

AccumGrad.m			Accumulation gradient in St. Elias continental side
BalanceDesign2.m		Second iteration of subsampling all measurement locations
BalanceDesign_BasicRegression.m	Final iteration of subsampling all measurement locations - 
				used for Paper II
BalanceDesign.m			Second iteration of subsampling all measurement locations
CentrelineDistance.m		Calculates distance from centreline
Import_Density.m		Initial data processing for density values
Import_SWE2.m			Second data processing for SWE data
Import_SWE.m			Initial data processing to obtain SWE values
Import_Topo.m			Initial data processing for topographic parameters
Import_Transect.m		Initial data processing for snow depth from transects
Import_Zigzag.m			Initial data processing for zigzag surveys
IWVT.m				Intergrated water vapour transport - unfinihsed work
Kriging.m			Kriging scrip with all the data set up and plotting
LargestGridSize.m		Changing sizes of DEM grid cells - unfinihsed work
MAIN.m				Calls all the 'Import' scripts
MeasurementLocations.m		Calculating the locations of snow depth measurements	
OPTIONS.m			Global variables and data processing options
Plot_Density.m			Plotting density
Plot_Maps.m			Plotting maps of the glacier (mostly initial stuff)
Plot_PaperII.m			Plotting for paper II (self contained)
Plot_PaperI.m			Plotting for paper I (self contained)
Plot_PresentationWB.m		Plotting for presentations (probably self contained)
Plot_SubsetAnimation.m		Plotting for initial stages of subset work
Plot_SWE.m			Plotting initial SWE data
Plot_TopoRegression.m		Plotting topographic parameters and their distributions
Plot_Zigzag.m			Plotting zigzag data
SmoothingTopo.m			Obtaining best smoothing for curvature
SnowDepth_Transect.m		Data processing and investigation for transect data
TestInterpMethodRMSE.m		RMSE and interpolation work
TopoRegression.m		Majority of work for topographic regressions (MLR and BMA)
Transferability.m		Transfer of regression coffficients to other glacier
VariogramCrossValidation.m	Cross validation with variogram data - unfinished
Variograms.m			Variogram work (initial)
WSMBdistribution.m		Monte Carlo analysis to investigate winter balance uncertainity
ZigzagRemoval.m			Exclude zigzag data 


------------------------------------------
.mat FILES WITH SAVED DATA

Basic.mat			Basic data about the glaciers
Full.mat			Full interpolation and SWE input data - used a lot for 
				plotting because this is the "best" estimate
Kriging.mat			Kriging winter balance
LR_SK_RK.mat			Interpolated winter balance
MLRrun.mat			MLR with different number of cross validation
Patterns.mat			Interpolation of subsets of WB input data found using basic 					linear regression (Paper II data)
PatternsTemp_HighNoise.mat	Subset work with WB inpertolated using LR with cross 
				validation and model averaging - high noise
PatternsTemp_LowNoise.mat	Subset work with WB inpertolated using LR with cross 
				validation and model averaging - low noise introduced
TopoBMS_MLR.mat			Comaprison of LR with Bayseian model averaging and multiple 
				linear regression
TopoSWE.mat			Topographic parameters and SWE input data
Variogram.mat			Variogram runs
varWSMB.mat			Monte Carlo runs for winter balance uncertainity
WSMBDistribution.mat		Monte Carlo runs for winter balance uncertainity

------------------------------------------
INPUT DATA

centrelinepoints2.csv
centrelinepoints.csv
d100_h0.xlsx
d200_h0.xlsx
d300_h0.xlsx
EmptyElevData.xlsx
FieldDataRevisedAP.xlsx
GlacierWP_UTM.xlsx
Narrow.jpg
sampling_d300_h0.xlsx
SPOTDEM_transectElev.csv
SPOT_SampledTopoParams_noZZ.csv
SPOT_TopoParams_noZZ.xlsx
SPOT_TopoParams.xlsx
summary_densitydata.xlsx
TransectMeasurementLocations.xlsx
VIWT_Oct-May_LowRes.nc
VIWT_Oct-May.nc
zigzag_corners_utm.xls

------------------------------------------
FOLDERS

BMS				Folder for doing Bayesian model averaging through R
Kriging				Folder for doing Dice Kriging through R
