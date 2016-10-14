%% ImSAnE: Tutorial with pointCloudDetector and Spherelike Fitter
% 
% Sample data: Drosophila embryo
%
% In this tutorial we detect the embryos apical surface using FIJI. The 
% point Cloud plugin (http://rsb.info.nih.gov/ij/plugins/point-cloud/) 
% is then used to export the point Cloud in a txt file, which we read using 
% the pointCloudDetector. The pointCloudDetector does not require fiji, any
% other segmentation tool that produces delimited text files will do. We
% simply chose FIJI for accessibility and popularity among biologists. 
% Next we fit a spherelike surface, generate a surface of interest
% and pullback the image data to various charts.
%
% This an example of how to detect and fit what we call spherelike surfaces,
% where the surface can be described as a slowly varying function of a
% preferred axis, which we call the z-axis. 
%
% Here is how to access the documentation of detectors and fitters that are
% introduced in this tutorial

% before running ImSAnE, we extract a pointCloud of the dataset in fiji: 

% load the data into fiji (http://fiji.sc/Fiji)
%   in the menu bar hit:
%       -> File -> Open
%   Navigate to the folder conaining the data, select the data and hit
%   open. Once fiji is done reading, you should see a new window containing
%   the data. We will now attempt to segment this data set. The goal is to
%   extract the edge of the drosophila embryo, that describes the surface
%   we are interseted in. 
% 
%   Step 1: Inspect the noise level. In this data, there may be some haze
%   around the embryo, which may throw of trhesholding algorithms. We can
%   get rid of this by subtracting the background. Here subtracting a fixed
%   number of 10 will be fine (, but we could also use the subtract
%   background function from fiji under the Process tab.)
%       -> Process -> Math -> Subtract
%   Enter 10 in the dialogue field and hit OK. You should see that the
%   haze is now reduced and overall, the image got a little bit dimmer. 

%   Step 2: Segment the image by creating a binary mask. This will convert
%   the stack into binary. In Fiji, there are plenty of algorithms
%   implemented to achieve a segmentation. In the menubar hit:
%   
%       -> Process -> Binary -> Make Binary 
%   To reproduce the data in this tutorial, in the Method section choose 
%   Huang, Background Dark and leave all fields un-ticked. Hit Ok.
%   You should now have the embryo (i.e. foreground) in black and the
%   surrounding background in white. As you scroll through the stack, you
%   may observe ocassional islands of background in the embryo, which we
%   have to get rid of next.
%
%   Step 3: Extract the edge of the binary mask. First we fill holes: 
%       -> Process -> Binary -> Fill Holes 
%   This should fill holes that the binarization algorithm did not consider
%   as foreground. Next extract the edge 
%       -> Process -> Find Edges 
%   This should produce an almost exclusively white image, which the
%   perimter of the embryo in black. Next we invert the lookup table 

%   Step 4: Invert the lookup table and export the point Cloud using the
%       point Cloud plugin. In the menubar hit:
%   
%       -> Image -> Lookup Tables -> Invert LUT
%   This should produce the inverse image, i.e. largely black with a thin
%   white line at the edge of the embryo. Now we are ready to export the
%   point cloud as a txt file. In the menubar hit: 
%
%       -> Plugins -> Point Cloud
%   In the dialog that opens up, tick 'output surface points only' and
%   'output as single file'. Hit Ok, and leave the header blank. Save the
%   data to your data Folder and note the name of the file. 

%   If all went well, we can now read an externally produced point Cloud 
%   into ImSAnE. 

%%
% 
doc surfaceDetection.pointCloudDetector


%% Initialize ImSAnE project
%
% We start by clearing the memory and closing all figures.
%
clear all; close all;
%%
%
% Setup a working directory for the project, where extracted surfaces,
% metadata and debugging output will be stored. Also specify the directory
% containing the data. 
%
[scriptPath,~,ext] = fileparts(matlab.desktop.editor.getActiveFilename);

dataDir    = fullfile(scriptPath, 'rawData');
projectDir = fullfile(scriptPath, 'projectFiles');


xp = project.Experiment(projectDir, dataDir);
cd(dataDir)

fileMeta                 = struct();
fileMeta.dataDir         = dataDir;
fileMeta.filenameFormat  = 'Time000005_8bit_bin2.tif'; % for full data sample use Time000000.tif
fileMeta.nChannels       = 1; % number of channels must be specified for matlab reader to function.
fileMeta.timePoints      = [5]; % for full data sample use 0;
fileMeta.stackResolution = [.52 .52 .52]; 
fileMeta.swapZT          = 0; % for full data sample use 1;

expMeta                  = struct();
expMeta.channelsUsed     = 1;
expMeta.channelColor     = 1;
expMeta.description      = 'Drosophila melanogaster embryo with GAP43 label';
expMeta.dynamicSurface   = 0;
expMeta.jitterCorrection = 0; % 1: Correct for sample translation
expMeta.fitTime          = fileMeta.timePoints(end); 
expMeta.detectorType     = 'surfaceDetection.pointCloudDetector';
expMeta.fitterType       = 'surfaceFitting.spherelikeFitter'; 
%

xp.setFileMeta(fileMeta);
xp.setExpMeta(expMeta);
%

xp.initNew();

%% Load data for surface detection and rescale to unit aspect ratio

xp.loadTime(5);

xp.rescaleStackToUnitAspect();

%% Detect the surface

myDetectOpts = struct('ssfactor', 1,'rmRadialOutliers', 1.5,...
         'fileName',[xp.fileMeta.dataDir,'/Time000005_8bit_bin2-pointCloud.txt'],...
         'zdim',3,'readerType','dlm','yInvert',1); 
xp.setDetectOptions(myDetectOpts);

%%

xp.detectSurface();


%% Inspect the point cloud in a cross section

inspectOptions= struct('dimension', 'x', 'value', 301, 'pointCloud', 'b');
xp.detector.inspectQuality(inspectOptions, xp.stack);

%% Inspect pointcloud in 3d
%
% Plot the points on the surface as a three dimensional point cloud.
% The subsampling factor reduces the number of points shown. 
%
ssfactor = 240;
xp.detector.pointCloud.inspect(ssfactor);

%% Fit the surface coarsly to prepare estimate of sample orientation

fitOptions = struct('R',6,'X',4,'Y',4,'e',2,'phase',0,'path',...
                    fullfile(projectDir, 'debugOutput'));
xp.setFitOptions(fitOptions);
xp.fitSurface();
%
% Determine sample orientation and perform fine fit
xp.determineROIfromFit();
fitOptions = struct('R',8,'X',3,'Y',3,'e',1,'phase',0,'shift',0);
xp.setFitOptions(fitOptions);
xp.fitSurface();
%
%% Inspect the fit in a cross section

zval = 400;
inspectOptions = struct('dimension','x','value',zval,'pointCloud','b');
xp.fitter.inspectQuality(inspectOptions, xp.detector, xp.stack);

%% Normally evolve the fitted surface by shift in pixels
%
% Depending on the part of the dataset we want to extract, we may want to 
% displace the surface by shift in the direction of the surface normal.
% Positive values move inwards, negative values outwards.
shift = 25;
xp.normallyEvolve(shift);

%% Set desired charts and generate generate pullbacks
%
% Now we have the surface in a region of the data we would like to extract,
% and next we want to pullback the data to the various available charts. 
% By setting the desired charts option of the fitter, we specify which charts to 
% compute.
%
% Implemented charts are
% 
% * 'cylinder1_index', fundamental cylinder chart
% * 'cylinder2_index', fundamental cylinder chart rotated by 180 degrees
% * 'cylinder1_proper', modification of cylinder chart, measuring distance on the surface      
% * 'cylinder2_proper', modification of cylinder chart, measuring distance on the surface rotated by 180 degrees     
% * 'polarLowerZ', chart covering lower pole,          
% * 'polarUpperZ', chart covering upper pole,             
% * 'anteriorEquidistant', chart covering anterior pole measuring distance on the surface
% * 'posteriorEquidistant', chart covering posterior pole measuring distance on the surface 
%
xp.fitter.setDesiredChart('cylinder1_proper',0);
xp.fitter.setDesiredChart('cylinder2_proper',0);
xp.fitter.setDesiredChart('polarLowerZ',0);
xp.fitter.setDesiredChart('polarUpperZ',0);
xp.fitter.setDesiredChart('anteriorEquidistant',0);
xp.fitter.setDesiredChart('posteriorEquidistant',0);
xp.generateSOI();

%% Pullback the stack to the desired charts
%
% Pass the region of interest and the current time to pull back the stack
% in the desired charts. This generates the data fields containing the
% pullback.
%
onionOpts = struct('nLayers', 3, 'layerDistance', 1, 'sigma', 20,'makeIP','both');
xp.SOI.pullbackStack(xp.stack, xp.currentROI, xp.currentTime, onionOpts);

%%  visualize the pullback on cylinder1_proper 
%
% Now we extract the data field from the surface of interest at the current
% time, which is the time of the fit.
%
data = xp.SOI.getField('data_MIP');
data = data(xp.tIdx(xp.currentTime));
%%
%
% then we specify, which chart over which patch we would like to use to 
% inspect the pullback. The following combinations of patch and transformation 
% are currently possible:
%
% * 'cylinder1_index', 'cylinder1'
% * 'cylinder1_index', 'cylinder1_proper'
% * 'cylinder2_index', 'cylinder2'
% * 'cylinder2_index', 'cylinder2_proper'
% * 'polarLowerZ_index', 'polarLowerZ'
% * 'polarUpperZ_index', 'polarUpperZ'
% * 'posteriorEquidistant_index', 'posteriorEquidistant'
% * 'anteriorEquidistant_index', 'anteriorEquidistant'
% 
% Remember that the chart needs to be specified as a desired chart, if you
% want to inspect the pullback to it. 
%
% This is an example that inspects the cylinder1_proper transform over the
% cylinder1_index patch. 
%
patchName     = 'cylinder1_index';
transformName = 'cylinder1';
pb = data.getPatch(patchName).getTransform(transformName).apply{1};
figure, imshow(pb',[],'InitialMagnification',66)

%% Save the surface of interest to disc
%
% Here we save the SOI using SOI.save. We set the following options:
%
% * dir:            The directory to save the SOI to.
% * imwriteOptions: Pullbacks are saved to image files using imwrite, we
% can pass options to change file format, compression etc. For example we
% could change this option to
% imwriteOptions = {'jp2', 'Mode', 'lossless'}; 
% * make8bit:       Often absolute intensities don't matter and 8 bit offers
% a large enough dynamic range. This options rescales the lookup table and
% converts to 8 bit before saving.
%
imwriteOptions = {'tif'};
saveDir = fullfile(projectDir, 'embryo');

options = struct('dir',saveDir,'imwriteOptions',{imwriteOptions},...
                    'make8bit',false);
xp.SOI.save(options) 



%% Where is the data? 
% 
% In the above specified saveDir, there is a fields folder, containg the 
% pullbacks in data, the metric and embedding for each desired chart as
% well as charts and transition maps in the atlas - each as a tif image for
% convenient inspection with other software such as Fiji. 

%% Load a SOI from disc
%
% All metadata is saved in SOI.xml. Pullbacks, geometry and any other data
% defined on the surface are saved to image files in subdirectories 
% reflecting the structure of patches and coordinates explained earlier in 
% this tutorial. We can reload a surface of interest with
% SOI.load(directory)
%
