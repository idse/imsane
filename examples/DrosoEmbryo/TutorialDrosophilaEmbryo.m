%% ImSAnE: Drosophila melanogaster embryo tutorial
%
% In this tutorial we detect the embryos apical surface during
% cellularization, fit a spherelike surface, generate a surface of interest
% and pullback the image data to various charts.
%
% This an example of how to detect and fit what we call spherelike surfaces,
% where the surface can be described as a slowly varying function of a
% preferred axis, which we call the z-axis. It also shows how to choose 
% from two different detectors and describes the options of loading an 
% externally provided point cloud. 
%
% After running the script, measurements in ImSAnE are explained and
% compared to uncorrected measurements. 
%
% Here is how to access the documentation of detectors and fitters that are
% introduced in this tutorial

%%
% 
% doc surfaceDetection.fastCylinderDetector
% doc surfaceDetection.radialEdgeDetector
%
% doc surfaceFitting.spherelikeFitter


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

dataDir = fullfile(scriptPath, 'rawData');
projectDir = fullfile(scriptPath, 'projectFiles');

%%
%
% Start by creating an experiment object, optionally pass on the project
% directory (otherwise it will ask), and change into the directory of the
% data. This serves as a frontend for data loading, detection, fitting etc.
%
xp = project.Experiment(projectDir, dataDir);
%%
%
% optionally switch working directory into the data directory.
%
cd(dataDir)

%% Set file and experiment metadata
%
% Set required additional information on the files
% 
% We assume one individual image stack for each time point, labeled by time.
% To be able to load the stack, we need to tell the project where the data
% is, what convention is assumed for the file names, available time points
% and the stack resolution. Options to modules in ImSAnE are organised in
% matlab structures, that is a pair of field name and value are provided
% for each option. 
%
% The following file metadata information is required   
%
% * 'directory'       , the project directory (full path)
% * 'dataDir'         , the data directory (full path)
% * 'filenameFormat'  , fprintf type format spec of file name
% * 'timePoints'      , list of times available stored as a vector
% * 'stackResolution' , stackresuolution in microns, e.g. [.25 .25 1]
%
%


%% 
%
% The following file metadata information is optional
%
% * 'imageSpace'      , bit depth of image, such as uint16 etc, defined in Stack class
% * 'stackSize'       , size of stack in pixels per dimension [xSize ySize zSize]
% * 'swapZT'          , for some datasets time is the third dimension, and z 
% the fourth. Swap = 1 if this is the case. 
%
% This tutorial uses a SPIM dataset, that is 2 fold downsampled, for
% speedup and memory requirements. There is also a full size data file
% included, however, it is strongly recommended to have at least 6 GB of
% free memory when trying this dataset. Moreover, some fileMeta and
% detector options should be changed, where indicated. 
fileMeta                 = struct();
fileMeta.dataDir         = dataDir;
fileMeta.filenameFormat  = 'Time%06d_8bit_bin2.tif'; % for full data sample use Time000000.tif
fileMeta.timePoints      = [5]; % for full data sample use 0;
fileMeta.stackResolution = [.5 .5 .5]; 
fileMeta.swapZT          = 0; % for full data sample use 1;
fileMeta.nChannels       = 1;

%%
%
% Set required additional information on the experiment. A verbal data set
% description, Jitter correct by translating the sample, which time point to 
% use for fitting, etc.
%
% The following project metadata information is required 
%
% * 'channelsUsed'   , the channels used, e.g. [1 3] for RGB
% * 'channelColor'   , mapping from element in channels used to RGB = 123
% * 'dynamicSurface' , Not implmented yet, future plan: boolean, false: static surface
% * 'detectorType'   , name of detector class, e.g. radialEdgeDetector
%                        (user thresholded), fastCylinderDetector
% * 'fitterType'     , name of fitter class, here spherelikeFitter
% * 'fitTime'        , time point used for fit
%
% The following project metadata information is optional 
%
% * 'description'     , string describing the data set set experiments metadata, 
%                                such as a description, and if the surface is dynamic,
%                                or requires drift correction of the sample.
% * 'jitterCorrection', Boolean, false: No fft based jitter correction 
%
expMeta                  = struct();
expMeta.channelsUsed     = 1;
expMeta.channelColor     = 1;
expMeta.description      = 'GAP43-mCherry labeled fruitfly embryo';
expMeta.dynamicSurface   = 0;
expMeta.jitterCorrection = 0; % 1: Correct for sample translation
expMeta.fitTime          = fileMeta.timePoints(1); 
expMeta.detectorType     = 'surfaceDetection.fastCylinderDetector';%
expMeta.fitterType       = 'surfaceFitting.spherelikeFitter'; 
%%
% Now set the meta data in the experiment.
%
xp.setFileMeta(fileMeta);
xp.setExpMeta(expMeta);
%%
%
% For a new project call initNew()
% This reads the stack size from the first available time point, then
% initialize fitter and detector and create fitOptions and detectOptions
% based on their defaults.
%
xp.initNew();

%% Load data for surface detection and rescale to unit aspect ratio
%
% Load a timepoint from the data.
% loadTime sets currentTime, loads the stack and resets detector and fitter
% with the appropriate options for that time. Optionally rescale the stack 
% to unit aspect ratio if desired.
%
xp.loadTime(5);

%%
%
% The stack used in this tutorial already has unit aspect ratio,
% therefore rescaleStackToUnitAspect will do nothing. However, most
% datasets won't, which will cause detectors to work less well. Then 
% rescaling axes to unit aspect ratio will be useful.
%
xp.rescaleStackToUnitAspect();

%% Detect the surface
%
% Before detecting the surface we set the options of the used detector.
% Changing the detectOptions resets the detector so that one cannot have
% for example a detected pointcloud that was detected with other parameters
% than the current one. 
%
% Default mode of this tutorial is the fast cylinder detector. Below are
% the options to set for the radialEdgeDetector. Don't forget to switch the
% detector type in the expMeta structure, when trying this out.
%
%% 
%
% options for fast cylinder detector
% Standard filters are applied to the images for point cloud detection. The
% options used are 
%
% * 'channel'             , vector of integers, specify which channel(s) to use for detection 
% * 'sigma'               , standard deviation of a gaussian filter in pixels 
% * 'ssfactor'            , integer, specifying the degree of data subsampling. 
% * 'nBins'               , number of radial bins to determine the point cloud in
% * 'rmRadialOutliers'    , remove radial outliers, low means stringent, high sloppy, 0 is off.
% * 'rmIntensityOutliers' , remove outliers based on intensity. 
% * 'zDim'                , specify long axis in data.
%
%%
%
% Comment below, when using radial edge detector 
% for full data sample change ssfactor to 4;
myDetectOpts = struct('channel', 1, 'sigma', 1, 'ssfactor', 2, ...
    'nBins', 120,'rmRadialOutliers', 2, 'rmIntensityOutliers',2,...
    'zDim', 2);  
%%
%
% options for radial edge detector
%
% * 'thresh' , determine cutoff for foreground detection
% * 'amin'   , minimal number of connected pixels to be foreground
% * 'bgdisc' , diameter of rolling ball filter for background estimation
% * 'dildisc', diameter of dilation disc for morphological opening and closing.
% * 'sp'     , permute stack dimensions before detection.
%
%%
%
%  Uncomment below, when using radial edge detector
%
%myDetectOpts = struct('channel', 1, 'sigma',1, 'ssfactor', 4,...
%   'rmRadialOutliers', 1.2, 'rmIntensityOutliers',2, 'thresh',75,...
%   'amin',5,'bgdisc',0,'dildisc',20,'sp',[1 3 2]);  

%%
%
% set the detect options in the project
%
xp.setDetectOptions(myDetectOpts);
%%
%
% calling detectSurface runs the surface detector and creates
% the detector.pointCloud object
%
xp.detectSurface();

%% Load a point cloud from disc
%
% if a point cloud was detected by other means, it can be loaded into the
% detector. Here we assume it was stored as a NPoints x 3 matrix, called
% points in matlab format. 
%
% uncomment the lines below, when loading a point cloud. 
%
% load('PointCloud.mat')
%
% pc = surfaceDetection.PointCloud(points);
%
% xp.detector.setPointCloud(pc);
%
% xp.detector.pointCloud.ROI.setRanges([1 ceil(max(points(:,1)))],... 
%     [1 ceil(max(points(:,2)))],[1 ceil(max(points(:,3)))]);
%

%% Inspect the point cloud in a cross section
%
% inspect point cloud over a cross section inthe  data. Dimensions 
% are x,y or z and the value has to be within the corresponding axis range.
%
inspectOptions= struct('dimension', 'z', 'value', 200, 'pointCloud', 'b');
xp.detector.inspectQuality(inspectOptions, xp.stack);

%% optionally save the quality inspection to disc
%
%storeQualityOptions = struct('inspectOptions',inspectOptions,...
% 'range',1:2:200,'outName',fullfile(projectDir,...
% 'debugOutput/InspectionDetector.tif'),'export','true','closeFig','true');
%xp.detector.storeQualityInspection(storeQualityOptions, xp.stack)

%% Inspect pointcloud in 3d
%
% Plot the points on the surface as a three dimensional point cloud.
% The subsampling factor reduces the number of points shown. 
%
ssfactor = 6;
xp.detector.pointCloud.inspect(ssfactor);

%% Fit the surface coarsly to prepare estimate of sample orientation
%
% We fit the surface in a cylindrical basis, with radius, eccentricity, centre 
% of mass and ellipse orientation as slowly varying polynomials of the z
% axis. 
%
% The fitter behaves similar to the detector, but the fitting of the
% surface is done by calling the function fitSurface() in the experiment
% class. This function calls fitter.fitSurface(detector.pointCloud) but
% also stores the resulting fittedParam in experiment class
%
% Initially we fit the surface with a coarse set of values, to obtain the
% orientation of the sample. Then we determine the frame with the samples
% symmetry axis and fit again with finer fit options.
% 
% The options of the fitter are 
% 
% * 'R'     , degree of polynomial fit for radius
% * 'X'     , degree of polynomial fit for centre of mass 1st coordinate
% * 'Y'     , degree of polynomial fit for centre of mass 2nd coordinate
% * 'e'     , degree of polynomial fit for ellipse eccentricity
% * 'phase' , degree of polynomial fit for ellipse orientation
% * 'path'  , directory name where to save debug output
%
fitOptions = struct('R',6,'X',4,'Y',4,'e',2,'phase',0,'path',...
                    fullfile(projectDir, 'debugOutput'));
xp.setFitOptions(fitOptions);
xp.fitSurface();
%
%% Determine sample orientation and perform fine fit
%
% Based on the carse fit, we determine the sample orientation, which we use 
% to rotate the detected point cloud, such that the samples long axis 
% coincides with the z-axis. Then we perform a second fit with finer chosen 
% fitOptions
%
xp.determineROIfromFit();
fitOptions = struct('R',8,'X',3,'Y',3,'e',3,'phase',0,'shift',0);
xp.setFitOptions(fitOptions);
xp.fitSurface();
%
%% Inspect the fit in a cross section
%
% inspect fit and point cloud over a cross section inthe  data. Dimensions 
% are x,y or z and the value has to be within the corresponding axis range.
%
zval = 211;
inspectOptions = struct('dimension','z','value',zval,'pointCloud','b');
xp.fitter.inspectQuality(inspectOptions, xp.detector, xp.stack);

%% Optionally save the quality inspection to disc
%
%storeQualityOptions = struct('inspectOptions',inspectOptions,...
% 'range',1:10,'outName',fullfile(projectDir,...
% 'debugOutput/InspectionFitter.tif'),'export','true','closeFig','true');
%xp.fitter.storeQualityInspection(storeQualityOptions, xp.detector, ...
%    xp.stack)

%% Inspect the bounding box on the sample
%
% draw the bounding box around the sample obtained by determineROIfromFit 
% in 3 different cross-sections through the centre of mass of the data.
%
scale = 100;
subplot(1, 3, 1);
imshow(xp.stack.getSlice('x', xp.currentROI.originp(1)));
view_sec = 'yz';
xp.detector.pointCloud.ROI.drawAxes(view_sec, scale);
xp.detector.pointCloud.ROI.drawROI(view_sec);
subplot(1, 3, 2);
imshow(xp.stack.getSlice('y', xp.currentROI.originp(2)));
view_sec = 'zx';
xp.detector.pointCloud.ROI.drawAxes(view_sec, scale);
xp.detector.pointCloud.ROI.drawROI(view_sec);
subplot(1, 3, 3);
imshow(xp.stack.getSlice('z', xp.currentROI.originp(3)));
view_sec = 'xy';
xp.detector.pointCloud.ROI.drawAxes(view_sec, scale);
xp.detector.pointCloud.ROI.drawROI(view_sec);

%% Normally evolve the fitted surface by shift in pixels
%
% Depending on the part of the dataset we want to extract, we may want to 
% displace the surface by shift in the direction of the surface normal.
% Positive values move inwards, negative values outwards.
shift = 15;
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
xp.fitter.setDesiredChart('cylinder1_proper',1);
xp.fitter.setDesiredChart('cylinder2_proper',0);
xp.fitter.setDesiredChart('polarLowerZ',0);
xp.fitter.setDesiredChart('polarUpperZ',0);
xp.fitter.setDesiredChart('anteriorEquidistant',1);
xp.fitter.setDesiredChart('posteriorEquidistant',1);
xp.generateSOI();
xp.SOI.NCalcInducedMetric(); % calculate the induced metric in all fundamental charts.

%% Pullback the stack to the desired charts
%
% Pass the region of interest and the current time to pull back the stack
% in the desired charts. This generates the data fields containing the
% pullback.

%% Pullback the stack to the desired charts
%
% Pass the region of interest and the current time to pull back the stack
% in the desired charts. This generates the data fields containing the
% pullback.
%

onionOpts = struct('nLayers', 3, 'layerDistance', 2, 'sigma', 20,'makeIP' ,'MIP');
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

%% Visualize the pullback on cylinder2
%
% Here we inspect the cylinder2 transform over the cylinder2_index patch.
%
patchName = 'cylinder2_index';
transformName = 'cylinder2';
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


%% measuring on the maps and in 3D

% Example 1: Length measurements
% Measure length of a line in cylinder 1 and posteriorEquidistant, and cross check with 3D. 
% Define a curve in the cylinder1 chart. 

% begin by loading the two pullbacks for visual comparisons later.
patchName     = 'posteriorEquidistant_index';
transformName = 'posteriorEquidistant';
pbPLZ2 = data.getPatch(patchName).getTransform(transformName).apply{1};

patchName     = 'cylinder1_index';
transformName = 'cylinder1';
pb = data.getPatch(patchName).getTransform(transformName).apply{1};


% Now define a curve. The format we assume Nx2 surface coordinates along curve

zBounds   = [840,940];
phiBounds = [200,400];
nPoints   = 100;
curveC1 = [linspace(zBounds(1),zBounds(2),nPoints);linspace(phiBounds(1),phiBounds(2),nPoints)]';
% curveC1 defines a line connecting the boundaries in z and phi by nPoints.
% For simplicity this is a straight line, but curved lines are possible, too. 


% Calculate the embedding of the curve in 3D, to measure the real length of
% the curve on the embryo.
eGrids = xp.SOI.embedding.getPatch(patchName).apply();
curve3D = [interp2(eGrids{1},curveC1(:,1),curveC1(:,2)),...
           interp2(eGrids{2},curveC1(:,1),curveC1(:,2)),...
           interp2(eGrids{3},curveC1(:,1),curveC1(:,2))];
       

% use transition maps to calculate the curve in a pole chart

cyl1Inv = xp.SOI.atlas.getChart('cylinder1').getInverse();
PLZ2Inv = xp.SOI.atlas.getChart('posteriorEquidistant').getInverse();
% 
% Calculate the transition map from cylinder 1 into the pole
tmap1PLZ = xp.SOI.atlas.getTransitionMap('cylinder1','posteriorEquidistant');
tmap1PLZ = PLZ2Inv.compose(tmap1PLZ).apply();

% now compute the curve in PLZ using tmap1PLZ; 
curvePLZ      = 0*curveC1;
curvePLZ(:,1) = size(pbPLZ2,1)-interp2(tmap1PLZ{1},curveC1(:,1),curveC1(:,2));
curvePLZ(:,2) = size(pbPLZ2,2)-interp2(tmap1PLZ{2},curveC1(:,1),curveC1(:,2));      
       
%% 
% Now we have the same curve in the cylinder1, in one of the pole maps and 
% in 3D. Lets inspect them

figure, 

% inspect the curve on cylinder 1
subplot(1,3,1)
imshow(pb,[])
hold on 
plot(curveC1(:,1),curveC1(:,2),'.')
% inspect the curve on the polar chart
subplot(1,3,2)
imshow(pbPLZ2,[])
hold on 
plot(curvePLZ(:,1),curvePLZ(:,2),'.')

% inspect the curve in 3D
subplot(1,3,3)
surf(double(eGrids{1}),double(eGrids{2}),double(eGrids{3}),double(pb))
axis equal, shading interp, colormap gray
hold on 
plot3(curve3D(:,1),curve3D(:,2),curve3D(:,3),'.','LineWidth',5), axis equal
view([-166 28])



%% Calculate curve lengths; 

       
% Now calculate the derivative of the curve to get the step size.       
curve3Dshift = circshift(curve3D, [1 0]);
dcurve3D = curve3Dshift - curve3D;
% if the curve is not closed, the first entry will be the
% distance between beginning and end, Remove this contribution
% by setting it to 0.
dcurve3D(1,:) = 0*dcurve3D(1,:);        
     
% the length in 3D is obtained by summing length of all the individual steps. 
l3D = sum(sqrt(sum(dcurve3D.^2,2)));


%%
% Measure the lengths in cylinder1 and in polar chart. Compare to
% 3D result

% first argument is time, second the curve, third the chart Name;
% the output is the correct length, taking curvature into account. 
lproperC1 = xp.SOI.properLength(xp.currentTime, curveC1, 'cylinder1');

% Compare now to uncorrected curve length in C1
% Just as we did the measurement in 3D, we simply add up step sizes in the
% 2D projection. 
curveC1shift = circshift(curveC1, [1 0]);
dcurveC1 = curveC1shift - curveC1;
% if the curve is not closed, the first entry will be the
% distance between beginning and end, Remove this contribution
% by setting it to 0.
dcurveC1(1,:) = 0*dcurveC1(1,:);        
       
lC1 = sum(sqrt(sum(dcurveC1.^2,2)));

%% repeat the same steps for the pole chart


lproperPLZ = xp.SOI.properLength(xp.currentTime, curvePLZ, 'posteriorEquidistant');

curvePLZshift = circshift(curvePLZ, [1 0]);
dcurvePLZ = curvePLZshift - curvePLZ;
% if the curve is not closed, the first entry will be the
% distance between beginning and end, Remove this contribution
% by setting it to 0.
dcurvePLZ(1,:) = 0*dcurvePLZ(1,:);        
       
lPLZ = sum(sqrt(sum(dcurvePLZ.^2,2)));

% compare your measurement results;
fprintf('\n')
disp('Measruement results');
fprintf('%f curve length in cylinder1 before metric correction compared to 3D. \n', lC1./l3D)
fprintf('%f curve length in cylinder1 after metric correction compared to 3D. \n', lproperC1/l3D)
fprintf('%f curve length in posteriorEquidistant before metric correction compared to 3D. \n', lPLZ/l3D)
fprintf('%f curve length in posteriorEquidistant after metric correction compared to 3D. \n', lproperPLZ/l3D)




%%  Area

% Example 2: Area measurements
% Measure some area in 2 charts, and compare to distortion corrected
% results; 

patchNamePole     = 'anteriorEquidistant_index';
transformNamePole = 'anteriorEquidistant';
pbPUZ2 = data.getPatch(patchNamePole).getTransform(transformNamePole).apply{1};

% first define a polygon in the cylinder1 chart. 

patchName     = 'cylinder1_index';
transformName = 'cylinder1';
pb = data.getPatch(patchName).getTransform(transformName).apply{1};

% We assume the area is described by a polyon stored as a Nx2 matrix.

% This is an example polyon. 
polyC1 = [588   430
          589   542
          493   541
          413   441
          484   361];
polyC1(:,1) = polyC1(:,1) -200;      
%%      
% To create a custom polygon, use matlabs impoly. Make sure the polygon is 
% not self intersecting in any of the charts.    

imshow(pb,[])
h = impoly;
polyC1 = wait(h);

%% Refine the polygon to get a more accurate description of the enclosed area; 
% Accurate description of the boundary of the area is particularly
% important when using transition maps from one chart into the next as
% lines in general get deformed under the transition map.
CoM = mean(polyC1);
x = polyC1(:,1)-CoM(1);
y = polyC1(:,2)-CoM(2);
phi = atan2(y,x);
r = sqrt(x.^2+y.^2);
phiDense = linspace(-pi,pi,50);
rDense = interp1([phi-2*pi;phi;phi+2*pi],[r;r;r],phiDense,'linear');
x = rDense.*cos(phiDense)+CoM(1);
y = rDense.*sin(phiDense)+CoM(2);
polyC1 = [x',y'];

%%
% Just as in case of length measurements,  use transition maps, to inspect 
% the same polygon in one of the polar maps.

% Calculate the transition map from cylinder 1 into the pole
tmap1PUZ = xp.SOI.atlas.getTransitionMap('cylinder1',transformNamePole);
PUZ2Inv  = xp.SOI.atlas.getChart(transformNamePole).getInverse();
tmap1PUZ = PUZ2Inv.compose(tmap1PUZ).apply();

% now compute the curve in PLZ using tmap1PLZ; 
polyPLZ      = 0*polyC1;
polyPLZ(:,1) = interp2(tmap1PUZ{1},polyC1(:,1),polyC1(:,2));
polyPLZ(:,2) = interp2(tmap1PUZ{2},polyC1(:,1),polyC1(:,2)); 

%% now inspect the polygon in both maps

figure, 

subplot(1,2,1)
imshow(pb,[])
hold on 
plot(polyC1([1:end,1],1),polyC1([1:end,1],2),'.-')

subplot(1,2,2)
imshow(pbPUZ2,[])
hold on 
plot(polyPLZ([1:end,1],1),polyPLZ([1:end,1],2),'.-')

%% now calculate the area under the polyon uncorrected 

[X,Y] = meshgrid(1:size(pb,2),1:size(pb,1));
mask = inpolygon(X,Y,polyC1(:,1),polyC1(:,2));
AC1 = sum(mask(:));

% and corrected;
AproperC1 = xp.SOI.properArea(xp.currentTime, mask, 'cylinder1');

% repeat the same measurements for the pole map

[X,Y] = meshgrid(1:size(pbPLZ2,2),1:size(pbPLZ2,1));
mask = inpolygon(X,Y,polyPLZ(:,1),polyPLZ(:,2));
APLZ = sum(mask(:));

% and corrected;
AproperPLZ = xp.SOI.properArea(xp.currentTime, mask, transformNamePole);

% compare the outcomes: 

fprintf('\n')
disp('Measruement results');
fprintf('%f polygon area in cylinder1 before metric correction. \n', AC1)
fprintf('%f polygon area in cylinder1 after metric correction. \n', AproperC1)
fprintf('%f polygon area in posteriorEquidistant before metric correction. \n', APLZ)
fprintf('%f polygon area in posteriorEquidistant after metric correction. \n', AproperPLZ)

AproperC1/AproperPLZ
