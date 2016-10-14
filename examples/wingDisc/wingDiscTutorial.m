%% ImSAnE Wing Disc Tutorial
%
% In this tutorial we detect the apical surface of the disc proper cells
% and the peripodial cells on one side of a wing disc two hours after
% puparium formation. We read out actin and E-cadherin on these surfaces.
%
% This an example of how to detect and fit what we call planar surfaces,
% where the surface can be described as a height for every x,y.
% It also demonstrates using multiple channels and making multiple
% surfaces. Finally, because finding folds as regions of high curvature is
% part of detecting the peripodial surface, this example shows how to
% calculate the surface metric and curvature.
%
% Note that ImSAnE is fully documented and additional information about
% available properties, methods and options can be found using the matlab
% documentation system. 
% 
% For example type:

%%
doc surfaceDetection.planarDetector

%% Initialize the project
%
% We start by creating an experiment object, which holds this metadata and 
% provides a frontend for a number of tasks such as loading the raw data.
% To construct the experiment object we need to pass dataDir, the directory 
% containing the raw data and projectDir, the directory where results of 
% the script will be saved.

clear all; close all;

[scriptPath,~,ext] = fileparts(matlab.desktop.editor.getActiveFilename);

dataDir = fullfile(scriptPath, 'rawData');
projectDir = fullfile(scriptPath, 'projectFiles');

xp = project.Experiment(projectDir, dataDir);

%%
% Next we set the metadata pertaining to the raw data files in the structure
% fileMeta. ImSAnE assumes that timepoints are saved individually and that 
% filenames in a timeseries are identical up to an integer specifying the 
% timepoint. Therefore we have
%
% * filenameFormat:               Filename, %u in the position of the integer
% * timePoints = [t1, t2, ..] :   List of times available. In this example 
% we have just a single time point 0.
% * stackResolution :             Stack resolution in micron.

fileMeta = struct();
fileMeta.dataDir = dataDir;
fileMeta.filenameFormat = 'wing_fixed_2hAP.ome.tif';
fileMeta.timePoints = [0];
fileMeta.stackResolution = [0.26 0.26 0.76]; 
fileMeta.rawStackSize = [800 800 50 2];
fileMeta.nChannels = 2;

%% 
% In the structure expMeta we set general parameters for the surface
% analysis we will do. 
%
% * channelsUsed:   Which channels do we need.
% * channelColor:   Assign color to the channels, RGB = [1 2 3].
%                   In this example the first channel is E-cadherin, and 
%                   the second is actin. We want these in green and red,
%                   respectively.
% * dynamicSurface: Does the surface shape change with time?  
%                   For a single time point this is false. True is not yet
%                   supported.
% * jitterCorrection:   Not needed here.
% * detectorType:       Which type of surface detector will be used.
%                       We will look at only one side of the wing, which is
%                       a planar surface so we use planarDetector.
% * fitterType:         Which type of fitter will be used.
%                       We fit planar surfaces using Thin Plate Spline:
%                       tpsFitter.

expMeta = struct();
expMeta.description = 'Fixed wing 2hAPF [Ecad, actin, widefield]';
expMeta.channelsUsed = [1 2];
expMeta.channelColor = [2 1 3];
expMeta.dynamicSurface = false;
expMeta.jitterCorrection = false;
expMeta.detectorType = 'surfaceDetection.planarEdgeDetector';
expMeta.fitterType = 'surfaceFitting.tpsFitter';

xp.setFileMeta(fileMeta);
xp.setExpMeta(expMeta);

%%
% Finally we call initNew(), which reads the stack size from the first 
% available time point, then initializes fitter and detector and creates 
% fitOptions and detectOptions based on their defaults.

xp.initNew();

%% Load a time point from the data
%
% Now that we have set up the project, we can load a time point.
% loadTime sets xp.currentTime, loads the stack into xp.stack 
% and resets detector and fitter with the default options for that time.
%
% We rescale to unit aspect ratio to make detection work better later and
% visualize in the right proportions.

xp.loadTime(0);
xp.rescaleStackToUnitAspect();

%% 
% xp.stack is not an array but a Stack object.
% The easy way to look at a slice through the data is using getSlice.

imshow(xp.stack.getSlice('x', 300), []);

%% Mask data
%
% Detection speed depends on stack size. We can speed it up by some basic
% masking. Masking out unwanted information can also improve detection and
% fitting.
% Here this is not really necessary and we only do it to demonstrate the 
% functionality, as this dataset has already been cropped to meet maximal
% download size requirements.
% Rather than creating 3D masks, we create 2D masks that are applied to
% each cross section. In this case we simply cut out the bottom half of the
% stack because we want to detect the top of the wing. We are doing this in
% cross section normal to the y direction. 

projectionMask = cell([3 1]);
projectionMask{2} = true([xp.stack.imageSize(1) xp.stack.imageSize(3)]);
projectionMask{2}(:, 130:end) = false;

xp.stack.setProjectionMask(projectionMask);

imshow(xp.stack.getSlice('y', 500), []);

%% Detect the surface
%
% planarDetector.detectSurface detects the surface as the position of the 
% maximal Gaussian z-derivative in some direction, i.e. the position of the
% largest intensity jump along some direction and smoothened over some
% scale.
%
% A number of detection options directly affect detection:
%
% * sigma :     Width of the Gaussian z-derivative.
% * channels :  Channels (summed) to use for detection.
% * zdir :      Dimension corresponding to z, minus flips direction.
% Flipping the direction can sometimes improve detection.
%
% Then there are options which filter the result and can be modified
% without redetecting:
%
% * maxIthresh:     Throw out points with MIP dimmer than this.
% * summedIthresh:  Throw out points with SIP dimmer than this.
% * sigZoutliers:   Remove height outliers after all other masks.
% * scaleZoutliers: Spatial scale of outlier removal.
%
% scaleZoutliers is the linear size of a region over which the
% distribution of height is computed, sigZoutliers is then a cutoff in
% units of standard deviation of this distribution to remove misdetected
% points far above or below the other points in the region.

detectOptions = xp.detector.defaultOptions;
detectOptions.sigma = 3.5;
detectOptions.zdir = -3;
detectOptions.maxIthresh = 0.02;  

% Calling detectSurface runs the surface detector and creates the point
% cloud in detector.pointCloud.

xp.setDetectOptions(detectOptions);
xp.detectSurface();

%%
% Different from the other detectors, the detected surface is
% represented not only by a PointCloud object but also by an image
% surfaceMatrix, containing z values for each xy.
% Looking at this height map masked by the filters specified in
% detectOptions one can judge how well the surface was detected.

imshow(xp.detector.mask.*xp.detector.surfaceMatrix, [],...
                                            'InitialMagnification', 40);

%% 
% One can then find better filter parameters without redetecting the
% surface by changing the second block of options in detectOptions and 
% calling resetMask and applyMasks. 

xp.detector.resetMask();

detectOptions.maxIthresh = 0.1; 
detectOptions.sigZoutliers = 2; 
detectOptions.scaleZoutliers = 30; 

xp.detector.setOptions(detectOptions);    
xp.detector.applyMasks();

imshow(xp.detector.mask.*xp.detector.surfaceMatrix, [],...
                                            'InitialMagnification', 40);

%%
% We can also inspect a point cloud cross section over the data with
% detector.inspectQuality. In the pointCloud option, 'c' specifies the 
% color cyan.

inspectOptions= struct('dimension', 'x', 'value', 620, 'pointCloud', 'c');
xp.detector.inspectQuality(inspectOptions, xp.stack);

%% 
% Or we can look at the point cloud in 3d, with some subsampling factor.
ssfactor = 50;
xp.detector.pointCloud.inspect(ssfactor);

%% Fit the surface for the disc proper cells
%
% By detecting the largest intensity jump along z for each x,y in the
% E-cad channel and filtering out local outliers we have found the apical
% surface of the disc proper cells. We can now fit a smooth surface
% representation to that.
%
% tpsFitter fits the pointcloud using a thin plate spline fit. It has the
% following options:
%
% * gridSize:     Size of grid on which to generate fitted surface
%               default [50 50], full size takes long.
% * smoothing:    TPS smoothing parameter (default 1000).

fitOptions = struct('smoothing', 500, 'gridSize', [100 100]);
xp.setFitOptions(fitOptions);
xp.fitSurface();

%%
% We can visualize the result on a cross section with
% fitter.inspectQuality.

xp.fitter.inspectQuality(inspectOptions, xp.detector, xp.stack);

%%
% The detector picks up the edge of the E-cad signal but the best read out
% goes solidly through it so we want to move the surface down a little. For
% this we use zEvolve, with a shift specified in pixels.

shift = 12;
xp.zEvolve(shift);

xp.fitter.inspectQuality(inspectOptions, xp.detector, xp.stack);

%%
% We now generate the Surface Of Interest. The charts to be generated are 
% specified in xp.fitter.charts. In this case there is only one, called
% 'xy'. 

xp.generateSOI();

%% Pull back the data to the surface
% 
% We pull back the data to the SOI using pullbackStack.

xp.SOI.pullbackStack(xp.stack, xp.currentROI, xp.currentTime);

%%
% To look at the pullback, we call the data field of the SOI at the right
% time, and get a particular patch from that with getPatch. A patch is a 
% part of a surface. In this case, there is only one called xy_index.
% Then we get the data in some patch in a particular coordinate system with
% getTransform. In this case there is only one coordinate system: xy.
% What we get is an object not only holding the image data but also
% metadata and methods to manipulate it. The actual data is obtained by
% calling the method apply. This returns a cell array with entries for each
% channel.

% xp.tIdx converts the time into an index in a list of time points
tidx = xp.tIdx(0);

% the first channel is Ecad
channel = 1;

discProperPatch = xp.SOI.data(tidx).getPatch('xy_index');
discProperImage = discProperPatch.getTransform('xy').apply{channel};
figure, imshow(discProperImage, [], 'InitialMagnification', 50);

%% Save the result
%
% Finally we save the SOI using SOI.save. We set the following options:
%
% * dir:            The directory to save the SOI to.
% * imwriteOptions: Pullbacks are saved to image files using imwrite, we
% can pass options to change file format, compression etc. For example we
% could change this option to
% imwriteOptions = {'jp2', 'Mode', 'lossless'}; 
% * make8bit:       Often absolute intensities don't matter and 8 bit offers
% a large enough dynamic range. This options rescales the lookup table and
% converts to 8 bit before saving.

imwriteOptions = {'tif', 'Compression', 'deflate'};
savedir = fullfile(scriptPath, 'discProperApicalSOI');

options = struct(   'dir',              savedir,...
                    'imwriteOptions',   {imwriteOptions},...
                    'make8bit',         true);
xp.SOI.save(options)

%%
% All metadata is saved in SOI.xml. Pullbacks, geometry and any other data
% defined on the surface are saved to image files in subdirectories 
% reflecting the structure of patches and coordinates explained earlier in 
% this tutorial. We can reload a surface of interest with
% SOI.load(directory)


%% Get the peripodial cells
%
% We also want to find the peripodial surface. Since these are lying on top
% of the disc proper cells except in the folds, the strategy is to detect
% the folds as regions of high curvature in the disc proper surface, mask
% out the folds from the detected point cloud and then fit a new surface
% without folds that when moved up captures the peripodial cells.

% first store columar SOI before we overwrite it with peripodial
columnarSOI = xp.SOI;

%%
% We start by computing the mean curvature (see supplementary text).
% This requires computing the metric first. We smoothen
% with a Gaussian filter on the scale of cells (about 10 pixels).

xp.SOI.NCalcInducedMetric('xy');
xp.SOI.NCalcCurvature('xy');

H = xp.SOI.getField('curvature').getPatch('xy_index').trace().apply{1};

sigma = 10;
H = mat2gray(imfilter(H, fspecial('gaussian', 3*sigma, sigma)));
imshow(H, [], 'InitialMagnification', 20)

%%
% By thresholding the mean curvature we create a mask that excludes high 
% mean curvature to get rid of the folds.

highCurv =  H < 0.4;
curvMask = ~highCurv;
curvMask = imerode(curvMask, strel('disk', 45));

xp.detector.setManualMask(curvMask);
xp.detector.applyMasks();

%%
% Looking at the fold-masked point cloud in a cross section we see that the
% mask works well.

xp.detector.inspectQuality(inspectOptions, xp.stack);

%%
% Fitting to this masked point cloud with we set smoothing higher because
% we are fitting a smoother surface (the peripodial cells don't fold).

fitOptions = struct('smoothing', 2000, 'gridSize', [100 100]);
xp.setFitOptions(fitOptions);
xp.fitSurface();

%%
% As before we shift the surface, this time up.

shift = -2;
xp.zEvolve(shift);

%%
inspectOptions= struct('dimension', 'x', 'value', 620, 'pointCloud', 'c');
xp.fitter.inspectQuality(inspectOptions, xp.detector, xp.stack);
hold on
plot(columnarSOI.embedding.patches{1}.apply{3}(:,inspectOptions.value),'r','LineWidth',2)
hold off

%%
% We generate the SOI again, pull back the data and look at it.

xp.generateSOI();
xp.SOI.pullbackStack(xp.stack, xp.currentROI, xp.currentTime);

perip = xp.SOI.data(1).patches{1}.getTransform('xy').apply{1};
figure, imshow(perip, [], 'InitialMagnification', 50);           

%%
pbcolor = cat(3, mat2gray(columnarSOI.data.patches{1}.apply{2}),mat2gray(columnarSOI.data.patches{1}.apply{1}),0*mat2gray(columnarSOI.data.patches{1}.apply{1}));
imshow(pbcolor,[])
imwrite(pbcolor, '/Users/idse/columnar.tif');

%%
pbcolor = cat(3, mat2gray(xp.SOI.data.patches{1}.apply{2}),mat2gray(xp.SOI.data.patches{1}.apply{1}),0*mat2gray(xp.SOI.data.patches{1}.apply{1}));
imshow(pbcolor,[])
imwrite(pbcolor, '/Users/idse/peripodial.tif');

%%
% Finally we save it to a different directory from before, because it is a
% different surface.

savedir = fullfile(scriptPath, 'peripodialSOI');
options = struct(   'dir',              savedir,...
                    'imwriteOptions',   {imwriteOptions},...
                    'make8bit',         true);
xp.SOI.save(options)


