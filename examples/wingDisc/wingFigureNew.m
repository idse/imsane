clear all;
close all;

%[testScriptPath,~,ext] = fileparts(matlab.desktop.editor.getActiveFilename);
%run(fullfile(testScriptPath, 'header'))

projectDir = '/Users/idse/WingNew';

% Start by creating an experiment object, optionally pass on the project
% directory (otherwise it will ask)
dataDir = projectDir;
xp = project.Experiment(projectDir, dataDir);

%%
%--------------------------------------------------------------------------
% initialize new project
%--------------------------------------------------------------------------

% If the experiment was stored, it will load. Otherwise one has to set file
% metadata and experiment metadata first. 
% One can either pass the appropriate structures to their setter or call
% the setters without arguments and be asked (NOT implemented yet)

fileMeta = struct();
fileMeta.dataDir = '/Users/idse/Dropbox/flows/flows_shared/data/summerschool Eaton/';
fileMeta.filenameFormat = '2013.08.04_fixed_2hAPF.lif - Series005_cropped_t%u.ome.tif';
fileMeta.timePoints = [0];
fileMeta.stackDimPermute = [1 2 3];
%fileMeta.stackResolution = [0.1 0.1 1]; % stack resolution in micron
fileMeta.stackResolution = [0.259 0.259 0.755]; % stack resolution in micron

expMeta = struct();
expMeta.description = 'Qbio school Eaton group, fixed 2hAPF';
expMeta.channelsUsed = [1 2];
expMeta.channelColor = [2 1 3];
expMeta.dynamicSurface = true;
expMeta.jitterCorrection = false;

expMeta.detectorType = 'surfaceDetection.planarEdgeDetector';
expMeta.fitterType = 'surfaceFitting.tpsFitter';

xp.setFileMeta(fileMeta);
xp.setExpMeta(expMeta);

% For a new project one also has to call initNew()
% This reads the stack size from the first available time point, then
% initialize fitter and detector and create fitOptions and detectOptions
% based on their defaults.
xp.initNew();

%%
%--------------------------------------------------------------------------
% load data
%--------------------------------------------------------------------------
% Now we can load a timepoint from the data.
% loadTime sets currentTime, loads the stack and resets detector and fitter
% with the appropriate options for that time.

xp.loadTime(0);

xp.rescaleStackToUnitAspect();

%% create a mask for the data

% dimensions = 3;
% xp.stack.generateProjectionMask(1);

projectionMask = cell([3 1]);

% imageSize is xyz, but in projectionMask we're using matlab yxz indexing
projectionMask{2} = true([xp.stack.imageSize(1) xp.stack.imageSize(3)]);
size(projectionMask{2})
projectionMask{2}(:, 130:end) = false;
size(projectionMask{2})

xp.stack.setProjectionMask(projectionMask);

imshow(projectionMask{2},[])

%%
%--------------------------------------------------------------------------
% detect the surface
%--------------------------------------------------------------------------
% Before detecting the surface we can set the right detect options.
% Changing the detectOptions resets the detector so that one cannot have
% for example a detected pointcloud that was detected with other parameters
% than the current one.

% detectOpts = struct('sigma', 3.5, 'channels', 1,...
%             'sigZthresh', 0.5, 'sigZscale', 3,...
%             'maxIthresh', 0.02, 'summedIthresh', 0,...
%             'sigZoutliers', 1, 'scaleZoutliers', 3, 'zdir', -3); 

detectOptions = struct(  'sigma', 3.5, 'channels', 1, 'zdir', -3,...
                        'maxIthresh', 0.02, 'summedIthresh', 0,...
                        'sigZoutliers', 1, 'scaleZoutliers', 3); 
                    
% calling detectSurface runs the surface detector and creates
% the detector.pointCloud object

xp.setDetectOptions(detectOptions);
xp.detectSurface();

imshow(xp.detector.mask.*xp.detector.surfaceMatrix, [])

%% to find good parameters without redetecting the surface

xp.detector.resetMask();
myDetectOpts = struct('sigma', 3.5, 'channels', 1,...
            'sigZthresh', 0.1, 'sigZscale', 3,...
            'maxIthresh', 0.1, 'summedIthresh', 0,...
            'sigZoutliers', 2, 'scaleZoutliers', 30, 'zdir', -3); 

xp.detector.setOptions(myDetectOpts);    
xp.detector.applyMasks();

imshow(xp.detector.mask.*xp.detector.surfaceMatrix, [])

%%
% inspect point cloud cross section over data
inspectOptions= struct('dimension', 'y', 'value', 730, 'pointCloud', 'b');
xp.detector.inspectQuality(inspectOptions, xp.stack);

%%
% inspect point cloud cross section over data
inspectOptions= struct('dimension', 'x', 'value', 250, 'pointCloud', 'b');
xp.detector.inspectQuality(inspectOptions, xp.stack);

%% 
% inspect pointcloud in 3d
ssfactor = 10;
figure, xp.detector.pointCloud.inspect(ssfactor);

%%
%--------------------------------------------------------------------------
% fit the surface for the columnar cells
%--------------------------------------------------------------------------

% make a ROI without aligning the pointCloud using determineRanges
margin = 0;
xp.detector.pointCloud.determineRanges(margin);
xp.updateROI();

%%

% then fit
fitOptions = struct('smoothing', 0, 'gridSize', [100 100]);
xp.setFitOptions(fitOptions);
xp.fitSurface();

%%
inspectOptions= struct('dimension', 'y', 'value', 700, 'pointCloud','b');
xp.fitter.inspectQuality(inspectOptions, xp.detector, xp.stack);

%%
% normally evolve (same code as above to inspect)
% takes the fitted surface and normally evolves it inward (or outward with
% negative numbers)
% subsampling is useful to speed up normal evolution of fullsize fit 

shift = 14;
xp.zEvolve(shift);

%%
% now generate the surface of interest
% the charts to be generated are specified in xp.fitter.charts

xp.generateSOI();

% generate pullbacks
xp.SOI.pullbackStack(xp.stack, xp.currentROI, xp.currentTime);
%xp.save();

%%
imstack = xp.stack.image.apply{1}(xp.ROI.ypRange(1):xp.ROI.ypRange(2),...
                                  xp.ROI.xpRange(1):xp.ROI.xpRange(2),...
                                  xp.ROI.zpRange(1):xp.ROI.zpRange(2));
imshow(max(imstack(1:800, 1:800, :), [], 3), [])

%%
columnar = xp.SOI.data(1).patches{1}.getTransform('xy').apply{1}(1:800, 1:800);
Zcolumn = double(xp.fitter.fittedPoints{3}(1:800, 1:800));
imshow(cat(3, mat2gray(max(imstack(1:800, 1:800, :), [], 3)),...
        mat2gray(columnar), 0*mat2gray(Zcolumn)), [])

%%
% WRITE columnar data and embedding

imwrite(columnar, 'columnar.tif');

% writing floats to tiff is a bitch
t = Tiff('columnEmb.tif', 'w'); 
tagstruct.ImageLength = size(Zcolumn, 1); 
tagstruct.ImageWidth = size(Zcolumn, 2); 
tagstruct.Compression = Tiff.Compression.None; 
tagstruct.SampleFormat = Tiff.SampleFormat.IEEEFP; 
tagstruct.Photometric = Tiff.Photometric.MinIsBlack; 
tagstruct.BitsPerSample = 32; %info.BitsPerSample; % 32; 
tagstruct.SamplesPerPixel = 1; %info.SamplesPerPixel; % 1; 
tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky; 
t.setTag(tagstruct); 
t.write(single(Zcolumn)); 
t.close();

%%
%--------------------------------------------------------------------------
% fit the surface for the peripodial cells
%--------------------------------------------------------------------------

% then fit
fitOptions = struct('smoothing', 1000, 'gridSize', [100 100]);
xp.setFitOptions(fitOptions);
xp.fitSurface();

%%

inspectOptions= struct('dimension', 'x', 'value', 550, 'pointCloud','b');
xp.fitter.inspectQuality(inspectOptions, xp.detector, xp.stack);

%%
inspectOptions= struct('dimension', 'y', 'value', 450, 'pointCloud','b');
xp.fitter.inspectQuality(inspectOptions, xp.detector, xp.stack);

%%

% compute the Gaussian and mean curvature
X = xp.fitter.fittedPoints;
sigma = 10;
[K,H,Pmax,Pmin] = mySurfature(X{1},X{2},X{3}, sigma);
imshow(H, [])

%%

% create a mask that excludes high mean curvature to get rid of the folds
highCurv = zeros(size(H));
highCurv(10:(end-10), 10:(end-10)) = mat2gray(abs(H(10:(end-10), 10:(end-10)))) > 0.25;

curvMask = true(size(xp.detector.mask));
curvMask(   xp.ROI.ypRange(1):xp.ROI.ypRange(2),...
            xp.ROI.xpRange(1):xp.ROI.xpRange(2)) = ~highCurv;
curvMask = imerode(curvMask, strel('disk', 45));
imshow(curvMask);

%%
%curvMask = false(size(xp.detector.mask));
xp.detector.setManualMask(curvMask);
xp.detector.applyMasks();

%%
% inspect point cloud cross section over data
inspectOptions= struct('dimension', 'y', 'value', 680, 'pointCloud', 'b');
xp.detector.inspectQuality(inspectOptions, xp.stack);

%%
% then fit again (slightly higher smoothing)
fitOptions = struct('smoothing', 2000, 'gridSize', [100 100]);
xp.setFitOptions(fitOptions);
xp.fitSurface();

%%

inspectOptions= struct('dimension', 'y', 'value', 800, 'pointCloud','b');
xp.fitter.inspectQuality(inspectOptions, xp.detector, xp.stack);


%% 3D inspection of fit

xp.fitter.inspectTPS();

%%
% normally evolve (same code as above to inspect)
% takes the fitted surface and normally evolves it inward (or outward with
% negative numbers)
% subsampling is useful to speed up normal evolution of fullsize fit 

shift = -2;
xp.zEvolve(shift);

%%
inspectOptions= struct('dimension', 'y', 'value', 700, 'pointCloud','b');
xp.fitter.inspectQuality(inspectOptions, xp.detector, xp.stack);

%%
% now generate the surface of interest
% the charts to be generated are specified in xp.fitter.charts

xp.generateSOI();

% generate pullbacks
xp.SOI.pullbackStack(xp.stack, xp.currentROI, xp.currentTime);
%xp.save();

%%

perip = xp.SOI.data(1).patches{1}.getTransform('xy').apply{1}(1:800, 1:800);
imshow(perip, []);           
%%
imwrite(perip, 'peripodial.tif');


%%

Zperip = double(xp.fitter.embedding.apply{3}(1:800, 1:800));
imshow(Zperip, []);

% writing floats to tiff is a bitch
t = Tiff('peripEmb.tif', 'w'); 
tagstruct.ImageLength = size(Zcolumn, 1); 
tagstruct.ImageWidth = size(Zcolumn, 2); 
tagstruct.Compression = Tiff.Compression.None; 
tagstruct.SampleFormat = Tiff.SampleFormat.IEEEFP; 
tagstruct.Photometric = Tiff.Photometric.MinIsBlack; 
tagstruct.BitsPerSample = 32; %info.BitsPerSample; % 32; 
tagstruct.SamplesPerPixel = 1; %info.SamplesPerPixel; % 1; 
tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky; 
t.setTag(tagstruct); 
t.write(single(Zperip)); 
t.close();

%%

MIP = max(xp.stack.image.apply{1}(:,:,1:100), [], 3);
SIP = sum(xp.stack.image.apply{1}(:,:,1:100), 3);
MIP = MIP(1:800, 1:800);
SIP = SIP(1:800, 1:800);

%%
imshow(MIP, [])
imwrite(MIP, 'MIP.tif');
imwrite(uint16(SIP), 'SIP.tif');

