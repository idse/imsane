%% ImSAnE: Drosophila melanogaster embryo tutorial
%
% In this tutorial we detect the embryos apical surface during
% cellularization, fit a spherelike surface, generate a surface of interest
% and pullback the image data to various charts.
%
% This an example of how to detect and fit what we call spherelike surfaces,
% where the surface can be described as a slowly varying function of a
% preferred axis, which we call the z-axis. It also demonstrates using
% batch processing of multiple time points with jitter correction. Moreover
% it shows how to choose from two different detectors and describes the 
% options of loading an externally provided point cloud.
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
%
% Setup a working directory for the project, where extracted surfaces,
% metadata and debugging output will be stored. Also specify the directory
% containing the data. 
%
[scriptPath,~,ext] = fileparts(matlab.desktop.editor.getActiveFilename);

%dataDir = fullfile(scriptPath, 'rawData');
dataDir = '/Users/idse/Dropbox/flows/flows_shared/data/zebrafish/huisken/video9/';
%'/Users/idse/Dropbox/flows/flows_shared/codes/ImSAnE_NatureMethods/examples/DrosoEmbryo/rawData/';
%projectDir = fullfile(scriptPath, 'projectFiles');
projectDir = dataDir;
%
% Start by creating an experiment object, optionally pass on the project
% directory (otherwise it will ask), and change into the directory of the
% data. This serves as a frontend for data loading, detection, fitting etc.
%
xp = project.Experiment(projectDir, dataDir);

% optionally switch working directory into the data directory.
%
cd(dataDir)

% Set file and experiment metadata
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
fileMeta.filenameFormat  = 'VSilent_T%d_G.tif';%'Time%06d_bin2.tif'; % for full data sample use Time000000.tif
fileMeta.filenameFormat  = 'VSilent_T%d_merged.ome.tif';
fileMeta.timePoints      = 1:129;%67:68;%1:129; % for full data sample use 0;
fileMeta.stackResolution = [.45 .45 1];%[.5 .5 .5]; 
fileMeta.swapZT          = 0; % for full data sample use 1;

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
expMeta.channelsUsed     = [1 2];
expMeta.channelColor     = [1 2];
expMeta.description      = 'Beating Zebrafish Heart';
expMeta.dynamicSurface   = 1;
expMeta.jitterCorrection = 0; % 1: Correct for sample translation
expMeta.fitTime          = 60; 
expMeta.detectorType     = 'surfaceDetection.fastCylinderDetector';%
%expMeta.fitterType       = 'surfaceFitting.cylinderMeshWrapper';
expMeta.fitterType       = 'surfaceFitting.meshWrapper'; 

% Now set the meta data in the experiment.
%
xp.setFileMeta(fileMeta);
xp.setExpMeta(expMeta);

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
xp.loadTime(fileMeta.timePoints(1)); %60

%
% The stack used in this tutorial already has unit aspect ratio,
% therefore rescaleStackToUnitAspect will do nothing. However, most
% datasets won't, which will cause detectors to work less well. Then 
% rescaling axes to unit aspect ratio will be useful.
%
xp.rescaleStackToUnitAspect();

%% Read surface mesh produced by meshlab

%outputMesh = ['/Users/idse/Dropbox/flows/flows_shared/data/zebrafish/huisken/video9/objsmooth2/heartPC' num2str(xp.currentTime) '_ilastik_mesh.ply'];
%outputMesh = ['/Users/idse/Dropbox/flows/flows_shared/data/zebrafish/huisken/video9/objIlastikNew/heartPC' num2str(xp.currentTime) '_ilastik_mesh.ply'];
%outputMesh = ['/Users/idse/Dropbox/flows/flows_shared/data/zebrafish/huisken/video9/objIlastikNewNewSmooth2/heartPC' num2str(xp.currentTime) '_ilastik_mesh.ply'];
%outputMeshFormat = '/Users/idse/Dropbox/flows/flows_shared/data/zebrafish/huisken/video9/objIlastikFinal/heartPC%d_ilastik_mesh.ply';
outputMeshFormat = '/Users/idse/Dropbox/flows/flows_shared/data/zebrafish/huisken/video9/objIlastikFinalSmooth2Mar23/heartPC%d_ilastik_mesh.ply';
outputMesh = sprintf(outputMeshFormat, xp.currentTime);
meshPermutation = [2,1,3]; % dumb stuff

tic
%mesh = readObj(outputMesh);
%mesh.f = mesh.f.v;
%mesh = read_obj_mod('/Users/idse/Dropbox/flows/flows_shared/data/zebrafish/huisken/video9/objsmooth2/heartPC1_ilastik_mesh.obj');
mesh = read_ply_mod(outputMesh);
toc
mesh.v = mesh.v(:,meshPermutation); 
mesh.vn = mesh.vn(:,meshPermutation);

% 
% C = zeros([size(mesh.v,1) 1]);
% trimesh(mesh.f, mesh.v(:,1), mesh.v(:,2), mesh.v(:,3), C, 'CDataMapping', 'direct')
% axis equal

%% set seeds points

seedsXinit = [127 338 106; 311 290 110; 332 326 252; 238 436 184; 127 205 273];
seedsXinit = [127 338 106; 311 290 110; 332 326 252; 307 431 185; 127 205 273];
%seedsXinit = [127 338 106; 311 290 110];
seeds = pointMatch(seedsXinit, mesh.v);

nSeeds = numel(seeds);
nVorSeeds = 3; %3
nDiskSeeds = nSeeds - nVorSeeds;

diskSeeds = seeds(nVorSeeds+1:end);
VorSeeds = seeds(1:nVorSeeds);
% round of improvement in the hope of removing cusps
% for i = 1:10
%     VorSeeds = perform_lloyd_mesh(mesh.v, mesh.f, VorSeeds);
% end

% % we need an extra row of zeros now
% diskSeeds = [diskSeeds, zeros(size(diskSeeds))];
% VorSeeds = [VorSeeds, zeros(size(VorSeeds))];

% keep track of seeds propagation
seedsX = seedsXinit;
seedsXall = {};
seedsXall{1} = seedsX;
seedsCounter = 2;

%%
% full atlas with overlap for everything else

fitOptions = struct('chartSeeds', VorSeeds, 'transitionWidth', 30,...
                    'diskSeeds', diskSeeds, 'diskRadius', 150, 'makeTMaps', true);
xp.setFitOptions(fitOptions);
xp.fitSurface(mesh);

%% visualize 

% whole mesh
xp.fitter.inspectMesh(1:3);

%% Inspect the fit in a cross section
%
% inspect fit and point cloud over a cross section inthe  data. Dimensions 
% are x,y or z and the value has to be within the corresponding axis range.
%
clf
zval = 200;
inspectOptions = struct('dimension','x','value',zval,'pointCloud','b', 'noalign', true);
xp.fitter.inspectQuality(inspectOptions, xp.detector, xp.stack);

%%
xp.fitter.setDesiredChart('conformal', true);
xp.fitter.setDesiredChart('exponential', false);
xp.fitter.setFitDomain(xp.stack.image.domain);
xp.generateSOI();

% %% DIRECT MESHUREMENT :)
% 
% n = mesh.vn;
% for i = 1:3
%     n(:,i) = n(:,i)./sqrt(sum(mesh.vn.^2,2));
% end
% meshArea(mesh.v + 20*n,mesh.f)
% meshArea(mesh.v,mesh.f)
% meshArea(mesh.v - 20*n,mesh.f)

%% USING TMAPS FOR PARTITION

cover = 1:nVorSeeds; 
margin = round(fitOptions.transitionWidth/2);
partitionMask = xp.SOI.partition(xp.currentTime, cover, margin);

% ti = 1;
% gti = ti;
% this = xp.SOI;
% 
% if isempty(this.g(gti).patches), this.NCalcInducedMetric(gti); end
% 
% A = 0;
% 
% for i = 1:nVorSeeds
% 
%     curChart = this.atlas(gti).charts{i};
%     gpatch = this.g(gti).patches{i};
%     detg = gpatch.determinant();
%     sqrtdetggrid = sqrt(detg.apply{1});
% 
%     A = A + sum(sqrtdetggrid(partitionMask{i}));
% end
% 
% A

%% Pullback the stack to the desired charts
%
% Pass the region of interest and the current time to pull back the stack
% in the desired charts. This generates the data fields containing the
% pullback.
%
% 21, 2, 20
onionOpts = struct('nLayers', 41, 'layerDistance', 1, 'sigma', 1,...
                    'makeIP', 'both', 'IPonly', true);
xp.SOI.pullbackStack(xp.stack, [], xp.currentTime, onionOpts);

%% SUM INTENSITY FROM ONION PATCHES NAIVE AND PROPER

ti = 1;
gti = ti;
this = xp.SOI;

totI1 = zeros([1 xp.SOI.nTimePoints]);
totI2 = zeros([1 xp.SOI.nTimePoints]);
totI1prop = zeros([1 xp.SOI.nTimePoints]);
totI2prop = zeros([1 xp.SOI.nTimePoints]);

partsI1 = zeros([nVorSeeds xp.SOI.nTimePoints]);
partsI2 = zeros([nVorSeeds xp.SOI.nTimePoints]);
partsI1prop = zeros([nVorSeeds xp.SOI.nTimePoints]);
partsI2prop = zeros([nVorSeeds xp.SOI.nTimePoints]);

totI13D = zeros([1 xp.SOI.nTimePoints]);
totI23D = zeros([1 xp.SOI.nTimePoints]);
A = zeros([1 xp.SOI.nTimePoints]);

if isempty(this.g(gti).patches), this.NCalcInducedMetric(gti); end

for i = 1:nVorSeeds

    X = xp.SOI.embedding(xp.currentTime).patches{i}.apply;
    geom = GaussGeometry(X{1},X{2},X{3},0);
    sqrtdetggrid = geom.dA;
%     curChart = this.atlas(gti).charts{i};
%     gpatch = this.g(gti).patches{i};
%     detg = gpatch.determinant();
%     sqrtdetggrid = sqrt(detg.apply{1});
    
    SIP = this.getField('data_SIP');
    for j = 1:2
        SIPgrids{j} = partitionMask{i}.*SIP(ti).patches{i}.apply{j};
    end

    partsI1(i,ti) = sum(SIPgrids{1}(:));
    partsI2(i,ti) = sum(SIPgrids{2}(:));
    
    totI1(ti) = totI1(ti) + partsI1(i,ti);
    totI2(ti) = totI2(ti) + partsI2(i,ti);

    A(ti) = A(ti) + sum(sqrtdetggrid(partitionMask{i}));
    
    ch1 = double(SIPgrids{1}).*sqrtdetggrid;
    ch2 = double(SIPgrids{2}).*sqrtdetggrid;

    partsI1prop(i, ti) = sum(ch1(~isnan(ch1)));
    partsI2prop(i, ti) = sum(ch2(~isnan(ch2)));
    
   	totI1prop(ti) = totI1prop(ti) + partsI1prop(i, ti);
    totI2prop(ti) = totI2prop(ti) + partsI2prop(i, ti);
end

% totI1
% totI2
A(ti)
totI1prop(ti)
totI2prop(ti)

%% CREATE 3D MASK

tic
SOImask = SOImask3D(mesh, onionOpts, xp.stack.imageSize([2 1 3]));
toc

%%
totI13D(ti) = sum(xp.stack.image.apply{1}(SOImask));
totI23D(ti) = sum(xp.stack.image.apply{2}(SOImask));

totI13D(ti)./totI1(ti)
totI13D(ti)./totI1prop(ti)
totI23D(ti)./totI2(ti)
totI23D(ti)./totI2prop(ti)

% %% AREA AND VOLUME DIRECT
% 
% A = 0;
% V = 0;
% 
% for i = 1:3
%     
% 	X = xp.SOI.embedding(xp.currentTime).patches{i}.apply;
%     
%     geom = GaussGeometry(X{1},X{2},X{3},0);
%     N = geom.N;
%     dA = geom.dA;
%     
%     dV = dA.*(X{1}.*N{1} + X{2}.*N{2} + X{3}.*N{3});
%     
%     V = V + sum(dV(partitionMask{i}))/3;
%     A = A + sum(dA(partitionMask{i}));
% end
% V measured in meshlab = 9695943

%%
stackAll = cat(4, xp.stack.image.apply{1}, xp.stack.image.apply{2}, SOImask);

zpos = 190;
imshow(cat(3,mat2gray(stackAll(:,:,zpos,1)),mat2gray(stackAll(:,:,zpos,2)),mat2gray(stackAll(:,:,zpos,3))),[])
cmap = [1 0 0; 0 1 0; 0 1 1];
ztol = 4;
for pi = 1:3
for li = 41%[1 21 41]
    
    X = embeddingStack{pi}(:,:,li,1);
    Y = embeddingStack{pi}(:,:,li,2);
    Z =  embeddingStack{pi}(:,:,li,3);
    
    ptsIdx = Z(:) < zpos + ztol & Z(:) > zpos - ztol;
    
    hold on
    scatter(X(ptsIdx),Y(ptsIdx), '.', 'MarkerEdgeColor', cmap(pi,:))
    hold off
end
end

% so patchwise normally evolve seems to work pretty well but the sampling
% density is just small in some regions
% maybe for conformal 4 we can get some agreement

%%
%----------------------------------------------------------------------
% loop for video
%----------------------------------------------------------------------

imwriteOptions = {'tif'};
saveDir = fullfile(projectDir, 'conformalHeart');
saveDir = '/Users/idse/conformalHeartSIP2';

options = struct('dir',saveDir,'imwriteOptions',{imwriteOptions},...
                    'make8bit',false);

% seeds positions 
seedsX = mesh.v(seeds,:);
% rotSeeds = cat(1, xp.fitter.fitOptions.chartSeeds(:,2), xp.fitter.fitOptions.diskSeeds(:,2));
% rotX = mesh.v(rotSeeds,:);
seedsXall{seedsCounter} = seedsX;
seedsCounter = seedsCounter+1;

% point in parametrization to hold fixed next round
fixedPtU = {};
fixedPtX = {};
for k = 1:nSeeds
    subm = xp.fitter.fittedParam.submeshes{k};
    fixedPtUtmp = sqrt(mean(subm.u{1}.^2)/2);
    fixedPtIdx = pointMatch(fixedPtUtmp, subm.u{1});
    fixedPtX{k} = subm.v(fixedPtIdx,:);
    fixedPtU{k} = subm.u{1}(fixedPtIdx,:);
end

% for debugging: all fittedParam
fittedParam = {};
fittedParam{fileMeta.timePoints(1)} = xp.fitter.fittedParam;
%%
% THE LOOP, 25, 30, 46, 57, 110 is a problem
for t = fileMeta.timePoints(111:end) %:end
    
    tid = tic;
    % raw data loading
    xp.loadTime(t);
    xp.rescaleStackToUnitAspect();
    
    % Read surface mesh produced by meshlab
    outputMesh = sprintf(outputMeshFormat, t);
    mesh = read_ply_mod(outputMesh);
    mesh.v = mesh.v(:,meshPermutation); 
    mesh.vn = mesh.vn(:,meshPermutation);
    
    % seeds for chart centers
    seeds = pointMatch(seedsX, mesh.v);
    diskSeeds = seeds(nVorSeeds+1:end);
    VorSeeds = seeds(1:nVorSeeds);

    % initialize fitter with overlap + disk
    fitOptions.chartSeeds = VorSeeds;
    fitOptions.diskSeeds = diskSeeds;
    fitOptions.fixedPtU = fixedPtU;
    fitOptions.fixedPtX = fixedPtX;
    xp.setFitOptions(fitOptions);
    xp.fitSurface(mesh); 

    % for debugging
    fittedParam{t} = xp.fitter.fittedParam;
    
    % seeds positions for the next iteration (rotX may come from popSOI)
    seedsX = mesh.v(seeds,:);
    seedsXall{seedsCounter} = seedsX;
    seedsCounter = seedsCounter+1;
%     rotSeeds = cat(1, xp.fitter.fitOptions.chartSeeds(:,2), xp.fitter.fitOptions.diskSeeds(:,2));
%     rotX = mesh.v(rotSeeds,:);

    % for next time, propagate fixedPtX
    for k = 1:nSeeds
        V = xp.fitter.fittedParam.submeshes{k}.v;
        fixedPtIdx = pointMatch(fixedPtX{k}, V);
        fixedPtX{k} = V(fixedPtIdx,:);
    end
        
    % populate SOI
    xp.fitter.setDesiredChart('conformal', true);
    xp.fitter.setDesiredChart('exponential', false);
    xp.fitter.setFitDomain(xp.stack.image.domain);
    xp.fitter.populateSOI(xp.SOI, xp.currentTime);

    % Pullback the stack to the desired charts
    xp.SOI.pullbackStack(xp.stack, [], xp.currentTime, onionOpts);
    
    % INTENSTITY MEASUREMENTS
    
    disp('2D intensity measurement');
    cover = 1:nVorSeeds; 
    margin = round(fitOptions.transitionWidth/2);
    partitionMask = xp.SOI.partition(t, cover, margin);
    
    for i = 1:nVorSeeds

        X = xp.SOI.embedding(t).patches{i}.apply;
        geom = GaussGeometry(X{1},X{2},X{3},0);
        sqrtdetggrid = geom.dA;

        SIP = this.getField('data_SIP');
        for j = 1:2
            SIPgrids{j} = partitionMask{i}.*SIP(t).patches{i}.apply{j};
        end
        
        partsI1(i,t) = sum(SIPgrids{1}(:));
        partsI2(i,t) = sum(SIPgrids{2}(:));

        totI1(t) = totI1(t) + partsI1(i,t);
        totI2(t) = totI2(t) + partsI2(i,t);

        A(t) = A(t) + sum(sqrtdetggrid(partitionMask{i}));

        ch1 = double(SIPgrids{1}).*sqrtdetggrid;
        ch2 = double(SIPgrids{2}).*sqrtdetggrid;

        partsI1prop(i, t) = sum(ch1(~isnan(ch1)));
        partsI2prop(i, t) = sum(ch2(~isnan(ch2)));

        totI1prop(t) = totI1prop(t) + partsI1prop(i, t);
        totI2prop(t) = totI2prop(t) + partsI2prop(i, t);
    end

    disp('3D intensity masking');
    tidmask = tic;
    SOImask = SOImask3D(mesh, onionOpts, xp.stack.imageSize([2 1 3]));
    toc(tidmask)
    totI13D(t) = sum(xp.stack.image.apply{1}(SOImask));
    totI23D(t) = sum(xp.stack.image.apply{2}(SOImask));
    
    save('intensityMeasurements', 'totI1', 'totI2', 'totI1prop', 'totI2prop', 'totI13D', 'totI23D', 'partsI1', 'partsI2', 'partsI1prop', 'partsI2prop');
    
    % Save
    if rem(t,20) == 0
        xp.SOI.save(options)
    end
    
    % save fittedParam
    fpfname = sprintf('fittedParam%d',t);
    fittedParamT = xp.fitter.fittedParam;
    save(fpfname, 'fittedParamT');

    disp('HERE COMES THE TOC TOC TOC TOC TOC TOC TOC TOC TOC TOC TOC TOC TOC TOC TOC TOC TOC TOC TOC TOC TOC TOC TOC');
    toc(tid)
end

%%
for t = 1:2
    t
    [totI13D(t)./totI1(t), totI23D(t)./totI2(t)]
    [totI13D(t)./totI1prop(t), totI23D(t)./totI2prop(t)]
end

%% INTENSITY GRAPH OVER TIME

% RESOLUTION : 3ms

cd('/Users/idse/Dropbox/flows/flows_shared/data/zebrafish/huisken/video9/')
load('intensityMeasurements')

cmap = [[0.8 0 0.4]; 
        0.8 1 0;
        0 1 0.4;
        1 0.5 0];
    
x = 3*(1:129);

scrsz = get(groot,'ScreenSize');
figure('Position',[1 scrsz(4)/2 scrsz(3)/2 scrsz(4)/2])

hold on

for i = 3:-1:1
    
    y =  sum(partsI1prop(1:i,:),1);
    %y =  sum(partsI1(1:i,:),1);
    y(y>0) = y(y>0)./totI13D(y>0);
    plot(x(y>0),y(y>0), 'k', 'LineWidth', 2)
    h = fill([x(y>0) x(end) x(1)], [y(y>0) 0 0], cmap(i,:),'EdgeColor','k', 'LineWidth', 2); 
    
%     y =  sum(partsI1(1:i,:),1);
%     y(y>0) = y(y>0)./totI13D(y>0);
%     plot(x(y>0),y(y>0), '-k', 'LineWidth', 4)
%     plot(x(y>0),y(y>0), '-','Color',cmap(i,:), 'LineWidth', 2)
end

axis([x(1) x(end) 0 1.5]);
fsize = 26;
% xlabel('time','fontweight','bold','fontsize',fsize);
% ylabel('relative intensity', 'fontweight','bold','fontsize',fsize);
set(gca,'color','none')
set(gca, 'LineWidth', 2);
set(gca,'FontSize', fsize)
set(gca,'FontWeight', 'bold')
pbaspect([2 1 1]);
hold off

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
saveDir = fullfile(projectDir, 'conformalHeartBLA');
saveDir = '/Users/idse/test';
options = struct('dir',saveDir,'imwriteOptions',{imwriteOptions},...
                    'make8bit',false);
xp.SOI.save(options)

%%

saveDir = '/Users/idse/Dropbox/flows/flows_shared/data/zebrafish/huisken/video9/conformalHeartMar23';
tic
fullSOI = surfaceAnalysis.SurfaceOfInterest(saveDir);
toc

%%
onionOpts = struct('nLayers', 41, 'layerDistance', 1, 'sigma', 1,...
                    'makeIP', 'SIP', 'IPonly', false);
                
fullSOI.makeIP(1, onionOpts)

% MAKE IT READ SOI.nLayers property?
