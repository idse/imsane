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

%%
profile on

dN = 20;

% addpath(genpath('/Users/idse/Downloads/Mesh_voxelisation'))

normal = mesh.vn./repmat(sqrt(sum(mesh.vn.^2,2)),[1 3]);
Vouter = mesh.v + dN*normal;
Vinner = mesh.v - dN*normal;

meshFV = struct();
meshFV.faces = mesh.f;

meshFV.vertices = Vouter;
maxes = round(max(meshFV.vertices));
mins = round(min(meshFV.vertices));
sizes = maxes - mins;
SOImaskOuter = VOXELISE(sizes(1),sizes(2),sizes(3),meshFV);

meshFV.vertices = Vinner;
maxes = round(max(meshFV.vertices));
mins = round(min(meshFV.vertices));
sizes = maxes - mins;
SOImaskInner = VOXELISE(sizes(1),sizes(2),sizes(3),meshFV);

offset = round(min(Vinner)) - round(min(Vouter));
bla = SOImaskOuter(offset(1):offset(1)+sizes(1)-1, offset(2):offset(2)+sizes(2)-1,...
                offset(3):offset(3)+sizes(3)-1);
SOImaskActual = SOImaskOuter;
SOImaskActual(offset(1):offset(1)+sizes(1)-1, offset(2):offset(2)+sizes(2)-1,...
                offset(3):offset(3)+sizes(3)-1) = bla - SOImaskInner;
profile viewer

%%

dN = 23;
%addpath(genpath('/Users/idse/Downloads/Mesh_voxelisation'))

normal = mesh.vn./repmat(sqrt(sum(mesh.vn.^2,2)),[1 3]);
Vouter = mesh.v + dN*normal;

meshFV = struct();
meshFV.faces = mesh.f;

meshFV.vertices = Vouter;
maxes = round(max(meshFV.vertices));
mins = round(min(meshFV.vertices));
sizes = maxes - mins;
SOImaskOut = VOXELISE(sizes(1),sizes(2),sizes(3),meshFV);

V = Vouter;

SOImask = {};

for i = 1:15
i
    V = V - 3*normal;

    meshFV.vertices = V;
    maxes = round(max(meshFV.vertices));
    mins = round(min(meshFV.vertices));
    sizes = maxes - mins;
    SOImaskIn = VOXELISE(sizes(1),sizes(2),sizes(3),meshFV);

    offset = round(min(V)) - round(min(Vouter));

    bla = SOImaskOut(offset(1):offset(1)+sizes(1)-1, offset(2):offset(2)+sizes(2)-1,...
                    offset(3):offset(3)+sizes(3)-1);
    SOImask{i} = SOImaskOut;
    SOImask{i}(offset(1):offset(1)+sizes(1)-1, offset(2):offset(2)+sizes(2)-1,...
                    offset(3):offset(3)+sizes(3)-1) = bla - SOImaskIn;
end

%%
greys = [0 1 2 3 4 3 2 1 0]/4;
greys = [2 1 1 1 1 1 2 2]/2;
L{1} = greys(1)*SOImask{1};
SOImaskTot  = L{1};
for i = 1:6
  
    SOImaskTot = SOImaskTot + greys(i+1)*(SOImask{2*i+1} - SOImask{2*i});
end
%%
i = 3;
markedLayer = SOImask{2*i+1} - SOImask{2*i};

%%
imshow(SOImaskTot(:,:,200))

%%

% make a mask made of "grid lines"
i = 3;
middleLayerMask = false(size(SOImask{1}));
for y = 50:50:350
    middleLayerMask(y:y+2,:,:) = SOImask{1}(y:y+2,:,:);
    %middleLayerMask(y:y+2,:,:) = SOImask{2*i+1}(y:y+2,:,:) - SOImask{2*i}(y:y+2,:,:);
end
z = 150;
middleLayerMask(:,:,z) = SOImask{1}(:,:,z);
%middleLayerMask(:,:,z) = SOImask{2*i+1}(:,:,z) - SOImask{2*i}(:,:,z);

imshow(SOImaskTot(:,:,200) + middleLayerMask(:,:,200))

%%

stacks = xp.stack.image.apply;

fname = fullfile(dataDir, 'SOIOnionMarkedLayer.tif');
if exist(fname)
    delete(fname);
end

zpos = 209;
Ilims = {};
for ch = 1:2
    Ilims{ch} = double([1.05*min(stacks{ch}(:)) 0.8*max(stacks{ch}(:))]);
end

Ilims = {[300 3500], [300 2500]};

meshFV.vertices = Vouter;
maxes = round(max(meshFV.vertices));
mins = round(min(meshFV.vertices));
sizes = maxes - mins;

for zpos = 1:xp.stack.imageSize(3) %218
    
    R = mat2gray(stacks{1}(:,:,zpos), Ilims{1});
    G = mat2gray(stacks{2}(:,:,zpos), Ilims{2});

    solidMaskSl = false(size(R));
    layer = zeros(size(R));
    if zpos > mins(3) && zpos < maxes(3)
        solidMaskSl(mins(2):mins(2)+sizes(2)-1, mins(1):mins(1)+sizes(1)-1) = SOImask{end}(:,:,zpos-mins(3))';
        %maskSl(mins(2):mins(2)+sizes(2)-1, mins(1):mins(1)+sizes(1)-1) = SOImask{end}(:,:,zpos-mins(3))';
        layer(mins(2):mins(2)+sizes(2)-1, mins(1):mins(1)+sizes(1)-1) = middleLayerMask(:,:,zpos-mins(3))';
    end
    
    alpha = 1;
    
    % add the grid of the outer layer
    R = R + alpha*layer;
    B = G;
    G = G + alpha*layer;
    
    % replace part by the onion mask image
    if zpos > mins(3)+150
        layeredMaskSl = zeros(size(R));
        markedSl = zeros(size(R));
        if zpos > mins(3) && zpos < maxes(3)
            layeredMaskSl(mins(2):mins(2)+sizes(2)-1, mins(1):mins(1)+sizes(1)-1) = SOImaskTot(:,:,zpos-mins(3))';
            markedSl(mins(2):mins(2)+sizes(2)-1, mins(1):mins(1)+sizes(1)-1) = markedLayer(:,:,zpos-mins(3))';
        end
        xpos = mins(1)+248;
        R(:,xpos:end) = alpha*layeredMaskSl(:,xpos:end) + markedSl(:,xpos:end);
        G(:,xpos:end) = alpha*layeredMaskSl(:,xpos:end) + markedSl(:,xpos:end);
        B(:,xpos:end) = markedSl(:,xpos:end);
        
    end
    
    % cut out
    if zpos > mins(3) + 152
        R(:,xpos+2:end) = 0;
        G(:,xpos+2:end) = 0;
        B(:,xpos+2:end) = 0;
    end
    colim = cat(3, R.*solidMaskSl, G.*solidMaskSl, B.*solidMaskSl);
    %colim = cat(3, R, G, B);
    
    imwrite(colim,fname,'tiff','Compression','none','WriteMode','append');
end

    %colim(1:280,1:200,1:3) = 0;
    figure,
    imshow(colim)

% ztol = 1;
% vinslice = mesh.v(:,3) < zpos + ztol & mesh.v(:,3) > zpos - ztol;
% hold on
% scatter(mesh.v(vinslice,1),mesh.v(vinslice,2),mesh.v(vinslice,3),'.g');
% scatter(Vinner(vinslice,1),Vinner(vinslice,2),Vinner(vinslice,3),'.g');
% scatter(Vouter(vinslice,1),Vinner(vinslice,2),Vinner(vinslice,3),'.b');
% hold off

%% set seeds points

seedsXinit = [127 338 106; 311 290 110; 332 326 252; 238 436 184; 127 205 273];
seedsXinit = [127 338 106; 311 290 110; 332 326 252; 307 431 185]%; 127 205 273];
seeds = pointMatch(seedsXinit, mesh.v);

nSeeds = numel(seeds);
nVorSeeds = 3;
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
% Voronoi cells with zero overlap for visualization

fitOptions = struct('chartSeeds', VorSeeds, 'transitionWidth', 0);
xp.setFitOptions(fitOptions);
xp.fitSurface(mesh);
noOverlap = xp.fitter.fittedParam;

%%
% full atlas with overlap for everything else

fitOptions = struct('chartSeeds', VorSeeds, 'transitionWidth', 50,...
                    'diskSeeds', diskSeeds, 'diskRadius', 150, 'makeTMaps', false);
xp.setFitOptions(fitOptions);
xp.fitSurface(mesh);

% %%
% fname = sprintf('fittedParam%d',xp.currentTime);
% fittedParam = xp.fitter.fittedParam;
% save(fname, 'fittedParam');

%% visualize 

% whole mesh
xp.fitter.inspectMesh(1:3);


%%
% whole heartedly 

this = xp.fitter;

% color different regions of the surface
v = this.fittedParam.mesh.v;
C = 0*v(:,1);
for i = 1:nVorSeeds
    C(noOverlap.submVidx{i} & ~this.fittedParam.submVidx{4}) = i;
end
for j = nVorSeeds+1:nSeeds
    C(this.fittedParam.submVidx{j}) = j;
end

trisurf(this.fittedParam.mesh.f, v(:,1), v(:,2), v(:,3), C, 'CDataMapping','direct')
shading interp 
cmap = hsv(numel(this.fittedParam.submeshes));
cmap = copper(numel(this.fittedParam.submeshes));
colormap(cmap);
axis equal

% add the chart boundaries 
hold on
for i = 1:numel(this.fittedParam.submeshes)
    v = this.fittedParam.submeshes{i}.v;
    bidx = this.fittedParam.submeshes{i}.b{1};
    plot3(v(bidx,1),v(bidx,2),v(bidx,3), '--', 'Color', cmap(i,:), 'LineWidth', 2);
end
hold off
camlight;
axis off;
viewdir = [1 0 1];
view(viewdir);

% seeds and normals
hold on
for i = 1:nSeeds
    v = this.fittedParam.mesh.v;
    vn = this.fittedParam.mesh.vn;
    scatter3(v(seeds(i),1), v(seeds(i),2),v(seeds(i),3),'filled', '.r')
    quiver3(v(seeds(i),1), v(seeds(i),2),v(seeds(i),3),...
            vn(seeds(i),1), vn(seeds(i),2),vn(seeds(i),3),100)
end
hold off

%%
%---------------------------------------------------------------------
% broken heart
%---------------------------------------------------------------------

cmap = hsv(numel(this.fittedParam.submeshes));
%cmap = copper(numel(this.fittedParam.submeshes));
separation = 30;

% voronoi only colors
Cvor = 0*v(:,1);
for i = 1:nVorSeeds
    Cvor(noOverlap.submVidx{i}) = i;
end

clf;
hold on;

% the Voronoi cells
for i = 1:nVorSeeds

    fittedParam = xp.fitter.fittedParam; % noOverlap
    submVidx = fittedParam.submVidx{i};
    v = fittedParam.submeshes{i}.v;
    f = fittedParam.submeshes{i}.f;

    d = mean(fittedParam.submeshes{i}.vn(:,[1 2 3]));
    d = d./norm(d);
    d = d*separation;

    for j = 1:3
        v(:,j) = v(:,j) +  d(j);
    end    
    trisurf(f, v(:,1), v(:,2),v(:,3), 0*v(:,1) + i, 'CDataMapping','direct')
end

% the stitches
for i = 1:nVorSeeds

    v = xp.fitter.fittedParam.submeshes{i}.v;
    f = xp.fitter.fittedParam.submeshes{i}.f;
    
    d = mean(xp.fitter.fittedParam.submeshes{i}.vn);
    d = d./norm(d);
    d = d*separation;

    for j = 1:3
        v(:,j) = v(:,j) +  d(j);
    end    

    bidx = xp.fitter.fittedParam.submeshes{i}.b{1};
    plot3(v(bidx,1),v(bidx,2),v(bidx,3), '--', 'Color', cmap(i,:), 'LineWidth', 2);
    
    % disk stitches on other regions
    for j = nVorSeeds + 1:nSeeds;
        
        if this.fittedParam.intersects(i,j)
            
            vdisk = xp.fitter.fittedParam.submeshes{j}.v;
            for k = 1:3
                vdisk(:,k) = vdisk(:,k) +  d(k);
            end  
            % from full mesh to submesh index
            submVidxj = this.fittedParam.submVidx{j};
            old2newj = zeros(size(submVidxj));
            old2newj(submVidxj) = 1:sum(submVidxj);

            % indices in full mesh of intersection
            isectIdx = this.fittedParam.intersectIdx{i,j};

            % indeces in j of points that are also contained in i
            isectIdxj = false([sum(submVidxj) 1]);
            isectIdxj(old2newj(isectIdx)) = true;

            % create a mask on the boundary curve of the disk
            bidxj = xp.fitter.fittedParam.submeshes{j}.b{1};
            bidxMask = true(size(bidxj));
            for l = find(~isectIdxj)'
                bidxMask(bidxj == l) = false;
            end

            % we need the connected components to plot as separate curves
            CC = bwconncomp(bidxMask);
            for cci = 1:CC.NumObjects
                bidxjp = bidxj(CC.PixelIdxList{cci});
                plot3(vdisk(bidxjp,1),vdisk(bidxjp,2),vdisk(bidxjp,3), '--', 'Color', cmap(j,:), 'LineWidth', 2);
            end
        end
    end
end

% add the disks
for j = nVorSeeds + 1:nSeeds;
    VorSeeds = xp.fitter.fitOptions.diskSeeds(j - nVorSeeds);
    submVidx = xp.fitter.fittedParam.submVidx{j};
    v = xp.fitter.fittedParam.submeshes{j}.v;
    f = xp.fitter.fittedParam.submeshes{j}.f;
    d = mean(xp.fitter.fittedParam.submeshes{j}.vn);
    d = d./norm(d);
    d = 4*d*separation;
    for k = 1:3
        v(:,k) = v(:,k) +  d(k);
    end    
    trisurf(f, v(:,1), v(:,2),v(:,3), C(submVidx,:), 'CDataMapping','direct')
end
axis equal;
shading interp
colormap(cmap)
camlight
view(viewdir);
axis off;
hold off;

%%
% submesh
xp.fitter.inspectMesh();
view([0 0 1]);

%% Inspect the fit in a cross section
%
% inspect fit and point cloud over a cross section inthe  data. Dimensions 
% are x,y or z and the value has to be within the corresponding axis range.
%
clf
zval = 200;
inspectOptions = struct('dimension','x','value',zval,'pointCloud','b', 'noalign', true);
xp.fitter.inspectQuality(inspectOptions, xp.detector, xp.stack);


%% normally evolve

shift = 2;
xp.normallyEvolve(shift);

%%
xp.fitter.setDesiredChart('conformal', true);
xp.fitter.setDesiredChart('exponential', false);
xp.fitter.setFitDomain(xp.stack.image.domain);
xp.generateSOI();

%% Pullback the stack to the desired charts
%
% Pass the region of interest and the current time to pull back the stack
% in the desired charts. This generates the data fields containing the
% pullback.
%
% 21, 2, 20
% 21, 2, 1
onionOpts = struct('nLayers', 41, 'layerDistance', 1, 'sigma', 1,...
                    'makeIP', 'both', 'IPonly', false);
xp.SOI.pullbackStack(xp.stack, [], xp.currentTime, onionOpts);

%%
%----------------------------------------------------------------------
% exploding heart with texture map
%----------------------------------------------------------------------

viewdir = ([1 1 1]);
az = 180; el = -90;
az = 115; el = -30;
az = 130; el = -60;
az = 110; el = -60;
az = 126; el = -68;
fittedParam = xp.fitter.fittedParam;
%cmap = hsv(numel(fittedParam.submeshes));
cmap = [[0.8 0 0.4]; 
        0.8 1 0;
        0 1 0.4;
        1 0.5 0];

tidx = xp.tIdx(xp.currentTime);
emb = xp.SOI.embedding(tidx);

separation = 140*[1 1 1 1.7 2];
nSeeds = 4;
extra = zeros([nSeeds 3]);
% extra(1,:) = 0*[0  3 1];
% extra(2,:) = 0*[0 -1 0];
% extra(3,:) = 0*[1 0 0];
% extra(4,:) = 0*[0 1 -4];

figure
clf
set(gcf,'color','w');

hold on;


%MIP = xp.SOI.getField('data');
%adjust = [0.3 0.6]; %for data

MIP = xp.SOI.getField('data_SIP');
adjust = [0.2 0.8]; %for SIP

%adjust = [0.5 0.8]; %for MIP


% texture 
C = {};
Cc = {};
Cmin = {};
Cmax = {};

for i = 1:nSeeds
    for ci = 1:numel(xp.expMeta.channelsUsed)
        im = MIP(tidx).patches{i}.apply{ci};
        if i == 1
            Cmin{ci} = double(min(im(:)));
            Cmax{ci} = adjust(ci)*double(max(im(:)));
        end
        Cc{i,ci} = mat2gray(im, [Cmin{ci} Cmax{ci}]);
    end
    if size(Cc,2) > 1
        C{i} = cat(3,Cc{i,1},Cc{i,2},Cc{i,2});
    else
        C{i} = Cc{i};
    end
end

% the pieces
ambstr = [0.1 0.4 0 0.1];
diffstr = [2 2 1 3];
for i = 1:nSeeds
    X = emb.patches{i}.apply;
    d = mean(xp.fitter.fittedParam.submeshes{i}.vn);
    d = d./norm(d);
    d = separation(i)*d;
    d = d+extra(i,:);
    
    hpatch = surf(X{1} + d(1),X{2} + d(2),X{3} + d(3), C{i}, 'FaceColor','texturemap');
    set(hpatch,'AmbientStrength',ambstr(i),...
        'SpecularStrength',1,...
      'DiffuseStrength',diffstr(i));
    %camlight('headlight');   
end

axis equal;
shading flat
colormap gray
%hcam = camlight('headlight');
light('style', 'infinite')
xl = 300;%300
yl = 400;%400
zl = -100;%-100
light('Position',[xl yl zl])
%scatter3(xl,yl,zl,'g');
%view(viewdir);
view(az,el);
%view([0 1 0]);
axis off;

% the stitches
for i = 1:nSeeds

    v = fittedParam.submeshes{i}.v;
    f = fittedParam.submeshes{i}.f;

    vdisk = fittedParam.submeshes{4}.v;
    ndisk = fittedParam.submeshes{4}.vn;
    
    d = mean(fittedParam.submeshes{i}.vn);
    d = d./norm(d);
    
    d = d*(separation(i));

    d = d+extra(i,:);
    
    for j = 1:3
        v(:,j) = v(:,j) +  d(j);
        if i == 3
            vdisk(:,j) = vdisk(:,j) +  d(j) + 2*d(j)./norm(d);
        else 
            vdisk(:,j) = vdisk(:,j) +  d(j);
        end
    end    

    bdryss = 3;
    lw = 5;
    st = '-';
    bidx = [fittedParam.submeshes{i}.b{1}(1:bdryss:end) fittedParam.submeshes{i}.b{1}(end)];
    bidxc = [bidx bidx(1)];
    
    % disk stitches on other regions
    j = 4;
    if fittedParam.intersects(i,j)

        % from full mesh to submesh index
        submVidxj = fittedParam.submVidx{j};
        old2newj = zeros(size(submVidxj));
        old2newj(submVidxj) = 1:sum(submVidxj);

        % indices in full mesh of intersection
        isectIdx = fittedParam.intersectIdx{i,j};

        % indeces in j of points that are also contained in i
        isectIdxj = false([sum(submVidxj) 1]);
        isectIdxj(old2newj(isectIdx)) = true;
        
        % create a mask on the boundary curve of the disk
        bidxj = fittedParam.submeshes{4}.b{1}(1:end);
        bidxj = [bidxj bidxj(1)];
        bidxMask = true(size(bidxj));
        for l = find(~isectIdxj)'
            bidxMask(bidxj == l) = false;
        end
        
        % we need the connected components to plot as separate curves
        CC = bwconncomp(bidxMask);
        for cci = 1:CC.NumObjects
            bidxjp = bidxj(CC.PixelIdxList{cci});
            vdiskb = vdisk(bidxjp,:);
            ndiskb = ndisk(bidxjp,:);
            vdiskb = vdiskb + 2*ndiskb;
            plot3(vdiskb([1:bdryss:end end],1),vdiskb([1:bdryss:end end],2),vdiskb([1:bdryss:end end],3), st, 'Color', cmap(4,:), 'LineWidth', lw);
        end
    end
    
    plot3(v(bidxc,1),v(bidxc,2),v(bidxc,3), st, 'Color', cmap(i,:), 'LineWidth', lw);
end
% 
% for i = 1:nSeeds
%     v = fittedParam.mesh.v;
%     vn = fittedParam.mesh.vn;
%     scatter3(v(seeds(i,1),1), v(seeds(i,1),2),v(seeds(i,1),3),'filled', '.r')
% %     quiver3(v(seeds(i,1),1), v(seeds(i,1),2),v(seeds(i,1),3),...
% %             vn(seeds(i,1),1), vn(seeds(i,1),2),vn(seeds(i,1),3),100)
% end

hold off;
%set(gca,'XLim',[0 550])
%set(gca,'ZLim',[0 500])

%%
%----------------------------------------------------------------------
% just the disk 3D
%----------------------------------------------------------------------

viewdir = ([1 1 1]);
az = 180; el = -90;
az = 115; el = -30;
az = 130; el = -60;
az = 110; el = -60;
az = -180; el = -46;
fittedParam = xp.fitter.fittedParam;
%cmap = hsv(numel(fittedParam.submeshes));
cmap = [[0.8 0 0.4]; 
        0.8 1 0;
        0 1 0.4;
        1 0.5 0];

tidx = xp.tIdx(xp.currentTime);
emb = xp.SOI.embedding(tidx);

separation = 140*[1 1 1 1.7 2];
nSeeds = 4;
extra = zeros([nSeeds 3]);
% extra(1,:) = 0*[0  3 1];
% extra(2,:) = 0*[0 -1 0];
% extra(3,:) = 0*[1 0 0];
% extra(4,:) = 0*[0 1 -4];

figure
clf
set(gcf,'color','w');

hold on;

MIP = xp.SOI.getField('data');

% texture 
C = {};
Cc = {};
Cmin = {};
Cmax = {};
%adjust = [0.2 0.8]; %adjust = [0.5 0.8]; for MIP
for i = 1:nSeeds
    for ci = 1:numel(xp.expMeta.channelsUsed)
        im = MIP(tidx).patches{i}.apply{ci};
        if i == 1
            Cmin{ci} = double(min(im(:)));
            Cmax{ci} = adjust(ci)*double(max(im(:)));
        end
        Cc{i,ci} = mat2gray(im, [Cmin{ci} Cmax{ci}]);
    end
    if size(Cc,2) > 1
        C{i} = cat(3,Cc{i,1},Cc{i,2},Cc{i,2});
    else
        C{i} = Cc{i};
    end
end

% the pieces
ambstr = [0.1 0.4 0 0.4];
diffstr = [2 2 1 1];
for i = 4%1:nSeeds
    X = emb.patches{i}.apply;
    d = mean(xp.fitter.fittedParam.submeshes{i}.vn);
    d = d./norm(d);
    d = separation(i)*d;
    d = d+extra(i,:);
    
    hpatch = surf(X{1} + d(1),X{2} + d(2),X{3} + d(3), C{i}, 'FaceColor','texturemap');
    set(hpatch,'AmbientStrength',ambstr(i),...
        'SpecularStrength',1,...
      'DiffuseStrength',diffstr(i));
    %camlight('headlight');   
end

axis equal;
shading flat
colormap gray
%hcam = camlight('headlight');
light('style', 'infinite')
xl = 300;%300
yl = 400;%400
zl = -100;%-100
light('Position',[xl yl zl])
%scatter3(xl,yl,zl,'g');
%view(viewdir);
view(az,el);
%view([0 1 0]);
axis off;

% the stitches
for i = 4%1:nSeeds

    v = fittedParam.submeshes{i}.v;
    f = fittedParam.submeshes{i}.f;

    vdisk = fittedParam.submeshes{4}.v;
    ndisk = fittedParam.submeshes{4}.vn;
    
    d = mean(fittedParam.submeshes{i}.vn);
    d = d./norm(d);
    
    d = d*(separation(i));

    d = d+extra(i,:);
    
    for j = 1:3
        v(:,j) = v(:,j) +  d(j);
        if i == 3
            vdisk(:,j) = vdisk(:,j) +  d(j) + 2*d(j)./norm(d);
        else 
            vdisk(:,j) = vdisk(:,j) +  d(j);
        end
    end    

    bdryss = 3;
    lw = 5;
    st = '-';
    bidx = [fittedParam.submeshes{i}.b{1}(1:bdryss:end) fittedParam.submeshes{i}.b{1}(end)];
    bidxc = [bidx bidx(1)];
    
    % disk stitches on other regions
    j = 4;
    if fittedParam.intersects(i,j)

        % from full mesh to submesh index
        submVidxj = fittedParam.submVidx{j};
        old2newj = zeros(size(submVidxj));
        old2newj(submVidxj) = 1:sum(submVidxj);

        % indices in full mesh of intersection
        isectIdx = fittedParam.intersectIdx{i,j};

        % indeces in j of points that are also contained in i
        isectIdxj = false([sum(submVidxj) 1]);
        isectIdxj(old2newj(isectIdx)) = true;
        
        % create a mask on the boundary curve of the disk
        bidxj = fittedParam.submeshes{4}.b{1}(1:end);
        bidxj = [bidxj bidxj(1)];
        bidxMask = true(size(bidxj));
        for l = find(~isectIdxj)'
            bidxMask(bidxj == l) = false;
        end
        
        % we need the connected components to plot as separate curves
        CC = bwconncomp(bidxMask);
        for cci = 1:CC.NumObjects
            bidxjp = bidxj(CC.PixelIdxList{cci});
            vdiskb = vdisk(bidxjp,:);
            ndiskb = ndisk(bidxjp,:);
            vdiskb = vdiskb + 2*ndiskb;
            plot3(vdiskb([1:bdryss:end end],1),vdiskb([1:bdryss:end end],2),vdiskb([1:bdryss:end end],3), st, 'Color', cmap(4,:), 'LineWidth', lw);
        end
    end
    
    plot3(v(bidxc,1),v(bidxc,2),v(bidxc,3), st, 'Color', cmap(i,:), 'LineWidth', lw);
end
% 
% for i = 1:nSeeds
%     v = fittedParam.mesh.v;
%     vn = fittedParam.mesh.vn;
%     scatter3(v(seeds(i,1),1), v(seeds(i,1),2),v(seeds(i,1),3),'filled', '.r')
% %     quiver3(v(seeds(i,1),1), v(seeds(i,1),2),v(seeds(i,1),3),...
% %             vn(seeds(i,1),1), vn(seeds(i,1),2),vn(seeds(i,1),3),100)
% end

hold off;
%set(gca,'XLim',[0 550])
%set(gca,'ZLim',[0 500])

%%
%----------------------------------------------------------------------
% atlas 
%----------------------------------------------------------------------

tidx = xp.tIdx(xp.currentTime);
data = xp.SOI.getField('data_SIP');
data = data(tidx);
bdryss = 5;


% BLAAA
MIP = xp.SOI.getField('data_SIP');

maxadjust = [0.4 0.9]; %for SIP
gamma = [1 1];

% texture 
C = {};
Cc = {};
Cmin = {};
Cmax = {};

for i = 1:nSeeds
    for ci = 1:numel(xp.expMeta.channelsUsed)
        im = MIP(tidx).patches{i}.apply{ci};
        if i == 1
            Cmin{ci} = double(min(im(:)));
            Cmax{ci} = adjust(ci)*double(max(im(:)));
        end
        im = mat2gray(im);
        Cc{i,ci} = imadjust(im, [0 maxadjust(ci)], [0 1], gamma(ci));
    end
    if size(Cc,2) > 1
        C{i} = cat(3,Cc{i,1},Cc{i,2},Cc{i,2});
    else
        C{i} = Cc{i};
    end
end

for i = 4%1:nSeeds;

%    addpath('/Users/idse/Dropbox/Ken/lattice_analysis/lattice_analysis_Aug28_2014/functions/')
%    subplot_tight(2,3,i);
    figure,

    type = 'conformal';
    patchName     = [type '_' num2str(i) '_index'];
    transformName = [type '_' num2str(i)];
    
    im = C{i};
    im(im==0) = 1;
    
    % could be commented out for old style boundaries
    im1 = data(tidx).patches{i}.apply{1};
    mask = imopen(im1 > 0, strel('disk', 20));
    mask = imerode(mask, strel('disk', 2));
    for c = 1:3
        im(:,:,c) = mask.*im(:,:,c);
    end
    im(im==0) = 1;
    % up to here
    
    im = permute(im, [2 1 3]);
    imshow(im,[],'InitialMagnification',66)

    % the stitches
    tidx = xp.tIdx(xp.currentTime);
    u = xp.fitter.fittedParam.submeshes{i}.u{1};
    f = xp.fitter.fittedParam.submeshes{i}.f;
    
    bidx = xp.fitter.fittedParam.submeshes{i}.b{1}(1:bdryss:end);
    
    curChart = xp.SOI.atlas(tidx).getChart(transformName);
    ub = u(bidx,1);
    vb = u(bidx,2);
    umin = curChart.image.boundary{1}(1);
    vmin = curChart.image.boundary{2}(1);
    ubp = (ub-umin)./curChart.image.stepSize(1);
    vbp = (vb-vmin)./curChart.image.stepSize(2);
    
    hold on
    
    % disk stitches on other regions
    j = 4;
    
    if fittedParam.intersects(i,j)

        % from full mesh to submesh index
        submVidxj = fittedParam.submVidx{j};
        old2newj = zeros(size(submVidxj));
        old2newj(submVidxj) = 1:sum(submVidxj);

        submVidxi = fittedParam.submVidx{i};
        old2newi = zeros(size(submVidxi));
        old2newi(submVidxi) = 1:sum(submVidxi);

        % indices in full mesh of intersection
        isectIdx = fittedParam.intersectIdx{i,j};

        % indeces in j of points that are also contained in i
        isectIdxj = false([sum(submVidxj) 1]);
        isectIdxj(old2newj(isectIdx)) = true;
        
        % create a mask on the boundary curve of the disk
        bidxj = fittedParam.submeshes{4}.b{1}(1:end);
        bidxj = [bidxj bidxj(1)];
        bidxMask = true(size(bidxj));
        for l = find(~isectIdxj)'
            bidxMask(bidxj == l) = false;
        end
        
        % we need the connected components to plot as separate curves
        CC = bwconncomp(bidxMask);
        for cci = 1:CC.NumObjects
            
            bidxjp = bidxj(CC.PixelIdxList{cci});
            
            new2oldj = find(submVidxj)';
            bidxIni = old2newi(new2oldj(bidxjp));
            
            udiskb = (u(bidxIni,1) - umin)/curChart.image.stepSize(1);
            vdiskb = (u(bidxIni,2) - vmin)/curChart.image.stepSize(2);
            udiskb = udiskb([1:bdryss:end end]);
            vdiskb = vdiskb([1:bdryss:end end]);
            plot(vdiskb, udiskb, st, 'Color', cmap(4,:), 'LineWidth', lw);
        end
    end
    
    B = bwboundaries(mask);
    hold on
    xpos = 110;
    plot(B{1}(:,1),B{1}(:,2),'Color',cmap(i,:),'LineWidth', lw)
    startdash = min(find(mask(xpos,:)));
    enddash = max(find(mask(xpos,:)));
    x = startdash:enddash;
    if i == 4
        plot(x,xpos*ones([1 length(x)]),'--','Color',cmap(i,:),'LineWidth',lw);
    end
    hold off

%     plot([ubp' ubp(1)],[vbp' vbp(1)], st, 'Color', cmap(i,:), 'LineWidth', lw);
    
    hold off
end

%%
%----------------------------------------------------------------------
% loop for video
%----------------------------------------------------------------------

imwriteOptions = {'tif'};
%saveDir = fullfile(projectDir, 'conformalHeart');
saveDir = '/Users/idse/conformalHeartSIP';

options = struct('dir',saveDir,'imwriteOptions',{imwriteOptions},...
                    'make8bit',false);

% seeds positions 
seedsX = mesh.v(seeds,:);
% rotSeeds = cat(1, xp.fitter.fitOptions.chartSeeds(:,2), xp.fitter.fitOptions.diskSeeds(:,2));
% rotX = mesh.v(rotSeeds,:);
seedsXall{seedsCounter} = seedsX;
seedsCounter = seedsCounter+1;

nSeeds = numel(seeds);

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

% 
% % starting position for seeds
% seedsX = [127 338 106; 311 290 110; 332 326 252;...
%             238 436 184; 127 205 273];
% nVorSeeds = 3;
% rotX = [];

for t = fileMeta.timePoints(121:end) %:end
    
    tic
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

%     % extra row of zeros for rotation seeds!!
%     diskSeeds = [diskSeeds, zeros(size(diskSeeds))];
%     VorSeeds = [VorSeeds, zeros(size(VorSeeds))];
%     
%     % seeds for chart rotation
%     if ~isempty(rotX)
%         rotSeeds = pointMatch(rotX, mesh.v);
%         diskSeeds(:,2) = rotSeeds(nVorSeeds+1:end);
%         VorSeeds(:,2) = rotSeeds(1:nVorSeeds);
%     end

% %________
%     fixedPtIdx = pointMatch(fixedPtX, subm{t}.v);
%     constrVidx = [subm{t}.seed fixedPtIdx];
%     constrU = [0 0; fixedPtU];
%     
%     % seed for the next one
%     if t < tmax
%         seed = pointMatch(seedX, mesh{t+1}.v);
%     end

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
    
%     % Save
%     if rem(t,20) == 0
%         xp.SOI.save(options)
%     end
    
    % save fittedParam
    fpfname = sprintf('fittedParam%d',t);
    fittedParamT = xp.fitter.fittedParam;
    save(fpfname, 'fittedParamT');

    disp('HERE COMES THE TOC TOC TOC TOC TOC TOC TOC TOC TOC TOC TOC TOC TOC TOC TOC TOC TOC TOC TOC TOC TOC TOC TOC');
    toc
end

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
saveDir = '/Users/idse/conformalHeartBla2';
options = struct('dir',saveDir,'imwriteOptions',{imwriteOptions},...
                    'make8bit',false);
xp.SOI.save(options)

%%
%----------------------------------------------------------------------
% conformal pringles
%----------------------------------------------------------------------

viewdir = ([0 0.1 1]);
az = 69; el = 43;
az = -200; el = -66;

cmap = [[0.8 0 0.4]; 
        0.8 1 0;
        0 1 0.4;
        1 0.5 0];

fittedParam = xp.fitter.fittedParam;

tidx = xp.tIdx(xp.currentTime);
emb = xp.SOI.embedding(tidx);

separation = 100*[1 1 1 1];

% texture 
C = {};
Cc = {};
Cmin = {};
Cmax = {};
adjust = [0.5 0.8];

nLayers = onionOpts.nLayers;%21;
halfLayers = (nLayers - 1)/2;

cstack = [];

for li = 1:nLayers
    
    idx = li - halfLayers - 1;
   
    % convert index to field name
    if idx < 0
        fieldName = ['data_layer_m' num2str(-idx)];
    elseif idx > 0
        fieldName = ['data_layer_p' num2str(idx)];
    else
        fieldName = 'data';
    end
    for ci = 1:numel(xp.expMeta.channelsUsed)
        
        MIP = xp.SOI.getField(fieldName);
        im = MIP(tidx).patches{4}.apply{ci};
        if li == 1
            Cmin{ci} = double(min(im(:)));
            Cmax{ci} = adjust(ci)*double(max(im(:)));
        end
        Cc{li,ci} = mat2gray(im, [Cmin{ci} Cmax{ci}]);
    end
    if size(Cc,2) > 1
        C{li} = cat(3,Cc{li,1},Cc{li,2},Cc{li,2});
    else
        C{li} = Cc{li};
    end
    
    cstack = cat(4, cstack, C{li});
end

cstack = permute(cstack, [1 2 4 3]);

figure
hold on;

% the stitches 2D
i =4;
tidx = xp.tIdx(xp.currentTime);
u = xp.fitter.fittedParam.submeshes{i}.u{1};
f = xp.fitter.fittedParam.submeshes{i}.f;
bidx = xp.fitter.fittedParam.submeshes{i}.b{1}(1:bdryss:end);

curChart = xp.SOI.atlas(tidx).getChart(transformName);
ub = u(bidx,1);
vb = u(bidx,2);
umin = 0;%curChart.image.boundary{1}(1);
vmin = 0;%curChart.image.boundary{2}(1);
ubp = (ub-umin)./curChart.image.stepSize(1);
vbp = (vb-vmin)./curChart.image.stepSize(2);

% the pieces
for li = 1:6:nLayers
% % 3D PRINGLE
%     X = emb.patches{4}.apply;
%     d = mean(xp.fitter.fittedParam.submeshes{4}.vn);
%     d = d./norm(d);
%     d = 40*(li-1)*d - 20*li*[1 0 0];
%     surf(X{1} + d(1),X{2} + d(2),X{3} + d(3), C{li}, 'FaceColor','texturemap');
% 2D STACK
    d = 15*li*[1 -1];
    surf(xp.SOI.atlas(1).charts{4}.apply{1} + d(1),...
            xp.SOI.atlas(1).charts{4}.apply{2} + d(2),...
            0*X{3} + 20*li,...
            C{li}, 'FaceColor','texturemap');
        
    plot3([ubp' ubp(1)] + d(1),[vbp' vbp(1)] + d(2), 0*[vbp' vbp(1)] + 20*li, st, 'Color', cmap(i,:), 'LineWidth', lw);
end

% X0 = mean(X{1}(~isnan(X{1})));
% Y0 = mean(X{2}(~isnan(X{2})));
% Z0 = mean(X{3}(~isnan(X{3})));
% 
%     li = 10;
%     d = mean(xp.fitter.fittedParam.submeshes{4}.vn);
%     d = d./norm(d);
%     d = 40*(li-1)*d - 20*li*[1 0 0];
%     %scatter(X0 + d(1) +200,Y0 + d(2),Z0 + d(3), '.b');
%     quiver3(X0 + d(1) +200,Y0 + d(2),Z0 + d(3),...
%              d(1), d(2), d(3), 1, 'LineWidth', 2, 'Color', 'red')
%      quiver3(X0 + d(1) +200,Y0 + d(2),Z0 + d(3),...
%          -d(1), -d(2), -d(3), 1, 'LineWidth', 2, 'Color', 'red')

axis equal;
shading flat
colormap gray
%camlight
%view(viewdir);
az = 161;
el = -28;
az = -144;
el = 24;
view(az,el);
%view([0 0 1]);
axis off;

hold off;

%% Xsection

xpos = 110;
R = flipud(imadjust(squeeze(cstack(:,xpos,:,1)), [0.1 0.8], [0 1], 0.9)');
G = flipud(imadjust(squeeze(cstack(:,xpos,:,2)), [0.1 1], [0 1], 1)');
mask = false(size(R));
bla = find(sum(R,1) > 0);
trimdist=4;
bla = bla(trimdist:end-trimdist+1);
% mask(:, bla) = true;
% R= mask.*R;
% G= mask.*G;
R = R(:,bla);
G = G(:,bla);
R(21,:)=1;
G(21,:)=1;
%figure, 
imshow(cat(3,R,G,G), [])

%% make time series figure

%data = xp.SOI.getField('data_MIP');
data = fullSOI.getField('data_MIP');
tindices = [1 3 5 7];
tindices = 1:5:30;
%tindices = [1 19 49 97];
tindices = [19 65 110];
gamma = [0.9 1 0.9];
highR = [0.4, 0.7, 0.45];
highC = [1, 0.9, 1];

%data = xp.SOI.getField('data_SIP');
%tindices = 1;

cmap = [[0.8 0 0.4]; 
        0.8 1 0;
        0 1 0.4;
        1 0.5 0];
lw = 3;

load('/Users/idse/Dropbox/flows/flows_shared/data/zebrafish/huisken/video9/ilt/TracksCoMove.mat')
load('/Users/idse/Dropbox/flows/flows_shared/data/zebrafish/huisken/video9/adjustCanvasOffset')
addpath('/Users/idse/Dropbox/Ken/lattice_analysis/lattice_analysis_Aug28_2014/functions/')
%tindices = 4;
pi = 4;

figure,

tidx = 49;
Ilim = {};
i =1;
Ilim{i} = double([min(data(tidx).patches{pi}.apply{i}(:)) max(data(tidx).patches{pi}.apply{i}(:))]);
tidx = 1;
i = 2;
Ilim{i} = double([min(data(tidx).patches{pi}.apply{i}(:)) 1.5*max(data(tidx).patches{pi}.apply{i}(:))]);

for i = 1:numel(tindices)

    tidx = tindices(i);
    subplot_tight(3,2*numel(tindices),[2*i-1, 2*i], [0 0]);
    %figure,

    im1 = data(tidx).patches{pi}.apply{1};
    mask = imopen(im1 > 0, strel('disk', 20));
    mask = imerode(mask, strel('disk', 2));
    %mask = true(size(im1));
    im2 = data(tidx).patches{pi}.apply{2};
    mask1 = im1 > 0;
    im1(im1==0) = max(im1(:));
    im2(im2==0) = max(im2(:));
    
    sigma = 5;
    im1 = mat2gray(im1);
    bg = imerode(im1, strel('disk', 14));
    im1 = im1 -  bg;
    
    im2 = mat2gray(im2);
    im2 = imfilter(im2,fspecial('Gauss', 6, 2));
    
    im1 = imadjust(mat2gray(im1), [0 highR(i)], [0 1], gamma(i));
    im2 = imadjust(mat2gray(im2), [0 highC(i)], [0 1], gamma(i));
    
    colim = cat(3, im1.*mask + ~mask, im2.*mask + ~mask, im2.*mask+ ~mask);

    imshow(colim);
    %title(['t = ' num2str(xp.fileMeta.timePoints(tidx))],...
    %                        'fontweight','bold','fontsize',16)
    B = bwboundaries(mask);
    hold on
    plot(B{1}(:,2),B{1}(:,1),'Color',cmap(4,:),'LineWidth', lw)
    hold off
    
    hold on 
    trackc = [];
    for trackidx = 1:numel(TracksCoMove)
        trackc = cat(1, trackc, TracksCoMove(trackidx).centres(tidx,:));
    end
    % order the points by hand, cause I'm a biologist
    bla = trackc(3,:);
    trackc(3,:) = trackc(5,:);
    trackc(5,:) = bla;
    trackc = cat(1,trackc,trackc(1,:));
    %plot(trackc(:,1), trackc(:,2),'-g');
    plot(trackc(:,1) - stDall(tidx,2), trackc(:,2)- stDall(tidx,1),'--',...
                        'Color', 0.8*[1 1 1], 'LineWidth', lw);
    hold off
    
%     
% 	u = fittedParam{fileMeta.timePoints(tidx)}.submeshes{pi}.u{1};
%     
%     cs = fittedParam{fileMeta.timePoints(tidx)}.submeshes{pi}.seed;
%     
%     hold on
%     umin = xp.SOI.atlas(tidx).charts{pi}.image.boundary{1}(1);
%     vmin = xp.SOI.atlas(tidx).charts{pi}.image.boundary{2}(1);
%     scatter(u(cs,1) - umin, u(cs,2) - vmin, 'r')
%     scatter(fixedPtU{pi}(1)-umin,fixedPtU{pi}(2)-vmin, 'g');
%     hold off
end

subplot_tight(3,2*numel(tindices),7:12, [0.01 0.05]);
%subplot(2, 3, [4 5 6]);
times = (1:tmax)*3;
yrange = 0:0.1:1.6;

hold on 
for i = 1:3
    plot(3*tindices(i)*ones(size(yrange)), yrange, 'LineWidth', 2, 'Color', cmap(4,:))
    plot(times, Aproper./max(Aproper), 'LineWidth', 4)
end
hold off

fsize = 14;

%# capture handle to current figure and axis
hFig = gcf;
hAx1 = gca;

%# create a second transparent axis, as a copy of the first
hAx2 = copyobj(hAx1,hFig);
delete( get(hAx2,'Children') )
set(hAx2, 'LineWidth', 2);
set(hAx2,'FontSize', fsize)
set(hAx2,'FontWeight', 'bold')
axis(hAx2, [1 times(end) yrange(1) yrange(end)]);
pbaspect(hAx2, [4 1 1]);
set(hAx2, 'Color','none', 'Box','off', ...
    'XGrid','off', 'YGrid','off')

%# show grid-lines of first axis, style them as desired,
%# but hide its tick marks and axis labels
set(hAx1, 'XColor',[0.9 0.9 0.9], 'YColor',0.8*[1 1 1], ...
    'XGrid','off', 'YGrid','on', 'MinorGridLineStyle','-', ...
    'XTickLabel',[], 'YTickLabel',[], 'Box', 'off');
xlabel(hAx1, ''), ylabel(hAx1, ''), title(hAx1, '')
set(hAx1, 'LineWidth', 2);
set(hAx1,'FontSize', fsize)
set(hAx1,'FontWeight', 'bold')
axis(hAx1, [1 times(end) yrange(1) yrange(end)]);
pbaspect(hAx1, [4 1 1]);
set(hAx1,'GridLineStyle','-')

%# link the two axes to share the same limits on pan/zoom
linkaxes([hAx1 hAx2], 'xy');


% INTENSITY GRAPH OVER TIME

% RESOLUTION : 3ms

cd('/Users/idse/Dropbox/flows/flows_shared/data/zebrafish/huisken/video9/')
load('intensityMeasurements')

cmap = [[0.8 0 0.4]; 
        0.8 1 0;
        0 1 0.4;
        1 0.5 0];
    
x = 3*(1:129);

%scrsz = get(groot,'ScreenSize');
%figure('Position',[1 scrsz(4)/2 scrsz(3)/2 scrsz(4)/2])

subplot_tight(3,2*numel(tindices),13:15, [0.01 0.05]);

hold on

for i = 3:-1:1
    
    y =  sum(partsI1prop(1:i,:),1);
    %y =  sum(partsI1(1:i,:),1);
    y(y>0) = y(y>0)./totI13D(y>0);
    plot(x(y>0),y(y>0), 'k', 'LineWidth', 2)
    h = fill([x(y>0) x(end) x(1)], [y(y>0) 0 0], cmap(i,:),'EdgeColor','k', 'LineWidth', 2); 
end

axis([x(1) x(end) 0 1.6]);
% xlabel('time','fontweight','bold','fontsize',fsize);
% ylabel('relative intensity', 'fontweight','bold','fontsize',fsize);
set(gca,'color','none')
set(gca, 'LineWidth', 2);
set(gca,'FontSize', fsize)
set(gca,'FontWeight', 'bold')
pbaspect([2 1 1]);
hold off
h = get(gca,'Children');
set(gca,'Children',[h(1) h(2) h(3) h(4) h(5) h(6)])

subplot_tight(3,2*numel(tindices),16:18, [0.01 0.05]);

hold on

for i = 3:-1:1
    
    %y =  sum(partsI1prop(1:i,:),1);
    y =  sum(partsI1(1:i,:),1);
    y(y>0) = y(y>0)./totI13D(y>0);
    plot(x(y>0),y(y>0), 'k', 'LineWidth', 2)
    h = fill([x(y>0) x(end) x(1)], [y(y>0) 0 0], cmap(i,:),'EdgeColor','k', 'LineWidth', 2); 
end

axis([x(1) x(end) 0 1.6]);
% xlabel('time','fontweight','bold','fontsize',fsize);
% ylabel('relative intensity', 'fontweight','bold','fontsize',fsize);
set(gca,'color','none')
set(gca, 'LineWidth', 2);
set(gca,'FontSize', fsize)
set(gca,'FontWeight', 'bold')
pbaspect([2 1 1]);
hold off

%%
% graph of area over time
load('/Users/idse/Dropbox/flows/flows_shared/data/zebrafish/huisken/video9/ilt/TracksCoMove.mat')
tmax = 129;
Aproper = zeros([1 tmax]);
Anaive = zeros([1 tmax]);
chidx = 4;

for tidx = 1:tmax
    
    trackc = [];
    for trackidx = 1:numel(TracksCoMove)
        trackc = cat(1, trackc, TracksCoMove(trackidx).centres(tidx,:));
    end
    % order the points by hand, cause I'm a biologist
    bla = trackc(3,:);
    trackc(3,:) = trackc(5,:);
    trackc(5,:) = bla;
    trackc = cat(1,trackc,trackc(1,:));
    
    chart = fullSOI.atlas(tidx).charts{chidx};
    
    mask = poly2mask(trackc(:,1) - stDall(tidx,2), trackc(:,2)- stDall(tidx,1),...
        chart.domain.gridSize(2),  chart.domain.gridSize(1));
    
    X = fullSOI.embedding(tidx).patches{chidx}.apply;
    geom = GaussGeometry(X{1}, X{2}, X{3}, 0);
    
    Aproper(tidx) = sum(geom.dA(mask));
    Anaive(tidx) = sum(mask(:));
end

%%
times = (1:tmax)*3;
plot(times, Aproper, 'LineWidth', 2)
fsize = 26;
set(gca, 'LineWidth', 2);
set(gca,'FontSize', fsize)
set(gca,'FontWeight', 'bold')
axis([1 times(end) 0 1.5*max(Aproper)]);
% hold on
% plot(times, Anaive)
% hold off


%%
fname =  '/Users/idse/Dropbox/flows/flows_shared/data/zebrafish/huisken/video9/conformalHeartMar23/fields/data_MIP/conformal_4_index/conformal_4/adjust1.tif';
info = imfinfo(fname);
num_images = numel(info);
tidx = 19;
A = imread(fname, tidx);
figure, 
imshow(A>0,[])
    
hold on 
trackc = [];
for trackidx = 1:numel(TracksCoMove)
    trackc = cat(1, trackc, TracksCoMove(trackidx).centres(tidx,:));
end
% order the points by hand, cause I'm a biologist
bla = trackc(3,:);
trackc(3,:) = trackc(5,:);
trackc(5,:) = bla;
trackc = cat(1,trackc,trackc(1,:));
plot(trackc(:,1), trackc(:,2),'-g');
%plot(trackc(:,1) - stDall(tidx,1), trackc(:,2)- stDall(tidx,2),'-g');
hold off

%%

saveDir = '/Users/idse/Dropbox/flows/flows_shared/data/zebrafish/huisken/video9/conformalHeartMar23';
%saveDir = '/Users/idse/conformalHeartSIP2';
tic
fullSOI = surfaceAnalysis.SurfaceOfInterest(saveDir);
toc


