%% ImSAnE: ZEBRAFISH BEATING HEART TUTORIAL

% In this tutorial we analyze a beating zebrafish heart. This involves
% extraction of a dynamic surface based on surface meshing, surface 
% partitioning and subsequent parametrization. We will also learn how to
% reorganize multi-layer data in a stack, and measure local tissue
% deformation by following a select group of nuclei. 
%
% Additional software required: Meshlab (http://meshlab.sourceforge.net). 
%
% The tutorial starts by defining the experiment, and commences with
% loading a previously generated pointCloud, that is turned into a mesh
% using meshlab. 
%
% This is the most advanced tutorial, it is recommended to work through the
% basic ones on the drosophila embryo or wing disc first.  

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

%dataDir = '/Volumes/data/shared/Dropbox/tissueflows/data/zebrafish/huisken/video9/';
dataDir    = fullfile(scriptPath, 'rawData');
projectDir = fullfile(scriptPath, 'projectFiles');

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
fileMeta                 = struct();
fileMeta.dataDir         = dataDir;
fileMeta.filenameFormat  = 'VSilent_T%d_merged.ome.tif';
fileMeta.timePoints      = 1:129;
fileMeta.stackResolution = [.45 .45 1];
fileMeta.swapZT          = 0; 

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
expMeta.jitterCorrection = 0; 
expMeta.fitTime          = 1; 
% the detector type is irrelevant here, as we are using an external
% detector
expMeta.detectorType     = 'surfaceDetection.fastCylinderDetector';
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
xp.loadTime(fileMeta.timePoints(1)); 
xp.rescaleStackToUnitAspect();


%% load a point cloud from disk and turn it into a mesh


% Specify the name format of the point cloud and read time point 1 saved in OBJ.  
PCfileFormat = fullfile(projectDir, 'objFiles', 'heartPC%d_ilastik.obj');
PCfile =  sprintf(PCfileFormat,1); 

% Specify the name of the output mesh format and save 
outputMeshFormat = fullfile(projectDir, 'objFiles', 'heartPC%d_ilastik_mesh.ply');
outputMesh = sprintf(outputMeshFormat,1);

% Specify the name and location of the meshlab script on disc. If there is
% no meshlab script, generate one using meshlab as follows:
%
%   Disclaimer: Meshlab is still under development and we noticed
%   differences across versions as well as operating systems. The below
%   description is based on v1.3.3 under MacOS. Generation of the script
%   should be very similar in the other operating systems, but differences
%   can occur across versions. We successfully tested: Win 7: meshlab 1.3.2
%   Ubuntu 12.04: meshlab 1.3.3Beta, MacOS 10.10: meshlab 1.3.3
%   
%   
%
%   Startup meshlab. (v1.3.3 on MacOS 10.10) 
%   In meshalb, import a point cloud 
%       -> File -> import mesh 
%   Select a point cloud on your disc in the dialog that appears.  
%
%   Next, subsample the point cloud by 
%       -> Filters -> sampling -> Poisson-disk Sampling 
%   and tick Base Mesh Subsampling. Default values for the rest are okay. 
%
%   Select the sampled layer under 
%       -> View -> Show Layer Dialog.
%
%   Now compute the normals of the point cloud by 
%       -> Filters -> Point Set -> Compute normals for point set.
%   Default values are okay.
%   
%   Now reconstruct the surface using poisson surface reconstruction at
%   default values
%       -> Filters -> Point Set -> Surface Reconstruction: poisson
%
%   Using different rendering modes, inspect the quality of the
%   reconstructed surface
%       -> Render -> Render Mode -> Wireframe
%  
%   To take on a convention on the orientation of faces (normals pointing
%   in the same direction), use 
%       -> Filters -> Normals, Curvatures and Orientation -> Invert Faces
%       Orientation
%   Untick Force Flip. Now all the faces should be oriented outwards. 
% 
%   To obtain a finer mesh, you can refine the triangles using 
%       -> Filters -> Remeshing, Simplification and Reconstruction ->
%       Subdivision surfaces loop
%   The default values should produce a fine enough mesh for our purposes. 
% 
%   Finally inspect and save the current filter script. 
%       -> Filters -> Show current filter script
%   Following the steps above, you should see the following list of
%   entries: 
%       Poisson-disk Sampling
%       Compute normals for point sets
%       Surface Reconstruction: Poisson
%       Invert Faces Orientation 
%       Subdivision Surfaces: Loop
%
%   Hit save script and navigate to the project directory (projectDir). 
%   Note the file name that you saved it as, e.g. heartFinal1.mlx and enter 
%   it below as a string in the variable myMeshlabScriptName;

myMeshlabScriptName = 'heartFinal1.mlx';
meshlabScript = fullfile(projectDir, myMeshlabScriptName);


% Now we create the mesh using meshlabserver and the just created script.
% Matlab has conflicting library defintions, which we resolved
% under Ubuntu and Windows, not in Mac Os. Here we will describe how to
% batch create meshes based on point clouds stored in a single folder using
% the clipboard and mac terminal. Windows and ubuntu useres can use the
% system command instead. 
% 
% Copy the command sctring into the clipboard. The syntax is -i for input,
% -o for output, -s for script name -om for options, here specifying that
% we whish to obtain the normals. 
command = ['meshlabserver -i ' PCfile ' -o ' outputMesh ' -s ' meshlabScript ' -om vn'];
clipboard('copy', command);
% 
% Open up the terminal, by going to the Finder, 
%       -> Go -> Utilities -> Terminal
%
% in the terminal, hit Command(apple) + v button. If all goes well, the
% meshlab server will start, producing the outputMesh based on the input
% point cloud file. 


%% optionally for MAC: For dynamic surfaces, we may whish to loop over all time points

% bash script that can be called on the command line, that will process all
% point clouds in a folder. 

% STILL TO DO! 


%% Read surface mesh produced by meshlab
% 
meshPermutation = [2,1,3]; % There is an axis permutation in the pointCloud, 
                           % that we need to take care of to have the
                           % surface axes match the data axes. 
                           
% read the output mesh, and crate a mesh structure. 
mesh = read_ply_mod(outputMesh);
mesh.v = mesh.v(:,meshPermutation); 
mesh.vn = mesh.vn(:,meshPermutation);

%% visualize mesh that was created by meshlab
% last argument in trimesh color codes vertices for y-position

trimesh(mesh.f, mesh.v(:,1), mesh.v(:,2), mesh.v(:,3), mesh.v(:,1));
colormap gray
axis equal

%% set parameters for atlas creation

% we build an atlas by dividing the mesh into overlapping submeshes and
% mapping those to the plane
% 
% ImSAnE creates two types of submeshes, both of which are created with a
% seed point on the surface as input.
%
% The first type is the Voronoi submesh:
% Given a number of seedpoints, the surface is divided in "geodesic Voronoi
% cells", meaning the submeshes are defined by all the points being closer
% to one seedpoint than any other, with the distance measured along the 
% surface. An overlap width between submeshes can be chosen by the user.
% By construction the Voronoi submeshes always cover the entire surface.
%
% The second submesh type is the disk submesh:
% This consist of all the points closer to the seedpoint than some user
% specified distance. Such a submesh tends to be less distorted because
% it has a very regular shape and can be chosen to be small, but it is not
% garantueed that choosing some disklike submeshes gives us complete
% coverage of the surface. Therefore, we may want to do both: create a
% Voronoi partition of the surface and create a disklike submesh around
% some point that interest us, such as the valve.

% Here we specify the rough XYZ position of where we want the seedpoint to
% be. These can be obtained by for example using the data cursur in the
% figure panel. 

VorSeedsXinit = [127 338 106; 311 290 110; 332 326 252]; 
diskSeedsXinit = [307 431 185];

% We then find the index into the mesh of the vertex closest
% to these positions.

diskSeeds = pointMatch(diskSeedsXinit, mesh.v);
VorSeeds = pointMatch(VorSeedsXinit, mesh.v);

% % Alternatively, we can just pick random points
% nVorSeeds = 3;
% nDiskSeeds = 1;
% VorSeeds = floor(rand([nVorSeeds 1])*size(mesh.v,1))+1;
% diskSeeds = floor(rand([nDiskSeeds 1])*size(mesh.v,1))+1;


% meshWrapper, as the name suggests, is a surface fitter class but rather
% than constructing the surface representation by fitting, it is a wrapper
% to work with externally created meshes
%
% we pass the mesh to fitSurface and that creates the submeshes that will
% be used for atlas creation
%
% The fitOptions structure should contain the following fields:
%
% VorSeeds:         seeds (vertex index) for the geodesic 
%                   Voronoi cells
% transitionWidth:  half-width of overlap between Voronoi cells
% diskSeeds:        center (vertex index) for disk like submeshes
% diskRadius:       radius of disk like submeshes
% makeTMaps:        boolean, should transition maps be made?

fitOptions = struct('VorSeeds', VorSeeds, 'transitionWidth', 50,...
                    'diskSeeds', diskSeeds, 'diskRadius', 150, 'makeTMaps', false);
xp.setFitOptions(fitOptions);
xp.fitSurface(mesh);

%%
% we can visualize a subset of the overlapping regions (indexing is such
% that disk regions follow Voronoi regions, so 1:3 shows just the Voronoi
% regions in this example
xp.fitter.inspectMesh(1:3);
view([0 0 1]);

%%
% we can also display a single submesh
xp.fitter.inspectMesh(2);

%% Inspect the fit in a cross section
%
% inspect fit and point cloud over a cross section in the data. Dimensions 
% are x,y or z and the value has to be within the corresponding axis range.

inspectOptions = struct('dimension','x', 'value', 300);
xp.fitter.inspectQuality(inspectOptions, xp.detector, xp.stack);

%% Generate the atlas

% We now generate the Surface Of Interest object, which involves building
% the atlas. The only type of chart currently implemented for meshes is the 
% conformal map described in the manuscript, which preserves shape but
% locally distorts size.

xp.generateSOI();

%% Create maps of the surface data
% 
% We call the process of creating the 2D maps of the surface data "pulling
% back the data" (in accordance with standard mathematical terminology)
% Calling the function SurfaceOfInterest.pullbackStack generates the Field 
% objects in SOI.fields that contain these pullbacks.
% The function is called as:
%
% pullbackStack(stack, ROI, time, onionOpts)
%
% Where the arguments are:
%
% stack:    Stack object containing the data
% ROI:      RegionOfInterest object containing affine
%           transformation from embedding coordinates to stack
%           ROI can be left empty: ROI = []
% time:     the time of the pullback, for storing in the SOI
%
% onionOpts: multi-layer pullback options, structure with fields
%   - nLayers:          number of layers
%   - layerDistance:    spacing of layers
%   - sigma:            smoothing of surface for layer
%                       displacement
%   - makeIP:           intensity projection: 'MIP', 'SIP',
%                       'both'
%   - IPonly:           only IP is stored as a field
%                       WARNING: this option will delete
%                       previously stored other layers if they
%                       were made
            
onionOpts = struct('nLayers', 41, 'layerDistance', 1, 'sigma', 1,...
                    'makeIP', 'both', 'IPonly', false);
 
xp.SOI.pullbackStack(xp.stack, [], xp.currentTime, onionOpts);

%% Visualize the atlas in 2D

% the summed intensity projection of each region at the current time
dataField = xp.SOI.getField('data_SIP');
tidx = xp.tIdx(xp.currentTime);
data = dataField(tidx);

% NOTE: For more information on the organization of the data structures, 
% see the supplemental information of the manuscript.

% the color version of the map (pullback) to the plane, stored to be used 
% as texture for the 3D renderinging the next block
color = {};

figure,
for i = 1:numel(data.patches)
    
    % the two channels
    R = mat2gray(data.patches{i}.apply{1});
	G = mat2gray(data.patches{i}.apply{2});
    
    % a little bit of manual adjustment of the lookup table
    R = imadjust(R, [0 0.5]);
    G = imadjust(G, [0 0.8]);
    
    % concatenate to make color image
    color{i} = cat(3,R,G,G);
    
    % make the background white 
    color{i}(color{i}==0) = 1;
    

    imwrite(permute(color{i}, [2 1 3]),['forweb' num2str(i) '.tif'])
    % show the map
    subplot(ceil(numel(data.patches)/2), 2, i)
    imshow(permute(color{i}, [2 1 3]),[],'InitialMagnification',66)
end

%% Visualize the atlas in 3D

% we take the submeshes and displace them along the average normal by a 
% distance set by separation
separation = 140;

figure
hold on;

% the pieces
for i = 1:numel(data.patches)
    
    % the 3D coordinates of a region
    X = xp.SOI.embedding(tidx).patches{i}.apply;
    
    % the mean normal direction
    d = mean(xp.fitter.fittedParam.submeshes{i}.vn);
    d = d./norm(d);
    d = separation*d;
    
    % show each region in 3D displaced along its normal by separation
    surf(X{1} + d(1),X{2} + d(2),X{3} + d(3), color{i}, 'FaceColor','texturemap');  
end

axis equal;
axis off;
shading flat
set(gcf,'color','w');
view(60, -60); % set viewing angle

hold off;

%% Loop for dynamic atlas generation

% For processing a time series we basically repeat the above steps in a
% loop. What is new is that we now constrain the mapping of a number of 
% fixed points on the surface to ensure temporal continuity in the maps.
%
% The use can safely ignore this, but we describe it briefly for
% completeness:
% In each timestep t, for each chart k, we match the seedPoint and a 2nd 
% point in 3D (fixedPtX) between two consecutive timepoints and then demand 
% they are mapped to the same 2D value (fixedPtU). We pass these values
% through fitOptions

% keep track of seed position propagation in time
VorSeedsX{1} = mesh.v(VorSeeds,:);
diskSeedsX{1} = mesh.v(diskSeeds,:);

% point in parametrization to hold fixed next round
fixedPtU = {};
fixedPtX = {};
for k = 1:numel(xp.fitter.fittedParam.submeshes)
    subm = xp.fitter.fittedParam.submeshes{k};
    fixedPtUtmp = sqrt(mean(subm.u{1}.^2)/2);
    fixedPtIdx = pointMatch(fixedPtUtmp, subm.u{1});
    fixedPtX{k} = subm.v(fixedPtIdx,:);
    fixedPtU{k} = subm.u{1}(fixedPtIdx,:);
end

for t = fileMeta.timePoints(2:end) 
    
    tidx = xp.tIdx(t);
    
    % raw data loading
    xp.loadTime(t);
    xp.rescaleStackToUnitAspect();
    
    % Read surface mesh produced by meshlab
    outputMesh = sprintf(outputMeshFormat, t);
    mesh = read_ply_mod(outputMesh);
    mesh.v = mesh.v(:,meshPermutation); 
    mesh.vn = mesh.vn(:,meshPermutation);

    % initialize fitter with overlap + disk
    fitOptions.VorSeeds = pointMatch(VorSeedsX{tidx-1}, mesh.v);
    fitOptions.diskSeeds = pointMatch(diskSeedsX{tidx-1}, mesh.v);
    fitOptions.fixedPtU = fixedPtU;
    fitOptions.fixedPtX = fixedPtX;
    xp.setFitOptions(fitOptions);
    xp.fitSurface(mesh); 
    
    % seeds positions for the next iteration
    VorSeedsX{tidx} = mesh.v(fitOptions.VorSeeds,:);
    diskSeedsX{tidx} = mesh.v(fitOptions.diskSeeds,:);

    % for next time, propagate fixedPtX
    for k = 1:numel(xp.fitter.fittedParam.submeshes)
        V = xp.fitter.fittedParam.submeshes{k}.v;
        fixedPtIdx = pointMatch(fixedPtX{k}, V);
        fixedPtX{k} = V(fixedPtIdx,:);
    end

    % populate SOI
    xp.fitter.populateSOI(xp.SOI, xp.currentTime);

    % Pullback the stack to the desired charts
    xp.SOI.pullbackStack(xp.stack, [], xp.currentTime, onionOpts);
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
saveDir = fullfile(projectDir, 'heartSOI');
options = struct('dir',saveDir,'imwriteOptions',{imwriteOptions},...
                    'make8bit',false);
xp.SOI.save(options)

%% multi layer maps as regular image stacks

% returns a cell array of regular image stack indexed by the patch (region)
% index 

stacks = xp.SOI.multilayer2stack(xp.currentTime);

% to save these stacks to tiff, use:
% stacks = xp.SOI.multilayer2stack(xp.currentTime, projectDir);

% visualize cross section:

pi = 4;     % patch index
xval = 100; % section index
% convert to RGB image with adjusted lookup table and right orientation
R = flipud(imadjust(squeeze(stacks{pi}(:,xval,:,1)))');
G = flipud(imadjust(squeeze(stacks{pi}(:,xval,:,2)))');
imshow(cat(3,R,G,G));

%% measure area spanned  by tracked nuclei

% we will use the maximum intensity projected layers, to inspect
% trajectories of cells on the conformal_4 patch.
data = xp.SOI.getField('data_MIP');
pi = 4; % tracks shown in the paper are from conformal patch 4. 
% Load stored tracks from file.
load('/Volumes/data/shared/Dropbox/tissueflows/data/zebrafish/huisken/video9/ilt/TracksCoMoveAdj.mat')

figure,

for t = fileMeta.timePoints

    tidx = xp.tIdx(t);
    
    % load the color images and display the tracks over the data
    
    im1 = imadjust(mat2gray(data(tidx).patches{pi}.apply{1}));
    im2 = imadjust(mat2gray(data(tidx).patches{pi}.apply{2}));
    colim = cat(3, im1, im2, im2);
    % show the image
    imshow(colim);
    
    hold on 
    trackc = [];
    % tracks are stored in a structure TracksCoMove. Each element contains
    % the centres of the nucleus as NTimes*2 array.
    for trackidx = 1:numel(TracksCoMove)
        trackc = cat(1, trackc, TracksCoMove(trackidx).centres(tidx,:));
    end
    % plot the tracks;
    scatter(trackc(:,1), trackc(:,2), 'ko');
    hold off
    
    % Now we compute the area enclosed by the tracked nuclei.
    
    % specify patch and transform name for proper area calculation.
    patchName     = 'conformal_4_index';
    transformName = 'conformal_4';

    % generate a mask image, that is 1 in the polygon spanned by the nuclei
    % and 0 elsewhere. 
    [X,Y] = meshgrid(1:size(im1,2),1:size(im1,1));
    % Tracks can move and therefore need to be sorted at every time point;
    clockSorted = ClockWiseSort(trackc);
    mask = inpolygon(X,Y,clockSorted(:,1),clockSorted(:,2));

    % is this a bug?? I needed to calculate the induced metric manually.
    xp.SOI.NCalcInducedMetric(tidx);
    
    % direct
    AC1(t) = sum(mask(:));
    % corrected measurement
    AproperC1(t) = xp.SOI.properArea(tidx, mask, transformName);
    
    
end

