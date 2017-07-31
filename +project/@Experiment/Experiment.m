classdef Experiment < handle_light
    % master class of each dataset acting as interface to detector, fitter
    % and SOI
    % a project should have its own directory
    
    %---------------------------------------------------------------------
    % license
    %---------------------------------------------------------------------

    % Copyright 2015 Idse Heemskerk and Sebastian Streichan
    %
    % This file is part of ImSAnE.
    % 
    % ImSAnE is free software: you can redistribute it and/or modify
    % it under the terms of the GNU General Public License as published by
    % the Free Software Foundation, either version 3 of the License, or
    % (at your option) any later version.
    % 
    % ImSAnE is distributed in the hope that it will be useful,
    % but WITHOUT ANY WARRANTY; without even the implied warranty of
    % MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    % GNU General Public License for more details.
    % 
    % You should have received a copy of the GNU General Public License
    % along with ImSAnE.  If not, see <http://www.gnu.org/licenses/>.
    %
    % We make an explicit exception giving permission to link this code
    % with GPL-incompatible matlab facilities.
    
    %---------------------------------------------------------------------
    % properties
    %---------------------------------------------------------------------
    
    properties (SetAccess = protected)

        % struct containing file names and location
        fileMeta = struct(...       
            'directory',    [],...      % project directory (full path)
            'dataDir',      [],...      % data directory (full path)
            'filenameFormat', [],...    % fprintf type format spec of file name
            'timePoints', 	[],...      % list of times available
            'series',   1,...
            'stackSize', [],...         % [xSize ySize zSize]
            'nChannels', [],...         % [xSize ySize zSize cSize]
            'imageSpace', [],...        % uint16 etc, defined in Stack class
            'stackResolution',[],...
            'swapZT',0);
        
        % struct containing imsane metadata
        expMeta = struct(...
            'description', [],...       % string describing the data set
            'channelsUsed',[],...       % the channels used, e.g. [1 3] for RG
            'channelColor',[],...       % mapping from element in channels used to RGB = 123
            'dynamicSurface', [],...    % Not implmented yet, future plan: boolean, false: static surface
            'detectorType', [],...      % name of detector class
            'fitterType', [],...        % name of fitter class
            'fitTime',[],...            % time point of fit for static surfaces
            'jitterCorrection',0);     % Boolean, false: No fft based jitter correction       
        
        
        detectOptions % option structure arrays of time for detector
        fitOptions    % option structure arrays of time for fitter  
        SOIOptions    % option structure arrays of time for SOI
        
        
        ROI           % region of interest object
        fittedParam   % parameters fitted through fitter
        
        
        currentTime   % the current time point
        stack         % the currently loaded stack
        detector      % the currently loaded detector 
        fitter        % the currently loaded fitter

        SOI           % the surface of interest. It has normal shift and time dpendence in its fields
        
    end
    
    properties (Dependent)
        
        metadataFile; % filename of meta data file
        nTimePoints;  % number of time points
        optionIndex;  % acces metadata structure element
        optionLength; % get length of metadata structure 
        currentROI;   % current region of interest
    end
    
    %---------------------------------------------------------------------
    % public methods
    %---------------------------------------------------------------------

    methods
        
        %------------------------------------------------------
        % constructor
        %------------------------------------------------------
        
        function this = Experiment(varargin)
            % constructor 
            %
            % Experiment()
            % Experiment(projectDir)
            % Experiment(projectDir,dataDir);
            
            % experiment directory
            if nargin > 0 && nargin < 2
                this.fileMeta.directory = varargin{1};
            elseif nargin == 2
                this.fileMeta.directory = varargin{1};
                this.fileMeta.dataDir   = varargin{2};
            else
                this.fileMeta.directory = uigetdir;
            end
            
            % debugging info should end up in a subdir
            global debugOutput;
            debugOutput.dir = fullfile(this.fileMeta.directory, 'debugOutput');
            if ~exist(debugOutput.dir, 'dir'), mkdir(debugOutput.dir); end;
                        
            if ~exist(this.metadataFile, 'file') 
                
                % disp('No metadata file found, call setFileMeta() and setExpMeta()');
            else
            
                % else if metadata is found
                disp(['Loading metadata from file' this.metadataFile]);
                meta = load(this.metadataFile);
                % If dataDir passed: overwrite the dataDir in fileMeta 
                meta.fileMeta.dataDir = this.fileMeta.dataDir; 
                this.setFileMeta(meta.fileMeta);
                this.setExpMeta(meta.expMeta);
                
              
                
                % set current time to first time
                if ~isempty(this.expMeta.fitTime)
                    ind = this.tIdx(this.expMeta.fitTime);
                else
                    ind = 1;
                end
                this.currentTime = this.fileMeta.timePoints(ind);
                disp('In Experiment, set the current time based on the saved file meta');
                
                
                % no setter at this level because the fields of these
                % structures are not universal
                this.detectOptions = meta.detectOptions;
                this.fitOptions = meta.fitOptions;
                this.SOIOptions = meta.SOIOptions;
                
                this.ROI = meta.ROI;
                this.fittedParam = meta.fittedParam;
                
                % reset fitter and detector create new fitter detector
                % object and set their options to those in this class
                this.resetDetector();
                this.resetFitter();
            end
        end
        
        %------------------------------------------------------
        % initialize new project (called after setting expMeta, fileMeta)
        %------------------------------------------------------       
                
        function initNew(this)
            % Initialize a new project
            % 
            % initNew()
            %
            % Given the experiment and fileMeta generate a new project
            % with settings according to fileMeta and expMeta.
            % 
            % 
            % see also setFileMeta and setExpMeta
            
            if isempty(this.expMeta) || isempty(this.fileMeta)
                error('Set expMeta and fileMeta before initNew');
            end 

            % set current time to first time
            this.currentTime = this.fileMeta.timePoints(1);
            
            % channels used default is one
            if isempty(this.expMeta.channelsUsed)
                this.expMeta.channelsUsed = 1;
            end
            
            if isempty(this.expMeta.channelColor)...
                        && numel(this.expMeta.channelsUsed) > 1
                this.expMeta.channelColor = [1 2 3];
            end

            % load the stackSize 
            hasbf = exist('bioformats_package.jar','file') ||...
                    exist('loci_tools.jar','file');
            justMeta = 1;
            if hasbf 
                disp('using bioformats to get metadata');
                this.loadStackBioformats(justMeta);
            else
                warning('did not find bioformats, trying without');
                assert(~isempty(this.fileMeta.nChannels),...
                    'fileMeta.nChannels needs to be specified');
                this.loadStack(justMeta);
            end
            
            % initialize detectOptions array to defaul value of
            % detector
            this.setDetectOptions();
            
            % same for fitOptions and SOIOptions
            this.setFitOptions();
            this.setSOIOptions();
        end
        
        %------------------------------------------------------
        % save current project to disc
        %------------------------------------------------------
        
        function save(this)
            % SAVE saves project settings.
            %
            % save()
            %
            % Stores meta data and options of fitters as well as detectors
            % to disc. The file format is a matlab file in the working
            % directory.
            
            fileMeta = this.fileMeta;
            expMeta  = this.expMeta;
            
            detectOptions = this.detectOptions;
            fitOptions = this.fitOptions;
            SOIOptions = this.SOIOptions;
            
            ROI = this.ROI;
            fittedParam = this.fittedParam;

            save(this.metadataFile, 'fileMeta', 'expMeta',...
                'detectOptions', 'fitOptions', 'SOIOptions', 'ROI',...
                'fittedParam');
            
            disp('Project has been saved');
        end
        
        %------------------------------------------------------
        % 3D phase: load data, fit surface
        %------------------------------------------------------
        
        function loadTime(this, t)
            % set the current time and load stack
            %
            % loadTime(t)
            %
            % sets the current time to t, loads the stack 
            % and creates the detector and fitter for that time
            if ~ismember(t, this.fileMeta.timePoints)
                error(['The specified time point ' num2str(t)... 
                    ' is not in the list of available time points.']);
            end
            this.currentTime = t;
            this.expMeta.fitTime = t;
            this.loadStack();
            
            if this.expMeta.dynamicSurface
                this.resetDetector();
                this.resetFitter();
            end
        end
        
        %------------------------------------------------------
        
        function loadStack(this, varargin)
            % Load stack from disc into project.
            %
            % loadStack()
            % loadStack(justMeta)
            %
            % Based on the metadata, load the stack at the current time
            % point using bioformats library. 
            % 
            % see also loadTime
            
            % use bioformats if posssible
            hasbf = exist('bioformats_package.jar','file') ||...
                    exist('loci_tools.jar','file');
            if hasbf
                this.loadStackBioformats(varargin{:});
                return;
            end
            
            fileName = sprintf(this.fileMeta.filenameFormat, this.currentTime);
            fullFileName = fullfile(this.fileMeta.dataDir, fileName);
            tmp = imfinfo(fullFileName);
            nImages = numel(tmp);
            tmp = tmp(1);
            isRGB = strcmp(tmp.ColorType,'truecolor');
            
            if nargin==2
                
                xSize = tmp.Width;
                ySize = tmp.Height;
                if isRGB
                    zSize = nImages;
                    assert(this.fileMeta.nChannels==3,...
                        'your data is RGB, fileMeta.nChannels should equal 3');
                else
                    zSize = nImages / this.fileMeta.nChannels;
                end
                this.fileMeta.stackSize = [xSize ySize zSize];
                    
                return
            end

            xSize = this.fileMeta.stackSize(1);
            ySize = this.fileMeta.stackSize(2);
            zSize = this.fileMeta.stackSize(3);
            nChannels = this.fileMeta.nChannels;
            
            if strcmp(tmp.ColorType,'grayscale')
                assert(nImages == zSize*nChannels,...
                     'fileMeta.nChannels is not consistent with data');
            elseif isRGB
                assert(nImages == zSize,...
                     'fileMeta.nChannels is not consistent with data');
            end
            
            % number of channels used
            nChannelsUsed = numel(this.expMeta.channelsUsed);
            
            % read the data
            ticID = tic;
            data = zeros([ySize xSize zSize nChannelsUsed], 'uint16');
            
            for i = 1:nImages
                
                im = imread(fullFileName, i);
                
                if isRGB
                    data(:,:, i, :) = im(:,:,this.expMeta.channelsUsed);	
                else
                    zidx = 1 + floor((i-1) / nChannels);
                    cidx = 1 + mod(i-1, nChannels);
                    if sum(this.expMeta.channelsUsed == cidx)
                        data(:,:, zidx, this.expMeta.channelsUsed == cidx) = im;
                    end
                end
                
                % progress indicator
                debugMsg(1,'.');
                if rem(i,80) == 0
                    debugMsg(1,'\n');
                end
            end

            debugMsg(1,'\n');
            dt = toc(ticID);
            debugMsg(1,['dt = ' num2str(dt) '\n']);

            % store data in stack object
            
            this.stack = surfaceDetection.Stack(data, this.fileMeta.stackResolution);
            this.stack.setChannelColor(this.expMeta.channelColor);
            this.stack.setDescription(this.expMeta.description);
            
            % store imageSpace in metadata
            this.fileMeta.imageSpace = this.stack.image.image;
        end
        
        function loadStackBioformats(this, varargin)
            % Load stack from disc into project.
            %
            % loadStack()
            % loadStack(justMeta)
            %
            % Based on the metadata, load the stack at the current time
            % point using bioformats library. 
            % 
            % see also loadTime
            
            if nargin==2
                justMeta = varargin{1};
            else
                justMeta = 0;
            end
            
            fileName = sprintf(this.fileMeta.filenameFormat, this.currentTime);
            fullFileName = fullfile(this.fileMeta.dataDir, fileName);
            
            % load the Bio-Formats library into the MATLAB environment
            autoloadBioFormats = 1;
            status = bfCheckJavaPath(autoloadBioFormats);
            assert(status, ['Missing Bio-Formats library. Either add loci_tools.jar '...
                'to the static Java path or add it to the Matlab path.']);

            debugMsg(1, ['Using bioformats version ' char(loci.formats.FormatTools.VERSION) '\n']);

            
            debugLevel = getpref('ImSAnE', 'msgLevel');
            
            % initialize logging
            if debugLevel > 0
                loci.common.DebugTools.enableLogging('INFO');
            else
                loci.common.DebugTools.enableLogging('OFF');
            end
            
            fullFileName
            if ~exist(fullFileName,'file')
                error('raw data file specified through dataDir and fileMeta.filenameFormat does not exist');
            end
            
            % Get the channel filler
            debugMsg(1, ['loading file : ' fullFileName '\n']);
            r = bfGetReader(fullFileName);
            r.setSeries(this.fileMeta.series-1);
            
            % corrupted metadata can give the wrong stack resolution and
            % cause problems, we should set resolution by hand

            if this.fileMeta.swapZT == 0
                stackSize = [r.getSizeX(), r.getSizeY(), r.getSizeZ(), r.getSizeT()];
            else
                stackSize = [r.getSizeX(), r.getSizeY(), r.getSizeT(), r.getSizeZ()];
            end
            debugMsg(2, ['stack size (xyzt) ' num2str(stackSize) '\n']);

            xSize = stackSize(1);
            ySize = stackSize(2);
            zSize = stackSize(3);
            
            % number of channels
            nChannels = r.getSizeC();
            
            nTimePts = stackSize(4);
            
            % update stack size in metadata 
            this.fileMeta.stackSize = [xSize, ySize, zSize];
            this.fileMeta.nChannels = nChannels;

            % if only reading stack size, stop here-----------------------
            if justMeta, return; end;
            
            nChannelsUsed = numel(this.expMeta.channelsUsed);

            if isempty(nChannelsUsed)
                error('nChannelsUsed is empty');
            end
            if any(this.expMeta.channelsUsed > nChannels)
                error('channelsUsed specifies channel larger than number of channels');
            end
            
            debugMsg(1, ['Loading stack for time: ' num2str(this.currentTime) '\n']);
                      
            % this is a stupid workaround
            % basically, even though you save each time point in a separate
            % file with the bioformats exporter, the metadata tells the
            % reader that it is part of a video and it will try to read all
            % the files as a single series
            % this is just a warning to the user, in the loop below is a
            % conditional to prevent it from reading all times
            if nTimePts > 1 
                debugMsg(1, 'More than one time point in metadata!\n');
            end
            
            % read the data
            ticID = tic;
            data = zeros([ySize xSize zSize nChannelsUsed], 'uint16');
            
            for i = 1:r.getImageCount()

                ZCTidx = r.getZCTCoords(i-1) + 1;
                
                % in the fused embryo data coming out of the python script,
                % Z and T are swaped. In general this isn't the case, thus
                % introduce a file metaField swapZT
                if this.fileMeta.swapZT == 0
                    zidx = ZCTidx(1);
                    tidx = ZCTidx(3);
                else 
                    zidx = ZCTidx(3);
                    tidx = ZCTidx(1);
                end
                cidx = ZCTidx(2);

                % see above: if there is only one timepoint all the planes
                % should be read, if there are multiple timepoints, only
                % the correct time shouldbe read
                if nTimePts == 1 || (nTimePts > 1 && this.currentTime == tidx-1)
                    
                    debugMsg(1,'.');
                    if rem(i,80) == 0
                        debugMsg(1,'\n');
                    end

                    dataCidx = find(this.expMeta.channelsUsed == cidx);
                    if ~isempty(dataCidx)
                        data(:,:, zidx, dataCidx) = bfGetPlane(r, i);
                    end
                end
            end

            debugMsg(1,'\n');
            dt = toc(ticID);
            debugMsg(1,['dt = ' num2str(dt) '\n']);

            % store data in stack object
            this.stack = surfaceDetection.Stack(data, this.fileMeta.stackResolution);
            this.stack.setChannelColor(this.expMeta.channelColor);
            this.stack.setDescription(this.expMeta.description);
            
            % store imageSpace in metadata
            this.fileMeta.imageSpace = this.stack.image.image;
        end
        
        %------------------------------------------------------
        
        function rescaleStackToUnitAspect(this)
            % Rescale the stack dimensions such that the resulting aspect
            % ratio is 1.
            % 
            % rescaleStackToUnitAspect()
            
            debugMsg(1,'Rescaling stack to unit aspect ratio\n');
            this.stack = this.stack.rescaleToUnitAspect();
            this.fileMeta.stackSize = this.stack.imageSize;
            
            % we will not store the rescaled stack we won't update the
            % experiment level resolution
            %this.fileMeta.stackResolution = this.stack.resolution;
        end
        
        %------------------------------------------------------
        % convert time to timePoint array index
        %------------------------------------------------------
        
        function idx = tIdx(this, t)
            % convert time to timePoint array index
            %
            % tIdx(t)
            %
            
            if ~ismember(t, this.fileMeta.timePoints)
                error(['The specified time point ' num2str(t)... 
                    ' is not in the list of available time points.']);
            end
            
            idx = find(this.fileMeta.timePoints == t);
        end
        
        
        %------------------------------------------------------
        % Setters for fileMeta and expMeta and fields
        %------------------------------------------------------       
        
        function setFileMeta(this, varargin)
            % set file metadata structure
            %
            % setFileMeta(fileMeta)
            % setFileMeta()
            %
            % User Input required without arguments
            %
            % Options are:    
            % 'directory'         project directory (full path)
            % 'dataDir'           data directory (full path)
            % 'filenameFormat'    fprintf type format spec of file name
            % 'timePoints' 	      Vector of times available
            % 'stackSize'         [xSize ySize zSize]
            % 'imageSpace'        uint16 etc, defined in Stack class
            % 'stackResolution'   Vector containing resolution of pixels
            % 'swapZT'            In some image files, time and z are
            %               swapped
            
            debugMsg(1, 'Experiment.setFileMeta()\n');
            
            if nargin == 2
                fMeta = varargin{1};
                
                if isfield(fMeta, 'dataDir')...
                            && exist(fMeta.dataDir, 'dir') == 0
                    msg = 'Select directory where data is located';
                    disp(msg);
                    % Default search path is the home directory;
                    fMeta.dataDir = uigetdir('~/', msg);
                    %error('data directory specified does not exist');
                end

            else
                fMeta = struct();
        
                msg = 'Select directory where data is located';
                disp(msg);
                fMeta.dataDir = uitgetdir(this.directory, msg);
                
                disp('WARNING: setFileMeta without arguments is unfinished');
            end
            
            % now actually set the fileMeta structure
            fields = fieldnames(this.fileMeta);

            for i = 1:numel(fields)
                if isfield(fMeta, fields{i})
                    this.fileMeta.(fields{i}) = fMeta.(fields{i});
                else
                    debugMsg(2, ['fileMeta.' fields{i} ' not provided, stays the same\n']);
                end
            end
        end
        
        %------------------------------------------------------
        
        function setExpMeta(this, varargin)
            % set experimental metadata structure
            %
            % setExpMeta(expMeta)
            % setExpMeta()
            %
            % User Input required without arguments
            %
            % Options are:    
            % 'description'      string describing the data set
            % 'channelsUsed'     the channels used, e.g. [1 3] for RG
            % 'channelColor'     mapping from element in channels used to RGB = 123
            % 'dynamicSurface'   Not implmented yet, future plan: boolean, false: static surface
            % 'detectorType'     name of detector class
            % 'fitterType'       name of fitter class
            % 'fitTime'          time point of fit for static surfaces
            % 'jitterCorrection' Boolean, false: No fft based translation
            %               estimation of images 
            
            debugMsg(1, 'Experiment.setExpMeta()\n');
            
            if nargin == 2
                
                xMeta = varargin{1};
                
                % check detector type
                if isfield(xMeta, 'detectorType') &&... 
                                    exist(xMeta.detectorType, 'class') == 0
                    error('detector type specified does not exist');
                end
                % check fitter type
                if isfield(xMeta, 'fitterType') &&... 
                                    exist(xMeta.fitterType, 'class') == 0
                    error(['fitter type specified: ' xMeta.fitterType ' does not exist']);
                end
                
            else
                xMeta = struct();
                disp('WARNING: setExpMeta without arguments is unfinished');
            end

            % update provided fields (this way fields cannot be added or
            % removed from the struct)
            fields = fieldnames(this.expMeta);
            
            for i = 1:numel(fields)
                if isfield(xMeta, fields{i})
                    this.expMeta.(fields{i}) = xMeta.(fields{i});
                else
                    debugMsg(1, ['expMeta.' fields{i} ' not provided, stays the same\n']);
                end
            end
        end
        
        %------------------------------------------------------
        % detector wrapper
        %------------------------------------------------------       
        
        function detectSurface(this, seed)
            % Run surface detector and create a pointCloud object
            %
            % detectSurface()
            % detectSurface(seed)
            %
            % calling detectSurface runs the surface detector and creates
            % the detector.pointCloud object
            
            if nargin == 1
                seed = [];
            end
            
            if ~isempty(seed)
                this.detector.detectSurface(this.stack, seed);
            else
                this.detector.detectSurface(this.stack);
            end
            
            PCROI = this.detector.pointCloud.ROI;
            
            % set/update the ROI (alignment) at the experiment level
            if isempty(this.ROI)
                if this.expMeta.jitterCorrection
                    newROIarray(1:this.nTimePoints) = struct(PCROI.irrep);
                else
                    newROIarray(1:this.optionLength()) = struct(PCROI.irrep);
                end
                this.ROI = newROIarray;
            end
                
            this.ROI(this.tIdx(this.currentTime)) = PCROI.irrep;     
        end
        
        %------------------------------------------------------
        
        function setDetectOptions(this, varargin)
        	% set options of detector
            %
            % setDetectOptions()
            % setDetectOptions(detectOptions)
            % setDetectOptions(detectOptions, allTimes)
            % 
            % detectOptions :  structure
            % allTimes:     boolean, only this time or all times
            %
            % sets options for current time,
            % without arguments initializes to detector default
            % resets the detector object
            % 
            % see also surfaceDetection.surfaceDetector
            
            debugMsg(1, 'Experiment.setDetectOptions()\n');
            
            if nargin >= 2
                this.detectOptions(this.optionIndex);
                this.detectOptions(this.optionIndex) = varargin{1};
                this.resetDetector();
            else
                % initialize detectOptions array to defaul value of
                % detector
                this.resetDetector();

                newOpts(this.optionLength()) = struct(this.detector.options);
                this.detectOptions = newOpts;
                for i = 1:this.optionLength()
                    this.detectOptions(i) = this.detector.options;
                end
            end
            
            if nargin == 3 && varargin{2}
                
                disp('Setting detect options for all times');
                
                for t = 1:this.nTimePoints
                    this.detectOptions(t) = varargin{1};
                end
            end
        end
        
        %------------------------------------------------------
        
        function resetDetector(this)
            % resets the detector
            %
            % resetDetector()

            debugMsg(1, 'Experiment.resetDetector():\n');
            
            this.detector = eval(this.expMeta.detectorType);
            
            % if not empty, set the options
            if ~isempty(this.detectOptions)
                this.detector.setOptions(this.detectOptions(this.optionIndex));
                debugMsg(2, 'resetDetector: setting detectOptions from experiment\n');
            end
        end
        
        %------------------------------------------------------
        % fitter wrapper
        %------------------------------------------------------       
        
        function fitSurface(this, varargin)
            % Fit the surface
            % 
            % fitSurface()
            % fitSurface(pointCloud)
            % fitSurface(mesh)
            %
            % pointCloud:   PointCloud object
            % mesh:     structure with fields v, f, vn
            %           respectively 3xN vertices, 3xM faces, 3xN vertex
            %           normals
            % 
            % fit the surface and store the resulting fitted parameters 
            % boolean align whether to use alignment
            
            if isempty(this.fitter)
                error('initialize fitter befor running fitSurface');
            end

            % call fitter.fitSurface
            if nargin == 1
                if ~isempty(this.detector.pointCloud)
                    this.fitter.fitSurface(this.detector.pointCloud);
                    this.fitter.setFitDomain(this.stack.image.domain);
                else
                    error('pass point cloud or run detector first');
                end
            else
                this.fitter.fitSurface(varargin{1});
                this.fitter.setFitDomain(this.stack.image.domain);
            end 
            
            % if shift is non-zero, also b evolve
            shift = this.fitter.fitOptions.shift;
            
            if shift ~= 0
                
                this.fitter.normallyEvolve(shift);
                
                % this used to say:
%                 % this nonsense is because normally evolve assumes that the
%                 % shift in fit options is a previous shift and adds its
%                 % shift to that
%                 fitOpts = this.fitter.fitOptions;
%                 shift = fitOpts.shift;
%                 fitOpts.shift = 0;
%                 this.fitter.setFitOptions(fitOpts);
%                 this.fitter.normallyEvolve(shift, fitOpts.normEvolveSS);
                % the last line
                % because the spherelikeFitter takes ss as an argument, but
                % it shouldn't because it can read it from fitoptions so
                % that should be changed
                % (alternative: make it work incorrectly in tpsFitter for
                % compatibility)
            end
            
            % extract fittedParam and store in this.fittedParam
            % create the appropriate structure array first if needed
            if isempty(this.fittedParam)
                
                newFittedParam(this.optionLength()) = struct(this.fitter.fittedParam);
                this.fittedParam = newFittedParam;
            end
            %this.fittedParam(this.optionIndex) = this.fitter.fittedParam;
                
            % fitSurface can change chartResolution based on fit
            this.SOIOptions(this.optionIndex, :) = this.fitter.charts;
            this.fitOptions(this.optionIndex) = this.fitter.fitOptions;
        end
        
        %------------------------------------------------------
        
        function normallyEvolve(this, shift, varargin)
            % normally evolve and update experiment fitted parameters
            %
            % normallyEvolve(shift)
            % normallyEvolve(shift,subSampling)
            % 
            % shift determines how much the surface is normally evolved.
            % Negative is out, positive inwards. subSampling is an integer
            % specifying the degree of subsampling along all directions
            % 
            % see also surfaceFitting.spherelikeFitter.normallyEvolve
            
            this.fitter.normallyEvolve(shift, varargin{:}); 
            
            %this.fittedParam(this.optionIndex) = this.fitter.fittedParam;
            this.SOIOptions(this.optionIndex, :) = this.fitter.charts;
            this.fitOptions(this.optionIndex) = this.fitter.fitOptions;
        end
        
                
        function zEvolve(this, shift)
            % evolve along z and update experiment fitted parameters
            %
            % zEvolve(shift)
            %
            % see also surfaceFitting.tpsFitter.zEvolve
            
            this.fitter.zEvolve(shift); 
            
            this.fittedParam(this.optionIndex) = this.fitter.fittedParam;
            this.SOIOptions(this.optionIndex, :) = this.fitter.charts;
            this.fitOptions(this.optionIndex) = this.fitter.fitOptions;
        end
 
        %------------------------------------------------------
        
          
        function setFitOptions(this, varargin)
            % set options of the fitter
            %
        	% setFitOptions()
            % setFitOptions(fitOptions)
            % setFitOptions(fitOptions, allTimes)
            % 
            % fitOptions :  structure
            % allTimes:     boolean, only this time or all times
            %
            % sets options for current time,
            % without arguments initializes for all times to fitter default
            % resets the fitter
            %
            % see also surfaceFitting.tpsFitter and
            % surfaceFitting.spherelikeFitter
            
            debugMsg(1, 'Experiment.setFitOptions()\n');
            
            if nargin >= 2
                
                % set options of the fitter object
                this.fitter.setFitOptions(varargin{1});
                
                % set the options in experiment equal to those of the
                % fitter
                this.fitOptions(this.optionIndex) = this.fitter.fitOptions;
                this.resetFitter();
                
            else
                % initialize fitOptions array to defaul value of fitter
                this.resetFitter();

                newOpts(this.optionLength()) = struct(this.fitter.fitOptions);
                this.fitOptions = newOpts;
                for i = 1:this.optionLength()
                    this.fitOptions(i) = this.fitter.fitOptions;
                end
            end
            
            % if all times
            if nargin == 3 && varargin{2}
                
                disp('Setting fit options for all times');
                
                for i = 1:this.nTimePoints
                    this.fitOptions(i) = this.fitter.fitOptions;
                end
            end
        end
        
        %------------------------------------------------------
        
        function setSOIOptions(this, varargin)
        	% Set options of the Surface of Interest
            %
            % setSOIOptions()
            % setSOIOptions(SOIOptions)
            % sets options for current time,
            % without arguments initializes to fitter default
            % resets the fitter
            %
            % see also surfaceAnalysis.surfaceOfInterest
            
            debugMsg(1, 'Experiment.setSOIOptions()\n');
            
            if nargin == 2
                this.SOIOptions = varargin{1};
                this.resetFitter();
            else
                % initialize SOIOptions array to defaul value of fitter
                this.resetFitter();

                newOpts(this.optionLength(),:) = struct(this.fitter.charts);
                this.SOIOptions = newOpts;
                for i = 1:this.optionLength()
                    this.SOIOptions(i,:) = this.fitter.charts;
                end
            end
        end

        %------------------------------------------------------
        
        function resetFitter(this)
            % resets the fitter 
            %
            % resetFitter()
            % 
            % create a new object with the same options and alignment

            
            if isempty(this.fileMeta.stackSize) || isempty(this.expMeta.fitterType)
                error('Set fileMeta.stackSize and expMeta.fitterType before resetFitter()');
            end
            
            debugMsg(1, 'Experiment.resetFitter():\n');

            this.fitter = eval([this.expMeta.fitterType '()']);
            
            % if not empty, set the options
            if ~isempty(this.fitOptions)
                debugMsg(2, 'resetFitter: setting fitOptions from experiment\n');
                this.fitter.setFitOptions(this.fitOptions(this.optionIndex));
                
            end
            if ~isempty(this.SOIOptions)
                debugMsg(2, 'resetFitter: setting SOIOptions from experiment\n');
                
                charts = this.SOIOptions(this.optionIndex, :);
                
                for i = 1: size(this.SOIOptions, 2)
                    this.fitter.setChartResolution(charts(i).name, charts(i).stepSize);
                    this.fitter.setDesiredChart(charts(i).name, charts(i).desired);
                end
                
            end
            
            
            % restore fitDomain from stored alignment
            if ~isempty(this.currentROI)
                bdry = {round(this.currentROI.xpRange),...
                    round(this.currentROI.ypRange), round(this.currentROI.zpRange)};
                embeddingSpace = diffgeometry.FiniteSet('alignedEmbeddingSpace',...
                                                            bdry, [1 1 1]);
                this.fitter.setFitDomain(embeddingSpace);
            end
        end

        %------------------------------------------------------
        % ROI stuff
        %------------------------------------------------------       
        
        function updateROI(this)
            % sync experiment ROI with detector ROI
            %
            % updateROI()
            % 
            % see also surfaceDetection.pointCloud
            
            tIdx =  this.tIdx(this.currentTime);
            this.ROI(tIdx) = this.detector.pointCloud.ROI.irrep;
        end
        
        function determineROIfromFit(this, varargin)
            % use existing fit to determine ROI
            %
            % determineROIfromFit()
            % determineROIfromFit(pc)
            %
            % This method is mainly used with spherelikeFitter and can't be
            % used without having fitted the data first.
            % As point clouds are oftentimes inhomogenous and noisy, we
            % determine the ROI around the SOI based on an initial fit.
            % Based on the fit, we determine the orientation of the
            % surface, to then generate a more accurate fit to the data
            % such that the major axis coincides with the z axis
            % 
            % see also surfaceFitting.spherelikeFitter
            debugMsg(1,'Experiment.determineROI()');
                            
            if isempty(this.fitter.fittedPoints)
               error('Cannot call setROI without arguments without fitting first'); 
            end
            
            % generate the alignment based on this fit
            embedding = this.fitter.fittedPoints;

            fitPoints = [embedding{1}(:), embedding{2}(:), embedding{3}(:)];
            fitPointCloud = surfaceDetection.PointCloud(fitPoints);

            margin = 10;
            fitPointCloud.determineROI(margin);
            relativeROI = fitPointCloud.ROI;

            if nargin == 1
                PC = this.detector.pointCloud;
            else
                PC = varargin{1};
            end

            PC.setROI(relativeROI);
            
            % update the ROI for all timepoints if its a dynamic one. 
            if this.expMeta.dynamicSurface == 1
                tIdx =  this.tIdx(this.currentTime);
                this.ROI(tIdx) = PC.ROI.irrep;
            elseif this.expMeta.jitterCorrection == 1
                for i = 1 : length(this.ROI)
                    this.ROI(i) = PC.ROI.irrep;
                end
            else
                this.ROI = struct('rotation',[],'translation',[],...
                    'xpRange',[],'ypRange',[],'zpRange',[]);
                this.ROI(1) = PC.ROI.irrep;
            end
        end
        
        %------------------------------------------------------
        % JitterCorrectROI
        %------------------------------------------------------ 
        
        function jitterCorrectROI(this,ROIT,TJitter)
            % use jitter correction to update ROI. 
            % 
            % jitterCorrectROI(ROIT,TJitter); 
            %
            % ROIT is the current ROI.  
            % TJitter is a translation obtained by 
            % see also surfaceDetection.stack.determineRelativeShift.
            % ROIT will be translated by TJitter. 

            T0 = ROIT.translation;
            
            ROIout = surfaceDetection.RegionOfInterest(ROIT.rotation, TJitter*T0);
            ROIout.setRanges(ROIT.xpRange, ROIT.ypRange, ROIT.zpRange);
           
            tIdx =  this.tIdx(this.currentTime);
            this.ROI(tIdx) = ROIout.irrep;
        end
      
        %------------------------------------------------------
        
        function batchProcess(this, timePoints, options)
            % Loop through all the steps upto pullback stack for a series of
            % timepoints
            % 
            % batchProcess(timePoints)
            % batchProcess(timePoints, options)
            %
            % options: structure with fields
            %   - unitAspect : boolean
            % 
            % Loop through all the batch processing time points and
            % depending on the options provided extract the surface in the
            % data or jitter correct(translate) the surface and then
            % extract the image data.
            % 
            % see also surfaceDetection.stack
            
            if nargin < 3
                options = struct('unitAspect', true);
            end
            
            if ~this.expMeta.dynamicSurface
                
                if this.expMeta.jitterCorrection == 1
                    % do jitter correction 
                    
                    % prepare jitter correction; Subsample current time;
                    FourierDownSample = 4;
                    stack_ss = this.stack.reduceStack(FourierDownSample);
                    % get the current ROI
                    ROIT = this.currentROI;
                    
                    for t = 1 : length(timePoints)
                        
                        % load current time point
                        this.loadTime(timePoints(t));                 
                        this.rescaleStackToUnitAspect();
                        
                        % compare previous subsampled with current
                        TJitter = this.stack.determineRelativeShift(stack_ss,...
                            FourierDownSample);

                        % Translate the current ROI by TJitter.
                        this.jitterCorrectROI(ROIT,TJitter);
                        
                        % Sub-sample current stack, and overwrite old
                        % sub-sampled stack;
                        stack_ss = this.stack.reduceStack(FourierDownSample);

                        this.SOI.pullbackStack(this.stack, this.currentROI, this.currentTime);
                        ROIT = this.currentROI;
                    end
                    
                else
                
                    for t = 1:length(timePoints)

                        this.loadTime(timePoints(t));
                        this.rescaleStackToUnitAspect();
                        this.SOI.pullbackStack(this.stack, this.currentROI, this.currentTime);
                    end  
                end
                
            else
                seed = [];
                
                for t = 1:length(timePoints)

                    debugMsg(1, ['Processing time point ' num2str(timePoints(t)) '\n']);
                    
                    % load stack and rescale
                    this.loadTime(timePoints(t));
                    if options.unitAspect
                        this.rescaleStackToUnitAspect();
                    end

                    % surface detection
                    this.setDetectOptions(this.detectOptions(t));
                    this.detectSurface(seed);

                    % surface fitting
                    this.setFitOptions(this.fitOptions(t));
                    this.fitSurface();

                    shift = this.fitOptions(t).shift;
                    if shift > 0
                        % THIS WON'T WORK FOR THE NON-PLANAR
                        this.zEvolve(shift);
                    end

                    % SOI and pullbacks
                    this.fitter.populateSOI(this.SOI, this.currentTime);
                    this.SOI.pullbackStack(this.stack, this.currentROI, this.currentTime);
                    
                    if options.seeded 
                        seed = this.fitter.generateSeed();
                    end 
                end  
            end 
        end
        
        %------------------------------------------------------
        % 2D phase: generate SOI, load pullbacks from file
        %------------------------------------------------------
        
        function generateSOI(this)
            % Generate a SOI object after fitting the point cloud.
            %
            % generateSOI()
            % 
            % We do not load fitted parameters because in general
            % that would  risk having fitted parameters that do not
            % match fit options, but when reloading from file they will
            % match and we will set the fitted parameters
            %
            % see also surfaceAnalysis.surfaceOfInterest and
            % surfaceFitting.spherelikeFitter.populateSOI

            if isempty(this.fitter)
                error('initialize fitter before generateSOI');
            end
            
            if isempty(this.fittedParam)
                error('fittedParam empty: fitSurface first');
            %else
            %    this.fitter.setFittedParameters(this.fittedParam(this.optionIndex));
            end

            embeddingSpace = this.fitter.fitDomain;
            if isempty(embeddingSpace)
                debugMsg(1, 'generateSOI: fitter.fitDomain empty, assuming stack domain\n')
                this.fitter.setFitDomain(this.stack.image.domain);
                embeddingSpace = this.fitter.fitDomain;
            end

            dataSpace = this.fileMeta.imageSpace;
            
            % construct surface of interest object;
            this.SOI = surfaceAnalysis.SurfaceOfInterest(dataSpace,...
                embeddingSpace, this.fileMeta.timePoints,...
                this.expMeta.dynamicSurface);

            % populate SOI
            if this.expMeta.dynamicSurface
                this.fitter.populateSOI(this.SOI, this.currentTime);
            else
                this.fitter.populateSOI(this.SOI);
            end
        end
        
        %------------------------------------------------------
        % dependent variables 
        %------------------------------------------------------       
        
        function metadataFile = get.metadataFile(this)
            % specify the name of the metadata file
            metadataFile = fullfile(this.fileMeta.directory, 'ImSAnE_MetaData.mat');
        end
        
        function nTimePoints = get.nTimePoints(this)
            % get the number of time points
            nTimePoints = numel(this.fileMeta.timePoints);
        end
        
        function idx = get.optionIndex(this)
            % convenient way to acces the right metadata structure element
            % depending on whether the surface is dynamic
            if this.expMeta.dynamicSurface %|| this.expMeta.jitterCorrection
                for idx = 1:this.nTimePoints
                    if this.fileMeta.timePoints(idx) == this.currentTime
                        return; 
                    end
                end
            else
                idx = 1;
            end
        end
        
        function L = get.optionLength(this)
            % convenient way to get length of metadata structure 
            % depending on whether the surface is dynamic
            if this.expMeta.dynamicSurface %|| this.expMeta.jitterCorrection
                L = this.nTimePoints;
            else
                L = 1;
            end
        end
        
        function ROI = get.currentROI(this)
            % get current alignment as an object

            if isempty(this.ROI)
                ROI = [];
            else
                if this.expMeta.dynamicSurface == 1 || this.expMeta.jitterCorrection == 1
                    irrep = this.ROI(this.tIdx(this.currentTime));
                    ROI = surfaceDetection.RegionOfInterest(irrep.rotation, irrep.translation);
                    ROI.setRanges(irrep.xpRange, irrep.ypRange, irrep.zpRange);
                else
                    irrep = this.ROI(1);
                    ROI = surfaceDetection.RegionOfInterest(irrep.rotation, irrep.translation);
                    ROI.setRanges(irrep.xpRange, irrep.ypRange, irrep.zpRange);
                end
            end
        end
    end
end