classdef Stack < handle_light
    %stackHolder Raw data and metadata, transform after ROI and axis selection
       
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
        
        % set at or directly after construction
        
        image                       % Map object representing the image
        resolution                  % vector of optical resolutions
        description                 % description of data (date, labels, ..)                      
        channelColor                % the color RGB = 123 of each image channel
        
        % set by setIntensityLimits
        
        Ilim                        % cell array of intensity limits
        
        % set by setProjectionMask
        
        projectionMask              %optional projection mask (3 2D masks) 
    end
    
    %---------------------------------------------------------------------
    % dependent and constant properties
    %---------------------------------------------------------------------
    
    properties (Dependent = true)
        
        aspect;       	% ratio of z-resolution to xy-resolution
        imageSize;      % size of the image [x y z] (actual data is [y x z])
    end
    
    properties (Constant = true)
        
        % defines the possible stack types as FiniteSet (e.g. grayscale_16bit)
        imageSpaces = struct(...
            'binary', diffgeometry.FiniteSet('binary', {[0 1]}, 1),... 
            'grayscale_16bit', diffgeometry.FiniteSet('grayscale_16bit',...
                                                {[0 2^(16)-1]}, 1),...
            'grayscale_8bit', diffgeometry.FiniteSet('grayscale_8bit',...
                                                        {[0 255]}, 1),...
            'RGB_8bit', diffgeometry.FiniteSet('RGB_8bit',...
                                {[0 255], [0 255], [0 255]}, [1 1 1]),...
            'RGB_16bit', diffgeometry.FiniteSet('RGB_16bit',...
                    {[0 2^(16)-1], [0 2^(16)-1], [0 2^(16)-1]}, [1 1 1]),...
            'twoColor_8bit', diffgeometry.FiniteSet('twoColor_8bit',...
                                {[0 255], [0 255]}, [1 1]),...
            'twoColor_16bit', diffgeometry.FiniteSet('twoColor_16bit',...
                                {[0 2^(16)-1], [0 2^(16)-1]}, [1 1])...
            );
    end
    
	%---------------------------------------------------------------------
    % protected methods
    %---------------------------------------------------------------------
    
    methods (Access = protected, Static)
       
        function [type nChannels] = detectStackType(stack)
            % DETECTSTACKTYPE Make a best guess of the stack type
            %
            % [type nChannels] = detectStackType(stack)
            %
            % stack:        4D stack with the 4 dimension color (or 3D for
            %               greyscale) or a cell array, e.g. {R,G,B}
            % type:         an entry in imageSpaces
            % nChannels:    number of channels
            
            % determine number of channels
            if isnumeric(stack) && ndims(stack) == 3 
                nChannels = 1;
                
            elseif isnumeric(stack) && ndims(stack) == 4
                nChannels = size(stack,4);
                
            elseif iscell(stack)
                nChannels = length(stack);
                
            else
                debugMsg(2, 'cannot determine number of channels. Assume single channel');
                nChannels = 1;
            end
            
            % determine channel type
            if isnumeric(stack)
                stackClass = class(stack);
                
            elseif iscell(stack)
                stackClass = class(stack{1});
            end
            
            % assign type
            if nChannels == 1
                channelPrefix = 'grayscale';
            elseif nChannels == 2
                channelPrefix = 'twoColor';
            elseif nChannels == 3
                channelPrefix = 'RGB';
            else
                error('more than 3 channels?');
            end
            
            if strcmp(stackClass, 'uint8')
                type = [channelPrefix '_8bit'];
            elseif strcmp(stackClass, 'uint16');
                type = [channelPrefix '_16bit'];
            else
                error('unknown stack type');
            end
        end 
    end
    
    %---------------------------------------------------------------------
    % public methods
    %---------------------------------------------------------------------
    
    methods
        
        % ------------------------------------------------------
        % constructor
        % ------------------------------------------------------
                
        function this = Stack(stack, varargin)
            % STACK Create a new Stack
            %
            % Stack(stack)
            % Stack(stack, resolution)
            % Stack(stack, resolution, description, type)
            %
            % stack:        4D stack with the 4 dimension color (or 3D for
            %               greyscale) or a cell array, e.g. {R,G,B}
            % type:         an entry in imageSpaces
            % resolution:   vector of resolutions
            % description:  description of data (date, labels, ..) 
            
            % for now we allow for making an empty stack to retrieve the
            % constant image spaces
            if isempty(stack)
                return;
            end

            % if stack type was passed use it, else detect it
            if nargin < 4
                [type nChannels] = this.detectStackType(stack);
            else
                type = varargin{3};
                if ~isfield(this.imageSpaces, type);
                    error(['unknown stack type ' type ' specified']);
                end
            end
            imageSpace = this.imageSpaces.(type);
            
            % set default channel to color mapping for color image
            if nChannels > 1
                this.channelColor = [1 2 3];
            end
            
            % if an array with multiple channels was passed, convert to
            % cell array containing the individual channels
            if isnumeric(stack)
                if imageSpace.dimension > 1 
                    grids = cell([1 imageSpace.dimension]);
                    for i = 1:imageSpace.dimension
                        grids{i} = stack(:,:,:,i);
                    end
                else
                    grids = {stack};
                end
            else
                grids = stack;
            end
            
            % next construct the image domain
            % we work with xyz, matlab yxz
            bdry = {[1 size(grids{1}, 2)] [1 size(grids{1}, 1)] [1 size(grids{1}, 3)]};
            imageDomain = diffgeometry.FiniteSet('sampleVolume', bdry, [1 1 1]);
            
            % now construct a map object holding the stack
            this.image = diffgeometry.Map(imageDomain, imageSpace, grids);
            
            % find intensity limits
            this.setIntensityLimits();
            
            % set metadata
            if nargin == 2
                this.resolution = varargin{1};
            elseif nargin >= 3
                this.resolution = varargin{1};
                this.description = varargin{2};
            else
                this.resolution = [1 1 1];
                disp(['No resolution provided for stackHolder,'...
                    'default [1 1 1]']);
            end
            
            % set projectionMask
            this.projectionMask = cell([3 1]);
        end 
        
        % ------------------------------------------------------
        % stack cropping, reduction, scaling
        % ------------------------------------------------------
        
        function redStack = reduceStack(this, ssfactor, average)
            % REDUCESTACK Generate down sampled Stack 
            %
            % redStack = reduceStack(ssfactor)
            % redStack = reduceStack(ssfactor, average)
            %
            % ssfactor: subsampling factor
            % average:  boolean, average over bins instead of just
            %           subsampling
            
            if nargin == 2
                average = false;
            end
            
            image   = this.image.apply();
            data = image{1}; 
            for i = 2 : length(image)
                data(:,:,:,i) = image{i};
            end
            
            if average
                
                ss = 4;
                ysize = size(data,1);
                xsize = size(data,2);
                zsize = size(data,3);
                
                newData = double(0*data(1:ss:ysize,1:ss:xsize,1:ss:zsize,:));
                weight = zeros(size(newData));
                for k = 1:ss
                    for j = 1:ss
                        for i = 1:ss
                            yidx = i:ss:ysize;
                            xidx = j:ss:xsize;
                            zidx = k:ss:zsize;
                            I = numel(yidx);
                            J = numel(xidx);
                            K = numel(zidx);
                            newData(1:I,1:J,1:K,:) = newData(1:I,1:J,1:K,:)...
                                        + double(data(yidx, xidx, zidx,:));
                            weight(1:I,1:J,1:K,:) = weight(1:I,1:J,1:K,:) + ones([I J K size(weight,4)]);
                        end
                    end
                end
                newData = newData./weight;
                newData = cast(newData, class(data));
            else
                newData = data(1:ssfactor:end,1:ssfactor:end,1:ssfactor:end,:);
            end
            
            redStack = surfaceDetection.Stack(newData,this.resolution/4, ...
                [this.description,' subsampled by factor ', num2str(ssfactor)]);
        end

        function bbox = getMaskBoundingBox(this)
            % CROPSTACK Generate a cropped stack 
            %
            % use generateProjectionMask or setProjectionMask 
            % to set the mask
            %
            % bbox = getMaskBoundingBox()
            %
            % bbox : bounding box [xmin xmax ymin ymax zmin zmax]
            
            stats = {};
            bb = {};
            for d = 1:3
                stats{d} = regionprops(this.projectionMask{d}, 'BoundingBox');
                bb{d} = floor([stats{d}.BoundingBox(1:2), stats{d}.BoundingBox(3:4) +  stats{d}.BoundingBox(1:2)]);
            end
            
            xmin = max(1, min(bb{1}(2), bb{3}(1)));
            ymin = max(1, min(bb{2}(2), bb{3}(2)));

            xmax = max(bb{1}(4), bb{3}(3));
            ymax = max(bb{2}(4), bb{3}(4));

            zmin = max(1, min(bb{2}(1), bb{1}(1)));
            zmax = max(1, min(bb{2}(3), bb{1}(3)));
            
            bbox = [xmin xmax ymin ymax zmin zmax];
        end
            
        function croppedStack = cropStack(this, bbox)
            % CROPSTACK Generate a cropped stack 
            %
            % use generateProjectionMask or setProjectionMask 
            % to set the cropping region
            %
            % croppedStack = cropStack()
            % croppedStack = cropStack(bbox)
            %
            % bbox : bounding box [xmin xmax ymin ymax zmin zmax]
            
            if nargin == 1
                bb = this.getMaskBoundingBox();
            else
                bb= bbox;
            end
            
            stackCell = this.image.apply;
            bb
            size(stackCell{1})
            cropped = stackCell{1}(bb(3):bb(4), bb(1):bb(2), bb(5):bb(6));

            for c = 2:length(stackCell)
                cropped(:,:,:,c) = stackCell{c}(bb(3):bb(4), bb(1):bb(2), bb(5):bb(6));
            end
            
            croppedStack = surfaceDetection.Stack(cropped,this.resolution, ...
                [this.description,' cropped']);
        end
        
        % ------------------------------------------------------
        % jitter correction: determine relative shift
        % ------------------------------------------------------
        
        function T = determineRelativeShift(this, stack_ss, ssfactor)
            % DETERMINERELATIVESHIFT Compute the shift (a translation) 
            % between self and a provided reference stack (subsampled by ssfactor).
            %
            % T = determineRelativeShift(stack_ss, ssfactor)
            %
            % stack_ss: subsampled reference stack
            % ssfactor: subsampling factor
            % T:        affine matrix containing relative translation
            
            % 1) subsample stack;
            self_ss = this.reduceStack(ssfactor);

            % 2) call xcorr3fft with self_ss and stack_ss
            debugMsg(2,'To Do: Implement fft for multi color data\n');
            shift = xcorr3fft(self_ss.image.apply{1},stack_ss.image.apply{1});  
            
            % 3) Return translation; 
            T = eye(4);
            T(1:3,end) = -shift;
        end
        
        % ------------------------------------------------------
        % metadata setters
        % ------------------------------------------------------
        
        function setResolution(this, resolution)
            % SETRESOLUTION Resolution setter
            %
            % setResolution(resolution)
            % 
            % resolution:   vector of resolutions
            
            this.resolution = resolution;
        end
       
        function setDescription(sH, description)
            % SETDESCRIPTION Set the string that describes the data set
            %
            % setDescription(description)
            
            sH.description = description;
        end
        
        function setChannelColor(this, channelColor)
            % SETCHANNELCOLOR Set mapping from channel to color
            %
            % setChannelColor(channelColor) 
            %
            % e.g channelColor = [3 1 2] for first channel B, second R, third G
            
            this.channelColor = channelColor;
        end
        
        % ------------------------------------------------------
        % set intensity limits
        % ------------------------------------------------------
        
        function setIntensityLimits(this, varargin)
            % SETINTENSITYLIMITS Compute intensity limits for image display. 
            %
            % setIntenstiyLimits()
            % setIntensityLimits(Ilim)
            %
            % Ilim: cell array of intensity limits for each channel
            %
            % If Ilim is not provided it is computed, based on the 
            % cumulative histogram of image intensity values 
            % compute the upper 95 and lower 10 percent limits of the histogram. 
            
            % TODO: deal with double, make percent cutoff an argument?
            
            debugMsg(1, 'Stack.setIntensityLimits()\n');
            
            % if limits are provided by the user
            if nargin == 2
                
                Ilim = varargin{1};
                if iscell(Ilim) && length(Ilim) == this.image.image.dimension
                    this.Ilim = Ilim;
                else
                    error('Ilim should be cell array whose size matches the number of channels');
                end
            
            % else detect the limits
            else
                data = this.image.apply();
                
                % for each channel
                for i = 1:this.image.image.dimension
                    
                    ssdata = data{i}(1:10:end, 1:10:end, 1:10:end);
                    
                    if isa(ssdata, 'uint16')
                        
                        debugMsg(2, 'uint16\n');
                        % avoid NaNs
                        n = histc(ssdata(~isnan(ssdata)),0:1:2^16);
                        
                    elseif isa(ssdata,'uint8')
                        
                        debugMsg(2, 'uint8\n');
                        % avoid NaNs
                        n = histc(ssdata(~isnan(ssdata)),0:1:2^8);
                    else
                        error(['data class ' class(data{i}) ' not uint16,'... 
                            'uint8, convert data to this type or set limits manually']);
                    end
                    
                    % check if cumulative sum is somwhere less than 10 percent of
                    % summed intensity. If not, set lower boundary to image min.. 
                    if ~isempty(find(cumsum(n)<.1*sum(n),1,'last'))
                        this.Ilim{i} = double([find(cumsum(n)<.1*sum(n),1,'last'),...
                            find(cumsum(n)>.99*sum(n),1,'first')]);
                    else
                        this.Ilim{i} = double([min(ssdata(~isnan(ssdata))),...
                            find(cumsum(n)>.99*sum(n),1,'first')]);
                    end
                end
            end
        end
        
        % ------------------------------------------------------
        % transformation
        % ------------------------------------------------------
        
        function newStack = rescaleToUnitAspect(this)
            % RESCALETOUNITASPECT Rescale the stack to unit aspect ratio
            %
            % newStack = rescaleToUnitAspect() 
            % 
            % assume x and y resolution the same, only scales z
            
            if this.aspect == 1
                newStack = this;
                return;
            end
            
            % prepare the perumutation of stacks, such that the low
            % resolution axis becomes the z axis;
            stacks     = this.image.apply;
            resolution = this.resolution;
            resolution = resolution([2,1,3]);
            [~,ii]     = sort(resolution);
            
            for i = 1:numel(stacks)
                
                debugMsg(2,['Channel ' num2str(i) '\n']);
                
% OLD WAY:
%                 curr = stacks{i};
%                 % permute;
%                 curr = permute(curr,ii);
%                 newnslices = round(size(curr,3)*this.aspect);
%                 scaled = zeros([size(curr,1) size(curr,2) newnslices], class(curr));
%                 for j=1:size(curr,1)
%                     debugMsg(1, '.');
%                     if rem(j,80) == 0
%                         debugMsg(1, '\n');
%                     end
%                     scaled(j,:,:) = imresize(squeeze(curr(j,:,:)),...
%                                                 [size(curr,2) newnslices]);
%                 end
%                 debugMsg(1,'\n');
%                 % permute back;
%                 scaled = ipermute(scaled,ii);

                tic
                
                % stack and stack differential
                stack = single(this.image.apply{i});
                dstack = stack(:,:,2:end) - stack(:,:,1:end-1);
                dstack = cat(3, dstack, 0*dstack(:,:,end));

                % refined zvalues
                nSlices = this.imageSize(3);
                zvalref = 1:1/this.aspect:nSlices;
                base = floor(zvalref);
                rel = zvalref - base;

                % rescale
                scaled = zeros([size(stack,1), size(stack,2), numel(zvalref)], class(this.image.apply{1}));
                for j = 1:numel(rel)
                    scaled(:,:,j) = stack(:,:,base(j)) + dstack(:,:,base(j))*rel(j);
                end
                
                stacks{i} = cast(scaled,class(this.image.apply{1}));
                
                toc
            end
            
            xres = this.resolution(1);
            newResolution = [xres xres xres];
            resolution = min(resolution)*ones(size(resolution));
            newStack = surfaceDetection.Stack(stacks, resolution);
            newStack.setChannelColor(this.channelColor);
            newStack.setDescription(this.description);
        end
        
        % ------------------------------------------------------
        % getters for dependent variables
        % ------------------------------------------------------

        function aspect = get.aspect(this)
            % get aspect ratio of stack resolution
            res    = this.resolution; 
            aspect = max(res)/min(res);
        end
        
        function imageSize = get.imageSize(this)
            % determine size of image grid
            imageSize = this.image.domain.gridSize;
        end
        
        % ------------------------------------------------------
        % getSlice
        % ------------------------------------------------------
        
        function slice = getSlice(this, dim, val)
            % GETSLICE Return a slice (color image) through the stack
            %
            % slice = getSlice(dim, val)
            %
            % dim:      dimension, accepts x,y,z or 1, 2, 3
            % val:      position along the dimension 
            % slice:    image array
            
            if ~isnumeric(dim)
                if strcmp(dim, 'x')
                    dim = 1;
                    dm = 2;
                elseif strcmp(dim, 'y')
                    dim = 2;
                    dm = 1;
                else
                    dim = 3;
                    dm = 3;
                end
            end
            
            bdry = this.image.domain.boundary;
            bdry{dim} = [val val];
            sliceDom = diffgeometry.FiniteSet('slice', bdry, [1 1 1]);
            
            if ~this.image.domain.hasSubset(sliceDom)
                error('slice index out of bounds');
            end
            
            sliceComponents = this.image.apply(sliceDom);
                        
            % projection mask along this direction (used in wing)
            
            mask = this.projectionMask{dm};
            if isempty(mask)
                mask = true(size(squeeze(sliceComponents{1})));
            end
            
            if length(sliceComponents) == 1
                slice = mat2gray(squeeze(sliceComponents{1}), this.Ilim{1}).*mask;
            else
                slice = zeros([size(squeeze(sliceComponents{1})) 3]);
                for i = 1:length(sliceComponents)
                    slice(:,:,this.channelColor(i)) = mat2gray(squeeze(sliceComponents{i}), this.Ilim{i}).*mask;
                end
            end
            
            % keep it left handed
            if dim == 1 
                slice = permute(slice, [2 1 3]);
            end
        end
        
        % ------------------------------------------------------
        % define mask
        % ------------------------------------------------------
        
        function generateProjectionMask(this, dimensions)
            % Define mask by polygon selection on max intensity projections 
            %
            % generateProjectionMask(dimensions)
            %
            % dimensions:   vector of directions, matlab style, e.g. [2 3]
            %               for [x z]
            %
            % Each slice along a direction is masked like the MIP, this is
            % coarse but much simpler than defining a truly 3D mask.
            
            if this.image.image.dimension > 1
                ncolors = 3;
            else
                ncolors = 1;
            end
            
            MIP = {};
            % imageSize is xyz, MIP also
            MIP{1} = zeros([this.imageSize(1) this.imageSize(3) ncolors]); 
            MIP{2} = zeros([this.imageSize(2) this.imageSize(3) ncolors]);
            MIP{3} = zeros([this.imageSize(2) this.imageSize(1) ncolors]);

            % create maximal intensity projections along all directions
            for d = dimensions
                for c = 1:this.image.image.dimension
                    data = this.image.apply;
                    MIP{d}(:,:, this.channelColor(c)) = mat2gray(max(data{c}, [], d));
                end
            end
            
            % create masks on the MIPs
            for d = 1:3
                if any(dimensions==d)
                    
                    xsize = size(MIP{d},2);
                    ysize = size(MIP{d},1);
                    
                    figure, imshow(MIP{d},[]);
                    rect = getrect;
                    xmin = max(1, rect(1));
                    ymin = max(1, rect(2));
                    xmax = min(xsize, rect(1)+rect(3));
                    ymax = min(ysize, rect(2)+rect(4));
                    [xmin xmax ymin ymax]
                    
                    this.projectionMask{d} = false([ysize xsize]);
                    this.projectionMask{d}(ymin:ymax,xmin:xmax) = true;
                    
                    % OLD:
%                     this.projectionMask{d} = roipoly(MIP{d});
                    close;
                else
                    this.projectionMask{d} = [];
                end
            end
        end
        
        function setProjectionMask(this, projectionMask)
            % Set projection masks
            %
            % setProcjectionMask(projectionMask)
            %
            % projectionMask:       cell array of 3 2D masks {x,y,z}
            % 
            % Each slice along a direction is masked with the 2D mask, 
            % this is coarse but much simpler than defining a truly 3D mask.
            
            if numel(projectionMask) ~= 3
                error('projectionMask should be cell array with 3 2D masks');
            end
             
            if ~isempty(projectionMask{2}) &&...
                ~all(size(projectionMask{2}) == [this.imageSize(1) this.imageSize(3)])
                error('projectionMask{2} has wrong dimensions');
            end

            if ~isempty(projectionMask{1}) &&...
                ~all(size(projectionMask{1}) == [this.imageSize(2) this.imageSize(3)])
                error('projectionMask{1} has wrong dimensions');
            end

            if ~isempty(projectionMask{3}) &&...
                ~all(size(projectionMask{3}) == [this.imageSize(2) this.imageSize(1)])
                error('projectionMask{3} has wrong dimensions');
            end
            
            this.projectionMask = projectionMask;
        end

        function mask = getMask(this)
            % GETMASK Create 3D mask from projection masks
            %
            % mask = getMask()
            %
            % mask:     3D mask (binary stack)
            %
            % Takes logical AND of the different projected masks

            mask = true([this.imageSize(2) this.imageSize(1) this.imageSize(3)]);
            
            % imageSize is xyz, but in projectionMask we're using matlab yxz indexing
            
            if ~isempty(this.projectionMask{3})
                mask = mask &...
                    repmat(this.projectionMask{3}, [1 1 this.imageSize(3)]);
            end
            
            if ~isempty(this.projectionMask{1})
                projMask = reshape(this.projectionMask{1},...
                                [this.imageSize(2) 1 this.imageSize(3)]);
                mask = mask & repmat(projMask, [1 this.imageSize(1) 1]);
            end
            
            if ~isempty(this.projectionMask{2})
                projMask = reshape(this.projectionMask{2},...
                                [1 this.imageSize(1) this.imageSize(3)]);
                mask = mask & repmat(projMask, [this.imageSize(2) 1 1]);
            end
        end
    end
    
end
