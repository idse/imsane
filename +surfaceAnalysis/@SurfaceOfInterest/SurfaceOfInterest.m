classdef SurfaceOfInterest < diffgeometry.Manifold2D
    % SurfaceOfInterest extends Manifold, adding methods that go beyond the
    % basic mathematical structure
   
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
    % dependent properties
    %---------------------------------------------------------------------
    
    properties (SetAccess = protected, Dependent = true)
        
        data               % data pullback
    end
    
    %---------------------------------------------------------------------
    % properties
    %---------------------------------------------------------------------
    
    properties (SetAccess = protected)
        
        nLayers            % number of layers
    end
    
    %---------------------------------------------------------------------
    % public methods
    %---------------------------------------------------------------------
    
    methods
        
        %------------------------------------------------------
        % constructor
        %------------------------------------------------------
        
        function this = SurfaceOfInterest(varargin)
            % SURFACEOFINTEREST Create a new SOI
            % 
            % SurfaceOfInterest(dataSpace, embeddingSpace, timePoints, dynamic)
            % SurfaceOfInterest(datadir)
            %
            % embeddingSpace:   FiniteSet
            % dataSpace:        space of image data
            % timePoints:       vector of integer timepoints for which fields
            %                   should be defined
            % dynamic:          if false, the geometry is time independent
            % datadir:          directory holding a saved manifold
            
            if nargin == 4
                dataSpace = varargin{1};
                varargin = varargin(2:4);
            end
            
            this = this@diffgeometry.Manifold2D(varargin{:});
            
            % there should always be a data field
            if nargin == 1 && all(this.getField('data')==0)
                
                dataDir = varargin{1};
                DOMnode = xmlread(fullfile(dataDir, 'SOI.xml'));
                node = xpathNode(DOMnode, 'SOI/fields/Field[name=''data'']/targetSpace/Set');
                dataSpace = diffgeometry.FiniteSet(node);
                
                this.createField('data', 'diffgeometry.TensorPatch', dataSpace, true);
                
            elseif nargin == 4 
                
                this.createField('data', 'diffgeometry.TensorPatch', dataSpace, true);
            end
            
            % determine nLayers indirectly
            fieldNames = this.getFieldNames;
            N = 1;
            for i = 1:numel(fieldNames)
                
                if strfind(fieldNames{i},'layer')
                    
                    lN = str2double(fieldNames{i}(13:end));
                    N = max(N,lN);
                end
            end
            this.nLayers = 2*N+1;
        end
        
        %------------------------------------------------------
        % pullbackStack
        %------------------------------------------------------
        
        function pullbackStack(this, stack, ROI, time, onionOpts)
            % PULLBACKSTACK Pull the stack data back to the surface
            %
            % pullbackStack(stack, ROI, time, onionOpts)
            %
            % stack:    Stack object containing the data
            % ROI:      RegionOfInterest object containing affine
            %           transformation from embedding coordinates to stack
            %           ROI can be left empty: ROI = []
            % time:     the time of the pullback, for storing in the SOI
            % onionOpts multi-layer pullback options, struct w fields
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

            % We have:
            % - a chart phi: M -> R^2
            % - an embedding X: M -> R^3
            % - a transformation T : R^3 -> R^3'
            % - an image Im: R^3' -> I 
            % then produce  Im \circ T \circ X \circ \phi^{-1}: R^2 -> I
            % T is the inverse of the transformation we used to
            % transform the raw stack to the stack we used for surface
            % detection
            
            % time index for data
            ti = this.tIdx(time);
            
            % initialize onion options
            if nargin == 4
                onionOpts = struct('nLayers', 1, 'layerDistance', 1,...
                            'sigma', 1, 'makeIP', false, 'IPonly', false,...
                            'zEvolve',false);
            else
                assert(isfield(onionOpts, 'nLayers') && rem(onionOpts.nLayers, 2),...
                'onionOptions should have field nLayers with odd integer value')
                if ~isfield(onionOpts, 'layerDistance')
                    onionOpts.layerDistance = 1;
                end
                if ~isfield(onionOpts, 'sigma'); onionOpts.sigma = 1; end
                if ~isfield(onionOpts, 'makeIP'); onionOpts.makeIP = 0; end
                if ~isfield(onionOpts, 'IPonly'); onionOpts.IPonly = 0; end
                if ~isfield(onionOpts, 'zEvolve'); onionOpts.zEvolve = 0; end
            end

            % store number of layers
            this.nLayers = onionOpts.nLayers;
            halfLayers = (this.nLayers - 1)/2;

            % create fields for onion layers
            for li = 1:halfLayers;

                name = 'data_layer_';
                PT = this.data.patchClass;
                TS = this.data.targetSpace;
                    
                if this.getField([name 'p' num2str(li)]) == 0 
                    this.createField([name 'p' num2str(li)], PT, TS, true);
                    this.createField([name 'm' num2str(li)], PT, TS, true);
                end
            end
                    
            % time index for geometry: can be different because for a
            % static geometry the index should be one for every time but
            % the data will still be time dependent
            if this.dynamic == false,   gti = 1;
            else                        gti = this.tIdx(time); end
            
            %-------------------------------------------------------
            % first transform the embedding if necessary: T \circ X
            %-------------------------------------------------------

            if ~isempty(ROI)
                % TODO: more precise input testing
                                
                % define map from the aligned space to the stack space
                imageVolume = stack.image.domain;
                alignedImageVolume = this.embedding.targetSpace;
                T = ROI.translation;
                R = ROI.rotation;
                %TJitter = alignment.TJitter;
                stackTransform = diffgeometry.AffineMap(alignedImageVolume,...
                                                        imageVolume, inv(T)*inv(R));
                
                % transformed embedding, has imageVolume as target
                embedding = diffgeometry.Field(this.embedding(gti).name,...
                    this.embedding(gti).patchClass, imageVolume,...
                    this.embedding(gti).topology);
                
                for i = 1:length(this.embedding(gti).patches)
                    newPatch = stackTransform.compose(this.embedding(gti).patches{i});
                    embedding.addPatch(newPatch);
                end
            else
                embedding = this.embedding(gti);
            end
            
            %-------------------------------------------------------
            % then pull back the image in each chart
            %-------------------------------------------------------
            
            for i=1:length(this.atlas(gti).charts)
                                
                curChart = this.atlas(gti).charts{i};
                curDom = curChart.domain;

                debugMsg(2, ['Pulling back stack to chart ' curChart.image.name '\n']);
                
                % (T \circ X) \circ \phi^{-1}
                EmbPatch = embedding.getPatch(curDom.name);
                chartInv = curChart.getInverse();
                chartEmb = EmbPatch.compose(chartInv);

                %---- mask interpolation outside data ---------------------

                % boundary curve: image of index set boundary
                ub = [curChart.apply{1}(1,:) curChart.apply{1}(:,end)'...
                    fliplr(curChart.apply{1}(end,:)) fliplr(curChart.apply{1}(:,1)')];
                vb = [curChart.apply{2}(1,:) curChart.apply{2}(:,end)'...
                    fliplr(curChart.apply{2}(end,:)) fliplr(curChart.apply{2}(:,1)')];
                
                umin = curChart.image.boundary{1}(1);
                vmin = curChart.image.boundary{2}(1);

                ub = (ub-umin)./curChart.image.stepSize(1);
                vb = (vb-vmin)./curChart.image.stepSize(2);

                usize = curChart.image.gridSize(1);
                vsize = curChart.image.gridSize(2);

                mask = poly2mask(double(ub), double(vb), vsize, usize);

                %---- create onion layers ---------------------
                
                projection = {};

                X = chartEmb.apply;
                dN = onionOpts.layerDistance;

                % compute normal
                [Nx,Ny,Nz] = GaussNormal(X{1},X{2},X{3},onionOpts.sigma);

                % normal displacement unit
                if onionOpts.zEvolve
                    debugMsg(2,'zEvolve\n');
                    dX = {0,0,dN};
                else
                    dX = {dN*Nx,dN*Ny,dN*Nz}; 
                end
                
                for ni = 1:numel(dX); dX{ni}(isnan(dX{ni})) = 0; end

                % now loop through the layers
                for li = 1:this.nLayers

                    idx = li - halfLayers - 1;

                    % normally evolved embedding
                    def = {X{1} + idx*dX{1}, X{2} + idx*dX{2}, X{3} + idx*dX{3}};
                    
                    %figure, imshow(def{3} > stack.image.domain.boundary{3}(2));
                    %[min(def{3}(:)) max(def{3}(:))]
                    
                    embLayer = diffgeometry.CoordinateMap(chartEmb.domain,...
                                            chartEmb.image, def);

                    % Im \circ (T \circ X \circ \phi^{-1})
                    projection = stack.image.compose(embLayer);
                    grids = projection.apply;
                    
                    % go through channels to mask
                    for ci = 1:numel(grids)
                        grids{ci}(~mask) = NaN;
                    end

                    % to deal with multi-channel pullbacks
                    if numel(grids) > 1
                        type = [0]; % non-tensorial index: multiple scalars
                    else
                        type = []; % single scalar
                    end

                    % convert index to field name
                    idx = li - halfLayers - 1;
                    if idx < 0
                        fieldName = ['data_layer_m' num2str(-idx)];
                    elseif idx > 0
                        fieldName = ['data_layer_p' num2str(idx)];
                    else
                        fieldName = 'data';
                    end

                    debugMsg(2, ['pulling back ' fieldName '\n']);

                    dataField = this.getField(fieldName);
                    dataPatch = dataField(ti).getPatch(curChart.domain.name);

                    if dataPatch ~= 0

                        debugMsg(2, ['field already defined on ' curChart.domain.name...
                            ' adding it as transformation\n']);

                        dataPatch.setTransform(curChart.image.name, grids, type);

                    else 
                        dataPatch = diffgeometry.TensorPatch(curChart.domain,...
                            projection.image, grids, type,...
                            curChart.image.name, this.atlas(gti));

                        dataField(ti).addPatch(dataPatch);
                    end  
                end
            end

            % make MIP/SIP if wanted
            if onionOpts.makeIP 
                
                if strcmp(onionOpts.makeIP, 'both')

                    opts = onionOpts;
                    opts.makeIP = 'MIP';
                    this.makeIP(time, opts);
                    opts.makeIP = 'SIP';
                    this.makeIP(time, opts);
                else
                    this.makeIP(time, onionOpts);
                end
                
                if onionOpts.IPonly
                    for li = 1:halfLayers;
                        name = 'data_layer_';
                        this.removeField([name 'p' num2str(li)]);
                        this.removeField([name 'm' num2str(li)]);
                    end 
                end
            end
        end
        
        %------------------------------------------------------
        % make Intensity Projection
        %------------------------------------------------------
        
        function makeIP(this, time, opts)
            % make intensity projection of multiplayer pullback
            %
            % makeIP(time, type)
            % makeIP(time, onionOpts)
            %
            % type:         'SIP' or 'MIP'
            % onionOpts:    see pullbackStack
            %
            % if only type provided:
            % assumes unit spacing for intensity correction

            if isstruct(opts)
                type = opts.makeIP;
                dN  = opts.layerDistance;
            else
                type = opts;
                dN = 1;
            end
            
            if nargin == 2 || ~(strcmp(type,'MIP') || strcmp(type,'SIP'))
                error('specify correct type for intensity projection: SIP, MIP');
            end
            
            assert(this.nLayers > 1, 'need multiple layers to make MIP or SIP');
            halfLayers = (this.nLayers - 1)/2;
            
            ti = this.tIdx(time);
            if this.dynamic == false,   gti = 1;
            else                        gti = ti; end
            
            % make field if necessary
            fieldName = ['data_' type];
            if this.getField(fieldName) == 0 
                PT = this.data.patchClass;
                TS = this.data.targetSpace;
                this.createField(fieldName, PT, TS, true);
            end
            
            IP = this.getField(fieldName);
            
            % now loop through charts and make mips
            for i=1:length(this.atlas(gti).charts)
                                
                curChart = this.atlas(gti).charts{i};
                domName = curChart.domain.name;
                chName = curChart.image.name;
                
                IPgrids = this.data(ti).getPatch(domName).getTransform(chName).apply;
                
                % for SIP we weight intensity with area change relative to
                % the middle surface to get the correct sum
                % this change is approximately -2*sqrt(detg).*H*dN
                % (linearizing in the normal displacement)
                % we could do it exactly
                
                if strcmp(type, 'SIP')
                    
                    X = this.embedding(gti).getPatch(domName).apply;
                    geom = GaussGeometry(X{1},X{2},X{3},0);
                    deltaA = -2*geom.H.*geom.dA;
                    dA = geom.dA;
                    
                    % a second approximation is to replace local area
                    % change by an overall area change factor because it is
                    % less sensitive to noise
                    % this could also be improved
                    
                    deltaAtot = sum(deltaA(~isnan(deltaA)));
                    dAtot = sum(dA(~isnan(deltaA)));
                    corrFactor = deltaAtot/dAtot;
                end
                
                % loop through the layers
                for li = 1:halfLayers

                    pField = this.getField(['data_layer_p' num2str(li)]);
                    mField = this.getField(['data_layer_m' num2str(li)]);
                    
                    pGrids = pField(ti).getPatch(domName).getTransform(chName).apply;
                    mGrids = mField(ti).getPatch(domName).getTransform(chName).apply;
                    
                    for ci = 1:numel(IPgrids)
                        
                        if strcmp(type,'MIP')
                            
                            IPgrids{ci} = max(IPgrids{ci}, pGrids{ci});
                            IPgrids{ci} = max(IPgrids{ci}, mGrids{ci});
                            
                        elseif strcmp(type, 'SIP')
                            
                            % correct for change in surface area by
                            % normally evolve (see above)
                            w = dN*li*corrFactor;
                            
                            IPgrids{ci} = double(IPgrids{ci})...
                                            + (1+w)*double(pGrids{ci})...
                                            + (1-w)*double(mGrids{ci});
                        end
                    end
                end
                
                IPpatch = IP(ti).getPatch(domName);
                dataPatch = this.data(ti).getPatch(domName);
                
                if IPpatch ~= 0

                    debugMsg(2, ['field already defined on ' curChart.domain.name...
                        ' adding it as transformation\n']);

                    IPpatch.setTransform(curChart.image.name, IPgrids, dataPatch.type);

                else 
                    IPpatch = diffgeometry.TensorPatch(dataPatch.domain,...
                        dataPatch.image, IPgrids, dataPatch.type,...
                        curChart.image.name, this.atlas(gti));

                    IP(ti).addPatch(IPpatch);
                end  
            end
        end
        
        %------------------------------------------------------
        % make stack of maps
        %------------------------------------------------------
        
        function registeredMultiLayer(this, saveDir, channel, chart, margin, tmax)
            % registeredMultiLayer 
            %
            % saveDir:          video written to multipage tif files in this dir
            % channel:          channel used for registration
            % chart:            name of chart
            % margin:           outer margin to ignore in cross correlation
            % tmax:             index of final time point
            
            if nargin == 5
                tmax = this.nTimePoints;
            end

            fieldName = 'data';
            shift = this.registerLayer(fieldName, channel, chart, margin);
            shift = cumsum(shift);

            for t = 1:tmax 
                this.multilayer2stack(this.timePoints(t), saveDir, shift(t,:));
            end
        end
        
        function shift = layer2video(this, saveDir, fieldName, channel, chart, varargin)
            % layer to tif video, optionally with registration
            %
            % shift = layer2video(saveDir, fieldName, channel, chart, registerFrames, margin)
            % shift = layer2video(saveDir, fieldName, channel, chart, registerFrames, margin, tmax)
            % 
            % saveDir:          video written to multipage tif files in this dir
            % fieldName:        field to make a video from
            % channel:  
            % chart:            name of chart
            % registerFrames:   use cross correlation to register frames
            % margin:           outer margin to ignore in cross correlation
            % tmax:             index of final time point
            %
            % shift: Tx2 array of relative shifts between consecutive
            % times (T relative to T-1)

            if nargin >= 7
                tmax = this.nTimePoints;
                registerFrames = varargin{1};
                margin = varargin{2};
            end
            if nargin == 8
                tmax = varargin{3};
            end
            
            L = margin;

            fname = fullfile(saveDir, [fieldName '_' chart '.tif']);

            if registerFrames
                shift = this.registerLayer(fieldName, channel, chart, margin);
                shift = cumsum(shift);
            end
            
            for t = 1:tmax

                field = this.getField(fieldName);
                patch = field(t).getPatch([chart '_index']);

                if t == 1

                    im2 = patch.getTransform(chart).apply{channel};
                    imwrite(im2, fname,'Compression','none');
                    ySize = size(im2,1);
                    xSize = size(im2,2);
                else

                    im1 = im2;
                    im2 = patch.getTransform(chart).apply{channel};

                    if registerFrames
                        
                        % transform image
                        xform = [1  0  0; 0  1  0; -shift(t,2) -shift(t,1)  1];
                        tform_translate = maketform('affine',xform);
                        cb = im2;
                        im2 = imtransform(cb, tform_translate,...
                                        'XData', [1 xSize],...
                                        'YData', [1 ySize]);
                    end
                    imwrite(im2, fname, 'writemode', 'append','Compression','none');
                end
            end
        end
        
        function shift = registerLayer(this, fieldName, channel, chart, margin)
            % find displacement between consecutive time points in a layer
            %
            % shift = registerLayer(this, fieldName, channel, chart, margin)
            % 
            % fieldName:        field to make a video from
            % channel:  
            % chart:            name of chart
            % margin:           outer margin to ignore in cross correlation
            %
            % shift: Tx2 array of relative shifts between consecutive
            % times (T relative to T-1)
            
            L = margin;
            shift = zeros([this.nTimePoints 2]);
            
            for t = 1:this.nTimePoints

                field = this.getField(fieldName);
                patch = field(t).getPatch([chart '_index']);

                if t == 1
                    im2 = patch.getTransform(chart).apply{channel};
                    ySize = size(im2,1);
                    xSize = size(im2,2);
                else
                    im1 = im2;
                    im2 = patch.getTransform(chart).apply{channel};

                    % determine relative shift
                    [shiftx, shifty] = xcorr2fft(im1(L:ySize-L,L:xSize-L),im2(L:ySize-L,L:xSize-L));
                    shift(t,:) = [shiftx shifty];
                end
            end
        end
        
        function stack = multilayer2stack(this, t, saveDir, shift)
            % multilayer map of region to image stack
            %
            % stack = multilayer2stack(t)
            % stack = multilayer2stack(t, saveDir)
            % stack = multilayer2stack(t, saveDir, shift)
            % 
            % stack:    cell array indexed by patch
            % t:        time index
            % shift:    [shiftx shifty] to apply a translatation
            % saveDir:  if directory is passed, stacks are written to
            %           multipage tif files
            
            stack = {};
            
            % make an onion stack
            halfLayers = (this.nLayers - 1)/2;
            
            % for each patch
            for pi = 1:numel(this.data(1).patches)

                pdef = this.data(1).patches{pi}.phi{1};
                nCh = numel(this.data(1).patches{pi}.phi);

                stack{pi} = zeros([size(pdef) this.nLayers nCh], class(pdef));
                        
                for li = 1:this.nLayers

                    idx = li - halfLayers - 1;

                    % convert index to field name
                    if idx < 0
                        fieldName = ['data_layer_m' num2str(-idx)];
                    elseif idx > 0
                        fieldName = ['data_layer_p' num2str(idx)];
                    else
                        fieldName = 'data';
                    end

                    data = this.getField(fieldName);
                    if numel(data)>1
                        data = data(t);
                    end
                    
                    % for each channel
                    for ci = 1:nCh
                        
                        pb = data.patches{pi}.apply{ci};
                
                        if nargin == 4
                            disp('transform image');
                            
                            xform = [1  0  0; 0  1  0; -shift(2) -shift(1)  1];
                            tform_translate = maketform('affine',xform);
                            cb = pb;
                            pb = imtransform(cb, tform_translate,...
                                            'XData', [1 size(cb,2)],...
                                            'YData', [1 size(cb,1)]);
                        end
                        
                        stack{pi}(:,:,li,ci) = pb;
                    
                        % save if saveDir was passed, write to file
                        if nargin >= 3

                            chrtName = this.data(1).patches{pi}.chartName;
                            fname = fullfile(saveDir, [chrtName '_stack_C'...
                                num2str(ci) '_T' num2str(t) '.tif']);
                            if li > 1
                                imwrite(pb, fname, 'writemode', 'append');
                            else
                                imwrite(pb, fname);
                            end
                        end
                    end
                end
            end
        end
        
        %------------------------------------------------------
        % measurements
        %------------------------------------------------------
        
        % general remark about stepSize in all measurements:
        %   We implemented a stepSize in each chart, that is not
        %   necessaryily isotropic. For example, on the frut fly embryo, we
        %   have dz = 1 and dphi = 1/R_max. This is not taken into account 
        %   at the level of the metric (which is natural, as the metric
        %   assumes chart values and not pixels. Since mostly the images and
        %   positions in the image will be considered as unit step pixels
        %   (e.g. from segmentation on the pullbacks), we assume the
        %   following convention: 
        %   Scalars, vectors and tensors are assumed to be fed in pixels
        %   and in the measurement functions changed to the according chart
        %   values;
        
        function lproper = properLength(this, time, curve, chartName)
            % calculate proper length of curve
            %
            % lproper = properLength(time, curve, chartName)
            %
            % time:         integer time point
            % curve:        Nx2 surface coordinates along curve
            % chartName:    chart in which curve is specified
            % 
            %   We assume the following convention: 
            %   Scalars, vectors and tensors are assumed to be fed in pixels
            %   and in the measurement functions changed to the according chart
            %   values;
            
                
                        
            ti = this.tIdx(time);
            if this.dynamic == false,   gti = 1;
            else                        gti = ti; end
            
            chart     = this.atlas(gti).getChart(chartName);
            stepSize  = chart.image.stepSize;
            u = curve;
            ushift = circshift(u, [1 0]);
            du = ushift - u;
            % if the curve is not closed, the first entry will be the
            % distance between beginning and end, Remove this contribution
            % by setting it to 0.
            du(1,:) = 0*du(1,:); 
            
            % correction by stepSize;
            du = [du(:,1)*stepSize(1),du(:,2)*stepSize(2)];
            
            if isempty(this.g(gti).patches)
                this.NCalcInducedMetric(time);
            end
            
            domName = this.atlas(gti).getChart(chartName).domain.name;
            g = this.g(gti).getPatch(domName).getTransform(chartName);
            %gcurve = g.apply({u(:,1),u(:,2)});
            metric    = g.apply();
            gcurve{1,1} = interp2(metric{1,1},u(:,1),u(:,2));
            gcurve{1,2} = interp2(metric{1,2},u(:,1),u(:,2));
            gcurve{2,1} = interp2(metric{2,1},u(:,1),u(:,2));
            gcurve{2,2} = interp2(metric{2,2},u(:,1),u(:,2));

            dlsq = du(:,1).*(gcurve{1,1}.*du(:,1) + gcurve{1,2}.*du(:,2))...
                  +du(:,2).*(gcurve{2,1}.*du(:,1) + gcurve{2,2}.*du(:,2));
            lproper = sum(sqrt(dlsq));
        end
        
        function Aproper = properArea(this, time, mask, chartName)
            % calculate proper area of region
            % 
            % Aproper = properArea(time, mask, chartName)
            %
            % time:         integer time point
            % curve:        binary mask with size matching chart
            % chartName:    chart in which mask is specified
            % 
            %   We assume the following convention: 
            %   Scalars, vectors and tensors are assumed to be fed in pixels
            %   and in the measurement functions changed to the according chart
            %   values;
            
            
                        
            ti = this.tIdx(time);
            if this.dynamic == false,   gti = 1;
            else                        gti = ti; end
            
            
            chart    = this.atlas(gti).getChart(chartName);
            stepSize = chart.image.stepSize;
            domName  = chart.domain.name;
            g        = this.g(gti).getPatch(domName).getTransform(chartName);
            
            %detg = g.determinant().apply{1};
            metric    = g.apply();
            
            detg      = metric{1,1}.*metric{2,2}-metric{1,2}.*metric{2,1};
            detg(isnan(detg)) = 0;
            temp      = sqrt(detg(mask==1));
            Aproper   = sum(temp(:))*stepSize(1)*stepSize(2);
        end
        
        function Aproper = properAngle(this, time,position, vector1, vector2, chartName)
            % calculate proper area of region
            % 
            % Aproper = properAngle(time, mask, chartName)
            %
            % time:         integer time point
            % position:     Nx2 surface coordinates
            % vector1 :     Nx2 vector
            % chartName:    chart in which mask is specified
            % 
            %   We assume the following convention: 
            %   Scalars, vectors and tensors are assumed to be fed in pixels
            %   and in the measurement functions changed to the according chart
            %   values;
            
            ti = this.tIdx(time);
            if this.dynamic == false,   gti = 1;
            else                        gti = ti; end
            
            u = position;
            
            
            chart     = this.atlas(gti).getChart(chartName);
            stepSize  = chart.image.stepSize;
            domName   = chart.domain.name;
            
            g         = this.g(gti).getPatch(domName).getTransform(chartName);
            
            %gpos      = g.apply({u(:,1),u(:,2)});
            metric    = g.apply();
            gpos{1,1} = metric{1,1}(u(:,1),u(:,2));
            gpos{1,2} = metric{1,2}(u(:,1),u(:,2));
            gpos{2,1} = metric{2,1}(u(:,1),u(:,2));
            gpos{2,2} = metric{2,2}(u(:,1),u(:,2));
            
            vector1   = [vector1(:,1)*stepSize(1),vector1(:,2)*stepSize(2)];
            vector2   = [vector2(:,1)*stepSize(1),vector2(:,2)*stepSize(2)];
            
            normV1sq  = ( gpos{1,1}.*vector1(:,1).^2+ 2*gpos{1,2}.*vector1(:,1).*vector1(:,2) ...
                            + gpos{2,2}.*vector1(:,2).^2);
            normV2sq  = ( gpos{1,1}.*vector2(:,1).^2+ 2*gpos{1,2}.*vector2(:,1).*vector2(:,2) ...
                            + gpos{2,2}.*vector2(:,2).^2);
            
            innerProd = gpos{1,1}.*vector1(:,1).*vector2(:,1) + gpos{1,2}.*vector1(:,1).*vector2(:,2) ...
                     +  gpos{2,1}.*vector1(:,2).*vector2(:,1) + gpos{2,2}.*vector1(:,2).*vector2(:,2);  
                        
            % for reasons of accuracy we compute the sqrt only once;
            Aproper = real(acos(sign(innerProd).*sqrt(innerProd.^2/normV1sq/normV2sq))); 
        end
        
        
        function Iproper = integrate(this, time, fieldname, mask, chartName)
            % calculate proper area of region
            % 
            % Aproper = properAngle(time, mask, chartName)
            %
            % time:         integer time point
            % position:     Nx2 surface coordinates
            % vector1 :     Nx2 vector
            % chartName:    chart in which mask is specified
            % 
            %   We assume the following convention: 
            %   Scalars, vectors and tensors are assumed to be fed in pixels
            %   and in the measurement functions changed to the according chart
            %   values;
            
            ti = this.tIdx(time);
            if this.dynamic == false,   gti = 1;
            else                        gti = ti; end
            
            
            chart     = this.atlas(gti).getChart(chartName);
            stepSize  = chart.image.stepSize;
            domName   = chart.domain.name;
            
            % get the metric;
            g         = this.getField('metric');
            if isempty(g(gti).patches)
                this.NCalcInducedMetric(gti);
            end
            g         = g(gti).getPatch(domName).getTransform(chartName);
     
            metric    = g.apply();
            
            detg      = sqrt(metric{1,1}.*metric{2,2}-metric{1,2}.*metric{2,1});
           
            % get the field in the desired representation;
            field    = this.getField(fieldname);
            field    = field(gti);            
            fieldrep = field.getPatch(domName).getTransform(chartName).apply();
            % the field representation may have multiple components; 
            Iproper = zeros(1,length(fieldrep));
            for i = 1 : length(fieldrep)
                
                maskedrep  = double(fieldrep{i}).*double(mask);
                maskeddetg = double(detg).*double(mask);
                maskedrep(isnan(maskedrep)) = 0;
                maskeddetg(isnan(maskeddetg)) = 0;
                temp       = maskedrep.*maskeddetg;
                
                Iproper(i) = sum(temp(:))*stepSize(1)*stepSize(2);%./sum(maskeddetg(:));
            end
        end
        
        
        %------------------------------------------------------
        % load segmentation
        %------------------------------------------------------
        
        function loadSegmentation(this, dataDir, segDir)
            % LOADSEGMENTATION load segmented data into SOI
            %
            % loadSegmentation(dataDir, segDir)
            %
            % dataDir:  directory containing saved version of this SOI
            % segDir:   directory containing the segmentation, this should
            %           be structured as a the data subdirectory in fields,
            %           i.e. have subdirectories for patches etc with the
            %           same file names as the raw counterparts
            
            % basically reuses code from Manifold2D.load but loads the
            % segmented files with the data field metadata
            
            DOMnode = xmlread(fullfile(dataDir, 'SOI.xml'));
            fieldsNode = xpathNode(DOMnode, 'SOI/fields');
            fieldNode = fieldsNode.getFirstChild;

            while ~isempty(fieldNode)

                if  strcmp(fieldNode.getNodeName, 'Field') &&...
                    strcmp(char(xpathText(fieldNode, 'name')), 'data')

                    t = str2num(fieldNode.getAttribute('time'));
                    patchClass = char(xpathText(fieldNode, 'patchClass'));
                    
                    node = xpathNode(fieldNode, 'targetSpace/Set');
                    tspaceClass= str2func(char(node.getAttribute('subclass')));
                    targetSpace = tspaceClass(node);

                    % create the segmentation field
                    fDynamic = 1; 
                    this.createField('segmentation', patchClass, targetSpace, fDynamic);
                    segField = this.getField('segmentation');

                    % the atlas that needs to be passed to patch
                    % constructor
                    if this.dynamic, gti = this.tIdx(t);
                    else gti = 1; end
                    patlas = this.atlas(gti);
                    
                    % patches: contain the actual field data
                    patchesNode = xpathNode(fieldNode, 'patches');
                    patchNode = patchesNode.getFirstChild;
                    
                    while ~isempty(patchNode)

                        if strcmp(patchNode.getNodeName, 'Map')

                            patchConstructor = str2func(patchClass);
                            domName = xpathText(patchNode, 'domain/Set/name');
                            pdir = fullfile(segDir, domName);
                            patch = patchConstructor(patchNode, pdir, patlas);
                            
                            % if the data cannot be found, phi is left empty
                            if ~isempty(patch.phi)
                                segField(this.tIdx(t)).addPatch(patch);
                            end
                        end
                        
                        patchNode = patchNode.getNextSibling;
                    end
                end
                fieldNode = fieldNode.getNextSibling;
            end
        end
        
        %------------------------------------------------------
        % getters for dependent properties
        %------------------------------------------------------
        
        function data = get.data(this)
            % return the metric tensor
            data = this.getField('data');
        end
    end
 
end