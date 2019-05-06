classdef ilastikSegmentor < surfaceAnalysis.twoDSegmentor
    %ilastikSegmentor Segmentation of pullbacks using ilastik. 
    
    properties 
        classifierPath % Specifies path to ilastik classifier for data type
        referenceHist  % Reference histogram
        prediction     % Output prediction from Ilastik
        predictions    % Output structure with batch processed Ilastik predictions
        cellLabel      % Segmented cell interiors
        cellLabels     % Output structure with batch segmented cell interiors
        patches        % image patches 
        batchTimes     % time points for batch processing;
        lattices       % Output structure for batch processed Lattice tranforms of each patch;
    end

    % ------------------------------------------------------
    % public methods
	% ------------------------------------------------------
    
    methods
        
        % ------------------------------------------------------
        % constructor
        % ------------------------------------------------------
        
        function this = ilastikSegmentor(image, classifierPath)
            % Constructor
            % ilastikSegmentor(chartName, SOI, transformName, time)
            
            this      = this@surfaceAnalysis.twoDSegmentor(image);
            this.classifierPath = classifierPath;
        end

    end
    
   
    methods 
        
        % ------------------------------------------------------
        %  preProcessing for Ilastik segmentaion -> 8 bit tranform
        % ------------------------------------------------------
        
        function preProcess(this,h5namedummy)
            % preProcess 
            
            
            % THIS WILL CHANGE! In general, we will do image histogram equalisation. For now
            % we only map the 16 bit image to 8 bit. 
          
            
            %tempImage = this.image.apply();
            tempImage = this.image;
            % image histogram and cumulative sum. ; 
            for k = 1 : length(tempImage)
            
                n      = histc(tempImage{k}(:),0:2^16);
                nn     = cumsum(n);

                amin   = find(nn<=min(nn)*.1,1,'last');
                if isempty(amin)
                    amin = find(nn>min(nn),1,'first');
                end
                amax   = find(nn>=max(nn)*.999,1);
                scaledImage  = uint8(255*mat2gray(tempImage{k},[amin amax]));

                dset = zeros([1,size(scaledImage,2),size(scaledImage,1),1,1],'uint8');

                dset(1,:,:,1,1) = scaledImage';

                if exist(sprintf(h5namedummy,k),'file')~=0
                    delete(sprintf(h5namedummy,k))
                end
                h5create(sprintf(h5namedummy,k),'/volume/data/',[1 size(scaledImage,2) size(scaledImage,1) 1 1],'Datatype','uint8')

                h5write(sprintf(h5namedummy,k),'/volume/data/',dset);
            
            end
        end
        
        % ------------------------------------------------------
        %  batch preProcessing for Ilastik segmentaion -> 8 bit tranform
        % ------------------------------------------------------
        
        function batchPreProcess(this,time,xp,chartName,transformName)
            % batchPreProcess images for Ilastik segmentation
            % batchPreProcess(time,xp,chartName,transformName);
            
            for c = 1 : xp.expMeta.channelsUsed
                for t = time 

                    h5name = sprintf(xp.fileMeta.filenameFormat,t);
                    h5name = sprintf([h5name(1:end-4),'_t%u','_c%u_',transformName,'.h5'],t,c);
                    h5name = fullfile(this.classifierPath,h5name);


                    %this.image = xp.SOI.data(t).getPatch([chartName,'_index']).getTransform(transformName);                
                    %read the file from disc instead;
                    fName = fullfile(xp.fileMeta.directory,'fields','data',...
                        transformName,[xp.fileMeta.filenameFormat(1:end-4),'_t',...
                        num2str(t),'_c',num2str(c),'_',transformName,'.tif']);
                    if exist(fName,'file')>0
                        this.image = cell(1,1);
                        this.image{1} = imread(fName);
                        
                        this.preProcess(h5name);
                    else
                        this.image = xp.SOI.getField('data').getPatch([chartName,'_index']).getTransform(transformName).apply();
                        this.preProcess(h5name);
                    end
                end
            end
        end
        
        % ------------------------------------------------------
        %  import Ilastik prediction
        % ------------------------------------------------------
        
        function importPrediction(this,h5name,foregroundChannel)
            % importPrediction Import Ilastik prediction;
            % importPrediction(h5name)
            % importPrediction(h5name, foregroundChannel)
            
            if nargin < 3
                foregroundChannel = 2;
            end
            
            
            if exist(h5name,'file')~=2
                error(['Cant locate processed Ilastik file. Expecting it in ', h5name])
            end
            file = h5read(h5name,'/exported_data');
            
            
            this.prediction = squeeze(file(foregroundChannel,:,:))';
            
%             tsName    = 'prediction';
%             bd        = [0 1];
%             stepSize  = .01;
%             image     = diffgeometry.FiniteSet(tsName,{bd},stepSize);
%             
%             fieldName  = 'prediction';
%             patchClass =  'diffgeometry.TensorPatch';
%             Field(name, patchClass, targetSpace, topology)
            
        end
        
        
        % ------------------------------------------------------
        %  batch import Ilastik prediction
        % ------------------------------------------------------
        
        function batchImportPrediction(this,h5nameDummy,time,xp,foregroundChannel)
            % batchImportPrediction Import Ilastik prediction;
            % importPrediction(h5nameDummy)
            % importPrediction(h5nameDummy, time, foregroundChannel)
            
            this.batchTimes = time;
            
            if nargin < 5
                foregroundChannel = 2;
            end
            
            %npatches = length(xp.SOI.data(this.batchTimes(1)).patches);
            npatches = length(xp.SOI.data(1).patches);
            
            
            this.predictions = struct([]);
            
            for np = 1 : npatches
            
                %nreps = length(xp.SOI.data(this.batchTimes(1)).patches{np}.availableReps);
                nreps = length(xp.SOI.data(1).patches{np}.availableReps);
                
                for nr = 1 : nreps
                    % define the h5name; for this loop over all the chart and
                    % transform Names; 
                    %chartName       = xp.SOI.data(this.batchTimes(1)).patches{np}.chartName;
                    %transformName   = xp.SOI.data(this.batchTimes(1)).patches{np}.availableReps{nr};
                    chartName       = xp.SOI.data(1).patches{np}.chartName;
                    transformName   = xp.SOI.data(1).patches{np}.availableReps{nr};
                    this.predictions(1).(transformName) = '';
                    
                    
                    

                    % now loop over time; 
                    for t = 1:length(this.batchTimes)
                        
                        h5name = [h5nameDummy(1:end-4),'_t%u_c1_',transformName,'_Probabilities.h5'];
                        h5name = sprintf(h5name,t);
                        % import this prediction;
                        if exist(sprintf(h5name,this.batchTimes(t)),'file')
                            this.importPrediction(sprintf(h5name,time(t)),foregroundChannel);
                            this.predictions(t).(transformName) = this.prediction;
                        end
                    end
                end
            end
            
        end
        
        % ------------------------------------------------------
        %  segment Membranes
        % ------------------------------------------------------
        
        function segmentMembranes(this,prob_th)
            
            if nargin <2
                prob_th = .5;
            end
            
            seg   = this.prediction>prob_th;            
            label = bwlabeln(~seg);
            % remove cells (4-connected components) of less than 40 pixels 
            % extra white pixel outside to avoid disconnected regions on the edge
            cells = true([size(label,1)+2,size(label,2)+2]);
            cells(2:end-1,2:end-1) = bwareaopen(label, 100, 4);
            % skeletonize space between cells to get the segmented membrane
            membrane = bwmorph(~cells, 'skel', 'inf');
            % remove membrane that doesn't end in a vertex
            membrane = bwmorph(membrane, 'shrink', 'inf');
            membrane = bwmorph(membrane, 'clean');

            % force pixels on boundary to vanish. This results in a single
            % connected region of boundary (including the corners, where no cells
            % are)
            membrane(1,:)   = 0;
            membrane(end,:) = 0;
            membrane(:,end) = 0;
            membrane(:,1)   = 0;
            % remove the extra boundary columns and rows. 
            membrane       = membrane(2:end-1,2:end-1);
            cells          = ~membrane;
            this.cellLabel = bwlabeln(cells,4);
        end
       
        % ------------------------------------------------------
        %  batch segment Membranes
        % ------------------------------------------------------
        
        function batchSegmentMembranes(this,prob_th)
            
            if nargin <2
                prob_th = .5;
            end
            
            this.cellLabels = struct([]);
            
            fNames  = fieldnames(this.predictions);
            nFields = length(fNames);
            
            for f = 1 : nFields
                
                for t = 1 : length(this.batchTimes)
                    
                    this.prediction = this.predictions(t).(fNames{f});
                    if ~isempty(this.prediction)
                        this.segmentMembranes(prob_th);
                        this.cellLabels(t).(fNames{f}) = this.cellLabel;
                    end
                end
            end
        end
        
        % ------------------------------------------------------
        %  turn membrane segmentation into lattice structure
        % ------------------------------------------------------
        
        function g = membraneToLattice(this)
            
            % get vertices
            vertices = bwmorph(~(this.cellLabel>0), 'branchpoints');
            %imshow(cat(3, mat2gray(membrane), mat2gray(vertices), zeros(size(membrane))));

            % make a list of the positions of all the vertices
            [vertexY, vertexX] = ind2sub(size(vertices), find(vertices));

            % make a cell array tracking which vertices belong to which cells
            ncells =  max(this.cellLabel(:));
            cellVertices = cell([1 ncells]);

            % Extract cell properties, such as centre of mass from image;
            % Remember that corner images 
            sp = regionprops(this.cellLabel);
           % Centroids = cat(1,sp.Centroid);
            Areas     = cat(1,sp.Area);

           %in  = find(Areas<max(Areas));

            % for each vertex, add index to cells that are 8-connected to it
            for i=1:length(vertexX)
                belongsto = unique(this.cellLabel(vertexY(i)-1:vertexY(i)+1,vertexX(i)-1:vertexX(i)+1));
                for j = 2:length(belongsto)
                    cellVertices{belongsto(j)} = [cellVertices{belongsto(j)}, i];
                end
            end

            vertexPosition = [vertexX, vertexY];
            g = GLattConversion(cellVertices(Areas<max(Areas)), vertexPosition);
            
        end
        
        % ------------------------------------------------------
        %  batch process membrane segmentation into lattice structure
        % ------------------------------------------------------
        
        function batchMembraneToLattice(this)
            
            this.lattices = struct([]);
            
            fNames  = fieldnames(this.cellLabels);
            nFields = length(fNames);
            
            for f = 1 : nFields
                
                for t = 1 : length(this.batchTimes)
                    
                    this.cellLabel = this.cellLabels(t).(fNames{f});
                    if ~isempty(this.prediction)
                        
                        g = this.membraneToLattice();
                        this.lattices(t).(fNames{f}) = g;
                    end
                end
            end
        end
    end
end
