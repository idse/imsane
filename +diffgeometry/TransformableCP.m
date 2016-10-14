classdef TransformableCP < diffgeometry.ComponentPatch
    % ComponentPatch whose components obey a coordinate transformation law
    % In particular tensors and the Christoffel connection, but not a
    % Jacobian.
   
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

        chart;      % chart in which the components are expressed
    end
    
    properties (Access = protected)
        
        atlas;      % atlas handle to check that chart actually exists etc
        % coordReps - storing computed transformations as componentPatches
        % 2D cell array with rows of the form {chartname, ComponentPatch} 
        coordReps;  
    end
    
    properties (Dependent = true)

        chartName;      % name of component chart
        availableReps;  % representations (transformations) available
    end
    
    %---------------------------------------------------------------------
    % protected methods
    %---------------------------------------------------------------------
    
    methods (Access = protected)
        
        function repIndex = getRepIndex(this, chartName)
            % GETREPINDEX get row index in coordReps of representation 
            
            for repIndex  = 1: size(this.coordReps,1)
                if strcmp(this.coordReps{repIndex, 1}, chartName)
                    return
                end
            end
           repIndex = 0;
        end
    end
    
    %---------------------------------------------------------------------
    % abstract methods
    %---------------------------------------------------------------------
    
    methods (Abstract = true)
        
        % CALCULATETRANSFORM Calculate components in a new chart
        % using transformation law that depends on type of object
        %
        % abstract method needs to be implemented for different
        % transformation laws
        newTensor = calculateTransform(chartName)
    end
    
    %---------------------------------------------------------------------
    % public methods
    %---------------------------------------------------------------------

    methods
        
        % ------------------------------------------------------
        % constructor
        % ------------------------------------------------------
        
        function this = TransformableCP(varargin)
            % TRANSFORMABLECP Create transformable component patch
            %
            % TransformableCP(domain, image, def, type, chartName, atlas)
            % TransformableCP(domain, image, def, TransformableCP)
            % TransformableCP(XMLnode, dataDir, atlas)
            %
            % The second is for subclasses to have a contructor with the
            % same syntax, copying additional properties from the argument
            % map. Because Map.compose will call the constructor of any
            % subclass with the same arguments.
            % The third is for loading saved maps and is called from the
            % manifold2D.load
            %
            % See also Map ComponentPatch
            
            if nargin==4
                
                assert(isa(varargin{4},'diffgeometry.TransformableCP'),...
                    'when passing four arguments, fourth should be TransformableCP');
                
                chartName = varargin{4}.chartName;
                atlas = varargin{4}.atlas;
                superClassArgin = varargin;
                
            elseif nargin == 6
                
                chartName = varargin{5};
                atlas = varargin{6};
                superClassArgin = varargin(1:4);
                
            elseif nargin == 3
                
                XMLnode = varargin{1};
                dataDir = varargin{2};
                atlas = varargin{3};
                
                chartName = xpathText(XMLnode, 'chartName');
                debugMsg(2, ['\t' chartName '\n']);
                
                % for a TransformableCP the defining Map is actually stored
                % in chart_index/chart instead of chart_index so modify dir
                superClassArgin = varargin;
                superClassArgin{2} = fullfile(dataDir, chartName);
                
            else
                error('TransformableCP: wrong number of arguments for constructor');
            end
            
            % call parent constructor
            % ComponentPatch(XMLnode, dataDir, atlas)
            this = this@diffgeometry.ComponentPatch(superClassArgin{:});
            this.atlas = atlas;
            
            % check the domain of the chart
            if this.atlas.getChartIndex(chartName) == 0
                error('chart is not in in atlas!');
            else
                argChart = this.atlas.getChart(chartName);
            end

            if ~( strcmp(argChart.domain.name, this.domain.name) ||...
                    strcmp(argChart.image.name, this.domain.name) )
                error(['transformableCP is defined on ' this.domain.name...
                    ', its component chart on ' argChart.domain.name...
                    ' - this is incompatible']);
            else
                this.chart = argChart;
            end
            
            % self reference so that defining coordinate system shows up as
            % representation, but leave object entry empty because
            % form reason matlab seems to copy by value here (i.e. not
            % treat it as a pointer
            this.coordReps = {this.chart.image.name, []};    
            
            % now load the "non-fundamental" reps from disk
            if nargin == 3

                repsNode = xpathNode(XMLnode, 'coordReps');
                repNode = repsNode.getFirstChild;
                
                while ~isempty(repNode)
                    
                    if strcmp(repNode.getNodeName, 'Map')

                        % rep name
                        repName = char(xpathText(repNode, 'domain/Set/name'));
                        debugMsg(2, ['\t' repName '\n']);
                        
                        % load 
                        repDir = fullfile(dataDir, repName);
                        rep = diffgeometry.ComponentPatch(repNode, repDir, atlas);
                        this.setTransform(repName, rep.phi, rep.type);
                    end
                    repNode = repNode.getNextSibling;
                end
            end
            
        end        
        
        % ------------------------------------------------------
        % transformation
        % ------------------------------------------------------
        
        function newCP = getTransform(this, chartName)
            % GETTRANSFORM Get ComponentPatch object in the specified chart
            %
            % newCP = getTransform(chartName)
            %
            % The domain of the patch is an 'index' domain, i.e.
            % the patch is a map T: M -> TensorSpace,
            % but the transformations are from specific coordinate systems,
            % the represent T\circ \phi^-1 : R^2 -> TensorSpace'

            chart = this.atlas.getChart(chartName);
            if chart == 0 
                error(['TensorPatch.transform: chart ' chartName 'does not exist']);
            end

            repIdx = this.getRepIndex(chartName);
            
            % else: see if it's already there and return it
            if  repIdx > 1
                newCP = this.coordReps{repIdx, 2};
            
            % if the repIndex==1, the patch is defined on an index set but
            % we want to return it as defined in an actual chart, this is
            % just swapping the domain in this case
            % YES, THIS IS A LITTLE AWKWARD
            elseif repIdx == 1
                
                debugMsg(1, 'getTranform returning defining rep\n');
                preCP = this;
                newCP = diffgeometry.ComponentPatch(chart.image,...
                                    preCP.image, preCP.apply, preCP.type);
            
            % if not, actually transform
            else
                this.calculateTransform(chartName);
            end
        end
        
        function setTransform(this, chartName, def, type)
            % SETTRANSFORM Define a transformation of a tensor
            %
            % setTransform(chartName, def, type)
            %
            % Useful because we often want to directly calculate tensors in
            % different charts, rather than compute them from existing
            % representations. For example the pullback in proper cylinder
            % coordinates is made directly, not computed from the cylinder
            % pullback.
            %
            % The domain of the patch is an 'index' domain, i.e.
            % the patch is a map T: M -> TensorSpace,
            % but the transformations are from specific coordinate systems,
            % the represent T\circ \phi^-1 : R^2 -> TensorSpace'
            
            chart = this.atlas.getChart(chartName);
            assert(chart ~= 0, ['ChristoffelPatch.transform: chart '...
                                               chartName 'does not exist']);
            
            repIdx = this.getRepIndex(chartName);
            
            % make a new ComponentPatch object
            % a field is defined from the manifold , i.e. an index set, to
            % the field space but other representations of a tensorPatch
            % are defined from some chart, so we have to deal with the
            % patch domain differently
            if repIdx==1,   patchDom = chart.domain;
            else            patchDom = chart.image;     end

            % image is the same after transformation, e.g. RGB_8bit for a scalar
            newCP = diffgeometry.ComponentPatch(patchDom,...
                                                this.image, def, type);

            % thus CP does not contain information on the chart it is
            % defined in and we need to store this here
            % add it to the reps
            if repIdx == 0
                nreps = length(this.coordReps);
                this.coordReps{nreps+1, 1} = chartName;
                this.coordReps{nreps+1, 2} = newCP;
                
            elseif repIdx == 1
                disp('overwriting defining rep');
                this.type = newCP.type;
                this.phi = newCP.phi;
                
            elseif repIdx >1
                disp('transformation generated before, overwriting');
                chartIndex = this.getRepIndex(chartName);
                this.coordReps{chartIndex, 2} = newCP;
            end
        end
        
        % ------------------------------------------------------
        % dependent props
        % ------------------------------------------------------
        
        function chartName = get.chartName(this)
            chartName = this.chart.image.name;
        end
        
        function availableReps = get.availableReps(this)
            availableReps = {};
            for i = 1:size(this.coordReps,1)
                availableReps = [availableReps, this.coordReps{i,1}];
            end
        end
        
        %------------------------------------------------------
        % XML & save
        %------------------------------------------------------
        
        function objectNode = save(this, docNode, options, varargin)
            % SAVE: saves TransformableCP data to tiff and returns metadata
            % as XML node
            % overloads ComponentPatch.save to deal with multiple representations 
            %
            % objectNode = save(docNode, options)
            % objectNode = save(docNode, options, filenamePostfix)
            %
            % docNode:              document node to which objectNode should belong
            % objectNode:           node representing the object
            % options:              struct with the following fields
            %   - dir               directory to save to
            %   - imwriteOptions    options to pass to imwrite
            %
            % tiff filename is of form 'cmp_i_j_filenamePostfix.tif'
            % filenamePostfix can be used to append time for example
            %
            % See also Map.save
            
            assert(isfield(options, 'dir'), 'options.dir needs to be specified');
            dir = options.dir;
            
            % the defining rep can be saved by the superclass method
            CPopts = options;
            CPopts.dir = fullfile(dir, this.availableReps{1});
            objectNode = save@diffgeometry.ComponentPatch(this, docNode, CPopts, varargin{:});
            
            % store the chartName
            elem_node = docNode.createElement('chartName');
            text_node = docNode.createTextNode(num2str(this.chart.image.name));
            elem_node.appendChild(text_node);
            objectNode.appendChild(elem_node);
            
            % for each rep that is not the defining rep, store
            reps = this.availableReps;      
            repsNode = docNode.createElement('coordReps');
            
            for ri = 2: length(reps)

                rdir = fullfile(dir, reps{ri});
                if ~exist(rdir, 'dir'), mkdir(rdir); end

                % all the components to loop over
                repOpts = options;
                repOpts.dir = rdir;
                repNode = this.getTransform(reps{ri}).save(docNode, repOpts, varargin{:});
                repsNode.appendChild(repNode);
            end
            objectNode.appendChild(repsNode);
        end
    end
end