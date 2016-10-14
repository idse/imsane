classdef Atlas < handle_light
    % Collection of charts and the transition maps between them
    % such that together the charts cover the entire manifold (surface).
    %
    % Charts and transition maps are CoordinateMap objects.
       
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
        
        charts;         % cell array containing the charts
        topology;       % TopologicalSpace object containing sets on which charts are defined
        transitionMaps  % cell array containing transition maps between charts
    end

    %---------------------------------------------------------------------
    % public methods
    %---------------------------------------------------------------------
    
    methods
        
        %------------------------------------------------------
        % constructor
        %------------------------------------------------------
        
        function this = Atlas(topologicalSpace, varargin)
            % ATLAS Create an atlas for a topological space
            % 
            % Atlas(topologicalSpace)
            % Atlas(topologicalSpace, XMLnode, dataDir)
            
            % if nargin == 0, do nothing: this is some matlab bullshit
            % matlab requires calling the constructor without arguments
            % for object array creation

            if nargin == 0

                return;
                
            elseif nargin == 1

                this.charts = {};
                this.transitionMaps = {};
            
            elseif nargin == 3
                
                atlasNode = varargin{1};
                dataDir = varargin{2};
                
                %---------------------------------
                debugMsg(2, 'loading the charts\n')
                %---------------------------------
                
                chartsNode = xpathNode(atlasNode, 'charts');
                ncharts = 0;
                % chart grids are in atlas/charts/chartname/
                chartsdir = fullfile(dataDir, 'atlas', 'charts');
                chartNode = chartsNode.getFirstChild;
                
                while ~isempty(chartNode)
                    
                    if strcmp(chartNode.getNodeName, 'Map')
                        
                        ncharts = ncharts + 1;
                        
                        % chart name and dir
                        chartName = char(xpathText(chartNode, 'image/Set/name'));
                        debugMsg(2, ['\t' chartName '\n']);
                        chartdir = fullfile(chartsdir, chartName);
                        
                        % load subsampled chart
                        % subsampled chart has smaller phigrids than domain
                        % apply interpolates it
                        ssChart = diffgeometry.CoordinateMap(chartNode, chartdir, this);
                        % third argument (atlas) was [] for a while, no
                        % recollection of why
                        
                        % create fullsize chart
                        this.charts{ncharts} = diffgeometry.CoordinateMap(...
                                ssChart.domain, ssChart.image, ssChart.apply);
                    end
                    chartNode = chartNode.getNextSibling;
                end
                
                %------------------------------------------
                debugMsg(2, 'loading the transition maps\n')
                %------------------------------------------
                
                tmapsNode = xpathNode(atlasNode, 'transitionMaps');
                ntmaps = 0;
                tmapsdir = fullfile(dataDir, 'atlas', 'transitionMaps');
                tmNode = tmapsNode.getFirstChild;
                
                while ~isempty(tmNode)
                    
                    if strcmp(tmNode.getNodeName, 'Map')
                        
                        ntmaps = ntmaps + 1;
                        
                        % tmap name
                        tmapName = [char(xpathText(tmNode, 'domain/Set/name'))...
                            '_' char(xpathText(tmNode, 'image/Set/name'))];
                        tmdir = fullfile(tmapsdir, tmapName);
                        
                        % load subsampled tmap
                        ssTmap = diffgeometry.CoordinateMap(tmNode, tmdir, this);
                        % third argument (atlas) was [] for a while, no
                        % recollection of why
                        
                        % create fullsize tmap
                        this.transitionMaps{ntmaps} = diffgeometry.CoordinateMap(...
                                ssTmap.domain, ssTmap.image, ssTmap.apply);
                    end
                    tmNode = tmNode.getNextSibling;
                end
                
            else
                error('wrong number of arguments');
            end
            
            this.topology = topologicalSpace;
        end

        %------------------------------------------------------
        % add chart
        %------------------------------------------------------
        
        function addChart(this, chart)
            % ADDCHART Add a chart to the atlas
            % 
            % addChart(chart)
            % 
            % chart should be a CoordinateMap whose domain is in the topology
            
             % check if domain is in topology
            if this.topology.getSetID(chart.domain.name) == 0
                error('add set to topology before defining a field over it');
            end
            
            % check that image is not in topology
            if this.topology.getSetID(chart.image.name) > 0
                error('image is in topology, did you mean to add a transition map?');
            end
                
            chartIndex = this.getChartIndex(chart.image.name);

            % check if chart is already in atlas
            if chartIndex == 0
                n = length(this.charts);
                this.charts{n+1} = chart;
            else
                disp(['chart ' chart.image.name ' is already in the atlas, overwriting']);
                this.charts{chartIndex} = chart;
            end

        end

        %------------------------------------------------------
        % add transition map
        %------------------------------------------------------
        
        function addTransitionMap(this, tMap)
            % ADDTRANSITIONMAP Add transition map between sets,
            % by adding a map between coordinates defined on sets (chart image)
            %
            % addTransitionMap(tMap)
            % 
            % tMap should be a CoordinateMap whose domain and image match
            % the images of two distinct charts in the atlas
            
            domID =  this.getChartIndex(tMap.domain.name);
            imID = this.getChartIndex(tMap.image.name);
            
            % both charts have to exist
            if imID == 0 || domID == 0 
                error(['one of the coordinates (chart image) you are '...
                    'trying to map is not in the atlas']);
            end
            
            % check if transitionmap already exists
            tMapIndex = getTMapIndex(this, tMap.domain.name, tMap.image.name);
            if tMapIndex == 0
                n = length(this.transitionMaps);
                this.transitionMaps{n+1} = tMap;
            else
                disp(['transition map : ' tMap.domain.name...
                    ' -> ' tMap.image.name ' already exists, overwriting']);
                this.transitionMaps{tMapIndex} = tMap;
            end
        end
        
        %------------------------------------------------------
        % getters
        %------------------------------------------------------
        
        function chart = getChart(this, imageSetName)
            % GETCHART Get chart from the atlas using image name, 
            % which is a unique identifier
            %
            % chart = getChart(imageSetName)
            
            chartIndex = this.getChartIndex(imageSetName);
            if chartIndex == 0
                error(['Atlas.getChart: chart ' imageSetName ' is not in atlas']);
            else
               chart = this.charts{chartIndex}; 
            end
        end
        
        function chartIndex = getChartIndex(this, imageSetName)
            % GETCHARTINDEX Get chart index (in atlas.charts) using image name. 
            % If the chart doesn't exist return 0.
            %
            % chartIndex = getChartIndex(imageSetName)

            if isempty(this.charts)
                debugMsg(2, 'Atlas.getChartIndex: atlas contains no charts\n');
                chartIndex = 0;
            else
                for chartIndex=1:length(this.charts)
                    if strcmp(this.charts{chartIndex}.image.name, imageSetName)
                        return;
                    end
                end
                chartIndex = 0;
            end
        end
        
        function tMap = getTransitionMap(this, domName, imName)
            % GETTRANSITIONMAP Get transition map using image and domain names. 
            %
            % tMap = getTransitionMap(domName, imName)
            
            tMapIndex = this.getTMapIndex(domName, imName);
            if tMapIndex == 0
                error('transition map not in atlas');
            else
                tMap = this.transitionMaps{tMapIndex};
            end
        end
        
        function tMapIndex = getTMapIndex(this, domName, imName)
            % GETTMAPINDEX Get transition map index using image and domain names. 
            % If the transition map doesn't exist return 0.
            %
            % tMapIndex = getTMapIndex(domName, imName)
            
            for tMapIndex=1:length(this.transitionMaps)
                
                currDomName = this.transitionMaps{tMapIndex}.domain.name;
                currImName = this.transitionMaps{tMapIndex}.image.name;
                
                if strcmp(currDomName, domName) && strcmp(currImName, imName)
                    return;
                end
            end
            tMapIndex = 0;
        end
        
        %------------------------------------------------------
        % XML & save
        %------------------------------------------------------
        
        function objectNode = save(this, docNode, options, varargin)
            % SAVE: save atlas data to tiff and return metadata as XML node
            % metadata is returned by XMLnode and linked into SOI.xml
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
            % charts and transition maps are subsampled when saved to save
            % disk space, because they vary very smoothly little
            % information is lost
            
            if nargin == 4
                fnamePostfix = ['_' varargin{1}];
            else
                fnamePostfix = [];
            end
            
            assert(isfield(options, 'dir'), 'options.dir needs to be specified');
            dir = options.dir;
            
            objectNode = docNode.createElement('Atlas');    
            
            %----------------------
            % save charts
            %----------------------
            
            elem_node = docNode.createElement('charts');
            
            for ci = 1:length(this.charts)
                
                % save subsampled version of chart
                ssChart = this.charts{ci}.subsample(2);
                
                suboptions = options;
                suboptions.dir = fullfile(dir, 'charts', this.charts{ci}.image.name);
                subelem_node = ssChart.save(docNode, suboptions, fnamePostfix);
                elem_node.appendChild(subelem_node);
            end
            objectNode.appendChild(elem_node);
            
            %----------------------
            % save transition maps
            %----------------------
            
            elem_node = docNode.createElement('transitionMaps');
            
            for tmi = 1:length(this.transitionMaps)
                
                % save subsampled tMap
                ssTmap = this.transitionMaps{tmi}.subsample(2);
                
                tmname = [ssTmap.domain.name '_' ssTmap.image.name];
                suboptions = options;
                suboptions.dir = fullfile(dir, 'transitionMaps', tmname);
                subelem_node = ssTmap.save(docNode, suboptions, fnamePostfix);
                elem_node.appendChild(subelem_node);
            end
            objectNode.appendChild(elem_node);
        end
    end
end
