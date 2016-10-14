classdef Field < handle_light
    % Map defined over the whole manifold, which is a collection of Map 
    % objects defined on different sets in the topology, here referred to
    % as patches.
       
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
    
    properties (SetAccess = private)
        
        name;           % field name
        topology;       % topological space for surface of interest
        patches;        % maps over different sets in the topology
        targetSpace;    % space into which the field maps (image of patches)
        patchClass;     % type of Map for the patches
    end
    
    %---------------------------------------------------------------------
    % public methods
    %---------------------------------------------------------------------
    
    methods
        
        %------------------------------------------------------
        % constructor
        %------------------------------------------------------
        
        function this = Field(name, patchClass, targetSpace, topology)
            % FIELD Create new field
            %
            % Field(name, patchClass, targetSpace, topology)
            %
            % For example, for the embedding the patchClass is
            % diffgeometry.CoordinateMap and the targetspace is a 3D set
            % representing the embeddingSpace.
            % Topology is the topologicalSpace over which the field is
            % defined.
            
            % if nargin == 0, do nothing: this is some matlab bullshit
            % matlab requires calling the constructor without arguments
            % for object array creation
            if nargin == 0
                return;
            end
                
            % check input types (incomplete)
            if ~isa(topology, 'diffgeometry.TopologicalSpace')
                error('4th argument should be topological space');
            end
            
            this.name = name;
            this.topology = topology;
            this.patches = {};
            this.patchClass = patchClass;
            
            if isa(targetSpace, 'diffgeometry.Set')
                this.targetSpace = targetSpace;
            else
                error('target space should be a Set');
            end
        end
        
        %------------------------------------------------------
        % add map, get map
        %------------------------------------------------------
        
        function addPatch(this, patch)
            % ADDPATCH Add patch to the field
            %
            % addPatch(patch)
            %
            % patch is a Map whose domain has to be in topology
            
            % check if domain is in topology
            if this.topology.getSetID(patch.domain.name) == 0
                error(['add domain ' patch.domain.name...
                    ' to topology before defining a field over it']);
            end
            if ~isa(patch, this.patchClass)
                error(['patch is not of class ' this.patchClass ', but ' class(patch)]);
            end
            if ~strcmp(patch.image.name, this.targetSpace.name)
                error(['target space of patch "' patch.image.name...
                    '" does not match field target space "'...
                    this.targetSpace.name '"'] );
            end
            
            % whereas a chart has a unique image, an patch should have
            % a unique domain (only one patch per set in the topology)
            patchIndex = this.getPatchIndex(patch.domain.name);

            % check if chart is already in atlas
            if patchIndex == 0
                n = length(this.patches);
                this.patches{n+1} = patch;
            else
                disp('there is already an patch for this set, overwriting');
                this.patches{patchIndex} = patch;
            end
        end
        
        function patch = getPatch(this, domainName)
            % GETPATCH Get patch using the name domain name
            %
            % patch = getPatch(domainName)
            %
            % patch is a Map object
            
            patchIndex = this.getPatchIndex(domainName);
            if patchIndex ~= 0
                patch = this.patches{patchIndex}; 
            else
                patch = 0;
            end
        end
        
        function patchIndex = getPatchIndex(this, domainName)
            % GETPATCHINDEX Get patch index (in Field.patches) using domain name
            % If the map doesn't exist return 0.
            %
            % patchIndex = getPatchIndex(domainName)
            
            for patchIndex=1:length(this.patches)
                if strcmp(this.patches{patchIndex}.domain.name, domainName)
                    return;
                end
            end
            patchIndex = 0;
        end
        
        function checkDefinition(this)
            % CHECKDEFINITION List on which sets a patch has been defined.
            %
            % checkDefinition()
            
            for i=1:length(this.topology.sets)
                curr = this.topology.sets{i}.name;
                if this.getPatchIndex(curr) == 0
                    disp(['field not defined on ' curr]);
                else
                    disp(['field defined on ' curr]);
                end
            end
        end
        
        function setName(this, name)
            % SETNAME Set field name
            %
            % setName(name)
            
            disp(['renaming ' this.name ' to ' name]);
            this.name = name;
        end
        
        %------------------------------------------------------
        % save
        %------------------------------------------------------
        
        function objectNode = save(this, docNode, options, fnamePostfix)
            % SAVE: saves Field data to tiff and returns metadata as XML node
            % metadata is returned by XMLnode and linked into SOI.xml
            %
            % objectNode = save(docNode, options)
            % objectNode = save(docNode, options, filenamePostfix)
            %
            % docNode:      document node to which objectNode should belong
            % objectNode:   node representing the object
            % options:              struct with the following fields
            %   - dir               directory to save to
            %   - imwriteOptions    options to pass to imwrite
            % fnamePostfix: file name postfix, e.g. timestamp
            
            objectNode = docNode.createElement('Field');

            % field name
            elem_node = docNode.createElement('name');
            text_node = docNode.createTextNode(this.name);
            elem_node.appendChild(text_node);
            objectNode.appendChild(elem_node);

            % patchClass
            elem_node = docNode.createElement('patchClass');
            text_node = docNode.createTextNode(this.patchClass);
            elem_node.appendChild(text_node);
            objectNode.appendChild(elem_node);

            % targetSpace
            elem_node = docNode.createElement('targetSpace');
            setNode = this.targetSpace.save(docNode);
            elem_node.appendChild(setNode);
            objectNode.appendChild(elem_node);
            
            % patches
            % POINTERS to prevent redundancy?
            % map: if docNode/SOI/atlas has a set, just give name
            % similarly: field targetSpace
            % set attribute pointer='atlas/...'
            elem_node = docNode.createElement('patches');
            
            for pi = 1:length(this.patches)
                
                % create directory
                % yes: within timeloop, in case we don't have all
                % patches at all times
                pdir = fullfile(options.dir, this.patches{pi}.domain.name);
                if ~exist(pdir, 'dir'), mkdir(pdir); end 
                
                popts = options;
                popts.dir = pdir;

                % optional subsampling
                if isfield(options, 'subsample')
                    patch = this.patches{pi}.subsample(options.subsample);
                    popts = rmfield(popts, 'subsample');
                else
                    patch = this.patches{pi};
                end
                
                patchNode = patch.save(docNode, popts, fnamePostfix);
                elem_node.appendChild(patchNode);
            end
            objectNode.appendChild(elem_node);
        end
    end
end
