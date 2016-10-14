classdef (Abstract) surfaceDetector < handle_light
    %   Extract a surface point cloud from a data stack
    %
    %   This is an abstract handle class that needs to be implemented for
    %   specific cases such as the wing disc.
    
    %---------------------------------------------------------------------
    % license
    %---------------------------------------------------------------------

    % Copyright 2014 Idse Heemskerk and Sebastian Streichan
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
        pointCloud          % point cloud object
        options             % detector options
    end
    
    properties (Constant, Abstract)
        defaultOptions     % default detector options
    end
    
    %---------------------------------------------------------------------
    % public methods
    %---------------------------------------------------------------------    

    methods (Abstract)
        
        % Abtract surface detection method
        % seed has a form that depends on the implementation
        detectSurface(this, stack, seed); 
    end
    
    methods
    
        function this = surfaceDetector()
            % SURFACEDETECTOR Create a new SurfaceDetector
            
            this.options = this.defaultOptions;
        end

        % ------------------------------------------------------
        % set options
        % ------------------------------------------------------

        function setOptions(this, options)
            % SETOPTIONS Update provided detector options
            %
            % For options see detectSurface help for specific subclass
            
            fields = fieldnames(this.defaultOptions);

            for i = 1:numel(fields)
                if isfield(options, fields{i})
                    this.options.(fields{i}) = options.(fields{i});
                else
                    disp([fields(i) ' unchanged']);
                end
            end
        end
        
        % ------------------------------------------------------
        % Write Quality inspection to disc.
        % ------------------------------------------------------
        
        function storeQualityInspection(this, options, stack)
            % Write quality inspection output to file.
            % 
            % storeQualityInspection(options, stack)
            %
            % inspect quality and store to file. Pass on options, that
            % specify the range of planes in the corresponding dimension to
            % be inspected. Optionally specify visibility of figures,
            % export and the inspection options. 
            %
            % options:              struct with fields:
            % - inspectOptions:     options for inspectQuality
            % - range:              range of planes along dimension
            % - outName:            filename of stored figure
            % - export:             boolean write to file
            % - closeFig:           boolean close figure after
            %
            % See also inspectQuality
            
            % parse the options
            % Options contains all writeOptions as fields and a structure
            % with options for inspectQuality.
            if ~isfield(options,'range')   
                options.range = 1:stack.zpSize;
            end
            if ~isfield(options,'outName');
                options.outName = 'qualityPointCloud.tif';
            end
            if ~isfield(options,'export');  
                options.export = 'true'; 
            end
            if ~isfield(options,'closeFig')
                options.closeFig = 'true';
            end
            if ~isfield(options,'inspectOptions')
                options.inspectOptions = struct();
            end
            
            
            if strcmp(options.closeFig,'true')
                h = figure;
            end
            % loop through all optionally set values (default is full size)
            for value = options.range
                options.inspectOptions.value = value;
                this.inspectQuality(options.inspectOptions, stack)
                
                if strcmp(options.export,'true') 
                    if (value-options.range(1))>0
                        I = getframe(h);
                        imwrite(I.cdata, options.outName,'tiff','Compression','none','WriteMode','append');
                    else
                        I = getframe(h);
                        imwrite(I.cdata, options.outName,'tiff','Compression','none');
                    end
                end
            end

            if strcmp(options.closeFig,'true')
                close(h);
            end
        end
        
        % ------------------------------------------------------
        % inspect quality of single slize and display image.
        % ------------------------------------------------------
        
        function inspectQuality(this, options, stack)
            % Inspect quality of detection by overlaying a point cloud
            % through a slice in the data
            % 
            % inspectQuality(options, stack)
            %
            % inspect quality of single slice in dimension specified by options 
            % and display image.
            %
            % options:  dimension: x,y,z 
            %           value
            %           pointCloud: color (matlab convention, g, w etc)

            if ~isfield(options,'dimension'), options.dimension = 'z'; end
            
            if options.dimension == 'x' 
                ind = [1 2 3];
            elseif options.dimension == 'y'
                ind = [2 3 1];
            else
                ind = [3 1 2];
            end

            % display image    
            imshow(stack.getSlice(options.dimension, options.value)); 
            
            % overlay point cloud. 
            if isfield(options, 'pointCloud')
                
                lspec = ['*' options.pointCloud];
                PC = this.pointCloud.unalignedPoints;
                
                PCplane = PC(round(PC(:,ind(1))) == options.value, :);
                hold on 
                plot(PCplane(:,ind(2)), PCplane(:,ind(3)), lspec, 'MarkerSize', 2)    
                hold off
            end
        end
        
    end
   
end
