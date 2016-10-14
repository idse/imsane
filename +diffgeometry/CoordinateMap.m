classdef CoordinateMap < diffgeometry.Map
    % Map whose image represent coordinates
    % for example: embedding, chart, transition map.
    %
    % For coordinate maps we can calculate Jacobians,
    % bijective coordinate maps have an inverse.
       
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
    
    properties (Access = private)
        
        inverse % stores the inverse map so it doesn't have to be recomputed
    end
    
    %---------------------------------------------------------------------
    % public methods
    %---------------------------------------------------------------------
    
    methods
        
        %------------------------------------------------------
        % constructor
        %------------------------------------------------------
        
        function this = CoordinateMap(varargin)
            % CoordinateMap Create a map between two sets of coordinates
            %
            % CoordinateMap(domain, image, definition)
            % CoordinateMap(domain, image, definition, map)
            % CoordinateMap(XMLnode, dataDir, atlas)
            %
            % The second is for subclasses to have a contructor with the
            % same syntax, copying additional properties from the argument
            % map. Because Map.compose will call the constructor of any
            % subclass with the same arguments.
            % The third is for loading saved maps and is called from the
            % manifold2D.load
            %
            % See also Map.Map
            
            % just call parent constructor
            this = this@diffgeometry.Map(varargin{:});
            
            this.inverse = [];
        end

        %------------------------------------------------------
        % generate inverse: may not work!
        %------------------------------------------------------
        
        function mInv = getInverse(this)
            % getInverse Construct inverse map by interpolation of domain over image
            %
            % mInv = getInverse()
            
            if this.domain.dimension ~= this.image.dimension
                error(['for inverse to exist domain dimension should ',...
                        'equal image dimension']);
            end
            
            % if it was computed before, we don't have to recompute
            if ~isempty(this.inverse)
                mInv = this.inverse;
            else
                
            debugMsg(2, 'CoordinateMap: computing inverse');
            
            % make grids of everything 
            imageGrids = this.image.makeGrids();
            domainGrids = this.domain.makeGrids();
            phiGrids = this.apply();
            
            % initialize inverse grid and interpolant cell arrays
            phiInvGrids = cell(this.phiSize);
            interpolant = cell(this.phiSize);
            
            % inverse is interpolating the domain values over the image
            % values
            % interp2/3 require meshgrid like input so TriScatteredInterp
            % needs to be used here
            % TriScatteredInterp needs double inputs...
            
            if this.domain.dimension == 2
                
                for i = 1:2
                    interpolant{i} = TriScatteredInterp(...
                        double(phiGrids{1}(:)), double(phiGrids{2}(:)),...
                                                double(domainGrids{i}(:)));
                    phiInvGrids{i} = interpolant{i}(double(imageGrids{1}),...
                                                    double(imageGrids{2}));
                end

            elseif this.domain.dimension == 3
                
                warning('3D inverse: untested code');
                for i = 1:3
                    interpolant{i} = TriScatteredInterp(...
                        double(phiGrids{1}(:)), double(phiGrids{2}(:)),...
                        double(phiGrids{2}(:)), double(domainGrids{i}(:)));
                    phiInvGrids{i} = interpolant{i}(double(imageGrids{1}),...
                        double(imageGrids{2}), double(imageGrids{3}));
                end
            else
                error('dimension is not 2 or 3?');
            end

            mInv = diffgeometry.CoordinateMap(this.image, this.domain, phiInvGrids);
            this.inverse = mInv;
            
            end
        end
        
        function setInverse(this, mInv)
            % setInverse Set the inverse to prevent calculating it numerically
            %
            % setInverse(mInv)
            %
            % this is useful for analytic coordinate maps where we can
            % calculate the inverse by hand
            
            this.inverse = mInv;
        end
        
        %------------------------------------------------------
        % Jacobian
        %------------------------------------------------------
        
        function J = getJacobian(this, mp)
            % getJacobian Jacobian of this map over other map on same domain
            % 
            % J = getJacobian(mp)
            %
            % If m, mp are charts over the same set S with dimensions
            % labeled by i, j, we calculate J^i_j = \partial m^i \partial mp^j.
            
            if ~strcmp(this.domain.name, mp.domain.name)
                error('the maps are not defined on the same domain');
            end
            
            debugMsg(1, 'CoordinateMap.getJacobian()\n');
            
            % J^i_j = \partial u^i/ \partial u'^j
            % here we use the chainrule
            % J^i_j = sum_k (\partial u^i \partial x^k)*(\partial x^k / 
            % \partial u'^i)
            % where x^k indicates the indices and u^i the charts 
            
            % u_k ~ this.phi{k}, du_{kl} ~ \partial_l u_k
            % derivatives of this map with respect to matrix indices
            duGrids = cell([this.image.dimension this.domain.dimension]);
            uGrids = this.apply;
            
            for i = 1: this.image.dimension
                [duGrids{i,1}, duGrids{i,2}] = gradient(double(uGrids{i}));    
            end
            
            % derivatives of other map wrt matrix indices
            dupGrids = cell([mp.image.dimension mp.domain.dimension]);
            upGrids = mp.apply;
               
            for i = 1: mp.image.dimension
                [dupGrids{i,1}, dupGrids{i,2}] = gradient(double(upGrids{i}));    
            end
            
            % now make matrix objects so we can take advantage of the
            % methods defined for those
            
            debugMsg(1, 'TODO: Jacobian index type determination is primitive:\n');
            % if this is the embedding, the embedding index should not be
            % contractable with some tangent space index so we set its
            % index type to 0
            if this.image.dimension == 2
                type = [1 -1];
            elseif this.image.dimension == 3
                type = [0 -1];
            end
            
            du = diffgeometry.ComponentPatch(this.domain,...
                diffgeometry.Set('du', numel(duGrids)), duGrids, type);
          
            dup = diffgeometry.SquareMatPatch(this.domain,...
                diffgeometry.Set('dup', numel(dupGrids)), dupGrids, [1 -1]);
            
            dupInv = dup.inverse();
            
            J = du.contract(2,1, dupInv);
        end
        
    end

end