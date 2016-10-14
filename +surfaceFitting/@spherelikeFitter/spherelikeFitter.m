classdef spherelikeFitter < surfaceFitting.surfaceFitter
    %spherelikeFitter Fit spherelike point cloud using polynomials
    
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
    
    properties (SetAccess = protected)
        
        % structure collecting information about the charts.
        % the rows have the structure 
        % {name, description, stepsize, intersections, fundamental, desired}
        % name: cell array containing strings describing the chart names
        % description: Verbal description of charts.
        % stepsize: 2d vector of with stepsizes in the charts,
        % sequenced as in name.
        % intersections: lists indices of charts with overlapping domain
        % fundamental: boolean, meaning whether the chart defines a set
        % in the topology
        % desired: boolean indicating whether to produce it
        charts = struct(...
            'name', {   'cylinder1',        'cylinder2',...
                        'polarLowerZ',      'polarUpperZ',...
                        'cylinder1_proper', 'cylinder2_proper', ...
                        'anteriorEquidistant',      'posteriorEquidistant'},...
            'description', {'Cylinder coordinates z, phi',...
                            'Cylinder shifted by pi',...
                            'Polar projection in lower z planes',...
                            'Polar projection in upper z planes',...
                            'Proper cylinder coordinates',...
                            'Shifted proper cylinder coordinates',...
                            'Proper coordinates in lower z planes',...
                            'Proper coordinates in upper z planes',},...
            'stepSize', {[1 0.1], [1 0.1], [1 1], [1 1], [1 1], [1 1], [1 1], [1 1]},...
            'intersections', {[2 3 4], [1 3 4], [1 2], [1 2], [], [], [], []},...
            'fundamental', {1, 1, 1, 1, 0, 0, 1 ,1},...
            'desired', {1, 1, 0, 0, 1, 1, 1, 1});
    end

    % ------------------------------------------------------
    % public methods
	% ------------------------------------------------------
    
    methods
        
        % ------------------------------------------------------
        % constructor
        % ------------------------------------------------------
        
        function this = spherelikeFitter()
            % Constructor
            %
            % spherelikeFitter()
            
            this      = this@surfaceFitting.surfaceFitter();
            
            % initialize fitOptions
            this.fitOptions.R = 4;
            this.fitOptions.X = 4;
            this.fitOptions.Y = 4;
            this.fitOptions.e = 1;
            this.fitOptions.phase = 1;
            this.fitOptions.shift = 0;
            this.fitOptions.path = [];
            
            % initialize fittedParam
            % we have the option between using fourier and polynomial fit
            % of the phase. We may make this an option in the future.

            this.fittedParam = struct('e',[],'R',[],'pR',[],'SR',[],...
                'muR',[],'pX',[],'SX',[],'muX',[],'pY',[],'SY',[],'muY',[],...
                'pphase',[],'Sphase',[],'muphase',[],'pe', [],'Se',[],'mue',[],...
                'zmin',[],'zmax',[],'zRmaxUpper',[],'zRmaxLower',[],...
                'medianphase',[]);
        end

        % ------------------------------------------------------
        % fit spherelike surface to point cloud
        % ------------------------------------------------------
        
        function fitSurface(this, pointCloud)
            % Fit a surface to the point cloud.
            % 
            % fitSurface(pointCloud)
            %
            % pointCloud is a pointCloud object, see surfaceDetection for
            % details. 
            % call parameter fit to fit parametrization. 
            % set fitDomain based on pointCloud.ROI
            % ROI wasn't stored as a finiteSet to begin with because then
            % its orientation would be a map and we didnt want to
            % complicate things by using the diffgeometry machinery before
            % getting to the level of fitter
            bdry = {round(pointCloud.ROI.xpRange),...
                round(pointCloud.ROI.ypRange), round(pointCloud.ROI.zpRange)};
            this.fitDomain = diffgeometry.FiniteSet('alignedEmbeddingSpace',...
                                                            bdry, [1 1 1]);

            % fit surface using parametrization
            this.parameterFit(pointCloud, this.fitOptions);

            % adjust stepsize in charts based on fit
            z     = unique(pointCloud.points(:,3));
            rsq   = polyval(this.fittedParam.pR, z, this.fittedParam.SR,...
                            this.fittedParam.muR);

            c1step = this.getChartResolution('cylinder1');
            dz = c1step(1);
            dphi = dz/sqrt(max(rsq));
            stepSize = [dz dphi];
            this.setChartResolution('cylinder1', stepSize);
            this.setChartResolution('cylinder2', stepSize);
            
            % generate cylinder1 embedding to inspect the fit and normally
            % evolve
            this.fittedPoints = this.generateEmbedding('cylinder1').apply;
        end
     
        % ------------------------------------------------------
        %  Evolve surface along surface normal.
        % ------------------------------------------------------
        
        function normallyEvolve(this, shift, varargin)
            %   Evolve surface by shift in the direction of the surface
            %   normal
            % 
            %   normallyEvolve(shift)
            %   normallyEvolve(shift, ss) 
            %
            %   Evolve a surface represented by a point cloud in the 
            %   direction of its normal. 
            %   shift is a multiplicative factor to the surface normal. 
            %   Option ss is a subsampling factor of
            %   the fitted pointcloud for speed 

            debugMsg(1, 'fitter.normallyEvolve()\n');
            
            X = this.fittedPoints;
            
            % subsampling
            if nargin == 3
                ss = varargin{1};
                this.fitOptions.normEvolveSS = ss;
                debugMsg(2, ['subsampling fitted surface by factor ' num2str(ss) '\n']);
                for i = 1:3
                    X{i} = X{i}(1:ss:end, 1:ss:end);
                end
            end

            % update shift
            this.fitOptions.shift = this.fitOptions.shift + shift;
            
            % compute surface normal
            [Nx,Ny,Nz] = surfnorm(X{1},X{2},X{3});
            % point cloud in 3 x np matrix
            x  = [X{1}(:),X{2}(:),X{3}(:)];
            % normal in 3 x np matrix
            n  = [Nx(:),Ny(:),Nz(:)];
            % shift point cloud by normal.
            xs = x  + shift*n;
            xs = xs + (rand(size(xs,1),size(xs,2))-.5);
            % remove nans;
            for k = 1 : 3
                n(isnan(x(:,k)),:)   = 0;
                x(isnan(x(:,k)),:)   = 0;
                xs(isnan(xs(:,k)),:) = 0;
            end
            % No points outside the shifted domain should remain.
            xs(xs(:,3)> (max(x(:,3))-shift),:) = 0;
            xs(xs(:,3)< (min(x(:,3))+shift),:) = 0;
            xs = xs(sum(xs,2)~=0,:);
            
            outCloud = surfaceDetection.PointCloud(xs);
            
            xpRange = [min(xs(:,1)) max(xs(:,1))];
            ypRange = [min(xs(:,2)) max(xs(:,2))];
            zpRange = [min(xs(:,3)) max(xs(:,3))];
            
            outCloud.ROI.setRanges(xpRange, ypRange, zpRange);
            
            this.fitSurface(outCloud);
        end
    end
    
   
    methods (Access = private)
        
        % ------------------------------------------------------
        %  Parametric fit
        % ------------------------------------------------------
        
        function parameterFit(this, pointCloud, fitOptions)
            % Parameter fit to point cloud slice by slice.
            %
            % parameterFit(pointCloud, fitOptions)
            %
            % Fit parameters of parametrization to pointCloud slice by
            % slice along symmetry axis. This assumes that parameters are
            % a slowly varying function of symmetry axis and fits
            % polynomial of degree specified in fitOptions. We fit ellipses
            % (using ellipseFit.m) end store the following parameters: 
            % centre, radius, eccentricity and phase.
            
            % obtain fit parameters from point cloud.
            x = pointCloud.points(:,1);
            y = pointCloud.points(:,2);
            z = pointCloud.points(:,3);
            % loop through z, fit ellipses to points and extract the radii 
            % R, Rmax, the eccentricity, the phase (ellipse angle) and the centre xc,yc;

            debugMsg(2, 'Begin fitting spherelike surface\n'); 
            
            rz = round(z);

            % zmin can be negative, since we align the point cloud such that the major axis 
            % coincides with the third axis for final detection.  		
            zmin = min(rz(:));
            zmax = max(rz(:));
            % structure to keep fit information for the slices;
            points = struct('xc',[],'yc',[],'R',[],'Rmax',[],'e',[],'phase',[]);
	
            zind = zmin:zmax;
         
            for t = 1:length(zind)
                
                rzeqt = (rz==zind(t));
               
                % require a minimal number of data points for fit. 
                if sum(rzeqt) > 20
                    
                    xi = x(rzeqt);
                    yi = y(rzeqt); 
                    % fit ellipse to pointcloud xi,yi
                    a = ellipseFit(xi,yi);
                    % There is a very rare bug in ellipse fit (line 669, numel(A) =
                    % 0), then, we return NaN; 
                    % store results
                    if ~isnan(a(1))
                        points(t).xc    = a(1);
                        points(t).yc    = a(2);
                        points(t).R     = min(a(3:4));
                        points(t).Rmax  = max(a(3:4));
                        points(t).e     = sqrt(1-( min(a(3:4))/max(a(3:4)))^2);
                        points(t).phase = a(5);
                    else
                        points(t).xc    = 0;
                        points(t).yc    = 0;
                        points(t).R     = 0;
                        points(t).Rmax  = 0;
                        points(t).e     = 0;
                        points(t).phase = 0;
                    end
                else
                    points(t).xc    = 0;
                    points(t).yc    = 0;
                    points(t).R     = 0;
                    points(t).Rmax  = 0;
                    points(t).e     = 0;
                    points(t).phase = 0;
                end
            end
            X     = real(double(cat(1,points.xc)));
            Y     = real(double(cat(1,points.yc)));
            R     = real(double(cat(1,points.R))); 
            Rmax  = real(double(cat(1,points.Rmax)));
            e     = real(double(cat(1,points.e)));
            phase = real(double(cat(1,points.phase)));

            % remove too small radii or Centres from this z-range. 
            index   = find(R(1:length(zind))>0.1);%,find(X(1:length(zind))>1));
            %index = intersect(tmp,find(Y(1:length(zind))>1)); 
            zind  = zind(index);
            
            
            % For fitted ellipse, phase is only mod pi, and since we fit
            % sine, we group phases -pi/2 < phase < pi/2 by addition or 
            % subtraction of integer multiples of pi;
            phase3               = phase(index);
            phase3(phase3>pi/2)  = phase3(phase3> pi/2) -  ceil(phase3(phase3> pi/2)/pi)*pi; 
            phase3(phase3<-pi/2) = phase3(phase3<-pi/2) - floor(phase3(phase3<-pi/2)/pi)*pi;
            
            
            phase(index)         = phase3;

            
            % An approximately constant phase should be in the middle of
            % the interval pi/2, or funny things happen in the fit. Therefore 
            % store the median of the phase, which we have to correctly
            % feedback into the generation of the embedding. 
            medphase = median(phase);
            this.fittedParam.medianphase = medphase;
            phase = phase-medphase;
            phase(phase>pi/2) =  phase(phase>pi/2) - pi;
            
            
            % if radii are very similar for 50 percent of all data points   
            % remove phase. Remove eccentricity correction if they are even closer. 
            if mean(abs(Rmax(index)-R(index))./(Rmax(index)+R(index)) < .001)>.5
                
                phase = 0*phase;
                e     = 0*phase;
                debugMsg(2, 'Similar major and minor axis. Removing eccentricity and phase!\n');
            elseif mean(abs(Rmax(index)-R(index))./(Rmax(index)+R(index)) < .005)>.5
           
                phase = 0*phase;
                debugMsg(2, 'Similar major and minor axis. Removing phase!\n');
            end
           
            
            % fit a radius R as a poynomial of z;
            fitPolynomial(double(zind),R,index,fitOptions.R,this,2);
            % fit slope of y0
            fitPolynomial(double(zind),Y,index,fitOptions.Y,this);
            % fit slope of x0            
            fitPolynomial(double(zind),X,index,fitOptions.X,this);
            % fit slope of e 
            fitPolynomial(double(zind),e,index,fitOptions.e,this);
            % fit slope to phase
            fitPolynomial(double(zind),phase,index,fitOptions.phase,this);
            % fit a fourier series to the sine of the phase;
            
            %fitFourier(double(zind),phase,index,fitOptions.phase,this)
            fitPolynomial(double(zind),phase,index,fitOptions.phase,this);
            
          
            % find the z-values over which fits should be generated
            this.determineZLimits(pointCloud.ROI);
            
            %---------------------------------------
            % make plots to inspect quality of fits
            %---------------------------------------
            
            %global debugOutput;
            %if debugOutput.fitQualityPlots
            pref = getpref;
            if pref.ImSAnE.fitQualityPlots 
            
            
                rsq    = real(polyval(this.fittedParam.pR, double(zind), this.fittedParam.SR,...
                 this.fittedParam.muR));                         
                xf     =  real(polyval(this.fittedParam.pX, double(zind), this.fittedParam.SX,...
                                 this.fittedParam.muX));
                yf     =  real(polyval(this.fittedParam.pY, double(zind), this.fittedParam.SY,...
                                 this.fittedParam.muY));
                ef     = real(polyval(this.fittedParam.pe,double(zind),this.fittedParam.Se,...
                    this.fittedParam.mue));
                %phasef = polyval(this.fittedParam.pphase,double(zind),this.fittedParam.Sphase,...
                %    this.fittedParam.muphase);
                %phasef = real(fitFourierEval(this,double(zind)));
                phasef     = real(polyval(this.fittedParam.pphase,double(zind),this.fittedParam.Sphase,...
                    this.fittedParam.muphase));
                
                
                %FolderName = debugOutput.dir;
                FolderName = this.fitOptions.path;
                
                h = figure('visible','off');   	    		
                plot(zind,phase(index),'r.');
                hold on
                plot(zind,phasef,'.')  , drawnow   
                legend('Data','Fit','Location','Best');
                saveas(h,fullfile(FolderName,'PhaseFit.jpg'));
                close(h);

                h = figure('visible','off');
                plot(zind,e(index),'r.');
                hold on
                plot(zind,ef,'.')     , drawnow   
                legend('Data','Fit','Location','Best');
                saveas(h,fullfile(FolderName,'EccentricityFit.jpg'));
                close(h);

                h = figure('visible','off');
                plot(zind,R(index),'r.');
                hold on
                plot(zind,real(sqrt(rsq)),'.'), drawnow    
                %axis([floor(zind(1)/10)*10 ceil(zind(end)/10)*10 0 ceil(max(sqrt(rsq))/10)*10]) , drawnow   
                legend('Data','Fit','Location','Best');
                saveas(h,fullfile(FolderName,'RadiusFit.jpg'));
                close(h);

                h = figure('visible','off');
                plot(zind,X(index),'r.');
                hold on
                plot(zind,xf,'.')     , drawnow   
                legend('Data','Fit','Location','Best');
                saveas(h,fullfile(FolderName,'CentreXFit.jpg'));
                close(h);

                h = figure('visible','off');
                plot(zind,Y(index),'r.');
                hold on
                plot(zind,yf,'.')     , drawnow   
                legend('Data','Fit','Location','Best');
                saveas(h,fullfile(FolderName,'CentreYFit.jpg'));
                close(h);
                close all;
            end
            debugMsg(2, 'Done fitting spherelike surface\n');    

        end

        %---------------------------------------------------------------
        % Determine the z-values over which the fit should be generated
        %---------------------------------------------------------------
        
        function determineZLimits(this, ROI)
            % Determine upper and lower limits for fit along z-axis.
            %
            % determineZLimits(ROI)
            %
            % determine the minimal and maximal z, either as zero crossings
            % of the radius or minima close to the pointcloud edge
            % also determine the zs for which the radius is maximal, so we
            % know the z-ranges over which we can think of Z(x,y) as a
            % (single-valued) function
            
            
            
            % 1) Coursly determine boundary given by the stack/ROI and alignment
            % 2) Refine based on largest positive connected component of radius
            % squared   
            z = ROI.zpRange(1):ROI.zpRange(end);
            
            % based on rsq, narrow down the boundary;
            rsq   = polyval(this.fittedParam.pR, double(z), this.fittedParam.SR,...
                            this.fittedParam.muR); 
                    
            % the below identification of the boundary has potential issues, 
            % especially when the fit polynomial has multiple zero crossings.            
            % zmin  = z(find(rsq>0,1,'first'));
            % zmax  = z(find(rsq>0,1,'last'));

            % estimate of zmin, zmax is based on identification of 
            % connected regions with positive rsq. Pick the largest one.
            rsqgt0 = rsq>0; % find where the squared radius is greather than zero;
            % determine connected components
            CC = bwconncomp(rsqgt0);
            if CC.NumObjects>1
                lcc = zeros(1,CC.NumObjects); 
                for i = 1 : CC.NumObjects

                    lcc(i) = length(CC.PixelIdxList{i});
                end
                % take the boundaries of largest connected component for zmin and zmax
                zminIdx = min(CC.PixelIdxList{lcc == max(lcc)}); %find
                zmaxIdx = max(CC.PixelIdxList{lcc == max(lcc)}); %find
            else
                zminIdx = min(CC.PixelIdxList{1});
                zmaxIdx = max(CC.PixelIdxList{1});
            end
            
            % check for local minima between zminIdx and zmaxIdx;
            % identify middle of z-range
            zmedIdx = floor((zmaxIdx+zminIdx)/2);
   
            
            % Take rsq in lower and upper half. Determine minmum and overwrite zmin by it; 
            rsqlow  = rsq(zminIdx:zmedIdx);
            % find the minimum;
            [~,ii] = min(rsqlow);
            % overwrite zminIdx, by shifting it by ii-1 (ii = 1 corresponds
            % to zminIdx);
            zminIdx = zminIdx+ii-1; 
            
            % same as above for upper half;
            rsqhigh = rsq(zmedIdx+1:zmaxIdx);    
            % find the minimum;
            [~,ii] = min(rsqhigh);
            % overwrite zmaxIdx, by shifting it by ii-1 (ii = 1 corresponds
            % to zmedIdx+1*);
            zmaxIdx = zmedIdx+ii;
            
            this.fittedParam.zmin  = z(zminIdx);
            this.fittedParam.zmax  = z(zmaxIdx);
            
            % after estimating the boundary remove outside values in rsq, as they
            % are a source of potential trouble
            rsq(1:max(zminIdx-1,1))     = 0;
            rsq(min(zmaxIdx+1,end):end) = 0;
            
            % TODO: do this with local maxima for shape that have multiple
            % maxima in radius
            Rmax = max(rsq);
            zRmax = z(rsq == Rmax);
            
            % zRmax can have multiple Values
            this.fittedParam.zRmaxLower = min(zRmax);
            this.fittedParam.zRmaxUpper = max(zRmax);
        end
    
        % ------------------------------------------------------
        %  Fit a polynomial of order n to data, and store result in fit
        %  Parameters
        % ------------------------------------------------------
        
        function fitPolynomial(X,Y,index,n,sF,power)
            % Fit a polynomial to data. 
            %
            % fitPolynomial(X,Y,index,n,sF,power)
            %
            % Find polynomial fit p(X) of order n to Y using polyfit. 
            % Y can be the radius, or the two coordinates of ellipse
            % centers.
            % index: Subset used for fitting only on Y;
            % The result is stored in fittedParam. 
            if nargin < 6
                power = 1; 
            end
            [p,S,mu] = polyfit(X,(Y(index).^power)',n);
            sF.fittedParam.(['p',inputname(2)])  = p;
            sF.fittedParam.(['S',inputname(2)])  = S;
            sF.fittedParam.(['mu',inputname(2)]) = mu;
        end
        
        function fitFourier(X,Y,index,n,sF,power)
            % Fit a fourier series to data.
            %
            % fitFourier(X,Y,index,n,sF,power)
            %
            % Find fourier series expression for Y using Fseries.m.
            % index: Subset used for fitting only on Y;
            % The result is stored in fittedParam. 
            if nargin < 6
                power = 1; 
            end
            [a,b] = Fseries(X,(Y(index).^power)',n,true);
            sF.fittedParam.(['a',inputname(2)])  = a;
            sF.fittedParam.(['b',inputname(2)])  = b;
            scale = [min(X),max(X)];%pi*(2*(x-x1)/(x2-x1) -1);
            sF.fittedParam.(['scale',inputname(2)]) = scale;

        end
        
        function ffit = fitFourierEval(sF,xin)
            % Evaluate fourier series fit.
            %
            % fitFourierEval(sF,xin)
            %
            % Evaluate fourier fit on xin. 
            
            x = xin(:)';
            % coefficeients from fourier fit
            a     = sF.fittedParam.aphase;
            b     = sF.fittedParam.bphase;
            % boundaries used in initial fourier fit scaling to -pi,pi
            % interval;
            scale = sF.fittedParam.scalephase;
     
            % rescale the input data and produce the output fit;
            scale_x = pi*(2*(x-scale(1))/(scale(2)-scale(1)) -1);
            fc      = @(x,n) cos(n'*x);
            fs      = @(x,n) sin(n'*x);
            if length(a)>1
                ffit    = a(2:end)'*fc(scale_x,1:length(b)) + b'*fs(scale_x,1:length(b))+a(1)/2;
            else 
                ffit    = a(1)/2*ones(size(x));
            end
            
            ffit    = reshape(ffit',size(xin));
            
        end
        
    end
end
