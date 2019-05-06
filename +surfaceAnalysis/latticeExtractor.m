 classdef latticeExtractor < handle
    %latticeExtractor Extract lattice from segmentations and ` in embedding space  
    
    properties 
        eGrids   % Embedding grids
        SOI      % surface of interest
        seg      % segmentation; Object of ilastikSegmentor class
        currLatt % current lattice; 
        currLatt2
        cellId
        lattices % structure to hold the lattices passed. 
    end

    %---------------------------------------------------------------------
    % public methods
    %---------------------------------------------------------------------
    
    methods
        
        %------------------------------------------------------
        % constructor
        %------------------------------------------------------
        
        function this = latticeExtractor(SurfOfInterest,Segmentation)
            % Constructor 
            % latticeExtractor(SurfOfInterest,Segmentation) % Specify SOI
            % and segmentation classe;
            
            this.SOI    = SurfOfInterest;
            this.seg    = Segmentation;
            this.eGrids = struct();    
            fNames      = fieldnames(this.seg.cellLabels);
            
            for i = 1 : length(fNames)
 
                this.eGrids.(fNames{i}) = this.SOI.embedding.getPatch(...
                    this.SOI.atlas.getChart(fNames{i}).domain.name).compose(...
                    this.SOI.atlas.getChart(fNames{i}).getInverse()).apply();
            end
            
            
        end
        
        %------------------------------------------------------
        % extract lattice from Segmentation; 
        %------------------------------------------------------
        
        function extractLattice(this,fName,cellLabel,xp)
            % extract Lattice from Segmentation 
            % extractLattice(this,fName,cellLabel) provide field Name e.g.
            % cylinder1_proper.
            
            % get vertices
            vertices  = bwmorph(~(cellLabel>0), 'branchpoints');
            %imshow(cat(3, mat2gray(membrane), mat2gray(vertices), zeros(size(membrane))));

            % make a list of the positions of all the vertices
            [vertexY, vertexX] = ind2sub(size(vertices), find(vertices));

            % make a cell array tracking which vertices belong to which cells
            ncells =  max(cellLabel(:));
            cellVertices = cell([1 ncells]);

            % Extract cell properties, such as centre of mass from image;
            % Remember that corner images 
            sp = regionprops(cellLabel);
            % Centroids = cat(1,sp.Centroid);
            Areas     = cat(1,sp.Area);
            % for each vertex, add index to cells that are 8-connected to it
            for i=1:length(vertexX)
                belongsto = unique(cellLabel(vertexY(i)-1:vertexY(i)+1,vertexX(i)-1:vertexX(i)+1));
                for j = 2:length(belongsto)
                    cellVertices{belongsto(j)} = [cellVertices{belongsto(j)}, i];
                end
            end
            
            vertexPosition = [vertexX, vertexY];
            % in at least one example, the above routine produced a cell
            % with 2 bonds only. This can't be! 
            CellsGLatt = cellVertices(Areas<max(Areas));
            ltt  = zeros(size(CellsGLatt));
            lbs  = zeros(size(CellsGLatt)); 
            for i = 1 : length(CellsGLatt)
                
                ltt(i) =  length(CellsGLatt{i});
                if ~isempty(CellsGLatt{i}) && ltt(i)>2;
                    vl     =  vertexPosition(CellsGLatt{i}(1:end),:)-vertexPosition(CellsGLatt{i}([2:end,1]),:);
                    vl     =  sqrt(sum(vl.^2,2));
                    if any(vl>1000)
                        lbs(i) = 1;
                    end
                end
            end
            % exclude the case of cells with 2 vertices, we require at
            % least 3
            CellsGLatt = CellsGLatt(ltt>2);
            
            % exlude all cells with too long bonds; 
            lbs        = lbs(ltt>2);
            CellsGLatt = CellsGLatt(lbs == 0);
            
            
            this.currLatt = GLattConversion(CellsGLatt, vertexPosition);
            
            % Compute the verticex coordinates in the embedding space. 
            verts3 = this.currLatt.verts;
            for i = 1 : length(this.currLatt.verts)
                verts3(i,1) = this.eGrids.(fName){1}(this.currLatt.verts(i,2),this.currLatt.verts(i,1));
                verts3(i,2) = this.eGrids.(fName){2}(this.currLatt.verts(i,2),this.currLatt.verts(i,1));
                verts3(i,3) = this.eGrids.(fName){3}(this.currLatt.verts(i,2),this.currLatt.verts(i,1));
            end
            this.currLatt.verts3 = verts3;
           
            % Compute CoM of each cell in the lattices; 
            this.currLatt.CoM  = NaN*ones(length(this.currLatt.cells),3);
            this.currLatt.CoM3 = NaN*ones(length(this.currLatt.cells),3);
            this.currLatt.phi  = NaN*ones(length(this.currLatt.cells),1);
            for i = 1 : length(this.currLatt.cells)
                currCell = this.currLatt.cells{i};
                if length(currCell) > 3
                    this.currLatt.CoM(i,:)  = mean(this.currLatt.verts(this.currLatt.bonds(currCell),:));
                    this.currLatt.CoM3(i,:) = mean(this.currLatt.verts3(this.currLatt.bonds(currCell),:));
                    this.currLatt.phi(i)    = atan2(this.currLatt.CoM3(i,2),this.currLatt.CoM3(i,1)); 
                end
            end
            
            % compute cell area; 
            this.currLatt.area = zeros(size(this.currLatt.cells));
            g = this.SOI.g.getPatch([fName '_index']);
            if strcmp('cylinder1',fName)
            eGrids = this.SOI.embedding.getPatch([fName,'_index']).apply();
            end
            if isa(g,'diffgeometry.TensorPatch')

                gim = g.apply();
                detg = sqrt(gim{1,1}.*gim{2,2}-gim{1,2}.*gim{2,1});
                [X,Y] = meshgrid(1:size(detg,2),1:size(detg,1));

                for i = 1: length(this.currLatt.cells)
                    v  = this.currLatt.verts(this.currLatt.bonds(this.currLatt.cells{i}),:);
                    IN = inpolygon(X,Y,v(:,1),v(:,2));
                    this.currLatt.area(i) = sum(sum(1*double(IN)));
                end
            end
            
            % this will be an alternative when loading the metric from
            % disc;
%             nameD = fullfile(xp.fileMeta.directory,'fields','metric',fname);
%             name = fullfile(nameD,[xp.fileMeta.filenameFormat(1:end-4),...
%                 '_metric_t%u_c%u%u_',fname,'.mat']);
%             t = 1;
%             if ~exist(sprintf(name,t,1,1))
% 
%                 error('Require metric to compute cell orientation');
%             else
% 
%                 gim = cell(2,2);
%                 for ii = 1 : 2
%                     for jj = 1 :2
%                         s = load(sprintf(name,t,ii,jj));
%                         gim{ii,jj} = s.im;
%                     end
%                 end
%             end
            
            
            % compute cell eccentricity;
            % Verbal of the code: In embedding space, compute covariance matrix of
            % each cell, get the eigenvalues and from these derive the
            % eccentricity.
            this.currLatt.eccentricity = zeros(size(this.currLatt.cells));
            this.currLatt.appolarity   = zeros(size(this.currLatt.cells));
            this.currLatt.appolarityaxis  = zeros(size(this.currLatt.cells));
            this.currLatt.orientation  = zeros(size(this.currLatt.cells));
            for k = 1 : length(this.currLatt.cells)
                verts    = (this.currLatt.bonds(this.currLatt.cells{k}(:),1));
                % the 3d coordinates of the cell around the center 
                verts2D =  [this.currLatt.verts(verts,1)-this.currLatt.CoM(k,1),...
                    this.currLatt.verts(verts,2)-this.currLatt.CoM(k,2)];
                verts3D = [this.currLatt.verts3(verts,1)-this.currLatt.CoM3(k,1),...
                    this.currLatt.verts3(verts,2)-this.currLatt.CoM3(k,2),...
                    this.currLatt.verts3(verts,3)-this.currLatt.CoM3(k,3)];

                if ~any(isnan(verts3D(:))) && ~any(isinf(verts3D(:))) && ...
                        ~isnan(this.currLatt.area(k)) && ~isinf(this.currLatt.area(k))

                    CovM(1,1) = sum((verts3D(:,1)).^2);
                    CovM(2,2) = sum((verts3D(:,2)).^2);
                    CovM(3,3) = sum((verts3D(:,3)).^2);
                    CovM(1,2) = sum((verts3D(:,1)).*(verts3D(:,2)));
                    CovM(2,1) = CovM(1,2);
                    CovM(1,3) = sum((verts3D(:,1)).*(verts3D(:,3)));
                    CovM(3,1) = CovM(1,3);
                    CovM(2,3) = sum((verts3D(:,2)).*(verts3D(:,3)));
                    CovM(3,2) = CovM(2,3);

                    CovM = CovM./(this.currLatt.area(k)+0.000001);% divsion by small number!        

                    [~,D] = eig(CovM);

                    this.currLatt.eccentricity(k) = sqrt(1-(D(2,2)./D(3,3)).^1);
                    %this.currLatt.appolarity(k)   = abs([0,0,1]*V(:,3)./norm(V(:,3)));
                else
                    this.currLatt.eccentricity(k) = NaN;
                end 
                
                % compute the orientation of each cell w.r.t. the AP axis;
                % By our convention, this coincides with the z axis. 
                
                verts = (this.currLatt.bonds(this.currLatt.cells{k}(:),1));

                verts = [this.currLatt.verts(verts,1)-this.currLatt.CoM(k,1)...
                        ,this.currLatt.verts(verts,2)-this.currLatt.CoM(k,2)];
                
                    
                    
                    
                if ~any(isnan(verts2D(:))) && ~any(isinf(verts2D(:))) && ...
                        ~isnan(this.currLatt.area(k)) && ~isinf(this.currLatt.area(k)) && ...
                    isa(g,'diffgeometry.TensorPatch') && strcmp(fName,'cylinder1')
                    % conly compute for the cylinder1&2_index charts. No
                    % need to recompute in proper charts etc.
                    
%                     CovM(1,1) = sum((verts2D(:,1)).^2);
%                     CovM(2,2) = sum((verts2D(:,2)).^2);
%                     CovM(1,2) = sum((verts2D(:,1)).*(verts2D(:,2)));
%                     CovM(2,1) = CovM(1,2);
%                     
%                     CovM = CovM./(this.currLatt.area(k)+0.000001);% divsion by small number!        

                    CovM(1,1) = sum((verts3D(:,1)).^2);
                    CovM(2,2) = sum((verts3D(:,2)).^2);
                    CovM(3,3) = sum((verts3D(:,3)).^2);
                    CovM(1,2) = sum((verts3D(:,1)).*(verts3D(:,2)));
                    CovM(2,1) = CovM(1,2);
                    CovM(1,3) = sum((verts3D(:,1)).*(verts3D(:,3)));
                    CovM(3,1) = CovM(1,3);
                    CovM(2,3) = sum((verts3D(:,2)).*(verts3D(:,3)));
                    CovM(3,2) = CovM(2,3);

                    CovM = CovM./(this.currLatt.area(k)+0.000001);% divsion by small number!        
                    
                    % in this case, the vector is in 3D, so you need to
                    % calculate it in the plane;
                    [V,D] = eig(CovM);
                    vec = V(:,3);% the leading eigenvector;
                    vec = vec/norm(vec); 
                    % add it to the cell center; 
                    vec1 = [this.currLatt.CoM3(k,1);this.currLatt.CoM3(k,2);this.currLatt.CoM3(k,3)]+2*vec;
                    vec2 = [this.currLatt.CoM3(k,1);this.currLatt.CoM3(k,2);this.currLatt.CoM3(k,3)]-2*vec;

                    d = (eGrids{1}-vec1(1)).^2+(eGrids{2}-vec1(2)).^2+(eGrids{3}-vec1(3)).^2;
                    [~,indk] = min(d(:));    
                    [I1,J1] = ind2sub(size(eGrids{1}),indk);
                    
                    d = (eGrids{1}-vec2(1)).^2+(eGrids{2}-vec2(2)).^2+(eGrids{3}-vec2(3)).^2;
                    [~,indk] = min(d(:));    
                    [I2,J2] = ind2sub(size(eGrids{1}),indk);
                    
                    
                    
                    
                    
                    %[V,D] = eig(CovM);
                 
                    %vec1 = V(:,2); % orientation of long axis in cell;
                    
                    vec1 = [I2-I1,J2-J1];
                    % adjust for step size in phi;
                    stepSize = this.SOI.atlas.charts{1}.image.stepSize;
                    vec1(2) = vec1(2)*stepSize(2);
                    vec1(1) = vec1(1)*stepSize(1);
                    vec2 = [1,0];  %e_phi
                    vec3 = [0,1]; % e_z;
                    vec2(1) = vec2(1)*stepSize(1);
                    vec2(2) = vec2(2)*stepSize(2);
                    vec3(1) = vec3(1)*stepSize(2);
                    vec3(2) = vec3(2)*stepSize(1);
                    
                    %Evaluate orientation between those vectors at the
                    %centre of the cell; 
                    curPos = round([this.currLatt.CoM(k,1),this.currLatt.CoM(k,2)]);
                    
                    gcp = zeros(2,2);
                    for ii = 1:2 
                        for jj = 1:2
                            gcp(ii,jj) = gim{ii,jj}(curPos(2),curPos(1));
                        end
                    end
                    % compute the inner product
                    temp1 = gcp(1,1)*vec1(1)*vec2(1)+gcp(1,2)*vec1(1)*vec2(2)+...
                        gcp(2,1)*vec1(2)*vec2(1)+gcp(2,2)*vec1(2)*vec2(2);
                    % compute the norms; 
                    temp2 = (gcp(1,1)*vec1(1)*vec1(1)+gcp(1,2)*vec1(1)*vec1(2)+...
                        gcp(2,1)*vec1(2)*vec1(1)+gcp(2,2)*vec1(2)*vec1(2));
                    temp3 = (gcp(1,1)*vec2(1)*vec2(1)+gcp(1,2)*vec2(1)*vec2(2)+...
                        gcp(2,1)*vec2(2)*vec2(1)+gcp(2,2)*vec2(2)*vec2(2));
                    
                    cosT = sqrt(temp1.^2./temp2./temp3); % matlab does things internally in rad not deg!
                    
                    this.currLatt.orientation(k) = acos(cosT);
                    
                    %compute the inner product
                    temp1 = gcp(1,1)*vec1(1)*vec3(1)+gcp(1,2)*vec1(1)*vec3(2)+...
                        gcp(2,1)*vec1(2)*vec3(1)+gcp(2,2)*vec1(2)*vec3(2);
                    % compute the norms; 
                    temp2 = (gcp(1,1)*vec1(1)*vec1(1)+gcp(1,2)*vec1(1)*vec1(2)+...
                        gcp(2,1)*vec1(2)*vec1(1)+gcp(2,2)*vec1(2)*vec1(2));
                    temp3 = (gcp(1,1)*vec3(1)*vec3(1)+gcp(1,2)*vec3(1)*vec3(2)+...
                        gcp(2,1)*vec3(2)*vec3(1)+gcp(2,2)*vec3(2)*vec3(2));
                    % the inner product between the z axis and the vector;
                    eZeV = sqrt(temp1.^2./temp2./temp3);
                    % multiply by length the cell;
                    tempL = this.currLatt.CoM3;
                    tempL = max(tempL(:,3))-min(tempL(:,3));
                    this.currLatt.appolarityaxis(k) = sqrt(temp2)*eZeV;
                    
                    
                    
                else
                    this.currLatt.orientation(k) = NaN;
                    this.currLatt.appolarityaxis(k) = NaN;
                end
                
                verts = (this.currLatt.bonds(this.currLatt.cells{k}(:),1));

                verts = [this.currLatt.verts(verts,1)-this.currLatt.CoM(k,1)...
                        ,this.currLatt.verts(verts,2)-this.currLatt.CoM(k,2)];
                
                if ~any(isnan(verts(:))) && ~any(isinf(verts(:))) && ...
                        ~isnan(this.currLatt.area(k)) && ~isinf(this.currLatt.area(k))
                    % repeat this in 2D to get an idea of the orientation
                    % of the cell
                   

%                     CovM2(1,1) = sum((verts(:,1)).^2);
%                     CovM2(2,2) = sum((verts(:,2)).^2);
%                     CovM2(1,2) = sum((verts(:,1)).*(verts(:,2)));
%                     CovM2(2,1) = CovM2(1,2);
% 
%                     CovM2 = CovM2./(this.currLatt.area(k)++0.000001);
% 
%                     [V,D] = eig(CovM2);
% 
%                     temp = atan2(V(2,2),V(1,2));
%                     if temp < 0
%                         temp = temp+pi;
%                     end
%                    this.currLatt.appolarity(k)  = temp;
                    % compute the absolute length of th AP axis.
                    temp = this.currLatt.bonds(:,1);
                    temp = this.currLatt.verts(temp(~isnan(temp)),1);
                    temp = max(temp) - min(temp);
                    % and then compare the relative contribution. 
                    this.currLatt.appolarity(k)  = (max(verts(:,1))-min(verts(:,1)))/temp;
                else
                    this.currLatt.appolarity(k) = NaN;
                end
            end
            
            
        end
       

        % ------------------------------------------------------
        %  batch batch extract lattice structure from segmentation; 
        % ------------------------------------------------------
        
        function batchExtractLattice(this)
            
            this.lattices = struct([]);
            
            fNames  = fieldnames(this.seg.cellLabels);
            nFields = length(fNames);
            

            for t = 1 : length(this.seg.batchTimes)
                    
                for f = 1 : nFields
                    cellLabel = this.seg.cellLabels(t).(fNames{f});
                    %if ~isempty(this.seg.prediction)
                        
                        disp(['Extracting Time Point ',num2str(t),' Lattice : ', fNames{f}])
                        this.extractLattice(fNames{f},cellLabel);
                        this.lattices(t).(fNames{f}) = this.currLatt;
                    %end
                end
            end
        end
        
        % ------------------------------------------------------
        %  fuse Lattices from spherelike SOI;  
        % ------------------------------------------------------
        
        function fuseSphereLikeLattices(this,t)
            % fuseSphereLikeLattices
            
                
            lattice1 = this.lattices(t).cylinder1;
            lattice1 = this.normalOutwardOrientation(lattice1);
            
            %phi1 = atan2(lattice1.CoM3(:,2),lattice1.CoM3(:,1));
            %in1  = find(abs(phi1)<3*pi/2);
            maxPhi = max(lattice1.CoM(:,2));
            minPhi = min(lattice1.CoM(:,2));
            meanPhi = mean(lattice1.CoM(~isnan(lattice1.CoM(:,2)),2));
            in1  = intersect(find(lattice1.CoM(:,2)>meanPhi/2),find(lattice1.CoM(:,2)<3*meanPhi/2));
            
            
            
            % reduce the first lattice;
            this.subLattice(lattice1,in1);
            latt = this.currLatt;
            % and make normals of cells point outwards
            latt = this.normalOutwardOrientation(latt); 
            % find the boundary of latt;
            index = this.boundaryCells(latt);
            
            
            % get the first lattice to fuse to; 
            lattice2 = this.lattices(t).cylinder2;
            % and make normals of cells point outwards
            lattice2 = this.normalOutwardOrientation(lattice2);
            
            
            % match all boundary cells in latt to the second lattice; 
            dmax = 20;
            [ix,dist] = track3D(lattice2.CoM3(:,1),lattice2.CoM3(:,2),lattice2.CoM3(:,3),...
                latt.CoM3(index,1),latt.CoM3(index,2),latt.CoM3(index,3),dmax);

            
       
            [latt,fuseInLattice] = this.adaptCells(latt,lattice2,ix,index);
                      
       
            
            %inCand = find(abs(fuseInLattice.phi)<=(pi/2+1/180*pi));
            
            inCand = intersect(find(fuseInLattice.CoM(:,2)>meanPhi/2),find(fuseInLattice.CoM(:,2)<3*meanPhi/2));
            ix(dist>10) = 0;
%             outCand = ix(dist>10);
%             inCand = setdiff(inCand,outCand);
%             inCand = setdiff(inCand,ix);
        
            finLattice     = this.stitchCorrectedLattices(latt,fuseInLattice,index,ix,inCand);
            
            
            this.currLatt  = finLattice;
            this.currLatt2 = fuseInLattice;
            
            this.cellId = index;
            
            
% add the caps! If you don't want this, then remove from here, 

            lattice3 = this.lattices(t).anteriorEquidistant;
            number = -300;
            in3 = lattice3.CoM3(:,3)<=number;
            this.subLattice(lattice3,in3);
            latt3 = this.currLatt;
          
            index3 = this.boundaryCells(latt3);

            [ix3,dist] = track3D(finLattice.CoM3(:,1),finLattice.CoM3(:,2),finLattice.CoM3(:,3),...
                latt3.CoM3(index3,1),latt3.CoM3(index3,2),latt3.CoM3(index3,3),20);
             
         
            [latt3,fuseInLattice3] = this.adaptCells(latt3,finLattice,ix3,index3);
            disp('done adapting polar UpperZ to the fusedLattice');
             this.currLatt  = latt3;
             this.currLatt2 = fuseInLattice3;
%             
    
    
    
            inCand3 = find(fuseInLattice3.CoM3(:,3)>=number);
            %ix3(dist>10) = 0;

             finLattice2 = this.stitchCorrectedLattices(latt3,fuseInLattice3,index3,ix3,inCand3);
            
         
            
            
            
            
            this.currLatt  = finLattice2;
            this.currLatt2 = latt3;
            
            
            % add the second pole: 
            
             lattice4 = this.lattices(t).posteriorEquidistant;
             number = 120;
             in4 = lattice4.CoM3(:,3)>=number;
             this.subLattice(lattice4,in4);
             latt4 = this.currLatt;
           
             index4 = this.boundaryCells(latt4); % the boundary cells on the pole.
 
             [ix4,dist] = track3D(finLattice2.CoM3(:,1),finLattice2.CoM3(:,2),finLattice2.CoM3(:,3),...
                 latt4.CoM3(index4,1),latt4.CoM3(index4,2),latt4.CoM3(index4,3),20);
              
          
             [latt4,fuseInLattice4] = this.adaptCells(latt4,finLattice2,ix4,index4);
             disp('done adapting polar UpperZ to the fusedLattice');
              this.currLatt  = latt4;
              this.currLatt2 = fuseInLattice4;
              
% %             
    
    
    
            inCand4 = find(fuseInLattice4.CoM3(:,3)<=number);
            ix4(dist>10) = 0;
            
            
            % there is a parity transformation between the two lattices,
            % that prevents proper fusion. To overcome, we change the
            % handedness convention in every cell of latt4:
            
            lattice1 = latt4;
            lattice1_new = lattice1;
            for i = 1 : length(lattice1.cells)
            %     lattice1.cells{i} = lattice1.cells{i}(end:-1:1);
            % change handedNess of cell;
            currCell1 = lattice1.cells{i};
                [conjCurrCell1,conjCurrBonds1,currBonds1] = getConjCell(lattice1,currCell1);
                lattice1_new.cells{i} = conjCurrCell1(end:-1:1);
            end
            latt4 = lattice1_new;

             finLattice3 = this.stitchCorrectedLattices(latt4,fuseInLattice4,index4,ix4,inCand4);
                                                        
            
         
            
            
            
            
            this.currLatt  = finLattice3;
            this.currLatt2 = latt4;
            
            
            
            
% stop commenting here            
            
        end
        
        
        % ------------------------------------------------------
        %  adapt overlapping cells in fusable lattices;
        % ------------------------------------------------------
        
        function [outLattice1,outLattice2] = adaptCells(this,lattice1,lattice2,ix,index)
        
      
            % remember: lattice2 is matched to index in lattice. (ix has size of index)

            countCell1Bigger = 0;
            countCell2Bigger = 0;


            verts1 = lattice1.verts3;
            verts2 = lattice2.verts3;

            bondsExcluded1 = zeros(1,length(lattice1.bonds))>0;
            bondsExcluded2 = zeros(1,length(lattice2.bonds))>0;
            % these list the bonds that are already considered to be taken out; 
            for i = 1 : length(index)

                if ix(i)>0

                    currBondsCell1 = lattice1.cells{index(i)};
                    currBondsCell2 = lattice2.cells{ix(i)};

                    % three cases:1) cell1 coordination number is bigger
                    %             2) cell1 coordination number is smaller; 
                    %             3) cells are of same coordination number, do nothing

                    if sum(~bondsExcluded1(currBondsCell1)) > sum(~bondsExcluded2(currBondsCell2))

                        % Case 1)  We only need to make changes in cell 1; 

                        % determine by how much does the coordination number differ
                        difference       = sum(~bondsExcluded1(currBondsCell1)) - sum(~bondsExcluded2(currBondsCell2));
                        countCell1Bigger = countCell1Bigger + difference;

                        % as where almost doing the same things (just different inputs)
                        % generate a subfunction excludedBonds that takes care of the
                        % bond exclusion;
                        for k = 1 : difference
                            bondsExcluded1 = this.findExclusionBonds(bondsExcluded1,currBondsCell1,lattice1,verts1,index,i);
                        end
                        %error('stop there')
                    elseif sum(~bondsExcluded1(currBondsCell1)) < sum(~bondsExcluded2(currBondsCell2)) 

                        % Case 2)  We only need to make changes in cell 2; 

                        % determine by how much does the coordination number differ
                        difference       = sum(~bondsExcluded2(currBondsCell2)) - sum(~bondsExcluded1(currBondsCell1));
                        countCell2Bigger = countCell2Bigger +difference;
                        % repeat as above just replace the index from 1 -> 2 and index
                        % by ix;
                        for k = 1 : difference
                            bondsExcluded2 = this.findExclusionBonds(bondsExcluded2,currBondsCell2,lattice2,verts2,ix,i);
                        end

                    end
                end
            end
            % 
            % sanity check: Each run with cell1 bigger cell2 produces 2 bonds to be
            % excluded; 

            if countCell1Bigger*2 == length(find(bondsExcluded1==1))
                disp('Correct number of bonds excluded in lattice 1.');
            else
                disp('Cell1 bigger cell2: Number of bonds excluded does not match number of runs*2.');
            end

            if countCell2Bigger*2 == length(find(bondsExcluded2==1))
                disp('Correct number of bonds excluded in lattice 2.');
            else
                disp('Cell1 smaller cell2: Number of bonds excluded does not match number of runs*2.');
            end

            % do the consistency check, i.e. chack that cells have same coordination number now. 
            wrong = 0*index;
            for i = 1 : length(index)

                if ix(i)>0

                    currBondsCell1 = lattice1.cells{index(i)};
                    currBondsCell2 = lattice2.cells{ix(i)};

                    % three cases:1) cell1 coordination number is bigger
                    %             2) cell1 coordination number is smaller; 
                    %             3) cells are of same coordination number, do nothing

                    if sum(~bondsExcluded1(currBondsCell1)) ~= sum(~bondsExcluded2(currBondsCell2))
                       % error('wrong!') 
                        wrong(i) = sum(~bondsExcluded1(currBondsCell1))-sum(~bondsExcluded2(currBondsCell2));
                    end
                end
            end

            if sum(wrong) == 0

                disp('All matched cells have same coordination number');
            else
                mistakes = find(wrong);
                disp(['Mistakes in the coordination Number at index: ',num2str(mistakes)])
            end


            % exclude the bonds; 
            outLattice1 = this.excludeBonds(lattice1,bondsExcluded1);
            % Remove all parts that aren't used
            outLattice1 = this.indexedSubLattice(outLattice1,1:length(outLattice1.cells));
            outLattice2 = this.excludeBonds(lattice2,bondsExcluded2);
            % Remove all parts that aren't used
            outLattice2 = this.indexedSubLattice(outLattice2,1:length(outLattice2.cells));
            % and again some sanity checking; 
            wrong = 0*index;
            for i = 1 : length(index)

                if ix(i)>0

                    currBondsCell1 = outLattice1.cells{index(i)};
                    currBondsCell2 = outLattice2.cells{ix(i)};

                    % three cases:1) cell1 coordination number is bigger
                    %             2) cell1 coordination number is smaller; 
                    %             3) cells are of same coordination number, do nothing

                    if length(currBondsCell1) ~= length(currBondsCell2)
                       % error('wrong!') 
                        wrong(i) = sum(~bondsExcluded1(currBondsCell1))-sum(~bondsExcluded2(currBondsCell2));
                    end
                end
            end

            if sum(wrong) == 0
                disp('After Matching, all matched cells have same coordination number');
            else
                mistakes = find(wrong);
                disp(['After Matching, there are mistakes in the coordination Number at index: ',num2str(mistakes)])
            end
        end
        
 
        
        % ------------------------------------------------------
        %  stitch corrected lattices;
        % ------------------------------------------------------
        
        function finOutLattice = stitchCorrectedLattices(this,lattice1,lattice2,index,ix,inCand)

            
            % start by point matching the bonds from lattice2 to bonds in
            % lattice1;
            
            bDistMax = 95; % a parameter that checks the cummulative distance of all matched bonds;
            % make the above thing an input into the function;
            bonds2Ix = zeros(length(lattice2.bonds),1);
            dontMatchCells = 0*ix;
            for i = 1 : length(ix)

                if ix(i)>0

                    % current cells
                    currCell1 = lattice1.cells{index(i)};
                    currCell2 = lattice2.cells{ix(i)};

                     % conjugate current cell, conjugate bonds and current bonds;
                    [conjCurrCell1,conjCurrBonds1,currBonds1] = getConjCell(lattice1,currCell1);
                    [conjCurrCell2,conjCurrBonds2,currBonds2] = getConjCell(lattice2,currCell2);

                    % current Vertices;
                    currv11  = lattice1.verts3(currBonds1(:,1),:);
                    currv12  = lattice1.verts3(currBonds1(:,2),:);
                    currv21  = lattice2.verts3(currBonds2(:,1),:);
                    currv22  = lattice2.verts3(currBonds2(:,2),:);

                    % center of mass of each bond;
                    currCoM1 = (currv11+currv12)/2;
                    currCoM2 = (currv21+currv22)/2;

                    % assign bonds such that orientation is preserved in both cells and closest partners are chosen; 
                    [bondIx,bDist] = assignBonds(currCoM1(:,1),currCoM1(:,2),currCoM1(:,3),...
                         currCoM2(:,1),currCoM2(:,2),currCoM2(:,3));

                    if min(bDist)<bDistMax 
                        
                        % check if the current mapping doesn't violate prior mappings; 
                        % 1) A bond in lattice2 might already point to a bond in lattice1 or
                        % 2) A bond in lattice1 might already be mapped to a different bond in
                        %   lattice2
                        flag = 0;

                        if any(bonds2Ix(currCell2)~=0) || any(bonds2Ix(conjCurrCell2)~=0)



                            temp = find(bonds2Ix(currCell2)~=0);
                            if any( ( bonds2Ix(currCell2(temp)) - currCell1(bondIx(temp))' ) ~= 0)


                                flag = 1;
                                dontMatchCells(i) = 1;
                                disp('A bond in lattice2 is already mapped and now a different one is mapped to it! Violation!Make this a debug message!');
                            end
                        end
                        if any(ismember(currCell1(bondIx),bonds2Ix)) || any(ismember(conjCurrCell1(bondIx),bonds2Ix))


                            temp = find(ismember(currCell1(bondIx),bonds2Ix));
                            if any( ( bonds2Ix(currCell2(temp)) - currCell1(bondIx(temp))' ) ~= 0)

                                flag = 1;
                                dontMatchCells(i) = 1;
                                disp('A bond in lattice1 is already mapped; Make this a debug message!');
                            end
                        end

                        if flag == 0    
                            
                            % all the bonds in currCell2 map to the bonds in currCell1(bondIx);
                            bonds2Ix(currCell2,1)     = currCell1(bondIx);
                            % do the same for the conjugates; 
                            bonds2Ix(conjCurrCell2,1) = conjCurrCell1(bondIx);
                        end

                    else
                        dontMatchCells(i) = 1;
                    end
                end
            end
            
            % bonds2Ix is a list of length(lattice2.bonds) pointing into the
            % bonds in lattice1; (point matching index); 
            % fBondsIx is just finding where a partner bond was identified.
            fBonds2Ix = find(bonds2Ix~=0);
            
            
            %%%%%%
            % Actual lattice fusion 
            %%%%%%

            % generate the output lattice based on the first lattice; Add
            % lattice 2 into it. 
            outLattice = lattice1;

            ncells = length(lattice1.cells);
            nbonds = length(lattice1.bonds);
            nverts = length(lattice1.verts3);

%            cellIndex     = [1:ncells,unique(inCand)'+ncells]';    

            for i = 1 : length(lattice2.cells)
                temp = lattice2.cells{i}+nbonds;
                outLattice.cells{end+1} = temp; 
            end
            
            outLattice.CoM3(end+1:end+length(lattice2.cells),:)  = lattice2.CoM3;
            outLattice.CoM(end+1:end+length(lattice2.cells),:)   = lattice2.CoM;
            outLattice.phi(end+1:end+length(lattice2.cells),:)   = lattice2.phi;
            outLattice.bonds  = [outLattice.bonds;[lattice2.bonds(:,1)+nverts,lattice2.bonds(:,2)+nverts,lattice2.bonds(:,3)+ncells*(lattice2.bonds(:,3)~=0),lattice2.bonds(:,4)+ncells*(lattice2.bonds(:,4)~=0)]];
            outLattice.verts  = [outLattice.verts;lattice2.verts];
            outLattice.verts3 = [outLattice.verts3;lattice2.verts3];
            
            % %% Loop through all the matched bonds and replace entries of lattice2 by entries of lattice1;  
            % don't forget to add nbonds to the index into bonds2Ix, as
            % this now points into the outLattice;
            for i = 1 : length(fBonds2Ix)

                % the name of the bond that goes out 
                bondOutName     = fBonds2Ix(i);
                % get the conjugates, too! 
                conjBondOutName = getConjCell(lattice2,bondOutName);

                % the name of the bond that replaces it; 
                bondReplaceName = bonds2Ix(fBonds2Ix(i));
                % get the conjugates, too!
                conjBondReplaceName = getConjCell(lattice1,bondReplaceName);

                % the vertices of the old bonds must be replaced in all incident bonds;
                % first vertex appears in these bonds
                vertsOut1 =  find(lattice2.bonds(:,1) == lattice2.bonds(bondOutName,1));
                vertsOut2 =  find(lattice2.bonds(:,2) == lattice2.bonds(bondOutName,1));
                % second vertex appears in these bonds
                vertsOut3 =  find(lattice2.bonds(:,1) == lattice2.bonds(bondOutName,2));
                vertsOut4 =  find(lattice2.bonds(:,2) == lattice2.bonds(bondOutName,2));

      %          disp('We dont need the dontUpdateBondFlag in latticeExtractor anymore! Make this a debug message!');
      %          dontUpdateBondFlag = 0;
                % replace the vertex indices
                if any(outLattice.bonds(vertsOut1+nbonds,1) < nverts) && any(lattice1.bonds(bondReplaceName,1) ~= outLattice.bonds(vertsOut1+nbonds,1))
                    % if a bond in the lattice2 part of outlattice is from lattice1 gets updated by another bond in
                    % lattice1: fuse the two and save them as the prior
                    % changed one; 
                    
      %              dontUpdateBondFlag = 0;

                    % fuse these two vertices; 
                    v1   = lattice1.verts(lattice1.bonds(bondReplaceName,1),:);
                    v2   = lattice1.verts(unique(outLattice.bonds(vertsOut1+nbonds,1)),:);
                    v13D = lattice1.verts3(lattice1.bonds(bondReplaceName,1),:);
                    v23D = lattice1.verts3(unique(outLattice.bonds(vertsOut1+nbonds,1)),:);

                    av   = (v1+v2)/2;
                    av3D = (v13D+v23D)/2;
                    lattice1.verts3(lattice1.bonds(unique(outLattice.bonds(vertsOut1+nbonds,1)),1),:) = av3D;
                    lattice1.verts(lattice1.bonds(unique(outLattice.bonds(vertsOut1+nbonds,1)),1),:)  = av;

                    outLattice.bonds(outLattice.bonds(:,1) == lattice1.bonds(bondReplaceName,1),1) = unique(outLattice.bonds(vertsOut1+nbonds,1));
                    outLattice.bonds(outLattice.bonds(:,2) == lattice1.bonds(bondReplaceName,1),2) = unique(outLattice.bonds(vertsOut1+nbonds,1));

                end

      %          if any(outLattice.bonds(vertsOut2+nbonds,2) < nverts) && any(lattice1.bonds(bondReplaceName,1) ~= outLattice.bonds(vertsOut2+nbonds,2))
      %              dontUpdateBondFlag = 0;
      %          end

                if any(outLattice.bonds(vertsOut3+nbonds,1) < nverts) && any(lattice1.bonds(bondReplaceName,2) ~= outLattice.bonds(vertsOut3+nbonds,1))
      %              dontUpdateBondFlag = 0;        
                    % fuse these two vertices; 
                    v1   = lattice1.verts(lattice1.bonds(bondReplaceName,2),:);
                    v2   = lattice1.verts(unique(outLattice.bonds(vertsOut3+nbonds,1)),:);
                    v13D = lattice1.verts3(lattice1.bonds(bondReplaceName,2),:);
                    v23D = lattice1.verts3(unique(outLattice.bonds(vertsOut3+nbonds,1)),:);

                    av   = (v1+v2)/2;
                    av3D = (v13D+v23D)/2;
                    lattice1.verts3(lattice1.bonds(unique(outLattice.bonds(vertsOut3+nbonds,1)),2),:) = av3D;
                    lattice1.verts(lattice1.bonds(unique(outLattice.bonds(vertsOut3+nbonds,1)),2),:)  = av;
                    
                    outLattice.bonds(outLattice.bonds(:,1) == lattice1.bonds(bondReplaceName,2),1) = unique(outLattice.bonds(vertsOut3+nbonds,1));
                    outLattice.bonds(outLattice.bonds(:,2) == lattice1.bonds(bondReplaceName,2),2) = unique(outLattice.bonds(vertsOut3+nbonds,1));
                end

      %          if any(outLattice.bonds(vertsOut4+nbonds,2) < nverts) && any(lattice1.bonds(bondReplaceName,2) ~= outLattice.bonds(vertsOut4+nbonds,2))
      %              dontUpdateBondFlag = 0;
      %          end

      %          if dontUpdateBondFlag == 0;

                    outLattice.bonds(vertsOut1+nbonds,1) = lattice1.bonds(bondReplaceName,1);
                    outLattice.bonds(vertsOut2+nbonds,2) = lattice1.bonds(bondReplaceName,1);
                    outLattice.bonds(vertsOut3+nbonds,1) = lattice1.bonds(bondReplaceName,2);
                    outLattice.bonds(vertsOut4+nbonds,2) = lattice1.bonds(bondReplaceName,2);

      %          end
                % now update the neighborbood relations

                swap   = [2,1];

                cellNamesLattice2     = lattice2.bonds(bondOutName,3:4);
                cellNamesLattice1     = lattice1.bonds(bondReplaceName,3:4);

     %           if dontUpdateBondFlag == 0

                if ~any(cellNamesLattice1 == 0)

                    % case1: cellNamesLattice1 doesn't contain the zero cell; 
                    % in bonds of lattice2 replace all the cell names that are matched with
                    % cells in lattice1 by the cell name of lattice1; Remove the
                    % cellNamesLattice2 from the inCand list; Do the same for the
                    % conjugates; Accordingly adjust the old bonds from lattice2; 
                    inCand = setdiff(inCand,cellNamesLattice2');
                    mem    = find(ismember(cellNamesLattice2,ix));
                    mem    = mem(cellNamesLattice2(mem)~=0); 


                    for k = 1 : length(mem)
                        outLattice.bonds(bondOutName+nbonds,mem(k)+2) = cellNamesLattice1(mem(k));
                        % conjugate
                        outLattice.bonds(conjBondOutName+nbonds,swap(mem(k))+2) = cellNamesLattice1(mem(k));
                    end
                else

                    % case2: cellNamesLattice1 contains the zero cell;
                    % Update the bond in lattice1: Replace the 0 by the neighbor that does
                    % not match the non-zero lattice1 cell; add this neighbor to the inCand
                    % list; Do the same for the conjugates; Accordingly adjust the old bonds from lattice2; 

                    mem    = find(ismember(cellNamesLattice2,ix));
                    mem    = mem(cellNamesLattice2(mem)~=0); 


                    if cellNamesLattice2(setdiff(1:2,mem)) ~= 0
                        outLattice.bonds(bondReplaceName,swap(mem)+2) = cellNamesLattice2(setdiff(1:2,mem))+ncells;
                        % conjugate
                        outLattice.bonds(conjBondReplaceName,mem+2)   = cellNamesLattice2(setdiff(1:2,mem))+ncells;
                    end
                    % update the list of incandidates; add the new neighbor
                    % of cell1 is not matched with it, and needs to stay in
                    % the outLattice after pruning out the excess cells; 
                    inCand = unique([inCand;cellNamesLattice2(setdiff(1:2,mem))]);
                end

                % replace the bond names in lattice2 with the new bond name; 
                if cellNamesLattice2(1) ~= 0 

                    temp = outLattice.cells{cellNamesLattice2(1)+ncells};
                    temp(temp == bondOutName+nbonds)     = bondReplaceName;
                    temp(temp == conjBondOutName+nbonds) = conjBondReplaceName;
                    outLattice.cells{cellNamesLattice2(1)+ncells} = temp;
                end

                if cellNamesLattice2(2) ~= 0 

                    temp = outLattice.cells{cellNamesLattice2(2)+ncells};
                    temp(temp == bondOutName+nbonds)     = bondReplaceName;
                    temp(temp == conjBondOutName+nbonds) = conjBondReplaceName;
                    outLattice.cells{cellNamesLattice2(2)+ncells} = temp;
                end

       %         end
            end
            

            
            % update matched cells in all the bonds; 
            for i = 1 : length(ix)
                if ix(i)>0 

                    if dontMatchCells(i) == 0
                        % if the bonds are not matched, just replace the
                        % lattice2 cell name with lattice1 cell name (in
                        % index(i), which is the boundary of lattice1);
                        temp = outLattice.bonds(nbonds+1:end,3);
                        temp(temp == ix(i)+ncells) = index(i);
                        outLattice.bonds(nbonds+1:end,3) = temp;

                        temp = outLattice.bonds(nbonds+1:end,4);
                        temp(temp == ix(i)+ncells) = index(i);
                        outLattice.bonds(nbonds+1:end,4) = temp;
                    else

                        if any(outLattice.bonds(outLattice.cells{index(i)},:) - outLattice.bonds(outLattice.cells{ix(i)+ncells},:) ~= 0)

                            % generally, they should have a neighbor that is not an
                            % inCand any more; So make them correct 0 cell
                            % by removing that neighbor.

                            currBond2     = lattice2.cells{ix(i)};
                            conjCurrBond2 = getConjCell(lattice2,currBond2);

                            % since this cell will not be in the system anymore, make it a
                            % zero cell;
 
                            outLattice.bonds(currBond2+nbonds,3) = 0;
                            outLattice.bonds(conjCurrBond2+nbonds,4) = 0;
                            
                            inCand = setdiff(inCand,ix(i));
                        end
                    end
                end
            end
            
            % set the neighboring cells of not inCand members to zero;
            notInCand = setdiff(1:length(lattice2.cells),inCand);
            temp = outLattice.bonds(:,3);
            temp(ismember(temp,notInCand+ncells)) = 0;
            outLattice.bonds(:,3) = temp; 

            temp = outLattice.bonds(:,4);
            temp(ismember(temp,notInCand+ncells)) = 0;
            outLattice.bonds(:,4) = temp; 
                    
            
            finOutLattice = this.indexedSubLattice(outLattice,[1:ncells,inCand'+ncells]);

        
            % finally some consistency check. This should only be
            % excecuted during debuggin;
            disp('skipping consistency check of inOutLattice! Make this debugging only (line 705)');
%             bonds = finOutLattice.bonds;
%             conjBondIndex  = getConjCell(finOutLattice,1:length(bonds));
%             if length(unique(conjBondIndex)) ~= length(conjBondIndex)
%                 disp('Not every bond in finOutLattice has a unique conjugate');
%             else
%                 disp('All bonds in finOutLattice have a unique conjugate');
%             end
            
        end
        
        
        % ------------------------------------------------------
        %  resort cells with normal pointing outwards;
        % ------------------------------------------------------
        
        function outLattice = normalOutwardOrientation(this,inLattice)
            % Orient cells in a lattice, such that the normal points outwards;      

            % Assume that the cells have previously been sorted, so all it
            % comes down to is pick a direction; 

           % hold on 
            outLattice = inLattice;
            for i = 1 : length(inLattice.cells)

                currCell        = inLattice.cells{i};
                CoM3            = inLattice.CoM3(i,:);
                currBonds       = inLattice.bonds(currCell,:);
                currVerts       = inLattice.verts3(currBonds(:,1),:);
                % centre the current vertices; 
                centerCurrVerts = [currVerts(:,1)-CoM3(1),currVerts(:,2)-CoM3(2),currVerts(:,3)-CoM3(3)];

                normal = -this.crossProduct3(centerCurrVerts(2,:),centerCurrVerts(1,:));
                normal = normal/norm(normal);

               %	plot3(currVerts([1:end,1],3),currVerts([1:end,1],1),currVerts([1:end,1],2),'-');

               %if norm(normal+CoM3) > norm(CoM3)
                    %%disp('correct orientation');
                   % quiver3(CoM3(3),CoM3(1),CoM3(2),normal(3),normal(1),normal(2),100);

                if norm(normal+CoM3) < norm(CoM3)
                    %disp('Correcting, ');

                    conjCurrCell = getConjCell(inLattice,currCell);
                    conjCurrCell = conjCurrCell(end:-1:1);

                    inLattice.cells{i} = conjCurrCell;

                end


            end
        end

        function vOut = crossProduct3(this,vIn1,vIn2)

            vOut    = 0*vIn1;
            vOut(1) = vIn1(2)*vIn2(3) - vIn1(3)*vIn2(2);
            vOut(2) = vIn1(3)*vIn2(1) - vIn1(1)*vIn2(3);
            vOut(3) = vIn1(1)*vIn2(2) - vIn1(2)*vIn2(1);
        end
        
        
        
        % ------------------------------------------------------
        %  find bonds to be excluded;
        % ------------------------------------------------------
        
        function bondsExcluded1 = findExclusionBonds(this,bondsExcluded1,currBondsCell1,lattice1,verts1,index,RunNumber)
    
            % only take the bonds in the cells that are not excluded yet;
            % Like this we prevent that a cell, which was already shortened
            % falsely removes another bond.
            currCell1  = currBondsCell1(~bondsExcluded1(currBondsCell1));
           
            
            % get the index of the conjugate bonds of currCell1 & 2, also 
            % put out the bond result;  
            [conjCurrCell1,~,currBonds1] = getConjCell(lattice1,currCell1);
           
            
            % sort the bonds in the currCell1 by bond length; While the
            % shortest seems to be the first candidate, it still may
            % neighbor a cell that was already changed; This is
            % unfavorable. Take the next bigger bond; 
            bl1 = sum((verts1(currBonds1(:,1),:)-verts1(currBonds1(:,2),:)).^2,2);
            [~,sortedIndex] = sort(bl1);
            flag    = 0;
            counter = 1;
            
            
            while  ( flag == 0 ) && ( counter <= length(currCell1) )
                
                
                % test if the currCell1(sortedIndex(counter)) and
                % conjCurrCell1(sortedIndex(counter)) is not in a previously
                % computed cell; 
                currIndex = index(1:RunNumber-1);
                outCandidates = cat(2,lattice1.cells{currIndex(currIndex~=0)}); 
                
                % remove the excluded ones;
                outCandidates = outCandidates(~bondsExcluded1(outCandidates));
                if all(~ismember([currCell1(sortedIndex(counter)),conjCurrCell1(sortedIndex(counter))],outCandidates))
                     % neither the bond currCell1(sortedIndex(counter)) nor its conjugate are member of any of the outcands; 
                     % Exclude them now and exit the loop;
                     
                     bondsExcluded1(currCell1(sortedIndex(counter)))     = 1;
                     bondsExcluded1(conjCurrCell1(sortedIndex(counter))) = 1;
                     
                     % exit the loop;
                     flag = 1; 
                else
                    counter = counter+1;
                end
                if counter == length(currCell1) && flag == 0
                    disp('Went out possible bonds to be excluded. Not matching cell!');
                end
            end
        end
        
        
        % ------------------------------------------------------
        %  exclude bonds
        % ------------------------------------------------------
        
        function lattice1Out = excludeBonds(this,lattice1,bondsExcluded1)

            % bondsExcluded1 are the identified bonds which have to
            % be taken out of lattice1, to match coordination
            % numbers; 

            bondsOut1   = find(bondsExcluded1 == 1);

            lattice1Out = lattice1;

            for i = 1 : length(bondsOut1)

                % fid the conjugate bond of bondsOut1;
                [conjBondIndex1,conjCurrBonds1,currBonds1] = getConjCell(lattice1,bondsOut1(i));

                conjI = find(bondsOut1 == conjBondIndex1);

                % bondsOut1(i) and bondsOut1(conjI) have now to be replaced;

                vertices1     = lattice1.bonds(bondsOut1(i),1:2);
               % conjVertices1 = lattice1.bonds(bondsOut1(conjI),1:2);

                % find all bonds that share the vertex vertices1(1);
                v1Bonds = setdiff(unique([find(lattice1.bonds(:,1) == vertices1(1));find(lattice1.bonds(:,2) == vertices1(1))]),...
                    [bondsOut1(i),bondsOut1(conjI)]);
                v2Bonds = setdiff(unique([find(lattice1.bonds(:,1) == vertices1(2));find(lattice1.bonds(:,2) == vertices1(2))]),...
                    [bondsOut1(i),bondsOut1(conjI)]);

                % compute the lengths of the bonds;
                v1BondsLengths = sum( (lattice1.verts3(lattice1.bonds(v1Bonds,1),:)-lattice1.verts3(lattice1.bonds(v1Bonds,2),:) ).^2 ,2);
                v2BondsLengths = sum( (lattice1.verts3(lattice1.bonds(v2Bonds,1),:)-lattice1.verts3(lattice1.bonds(v2Bonds,2),:) ).^2 ,2);

                if mean(v1BondsLengths) < mean(v2BondsLengths)

                    % the average lengths of the bonds attached to vertices1(1) is
                    % shorter, so keep vertices1(2); This way the on average shorter bonds get longer;  
                    vertKeep = vertices1(2);
                    vertOut  = vertices1(1);

                    %temp = lattice1.bonds(v1Bonds,1:2); 
                    temp = lattice1Out.bonds(v1Bonds,1:2);
                    temp(temp== vertOut) = vertKeep;
                    lattice1Out.bonds(v1Bonds,1:2) = temp;
                else

                    % the average lengths of the bonds attached to vertices1(2) is
                    % shorter, so keep vertices1(1); This way the on average shorter bonds get longer;  
                    vertKeep = vertices1(1);
                    vertOut  = vertices1(2);

                    %temp = lattice1.bonds(v2Bonds,1:2); 
                    temp = lattice1Out.bonds(v2Bonds,1:2);
                    temp(temp== vertOut) = vertKeep;
                    lattice1Out.bonds(v2Bonds,1:2) = temp;
                end

                % don't forget to remove the bond from its neighboring cells
                correctCells = unique([currBonds1(3:4),conjCurrBonds1(3:4)]);
                correctCells = correctCells(correctCells~=0);
                for k = 1 : length(correctCells)

                    temp = lattice1Out.cells{correctCells(k)};
                    temp = temp(temp~=bondsOut1(i));
                    temp = temp(temp~=bondsOut1(conjI));
                    lattice1Out.cells{correctCells(k)} = temp;
                end

                lattice1 = lattice1Out;
            end



        end 


        
        % ------------------------------------------------------
        %  find sublattice of a lattice depending on index of cells;
        % ------------------------------------------------------
        
         function outLattice = indexedSubLattice(this,inLattice,cellIndex)


                % Use only indexed cells, unique indexed bonds and vertices; 
                bonds  = inLattice.bonds;
                cells  = inLattice.cells;
                verts  = inLattice.verts; 
                verts3 = inLattice.verts3;

                % we need the zero cell, since we don't keep it in our original lattice
                % structure, yet want the conjugate of each bond to be kept; Use the 
                % conjZeroCell to determine which of the zeroCell bonds should be used;
                % As some may be excluded;
                zeroCell     = find(inLattice.bonds(:,3) == 0)';
                zeroCell2    = find(inLattice.bonds(:,4) == 0)';
                zeroCell     = setdiff(zeroCell,zeroCell2);
                conjZeroCell = getConjCell(inLattice,zeroCell);

                % some of these bonds are duplicated; only use first occurence;
                [conjZeroCell,unIx] = unique(conjZeroCell,'first');
                zeroCell = zeroCell(unIx);

                % 1) Map current cell names in the bonds to new cellIndex;  
                newBonds = bonds;
                for i = 1 : length(cellIndex)

                    temp = bonds(:,3)==cellIndex(i);
                    newBonds(temp,3) = i;

                    temp = bonds(:,4)==cellIndex(i);
                    newBonds(temp,4) = i;
                end

                newBonds(~ismember(bonds(:,3),cellIndex),3) = 0;
                newBonds(~ismember(bonds(:,4),cellIndex),4) = 0;

               % inLattice.bonds = newBonds;

                % 2) Find bonds addressed by indexed cells; Take unique ones only! 
                usedBonds    = cat(2,cells{cellIndex});
                useZeroCell  = zeroCell(ismember(conjZeroCell,usedBonds));


                usedBonds = [usedBonds,useZeroCell];
                [usedBonds2,ixUsedBonds] = unique(usedBonds);

            %     unUsedBonds = setdiff(1:length(inLattice.bonds),usedBonds);
            %     
            %     conjBondIndex  = getConjCell(inLattice,usedBonds);
            %     conjBondIndex  = setdiff(conjBondIndex,usedBonds); 
            %     bonds(conjBondIndex,:)
            %     disp('these are conjugate bonds, but not used in the cells')
            %     pause
            %     
            %     
                bonds = newBonds(usedBonds2,:);

                % 3) replace all indexed cell entries by address of unique bonds;
                for i = 1 : length(cellIndex)
                    temp = cells{cellIndex(i)}; 
                    for k = 1 : length(temp)
                        temp(k) = find(usedBonds2 == temp(k));
                    end
                    cells{cellIndex(i)} = temp;
                end
                cells = cells(cellIndex);



                % 4) in the unique bonds, find the addressed vertices and keep only those. 
                usedVerts = [bonds(:,1);bonds(:,2)];
                [usedVerts2,~,C] = unique(usedVerts);
                verts      = verts(usedVerts2,:);
                verts3     = verts3(usedVerts2,:);

                % 5) replace indices to vertices in unique bonds by mapping; 
                bonds(:,1) = C(1:size(bonds,1));
                bonds(:,2) = C((size(bonds,1)+1):2*size(bonds,1));


                % theses bonds may still have cell entries that are not currently
                % in the cell index; Set them to the zero cell; Make sure they all
                % have a conjugate bond, to have a prober lattmin structure.

                %zeroBonds = find(~ismember(bonds(:,4),cellIndex));
        %         zeroBonds = find(bonds(:,4)> length(cellIndex));
        %         for i = 1 : length(zeroBonds)
        %             
        %             j = find(sum([bonds(:,2)-bonds(zeroBonds(i),1),bonds(:,1)-bonds(zeroBonds(i),2),...
        %                 bonds(:,4)-bonds(zeroBonds(i),3),bonds(:,3)-bonds(zeroBonds(i),4)].^2,2)==0);
        %             if isempty(j)
        %                 bonds(zeroBonds(i),4) = 0;
        %                 % generat the conjugate bond, as it is missing;
        %                 bonds(end+1,:) = [bonds(zeroBonds(i),2),bonds(zeroBonds(i),1),bonds(zeroBonds(i),4),bonds(zeroBonds(i),3)];
        %             else
        %                 bonds(zeroBonds(i),4) = 0;
        %                 bonds(j,3) = 0;
        %             end
        %         end
                zeroBonds = find(bonds(:,4) == 0);

                for i = 1 : length(zeroBonds)

                    j = find(sum([bonds(:,2)-bonds(zeroBonds(i),1),bonds(:,1)-bonds(zeroBonds(i),2),...
                        bonds(:,4)-bonds(zeroBonds(i),3),bonds(:,3)-bonds(zeroBonds(i),4)].^2,2)==0);
                    if isempty(j)

                        % generat the conjugate bond, as it is missing;
                        bonds(end+1,:) = [bonds(zeroBonds(i),2),bonds(zeroBonds(i),1),bonds(zeroBonds(i),4),bonds(zeroBonds(i),3)];
                    end
                end


                disp('Make proper sorting message of cell a debugging message!')
                %Check that cell is properly sorted, 
        %         for i = 1 : length(cells)
        % 
        %             temp = cells{i};
        %             tempV = bonds(temp,1:2);
        %             ii = [1:length(temp),1];
        %             for k = 1 : length(tempV)
        %                 if  tempV(ii(k),2)~=tempV(ii(k+1),1)
        %                     tempV
        %                     i
        %                 end
        %             end
        %         end

                % make sure all bonds have a counter bond?


                % 6) loop through all the cells, and sort them so they are oriented
                % outwards; 
                % disp('Cell orientation with normal pointing outwards needs to be implemented');




                outLattice = struct('cells',[],'bonds',[],'verts',[],'verts3',[],...
                        'CoM',[],'CoM3',[],'phi',[]);
                outLattice.cells  = cells;
                outLattice.bonds  = bonds;
                outLattice.verts  = verts;
                outLattice.verts3 = verts3;
                outLattice.CoM    = inLattice.CoM(cellIndex,:);
                outLattice.CoM3   = inLattice.CoM3(cellIndex,:);
                outLattice.phi    = inLattice.phi(cellIndex);

        end

        
        
        % ------------------------------------------------------
        %  find index of boundary cells in a lattice;
        % ------------------------------------------------------
        
        function index = boundaryCells(this,latt)

            index = zeros(size(latt(1).cells));
            for k = 1 : length(latt(1).cells)
                temp = latt(1).cells{k};
                bds  = latt(1).bonds(temp,3:4);
                if any(bds(:) == 0)
                   % cell is on the boundary;
                   index(k) = 1;
                end
            end
            index = find(index);
        end
        
        
        % ------------------------------------------------------
        %  extract subLattice based on index; 
        % ------------------------------------------------------
        
        function subLattice(this,lattice,index)

            % given a Lattmin lattice, and a set of indices of cells, generate a
            % sublattice containing only those cells; 

            verts2 = lattice.verts;
            verts3 = lattice.verts3;
            com    = lattice.CoM;
            com3   = lattice.CoM3;
            phi    = lattice.phi;
            bonds  = lattice.bonds;
            cells  = lattice.cells;

            Cells  = cells(index);
            CoM    = com(index,:);
            CoM3   = com3(index,:);
            phi    = phi(index);
            usedVerts = zeros(length(verts2),1);
            for k = 1 : length(Cells)

                temp = bonds(Cells{k},1);
                usedVerts(temp) = 1;
                Cells{k} = temp;
            end

            unVerts = find(usedVerts == 1);

            for k = 1 : length(Cells)
                temp = Cells{k};
                temp2 = 0*temp;
                for kk = 1 : length(temp)
                    temp2(kk) = find(unVerts == temp(kk)); 
                end

                Cells{k} = temp2';
                if length(unique(Cells{k}))~= length(Cells{k})
                    disp('Wrong Cell Length');
                end
            end
        %  
            g = GLattConversion(Cells, verts2(unVerts,1:3)); 
            g.verts  = verts2(unVerts,:);
            g.verts3 = verts3(unVerts,:);
            g.CoM    = CoM;
            g.CoM3   = CoM3;
            g.phi    = phi;
            this.currLatt = g;
        end
        
        
    end
end
