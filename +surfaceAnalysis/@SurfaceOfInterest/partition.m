function partitionMask = partition(this, time, chartIndices, margin)
    % For some set of charts generate masks to get the non-overlapping
    % parts
    %
    % partitionMask = partition(time, chartIndices, margin)
    %
    % partitionMask:    cell array of binary masks for each of the charts
    %                   specified in chartIndices
    %
    % time:             time point
    % chartIndices:     indices into SOI.atlas.charts 
    % margin:           how much to erode the overlapping masks by before
    %                   removing the intersection 
    %                   (with 0, the first chart will not be cropped and
    %                   overlap will be removed from the remaining chart
    %                   masks)

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

    if nargin == 3
        margin = 0;
    end

    embMasks = {};
    partitionMask = {};

    for i = 1:numel(chartIndices)

        % mask
        embMasks{i} = this.embedding(time).patches{chartIndices(i)}.apply{1} > 0;
        embMasks{i} = imerode(embMasks{i}, strel('disk', margin));

        mask = false(size(embMasks{i}));

        for j = 1:i-1

            chName1 = this.atlas(time).charts{chartIndices(j)}.image.name;
            %['conformal_' num2str(j)];
            chName2 = this.atlas(time).charts{chartIndices(i)}.image.name;
            dm = tmapMask(this, time, embMasks{j}, chName1, chName2);
            mask = mask | dm;
        end

        partitionMask{i} = embMasks{i}.*~mask > 0;
    end
end

function mask2 = tmapMask(SOI, time, mask1, chartName1, chartName2)

    tmap12 = SOI.atlas(time).getTransitionMap(chartName1, chartName2);
    tmap21 = SOI.atlas(time).getTransitionMap(chartName2, chartName1);
    
    tmap12mask = ~isnan(tmap12.apply{1});
   
    % determine the boundaries of the mask intersected with the tmap 
    [b,~,~,A] = bwboundaries(tmap12mask.*mask1);
    
    % this is some stupid stuff to get the inner boundary approximately
    % lying inside the region, rather than just outside
    b2 = bwboundaries(imerode(tmap12mask.*mask1,strel('square',3)));
    
    nBdries = numel(b);

    outerBdry = 1;
    for i = 1:nBdries

        enclosing_boundary = find(A(i,:));

        if isempty(enclosing_boundary)
            outerBdry = i;
        else
            b{i} = b2{i};
        end
    end
    
    % find linear indices of boundary
    bdryInd = {};
    for i = 1:nBdries
        bdryInd{i} = sub2ind(size(tmap12mask), b{i}(:,1), b{i}(:,2));
    end
    
    % map the boundary to the other chart using the tmap
    mappedB = {};
    for i = 1:nBdries
        u1 = tmap12.apply{1}(bdryInd{i}) - tmap21.domain.boundary{1}(1);
        u2 = tmap12.apply{2}(bdryInd{i}) - tmap21.domain.boundary{2}(1);
        mappedB{i} = [u1 u2];
    end

    % determine which is the outer boundary in the 2nd chart
    outBdry = 1;
    for bi = 2:nBdries
        if all(inpolygon(mappedB{1}(:,1), mappedB{1}(:,2),...
                         mappedB{bi}(:,1),mappedB{bi}(:,2)));
            outBdry = bi;
        end
    end

    tmap21mask = ~isnan(tmap21.apply{1});

    mask2 = poly2mask(mappedB{outBdry}(:,1),mappedB{outBdry}(:,2),...
                        size(tmap21mask,1), size(tmap21mask,2));
                                    
    for bi = setdiff(1:nBdries,outBdry)
        mask2 = mask2 - poly2mask(mappedB{bi}(:,1),mappedB{bi}(:,2),...
                        size(tmap21mask,1), size(tmap21mask,2));
    end
end