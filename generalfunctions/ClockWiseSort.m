function [ClockSorted,Index] = ClockWiseSort(NewVertices)

    NewCell = 1 : size(NewVertices,1);   
    NewCellLength = length(NewCell);
    %   first we need to find the center of mass of the vertices
    if NewCellLength == 1
        Center = NewVertices(NewCell,:);
    else
        Center = sum(NewVertices(NewCell,:))/NewCellLength;
    end
    %   Now take the NewVertices, and subtract the Center from each entry
    MovedVertices = NewVertices(NewCell,:);
    MovedVertices(:,2) = MovedVertices(:,2) - Center(2);
    MovedVertices(:,1) = MovedVertices(:,1) - Center(1);
    [~,Index] = sort(cart2pol(MovedVertices(:,1),MovedVertices(:,2))); 
    % Now we actually sort counter clock, just invert.
    Index = Index(length(Index):-1:1);
    ClockSorted = NewVertices(NewCell(Index),:); 
end