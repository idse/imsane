function result = fastLinearInterp(inputGrids,image,ss)


image = single(image);


imageInterp = zeros(size(image,1)*ss,size(image,2)*ss,'single');
% divSum      = ones(size(image,1)*ss,size(image,2)*ss,'single');

imageInterp(1:2:end,1:2:end)   = image;
% 
% interpolate colums;
imageInterp(1:2:end,2:2:end)   = imageInterp(1:2:end,2:2:end)+image*.5;
imageInterp(1:2:end,2:2:end-1) = imageInterp(1:2:end,2:2:end-1)+image(:,2:end)*.5;

% interpolate rows; 
imageInterp(2:2:end-2,:)     = (imageInterp(1:2:end-2,:) + imageInterp(3:2:end,:))*.5;

result = imageInterp(1:size(inputGrids{1},1),1:size(inputGrids{1},2));



end