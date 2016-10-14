function preprocessData(this, surfaceDetectionChannel, saveSmall, bbox)
% PREPROCESSDATA, crop/scale raw data, make small version for Ilastik
%
% REQUIRES BIOFORMATS
% 
% preprocessData(surfaceDetectionChannel, saveSmall, bbox)
%
% bbox: bounding box [xmin xmax ymin ymax zmin zmax]


% determine ROI for first and last time point

channelColor = this.expMeta.channelColor;
channelsUsed = this.expMeta.channelsUsed;

expMeta = this.expMeta;
expMeta.channelsUsed = surfaceDetectionChannel;
expMeta.channelColor = 1;

this.setExpMeta(expMeta);

if nargin == 3

    % load first time point and set mask
    this.loadTime(this.fileMeta.timePoints(1));
    this.rescaleStackToUnitAspect();

    yval = round(this.stack.imageSize(2)/2);

    close all;

    happy = 'No';
    while(strcmp(happy,'No'))

        this.stack.generateProjectionMask([1 2 3])

        imshow(this.stack.getSlice('y', yval), []);

        happy = questdlg('Happy with mask?', 'Happy check', 'Yes', 'No', 'No');
    end
    close
    drawnow
    
    bbi = this.stack.getMaskBoundingBox;

    % load last time point and set mask
    this.loadTime(this.fileMeta.timePoints(end));
    this.rescaleStackToUnitAspect();

    happy = 'No';
    while(strcmp(happy,'No'))

        this.stack.generateProjectionMask([1 2 3])

        imshow(this.stack.getSlice('y', yval), []);

        happy = questdlg('Happy with mask?', 'Happy check', 'Yes', 'No', 'No');
    end
    close
    drawnow

    bbf = this.stack.getMaskBoundingBox;

    % make a joint bounding box

    bbtot = bbi;
    bbtot([1 3 5]) = min(bbi([1 3 5]),bbf([1 3 5]));
    bbtot([2 4 6]) = max(bbi([2 4 6]),bbf([2 4 6]));
else
    bbtot = bbox;
end

save(fullfile(this.fileMeta.dataDir,'ROIbbox.mat'),'bbtot');

debugMsg(1,['Bounding box: ' num2str(bbtot) '\n']);

%--------------------------------
% BIG LOOP processing all timepts
%--------------------------------

expMeta.channelsUsed = channelsUsed;
expMeta.channelColor = channelColor;
this.setExpMeta(expMeta);

cropdir = fullfile(this.fileMeta.dataDir,'cropscale');
if ~exist(cropdir,'dir')
    mkdir(cropdir);
end

smalldir = fullfile(this.fileMeta.dataDir,'small');
if ~exist(smalldir,'dir')
    mkdir(smalldir);
end

delete(fullfile(smalldir,'*'));
delete(fullfile(cropdir,'*'));

ssfactor = 4;
average = true;

for ti = 1:numel(this.fileMeta.timePoints)
    
    t = this.fileMeta.timePoints(ti);
    
    this.loadTime(t);
    this.rescaleStackToUnitAspect();
    
    % crop
    croppedStack = this.stack.cropStack(bbtot);
    %imshow(croppedStack.getSlice('x',500));

    [~,barefn,~] = fileparts(sprintf(this.fileMeta.filenameFormat,t));
    outputPath = fullfile(cropdir, [barefn '_cropscale.ome.tif']);
    I = cat(4,croppedStack.image.apply{:});
    bfsave(I, outputPath);
    
    if saveSmall
        
        % small
        dch = croppedStack.image.apply{surfaceDetectionChannel};
        croppedOneChannel = surfaceDetection.Stack(dch);
        cropRedStack = croppedOneChannel.reduceStack(ssfactor, average);
        % imshow(cropRedStack.getSlice('x',125))

        [~,barefn,~] = fileparts(sprintf(this.fileMeta.filenameFormat,t));
        outputPath = fullfile(smalldir, [barefn '_small.ome.tif']);
        I = cat(4,cropRedStack.image.apply{:});
        bfsave(I, outputPath);
    end
end

end
