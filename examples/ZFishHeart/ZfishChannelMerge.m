[scriptPath,~,ext] = fileparts(matlab.desktop.editor.getActiveFilename);

dataDir = '/Users/idse/Dropbox/flows/flows_shared/data/zebrafish/huisken/video9/';
projectDir = dataDir;

xp = project.Experiment(projectDir, dataDir);
cd(dataDir)

fileMeta                 = struct();
fileMeta.dataDir         = dataDir;
fileMeta.filenameFormat  = 'VSilent_T%d_G.tif';%'Time%06d_bin2.tif'; % for full data sample use Time000000.tif
fileMeta.timePoints      = 1:129; % for full data sample use 0;
fileMeta.stackResolution = [.45 .45 1];%[.5 .5 .5]; 
fileMeta.swapZT          = 0; % for full data sample use 1;

expMeta                  = struct();
expMeta.channelsUsed     = 1;
expMeta.channelColor     = 1;
expMeta.description      = 'Beating Zebrafish Heart';
expMeta.dynamicSurface   = 1;
expMeta.jitterCorrection = 0; % 1: Correct for sample translation
expMeta.fitTime          = 60; 
expMeta.detectorType     = 'surfaceDetection.fastCylinderDetector';%
expMeta.fitterType       = 'surfaceFitting.meshWrapper'; 

xp.setFileMeta(fileMeta);
xp.setExpMeta(expMeta);


%% load, merge

channels = [1 2];

for t = fileMeta.timePoints

    mergedFName = sprintf('VSilent_T%d_merged.ome.tif', t);
    
    if ~exist(mergedFName, 'file')
        channelData = {};

        for c = channels

            if c == 1
                fileMeta.filenameFormat  = 'VSilent_T%d_G.tif';
            else
                fileMeta.filenameFormat  = 'VSilent_T%d_R.tif';
            end

            xp = project.Experiment(projectDir, dataDir);
            xp.setFileMeta(fileMeta);
            xp.setExpMeta(expMeta);
            xp.initNew();
            xp.loadTime(t);
            channelData{c} = xp.stack.image.apply{1};
        end

        stack4D = cat(4, channelData{:});
        bfsave(stack4D, mergedFName, 'dimensionOrder', 'XYZCT');
    end
end
