% EXAMPLE HOW TO LOAD SAVED SOI

% define the path to the SOI
ImSAnEpath = getpref('ImSAnE', 'path');
SOIdir = fullfile(ImSAnEpath, 'examples', 'wingDisc', 'discProperApicalSOI');
%SOIdir = uigetdir;

% load the SOI
SOI = surfaceAnalysis.SurfaceOfInterest(SOIdir);

%%
% display surface data
t = 1;
imshow(SOI.data(t).patches{1}.apply{1}, [])