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

%%
%-----------------------------------------
% set up ImSAnE path
%-----------------------------------------

% entries of matlab search path
pathentries = regexp(path, pathsep, 'split');

% find ones containing ImSAnE and remove from path
disp('removing old ImSAnE entries from path');
for i = 1:length(pathentries)
    if ~isempty(regexp(pathentries{i}, 'ImSAnE','once'))...
       || ~isempty(regexp(pathentries{i}, 'imsane','once'))
        rmpath(pathentries{i});
    end
end

% path of current script is new ImSAnE directory
[ImSAnEpath,~,ext] = fileparts(matlab.desktop.editor.getActiveFilename);
addpath(genpath(ImSAnEpath));
disp('added ImSAnE directory containing setup to path');

% storing the path in the startup directory
upath = userpath;
upath = upath(1:end-1);
savepath(fullfile(upath, 'pathdef.m'));

%  add path to settings for easy retrieval
setpref('ImSAnE', 'path', ImSAnEpath);

%%
%-----------------------------------------
% settings
%-----------------------------------------

% detail of the level of output messages 
% admittedly not very well implemented
% 1: function names
% 2: function details
% 3: map call
msgLevel = 2;
setpref('ImSAnE', 'msgLevel', msgLevel);

% store figures to check quality of surface fits
fitQualityPlots =  1;
setpref('ImSAnE', 'fitQualityPlots', fitQualityPlots);

%%
%-------------------------------------------------------
% compile mex code for mesh distance by Gabriel Peyre
%-------------------------------------------------------

cd(fullfile(ImSAnEpath, 'external', 'fast_marching'));
compile_mex

%% 
%-------------------------------------------------------
% For Linux and meshlab usage: Set environment variable
%-------------------------------------------------------

% In case using sytem command in combination with meshlab faild due to a 
% symbol lookup error: This may be caused by matlab system command using 
% the wrong libraries or doesnt find libcommon.so.1.  
% In this case, set the library path variable in matlab to point to the
% right library. In Ubuntu this is done by copying the line below into the command window: 

pathToLibCommon = '/usr/lib/x86_64-linux-gnu';
setenv('LD_LIBRARY_PATH',sprintf(genpath(pathToLibCommon),getenv('LD_LIBRARY_PATH')));

