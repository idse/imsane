function saveFloatTiff(im, filename)
    % SAVEFLOATTIFF Save a floating point matrix to 32 bit tiff
    % imwrite doesn't do this for some reason
    % we convert to single in the process
    %
    % saveFloatTiff(im, filename)
    %
    % im: Image
    % filename: Name of file

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
    
    t = Tiff(filename, 'w'); 
    tagstruct.ImageLength = size(im, 1); 
    tagstruct.ImageWidth = size(im, 2); 
    
    %tagstruct.Compression = Tiff.Compression.None; 
    tagstruct.Compression = Tiff.Compression.Deflate; 
    tagstruct.SampleFormat = Tiff.SampleFormat.IEEEFP; 
    tagstruct.Photometric = Tiff.Photometric.MinIsBlack; 
    tagstruct.BitsPerSample = 32; %info.BitsPerSample; % 32; 
    tagstruct.SamplesPerPixel = 1; %info.SamplesPerPixel; % 1; 
    tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky; 
    t.setTag(tagstruct); 

    t.write(single(im));
    t.close();
end