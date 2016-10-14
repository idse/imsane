function debugMsg(level, msg)
    % Set the debugging level 
    % 
    % debugMsg(level, msg)
    %
    % level 1: function names
    % level 2: function details
    % level 3: map call
    % msg: string containing debugging message
    
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
    
    debugLevel = getpref('ImSAnE', 'msgLevel');
    
    if level <= debugLevel
        
        for i = 2:level
            indent = '  ';
            msg = [indent msg];
        end
        fprintf(msg);
    end 
end