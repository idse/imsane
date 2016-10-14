function convolved = myconvn(stack, filter)
    % convn with periodic bcs
    %
    % myconvn(stack,filter)
    
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
    kersize = size(filter);
    padsize = ceil(kersize / 2) - 1;
    
    if padsize(1) > size(stack, 1) || padsize(2) > size(stack, 2) ||...
        ( numel(padsize) > 2 && padsize(3) > size(stack,3) )
        
        error('the kernel is larger than the stack!');
    end
    
    dim = 1;
    if padsize(dim) > 0

        padding = flipdim(stack(end-padsize(dim) + 1:end, :, :), dim);
        stack = cat(dim, stack, padding);
        padding = flipdim(stack(1: padsize(dim), :, :), dim);
        stack = cat(dim, padding, stack);
    end
    
    dim = 2;
    if padsize(dim) > 0 

        padding = flipdim(stack(:, end-padsize(dim) + 1:end, :), dim);
        stack = cat(dim, stack, padding);
        padding = flipdim(stack(:, 1: padsize(dim), :), dim);
        stack = cat(dim, padding, stack);
    end
    
    dim = 3;
    if numel(padsize) > 2 && padsize(dim) > 0

        padding = flipdim(stack(:, :, end-padsize(dim) + 1:end), dim);
        stack = cat(dim, stack, padding);
        padding = flipdim(stack(:, :, 1: padsize(dim)), dim);
        stack = cat(dim, padding, stack);
    end

    convolved = convn(stack, filter, 'valid');
    
end

