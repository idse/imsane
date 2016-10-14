function out = xpathText(node, path)
    % XPATHTEXT Shorthand for loading text from node in xml file
    %
    % xpathText(node, path)
    % 
    % if path doesn't exist, return empty string
    
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
    
    import javax.xml.xpath.*;
    factory = XPathFactory.newInstance;
    xpath = factory.newXPath;

    expression = xpath.compile(path);
    
    evalExp = expression.evaluate(node, XPathConstants.NODE);
    
    % if path doesn't exist, return empty string
    if isempty(evalExp)
        out = char([]);
    else
        out = char(evalExp.getTextContent);
    end
end
