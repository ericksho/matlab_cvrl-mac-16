%DIR2 List directory.
%   DIR2 directory_name lists the files in a directory. Pathnames and
%   wildcards may be used.  For example, DIR *.m lists all program files
%   in the current directory.
%
%   D = DIR2('directory_name') returns the results in an M-by-1
%   structure for MAC OSX with the fields: 
%       name    -- Filename
%       date    -- Modification date
%       bytes   -- Number of bytes allocated to the file
%       isdir   -- 1 if name is a directory and 0 if not
%       datenum -- Modification date as a MATLAB serial date number.
%                  This value is locale-dependent.
%
%   See also WHAT, CD, TYPE, DELETE, LS, RMDIR, MKDIR, DATENUM.

%   Copyright 1984-2010 The MathWorks, Inc.
%   $Revision: 5.12.4.6 $  $Date: 2010/06/15 01:38:35 $
%   Built-in function.

function [files] = dir2(path)
        
    if nargin ~= 1
        path = pwd;
    end
    
    files = dir(path);
    while strcmp(files(1).name,'.') || strcmp(files(1).name,'..') || strcmp(files(1).name,'.DS_Store')
        files(1) = [];
    end
    
    if nargout == 0
        names = [];
        for i = 1:numel(files)
            if i == 1
                names = files(i).name;
            end
            names = sprintf('%s \t %s',names,files(i).name);
        end
        clear files;
        
        disp(names);
    end
end