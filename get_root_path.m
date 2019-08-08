function rootPath=get_root_path()
%% rootPath=get_root_path()
%
%  Finds the location of the current file, specifically the absolute path, and returns the
%  folder
%
%   Inputs: None
%
%   Outputs:
%    - rootPath: String containing the absolute path to the folder where this function is
%       running. This should be found in the root directory of the DTF Git
%

currPath=mfilename('fullpath');
rootPath=fileparts(currPath);

end