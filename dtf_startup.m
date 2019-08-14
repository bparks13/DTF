CCC;

rootPath=get_root_path();

addpath(fullfile(rootPath,'Calculate'),...
    fullfile(rootPath,'Helper'),...
    fullfile(rootPath,'Load'),...
    fullfile(rootPath,'Output'));

if exist(fullfile(get_root_path,'Files'),'dir') == 0
    mkdir(fullfile(rootPath,'Files'));
end

CCC;