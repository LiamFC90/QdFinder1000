clc
%% Preamble
ver = '1.001';%build ID
fprintf('Initializing QdFinder: Version %s on %s\n',ver,string(computer));
!git -v

%% First time only

%Creates data folder struct
if not(isfolder("Data"))
    mkdir("Data")
end

%Creates Analysis folder struct
if not(isfolder("Analysis"))
    mkdir("Analysis")
end

%Creates sorting keys
if not(isfile("Library\sortingP_idxGuide_1emitter_backgrd_dc.mat")) || not(isfile("Library/sortingP_idxGuide_1emitter_backgrd_dc.mat"))
    makeSortingMatrix_1emitter_backgrd_dc
else
    fprintf('Found key: makeSortingMatrix_1emitter_backgrd_dc\n')
end
if not(isfile("Library\sortingP_idxGuide_2emitters_backgrd_dc.mat")) || not(isfile("Library/sortingP_idxGuide_2emitters_backgrd_dc.mat"))
    makeSortingMatrix_2emitters_backgrd_dc
else
    fprintf('Found key: makeSortingMatrix_2emitters_backgrd_dc\n')
end




%% Start Up
%Set up changes for OS. Set up supports Linux and Windows. No MacOS support
%Adds both folders to path so matlab can see them
fprintf('Adding resources to path: \n')

if isunix
    addpath(genpath("Data/"))
elseif ispc
    addpath(genpath("Data\"))
elseif ismac
    addpath(genpath("Data/"))
end
fprintf('...Data...\n')

if isunix
    addpath(genpath("Analysis/"))
elseif ispc
    addpath(genpath("Analysis\"))
elseif ismac
    addpath(genpath("Analysis/"))
end
fprintf('...Analysis...\n')

if not(isfile("Library\DirPath.mat"))
DirPath = cd;
fprintf('Saving directory information\n Main directory is at:\n%s\n',DirPath)
if ispc
save("Library\DirPath.mat","DirPath")
elseif isunix || ismac
save("Library/DirPath.mat","DirPath")
end
end

%% End, enter program
fprintf('Done\n\n\n')
fprintf('Execute "Main()" to begin analysis\n')
