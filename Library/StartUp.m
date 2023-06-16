clc
%% Preamble
%%Backend setup hidden from user
warning('off','MATLAB:MKDIR:DirectoryExists')
%%
printBreak(3)
cd
ver = '1.002';%build ID
fprintf('Initializing QdFinder: Version %s on %s\n',ver,string(computer));
if ispc
    !git -v
else
    !git --version
end

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
if not(isfile("Library\sortingP_idxGuide_1emitter_backgrd_dc.mat")) && not(isfile("Library/sortingP_idxGuide_1emitter_backgrd_dc.mat"))
    fprintf('Get single emitter sort key\n')
    makeSortingMatrix_1emitter_backgrd_dc
    fprintf('Found key: makeSortingMatrix_1emitter_backgrd_dc\n')
else
    fprintf('Found key: makeSortingMatrix_1emitter_backgrd_dc\n')
end
if not(isfile("Library\sortingP_idxGuide_2emitters_backgrd_dc.mat")) && not(isfile("Library/sortingP_idxGuide_2emitters_backgrd_dc.mat"))
    fprintf('Get multi emitter sort key\n')
    makeSortingMatrix_2emitters_backgrd_dc
    fprintf('Found key: makeSortingMatrix_2emitters_backgrd_dc\n')
else
    fprintf('Found key: makeSortingMatrix_2emitters_backgrd_dc\n')
end
printBreak



%% Start Up
%Set up changes for OS. Set up supports Linux and Windows. Might support
%MacOS. I dont have an .iso to check with.
%Adds both folders to path so matlab can see them
fprintf('Adding resources to path: \n')

if isunix
    addpath(genpath("Data/"))
elseif ispc
    addpath(genpath("Data\"))
elseif ismac
    addpath(genpath("Data/"))
end
fprintf('..Data..\n')

if isunix
    addpath(genpath("Analysis/"))
elseif ispc
    addpath(genpath("Analysis\"))
elseif ismac
    addpath(genpath("Analysis/"))
end
fprintf('..Analysis..\n')

if not(isfile("Library\DirPath.mat"))
DirPath = cd;
fprintf('Saving directory information\n Main directory is at:\n%s ..\n',DirPath)
if ispc
save("Library\DirPath.mat","DirPath")
elseif isunix || ismac
save("Library/DirPath.mat","DirPath")
end
end
printLine


%% End, enter program
fprintf('Done\n\n\n')
fprintf('Execute "Main()" to begin analysis\n')
printBreak
printLine
