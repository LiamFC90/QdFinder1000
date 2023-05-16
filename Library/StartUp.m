clc
%% Preamble
ver = '1.000';%build ID
fprintf('Initializing QdFinder: Version %s on %s\n',ver,string(computer));

%% First time only

%Creates data folder struct
if not(isfolder("Data"))
    mkdir("Data")
end

%Creates Analysis folder struct
if not(isfolder("Analysis"))
    mkdir("Analysis")
end

%% Start Up
%Set up changes for OS. Set up supports Linux and Windows. No MacOS support
%Adds both folders to path so matlab can see them
fprintf('Adding resources to path: \n')

if isunix
    addpath(genpath("Data/"))
elseif ispc
    addpath(genpath("Data\"))
end
fprintf('...Data...\n')

if isunix
    addpath(genpath("Analysis/"))
elseif ispc
    adddpath(genpath("Analysis\"))
end
fprintf('...Analysis...\n')

%% End, enter program
fprintf('Done\n\n\n')
fprintf('Execute "Main()" to begin analysis\n')
