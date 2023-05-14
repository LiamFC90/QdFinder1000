 clc


%% First time only

%Creates data folder struct
if not(isfolder("Data"))
    fprintf('\n->>>FirstTimeSetup\n->>>Creating "Data" folder\n')
    mkdir("Data")
    
end

%Creates Analysis folder struct
if not(isfolder("Analysis"))
    fprintf('->>>Creating "Analysis" folder\n')
    mkdir("Analysis")
    
end

%% Start Up

%Adds both folders to path so matlab can see them
fprintf('Adding resources to path: \n')

addpath(genpath("Data\"))
fprintf('\\Data\\\n')

addpath(genpath("Analysis\"))
fprintf('\\Analysis\\\n')

%% End, enter program
fprintf('Done\n\n\n')
fprintf('Execute "Main()" to begin analysis\n')
