    


%Creates data folder struct
if not(isfolder("Data"))
    fprintf('\n->>>FirstTimeSetup\n->>>Creating "Data" folder')
    mkdir("Data")
    
end

%Creates Analysis folder struct
if not(isfolder("Analysis"))
    fprintf('\n->>>Creating "Analysis" folder\n')
    mkdir("Analysis")
    
end

%Adds both new folders to path so matlab can see them
addpath("Data\",'-end')
addpath("Analysis\",'-end')
