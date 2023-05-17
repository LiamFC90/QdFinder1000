%% closes and clears vars and figures before closing

clearvars
clc
close all 
closeNoPrompt(matlab.desktop.editor.getAll)

%Clears file paths for respective OS
if isunix
    rmpath(genpath("Analysis/"))
    rmpath(genpath("Data/"))
elseif ispc
    rmpath(genpath("Analysis\"))
    rmpath(genpath("Data\"))
elseif ismac
    rmpath(genpath("Analysis/"))
    rmpath(genpath("Data/"))
end
fprintf('\n\nCleanedUp\n\n')