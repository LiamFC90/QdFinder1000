%closes and clears vars and figures before closing
clearvars
clc
close all 
closeNoPrompt(matlab.desktop.editor.getAll)

rmpath("Analysis/")
rmpath("Data/")
fprintf('\nRemoved all paths\n')