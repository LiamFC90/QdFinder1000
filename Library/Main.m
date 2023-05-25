 %% Begin Main program
   %should remain dominated by user input. Prefer all work done in function
   %workspace. This workspace should hold handles determined by params of
   %data file. These are then passed to function workspaces.
   %Run this script to execute the program
%% Begin Setup
close all
clearvars
% Set output config
%rerun startup to ensure that file has been parsed since git update
StartUp()
fprintf('...Executing Main()...\n')
load(which("DirPath.mat"))
cd(DirPath)
clear DirPath
configparams();





%% User Inputs
%Ask user what to do

printBreak
fprintf('What task should be started? Input [X] to select\n[A]nalysis\n[S]imulation\n')
taskselect = input('>>> ','s');
printBreak





%% Point program
%

if any(taskselect == 'a') || any(taskselect == 'A')
    [Mcounts,foldcent] = doAnikan();
elseif any(taskselect == 's') || any(taskselect == 'S')
    doSabine();
else
    error('Failed to assume task request. Try again with diff input [X]')
end
