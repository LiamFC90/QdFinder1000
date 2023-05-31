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
TalkToUser = 'engage';
while exist('TalkToUser','var')
printBreak(2)

fprintf('What task should be started? Input [X] to select\n[A]nalysis\n[S]imulation\n')
taskselect = input('>>> ','s');
printBreak





%% Point program
%

decidelogic = ninputStringLogic('A','S',taskselect);
if decidelogic == 1 %do analysis
    clear TalkToUser
    printLine(20)
    printBreak
    simlogic = ninputStringLogic('D','S',input('Load a [D]ata file or [S]imulation file?\n>>>','s'));
        if simlogic == 1 %do data file
            [Mcounts,foldcent] = doAnikan('rawptu');
        elseif simlogic ==2 %do sim file
            [Mcounts,foldcent] = doAnikan('sim');
        else
            warning('Failed to assume file type request. Restarting')
            TalkToUser = 'engage';
            printBreak(2)      
        end
elseif decidelogic == 2%do simulation
    clear TalkToUser
    printLine(20)
    printBreak(2)
    fprintf('Preparing Simulation:\n')
    doSabine();
else
    warning('Failed to assume task request [%s]. Try again with diff input [X]',taskselect)
end

end

printBreak
fprintf('Program complete, rerun Main to restart.\n')
printBreak







