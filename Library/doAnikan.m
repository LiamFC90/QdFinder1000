function [Mcounts,foldcent] = doAnikan(jobtype)
%DOANIKAN Sets up and requests an analysis job
%   Detailed explanation goes here
printLine(20)
printBreak(2)
fprintf('Analysis: jobtype is %s\n',jobtype)
printBreak
%% Gather basic information


%%%foldbasic = [bint tc Wh sigmai sigmaiy sigmaj sigmajy];

%% start job
if strcmp(jobtype,'rawptu')
    [fold,foldpath] = get_AnikanInput(jobtype);
    get_dtimeplot(cell2mat(foldpath(1)),jobtype)
    %show user decay before asking for tc selection
    [foldbasic] = get_AnikanBasicParams;
    get_B_and_DC_Params_6ch(cell2mat(foldpath(3)),cell2mat(foldpath(2)),foldbasic(2),foldbasic(3),0,fold(1),fold(2),foldbasic(1))
    [Mcounts,foldcent] = get_centroidPOS(fold, foldpath, foldbasic, 0);

elseif strcmp(jobtype,'sim')
    [fold,foldpath] = get_AnikanInput(jobtype);
    get_dtimeplot(cell2mat(foldpath(1)),jobtype)
    [foldbasic] = get_AnikanBasicParams;
    simtag = strcat('X_',string(fold(1)),'_Y_',string(fold(2)),'_Tau_',string(fold(3)),'_rate_',string(fold(4)));
    fold = [fold simtag];
    get_B_and_DC_Params_6ch(cell2mat(foldpath(3)),cell2mat(foldpath(2)),foldbasic(2),foldbasic(3),1,55555,simtag,foldbasic(1))
    [Mcounts,foldcent] = get_centroidPOS(fold, foldpath, foldbasic, 1);
    % $work$
end

%%% fold = [folderE, fileE, folderB, fileB, folderD, fileD, CF1, CF2, CF3, CF4, tc, bint,sigmai, sigmaiy, sigmaj, sigmajy];
%%% foldpath = {pathE, pathB, pathD};
%%%%OR FOR SIM, foldpath point to sim.
%%%% fold = [XposNM, YposNM, TauVal, ratevl, simtag];


% get_dtimeplot(pathE)
% get_B_and_DC_Params_6ch(pathD,pathB,tc,Wh,simulated,folderE,fileE,bint)
% [Mcounts,foldcent] = get_centroidPOS(fold, foldpath, 0);


end

