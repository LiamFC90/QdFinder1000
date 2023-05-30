function [Mcounts,foldcent] = get_centroidPOS(fold,pathfold,foldbasic,simulated)
%GET_CENTROIDPOS finds a single centroid position for the system. 
%   Detailed explanation goes here


%%Decode fold
folderE = fold(1);
fileE = fold(2);
folderB = fold(3);
fileB = fold(4);
folderD = fold(5);
fileD = fold(6);
cf1 = fold(7); 
cf2 = fold(8); 
cf3 = fold(9); 
cf4 = fold(10);
bint = foldbasic(1);
tc = foldbasic(2);
Wh = foldbasic(3);
sigmai = foldbasic(4);
sigmaiy = foldbasic(5);
sigmaj = foldbasic(6);
sigmajy = foldbasic(7);
pathE = string(pathfold(1));
pathB = string(pathfold(2));
pathD = string(pathfold(3));

%% Load BandD file
load(which(strcat('BandD_tc',string(tc),'.mat')));
N = bint/(T*10^-9);%number of periods in time bint_pos

%% Part 1: find a single centroid location

[output,~] =Read_PTU_V1_Barelli_fast(pathE);% this function must be in the same directory as this progra
MeasDesc_GlobalResolution = output.Headers.MeasDesc_GlobalResolution;%macrotime resolution in seconds
MeasDesc_Resolution = output.MeasDesc_Resolution; %microt resolution in seconds
dtime = output.ph_dtime.*MeasDesc_Resolution*(10^9); %microt in nanoseconds
sync = output.ph_sync.*MeasDesc_GlobalResolution; % in seconds
total = size(sync,1);%total number of photons in dataset
ch = output.ph_channel;
T = round(output.Headers.MeasDesc_GlobalResolution*10^9,-2);
if max(ch)==7
    ch = ch-2;
end
[sync, ch, dtime, total] = cutRepeats(sync,ch,dtime,total);
[sync, ch, dtime, total] = cutRepeats(sync,ch,dtime,total);
tmax = max(sync);

%1)find phiE
Mcounts = binPhotons_tc(tc,bint,sync,dtime,ch);
% figure; plot(Mcounts(:,1))
[mBcounts] = get_mBcounts(pathB,tc,bint,simulated); %mean of the total background counts in each channel per unit time  of bint
j = 1;
start = 1;
for i = 1:total
    if sync(i)>j*bint
        ddtime = dtime(start:i-1);
        tau(j) = mean(ddtime);
        tau_std(j) = std(ddtime);
        tau_leng(j) = length(ddtime);
        cch = ch(start:i-1);
        ssync = sync(start:i-1);
        Hist_E1 = get_Hist_tc(ddtime',tc,cch,2);
        Hist_E2 = get_Hist_tc(ddtime',tc,cch,3);
        Hist_E3 = get_Hist_tc(ddtime',tc,cch,4);
        Hist_E4 = get_Hist_tc(ddtime',tc,cch,5);

        hist_E1 = (Hist_E1(1,:)./(max(ssync)-min(ssync)))-hist_B1;
        hist_E2 = (Hist_E2(1,:)./(max(ssync)-min(ssync)))-hist_B2;
        hist_E3 = (Hist_E3(1,:)./(max(ssync)-min(ssync)))-hist_B3;
        hist_E4 = (Hist_E4(1,:)./(max(ssync)-min(ssync)))-hist_B4;

        fE1a = hist_E1(1)/sum(hist_E1);
        fE1b = hist_E1(2)/sum(hist_E1);
        fE2a = hist_E2(1)/sum(hist_E2);
        fE2b = hist_E2(2)/sum(hist_E2);
        fE3a = hist_E3(1)/sum(hist_E3);
        fE3b = hist_E3(2)/sum(hist_E3);
        fE4a = hist_E4(1)/sum(hist_E4);
        fE4b = hist_E4(2)/sum(hist_E4);
        A = [fE1a fE2a fE3a fE4a fE1b fE2b fE3b fE4b];

        fEa(j) = mean([fE1a fE2a fE3a fE4a]);
        fEb(j) = mean([fE1b fE2b fE3b fE4b]);

        %                 end
        j= j+1;
        start = i;
        clear Hist_E1 Hist_E2 Hist_E3 Hist_E4...
            hist_E1 hist_E2 hist_E3 hist_E4...
            fE1a fE1b fE2a fE2b fE3a fE3b fE4a fE4b...
            ddtime cch ssync A ic
    end
end


clearvars -except sync ch dtime total bint  N T tc bint Mcounts mBcounts...
    fD1a fD2a fD3a fD4a fD1b fD2b fD3b fD4b...
    fBa fBb fEa fEb...
    PB1 PB2 PB3 PB4...
    ianalyze ...
    fnum2 fname_B fname_dc E_tag...
    pfia_list ifi pfja_list positions R2 simulated...
    cf1 cf2 cf3 cf4 sigmai sigmaiy sigmaj sigmajy savePath Sim_tag wH Etag E_num cntfile dirlist fold pathfold foldbasic




thresholdM = 1000; %MUST HAVE THIS MANY COUNTS TO BE A BIN
sM = sum(Mcounts,2);
x0 = NaN.*ones(1,size(Mcounts,1));
y0 = NaN.*ones(1,size(Mcounts,1));
P0 = NaN.*ones(1,size(Mcounts,1));
LnL0 = NaN.*ones(1,size(Mcounts,1));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for icnt = 1:size(Mcounts,1)
    disp(strcat(num2str(icnt),' of ',num2str(size(Mcounts,1))))
    if sM(icnt)>=thresholdM
        testThresh(icnt) =1;

    end
    if sM(icnt)<=thresholdM
        disp('under counts threshold, see line 176/237')
        testThresh(icnt) = 0;
    end
    if testThresh(icnt)==1
        A1 = (sum(Mcounts(icnt,1:4))-sum(mBcounts(1:4)))/N;
        A2 = (sum(Mcounts(icnt,5:8))-sum(mBcounts(5:8)))/N;
        %pixel constants
        global x1
        global x2
        global x3
        global x4
        global y1
        global y2
        global y3
        global y4
        global w
        mag = 250;%magnification
        w = 100*1000/mag;
        cladding = 15*1000/mag;
        % disp('cladding is zero')
        %                 cladding = 0;

        W = w+cladding;
        leng = W/2; %with cladding
        x1 = 1*leng; y1 = 1*leng;
        x2 = -1*leng; y2 = 1*leng;
        x3 = -1*leng; y3 = -1*leng;
        x4 = 1*leng; y4 = -1*leng;

        %         disp('Getting a single centroid position.')
        if ispc
            load("Library\sortingP_idxGuide_1emitter_backgrd_dc.mat")
        else
            load("Library/sortingP_idxGuide_1emitter_backgrd_dc.mat")
        end


        if icnt == 1
            a = -150;
            b = 150;
            rx = (b-a).*rand(1) + a;
            ry = (b-a).*rand(1) + a;
            emit = 1e-2;
            xyL = [rx ry emit];
        end
        if icnt>1
            if isnan(x0(icnt-1))==0
                a = x0(icnt-1)-20;
                b = x0(icnt-1)+20;
                rx = (b-a).*rand(1) + a;
                a = y0(icnt-1)-20;
                b = y0(icnt-1)+20;
                ry = (b-a).*rand(1) + a;
                emit = Esolvedlist(icnt-1);
                xyL = [rx ry emit];
            end
        end
        if icnt>1
            if isnan(x0(icnt-1))==1
                a = -150;
                b = 150;
                rx = (b-a).*rand(1) + a;
                ry = (b-a).*rand(1) + a;
                emit = 1e-2;
                xyL = [rx ry emit];
            end
        end

        irun = 0;
        deltaL =1;
        while deltaL>0
            irun = irun+1;
            i =0;
            gAll = zeros(1,2);
            dl = 50;
            dN = 1e-3;

            for n1 = -1*dl+xyL(1):dl:1*dl+xyL(1)
                for n2 = -1*dl+xyL(2):dl:1*dl+xyL(2)
                    i = i+1;
                    gAll(i,:) = [n1 n2];
                    [lnL(i),Pi(i)] = getL_41emitterSolvedP0(Mcounts(icnt,:),N,n1,n2,...
                        cf1,cf2,cf3,cf4,...
                        sigmai,sigmaiy,...
                        fEa(icnt),fEb(icnt),...
                        fBa,fBb,...
                        PB1,PB2,PB3,PB4,...
                        fD1a,fD1b,fD2a,fD2b,...
                        fD3a,fD3b,fD4a,fD4b,...
                        key1a,key2a,key3a,key4a,...
                        key1b,key2b,key3b,key4b,...
                        keyMulti,keyNone,...
                        A1,A2);

                end
            end

            idx = find(lnL==max(lnL));
            if length(idx)>1
                testFlat=1;
                disp('it is flat')
                xyL = [NaN NaN];
                deltaL=-1;
                megtest = lnL;
                Pisolved = NaN;

            end
            if length(idx)==1
                xyL = gAll(idx,:);
                lnLL(irun) = lnL(idx);
                Pisolved = Pi(idx);
                clear dl n1 n2 n3 n4 i gAll L lnL idx P0solved
                if irun>1
                    deltaL = abs(lnLL(irun)-lnLL(irun-1));
                end
            end
            if irun>200
                xyL = [NaN NaN];
                deltaL=-1;
                disp('could not converge')
                Pisolved = NaN;

            end
        end
        if deltaL>-1
            deltaL =1;
        end
        while deltaL>0
            irun = irun+1;
            i =0;
            gAll = zeros(1,2);
            dl = 10;
            dN = 1e-3;

            for n1 = -1*dl+xyL(1):dl:1*dl+xyL(1)
                for n2 = -1*dl+xyL(2):dl:1*dl+xyL(2)
                    i = i+1;
                    gAll(i,:) = [n1 n2];
                    [lnL(i),Pi(i)] = getL_41emitterSolvedP0(Mcounts(icnt,:),N,n1,n2,...
                        cf1,cf2,cf3,cf4,...
                        sigmai,sigmaiy,...
                        fEa(icnt),fEb(icnt),...
                        fBa,fBb,...
                        PB1,PB2,PB3,PB4,...
                        fD1a,fD1b,fD2a,fD2b,...
                        fD3a,fD3b,fD4a,fD4b,...
                        key1a,key2a,key3a,key4a,...
                        key1b,key2b,key3b,key4b,...
                        keyMulti,keyNone,...
                        A1,A2);


                end
            end

            idx = find(lnL==max(lnL));
            if length(idx)>1
                testFlat=1;
                disp('it is flat')
                xyL = [NaN NaN];
                deltaL=-1;
                megtest = lnL;
                Pisolved = NaN;

            end
            if length(idx)==1
                xyL = gAll(idx,:);
                lnLL(irun) = lnL(idx);
                Pisolved = Pi(idx);
                clear dl n1 n2 n3 n4 i gAll L lnL idx P0solved
                if irun>1
                    deltaL = abs(lnLL(irun)-lnLL(irun-1));
                end
            end
            if irun>200
                xyL = [NaN NaN];
                deltaL=-1;
                disp('could not converge')
                Pisolved = NaN;

            end
        end
        if deltaL>-1
            deltaL =1;
        end
        while deltaL>0
            irun = irun+1;
            i =0;
            gAll = zeros(1,2);
            dl = 1;
            dN = 1e-3;

            for n1 = -1*dl+xyL(1):dl:1*dl+xyL(1)
                for n2 = -1*dl+xyL(2):dl:1*dl+xyL(2)
                    i = i+1;
                    gAll(i,:) = [n1 n2];
                    [lnL(i),Pi(i)] = getL_41emitterSolvedP0(Mcounts(icnt,:),N,n1,n2,...
                        cf1,cf2,cf3,cf4,...
                        sigmai,sigmaiy,...
                        fEa(icnt),fEb(icnt),...
                        fBa,fBb,...
                        PB1,PB2,PB3,PB4,...
                        fD1a,fD1b,fD2a,fD2b,...
                        fD3a,fD3b,fD4a,fD4b,...
                        key1a,key2a,key3a,key4a,...
                        key1b,key2b,key3b,key4b,...
                        keyMulti,keyNone,...
                        A1,A2);


                end
            end

            idx = find(lnL==max(lnL));
            if length(idx)>1
                testFlat=1;
                disp('it is flat')
                xyL = [NaN NaN];
                deltaL=-1;
                megtest = lnL;
                Pisolved = NaN;

            end
            if length(idx)==1
                xyL = gAll(idx,:);
                lnLL(irun) = lnL(idx);
                Pisolved = Pi(idx);
                clear dl n1 n2 n3 n4 i gAll L lnL idx P0solved
                if irun>1
                    deltaL = abs(lnLL(irun)-lnLL(irun-1));
                end
            end
            if irun>200
                xyL = [NaN NaN];
                deltaL=-1;
                disp('could not converge')
                Pisolved = NaN;

            end
        end
        if deltaL>-1
            deltaL =1;
        end
        while deltaL>0
            irun = irun+1;
            i =0;
            gAll = zeros(1,2);
            dl = .1;
            dN = 1e-4;

            for n1 = -1*dl+xyL(1):dl:1*dl+xyL(1)
                for n2 = -1*dl+xyL(2):dl:1*dl+xyL(2)
                    i = i+1;
                    gAll(i,:) = [n1 n2];
                    [lnL(i),Pi(i)] = getL_41emitterSolvedP0(Mcounts(icnt,:),N,n1,n2,...
                        cf1,cf2,cf3,cf4,...
                        sigmai,sigmaiy,...
                        fEa(icnt),fEb(icnt),...
                        fBa,fBb,...
                        PB1,PB2,PB3,PB4,...
                        fD1a,fD1b,fD2a,fD2b,...
                        fD3a,fD3b,fD4a,fD4b,...
                        key1a,key2a,key3a,key4a,...
                        key1b,key2b,key3b,key4b,...
                        keyMulti,keyNone,...
                        A1,A2);


                end
            end

            idx = find(lnL==max(lnL));
            if length(idx)>1
                testFlat=1;
                disp('it is flat')
                xyL = [NaN NaN];
                deltaL=-1;
                megtest = lnL;
                Pisolved = NaN;

            end
            if length(idx)==1
                xyL = gAll(idx,:);
                lnLL(irun) = lnL(idx);
                Pisolved = Pi(idx);
                clear dl n1 n2 n3 n4 i gAll L lnL idx P0solved
                if irun>1
                    deltaL = abs(lnLL(irun)-lnLL(irun-1));
                end
            end
            if irun>200
                xyL = [NaN NaN];
                deltaL=-1;
                disp('could not converge')
                Pisolved = NaN;

            end
        end
        if deltaL>-1
            deltaL =1;
        end
        x0(icnt) = xyL(1,1);
        y0(icnt) = xyL(1,2);
        
        Esolvedlist(icnt) = Pisolved;
        if isnan(x0(icnt))==0
            LnL0(icnt) = lnLL(irun);
        end
        clear n1 n2 n3 i gAll L lnL idx P0solved
    end

end
time = 0:bint:(icnt-1)*bint;
x0 = x0(1:length(time));
y0 = y0(1:length(time));
[dataX,~] = DriftFit2(time,x0);
[dataY,~] = DriftFit2(time,y0);
for n = 1:numel(icnt)
    x01(n) = x0(n) - dataX.p1*time(n);
    y01(n) = y0(n) - dataY.p1*time(n);
end
if all(x01 == x0) && all(y01 == y0)
    fprintf('\n line 445 test\n')
end
%%Decode fold 2
folderE = fold(1);
fileE = fold(2);
folderB = fold(3);
fileB = fold(4);
folderD = fold(5);
fileD = fold(6);
cf1 = fold(7); 
cf2 = fold(8); 
cf3 = fold(9); 
cf4 = fold(10);
bint = foldbasic(1);
tc = foldbasic(2);
Wh = foldbasic(3);
sigmai = foldbasic(4);
sigmaiy = foldbasic(5);
sigmaj = foldbasic(6);
sigmajy = foldbasic(7);
pathE = string(pathfold(1));
pathB = string(pathfold(2));
pathD = string(pathfold(3));
foldcent = [x0' y0' Esolvedlist' LnL0']; %new foldvar for CENT POS


figure('units','inch','position',[1,1,3,3]);
scatter(x0,y0,'k','filled'); hold on
grid on; box on;
set(gcf,'color','w');
xlabel('x (nm)');ylabel('y (nm)')
axis equal
title(strcat(string(folderE),'_',string(fileE),'_tc',string(tc),'ns_bint',string(bint),'s'),Interpreter="none")


if ispc
    savepath = strcat("Analysis\",string(folderE),'\',string(fileE));
    temppath = cd(savepath);
    save(strcat('Centroid_tc',string(tc),'ns_bint',string(bint),'s.mat'),"foldcent")
    exportgraphics(gcf,strcat('Centroid_tc',string(tc),'ns_bint',string(bint),'s.png'))
    cd(temppath)
elseif isunix || ismac
    savepath = strcat("Analysis/",string(folderE),'/',string(fileE));
    temppath = cd(savepath);
    save(strcat('Centroid_tc',string(tc),'ns_bint',string(bint),'s.mat'),"foldcent")
    exportgraphics(gcf,strcat('Centroid_tc',string(tc),'ns_bint',string(bint),'s.png'))
    cd(temppath)
end
fprintf('Found Centroid at:\nX = %.2f +- %.2f\nY = %.2f +- %.2f\n',mean(x0),std(x0),mean(y0),std(y0))
end

