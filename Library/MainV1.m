for tc = 15
clc
tc
close all
clearvars -except dataout tc
%clear all



% this program first finds one centroid with the pulsed analysis method
% where MLE searches for x and y (and the program uses area ratios to solve
% for Epsilon).

% next, it searches for two emitters, assumes phi1 = 1, and finds phi2 from
% tail fit to monoexponential. MLE searches for the two positions (and the
% program uses area ratios to solve for Epsilon1 and Epsilon 2).

set(0,'DefaultaxesFontName', 'Times new Roman')
set(0,'DefaultlegendFontName', 'Times new Roman')
set(groot,'defaultLineLineWidth',.8)
set(0,'DefaultaxesLineWidth', 1.2)
set(0,'DefaultaxesFontSize', 10)
set(0,'DefaultTextFontSize', 10)
set(0,'defaultfigurecolor',[1 1 1])
set(groot,'defaultAxesXGrid','on')
set(groot,'defaultAxesYGrid','on')


% dirlist = dir('D:\210325.sptw\selected210325\');
% dirlist = dir('D:\200629.sptw\selected200629\');
dirlist = dir('D:\210106.sptw\selected210106\');

% dirlist = dir('D:\210402.sptw\selected210402\');


clearvars -except cntfile dirlist tc dataout

simulated = 0;%1; %
%%%%%%tc = 30;%ns  CUT TIME PLAY AROUND WITH SINGLE EMMITER SHOULD NOT CHANGE ALOT
% pfia_list = 1;% pfia_list = 1-exp(-(1/10)*tc);%[0.80 0.9 0.95];
fittingPhi2 =1; % if fitting its fitting DECAY. MULTIPLE EMMITTERS BECOME DIFFICULT TO ASSIGN

bint = 10;%seconds   BIN TIME

wH = 0.96;%bin times in ns for the histogram of the decay for background and data

%cf1 = 1; cf2 = 1; cf3 =1; cf4 = 1;%1/6/21
cf1 = 1; cf2 = .8604; cf3 =.8455; cf4 = .8941;%1/6/21
sigmai= 180; sigmaiy = 180;%1/6/21
sigmaj= 180; sigmajy = 180;%1/6/21

% cf1 = 1; cf2 = 0.866; cf3 =0.855;cf4 = 0.910;%6.29.20
% sigmai= 180; sigmaiy = 187;%6.29.20
% sigmaj= 180; sigmajy = 187;%6.29.20

% cf1 = 1; cf2 = 0.86763; cf3 =0.85739;cf4 = 0.91090;%4/2/21
% sigmai= 174; sigmaiy = 172;%4/2/21
% sigmaj= 174; sigmajy = 172;%4/2/21

savePath = 'C:\Users\liamk\Desktop\QD_ver_4\RE__Alan_Van_Ordens_Zoom_Meeting\output220220';
if simulated ==0

    %folder = dirlist(1).folder;
    %fname = dirlist(cntfile).name;

    %Epath = dirlist(1).folder;
    Etag = 'QD_SingleDot';


    fname_E = which("210820_57.ptu");

    fname_dc = which("210820_94.ptu");
    fname_B = which("210820_88.ptu");


    %     fname_dc = strcat('D:\210325.sptw\210325tr_3.ptu');
    %     fname_B = strcat('D:\210325.sptw\210325tr_76.ptu');

    %     fname_B = strcat('D:\200629.sptw\200629tr_35.ptu');
    %     fname_dc = strcat('D:\200629.sptw\200629tr_34.ptu');

    %     fname_B = strcat('D:\210402.sptw\210128tr_46.ptu'); %AB
    %     fname_dc = strcat('D:\210402.sptw\210128tr_47.ptu');
    %% for when this was not in a big loop:

    %     Epath = 'D:\210106.sptw\';
    %     Etag = '210106tr_';
    %
    %     dc_num = 57; %number of dark count file
    %
    %
    %     B_num = 56;% no filter
    %     E_num = 51;% no filter
    %
    %     %     B_num = 54;% ABCD
    %     %     E_num = 53;% ABCD
    %     fname_dc = strcat(Epath,Etag,num2str(dc_num),'.ptu');
    %     fname_B = strcat(Epath,Etag,num2str(B_num),'.ptu');
    %     fname_E = strcat(Epath,Etag,num2str(E_num),'.ptu');
    %%
end

axisvar = 'Off';
posvar = 'Close';
LTcase = 'E2short';
Emvar = 'EprobEqual';
SimTimevar = '30';



if simulated ==1
    sigmai = 160;%in nm
    sigmaiy = 160;
    sigmaj = 160;%in nm
    sigmajy = 160;

    T=400;

    cf1 = 1; cf2 = 1; cf3 =1;cf4 = 1;%
    Sim_tag = strcat('simulate2E_Axis_',axisvar,'_Pos_',posvar,'_Case_',LTcase,'_Emission_',Emvar,'_SimTimeTest_',SimTimevar);
    fname_dc = strcat('C:\Users\liamk\Desktop\QDver1.2.3\RE__Alan_Van_Ordens_Zoom_Meeting\output220220\Sims\',Sim_tag,'\DC.mat');
    fname_B = strcat('C:\Users\liamk\Desktop\QDver1.2.3\RE__Alan_Van_Ordens_Zoom_Meeting\output220220\Sims\',Sim_tag,'\B.mat');
    fname_E = strcat('C:\Users\liamk\Desktop\QDver1.2.3\RE__Alan_Van_Ordens_Zoom_Meeting\output220220\Sims\',Sim_tag,'\E1E2.mat');

    %     fname_E = strcat('C:\Users\dunla\OneDrive\Desktop\programs210412\output210612\simulate210601_v3_E1E2.mat');

    %     fname_dc = 'C:\Users\dunla\OneDrive\Desktop\programs210412\output210612\simulate210601_v1_DC_plusIRF2.mat';
    %     fname_B = 'C:\Users\dunla\OneDrive\Desktop\programs210412\output210612\simulate210601_v1_Bkg_plusIRF2.mat';
    %     fname_E = strcat('C:\Users\dunla\OneDrive\Desktop\programs210412\output210612\simulate210601_v3_E1E2_plusIRF2.mat');
end

%program begins

% f = get_B_and_DC_Params_6ch_210408(fname_dc,fname_B,tc,'C:\Users\dunla\OneDrive\Desktop\programs210412\output\BandDC_params.mat',wH,simulated);
f = get_B_and_DC_Params_6ch_210612(fname_dc,fname_B,tc,'C:\Users\liamk\Desktop\QD_ver_4\RE__Alan_Van_Ordens_Zoom_Meeting\output220220\BandDC_params.mat',wH,simulated);

load('C:\Users\liamk\Desktop\QD_ver_4\RE__Alan_Van_Ordens_Zoom_Meeting\output220220\BandDC_params.mat');
N = bint/(T*10^-9);%number of periods in time bint_pos

[mBcounts] = get_mBcounts(fname_B,tc,bint,simulated); %mean of the total background counts in each channel per unit time  of bint

if simulated==0
    [output,txtout] =Read_PTU_V1_Barelli_fast(fname_E);% this function must be in the same directory as this progra
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
end
if simulated==1
    load(fname_E,'ch','dtime','sync')
    total = length(sync);
    ch = ch+1;

    T=400;

end

[sync, ch, dtime, total] = cutRepeats(sync,ch,dtime,total);
[sync, ch, dtime, total] = cutRepeats(sync,ch,dtime,total);

tmax = max(sync);

% if fittingPhi2 ==0
%     pfja_list = (1-exp(-(1/70)*tc))*ones(1,round(tmax/bint));
% end
clearvars -except sync ch dtime total N T tc bint ...
    fD1a fD2a fD3a fD4a fD1b fD2b fD3b fD4b...
    fBa fBb fEa fEb...
    PB1 PB2 PB3 PB4...
    ianalyze...
    fnum2 fname_B2 E_tag...
    hist_B1 hist_B2 hist_B3 hist_B4...
    ifi pfia_list positions pfja_list R2 simulated...
    cf1 cf2 cf3 cf4 sigmai sigmaiy sigmaj sigmajy...
    Sim_tag Etag savePath Epath fname_B fname_dc wH E_num cntfile dirlist
%% Part 1: find a single centroid location

%1)find phiE
Mcounts = binPhotons_tc(tc,bint,sync,dtime,ch);
% figure; plot(Mcounts(:,1))
[mBcounts] = get_mBcounts(fname_B,tc,bint,simulated); %mean of the total background counts in each channel per unit time  of bint
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
    cf1 cf2 cf3 cf4 sigmai sigmaiy sigmaj sigmajy savePath Sim_tag wH Etag E_num cntfile dirlist



thresholdM = 1000;
sM = sum(Mcounts,2);
x0 = NaN.*ones(1,size(Mcounts,1));
y0 = NaN.*ones(1,size(Mcounts,1));
P0 = NaN.*ones(1,size(Mcounts,1));
LnL0 = NaN.*ones(1,size(Mcounts,1));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for icnt = 6:9%1:size(Mcounts,1)
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
        load('C:\Users\liamk\Desktop\QD_ver_4\RE__Alan_Van_Ordens_Zoom_Meeting\output220220\sortingP_idxGuide_1emitter_backgrd_dc.mat');
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

load train
%sound(y,Fs)
figure('units','inch','position',[1,1,3,3]);
scatter(x0,y0,'k','filled'); hold on
grid on; box on;
set(gcf,'color','w');
xlabel('x (nm)');ylabel('y (nm)')
axis equal

%now look for two emitters

%1) Assume phi1 = 1; find phi2 by fitting the tail of the fluor decay.

f = get_B_and_DC_Params_6ch_210612(fname_dc,fname_B,tc,'C:\Users\liamk\Desktop\QD_ver_4\RE__Alan_Van_Ordens_Zoom_Meeting\output220220\BandDC_params.mat',wH,simulated);

% f = get_B_and_DC_Params_6ch_210408(fname_dc,fname_B,tc,'C:\Users\dunla\OneDrive\Desktop\programs210412\output\BandDC_params.mat',wH,simulated);
load('C:\Users\liamk\Desktop\QD_ver_4\RE__Alan_Van_Ordens_Zoom_Meeting\output220220\BandDC_params.mat');

pfia = 1; %Assume phi1 = 1
fittingPhi2 =1;   %WAS 1
if fittingPhi2 ==1
    disp(strcat('Getting phi2, the probability of arrivals from emitter 2 before the time tc = ',num2str(tc),' ns.'));
    disp('...')
    ihist = 1;
    leng = length(sync);
    n_bin = 1;
    i = 0;
    for iphoton = 1:leng
        i = i+1;
        if sync(iphoton)<bint*n_bin
            dt_n(i) = dtime(iphoton);
            ch_n(i) = ch(iphoton);
            Mt_n(i) = sync(iphoton);
        end
        if sync(iphoton)>=bint*n_bin
            [t_hist,y_hist] = get_histogram(dt_n,ch_n,T,0.96,2);
            yh1 = (y_hist./bint)-yB1;
            [t_hist,y_hist] = get_histogram(dt_n,ch_n,T,0.96,3);
            yh2 = (y_hist./bint)-yB2;
            [t_hist,y_hist] = get_histogram(dt_n,ch_n,T,0.96,4);
            yh3 = (y_hist./bint)-yB3;
            [t_hist,y_hist] = get_histogram(dt_n,ch_n,T,0.96,5);
            yh4 = (y_hist./bint)-yB4;

            yh = yh1+yh2+yh3+yh4;
            x = t_hist(tc:length(t_hist));%was 200
            y = log(yh(tc:length(t_hist)));
            [fitresult, gof] = LinearFit1(x, y);
            R2(n_bin) = -1*fitresult.m;
            bfit(n_bin)=fitresult.b;
            pfja_list(n_bin) = 1-exp(-R2(n_bin)*tc);


%             %             if floor(n_bin/10)==n_bin/10
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%FIGS%OF%DECAY%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
wantfigs = 1;
if wantfigs == 1
    figure('units','inch','position',[1,1,6,3]);
    subplot(1,2,1)
    plot(t_hist,log(yh),'k','linewidth',1);
    hold on
    plot(x,y,'k','linewidth',1);
    hold on
    plot(x,fitresult.m.*x+fitresult.b,'--r','linewidth',2)
    plot(x,fitresult.m.*x+fitresult.b,'--r','linewidth',2)
    grid on; box on
    xlabel(strcat('{\phi}_{\epsilon}_2 = ',num2str(pfja_list(n_bin)),', {\tau}_2=',num2str(1/R2(n_bin)),' ns'))
    xlim([0 400])
    ylim([1 inf])
    subplot(1,2,2)
    plot(t_hist,log(yh),'k','linewidth',1);
    hold on
    plot(x,y,'k','linewidth',1);
    hold on
    plot(x,fitresult.m.*x+fitresult.b,'--r','linewidth',2)
    grid on; box on
    xlabel(strcat('{\phi}_{\epsilon}_2 = ',num2str(pfja_list(n_bin)),', {\tau}_2=',num2str(1/R2(n_bin)),' ns'))
    xlim([0 3*tc])
    ylim([0 inf])
    aaatemp = 5;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             %             end

            clear x y fitresult gof t_hist y_hist
            n_bin=n_bin+1;%HELLO
            i=0;
            dt_n = 0;
            Mt_n=0;
            Ch_n=0;
            ihist=ihist+1;

        end
    end
    clear ihist leng n_bin i iphoton dt_n ch_n Mt_n t_hist y_hist x y
end
figure('units','inch','position',[1,1,3,3]);
scatter(x0,y0,'k','filled'); hold on
grid on; box on;
set(gcf,'color','w');
xlabel('x (nm)');ylabel('y (nm)')
axis equal

knowPhi = 0;   %WAS 0
if knowPhi==1
    disp('change known tau here!!!line664ish')
    Riknown = 1/5;  %WAS 1/25
    Rjknown = 1/122;  %WAS 1/50
    pfia = 1-exp(-Riknown*tc);%prob decay before tc for short 
    pfja = 1-exp(-Rjknown*tc);%prob decay before tc for long

    %     pfia_listFit = pfia_list;
    %     pfia_list = pfia+0*pfia_list;

    pfja_listFit = pfja_list;
    pfja_list = pfja+0*pfja_list;

    %     figure; plot(pfja_listFit,'k'); hold on; plot(pfja_list,'--k');
    %     plot(pfia_listFit,'r'); hold on; plot(pfia_list,'--r'); hold on

end
drawnow()

%2) Locate both emitters and P1 and P2

Results = NaN*ones(size(Mcounts,1),12);

for icnt = 6:9%1:size(pfja_list,2)
    fprintf('\nNow looking at bin: %d',icnt)
    badpoint = 0;
    testPhi2 = 1;
    testInit = 1;
    if pfja_list(icnt)<0.2
        testPhi2 = 0;
        fprintf('\nCan not find short lifetime. Moving to next bin')
        %testPhi2 = 1;
    end
    if isnan(x0(icnt))==1
        testInit = 0;
        fprintf('\nNew Error.Save Workspace and inputs for investigation. Moving to next bin\nicnt = %d',icnt)
        fprintf('\n%d',x0)
    end
    
    if testPhi2 ==1
        if testInit ==1

            load('C:\Users\liamk\Desktop\QD_ver_4\RE__Alan_Van_Ordens_Zoom_Meeting\output220220\sortingP_idxGuide_2emitters_backgrd_dc.mat');

            xyL0 = [x0(icnt) y0(icnt) x0(icnt) y0(icnt) pfia pfja_list(icnt)];

            %             if icnt ==1
            %                 init=1;
            %             end
            init =1;


            if init==1
                for iglobal = 1:100

                    a = x0(icnt)+10;
                    b = x0(icnt)-10;
                    rx = (b-a).*rand(1) + a;

                    a = y0(icnt)+10;
                    b = y0(icnt)-10;
                    ry = (b-a).*rand(1) + a;

                    a = x0(icnt)+10;
                    b = x0(icnt)-10;
                    rx2 = (b-a).*rand(1) + a;

                    a = y0(icnt)+10;
                    b = y0(icnt)-10;
                    ry2 = (b-a).*rand(1) + a;

                    gAll(iglobal,:) = [rx ry rx2 ry2 pfia pfja_list(icnt)];


                    [lnL(iglobal),P1(iglobal),P2(iglobal)] = getL_4p(Mcounts(icnt,:),N,rx,ry,rx2,ry2,...
                        sigmai,sigmaiy,pfia,...
                        sigmaj,sigmajy,pfja_list(icnt),...
                        fBa,fBb,...
                        PB1,PB2,PB3,PB4,...
                        fD1a,fD1b,fD2a,fD2b,...
                        fD3a,fD3b,fD4a,fD4b,...
                        key1a,key2a,key3a,key4a,...
                        key1b,key2b,key3b,key4b,...
                        keyMulti,keyNone,A1,A2,...
                        cf1,cf2,cf3,cf4);

                end
                if isnan(lnL)
                    lnL = 0;
                    disp('FAILED LINE 723; MOVING ON...')
                    badpoint = 1;
                end
                lnLmax = min(abs(lnL));
                row = find(lnL==-1*lnLmax);
                if length(row)>1
                    row = 1;
                    disp('Doh!!! multiple maxima found in global search, see line 494.')
                end
                xyL = gAll(row,:);
                xyLglobal = xyL
                clear L lnL gAll
            end
            %tempi = 1;

            irun = 0;
            deltaL =1;
            fprintf('\nFirst pass try')
            while deltaL>0
                irun = irun+1;
                
                fprintf('\n%d/31',irun)
                i =0;
                gAll = zeros(1,6);
                dl = 5;
                for n1 = -1*dl+xyL(1):dl:1*dl+xyL(1)
                    for n2 = -1*dl+xyL(2):dl:1*dl+xyL(2)%pos1:pos2:pos3 points of calc
                        for n3 = -1*dl+xyL(3):dl:1*dl+xyL(3)
                            for n4 =-1*dl+xyL(4):dl:1*dl+xyL(4)
                                for n5 = xyL(5)
                                    for n6 =xyL(6)
                                        i = i+1;
                                        gAll(i,:) = [n1 n2 n3 n4 n5 n6];
                                        [lnL(i),P1(i),P2(i)] = getL_4p(Mcounts(icnt,:),N,n1,n2,n3,n4,...
                                            sigmai,sigmaiy,n5,...
                                            sigmaj,sigmajy,n6,...
                                            fBa,fBb,...
                                            PB1,PB2,PB3,PB4,...
                                            fD1a,fD1b,fD2a,fD2b,...
                                            fD3a,fD3b,fD4a,fD4b,...
                                            key1a,key2a,key3a,key4a,...
                                            key1b,key2b,key3b,key4b,...
                                            keyMulti,keyNone,A1,A2,...
                                            cf1,cf2,cf3,cf4);                                       
                                    end
                                end
                            end
                        end
                    end
                end
                idx = find(lnL==max(lnL));
                if length(idx)>1
                    testFlat=1;
                    disp('it is flat')
                    deltaL=-1;
                elseif isnan(lnL)
                    lnL = 0;
                    disp('FAILED LINE 767, MOVING ON...')
                    deltaL=-1;
                    badpoint = 1;
                end
                if length(idx)==1
                    xyL = gAll(idx,:);
                    %lnLL(irun) = lnL(idx);%These 3 were commented out
                    %E1(irun) = P1(idx);
                    %E2(irun) = P2(idx);
                    clear dl n1 n2 n3 n4 n5 n6 n7 n8 i gAll L lnL idx P1 P2
                    if irun>1
                        deltaL = abs(lnLL(irun)-lnLL(irun-1));
                    end
                end
                if irun>30
                    disp('could not converge')
                    deltaL = -1;
                end

            end
            if deltaL>-1
                deltaL =1;
            end



            %             while deltaL>0
            %                 irun = irun+1;
            %                 i =0;
            %                 gAll = zeros(1,6);
            %                 dl = 1;
            %
            %                 for n1 = -1*dl+xyL(1):dl:1*dl+xyL(1)
            %                     for n2 = -1*dl+xyL(2):dl:1*dl+xyL(2)
            %                         for n3 = -1*dl+xyL(3):dl:1*dl+xyL(3)
            %                             for n4 =-1*dl+xyL(4):dl:1*dl+xyL(4)
            %                                 for n5 = xyL(5)
            %                                     for n6 =xyL(6)
            %                                         i = i+1;
            %                                         gAll(i,:) = [n1 n2 n3 n4 n5 n6];
            %                                         [lnL(i),P1(i),P2(i)] = getL_4p(Mcounts(icnt,:),N,n1,n2,n3,n4,...
            %                                             sigmai,sigmaiy,n5,...
            %                                             sigmaj,sigmajy,n6,...
            %                                             fBa,fBb,...
            %                                             PB1,PB2,PB3,PB4,...
            %                                             fD1a,fD1b,fD2a,fD2b,...
            %                                             fD3a,fD3b,fD4a,fD4b,...
            %                                             key1a,key2a,key3a,key4a,...
            %                                             key1b,key2b,key3b,key4b,...
            %                                             keyMulti,keyNone,A1,A2,...
            %                                             cf1,cf2,cf3,cf4);
            %                                     end
            %                                 end
            %                             end
            %                         end
            %                     end
            %                 end
            %                 idx = find(lnL==max(lnL));
            %                 if length(idx)>1
            %                     testFlat=1;
            %                     disp('it is flat')
            %                     deltaL=-1;
            %                 end
            %                 if length(idx)==1
            %                     xyL = gAll(idx,:);
            %                     %                 lnLL(irun) = lnL(idx);
            %                     %                 E1(irun) = P1(idx);
            %                     %                 E2(irun) = P2(idx);
            %                     clear dl n1 n2 n3 n4 n5 n6 n7 n8 i gAll L lnL idx P1 P2
            %                     if irun>1
            %                         deltaL = abs(lnLL(irun)-lnLL(irun-1));
            %                     end
            %                 end
            %                 if irun>30
            %                     disp('could not converge')
            %                     deltaL = -1;
            %                 end
            %
            %             end
            %             if deltaL>-1
            %                 deltaL =1;
            %             end

            %tempi = 1;
            fprintf('\nSecond pass try')
            while deltaL>0
                irun = irun+1;
                fprintf('\n%d/31',irun)
                i =0;
                gAll = zeros(1,6);
                dl = 1;

                for n1 = -1*dl+xyL(1):dl:1*dl+xyL(1)
                    for n2 = -1*dl+xyL(2):dl:1*dl+xyL(2)
                        for n3 = -1*dl+xyL(3):dl:1*dl+xyL(3)
                            for n4 =-1*dl+xyL(4):dl:1*dl+xyL(4)
                                for n5 = xyL(5)
                                    for n6 =xyL(6)
                                        i = i+1;
                                        gAll(i,:) = [n1 n2 n3 n4 n5 n6];
                                        [lnL(i),P1(i),P2(i)] = getL_4p(Mcounts(icnt,:),N,n1,n2,n3,n4,...
                                            sigmai,sigmaiy,n5,...
                                            sigmaj,sigmajy,n6,...
                                            fBa,fBb,...
                                            PB1,PB2,PB3,PB4,...
                                            fD1a,fD1b,fD2a,fD2b,...
                                            fD3a,fD3b,fD4a,fD4b,...
                                            key1a,key2a,key3a,key4a,...
                                            key1b,key2b,key3b,key4b,...
                                            keyMulti,keyNone,A1,A2,...
                                            cf1,cf2,cf3,cf4);
                                    end
                                end
                            end
                        end
                    end
                end
                idx = find(lnL==max(lnL));


                if length(idx)>1
                    testFlat=1;
                    disp('it is flat')
                    deltaL=-1;
                    %badpoint = 1;%new %toss this bin move on
                elseif isnan(lnL)
                    lnL = 0;
                    disp('FAILED LINE 888. MOVING ON...')
                    deltaL=-1
                    badpoint = 1;%toss and move on
                end
                if length(idx)==1
                    xyL = gAll(idx,:);
                    lnLL(irun) = lnL(idx);
                    E1(irun) = P1(idx);
                    E2(irun) = P2(idx);
                    clear dl n1 n2 n3 n4 n5 n6 n7 n8 i gAll L lnL idx P1 P2
                    if irun>1
                        deltaL = abs(lnLL(irun)-lnLL(irun-1));
                    end
                end
                if irun>30
                    disp('Could not converge')
                    deltaL = -1;
                    %badpoint = 1;%new%toss and move on
                end


            end
            if badpoint == 0
                fprintf('\nFinished bin %d',icnt)
            elseif badpoint == 1
                fprintf('\nFinished bin %d, bin is bad, moving on',icnt)
            end
            if badpoint == 0
                [nph1, nph2] = getNph2(Mcounts(icnt,:),cf1,cf2,cf3,cf4,mBcounts,sigmai);
                dummyvar(icnt,:) = [nph1, nph2];
                [error1a, error1b, error2a, error2b] = getPrecisionStandAlone2(Mcounts(icnt,:),mBcounts,xyL(1),xyL(2),xyL(3),xyL(4),nph1,nph2,cf1,cf2,cf3,cf4,sigmai,sigmaiy);
                dummyvar2(icnt,:) = [error1a, error1b, error2a, error2b];
                XYL(icnt,:) =xyL;
                EE1(icnt,1) = E1(irun);
                EE2(icnt,1) = E2(irun);
                logL(icnt,1) = lnLL(irun);
                xi = XYL(icnt,1);yi = XYL(icnt,2);
                xj = XYL(icnt,3);yj = XYL(icnt,4);
                R(icnt) = sqrt(((xi-xj)^2)+((yi-yj)^2));

                Results(icnt,:) = horzcat(tc,bint,XYL(icnt,:),EE1(icnt,1),EE2(icnt,1),logL(icnt,1),R(icnt));
                %disp('Found:')
                X1 = XYL(icnt,1);
                Y1 = XYL(icnt,2);
                e1 = EE1(icnt);
                X2 = XYL(icnt,3);
                Y2 = XYL(icnt,4);
                e2 = EE2(icnt);
                fprintf('\nFound E1 at \nX = %d\nY = %d\ne1 = %d\nAnd E2 at \nX = %d\nY = %d\ne2 = %d',XYL(icnt,1),XYL(icnt,2),EE1(icnt),XYL(icnt,3),XYL(icnt,4),EE2(icnt))

                if icnt>1
                    delta = abs([XYL(icnt-1,:)-XYL(icnt,:)])
                    if max(delta)>40
                        testDeltaDist(icnt)=0;
                    end
                end

                scatter(X1,Y1,'b','filled'); hold on
                scatter(X2,Y2,'r','filled'); hold on
                drawnow()
                fprintf('%d\n',icnt)
                clear xyL
            elseif badpoint == 1
                clear xyL
            end



        end










    end


end





for i = 1:size(XYL,1)
    Results(i,1)=tc;
    Results(i,2)=bint;
    Results(i,12)=R(i);
    Results(i,8)= pfja_list(i);
    Results(i,7) =pfia;
end

time = 0:bint:(icnt-1)*bint;
tempstop = 1;%fix E0 too!!!!!

X1 = Results(:,3);
Y1 = Results(:,4);
X2 = Results(:,5);
Y2 = Results(:,6);
E1 = Results(:,9);
E2 = Results(:,10);
E0 = Esolvedlist;
X0 = x0;
Y0 = y0;

X0 = X0(1:length(time));
X1 = X1(1:length(time));
X2 = X2(1:length(time));
Y0 = Y0(1:length(time));
Y1 = Y1(1:length(time));
Y2 = Y2(1:length(time));
E0 = E0(1:length(time));
E1 = E1(1:length(time));
E2 = E2(1:length(time));

figure('units','inch','position',[1,1,6.5,3]);
subplot(3,1,1)
plot(time,X0,'k'); hold on
plot(time,X1,'b'); hold on
plot(time,X2,'r'); hold on
ylabel('x (nm)')
legend('x0','x1','x2')

subplot(3,1,2)
plot(time,Y0,'k'); hold on
plot(time,Y1,'b'); hold on
plot(time,Y2,'r'); hold on
ylabel('y (nm)')
legend('y0','y1','y2')

subplot(3,1,3)
plot(time,E0,'k'); hold on
plot(time,E1,'b'); hold on
plot(time,E2,'r'); hold on
ylabel('{\epsilon}')
legend('{\epsilon}','{\epsilon}_1','{\epsilon}_2')

% scatter(x0,y0,'k','filled'); hold on
% grid on; box on;
% set(gcf,'color','w');
% xlabel('x (nm)');ylabel('y (nm)')
% axis equal


%% linear fit and subtraction

fitdrift =1;
if fitdrift==1
    [fitresultx, gof] = DriftFit2(time,X0);
    for i = 1:length(X0)
        X0f(i)=X0(i)-fitresultx.p1*time(i);
        X1f(i)=X1(i)-fitresultx.p1*time(i);
        X2f(i)=X2(i)-fitresultx.p1*time(i);

    end
    [fitresulty, gof] = DriftFit2(time,Y0);
    for i = 1:length(X0)
        Y0f(i)=Y0(i)-fitresulty.p1*time(i);
        Y1f(i)=Y1(i)-fitresulty.p1*time(i);
        Y2f(i)=Y2(i)-fitresulty.p1*time(i);
    end
end

figure('units','inch','position',[1,1,3,3]);
scatter(X0f,Y0f,'k','filled'); hold on
scatter(X1f,Y1f,'b','filled'); hold on
scatter(X2f,Y2f,'r','filled'); hold on

grid on; box on;
set(gcf,'color','w');
xlabel('x (nm)');ylabel('y (nm)')
axis equal
legend('Centroid','Short \tau','Long \tau')

n_0 = nnz(~isnan(X0f));
n_1 = nnz(~isnan(X1f));
n_2 = nnz(~isnan(X2f));


mnPos0 = [n_0 nanmean(X0f) nanstd(X0f) nanstd(X0f)/sqrt(n_0) nanmean(Y0f) nanstd(Y0f) nanstd(Y0f)/sqrt(n_0)]
mnPos1 = [n_1 nanmean(X1f) nanstd(X1f) nanstd(X1f)/sqrt(n_1) nanmean(Y1f) nanstd(Y1f) nanstd(Y1f)/sqrt(n_1)]
mnPos2 = [n_2 nanmean(X2f) nanstd(X2f) nanstd(X2f)/sqrt(n_2) nanmean(Y2f) nanstd(Y2f) nanstd(Y2f)/sqrt(n_2)]

vdiff = mnPos1-mnPos2;
vdiff(3) = sqrt(nanstd(X1)^2+nanstd(X2)^2);
vdiff(4) = mnPos1(4)+mnPos2(4);
vdiff(6) = sqrt(nanstd(Y1)^2+nanstd(Y2)^2);
vdiff(7) = mnPos1(7)+mnPos2(7);


mnPos = vertcat(mnPos0,mnPos1,mnPos2,vdiff);

savePath = 'C:\Users\liamk\Desktop\QD_ver_4\RE__Alan_Van_Ordens_Zoom_Meeting\output220220\SimsResults\';



if simulated==0
    %     XYname = strcat(savePath,'pos2centroids_',Etag,num2str(E_num),'_bint_',num2str(bint),'s_tc_',num2str(tc),'ns.mat');
    XYname = strcat(savePath,'pos2centroids_',Etag,'_bint_',num2str(bint),'s_tc_',num2str(tc),'ns.mat');

end
if simulated==1
    XYname = strcat(savePath,Sim_tag,'_bint_',num2str(bint),'s_tc_',num2str(tc),'ns.mat');
end

save(XYname);

resultdatatable = [nanmean(X0f) nanstd(X0f) nanmean(Y0f) nanstd(Y0f) nanmean(X1f) nanstd(X1f) nanmean(Y1f) nanstd(Y1f) nanmean(X2f) nanstd(X2f) nanmean(Y2f) nanstd(Y2f)]

load train
sound(y,Fs)
dataout(tc,1) = mnPos(4,2);
dataout(tc,2) = mnPos(4,4);
dataout(tc,3) = mnPos(4,5);
dataout(tc,4) = mnPos(4,7);
posstats(mnPos)
end