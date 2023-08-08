function [Mcounts,foldcent] = get_centroidPOS(fold,pathfold,foldbasic,simulated)
%GET_CENTROIDPOS finds a single centroid position for the system.
%   Detailed explanation goes here


%%Decode fold
if simulated == 0
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
end
if simulated == 1
    cf1 = 1;
    cf2 = 1;
    cf3 = 1;
    cf4 = 1;
end

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
if simulated == 0
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
end
if simulated ==1
    load(pathE,'dtime','ch','sync','total')
    ch = ch + 1;
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
    cf1 cf2 cf3 cf4 sigmai sigmaiy sigmaj sigmajy savePath Sim_tag wH Etag E_num cntfile dirlist fold pathfold foldbasic fileE folderE




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

if ispc
    
    load(strcat("Analysis\",string(folderE),"\",string(fileE),"\BandD_tc",string(tc),".mat"))
else
    load(strcat("Analysis/",string(FolderE),"/",string(fileE),"\BandD_tc",string(tc),".mat"))
end

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
            [fitresult, ~] = LinearFit1(x, y);
            R2(n_bin) = -1*fitresult.m;
            bfit(n_bin)=fitresult.b;
            pfja_list(n_bin) = 1-exp(-R2(n_bin)*tc);


%             %             if floor(n_bin/10)==n_bin/10
            % figure('units','inch','position',[1,1,6,3]);
            % subplot(1,2,1)
            % plot(t_hist,log(yh),'k','linewidth',1);
            % hold on
            % plot(x,y,'k','linewidth',1);
            % hold on
            % plot(x,fitresult.m.*x+fitresult.b,'--r','linewidth',2)
            % plot(x,fitresult.m.*x+fitresult.b,'--r','linewidth',2)
            % grid on; box on
            % xlabel(strcat('{\phi}_{\epsilon}_2 = ',num2str(pfja_list(n_bin)),', {\tau}_2=',num2str(1/R2(n_bin)),' ns'))
            % xlim([0 400])
            % ylim([1 inf])
            % subplot(1,2,2)
            % plot(t_hist,log(yh),'k','linewidth',1);
            % hold on
            % plot(x,y,'k','linewidth',1);
            % hold on
            % plot(x,fitresult.m.*x+fitresult.b,'--r','linewidth',2)
            % grid on; box on
            % xlabel(strcat('{\phi}_{\epsilon}_2 = ',num2str(pfja_list(n_bin)),', {\tau}_2=',num2str(1/R2(n_bin)),' ns'))
            % xlim([0 3*tc])
            % ylim([0 inf])
            % aaatemp = 5;
% 
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
    Rjknown = 1/50;  %WAS 1/50
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

for icnt = 1:size(pfja_list,2)
    badpoint = 0;
    testPhi2 = 1;
    testInit = 1;
    if pfja_list(icnt)<0.2
        testPhi2 = 0;
        %testPhi2 = 1;
    end
    if isnan(x0(icnt))==1
        testInit = 0;
    end
    
    if testPhi2 ==1
        if testInit ==1

            
            if ispc
                load("Library\sortingP_idxGuide_2emitters_backgrd_dc.mat")
            else
                load("Library/sortingP_idxGuide_2emitters_backgrd_dc.mat")
            end

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
            while deltaL>0
                irun = irun+1;
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
            while deltaL>0
                irun = irun+1;
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
                elseif isnan(lnL)
                    lnL = 0;
                    disp('FAILED LINE 888. MOVING ON...')
                    deltaL=-1
                    badpoint = 1;
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
                    disp('could not converge')
                    deltaL = -1;
                end


            end

            if badpoint == 0
                XYL(icnt,:) =xyL;
                EE1(icnt,1) = E1(irun);
                EE2(icnt,1) = E2(irun);
                logL(icnt,1) = lnLL(irun);
                xi = XYL(icnt,1);yi = XYL(icnt,2);
                xj = XYL(icnt,3);yj = XYL(icnt,4);
                R(icnt) = sqrt(((xi-xj)^2)+((yi-yj)^2));

                Results(icnt,:) = horzcat(tc,bint,XYL(icnt,:),EE1(icnt,1),EE2(icnt,1),logL(icnt,1),R(icnt));
                disp('Found:')
                X1 = XYL(icnt,1)
                Y1 = XYL(icnt,2)
                e1 = EE1(icnt)
                X2 = XYL(icnt,3)
                Y2 = XYL(icnt,4)
                e2 = EE2(icnt)

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


% if simulated == 0 %do drift correction if not sim data
%     time = 0:bint:(icnt-1)*bint;
%     x0 = x0(1:length(time));
%     y0 = y0(1:length(time));
%     [dataX,~] = DriftFit2(time,x0);
%     [dataY,~] = DriftFit2(time,y0);
%     for n = 1:icnt
%         x01(n) = x0(n) - dataX.p1*time(n);
%         y01(n) = y0(n) - dataY.p1*time(n);
%     end
%     if all(x01 == x0) && all(y01 == y0)
%         fprintf('\n line 445 DRIFT CHECK. IF THIS PRINTS, DRIFT CORRECT IS NOT WORKING! KNOWN ISSUE\n')
%     end
% end




%%Decode fold 2nd time around since it gets cleared.
if simulated == 0
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
end
if simulated == 1
    folderE = 55555;
    fileE = fold(5);
    cf1 = 1;
    cf2 = 1;
    cf3 = 1;
    cf4 = 1;
end

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

% x01 = x01(1:length(Esolvedlist));
% y01 = y01(1:length(Esolvedlist));
% LnL0 = LnL0(1:length(Esolvedlist));
% 
tempsize = numel(X0f);
Y0f = Y0f(1:tempsize);
Esolvedlist = Esolvedlist(1:tempsize);
LnL0 = LnL0(1:tempsize);
clear tempsize
foldcent = [X0f' Y0f' Esolvedlist' LnL0']; %new foldvar for CENT POS

figure('units','inch','position',[1,1,3,3]);
scatter(X0f,Y0f,'k','filled'); hold on
scatter(X1f,Y1f,'b','filled'); hold on
scatter(X2f,Y2f,'r','filled'); hold on
grid on; box on;
set(gcf,'color','w');
xlabel('x (nm)');ylabel('y (nm)')
axis equal
title(strcat(string(folderE),'_',string(fileE),'_tc',string(tc),'ns_bint',string(bint),'s'),Interpreter="none")
legend('Centroid','Short','Long');


if ispc
    savepath = strcat("Analysis\",string(folderE),'\',string(fileE));
    temppath = cd(savepath);
    save(strcat('Centroid_tc',string(tc),'ns_bint',string(bint),'s.mat'),"foldcent")
    exportgraphics(gcf,strcat('Centroid_tc',string(tc),'ns_bint',string(bint),'s.png'))
    %exportgraphics(hdecay,strcat('Hdecay_tc',string(tc),'ns_bint',string(bint),'s.png'))
    cd(temppath)
elseif isunix || ismac
    savepath = strcat("Analysis/",string(folderE),'/',string(fileE));
    temppath = cd(savepath);
    save(strcat('Centroid_tc',string(tc),'ns_bint',string(bint),'s.mat'),"foldcent")
    exportgraphics(gcf,strcat('Centroid_tc',string(tc),'ns_bint',string(bint),'s.png'))
    %exportgraphics(hdecay,strcat('Hdecay_tc',string(tc),'ns_bint',string(bint),'s.png'))
    cd(temppath)
end
printBreak
fprintf('Found Centroid at:\nX = %.2f +- %.2f\nY = %.2f +- %.2f\n',mean(X0f,'omitnan'),std(X0f,'omitnan'),mean(Y0f,'omitnan'),std(Y0f,'omitnan'))
printBreak
printLine
end

