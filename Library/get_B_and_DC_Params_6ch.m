 
function [] = get_B_and_DC_Params_6ch(fname_dc,fname_B,tc,wH,simulated,folderE,fileE,bint)

if simulated==0
[output,txtout] =Read_PTU_V1_Barelli_fast(fname_dc);% this function must be in the same directory as this progra
MeasDesc_GlobalResolution = output.Headers.MeasDesc_GlobalResolution;%macrotime resolution in seconds
MeasDesc_Resolution = output.MeasDesc_Resolution; %microt resolution in seconds
ddtime = output.ph_dtime.*MeasDesc_Resolution*(10^9); %microt in nanoseconds
ssync = output.ph_sync.*MeasDesc_GlobalResolution; % in seconds
total = size(ssync,1);%total number of photons in dataset
ch = output.ph_channel;
T = round(output.Headers.MeasDesc_GlobalResolution*10^9,-2);

end
if simulated==1
    load(fname_dc,'ch','dtime','sync')
    ssync = sync;
    ddtime = dtime;
    ch = ch+1;
    
    T=400
end

if max(ch)==7
    ch = ch-2;
end

bintDC = max(ssync)-0.1*max(ssync);

M_dc = binPhotons_tc(tc,bintDC,ssync,ddtime,ch);
M1_dc = mean(M_dc(:,1)+M_dc(:,5));
M2_dc = mean(M_dc(:,2)+M_dc(:,6));
M3_dc = mean(M_dc(:,3)+M_dc(:,7));
M4_dc = mean(M_dc(:,4)+M_dc(:,8));


RD1=mean(M1_dc)*10^-9/bintDC;
RD2=mean(M2_dc)*10^-9/bintDC;
RD3=mean(M3_dc)*10^-9/bintDC;
RD4=mean(M4_dc)*10^-9/bintDC;

fD1a = 1-exp(-RD1*tc);
fD1b = exp(-RD1*tc)-exp(-RD1*T);
fD2a = 1-exp(-RD2*tc);
fD2b = exp(-RD2*tc)-exp(-RD2*T);
fD3a = 1-exp(-RD3*tc);
fD3b = exp(-RD3*tc)-exp(-RD3*T);
fD4a = 1-exp(-RD4*tc);
fD4b = exp(-RD4*tc)-exp(-RD4*T);

Hist_dc1 = get_Hist_tc(ddtime',tc,ch,2);
Hist_dc2 = get_Hist_tc(ddtime',tc,ch,3);
Hist_dc3 = get_Hist_tc(ddtime',tc,ch,4);
Hist_dc4 = get_Hist_tc(ddtime',tc,ch,5);

hist_dc1 = Hist_dc1(1,:)./max(ssync)
hist_dc2 = Hist_dc2(1,:)./max(ssync);
hist_dc3 = Hist_dc3(1,:)./max(ssync);
hist_dc4 = Hist_dc4(1,:)./max(ssync);



clearvars -except T tc bint_E fname_dc fname_B fname_E...
    fD1a fD1b fD2a fD2b fD3a fD3b fD4a fD4b...
    RD1 RD2 RD3 RD4...
    hist_dc1 hist_dc2 hist_dc3 hist_dc4 tag fname_params wH simulated fileE folderE

if simulated==0
[output,txtout] =Read_PTU_V1_Barelli_fast(fname_B);% this function must be in the same directory as this progra
MeasDesc_GlobalResolution = output.Headers.MeasDesc_GlobalResolution;%macrotime resolution in seconds
MeasDesc_Resolution = output.MeasDesc_Resolution; %microt resolution in seconds
ddtime = output.ph_dtime.*MeasDesc_Resolution*(10^9); %microt in nanoseconds
ssync = output.ph_sync.*MeasDesc_GlobalResolution; % in seconds
total = size(ssync,1);%total number of photons in dataset
ch = output.ph_channel;
T = round(output.Headers.MeasDesc_GlobalResolution*10^9,-2);

end
if simulated==1
    load(fname_B,'ch','dtime','sync')
    ssync = sync;
    ddtime = dtime;
        ch = ch+1;
    T=400

end

if max(ch)==7
    ch = ch-2;
end


[xB,y1] = get_histogram(ddtime,ch,T,wH,2);
yB1 = y1./max(ssync);
[xB,y2] = get_histogram(ddtime,ch,T,wH,3);
yB2 = y2./max(ssync);
[xB,y3] = get_histogram(ddtime,ch,T,wH,4);
yB3 = y3./max(ssync);
[xB,y4] = get_histogram(ddtime,ch,T,wH,5);
yB4 = y4./max(ssync);

Hist_b1 = get_Hist_tc(ddtime',tc,ch,2)
Hist_b2 = get_Hist_tc(ddtime',tc,ch,3);
Hist_b3 = get_Hist_tc(ddtime',tc,ch,4);
Hist_b4 = get_Hist_tc(ddtime',tc,ch,5);
hist_B1 = (Hist_b1(1,:)./max(ssync))
hist_B2 = (Hist_b2(1,:)./max(ssync));
hist_B3 = (Hist_b3(1,:)./max(ssync));
hist_B4 = (Hist_b4(1,:)./max(ssync));
hist_b1 = (Hist_b1(1,:)./max(ssync))-hist_dc1
hist_b2 = (Hist_b2(1,:)./max(ssync))-hist_dc2;
hist_b3 = (Hist_b3(1,:)./max(ssync))-hist_dc3;
hist_b4 = (Hist_b4(1,:)./max(ssync))-hist_dc4;


test = min(min([hist_b1(1);hist_b2(1);hist_b3(1);hist_b4(1)]));
if test<0
    hist_B1=hist_dc1;
    hist_B2 = hist_dc2;
    hist_B3 = hist_dc3;
    hist_B4 = hist_dc4;
    PB1=1e-12;
    PB2=1e-12;
    PB3=1e-12;
    PB4=1e-12;
    fBa=1e-12;
    fBb=1e-12;
end
% test2 = min(min([hist_b1(2);hist_b2(2);hist_b3(2);hist_b4(2)]));
% if test2>0
%     hist_b1(2) = 0;
%     hist_b2(2) = 0;
%     hist_b3(2) = 0;
%     hist_b4(2) = 0;
% end




if test>0
    fB1a = hist_b1(1)/sum(hist_b1);
    fB1b = hist_b1(2)/sum(hist_b1);
    fB2a = hist_b2(1)/sum(hist_b2);
    fB2b = hist_b2(2)/sum(hist_b2);
    fB3a = hist_b3(1)/sum(hist_b3);
    fB3b = hist_b3(2)/sum(hist_b3);
    fB4a = hist_b4(1)/sum(hist_b4);
    fB4b = hist_b4(2)/sum(hist_b4);
    
    fBa = mean([fB1a fB2a fB3a fB4a]);
    fBb = mean([fB1b fB2b fB3b fB4b]);
    
    bintB = max(ssync)-0.1*max(ssync);
    Mcounts_background = binPhotons_tc(tc,bintB,ssync,ddtime,ch);
    if size(Mcounts_background,1)>1
    Mcounts_background = Mcounts_background(1:length(Mcounts_background)-1,:);
    end
    M_bk = mean(Mcounts_background,1)./bintB;
    m1a = M_bk(1,1);
    m1b = M_bk(1,5);
    m2a = M_bk(1,2);
    m2b = M_bk(1,6);
    m3a = M_bk(1,3);
    m3b = M_bk(1,7);
    m4a = M_bk(1,4);
    m4b = M_bk(1,8);
    
    N = 1/(T*10^-9);%number of periods in time t
    for simplex=1:1
        %% simplex analysis for pixel 1
        % get random starting conditions for probing paramter space
        %         a = .00001;
        %         b = .001;
        %         r1 = (b-a).*rand(1) + a;
        r1 = 1e-5;
        xyL = [r1];
        irun = 0;
        deltaL =1;
        while deltaL>0
            irun = irun+1;
            i =0;
            gAll = zeros(1,1);
            dl = 1e-4;
            for r1 = 0:dl:1*dl+xyL(1)
                i = i+1;
                gAll(i,:) = [r1];
                [L(i),lnL(i)] = getL_P_background(N,m1a,m1b,fBa,fBb,...
                    r1,...
                    fD1a,fD1b);
            end
            idx = find(lnL==max(lnL));
            xyL = gAll(idx,:);
            lnLL(irun) = lnL(idx);
            clear dl n1 n2 n3 n4 i gAll L lnL idx
            if irun>1
                deltaL = abs(lnLL(irun)-lnLL(irun-1));
            end
        end
        deltaL =1;
        while deltaL>0
            irun = irun+1;
            i =0;
            gAll = zeros(1,1);
            dl = .00001;
            for r1 = 0:dl:1*dl+xyL(1)
                i = i+1;
                gAll(i,:) = [r1];
                [L(i),lnL(i)] = getL_P_background(N,m1a,m1b,fBa,fBb,...
                    r1,...
                    fD1a,fD1b);
            end
            idx = find(lnL==max(lnL));
            xyL = gAll(idx,:);
            lnLL(irun) = lnL(idx);
            clear dl n1 n2 n3 n4 i gAll L lnL idx
            if irun>1
                deltaL = abs(lnLL(irun)-lnLL(irun-1));
            end
        end
        PB1 = xyL;
        clear xyL
        %% simplex analysis for pixel 2
        % get random starting conditions for probing paramter space
        a = .00001;
        b = .001;
        r1 = (b-a).*rand(1) + a;
        xyL = [r1];
        irun = 0;
        deltaL =1;
        while deltaL>0
            irun = irun+1;
            i =0;
            gAll = zeros(1,1);
            dl = .0001;
            for r1 = 0:dl:1*dl+xyL(1)
                i = i+1;
                gAll(i,:) = [r1];
                [L(i),lnL(i)] = getL_P_background(N,m2a,m2b,fBa,fBb,...
                    r1,...
                    fD2a,fD2b);
            end
            idx = find(lnL==max(lnL));
            xyL = gAll(idx,:);
            lnLL(irun) = lnL(idx);
            clear dl n1 n2 n3 n4 i gAll L lnL idx
            if irun>1
                deltaL = abs(lnLL(irun)-lnLL(irun-1));
            end
        end
        deltaL =1;
        while deltaL>0
            irun = irun+1;
            i =0;
            gAll = zeros(1,1);
            dl = .00001;
            for r1 = 0:dl:1*dl+xyL(1)
                i = i+1;
                gAll(i,:) = [r1];
                [L(i),lnL(i)] = getL_P_background(N,m2a,m2b,fBa,fBb,...
                    r1,...
                    fD2a,fD2b);
            end
            idx = find(lnL==max(lnL));
            xyL = gAll(idx,:);
            lnLL(irun) = lnL(idx);
            clear dl n1 n2 n3 n4 i gAll L lnL idx
            if irun>1
                deltaL = abs(lnLL(irun)-lnLL(irun-1));
            end
        end
        PB2 = xyL;
        clear xyL
        %% simplex analysis for pixel 3
        % get random starting conditions for probing paramter space
        a = .00001;
        b = .001;
        r1 = (b-a).*rand(1) + a;
        xyL = [r1];
        irun = 0;
        deltaL =1;
        while deltaL>0
            irun = irun+1;
            i =0;
            gAll = zeros(1,1);
            dl = .0001;
            for r1 = 0:dl:1*dl+xyL(1)
                i = i+1;
                gAll(i,:) = [r1];
                [L(i),lnL(i)] = getL_P_background(N,m3a,m3b,fBa,fBb,...
                    r1,...
                    fD3a,fD3b);
            end
            idx = find(lnL==max(lnL));
            xyL = gAll(idx,:);
            lnLL(irun) = lnL(idx);
            clear dl n1 n2 n3 n4 i gAll L lnL idx
            if irun>1
                deltaL = abs(lnLL(irun)-lnLL(irun-1));
            end
        end
        deltaL =1;
        while deltaL>0
            irun = irun+1;
            i =0;
            gAll = zeros(1,1);
            dl = .00001;
            for r1 = 0:dl:1*dl+xyL(1)
                i = i+1;
                gAll(i,:) = [r1];
                [L(i),lnL(i)] = getL_P_background(N,m3a,m3b,fBa,fBb,...
                    r1,...
                    fD3a,fD3b);
            end
            idx = find(lnL==max(lnL));
            xyL = gAll(idx,:);
            lnLL(irun) = lnL(idx);
            clear dl n1 n2 n3 n4 i gAll L lnL idx
            if irun>1
                deltaL = abs(lnLL(irun)-lnLL(irun-1));
            end
        end
        PB3 = xyL;
        clear xyL;
        %% simplex analysis for pixel 4
        % get random starting conditions for probing paramter space
        a = .00001;
        b = .001;
        r1 = (b-a).*rand(1) + a;
        xyL = [r1];
        irun = 0;
        deltaL =1;
        while deltaL>0
            irun = irun+1;
            i =0;
            gAll = zeros(1,1);
            dl = .0001;
            for r1 = 0:dl:1*dl+xyL(1)
                i = i+1;
                gAll(i,:) = [r1];
                [L(i),lnL(i)] = getL_P_background(N,m4a,m4b,fBa,fBb,...
                    r1,...
                    fD4a,fD4b);
            end
            idx = find(lnL==max(lnL));
            xyL = gAll(idx,:);
            lnLL(irun) = lnL(idx);
            clear dl n1 n2 n3 n4 i gAll L lnL idx
            if irun>1
                deltaL = abs(lnLL(irun)-lnLL(irun-1));
            end
        end
        deltaL =1;
        while deltaL>0
            irun = irun+1;
            i =0;
            gAll = zeros(1,1);
            dl = .00001;
            for r1 = 0:dl:1*dl+xyL(1)
                i = i+1;
                gAll(i,:) = [r1];
                [L(i),lnL(i)] = getL_P_background(N,m4a,m4b,fBa,fBb,...
                    r1,...
                    fD4a,fD4b);
            end
            idx = find(lnL==max(lnL));
            xyL = gAll(idx,:);
            lnLL(irun) = lnL(idx);
            clear dl n1 n2 n3 n4 i gAll L lnL idx
            if irun>1
                deltaL = abs(lnLL(irun)-lnLL(irun-1));
            end
        end
        PB4 = xyL;
        PBk = [PB1 PB2 PB3 PB4];
    end
    PB1 = mean(PBk(:,1));
    PB2 = mean(PBk(:,2));
    PB3 = mean(PBk(:,3));
    PB4 = mean(PBk(:,4));
end
clearvars -except T tc bint_E fname_dc fname_B fname_E...
    fD1a fD1b fD2a fD2b fD3a fD3b fD4a fD4b...
    RD1 RD2 RD3 RD4...
    hist_dc1 hist_dc2 hist_dc3 hist_dc4...
    PB1 PB2 PB3 PB4...
    fBa fBb...
    hist_B1 hist_B2 hist_B3 hist_B4 tag fname_params...
    yB1 yB2 yB3 yB4 xB simulated fileE folderE



if ispc
    mkdir(strcat("Analysis\",string(folderE),'\',string(fileE)))
    save(strcat("Analysis\",string(folderE),'\',string(fileE),"\BandD_tc",string(tc),".mat"))
    addpath(genpath("Analysis\"))
elseif ismac || isunix
    mkdir(strcat("Analysis/",string(folderE),'/',string(fileE)))
    save(strcat("Analysis/",string(folderE),'/',string(fileE),"/BandD_tc",string(tc),".mat"))
    addpath(genpath("Analysis/"))


end


%f = 1 old output logic
end


function f = binPhotons_tc(tc,bint,ssync,ddtime,ch)
sz = round(max(ssync)/bint);

i = 0;
j = 0;
n_bin = 1;
m4a = zeros(1,sz);
m4b = zeros(1,sz);
m5a = zeros(1,sz);
m5b = zeros(1,sz);
m6a = zeros(1,sz);
m6b = zeros(1,sz);
m7a = zeros(1,sz);
m7b = zeros(1,sz);
m2plus = zeros(1,sz);
dt_n = 0;
Mt_n=0;
ch_n=0;

leng = length(ssync);
for iphoton = 1:leng
    i = i+1;
    if ssync(iphoton)<bint*n_bin
        dt_n(i) = ddtime(iphoton);
        ch_n(i) = ch(iphoton);
        Mt_n(i) = ssync(iphoton);
    end
    if ssync(iphoton)>=bint*n_bin
        %organize photons from nth bin
        m2plus(n_bin) = length(Mt_n)-length(unique(Mt_n));
        for cnt = 1:length(dt_n)-3
            if Mt_n(cnt+1)==Mt_n(cnt)
                dt_n(cnt)=NaN;
                dt_n(cnt+1)=NaN;
                if Mt_n(cnt+2)==Mt_n(cnt)
                    dt_n(cnt+2)=NaN;
                    if Mt_n(cnt+3)==Mt_n(cnt)
                        dt_n(cnt+3)=NaN;
                    end
                end
            end
        end
        for cnt = 1:length(dt_n)
            if dt_n(cnt)<tc
                if ch_n(cnt)==2
                    m4a(n_bin) = m4a(n_bin)+1;
                end
                if ch_n(cnt)==3
                    m5a(n_bin) = m5a(n_bin)+1;
                end
                if ch_n(cnt)==4
                    m6a(n_bin) = m6a(n_bin)+1;
                end
                if ch_n(cnt)==5
                    m7a(n_bin) = m7a(n_bin)+1;
                end
            end
            
            if dt_n(cnt)>=tc
                if ch_n(cnt)==2
                    m4b(n_bin) = m4b(n_bin)+1;
                end
                if ch_n(cnt)==3
                    m5b(n_bin) = m5b(n_bin)+1;
                end
                if ch_n(cnt)==4
                    m6b(n_bin) = m6b(n_bin)+1;
                end
                if ch_n(cnt)==5
                    m7b(n_bin) = m7b(n_bin)+1;
                end
            end
        end
        
        n_bin=n_bin+1;
        i=0;
        dt_n = 0;
        Mt_n=0;
        Ch_n=0;
        
    end
end
f = horzcat(m4a',m5a',m6a',m7a',m4b',m5b',m6b',m7b',m2plus');
end


function f = get_Hist_tc(dt,tc,ch,chN)
j=0;

ddt(1)=NaN;
for i = 1:length(dt)
    if ch(i)==chN
        j = j+1;
        ddt(j)=dt(i);
    end
end
uhist= zeros(1,2);
if isnan(ddt(1))==1
    uhist(1,1)=0;
    uhist(1,2)=0;
end
if isnan(ddt(1))==0
    for mv = 1:length(ddt)
        if ddt(1,mv)>=tc
            uhist(1,2) = uhist(1,2)+1;
        end
        if ddt(1,mv)<tc
            uhist(1,1) = uhist(1,1)+1;
            
        end
    end
end
f = uhist;
end

function [L,lnL] = getL_P_background(N,m1a,m1b,pfBa,pfBb,...
    PB1,...
    pfD1a,pfD1b)

% background pixel 1 arrival
PB1_mat = [PB1; 1-PB1];
Bmat = [1;0];

%combine above with background arrival time
PfB_mat = [pfBa pfBb 1-pfBa-pfBb];
P1 = PB1_mat*PfB_mat;
P1 = reshape(P1',[],1);

Bmata=[1 0 0];
Bmatb = [0 1 0];

Bmata = Bmat*Bmata;
Bmatb = Bmat*Bmatb;
Bmata = reshape(Bmata',[],1);
Bmatb = reshape(Bmatb',[],1);

%detector 1 dark counts
pfD1 = [pfD1a pfD1b 1-pfD1a-pfD1b];
P2 = P1*pfD1;
P2 = reshape(P2',[],1);
PF = P2';

DC = ones(length(Bmata),1);
DCa = [1 0 0];
DCb = [0 1 0];
filler = [1 1 1];

DCa = DC*DCa;
DCb = DC*DCb;
DCa = reshape(DCa',[],1);
DCb = reshape(DCb',[],1);
Bmata = Bmata*filler;
Bmatb = Bmatb*filler;
Bmata = reshape(Bmata',[],1);
Bmatb = reshape(Bmatb',[],1);


%make key

%combine Bxx and dxx
leng = length(Bmata);
c1a = ones(leng,1);
c2a = c1a;
c3a = c1a;
c4a = c1a;
c1b = c1a;
c2b = c1a;
c3b = c1a;
c4b = c1a;

c1a = c1a.*Bmata+c1a.*DCa;
c1b = c1b.*Bmatb+c1b.*DCb;

for i = 1:leng
    if c1a(i)>0
        c1b(i)=0;
    end
end

Mob = horzcat(c1a,c1b);

Mob1 = Mob.*0;
for ic = 1:size(Mob,1)
    for jc = 1:size(Mob,2)
        if Mob(ic,jc)>=1
            Mob1(ic,jc)=1;
        end
    end
end

MobS = sum(Mob1,2);

key1a=zeros(leng,1);
key1b=zeros(leng,1);
keyNone=zeros(leng,1);

for i = 1:leng
    if MobS(i)==1
        if Mob1(i,1)==1
            key1a(i) = 1;
        end
        if Mob1(i,2)==1
            key1b(i) = 1;
        end
    end
    if MobS(i)==0
        keyNone(i) = 1;
    end
end

% sort for observables
%one photons detected
Pm1a = PF*key1a;
Pm1b = PF*key1b;
Pnone = PF*keyNone;


%% checkpoint
% Pnone=Pnone
% Ptest = 1-Pm1a-Pm1b


L = (Pm1a^m1a)*...
    (Pm1b^m1b)*...
    ((1-Pm1a-Pm1b)^(N-m1a-m1b));


lnL = m1a*log(Pm1a)+...
    m1b*log(Pm1b)+...
    (N-m1a-m1b)*log(1-Pm1a-Pm1b);

% L1 = 0;
end

function [t_hist,y_hist]=get_histogram(ddtime,ch,T,w,chN)
nbin= round(T/w);
y_hist = zeros(1,nbin);
t_hist = y_hist;
for i =1:nbin
    t_hist(1,i) = i*w;
end
for i = 1:nbin
    for j = 1:length(ddtime)
        if ch(j) == chN
        if ddtime(j)<=t_hist(i)
            if ddtime(j)>t_hist(i)-w
                y_hist(1,i) = y_hist(i)+1;
            end
        end
        end
    end
end

end
