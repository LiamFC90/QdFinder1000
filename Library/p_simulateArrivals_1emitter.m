


T = 400; %laser pulse ns
t = 5;%time in seconds that sim runs for WAS 5
sigmai = 180;%in nm PSF size
sigmaiy = 180;

mui =  -50;    %pos of emmiter
muiy = 100;

Ri = 1/50; %lifetimes of emmiter 1 off 50 ns
RB = 1/3;  %lifetime Background 1 off 3ns

Pi = .05;   %rate of emmision at ordinary levels EFFECTS PERFORMNCE WAS .05

RD1 = 5e-7; %rate of dark counts GOOD LEVEL NO TOUCHY   WAS 5e-7
RD2 = RD1;
RD3 = RD1;
RD4 = RD1;

PB1 = .005;  %prob of backround counts GOOD LEVEL NO TOUCHY   WAS .005
PB2 = PB1;
PB3 = PB1;
PB4 = PB1;

%program begins
Sync = NaN;
Dtime = NaN;
Ch = NaN;
N = t/(T*10^-9);%number of periods in time t
iT = 0;
dtime = NaN.*zeros(round(N*2),1);
Mt = NaN.*zeros(round(N*2),1);
ch = NaN.*zeros(round(N*2),1);

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
% cladding = 0;
W = w+cladding;
leng = W/2; %with cladding
x1 = 1*leng; y1 = 1*leng;
x2 = -1*leng; y2 = 1*leng;
x3 = -1*leng; y3 = -1*leng;
x4 = 1*leng; y4 = -1*leng;

for cnt = 1:N
    clear dtmat posmat
    dtmat = NaN.*zeros(1,10);
    posmat = NaN.*zeros(1,10);
    %dark count contribution
    dt = exprnd(1/RD1);
    if dt<T
        dtmat(6) = dt;
        posmat(6) = 1;
    end
    clear dt
    dt = exprnd(1/RD2);
    if dt<T
        dtmat(7) = dt;
        posmat(7)=2;
    end
    clear dt
    dt = exprnd(1/RD3);
    if dt<T
        dtmat(8) = dt;
        posmat(8)=3;
    end
    clear dt
    dt = exprnd(1/RD4);
    if dt<T
        dtmat(9) = dt;
        posmat(9)=4;
    end
    clear dt
    
    %emitter i
    ic = rand(1);
    if ic<Pi
        
        dt =exprnd(1/Ri);
        if dt<T
            dtmat(5)=dt;
            X = normrnd(mui,sigmai);
            Y = normrnd(muiy,sigmaiy);
            if abs(X-x1)<w/2
                if abs(Y-y1)<w/2
                    posmat(5) =1;
                end
            end
            if abs(X-x2)<w/2
                if abs(Y-y2)<w/2
                    posmat(5)=2;
                end
            end
            if abs(X-x3)<w/2
                if abs(Y-y3)<w/2
                    posmat(5)=3;
                end
            end
            if abs(X-x4)<w/2
                if abs(Y-y4)<w/2
                    posmat(5) = 4;
                end
            end
        end
        clear dt
    end
    

    
    %background
    ic = rand(1);
    if ic<PB1
        dt =exprnd(1/RB);
        if dt<T
            dtmat(1)=dt;
            posmat(1)=1;
        end
        clear dt
    end
    ic = rand(1);
    if ic<PB2
        dt =exprnd(1/RB);
        if dt<T
            dtmat(2)=dt;
            posmat(2)=2;
        end
        clear dt
    end
    ic = rand(1);
    if ic<PB3
        dt =exprnd(1/RB);     %Pixel Downtime
        if dt<T
            dtmat(3)=dt;
            posmat(3)=3;
        end
        clear dt
    end
    ic = rand(1);
    if ic<PB4
        dt =exprnd(1/RB);
        if dt<T
            dtmat(4)=dt;
            posmat(4)=4;
        end
        clear dt
    end
    
    test = isnan(max(posmat));
    if test==0%there is at least one count detected
        for i = 1:10
            if isnan(posmat(i))==0
                
                iT = iT+1;
                
                dtime(iT,1) = dtmat(i);
                Mt(iT,1) =cnt*T/(10^9);
                ch(iT,1) = posmat(i);
            end
        end
    end
end
Mt = Mt(~isnan(Mt));
dtime = dtime(~isnan(dtime));
ch = ch(~isnan(ch));

Sync = vertcat(Sync,Mt);
Dtime = vertcat(Dtime,dtime);
Ch = vertcat(Ch,ch);

clear Mt dtime ch


Sync = Sync(~isnan(Sync));
Dtime = Dtime(~isnan(Dtime));
Ch = Ch(~isnan(Ch));

[sync, ch, dtime, total] = cutRepeats(Sync,Ch,Dtime,length(Dtime));
[sync, ch, dtime, total] = cutRepeats(sync,ch,dtime,length(dtime));

[t_hist,y_hist] = get_histogram(dtime,ch,T,0.96,2);
figure; plot(t_hist,y_hist)

load train
sound(y,Fs)

save(saveName);
toc

function [ssync, cch, ddtime, ttotal] = cutRepeats(sync,ch,dtime,total)
%cuts out secondary arrivals that occur in the same pixel during the same
%pulse period because analysis assumes the dead time is equal to the pulse
%period and these sorts of events are impossible.
for i = 1:total-1
    if sync(i) == sync(i+1)
        if ch(i) == ch(i+1)
            if dtime(i)>dtime(i+1)
                sync(i) = NaN;
                ch(i) = NaN;
                dtime(i) = NaN;
            end
            if dtime(i)<dtime(i+1)
                sync(i+1) = NaN;
                ch(i+1) = NaN;
                dtime(i+1) = NaN;
            end
        end
    end
end
sync = sync(~isnan(sync));
dtime = dtime(~isnan(dtime));
ch = ch(~isnan(ch));
totalp = length(sync);
for i = 1:totalp-1
    if sync(i) == sync(i)+1
        if ch(i) == ch(i+1)
            if dtime(i)>dtime(i+1)
                sync(i) = NaN;
                ch(i) = NaN;
                dtime(i) = NaN;
            end
            if dtime(i)<dtime(i+1)
                sync(i+1) = NaN;
                ch(i+1) = NaN;
                dtime(i+1) = NaN;
            end
        end
    end
end
ssync = sync(~isnan(sync));
ddtime = dtime(~isnan(dtime));
cch = ch(~isnan(ch));
ttotal = length(ssync);
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
        if ch(j) ==chN
            if ddtime(j)<=t_hist(i)
                if ddtime(j)>t_hist(i)-w
                    y_hist(1,i) = y_hist(i)+1;
                end
            end
        end
    end
end

end
