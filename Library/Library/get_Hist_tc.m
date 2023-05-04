function f = get_Hist_tc(dt,tc,ch,chN)
%GET_HIST_TC Summary of this function goes here
%   Detailed explanation goes here
j=0;
for i = 1:length(dt)
    if ch(i)==chN
        j = j+1;
        ddt(j)=dt(i);
    end
end
uhist= zeros(1,2);
for mv = 1:length(ddt)
    if ddt(1,mv)>=tc
        uhist(1,2) = uhist(1,2)+1;
    end
    if ddt(1,mv)<tc
        uhist(1,1) = uhist(1,1)+1;

    end
end
f = uhist;
end

