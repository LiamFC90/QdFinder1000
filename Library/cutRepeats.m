function [ssync, cch, ddtime, ttotal] = cutRepeats(sync, ch, dtime, total)
%CUTREPEATS cuts out secondary arrivals that occur in the same pixel during
%the same pulse period because analysis assumes the dead time is equal to
%the pulse period and these sorts of events are impossible.

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

