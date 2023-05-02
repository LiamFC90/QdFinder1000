function f = binPhotons_tc(tc, bint, ssync, ddtime, ch)
%BINPHOTONS_TC Summary of this function goes here
%   Detailed explanation goes here
sz = round(max(ssync)/bint);

i = 0;
%j = 0;
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
Mt_n = 0;
ch_n = 0;

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

