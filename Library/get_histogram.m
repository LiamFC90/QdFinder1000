function [t_hist,y_hist] = get_histogram(ddtime,ch,T,w,chN)
%GET_HISTOGRAM Builds histogram of decays from ddtime
%   Should be used once per channel!!!! so an example would look like
%   get_histogram(dt_n,ch_n,T,.96,2). This looks at channel 2. .96 = w is
%   for 96ns which is a factor of 16ns so we don't get combing in the
%   histogram. 

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

