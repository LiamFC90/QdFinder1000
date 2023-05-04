function [mBcounts] = get_mBcounts(fname_B ,tc, bint, simulated)
%GET_MBCOUNTS Finds the expected mean backround and then creates a bin
%system identical to mcounts so background can be removed per bin. 
%   This method will break down if background radiation is not constant. 
if simulated ==0
    [output,~] =Read_PTU_V1_Barelli_fast(fname_B);% this function must be in the same directory as this progra
    MeasDesc_GlobalResolution = output.Headers.MeasDesc_GlobalResolution;%macrotime resolution in seconds
    MeasDesc_Resolution = output.MeasDesc_Resolution; %microt resolution in seconds
    ddtime = output.ph_dtime.*MeasDesc_Resolution*(10^9); %microt in nanoseconds
    ssync = output.ph_sync.*MeasDesc_GlobalResolution; % in seconds
    ch = output.ph_channel;
    total = length(ssync);
    if max(ch)==7
        ch = ch-2;
    end
    [ssync, ch, ddtime, ~] = cutRepeats(ssync,ch,ddtime,total);
end
if simulated==1
    load(fname_B,'ch','dtime','sync')
    ssync = sync;
    ddtime = dtime;
    ch = ch+1;

end
counts = binPhotons_tc(tc,bint,ssync,ddtime,ch);
% counts = counts(1:length(counts)-1,:);
mBcounts = mean(counts,1);
end

