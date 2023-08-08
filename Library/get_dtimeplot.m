
function [plotvar] = get_dtimeplot(filepath,jobtype)
%GET_DTIMEPLOT Prints a histogram that shows decays for the data file
%   Detailed explanation goes here
temp = 5;
if string(jobtype) == 'rawptu'
    [output, ~]= Read_PTU_V1_Barelli_fast(filepath);
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
elseif string(jobtype) == 'sim'
    load(filepath,"dtime")
end
   
figure;
histogram(dtime)
set(gca,'YScale','log')
title('Decay Params')
ylabel('freq')
xlabel('time (ns)')
drawnow
figure;
histogram(sync)
title('Intensity')
ylabel('n')
xlabel('time(s)')
end

