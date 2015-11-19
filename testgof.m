%%clear all;
clc; format shortg;
%%load('patient_vol_flow.mat');

patlen = length(patient);
pcount = 0;
skip = 10;
maxstage = 0;

for p=1:skip:patlen
    clc; disp([num2str(patlen - p) ' remaining...']);
    pat = patient(p);
    trialnum = 1;
    next = false;
    
    while trialnum <= 6
        data = patient(p).trial(trialnum);
        [~,volind] = max(data.vol);
        [~,flowind] = max(data.flow);
        vollen = length(data.vol);
        flowlen = length(data.flow);
        
        if volind ~= vollen && flowind ~= flowlen
            break;
        elseif trialnum <= 6
            trialnum = trialnum + 1;
        else
            next = true;
            break;
        end
    end
    
    if isnan(pat.pctEmph_Slicer) || isnan(pat.pctGasTrap_Slicer) ...
            || pat.finalGold > maxstage || next
        continue;
    else
        pcount = pcount + 1;
    end
    
    try
        temp = scale(data.vol);
        [res gof] = fit((0:(length(temp)-1))',temp, fittype('exp2'));
        patdata(pcount,1) = abs(sqrt(gof.rsquare));
        
        [~,index] = max(data.flow);
        temp = scale(data.flow(index:end));
        [res gof] = fit((0:(length(temp)-1))',temp, fittype('poly3'));
        patdata(pcount,2) = abs(sqrt(gof.rsquare));
    catch err
        disp(['Error on patient ' num2str(p)]);
        disp(err.message);
        pcount = pcount - 1;
    end
end