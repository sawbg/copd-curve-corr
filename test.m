%clear all;
clc; format shortg;
%load('patient_vol_flow.mat');

patlen = length(patient);
pcount = 0;
skip = 1;
maxstage = 1;

pcorr = @(x,y) 100*abs(corr(x,y));

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
    
    patdata(pcount,1) = pat.finalGold;
    patdata(pcount,2) = pat.pctEmph_Slicer;
    patdata(pcount,3) = pat.pctGasTrap_Slicer;
    
    % Some improved results when data scaled
    % EXPLORE FURTHER
    [patdata(pcount,4) patdata(pcount,5)] ...
        = max(diff(scale(data.vol),2));
    [patdata(pcount,6) patdata(pcount,7)] ...
        = max(diff(scale(data.flow),2));
    
    patdata(pcount,8) = norm([patdata(pcount,4) patdata(pcount,5)]);
    patdata(pcount,9) = norm([patdata(pcount,6) patdata(pcount,7)]);
    
    patdata(pcount,10) ...
        = geomean([patdata(pcount,4) patdata(pcount,5)]);
    patdata(pcount,11) ...
        = geomean([patdata(pcount,6) patdata(pcount,7)]);
    
    try
        [~,index] = max(data.vol);
        temp = scale(data.vol(index:end));
        res = fit((0:(length(temp)-1))',temp, fittype('exp1'));
        patdata(pcount,12) = res.a;
        patdata(pcount,13) = res.b;
        
        [~,index] = max(data.flow);
        temp = scale(data.flow(index:end));
        [res gof] = fit((0:(length(temp)-1))',temp, fittype('exp1'));
        patdata(pcount,14) = res.a;
        patdata(pcount,15) = res.b;
    catch err
        disp(['Error on patient ' num2str(p)]);
        disp(err.message);
        pcount = pcount - 1;
    end
end

fprintf(['Parameter\tEmp. Time-Volume\tEmp. Volume-Flow\t' ...
    'Gas Trap Time-Volume\tGas Trap Volume-Flow\n']);

fprintf('Index\t%.0f%%\t%.0f%%\t%.0f%%\t%.0f%%\n', ...
    pcorr(patdata(:,2),patdata(:,5)), ...
    pcorr(patdata(:,2),patdata(:,7)), ...
    pcorr(patdata(:,3),patdata(:,5)), ...
    pcorr(patdata(:,3),patdata(:,7)));

fprintf('Value\t%.0f%%\t%.0f%%\t%.0f%%\t%.0f%%\n', ...
    pcorr(patdata(:,2),patdata(:,4)), ...
    pcorr(patdata(:,2),patdata(:,6)), ...
    pcorr(patdata(:,3),patdata(:,4)), ...
    pcorr(patdata(:,3),patdata(:,6)));

fprintf('Norm\t%.0f%%\t%.0f%%\t%.0f%%\t%.0f%%\n', ...
    pcorr(patdata(:,2),patdata(:,8)), ...
    pcorr(patdata(:,2),patdata(:,9)), ...
    pcorr(patdata(:,3),patdata(:,8)), ...
    pcorr(patdata(:,3),patdata(:,9)));

fprintf('Mean\t%.0f%%\t%.0f%%\t%.0f%%\t%.0f%%\n', ...
    pcorr(patdata(:,2),patdata(:,10)), ...
    pcorr(patdata(:,2),patdata(:,11)), ...
    pcorr(patdata(:,3),patdata(:,10)), ...
    pcorr(patdata(:,3),patdata(:,11)));

fprintf('Exp a\t%.0f%%\t%.0f%%\t%.0f%%\t%.0f%%\n', ...
    pcorr(patdata(:,2),patdata(:,12)), ...
    pcorr(patdata(:,2),patdata(:,14)), ...
    pcorr(patdata(:,3),patdata(:,12)), ...
    pcorr(patdata(:,3),patdata(:,14)));

fprintf('Exp b\t%.0f%%\t%.0f%%\t%.0f%%\t%.0f%%\n', ...
    pcorr(patdata(:,2),patdata(:,13)), ...
    pcorr(patdata(:,2),patdata(:,15)), ...
    pcorr(patdata(:,3),patdata(:,13)), ...
    pcorr(patdata(:,3),patdata(:,15)));