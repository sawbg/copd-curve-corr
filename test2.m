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
    
    patdata(pcount,11) = pat.finalGold;
    patdata(pcount,12) = pat.pctEmph_Slicer;
    patdata(pcount,13) = pat.pctGasTrap_Slicer;
    
    try
        temp = scale(data.vol);
        [res gof] = fit((0:(length(temp)-1))',temp, fittype('exp2'));
        patdata(pcount,1) = res.a;
        patdata(pcount,2) = res.b;
        patdata(pcount,3) = res.c;
        patdata(pcount,4) = res.d;
        patdata(pcount,5) = abs(sqrt(gof.rsquare));
        
        [~,index] = max(data.flow);
        temp = scale(data.flow(index:end));
        [res gof] = fit((0:(length(temp)-1))',temp, fittype('poly3'));
        patdata(pcount,6) = res.p1;
        patdata(pcount,7) = res.p2;
        patdata(pcount,8) = res.p3;
        patdata(pcount,9) = res.p4;
        patdata(pcount,10) = abs(sqrt(gof.rsquare));
    catch err
        disp(['Error on patient ' num2str(p)]);
        disp(err.message);
        pcount = pcount - 1;
        pause;
    end
end
%%
fprintf(['Parameter\tEmp. Time-Volume\tEmp. Volume-Flow\t' ...
    'Gas Trap Time-Volume\tGas Trap Volume-Flow\n']);

fprintf('Coeff 1\t%.0f%%\t%.0f%%\t%.0f%%\t%.0f%%\n', ...
    pcorr(patdata(:,12),patdata(:,1)), ...
    pcorr(patdata(:,12),patdata(:,6)), ...
    pcorr(patdata(:,13),patdata(:,1)), ...
    pcorr(patdata(:,13),patdata(:,6)));

fprintf('Coeff 2\t%.0f%%\t%.0f%%\t%.0f%%\t%.0f%%\n', ...
    pcorr(patdata(:,12),patdata(:,2)), ...
    pcorr(patdata(:,12),patdata(:,7)), ...
    pcorr(patdata(:,13),patdata(:,2)), ...
    pcorr(patdata(:,13),patdata(:,7)));

fprintf('Coeff 3\t%.0f%%\t%.0f%%\t%.0f%%\t%.0f%%\n', ...
    pcorr(patdata(:,12),patdata(:,3)), ...
    pcorr(patdata(:,12),patdata(:,8)), ...
    pcorr(patdata(:,13),patdata(:,3)), ...
    pcorr(patdata(:,13),patdata(:,8)));

fprintf('Coeff 4\t%.0f%%\t%.0f%%\t%.0f%%\t%.0f%%\n', ...
    pcorr(patdata(:,12),patdata(:,4)), ...
    pcorr(patdata(:,12),patdata(:,9)), ...
    pcorr(patdata(:,13),patdata(:,4)), ...
    pcorr(patdata(:,13),patdata(:,9)));

fprintf('R-Value\t%.0f%%\t%.0f%%\t%.0f%%\t%.0f%%\n', ...
    pcorr(patdata(:,12),patdata(:,5)), ...
    pcorr(patdata(:,12),patdata(:,10)), ...
    pcorr(patdata(:,13),patdata(:,5)), ...
    pcorr(patdata(:,13),patdata(:,10)));