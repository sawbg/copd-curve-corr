clc; format shortg;
if exist('patient','var') == 0, load('patient_vol_flow.mat'); end

patlen = length(patient);
pcount = 0;
skip = 1;

%pcorr = @(x,y) 100*abs(corr(x,y));
pcorr = @(x,y) corr(x,y);

for maxstage=0:4
    disp(['MAX STAGE: ' num2str(maxstage)]);
    disp('');
    
    for p=1:skip:patlen
        if mod(patlen-p,1000) == 0
            disp([num2str(patlen - p) ' remaining...']);
        end
        
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
            temp = scale(data.vol);
            [res gof] = fit((0:(length(temp)-1))',temp, fittype('exp2'));
            patdata(pcount,12) = res.a;
            patdata(pcount,13) = res.b;
            patdata(pcount,14) = res.c;
            patdata(pcount,15) = res.d;
            patdata(pcount,16) = abs(sqrt(gof.rsquare));
            
            [~,index] = max(data.flow);
            temp = scale(data.flow(index:end));
            [res gof] = fit((0:(length(temp)-1))',temp, fittype('poly3'));
            patdata(pcount,17) = res.p1;
            patdata(pcount,18) = res.p2;
            patdata(pcount,19) = res.p3;
            patdata(pcount,20) = res.p4;
            patdata(pcount,21) = abs(sqrt(gof.rsquare));
        catch err
            disp(['Error on patient ' num2str(p)]);
            disp(err.message);
            pcount = pcount - 1;
            pause;
        end
    end
    
    disp('');
    
    fprintf(['Parameter\tEmp. Time-Volume\tEmp. Volume-Flow\t' ...
        'Gas Trap Time-Volume\tGas Trap Volume-Flow\n']);
    
    fprintf('Index\t%.3f\t%.3f\t%.3f\t%.3f\n', ...
        pcorr(patdata(:,2),patdata(:,5)), ...
        pcorr(patdata(:,2),patdata(:,7)), ...
        pcorr(patdata(:,3),patdata(:,5)), ...
        pcorr(patdata(:,3),patdata(:,7)));
    
    fprintf('Value\t%.3f\t%.3f\t%.3f\t%.3f\n', ...
        pcorr(patdata(:,2),patdata(:,4)), ...
        pcorr(patdata(:,2),patdata(:,6)), ...
        pcorr(patdata(:,3),patdata(:,4)), ...
        pcorr(patdata(:,3),patdata(:,6)));
    
    fprintf('Norm\t%.3f\t%.3f\t%.3f\t%.3f\n', ...
        pcorr(patdata(:,2),patdata(:,8)), ...
        pcorr(patdata(:,2),patdata(:,9)), ...
        pcorr(patdata(:,3),patdata(:,8)), ...
        pcorr(patdata(:,3),patdata(:,9)));
    
    fprintf('Mean\t%.3f\t%.3f\t%.3f\t%.3f\n', ...
        pcorr(patdata(:,2),patdata(:,10)), ...
        pcorr(patdata(:,2),patdata(:,11)), ...
        pcorr(patdata(:,3),patdata(:,10)), ...
        pcorr(patdata(:,3),patdata(:,11)));
    
    fprintf('Coeff a\t%.3f\t%.3f\t%.3f\t%.3f\n', ...
        pcorr(patdata(:,2),patdata(:,12)), ...
        pcorr(patdata(:,2),patdata(:,17)), ...
        pcorr(patdata(:,3),patdata(:,12)), ...
        pcorr(patdata(:,3),patdata(:,17)));
    
    fprintf('Coeff b\t%.3f\t%.3f\t%.3f\t%.3f\n', ...
        pcorr(patdata(:,2),patdata(:,13)), ...
        pcorr(patdata(:,2),patdata(:,18)), ...
        pcorr(patdata(:,3),patdata(:,13)), ...
        pcorr(patdata(:,3),patdata(:,18)));
    
    fprintf('Coeff c\t%.3f\t%.3f\t%.3f\t%.3f\n', ...
        pcorr(patdata(:,2),patdata(:,14)), ...
        pcorr(patdata(:,2),patdata(:,19)), ...
        pcorr(patdata(:,3),patdata(:,14)), ...
        pcorr(patdata(:,3),patdata(:,19)));
    
    fprintf('Coeff d\t%.3f\t%.3f\t%.3f\t%.3f\n', ...
        pcorr(patdata(:,2),patdata(:,15)), ...
        pcorr(patdata(:,2),patdata(:,20)), ...
        pcorr(patdata(:,3),patdata(:,15)), ...
        pcorr(patdata(:,3),patdata(:,20)));
end