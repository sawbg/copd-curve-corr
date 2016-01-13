% Make sure required variables are in workspace
if exist('pdcell','var') == 0 || exist('pdcolnames','var') == 0
    error('No data!');
end

% iterate through all five simulations (for five max stages)
for k = 1:length(pdcell)
    file = fopen(['pat_data_max_stage_' num2str(k-1)]);
    patdata = pdcell{k};
    patids = patdata{1};
    patdata = patdata{2};
    pdsize = size(patdata);
    
    % write column names
    for m = 1:length(pdcolnames)
        fprintf(file, '%s', pdcolnames{m});
        
        if m < length(pdcolnames)
            fprintf(file, ',');
        else
            fprintf(file, '\n');
        end
    end
    
    % write patient data
    for m = 1:pdsize(1)
        fprintf(file, '%s,', patids(m));
        
        % write individual piece of data
        for n = 1:pdsize(2)
            fprintf(file, '%f', patdata(m,n));
            
            if n < pdsize(2)
                fprintf(file, ',');
            else
                fprintf(file, '\n');
            end
        end
    end
    
    fclose(file);
end