thresh = 0.99;

timevolcorr = 0;
volflowcorr = 0;

for k=1:pcount
    if abs(patdata(k,1)) >= thresh, timevolcorr = timevolcorr + 1; end
    if abs(patdata(k,2)) >= thresh, volflowcorr = volflowcorr+ 1; end
end

timevolover = timevolcorr / pcount * 100
volflowover = volflowcorr / pcount * 100