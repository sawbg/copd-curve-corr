clear all;
loops = 100000;
nelm = 1000;
vals = zeros(1,loops);

for k=1:loops
    val(k) = corr(rand(nelm,1),rand(nelm,1));
    clc;
    disp(k);
end

val = sort(val);
hist(val,100);
average = mean(val)
deviation = std(val)
minimum = min(val)
maximum = max(val)