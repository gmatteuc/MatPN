function [p] = quantile_test(distr1, distr2, nboot, frac_samp,quant)

% calculates p value in quantile comparison with bootstrap.
% distr1: first data distribution
% distr2: second data distribution
% nboot: number of bootstrap repetitions (use 2000)
% frac_sample: fration of data used at every bootsrap repetition (use 0.5)
% quant: quantile to compare (probability)


theta = zeros(1,nboot);

for i = 1:nboot

idx1 = randperm(length(distr1));
idx2 = randperm(length(distr2));


s1 = distr1(idx1(1:round(frac_samp*end)));
s2 = distr2(idx2(1:round(frac_samp*end)));

theta(i) = qdiff(s1,s2,quant);

end

delta_vect = sort(theta);
% plot(delta_vect)

min_delta = sum(theta<0);
zero_delta = sum(theta==0);

p = (min_delta+0.5*zero_delta)/nboot;

end


