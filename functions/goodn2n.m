function [ n ] = goodn2n( area,goodn )

% load tuning analysis results and indexing
load('Indexing.mat')
% load neuron selection results
load(['selected_population_',area,'.mat'])
n=selectedsi{goodn}(1);

end

