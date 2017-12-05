function [ goodn ] = n2goodn( area, n )

% load tuning analysis results and indexing
load('Indexing.mat')
% load neuron selection results
load(['selected_population_',area,'.mat'])
nlist=cell2mat(selectedsi);
nlist=nlist(:,1);
goodn=find(n==nlist);

end

