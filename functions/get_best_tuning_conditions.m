function [ max_indval, best_SF, best_TF, sorted_idx ] = get_best_tuning_conditions( indval )

%[ max_indval, best_SF, best_TF, sorted_idx ] = get_best_tuning_conditions( indval )
%
%Find maximum tuning conditions for every neuron and sort them 
% according to tuning level

%-----------------------------------------------------------------------

load('Indexing.mat');

SI=indval;

% find maximum SI conditions (SF_ind e SF_ind ontain max SI frequencies for every neuron)
SF_ind=zeros(1,size(M,1));
TF_ind=zeros(1,size(M,1));
max_SI=zeros(1,size(M,1));
for nn=1:size(M,1)
temp_SI=squeeze(SI(nn,:,:,1))+0.0001*rand(size(squeeze(SI(nn,:,:,1))));
max_temp_SI=max(temp_SI(:)); 
idx_max_temp_SI=find(temp_SI(:)==max_temp_SI);
[SF_ind(nn),TF_ind(nn)]=ind2sub(size(temp_SI),idx_max_temp_SI);
max_SI(nn)=SI(nn,SF_ind(nn),TF_ind(nn),1);
end

% sort maximum SI conditions
[~,idx_max_SI] = sort(max_SI,'descend');

best_TF=TF_ind;
best_SF=SF_ind;
max_indval=max_SI;
sorted_idx=idx_max_SI;

end

