function [ tuning_matrix, tuning_matrix_error, pref_SF, pref_TF ] = get_frequency_tuning_PN( nn, stimidx, tuning_curve, sigperm_fr)

%[ tuning_matrix, tuning_matrix_error, pref_SF, pref_TF ] = get_frequency_tuning_PN( nn, stimulustype, tuning_curve, sigperm_fr)
%
%Comupute tuning curve and parameters (selectivity indices) for neuron nn 
%and stimulustype of spatial frequency SF and temporal frequency TF
%

%--------------------------------------------------------------------------

% set some pars
pars = set_pars_PN();
SF=pars.stimPars.SF;
TF=pars.stimPars.TF;
numSF=pars.stimPars.numSF;
numTF=pars.stimPars.numTF;
    
    % get desired tuning curve
    tuning_matrix=zeros(numSF,numTF);
    tuning_matrix_error=zeros(numSF,numTF);
    % slice for selected bitcodes
    for i=1:numSF
        for j=1:numTF
              [tuning_matrix(i,j),midx]=randmax(tuning_curve(:,nn,i,j,stimidx));
              tuning_matrix_error(i,j)=sigperm_fr(midx,nn,i,j,stimidx);
        end
    end
    
    % find preferred direction
    [~,tempidx]=randmax(tuning_matrix(:));
    [rowm, colm] = ind2sub(size(tuning_matrix),tempidx);
    pref_SF=SF(rowm);
    pref_TF=TF(colm);
    
end

