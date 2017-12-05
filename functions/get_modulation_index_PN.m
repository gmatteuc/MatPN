function [ modulation_index,modulation_index_bis ] = get_modulation_index_PN( nn, SF , TF , DIR, stimulustype, S)

%Comupute modulation index for neuron nn
%and stimulustype of spatial frequency SF and temporal frequency TF

%--------------------------------------------------------------------------

% set some pars
pars = set_pars_PN();
pre_time=pars.stimPars.pre_time;
stim_time=pars.stimPars.stim_time;
post_delay=pars.stimPars.post_delay;

% select bitcodes for the tuning curve
selected_idx = get_indexes_PN( SF, TF, DIR, stimulustype  );

% skip (output NaNs) parameters combinations for which ter is no bitcode
if isempty(selected_idx)
    modulation_index=0;
    modulation_index_bis=0;
else
    
    % get psth
    S_ts=[];
    for i=1:size(S.SPIKEtimes_bis,3)
        S_ts=[S_ts;S.SPIKEtimes_bis{ selected_idx, nn, i}-pre_time];
    end
    % take spikes falling into analysis window only
    S_ts=S_ts(S_ts<=(stim_time+post_delay));
    S_ts=S_ts(S_ts>=0);
    % set sampling frequency
    T_s=0.010;
    hedges=0:T_s:(stim_time+post_delay);
    % bin psth
    [psth]=hist(S_ts,hedges);
    
    % ---- compute F1z ----
    [ F1z, ~, ~, ~ ] = get_F1z( psth, TF, T_s );
    
    %  ---- compute F1/F0 ----
    [ F1F0, ~, ~, ~ ] = get_F1F0( psth, TF, T_s );
    
    % assign output
    modulation_index = F1z;
    modulation_index_bis = F1F0;
    
end

end

