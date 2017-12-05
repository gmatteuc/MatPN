function [ tuning_curve, c_tuning_curve,d_tuning_curve,...
    tuning_curve_matrix, pref_DIR,...
    OSI, DSI, DI, c_pref_DIR,...
    c_OSI, c_DSI, c_DI, ...
    d_pref_DIR, d_OSI, d_DSI, d_DI, ...
    sigperm_counts, muperm_counts,...
    response_counts, spontaneous_counts,...
    sigperm_fr, muperm_fr, response_fr,...
    spontaneous_fr, tcorr, sigtrials, tuning_matrix, tuning_matrix_error, ...
    pref_SF, pref_TF, c_tuning_matrix, c_tuning_matrix_error, ...
    c_pref_SF, c_pref_TF,...
    d_tuning_matrix, d_tuning_matrix_error, ...
    d_pref_SF, d_pref_TF, Zp, Zc, Rp, Rc, PI, ...
    c_Zp, c_Zc, c_Rp, c_Rc, c_PI, ...
    d_Zp, d_Zc, d_Rp, d_Rc, d_PI, ...
    modulation_index, modulation_index_bis, tuning_curve_z, basal_fr, basal_fr_std ] = analyze_tuning_PN( S )

% [ tuning_curve, c_tuning_curve,d_tuning_curve,...
%     tuning_curve_matrix, pref_DIR,...
%     OSI, DSI, DI, c_pref_DIR,...
%     c_OSI, c_DSI, c_DI, ...
%     d_pref_DIR, d_OSI, d_DSI, d_DI, ...
%     sigperm_counts, muperm_counts,...
%     response_counts, spontaneous_counts,...
%     sigperm_fr, muperm_fr, response_fr,...
%     spontaneous_fr, tcorr, sigtrials, tuning_matrix, tuning_matrix_error, ...
%     pref_SF, pref_TF, c_tuning_matrix, c_tuning_matrix_error, ...
%     c_pref_SF, c_pref_TF,...
%     d_tuning_matrix, d_tuning_matrix_error, ...
%     d_pref_SF, d_pref_TF, d_Zp, d_Zc, d_Rp, d_Rc, d_PI, ...
%     modulation_index, modulation_index_bis, d_tuning_curve_z, basal_fr, basal_fr_std ] = analyze_tuning_PN( S )
%
%
%Perform complete direction, spatial and tempral frequency tuning analysis
%on the input SPIKEMAT:
%
% tuning_curve(:,nn,i,j,k) tuning_curve_error(:,nn,i,j,k) tuning_matrix(:,:,nn,k)
% tuning_matrix_error(:,:,nn,k) pref_DIR(nn,i,j,k) pref_SF(nn,k) pref_TF(nn,k)
% OSI(nn,i,j,k) DSI(nn,i,j,k) DCIr(nn,i,j,k) DCIa(nn,i,j,k) OCIr(nn,i,j,k)
% OCIa(nn,i,j,k) diff_grpl(nn) Rp(nn) Rc(nn) modulation_index(nn,i,j)
%
% with nn= neuron number, i=SF index, j=TF index, k=stimulustype index
%--------------------------------------------------------------------------

% set some pars
pars = set_pars_PN();
SF=pars.stimPars.SF;
TF=pars.stimPars.TF;
DIR=pars.stimPars.DIR;
numSF=pars.stimPars.numSF;
numTF=pars.stimPars.numTF;
numDIR=pars.stimPars.numDIR;
numtrial=pars.stimPars.ntrials;
stimulustype=pars.stimPars.stimulustype;
numstimulustype=pars.stimPars.numstimulustype;

% initialize storage variables
tuning_curve=zeros(numDIR,size(S.SPIKEmean,2),numSF,numTF,numstimulustype);
c_tuning_curve=zeros(numDIR,size(S.SPIKEmean,2),numSF,numTF,numstimulustype);
d_tuning_curve=zeros(numDIR,size(S.SPIKEmean,2),numSF,numTF,numstimulustype);
tuning_curve_z=zeros(numDIR,size(S.SPIKEmean,2),numSF,numTF,numstimulustype);
basal_fr=zeros(size(S.SPIKEmean,2),1);
basal_fr_std=zeros(size(S.SPIKEmean,2),1);
tuning_curve_matrix=zeros(numDIR,numtrial,size(S.SPIKEmean,2),numSF,numTF,numstimulustype);
pref_DIR=zeros(size(S.SPIKEmean,2),numSF,numTF,numstimulustype);
c_pref_DIR=zeros(size(S.SPIKEmean,2),numSF,numTF,numstimulustype);
d_pref_DIR=zeros(size(S.SPIKEmean,2),numSF,numTF,numstimulustype);
OSI=zeros(size(S.SPIKEmean,2),numSF,numTF,numstimulustype);
DSI=zeros(size(S.SPIKEmean,2),numSF,numTF,numstimulustype);
DI=zeros(size(S.SPIKEmean,2),numSF,numTF,numstimulustype);
c_OSI=zeros(size(S.SPIKEmean,2),numSF,numTF,numstimulustype);
c_DSI=zeros(size(S.SPIKEmean,2),numSF,numTF,numstimulustype);
c_DI=zeros(size(S.SPIKEmean,2),numSF,numTF,numstimulustype);
d_OSI=zeros(size(S.SPIKEmean,2),numSF,numTF,numstimulustype);
d_DSI=zeros(size(S.SPIKEmean,2),numSF,numTF,numstimulustype);
d_DI=zeros(size(S.SPIKEmean,2),numSF,numTF,numstimulustype);
sigperm_counts=zeros(numDIR,numtrial,size(S.SPIKEmean,2),numSF,numTF,numstimulustype);
muperm_counts=zeros(numDIR,numtrial,size(S.SPIKEmean,2),numSF,numTF,numstimulustype);
response_counts=zeros(numDIR,numtrial,size(S.SPIKEmean,2),numSF,numTF,numstimulustype);
spontaneous_counts=zeros(numDIR,numtrial,size(S.SPIKEmean,2),numSF,numTF,numstimulustype);
sigperm_fr=zeros(numDIR,size(S.SPIKEmean,2),numSF,numTF,numstimulustype);
muperm_fr=zeros(numDIR,size(S.SPIKEmean,2),numSF,numTF,numstimulustype);
response_fr=zeros(numDIR,size(S.SPIKEmean,2),numSF,numTF,numstimulustype);
spontaneous_fr=zeros(numDIR,size(S.SPIKEmean,2),numSF,numTF,numstimulustype);
tcorr=zeros(size(S.SPIKEmean,2),numSF,numTF,numstimulustype);
sigtrials =zeros(numtrial,size(S.SPIKEmean,2),numSF,numTF,numstimulustype);
tuning_matrix=zeros(numSF,numTF,size(S.SPIKEmean,2));
tuning_matrix_error=zeros(numSF,numTF,size(S.SPIKEmean,2));
pref_SF=zeros(size(S.SPIKEmean,2),numstimulustype);
pref_TF=zeros(size(S.SPIKEmean,2),numstimulustype);
c_tuning_matrix=zeros(numSF,numTF,size(S.SPIKEmean,2));
c_tuning_matrix_error=zeros(numSF,numTF,size(S.SPIKEmean,2));
c_pref_SF=zeros(size(S.SPIKEmean,2),numstimulustype);
c_pref_TF=zeros(size(S.SPIKEmean,2),numstimulustype);
d_tuning_matrix=zeros(numSF,numTF,size(S.SPIKEmean,2));
d_tuning_matrix_error=zeros(numSF,numTF,size(S.SPIKEmean,2));
d_pref_SF=zeros(size(S.SPIKEmean,2),numstimulustype);
d_pref_TF=zeros(size(S.SPIKEmean,2),numstimulustype);
Rp=zeros(size(S.SPIKEmean,2),numSF,numTF);
Rc=zeros(size(S.SPIKEmean,2),numSF,numTF);
Zp=zeros(size(S.SPIKEmean,2),numSF,numTF);
Zc=zeros(size(S.SPIKEmean,2),numSF,numTF);
PI=zeros(size(S.SPIKEmean,2),numSF,numTF);
c_Rp=zeros(size(S.SPIKEmean,2),numSF,numTF);
c_Rc=zeros(size(S.SPIKEmean,2),numSF,numTF);
c_Zp=zeros(size(S.SPIKEmean,2),numSF,numTF);
c_Zc=zeros(size(S.SPIKEmean,2),numSF,numTF);
c_PI=zeros(size(S.SPIKEmean,2),numSF,numTF);
d_Rp=zeros(size(S.SPIKEmean,2),numSF,numTF);
d_Rc=zeros(size(S.SPIKEmean,2),numSF,numTF);
d_Zp=zeros(size(S.SPIKEmean,2),numSF,numTF);
d_Zc=zeros(size(S.SPIKEmean,2),numSF,numTF);
d_PI=zeros(size(S.SPIKEmean,2),numSF,numTF);
modulation_index=zeros(numDIR,size(S.SPIKEmean,2),numSF,numTF,numstimulustype);
modulation_index_bis=zeros(numDIR,size(S.SPIKEmean,2),numSF,numTF,numstimulustype);

% get direction tuning for all neurons
for nn=1:size(S.SPIKEmean,2)
    tic
    for i=1:numSF
        for j=1:numTF
            for k=1:numstimulustype
                
                [ tuning_curve(:,nn,i,j,k), c_tuning_curve(:,nn,i,j,k),...
                    tuning_curve_matrix(:,:,nn,i,j,k), pref_DIR(nn,i,j,k),...
                    OSI(nn,i,j,k), DSI(nn,i,j,k), DI(nn,i,j,k), c_pref_DIR(nn,i,j,k),...
                    c_OSI(nn,i,j,k), c_DSI(nn,i,j,k), c_DI(nn,i,j,k), ...
                    sigperm_counts(:,:,nn,i,j,k), muperm_counts(:,:,nn,i,j,k),...
                    response_counts(:,:,nn,i,j,k), spontaneous_counts(:,:,nn,i,j,k),...
                    sigperm_fr(:,nn,i,j,k), muperm_fr(:,nn,i,j,k), response_fr(:,nn,i,j,k),...
                    spontaneous_fr(:,nn,i,j,k), tcorr(nn,i,j,k), sigtrials(:,nn,i,j,k)...
                    ] = get_direction_tuning_PN_NEW( nn, SF(i), TF(j), stimulustype{k}, S);
                
            end
        end
    end
    toc
end

% compute new corrected tcs
for nn=1:size(response_fr,2)
    
    s_fr=spontaneous_fr(:,nn,:,:,:);
    basal_fr(nn)=mean(s_fr(:));
    basal_fr_std(nn)=std(s_fr(:));
    d_tuning_curve(:,nn,:,:,:) = max(0,response_fr(:,nn,:,:,:)-basal_fr(nn));
    tuning_curve_z(:,nn,:,:,:) = (response_fr(:,nn,:,:,:)-basal_fr(nn))/basal_fr_std(nn);

end
for nn=1:size(S.SPIKEmean,2)
    tic
    for i=1:numSF
        for j=1:numTF
            for k=1:numstimulustype
                
                % on "raw" tuning curves
                [ d_OSI(nn,i,j,k),d_DSI(nn,i,j,k),d_DI(nn,i,j,k),d_pref_DIR(nn,i,j,k),~  ] = compute_SIs( d_tuning_curve(:,nn,i,j,k) );
            
            end
        end
    end
    toc
end

% get frequency tuning for all neurons
for nn=1:size(S.SPIKEmean,2)
    tic
    for k=1:numstimulustype
        
        [ tuning_matrix(:,:,nn,k), tuning_matrix_error(:,:,nn,k), pref_SF(nn,k), pref_TF(nn,k) ] = get_frequency_tuning_PN( nn, k, tuning_curve, sigperm_fr);
        [ c_tuning_matrix(:,:,nn,k), c_tuning_matrix_error(:,:,nn,k), c_pref_SF(nn,k), c_pref_TF(nn,k) ] = get_frequency_tuning_PN( nn, k, c_tuning_curve, sigperm_fr);
        [ d_tuning_matrix(:,:,nn,k), d_tuning_matrix_error(:,:,nn,k), d_pref_SF(nn,k), d_pref_TF(nn,k) ] = get_frequency_tuning_PN( nn, k, c_tuning_curve, sigperm_fr);

    end
    toc
end

% get stimulus tuning and patterness analysis for all neurons
for nn=1:size(S.SPIKEmean,2)
    tic
    for i=1:numSF
        for j=1:numTF
            
            [ Zp(nn,i,j),Zc(nn,i,j),Rp(nn,i,j),Rc(nn,i,j),PI(nn,i,j) ] = get_plaid_analysis_PN( nn, SF(i), TF(j), tuning_curve);
            [ c_Zp(nn,i,j),c_Zc(nn,i,j),c_Rp(nn,i,j),c_Rc(nn,i,j),c_PI(nn,i,j) ] = get_plaid_analysis_PN( nn, SF(i), TF(j), c_tuning_curve);
            [ d_Zp(nn,i,j),d_Zc(nn,i,j),d_Rp(nn,i,j),d_Rc(nn,i,j),d_PI(nn,i,j) ] = get_plaid_analysis_PN( nn, SF(i), TF(j), c_tuning_curve);
            
        end
    end
    toc
end

% get modulation index for all neurons
for nn=1:size(S.SPIKEtimes,2)
    tic
    for i=1:numSF
        for j=1:numTF
            for h=1:numDIR
                for k=1:numstimulustype
                    
                    [ modulation_index(h,nn,i,j,k),modulation_index_bis(h,nn,i,j,k)] = get_modulation_index_PN( nn, SF(i), TF(j), DIR(h), stimulustype{k}, S);
                    
                end
            end
        end
    end
    toc
end

end