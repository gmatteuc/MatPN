close all
clear all
clc

% set some pars
pars = set_pars();
SF=pars.stimPars.SF;
TF=pars.stimPars.TF;
DIR=pars.stimPars.DIR;
numSF=pars.stimPars.numSF;
numTF=pars.stimPars.numTF;
numDIR=pars.stimPars.numDIR;
stimulustype=pars.stimPars.stimulustype;
numstimulustype=pars.stimPars.numstimulustype;

% load SPIKEMAT
S=load('SPIKEMAT_12_01_2016_b6.mat');

% initialize storage variables
tuning_curve=zeros(numDIR,size(S.SPIKEmean,2));
tuning_curve_error=zeros(numDIR,size(S.SPIKEmean,2));
tuning_matrix=zeros(numSF,numTF,size(S.SPIKEmean,2));
tuning_matrix_error=zeros(numSF,numTF,size(S.SPIKEmean,2));
pref_DIR=zeros(size(S.SPIKEmean,2),numSF,numTF,numstimulustype);
pref_SF=zeros(size(S.SPIKEmean,2),numstimulustype);
pref_TF=zeros(size(S.SPIKEmean,2),numstimulustype);
OSI=zeros(size(S.SPIKEmean,2),numSF,numTF,numstimulustype);
DSI=zeros(size(S.SPIKEmean,2),numSF,numTF,numstimulustype);
CIr=zeros(size(S.SPIKEmean,2),numSF,numTF,numstimulustype);
CIa=zeros(size(S.SPIKEmean,2),numSF,numTF,numstimulustype);


% get direction tuning for all neurons

for nn=1:size(S.SPIKEmean,2)
    tic
    for i=1:numSF
        for j=1:numTF
            for k=1:numstimulustype
                
                [ tuning_curve(:,nn,i,j,k), tuning_curve_error(:,nn,i,j,k), pref_DIR(nn,i,j,k), DSI(nn,i,j,k), OSI(nn,i,j,k), CIr(nn,i,j,k), CIa(nn,i,j,k) ] = get_direction_tuning( nn, SF(i), TF(j), stimulustype{k}, S);
                
            end
        end
    end
    toc
end

% get frequency tuning for all neurons
for nn=1:size(S.SPIKEmean,2)
    tic
    for k=1:numstimulustype
        
        [ tuning_matrix(:,:,nn,k), tuning_matrix_error(:,:,nn,k), pref_SF(nn,k), pref_TF(nn,k) ] = get_frequency_tuning( nn, stimulustype{k}, S);
        
    end
    toc
end
