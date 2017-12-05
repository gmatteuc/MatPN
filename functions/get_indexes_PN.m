function [ indexes ] = get_indexes_PN( SF, TF, DIR, stimulustype  )

%[ indexes ] = get_indexes( SF, TF, DIR, stimulustype  )
%
%Retrive indexes (i.e. bitcodes) corresponding to a particular set of
%direction, spatial and temporal, frequency condition.
%
%Input a nonnumeric variable (such as 'all') in the input field for wich
%you don't want to do any selection. Usueful to get indexes necessary to
%slice SPIKEMAT in order to obtain tuning curves.

%--------------------------------------------------------------------------

% set pars  
pars = set_pars_PN();

switch stimulustype

    % choose stimulus type and load bitcodes lookup table
    case 'grating'
        % cut needed rows
        bcode=pars.stimPars.gratings_bitcodes;
        
        % choose stimulus type and load bitcodes lookup table
    case 'plaid'
        % cut needed rows
        bcode=pars.stimPars.plaids_bitcodes;
        
end

% if a TF constraint has been specified reduce the current bitcodes set to compy
%-------------------------------------------------------------------------%
if isnumeric(TF)
    
    switch TF
        case pars.stimPars.TF(1)
            % cut needed rows
            inds1=1:numel(pars.stimPars.DIR)*numel(pars.stimPars.SF);
        case pars.stimPars.TF(2)
            % cut needed rows
            inds1=(1:numel(pars.stimPars.DIR)*numel(pars.stimPars.SF))+ ...
                numel(pars.stimPars.DIR)*numel(pars.stimPars.SF);
    end
    bcode1=bcode(inds1);
else
    bcode1=bcode;
end
%-------------------------------------------------------------------------%


% if a SF constraint has been specified reduce the current bitcodes set to compy
%-------------------------------------------------------------------------%
if isnumeric(SF)
    switch SF
        case pars.stimPars.SF(1)
            % cut needed rows
            inds2=1:numel(pars.stimPars.DIR);
        case pars.stimPars.SF(2)
            % cut needed rows
            inds2=(1:numel(pars.stimPars.DIR))+ ...
                numel(pars.stimPars.DIR);
    end
    bcode2=bcode1(inds2);
else
    bcode2=bcode1;
end
%-------------------------------------------------------------------------%

% if a DIR constraint has been specified reduce the current bitcodes set to compy
%-------------------------------------------------------------------------%
if isnumeric(DIR)
    inds3 = find(pars.stimPars.DIR==DIR);
    bcode3=bcode2(inds3);
else
    bcode3=bcode2;
end
%-------------------------------------------------------------------------%

% needed indexes
indexes=bcode3;

end

