function [ Stim ] = get_noise_movie_PN()

% set pars
pars = set_pars_PN();
stimuli_folder=pars.stim_folder;
addpath(stimuli_folder);

% load noise movie
movname='im_matrix_909_downsampled_MODULATED_LONG_10_medium.mat';
load(movname)
% rename noise movie storage variable
rvideo=S(:,:,1:1818,1:end); % NB: 1818 = number of rames in 1 min at 30 hz = size(S,3) = chunk length
clear S

% reshape noise movie
Stim=zeros(size(rvideo,3)*size(rvideo,4),size(rvideo,2)*size(rvideo,1));
for ll=1:size(rvideo,4) % loop over chunks
    for k=1:size(rvideo,3) % loop over frames
        % get current frame
        frm=rvideo(:,:,k,ll);
        % stack successive frames rowise
        Stim(k+(ll-1)*size(rvideo,3),:)=frm(:);
    end
end

end

