function [] = spike_triggered_average_gaussian_PN_MEDIUM( sind, bind )

% [] = spike_triggered_average_gaussian_PN_MEDIUM( sind, bind )
% perform sta analysis of given block in given session
% and save results in the appropriate folders
%--------------------------------------------------------------------------

% set pars
pars=set_pars_PN();

% retrive session name
load('Indexing.mat')
session_list=pars.listSessions(1,:);
current_session = session_list{sind};
ssidx = (M(:,1)==sind); %#ok<NODEF>
blocks = unique(M(ssidx,2));
current_block=blocks(bind);
session_name=[current_session,'_b',num2str(current_block)];

% reduce indexing matrix to this session only
M_red=M(M(:,1)==sind,:);
M_red=M_red(M_red(:,2)==current_block,:);

% add useful folders
home_folder=pars.processed_data_folder;
fprintf(['\nLoading ','NOISEMAT_',session_name,'.mat\n'])
data_file=[home_folder,filesep,'NOISEMAT_',session_name,'.mat'];
load(data_file', 'psth_by_frames')
code_folder=pars.code_folder;
addpath(code_folder);
out_folder=[home_folder,filesep,'STA_results_',session_name];
mkdir(out_folder);

% set max pre-spike lag
maxlag=pars.STA_depth;
fwidth=pars.STA_width;
fheight=pars.STA_height;
% set regularization parameter
ridgeparam = pars.STA_ridgeparam;

% loop over neurons
numstart=1;
Nneu=size(psth_by_frames,2);
for neuronum = numstart:Nneu
    
    % skip neuron if already analyzed
    fnam=[out_folder,'/neuron_',num2str(neuronum),'_medium_results'];
    if exist(fnam,'dir')==7
        fprintf(['neuron ',num2str(neuronum),' already analyzed, skipping... \n'])
        numstart=neuronum+1;
    else
        
        if M_red(neuronum,8)==0 % process good units only
            
            % initialize variables in witch store results at different lags
            Dstafr=zeros(fheight,fwidth,maxlag);
            Dwstafr=zeros(fheight,fwidth,maxlag);
            DZwstafr=zeros(fheight,fwidth,maxlag);
            DZstafr=zeros(fheight,fwidth,maxlag);
            Drawcov=zeros(fheight*fwidth,fheight*fwidth,maxlag);
            Dneuronum=neuronum; %#ok<NASGU>
            
            % load stimulus
            movname='im_matrix_909_downsampled_MODULATED_LONG_10_medium.mat';
            load(movname)
            rvideo=S; % once you used t change it before running if you wnt to exclude some trials
            clear S
            
            % reformat stimulus movie
            Stim=zeros(size(rvideo,3)*size(rvideo,4),size(rvideo,2)*size(rvideo,1));
            for ll=1:size(rvideo,4)
                for k=1:size(rvideo,3)
                    frm=rvideo(:,:,k,ll);
                    % frames are transformed to columns, succensive rows represent successive frames
                    Stim(k+(ll-1)*size(rvideo,3),:)=frm(:);
                end
            end
            
            % loop pre-spike lags (along the depth of the filter)
            for nlag = 1:maxlag
                
                if nlag==1
                    % create output folder and cd in
                    fnam=[out_folder,'/neuron_',num2str(neuronum),'_medium_results'];
                    mkdir(fnam);
                    oldfold=cd(fnam);
                else
                end
                
                % reformat spike data into a single spike per frame vector aligned with stimulus movie
                sp=zeros(1,size(psth_by_frames,1)*size(psth_by_frames,3));
                hsp=zeros(1,size(psth_by_frames,3));
                for mm=1:size(psth_by_frames,1)
                    for ff=1:size(psth_by_frames,3)
                        sp(ff+(mm-1)*size(psth_by_frames,3))=psth_by_frames(mm,neuronum,ff);
                        hsp(ff)=hsp(ff)+psth_by_frames(mm,neuronum,ff);
                    end
                end
                
                totspikes=sum(hsp);
                message=['\nAnalyzing neuron ',num2str(neuronum),' session ',session_name,'...\n'];
                fprintf(message)
                message=['total number of spikes = ',num2str(totspikes),'\n'];
                fprintf(message)
                % message=[num2str(fwidth*fheight*50),' spikes required for a good STC reconstruction according to Rust rule\n'];
                
                if nlag==1
                    Dtotspikes=totspikes; %#ok<NASGU>
                    % plot and save global PSTHs
                    f1=figure;
                    set(f1,'Position',[10,10,1500,1000]);
                    plot(hsp,'color',[0.1,0.3,0.9],'LineWidth',1.5)
                    hold on
                    plot(sgolayfilt(hsp, 2, 27),'k','LineWidth',3)
                    legend('Mean Noise PSTH','Mean Noise PSTH - Smoothed');
                    ylabel('spikes per frame');
                    xlabel('frame number');
                    title(['Neuron ',num2str(neuronum),': ',num2str(totspikes),' spikes before'])
                    hold off
                    filename1='global PSTH.jpg';
                    set(gcf, 'PaperPositionMode', 'auto');
                    saveas(f1, filename1);
                    close all
                else
                end
                
                %% perform STA analysis at current lag ---------------------------
                
                % set analysisis parameters
                CriticalSize = 1e8;
                n=1;
                % reshape sp ad appy lag
                sp = sp';
                sp = circshift(sp,-nlag);
                sp(end-nlag+1:end)=zeros(1,nlag);
                % standardize stimulus
                Stim=(Stim-mean(Stim(:)))*(std(Stim(:))).^-1;
                
                if nlag==1 && neuronum==numstart % do the work on the correlation matrix only once
                    
                    % compute stimulus ensemble covariance matrix
                    [~,~,~,rawcov] = simpleSTC( Stim,sp, n, CriticalSize);
                    
                    % perform Tikhonov regularized inversion of the covariance matrix
                    covInv = inv(rawcov+ridgeparam*eye(size(rawcov)));
                    covInvsqrt = sqrtm(covInv); %#ok<NASGU>
                    
                end
                
                % compute STA
                [sta] = simpleSTA(Stim, sp, n,CriticalSize);
                
                % whiten STA
                wsta = covInv*sta; %#ok<MINV>
                
                % reshape and visualize sta
                stafr=reshape(sta,fheight,fwidth);
                if nlag==3
                    f1=figure;
                    set(f1,'Position',[10,10,1500,1000]); imagesc(stafr); colormap('gray'); set(gca,'dataAspectRatio',[1 1 1]); colorbar; title('raw STA filter');
                    set(gcf, 'PaperPositionMode', 'auto');
                    filename1='raw STA.jpg';
                    saveas(f1, filename1);
                    close all
                else
                end
                
                % reshape and visualize wsta
                wstafr=reshape(wsta,fheight,fwidth);
                if nlag==3
                    f1=figure;
                    set(f1,'Position',[10,10,1500,1000]); imagesc(wstafr); colormap('gray'); set(gca,'dataAspectRatio',[1 1 1]); colorbar; title('whitened STA filter');
                    set(gcf, 'PaperPositionMode', 'auto');
                    filename1='whitened STA.jpg';
                    saveas(f1, filename1);
                    close all
                else
                end
                
                %% perform permutation test at current lag ---------------------------
                
                % initialize permuation test variables
                nperm=30;
                staT=zeros(size(sta,1),nperm);
                wstaT=zeros(size(wsta,1),nperm);
                
                % loop over permutations
                for j=1:nperm
                    
                    % randomly reshuffle spike timestamps in time
                    spP=sp(randperm(length(sp)));
                    
                    % perform ST analysis again
                    spP = circshift(spP,-nlag);
                    spP(end-nlag+1:end)=zeros(1,nlag);
                    [staP] = simpleSTA( Stim, spP, n, CriticalSize);
                    
                    % whiten permuted STA
                    wstaP = covInv*staP; %#ok<MINV>
                    staT(:,j)=staP;
                    wstaT(:,j)=wstaP;
                    
                    message=['permutation number ',num2str(j),' completed\n'];
                    fprintf(message)
                end
                
                % compute permutation means
                musta=mean(staT,2);
                muwsta=mean(wstaT,2);
                
                % compute permuatation standard deviations
                sigsta=std(staT,0,2);
                sigwsta=std(wstaT,0,2);
                
                %% permutation test results visualization
                
                % reshape permutation results
                mustafr=reshape(musta,fheight,fwidth);
                muwstafr=reshape(muwsta,fheight,fwidth);
                
                % compute Z-scored STA, whitened STA and STC
                Zstafr = (stafr-mustafr)./mean(sigsta(:));
                Zwstafr = (wstafr-muwstafr)./mean(sigwsta(:));
                
                % visualize Z-scored STA, whitened STA and STC
                Z1=figure; set(Z1,'Position',[10,10,1500,1000]); set(gca,'dataAspectRatio',[1 1 1]);
                imagesc(Zstafr,[-5,5]); colormap('gray'); colorbar; title('Z-scored STA filter'); set(gcf, 'PaperPositionMode', 'auto');
                Z2=figure; set(Z2,'Position',[10,10,1500,1000]); set(gca,'dataAspectRatio',[1 1 1]);
                imagesc(Zwstafr,[-5,5]); colormap('gray'); colorbar; title('Z-scored whitened STA filter'); set(gcf, 'PaperPositionMode', 'auto');
                filename1=['Z_scored_STA_at_lag',num2str(nlag),'.jpeg'];
                filename2=['Z_scored_wSTA_at_lag',num2str(nlag),'.jpeg'];
                saveas(Z1, filename1);
                saveas(Z2, filename2);
                close all
                
                %% compute significance maps ------------------------------
                
                N=3; % confidence level for the map in sigmas
                
                % significativity for wSTA
                signifmapwsta=zeros(size(wsta));
                for i=1:size(wsta,1)
                    if abs(wsta(i)-muwsta(i))>N*sigwsta(i)
                        signifmapwsta(i)=abs(abs(wsta(i))-abs(muwsta(i)));
                    end
                end
                signifmapwstafr=reshape(signifmapwsta,fheight,fwidth);
                
                % significativity for STA
                signifmapsta=zeros(size(sta));
                for i=1:size(sta,1)
                    if abs(sta(i)-musta(i))>N*sigsta(i)
                        signifmapsta(i)=abs(abs(sta(i))-abs(musta(i)));
                    end
                end
                signifmapstafr=reshape(signifmapsta,fheight,fwidth);
                
                %% plot global analysis results
                
                imag=figure;
                set(imag,'Position',[10,10,1500,1000]);
                % raw results
                subplot(3,2,1)
                imagesc(wstafr,[min(wstafr(:)),max(wstafr(:))]); colormap('paruly'); set(gca,'dataAspectRatio',[1 1 1]);
                title('reconstructed wSTA filter', 'FontSize', 7);
                subplot(3,2,2)
                imagesc(stafr,[min(stafr(:)),max(stafr(:))]); colormap('paruly'); set(gca,'dataAspectRatio',[1 1 1]);
                title('reconstructed STA filter', 'FontSize', 7);
                % signifmaps
                subplot(3,2,3)
                imagesc(signifmapwstafr,[min(signifmapwstafr(:)),max(signifmapwstafr(:))+0.0001]); colormap('paruly'); set(gca,'dataAspectRatio',[1 1 1]);
                title('signifmap wSTA filter', 'FontSize', 7);
                subplot(3,2,4)
                imagesc(signifmapstafr,[min(signifmapstafr(:)),max(signifmapstafr(:))+0.0001]); colormap('paruly'); set(gca,'dataAspectRatio',[1 1 1]);
                title('signifmap STA filter', 'FontSize', 7);
                % raw permutation means
                subplot(3,2,5)
                imagesc(muwstafr,[min(muwstafr(:)),max(muwstafr(:))]); colormap('paruly'); set(gca,'dataAspectRatio',[1 1 1]);
                title('wSTA permutation mean', 'FontSize', 7);
                subplot(3,2,6)
                imagesc(mustafr,[min(mustafr(:)),max(mustafr(:))]); colormap('paruly'); set(gca,'dataAspectRatio',[1 1 1]);
                title('STA permutation mean', 'FontSize', 7);
                filenam=['ST analysis result at lag',num2str(nlag),'.jpeg'];
                saveas(imag, filenam);
                close all
                
                %% save results
                
                % store results at current lag
                Dstafr(:,:,nlag)=stafr;
                Dwstafr(:,:,nlag)=wstafr;
                DZwstafr(:,:,nlag)=Zwstafr;
                DZstafr(:,:,nlag)=Zstafr;
                Drawcov(:,:,nlag)=rawcov;
                
                if nlag==maxlag
                    % save results
                    save(['medium_ST_results_neuron_',num2str(neuronum),'.mat'],'Dstafr','Dwstafr','DZwstafr','DZstafr','Drawcov','Dneuronum','Dtotspikes')
                    % cd back to original folder
                    cd(oldfold)
                else
                end
                
                close all
                
            end
        end
    end
end
end

