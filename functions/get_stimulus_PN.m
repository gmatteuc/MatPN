function stimul = get_stimulus_PN( SF, TF, DIR,stimulustype )

pars=set_pars_PN();
gcode=get_indexes_PN( SF, TF, DIR, stimulustype );
gpath=pars.stim_folder;

switch stimulustype
    case 'grating'
        
        gname=['0',num2str(gcode),'_grating_',num2str(DIR),'_',num2str(SF),'_',num2str(TF)];
        oldfold=cd([gpath,filesep,gname]);
        % stimulus + isi
        stimul=0.5*ones(1080/3,1920/3,25+27+25);
        for l=1:27
            if l>9
                stimul(:,:,25+l) = imresize(im2double(imread(['grating_',num2str(DIR),'_frame_0',num2str(l)],'bmp')),1/3);
            else
                stimul(:,:,25+l) = imresize(im2double(imread(['grating_',num2str(DIR),'_frame_00',num2str(l)],'bmp')),1/3);
            end
        end
        stimul=stimul-0.5;
        
    case 'plaid'
        
        gname=['0',num2str(num2str(gcode)),'_plaid_',num2str(DIR),'_',num2str(SF),'_',num2str(TF)];
        oldfold=cd([gpath,filesep,gname]);
        % stimulus + isi
        stimul=0.5*ones(1080/3,1920/3,25+27+25);
        for l=1:27
            if l>9
                stimul(:,:,25+l) = imresize(im2double(imread(['plaid_',num2str(DIR),'_frame_0',num2str(l)],'bmp')),1/3);
            else
                stimul(:,:,25+l) = imresize(im2double(imread(['plaid_',num2str(DIR),'_frame_00',num2str(l)],'bmp')),1/3);
            end
        end
        stimul=stimul-0.5;
        
end

cd(oldfold)

end

