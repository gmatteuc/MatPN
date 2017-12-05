function [ centres_simple, centres_complex, bool_aic ] = GLMcenters( MI_distr,replim,choosenmodel ) 

% ------------------------------------------------------------------%
% perform gaussian mixture model analysis and return gaussian centres

X=MI_distr;
if choosenmodel==1
    k=1;
    repcount=0;
    bool_conv=0;
    
    while bool_conv==0 && repcount<replim
        % fit mixture model
        gm{k} = fitgmdist(X,k);
        repcount=repcount + 1;
        if gm{k}.Converged == 1
            bool_conv=1;
        elseif gm{k}.Converged == 0 && repcount==replim
            fprintf(['\nWARNING!!! convergence not reached after all trials at k=',num2str(k),' complex_center =',num2str(min(gm{k}.mu),'%0.2f'),' simple_center =',num2str(max(gm{k}.mu),'%0.2f'),'\n'])
        else
            %             fprintf('WARNING! convergence not reached')
        end
        
    end
else
    for k = 1:2
        
        repcount=0;
        bool_conv=0;
        
        while bool_conv==0 && repcount<replim
            % fit mixture model
            gm{k} = fitgmdist(X,k);
            repcount=repcount + 1;
            if gm{k}.Converged == 1
                bool_conv=1;
            elseif gm{k}.Converged == 0 && repcount==replim
                fprintf(['\nWARNING!!! convergence not reached after all trials at k=',num2str(k),' complex_center =',num2str(min(gm{k}.mu),'%0.2f'),' simple_center =',num2str(max(gm{k}.mu),'%0.2f'),'\n'])
            else
                %             fprintf('WARNING! convergence not reached')
            end
            
        end
    end
end
centres_complex=min(gm{choosenmodel}.mu);
centres_simple=max(gm{choosenmodel}.mu);
if choosenmodel==2
    bool_aic=gm{1}.AIC>gm{2}.AIC;
elseif choosenmodel==1
    bool_aic=NaN;
end
% -----------------------------------------------%

end

