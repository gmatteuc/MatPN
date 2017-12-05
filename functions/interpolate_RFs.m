function [ RFs ] = interpolate_RFs( RFs )

%[ RFs ] = interpolate_RFs( RFs )
%
%Add to RFs datastructure the interpolated RFs

%--------------------------------------------------------------------------

RFs.DZwstafr_int=zeros((size(RFs.DZwstafr,1))*30,(size(RFs.DZwstafr,2))*30,size(RFs.DZwstafr,3),size(RFs.DZwstafr,4));
for nn=1:max(RFs.Dneuronum);
    
    % spline interpolation
    for jj=1:size(RFs.DZwstafr,3)
        % find index corresponding to current neuron number (nn)
        kk=find(RFs.Dneuronum==nn);
        [Xo,Yo] = meshgrid(1:1:size(RFs.DZwstafr,2), 1:1:size(RFs.DZwstafr,1));
        [Xp,Yp] = meshgrid(linspace(1,size(RFs.DZwstafr,2),(size(RFs.DZwstafr,2))*30), linspace(1,size(RFs.DZwstafr,1),(size(RFs.DZwstafr,1))*30));
        if isempty(kk)
        else
            fram=squeeze(RFs.DZwstafr(:,:,jj,kk));
            if sum(isnan(fram(:)))
            fram=zeros(size(RFs.DZwstafr,1),size(RFs.DZwstafr,2));
            else
            end
            RFs.DZwstafr_int(:,:,jj,kk)=interp2(Xo,Yo ,fram,Xp,Yp,'spline');
        end
    end
end

end

