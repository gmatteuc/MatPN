function channel = GetChannel(nn,SpikesChannels)

    %SpikesChannels = SPIKES.channel;
    %channels = find(cellfun(@(x) any(x==nn),SpikesChannels));
    try
    channel=SpikesChannels{nn};
    catch
    channel=SpikesChannels(nn);
    end
    
    if isempty(channel)
        display('Warning : no neurons in these channels !');        
    end
     

end