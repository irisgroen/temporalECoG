function plot_chan_ind = checkForHDgrid(channels)

HDind = contains(channels.group, 'HDgrid');
plot_chan_ind = [];
if any(HDind)
    plot_chan_ind{1} = find(~HDind);
    HDind = find(HDind); HDindm = length(HDind)/2;
    plot_chan_ind{2} = HDind(1:HDindm-1); 
    plot_chan_ind{3} = HDind(HDindm:end);         
else 
    plot_chan_ind{1} = 1:height(channels); 
end

end