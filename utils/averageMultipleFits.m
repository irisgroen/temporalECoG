function [m,se,channels] = averageMultipleFits(data, channels, fun)

if ~exist('fun', 'var') || isempty(fun)
    fun = @mean;
end

% if there are more than one dimension per channel, vectorize data
if ndims(data) == 3
    multiDimData = true;
    dataSz = size(data);
    data = reshape(data, [dataSz(1) * dataSz(2) dataSz(3)]);
else
    multiDimData = false;
end

% reshape data such that repeated fits are in last dimension
nAreas = length(unique(channels.index));
dataSz2 = size(data);
data = reshape(data, [dataSz2(1) nAreas dataSz2(2)/nAreas]);

% average across fits
m = fun(data,3,'omitnan');

% compute 68% confidence interval
se = squeeze(prctile(data,[15.87 84.13],3));

channels = channels(1:nAreas,:);

% reshape back to original shape
if multiDimData
    m = reshape(m, [dataSz(1) dataSz(2) nAreas]); 
    if ~isempty(se)
        se = reshape(se, [dataSz(1) dataSz(2) 2 nAreas]);
    end
end

end

