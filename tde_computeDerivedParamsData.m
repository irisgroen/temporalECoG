function [derivedPrm, derivedPrmNames] = tde_computeDerivedParamsData(data, t)
% Assumes data = samples x stim x channels
%
% 2022 Iris Groen

[nSamples, nStim, nChan] = size(data);
data = smooth(data(:),5);
data = reshape(data, [nSamples, nStim, nChan]);
derivedPrm = nan(3,nChan);
cond_ind = [6 6 1];

% Find the timepoint corresponding to the maximum magnitude
[M,I] = max(data,[],1);
I = squeeze(I);
T = repmat(t, [1 nStim nChan]);

t_max = T(I);
derivedPrm(1,:) = t_max(cond_ind(1),:); 

% Find the ratio between the asymptote and the maximum magnitude
t_offset = t == 0.5;
A = data(t_offset,:,:);
r_asymp = squeeze(A./M);
derivedPrm(2,:,:) = r_asymp(cond_ind(2),:);

% Find the full width at half max value
halfMax = (min(data) + max(data)) / 2;
data_thresh = data >= halfMax;
fwhmx = nan(1,nChan);
for jj = 1:nChan
    index1 = find(data_thresh(:,cond_ind(3),jj), 1, 'first');
    index2 = find(data_thresh(:,cond_ind(3),jj), 1, 'last');
    fwhmx(jj) = t(index2) - t(index1);
end

derivedPrm(3,:)    = fwhmx;

derivedPrmNames =  {'Time to peak', 'Ratio sustained/transient', 'Fwhm'};
end