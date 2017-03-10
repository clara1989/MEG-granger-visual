%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to apply spatial filters from mous_lcmv_parcellate to 
% sensor-level fourier output from ft_freqanalysis 
%
% Inputs:
% - freq : sensor-level freq data from ft_freqanalysis (output = 'fourier')
% - parcel_name = string of parcel required. MUST match name within the
% atlas
% - p = parcellated data computed using mous_lcmv_parcellate
%
% Output:
% - fourier_source = nrpttap*chan complex double (fourier)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [fourier_source] = extract_fourier_parcellation(freq,parcel_name,p)

nchan = numel(freq.label);
nfreq  = numel(freq.freq);
nrpttap = size(freq.fourierspctrm,1);

indx = find(ismember(p.label,parcel_name));

fourier_source = p.filter{indx,1}(1,:).*p.filter{indx,1}(2,:)*reshape(permute(freq.fourierspctrm,[2 3 1]), [nchan nfreq*nrpttap]);
fourier_source = reshape(fourier_source, [nfreq nrpttap]);
fourier_source = rot90(fourier_source);
end
