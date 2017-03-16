%%
% Function to extract VE timecourse from sensor level data using the
% spatial filter computed from mous_lcmv_parcellate
%
% Inputs:
% - data = preprocessed sensor level data with all trials
% - parcel_name = string of parcel required. MUST match name within the
% atlas
% - p = parcellated data computed using mous_lcmv_parcellate
%
%% 

function [virtualchannel_raw] = get_VE_from_parcellated_filter(data,parcel_name,p)

indx = find(ismember(p.label,parcel_name)); % find index of the required label
    
    virtualchannel_raw = [];
    label = strrep(parcel_name,'_',' '); %Remove underscores - better for plotting later
    virtualchannel_raw.label = {label};
    virtualchannel_raw.trialinfo = data.trialinfo;
    for i=1:(length(data.trialinfo))
        % note that this is the non-filtered raw data
        virtualchannel_raw.time{i}       = data.time{i};
        virtualchannel_raw.trial{i}(1,:) = p.filter{indx,1}(1,:).*p.filter{indx,1}(2,:)*data.trial{i}(:,:);
    end
end
    
    