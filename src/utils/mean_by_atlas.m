function [Network_value] = mean_by_atlas(Volume, Atlas, atlas_order)
% =========================================================================
% Function: mean_by_atlas
% -------------------------------------------------------------------------
% Description:
%   Calculates mean intensity values within atlas-defined regions.
%
% Inputs:
%   Volume      - Path to Nifti volume file (data)
%   Atlas       - Path to Nifti atlas file (labels)
%   atlas_order - Vector of label IDs to extract
%
% Output:
%   Network_value - [1 x N_labels] Mean value for each requested label
%
% Author: Qirui Zhang, Farber Institute for Neuroscience, Thomas Jefferson University
% Date: 12/30/2025
% =========================================================================

info = niftiinfo(Volume);
V = niftiread(info);
V_2d = reshape(V, [size(V,1)*size(V,2)*size(V,3), 1]);

info = niftiinfo(Atlas);
A = niftiread(info);
A_2d = reshape(A, [size(V,1)*size(V,2)*size(V,3), 1]);

% network_number = sort(unique(A_2d));
%
% for i = 2:length(network_number)
%
% index = find(A_2d == network_number(i));
% Network_value(:,i-1) = mean(V_2d(index));
% end

network_number = length(atlas_order);
Network_value = zeros(1, network_number);

for i = 1:network_number
  index = (A_2d == atlas_order(i));
  value = V_2d(index);
  % value1 = value(find(value ~=0 ));
  Network_value(1, i) = nanmean(value);
end
end