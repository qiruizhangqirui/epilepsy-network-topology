function parcel_to_volume(parcel_value, atlas, output_name, atlas_order)
% =========================================================================
% Function: parcel_to_volume
% -------------------------------------------------------------------------
% Description:
%   Maps parcel-wise values back to a Nifti volume based on an atlas.
%
% Inputs:
%   parcel_value - Vector of values to map
%   atlas        - Path to Nifti atlas file
%   output_name  - Path for output Nifti file
%   atlas_order  - Vector of label IDs corresponding to parcel_value entries
%
% Author: Qirui Zhang, Farber Institute for Neuroscience, Thomas Jefferson University
% Date: 12/30/2025
% =========================================================================

info = niftiinfo(atlas);
V = niftiread(info);
V_2d = reshape(V, [size(V,1)*size(V,2)*size(V,3), 1]);
V1_2d = zeros(size(V_2d));

% Map values
for i = 1:length(atlas_order)
    ind = (V_2d == atlas_order(i));
    V1_2d(ind) = parcel_value(i);
end

V1_3d = reshape(V1_2d, [size(V,1), size(V,2), size(V,3)]);

if contains(info.Datatype, 'int16')
    V1_3d = int16(V1_3d);
elseif contains(info.Datatype, 'single')
    V1_3d = single(V1_3d);
elseif contains(info.Datatype, 'uint8')
    info.Datatype = 'double';
end

niftiwrite(V1_3d, output_name, info);
end