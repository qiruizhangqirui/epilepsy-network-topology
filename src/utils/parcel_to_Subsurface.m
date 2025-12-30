function parcel_to_Subsurface(parcel_value, output_name, color_range, cmap)
% =========================================================================
% Function: parcel_to_Subsurface
% -------------------------------------------------------------------------
% Description:
%   Visualizes subcortical values and saves the figure as a TIFF image.
%   Uses 'plot_subcortical' (must be in path).
%
% Inputs:
%   parcel_value - Data vector for subcortical regions
%   output_name  - File path for the output image
%   color_range  - [min max] for colormap scaling
%   cmap         - Colormap to use
%
% Author: Qirui Zhang, Farber Institute for Neuroscience, Thomas Jefferson University
% Date: 12/30/2025
% =========================================================================

figure;
plot_subcortical(parcel_value, 'color_range', color_range, 'cmap', cmap, 'ventricles', 'False');

print(gcf, '-dtiff', '-r300', output_name);
close(gcf);

end