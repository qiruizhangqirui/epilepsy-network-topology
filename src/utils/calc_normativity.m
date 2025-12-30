function [net_max, normv, nonnormv] = calc_normativity(Dice_coeff)
% =========================================================================
% Function: calc_normativity
% -------------------------------------------------------------------------
% Description:
%   Compute normativity metrics from a Dice coefficient matrix.
%
% Inputs:
%   Dice_coeff - [nAtoms x nTargets] double matrix. 
%                nTargets can be networks/components (cortical) or subcortical groups.
%
% Outputs:
%   net_max    - [1 x nTargets] Network-wise (column-wise) maximum Dice
%   normv      - Scalar, mean of net_max (Network Normativity)
%   nonnormv   - Scalar, 1 - mean of atom-wise maximum Dice (Non-normativity)
%
% Author: Qirui Zhang, Farber Institute for Neuroscience, Thomas Jefferson University
% Date: 12/30/2025
% =========================================================================

% Handle empty input robustly
if isempty(Dice_coeff)
    net_max  = [];
    normv    = nan;
    nonnormv = nan;
    return
end

% Network-wise maximum Dice (1 x nTargets)
net_max = max(Dice_coeff, [], 1);

% Atom-wise maximum Dice (nAtoms x 1)
atom_max = max(Dice_coeff, [], 2);

% Metrics
normv    = mean(net_max,  'omitnan');
nonnormv = 1 - mean(atom_max, 'omitnan');

end
