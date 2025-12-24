function [net_max, normv, nonnormv] = calc_normativity(Dice_coeff)
% Compute normativity metrics from a Dice coefficient matrix.
%
% Input
%   Dice_coeff : [nAtoms x nTargets] double
%                nTargets can be networks/components (cortical) or subcortical groups.
%
% Output
%   net_max    : [1 x nTargets] network-wise (column-wise) maximum Dice
%   normv      : scalar, mean(net_max)
%   nonnormv   : scalar, 1 - mean(atom-wise maximum Dice)

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
