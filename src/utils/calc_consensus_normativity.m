function [ConsensusDice_S200, ConsensusDice_17, TemplateNamesOut, nParcelPerNet17] = ...
    calc_consensus_normativity(NetMaxCell, TemplateNamesIn, cons)

Nsub  = size(NetMaxCell, 1);
NPARC = cons.NPARC;

SumAll   = zeros(Nsub, NPARC);
CountVec = zeros(1, NPARC);

TemplateNamesOut = {};

for jj = 1:numel(TemplateNamesIn)
    atlas_name = TemplateNamesIn{jj};

    % Data_mat: [Nsub x Mnet]
    Data_mat = cell2mat(NetMaxCell(:, jj));
    Mnet = size(Data_mat, 2);

    if ~ismember(atlas_name, {'HCPICA','UKBICA'})
        % ----- Label atlas -----
        [lh_labels, rh_labels] = get_conte69_labels(atlas_name, cons.P);
        L = [double(lh_labels(:)); double(rh_labels(:))];

        ok = (cons.labels_all > 0) & (L > 0);
        counts = accumarray([L(ok), cons.labels_all(ok)], 1, [Mnet, NPARC]);

        Tproj = bsxfun(@rdivide, counts, max(cons.cnt_parcel, 1));
        Y_tmp = Data_mat * Tproj;

        SumAll   = SumAll + Y_tmp;
        CountVec = CountVec + 1;
        TemplateNamesOut{end+1} = atlas_name; %#ok<AGROW>

    else
        % ----- ICA atlas -----
        if strcmp(atlas_name, 'HCPICA')
            icadir = cons.HCPICA_dir;
        else
            icadir = cons.UKBICA_dir;
        end

        d = dir(fullfile(icadir, '*.mat'));
        if isempty(d)
            warning('Skipped %s: no ICA component .mat files found in %s.', atlas_name, icadir);
            continue
        end

        fnames = {d.name};
        nums = cellfun(@(s) str2double(regexp(s,'(\d+)(?=\.mat$)','match','once')), fnames);
        [~, idx_sort] = sort(nums, 'ascend', 'MissingPlacement','last');
        comp_files = fnames(idx_sort);
        comp_names = regexprep(comp_files, '\.mat$', '');

        nComp = min(Mnet, numel(comp_names));
        if nComp < 1
            warning('Skipped %s: component count mismatch.', atlas_name);
            continue
        end

        M = false(nComp, NPARC);
        for cc = 1:nComp
            comp_name = comp_names{cc};
            cons.P.(comp_name) = fullfile(icadir, [comp_name '.mat']); %#ok<STRNU>

            [lh_mask, rh_mask] = get_conte69_labels(comp_name, cons.P);
            mask = [double(lh_mask(:)) > 0; double(rh_mask(:)) > 0];

            ok = (cons.labels_all > 0) & mask;
            if any(ok)
                parcels_hit = unique(cons.labels_all(ok));
                parcels_hit(parcels_hit <= 0) = [];
                M(cc, parcels_hit) = true;
            end
        end

        Y_tmp = zeros(Nsub, NPARC);
        for p = 1:NPARC
            idx_comp = M(:, p);
            if any(idx_comp)
                vals = Data_mat(:, idx_comp);
                Y_tmp(:, p) = mean(vals, 2, 'omitnan');
            else
                Y_tmp(:, p) = 0;
            end
        end

        cover_vec = any(M, 1);

        SumAll   = SumAll + Y_tmp;
        CountVec = CountVec + double(cover_vec);
        TemplateNamesOut{end+1} = atlas_name; %#ok<AGROW>
    end
end

denom = max(CountVec, 1);
ConsensusDice_S200 = bsxfun(@rdivide, SumAll, denom);

% 200 -> 17
ConsensusDice_17 = nan(Nsub, 17);
nParcelPerNet17  = zeros(17, 1);

for k = 1:17
    colk = (cons.map200to17 == k);
    nParcelPerNet17(k) = sum(colk);
    if any(colk)
        ConsensusDice_17(:, k) = mean(ConsensusDice_S200(:, colk), 2, 'omitnan');
    end
end

end


function [lh_labels, rh_labels] = get_conte69_labels(Altas_name, P)
%GET_CONTE69_LABELS  Load vertex-wise labels in Conte69 (fs_LR_32k) space.
%
% This function returns left- and right-hemisphere vertex labels for a given
% cortical atlas or ICA component in Conte69 space.
%
% Behavior by atlas type:
%   - AS200Y17:
%       (1) Load the Schaefer-200 parcellation (lh_labels / rh_labels)
%       (2) Map parcel labels (1–200) to Yeo 17-network labels using the
%           provided assignment vector (labels of 0 remain 0)
%   - EG17 / MG360J12 / TY7:
%       Directly return the lh_labels / rh_labels stored in the atlas .mat file
%   - ICA components or other inputs:
%       The .mat file is assumed to contain binary or labeled vertex masks
%       (lh_labels / rh_labels) for that specific component
%
% Inputs
%   Altas_name : char or string
%       Atlas name (e.g., 'AS200Y17', 'EG17', 'MG360J12', 'TY7', or ICA component name)
%   P : struct
%       Structure containing full paths to atlas/component .mat files
%
% Outputs
%   lh_labels, rh_labels : double column vectors
%       Vertex-wise labels for left and right hemispheres in Conte69 space

    switch Altas_name

        case 'AS200Y17'
            % Load Schaefer-200 parcellation (vertex-wise ROI labels)
            S200 = load(P.AS200Y17_mat);
            lh_labels = double(S200.lh_labels(:));
            rh_labels = double(S200.rh_labels(:));

            % Load ROI-to-network mapping (200 parcels -> 17 networks)
            Smap = load(P.AS200Y17_assign_mat);
            mapping = double(Smap.mapping(:));   % length = 200

            % Apply mapping only to valid (>0) parcel labels
            idxL = lh_labels > 0;
            idxR = rh_labels > 0;
            lh_labels(idxL) = mapping(lh_labels(idxL));
            rh_labels(idxR) = mapping(rh_labels(idxR));

        case 'EG17'
            S = load(P.EG17_mat);
            lh_labels = double(S.lh_labels(:));
            rh_labels = double(S.rh_labels(:));

        case 'MG360J12'
            S = load(P.MG360J12_mat);
            lh_labels = double(S.lh_labels(:));
            rh_labels = double(S.rh_labels(:));

        case 'TY7'
            S = load(P.TY7_mat);
            lh_labels = double(S.lh_labels(:));
            rh_labels = double(S.rh_labels(:));

        otherwise
            % ICA component or other atlas:
            % The corresponding .mat file is assumed to store vertex-wise
            % masks or labels in lh_labels / rh_labels
            S = load(P.(Altas_name));
            lh_labels = double(S.lh_labels(:));
            rh_labels = double(S.rh_labels(:));
    end
end
