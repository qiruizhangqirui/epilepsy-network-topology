% =========================================================================
% Script: S1_Org_outputs.m
% -------------------------------------------------------------------------
% Description:
%   Organizes input data (NCT derivatives, Hubness maps) into a unified
%   structure and output directory.
%
%   - Reads subject list from Excel.
%   - Parses noise atom list.
%   - Aggregates correspondence results from CSVs.
%   - Computes subcortical dice coefficients.
%   - Computes clean/raw k-hubness maps and saves as .nii.gz.
%
% Inputs:
%   - data/Subjects.xlsx
%   - data/NCT_Derivatives/
%   - data/SPARK_Derivatives/
%   - data/Noise_atoms_new.txt
%   - data/NCT_atlases/
%
% Outputs:
%   - outputs/All_derivatives_struct.mat
%   - outputs/Hubness/*.nii.gz (compressed)
%
% Author: Qirui Zhang,  Farber Institute for Neuroscience, Department of Neurology, Thomas Jefferson University
% Email: qirui.zhang@jefferson.edu;fmrizhangqr@126.com
% Date: 12/24/2025
% =========================================================================

clear; clc;

%% =========================
% Part 1: Path Configuration
% =========================
% Determine project root (assumes src/pipeline/S1_Org_outputs.m)
this_file = mfilename('fullpath');

if ~isempty(this_file)
    % Case 1: script/function is being executed from a .m file
    project_root = fileparts(fileparts(fileparts(this_file)));
else
    % Case 2: executed from Command Window or interactive context
    project_root = fileparts(fileparts(pwd));
end

data_dir     = fullfile(project_root, 'data');
result_dir   = fullfile(project_root, 'outputs');

if ~exist(result_dir, 'dir'); mkdir(result_dir); end

% Add utils to path
addpath(genpath(fullfile(project_root, 'src', 'utils')));

pathSubjects = fullfile(data_dir, 'Subjects.xlsx');

% Anonymized outputs created by your previous scripts
nct_root     = fullfile(data_dir, 'NCT_Derivatives');        % <newID>/*.csv
atom_root    = fullfile(data_dir, 'SPARK_Derivatives');      % <newID>/ratom*.nii.gz

% Anonymized noise list
pathNoise    = fullfile(data_dir, 'Noise_atoms_new.txt');

%% =========================
% Part 2: Load Subject Metadata
% =========================
subjectTable = readtable(pathSubjects);
newIDs       = string(subjectTable.ID_new);
subjectIDs   = cellstr(newIDs);   % use anonymized IDs everywhere below

%% =========================
% Part 3: Parse Noise List
% =========================
noiseLines  = readlines(pathNoise);
noiseLines  = strtrim(noiseLines);
noiseLines  = erase(noiseLines, '.png');   % remove suffix

noiseAtomsStruct = struct();

for i = 1:length(noiseLines)
    full_name = noiseLines(i);

    % Split into: subject_id + image_name
    tokens = regexp(full_name, '^([^_]+)_(.+)$', 'tokens');

    if ~isempty(tokens)
        sID = tokens{1}{1};
        imgName = tokens{1}{2};

        if isfield(noiseAtomsStruct, sID)
            noiseAtomsStruct.(sID){end+1} = imgName;
        else
            noiseAtomsStruct.(sID) = {imgName};
        end
    end
end


%% =========================
% Part 4: Aggregate Data & Compute Dice
% =========================
All_derivatives_struct = struct();

for iSub = 1:length(subjectIDs)

    subjID = subjectIDs{iSub};

    % ---- Load NCT output CSVs (anonymized folder)
    NCT_dir  = fullfile(nct_root, subjID);
    NCT_list = dir(fullfile(NCT_dir, '*.csv'));

    for iFile = 1:length(NCT_list)
        atlas_name = extractBetween(NCT_list(iFile).name, 'correspondence_', '.csv');
        NCT_result = read_network_correspondence(fullfile(NCT_dir, NCT_list(iFile).name));
        All_derivatives_struct.(subjID).(atlas_name{1}) = NCT_result;
    end

    % Use the last loaded atlas to get atom names (assumes all atlases share the same Atoms list)
    atom_list = All_derivatives_struct.(subjID).(atlas_name{1}).Atoms;

    % ---- Map noisy atom filenames to indices in atom_list
    if ~isfield(noiseAtomsStruct, subjID)
        noise_idx = NaN;
    else
        Noise_Atom = noiseAtomsStruct.(subjID);
        noise_idx  = zeros(1, length(Noise_Atom));

        for i = 1:length(Noise_Atom)
            idx = find(strcmp(atom_list, Noise_Atom{i}));
            if ~isempty(idx)
                noise_idx(i) = idx;
            else
                noise_idx(i) = NaN;
            end
        end
    end
    All_derivatives_struct.(subjID).noise_idx = noise_idx;

    % ---- Subcortical dice (atoms vs subcortical atlas)
    atom_nii_dir  = dir(fullfile(atom_root, subjID, 'ratom*.nii*'));
    atom_nii_list = fullfile({atom_nii_dir.folder}, {atom_nii_dir.name});

    atomNiftiList  = atom_nii_list;
    subcorticalRef = fullfile(data_dir, 'Subcortical_8ROI_FLS2mm.nii');
    mask           = fullfile(data_dir, 'FSLMNI2mm.nii');
    % Compute Dice
    diceMatrix = image_cross_dice(atomNiftiList, subcorticalRef, mask, 1);

    All_derivatives_struct.(subjID).Subcortical.Atoms      = All_derivatives_struct.(subjID).AS200Y17.Atoms;
    All_derivatives_struct.(subjID).Subcortical.Networks   = {'Accumbens','Amygdala','Caudate','Hippocampus','Pallidum','Putamen','Thalamus','Cerebellum'};
    All_derivatives_struct.(subjID).Subcortical.Dice_coeff = diceMatrix;

end

save(fullfile(result_dir, 'All_derivatives_struct.mat'), 'All_derivatives_struct');

%% =========================
% Part 5: Hubness Map Generation
% =========================
indir = [atom_root filesep];
subj_list = dir(indir);

threshold = 0;

for i = 3:length(subj_list)

    subjID = subj_list(i).name;

    nii_list = dir(fullfile(indir, subjID, 'ratom*.nii*'));

    nii_paths = arrayfun(@(f) fullfile(f.folder, f.name), ...
        nii_list, 'UniformOutput', false);

    % -------------------------
    % Raw k-hubness
    % -------------------------
    info0 = niftiinfo(nii_paths{1});
    img0  = niftiread(info0);

    img_sum = double(img0 > threshold);
    for k = 2:length(nii_paths)
        img_k = niftiread(nii_paths{k});
        img_sum = img_sum + double(img_k > threshold);
    end

    img_sum = cast(img_sum, 'like', img0);
    info0.Datatype = class(img_sum);

    output_nii = fullfile(indir, subjID, [subjID '_khubness.nii']);
    niftiwrite(img_sum, output_nii, info0);

    % ---- compress to .nii.gz
    gzip(output_nii);
    delete(output_nii);

    % -------------------------
    % Cleaned k-hubness
    % -------------------------
    if ~isfield(noiseAtomsStruct, subjID)

        % No noise to remove -> same as raw
        output_clean = fullfile(indir, subjID, [subjID '_khubness_clean.nii']);
        niftiwrite(img_sum, output_clean, info0);
        gzip(output_clean);
        delete(output_clean);

    else
        Noise_Atom = noiseAtomsStruct.(subjID);
        remove_idx = false(length(nii_paths), 1);

        for j = 1:length(Noise_Atom)
            remove_idx = remove_idx | contains(nii_paths, Noise_Atom{j});
        end

        nii_paths_clean = nii_paths(~remove_idx);

        if ~isempty(nii_paths_clean)

            info0c = niftiinfo(nii_paths_clean{1});
            img0c  = niftiread(info0c);

            img_sum = double(img0c > threshold);
            for k = 2:length(nii_paths_clean)
                img_k = niftiread(nii_paths_clean{k});
                img_sum = img_sum + double(img_k > threshold);
            end

            img_sum = cast(img_sum, 'like', img0c);
            info0c.Datatype = class(img_sum);

            output_nii = fullfile(indir, subjID, [subjID '_khubness_clean.nii']);
            niftiwrite(img_sum, output_nii, info0c);

            gzip(output_nii);
            delete(output_nii);
        end
    end
end
