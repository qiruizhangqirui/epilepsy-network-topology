function result = read_network_correspondence(filename)
%READ_NETWORK_CORRESPONDENCE  Read network correspondence CSV and extract Dice coefficients and p-values
%
% Input:
%   filename - Path to a CSV file in which each cell is formatted as:
%              "value (p=0.1234)"
%
% Output:
%   result - Structure with fields:
%       .Networks    - Column names (networks)
%       .Atoms       - Row names (atoms / brain regions)
%       .Dice_coeff  - Numeric matrix [nAtoms x nNetworks]
%       .P_value     - Numeric matrix [nAtoms x nNetworks]

    % Read table with all entries forced to strings
    opts = detectImportOptions(filename, 'NumHeaderLines', 0);
    opts = setvartype(opts, 'char');
    data = readtable(filename, opts);

    % Extract network and atom labels
    result.Networks = data.Properties.VariableNames(2:end);
    result.Atoms    = data{:,1};

    % Initialize output matrices
    nRows = height(data);
    nCols = width(data) - 1;
    result.Dice_coeff = NaN(nRows, nCols);
    result.P_value    = NaN(nRows, nCols);

    % Regular expression to extract value and p-value
    % Matches strings like: "0.45 (p=0.012)"
    pattern = '([\d\.Ee+-]+)\s*\(p=([\d\.Ee+-]+)\)';

    % Parse each cell
    for i = 1:nRows
        for j = 1:nCols
            cellStr = strtrim(data{i, j+1}{1});
            tokens = regexp(cellStr, pattern, 'tokens');
            if ~isempty(tokens)
                result.Dice_coeff(i, j) = str2double(tokens{1}{1});
                result.P_value(i, j)    = str2double(tokens{1}{2});
            end
        end
    end
end
