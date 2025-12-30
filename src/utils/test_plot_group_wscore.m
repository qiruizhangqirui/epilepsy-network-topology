% Test script for plot_group_wscore
clear; clc;

% Create dummy data
rng(123);
N = 100;
w_scores = randn(N, 1);
group_names = {'GroupA', 'GroupB', 'GroupC'};
GroupTable = table(rand(N, 1) > 0.5, rand(N, 1) > 0.5, rand(N, 1) > 0.5, ...
    'VariableNames', group_names);

% Call the function with 4 outputs
try
    [w_vec, p_vec, labels, std_vec] = plot_group_wscore(w_scores, GroupTable, 'Test Plot');

    fprintf('Function call successful.\n');
    fprintf('Median W-scores: %s\n', mat2str(w_vec, 3));
    fprintf('Standard Deviations: %s\n', mat2str(std_vec, 3));

    if length(std_vec) == 3
        fprintf('Verification PASSED: std_vec has correct length.\n');
    else
        fprintf('Verification FAILED: std_vec has incorrect length.\n');
    end

    disp('Closing figure...');
    close(gcf);

catch ME
    fprintf('Verification FAILED with error: %s\n', ME.message);
end
