function [w_scores, Resid, Yhat] = compute_wscore_pcn(X_cov_no_const, Y_all, idx_HP, n_blocks)
% Compute w-scores using PCNtoolkit BLR (Python wrapper)
% Block-wise parallelization over Y columns
% n_blocks controls both number of blocks and number of parallel workers

arguments
    X_cov_no_const double
    Y_all          double
    idx_HP         logical
    n_blocks       (1,1) double {mustBeInteger, mustBePositive} = 3
end

N = size(Y_all,1);
if size(X_cov_no_const,1) ~= N
    error('X_cov_no_const and Y_all must have the same number of rows.');
end
if numel(idx_HP) ~= N
    error('idx_HP length must match the number of subjects.');
end

P = size(Y_all,2);

% If fewer variables than blocks, reduce blocks
n_blocks = min(n_blocks, P);

% Preallocate outputs
w_scores = nan(N, P);
Resid    = nan(N, P);
Yhat     = nan(N, P);

% Start/resize parallel pool if needed
if P > 1 && n_blocks > 1
    pool = gcp('nocreate');
    if isempty(pool) || pool.NumWorkers ~= n_blocks
        parpool('local', n_blocks);
    end
end

% Split columns into n_blocks approximately equal blocks
edges = round(linspace(1, P+1, n_blocks+1));
cols_blocks = cell(1, n_blocks);
for b = 1:n_blocks
    cols_blocks{b} = edges(b):(edges(b+1)-1);
end

% Per-worker initializer: import module + constant NumPy X/idx
pyConst = parallel.pool.Constant(@() local_init_pcn_wscore(X_cov_no_const, idx_HP));

% -------------------------
% Case 1: single variable
% -------------------------
if P == 1
    S = pyConst.Value;
    Y_np = py.numpy.array(Y_all);

    out = S.mod.pcn_wscore(S.X_np, Y_np, S.HP_np, ...
        pyargs('verbose', true, 'progress_every', int32(1), 'silence_internal', true));

    w_scores(:,1) = double(out{1});
    Resid(:,1)    = double(out{2});
    Yhat(:,1)     = double(out{3});
    return
end

% -------------------------
% Case 2: multiple variables (block-wise parallel)
% -------------------------
Z_blk  = cell(1, n_blocks);
R_blk  = cell(1, n_blocks);
Yh_blk = cell(1, n_blocks);

parfor b = 1:n_blocks
    cols = cols_blocks{b};
    if isempty(cols)
        Z_blk{b}  = [];
        R_blk{b}  = [];
        Yh_blk{b} = [];
        continue
    end

    S = pyConst.Value;

    Y_part = Y_all(:, cols);
    Y_np   = py.numpy.array(Y_part);

    out = S.mod.pcn_wscore(S.X_np, Y_np, S.HP_np, ...
        pyargs('verbose', true, 'progress_every', int32(1), 'silence_internal', true));

    Z_blk{b}  = double(out{1});
    R_blk{b}  = double(out{2});
    Yh_blk{b} = double(out{3});
end

% Merge blocks
for b = 1:n_blocks
    cols = cols_blocks{b};
    if isempty(cols), continue; end
    w_scores(:, cols) = Z_blk{b};
    Resid(:, cols)    = R_blk{b};
    Yhat(:, cols)     = Yh_blk{b};
end

end

% ---- helper: runs once per worker ----
function S = local_init_pcn_wscore(X_cov_no_const, idx_HP)
mod = py.importlib.import_module('pcn_wscore');
S.mod   = mod;
S.X_np  = py.numpy.array(X_cov_no_const);
S.HP_np = py.numpy.array(idx_HP);
end
