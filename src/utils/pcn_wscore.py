# -*- coding: utf-8 -*-
# =========================================================================
# Script: pcn_wscore.py
# -------------------------------------------------------------------------
# Description:
#   Utility functions for computing w-scores using PCNtoolkit.
#   Includes:
#     - pcn_wscore: Main function to fit normative model and compute deviation (w-score).
#     - Helper functions for data standardization and I/O suppression.
#
# Author: Qirui Zhang, Farber Institute for Neuroscience, Thomas Jefferson University
# Date: 12/30/2025
# =========================================================================

import numpy as np
import warnings
import logging
import io, contextlib

from pcntoolkit import NormativeModel, NormData
from pcntoolkit.regression_model.blr import BLR
from pcntoolkit.math_functions.basis_function import (
    BsplineBasisFunction,
    PolynomialBasisFunction,
)

EPS = 1e-12



def _to_2d(a):
    """Convert input to a 2D NumPy array."""
    a = np.asarray(a, dtype=float)
    a = np.squeeze(a)
    if a.ndim == 1:
        a = a[:, None]
    if a.ndim != 2:
        raise ValueError(f"Expected a 2D array, got shape {a.shape}")
    return a


def _standardize_hp(X, idx_hp):
    """Standardize covariates using statistics from healthy participants only."""
    X = np.asarray(X, dtype=float)
    hp = np.asarray(idx_hp, dtype=bool).reshape(-1)

    mu = X[hp].mean(axis=0)
    std = X[hp].std(axis=0, ddof=1)
    zero = std < EPS

    Xs = X - mu
    if np.any(~zero):
        Xs[:, ~zero] = Xs[:, ~zero] / std[~zero]
    if np.any(zero):
        Xs[:, zero] = 0.0

    return Xs, mu, std, zero


def _get_var(ds, keys):
    """Extract a variable from a pcntoolkit prediction object using multiple aliases."""
    try:
        var_names = list(ds.data_vars.keys())
    except Exception as e:
        raise TypeError(
            f"Prediction object does not contain data_vars: {type(ds)}; original error: {e}"
        )

    lower_map = {k.lower(): k for k in var_names}
    for k in keys:
        lk = k.lower()
        if lk in lower_map:
            return np.asarray(ds[lower_map[lk]].values, dtype=float)

    raise KeyError(f"Variables {keys} not found. Available keys: {var_names}")


@contextlib.contextmanager
def _silence_stdout_stderr(silence=True):
    """Optionally suppress stdout and stderr."""
    if not silence:
        yield
        return
    buf_out, buf_err = io.StringIO(), io.StringIO()
    with contextlib.redirect_stdout(buf_out), contextlib.redirect_stderr(buf_err):
        yield


def _set_pcntoolkit_log_level(level=logging.ERROR):
    """Set pcntoolkit-related loggers to a given level."""
    for name in ("pcntoolkit", "pcntoolkit.core", "pcntoolkit.normative"):
        logging.getLogger(name).setLevel(level)


def pcn_wscore(
    X,
    Y,
    idx_hp,
    verbose=True,
    progress_every=50,
    use_bspline=False,
    bspline_nknots=4,
    bspline_degree=3,
    poly_degree=3,
    suppress_warnings=True,
    silence_internal=True,   # key: suppress internal stdout/stderr
    **kwargs
):
    if suppress_warnings:
        warnings.filterwarnings("ignore")

    _set_pcntoolkit_log_level(logging.ERROR)

    X = _to_2d(X)
    Y = _to_2d(Y)
    idx_hp = np.asarray(idx_hp).astype(bool).reshape(-1)

    if Y.shape[0] != X.shape[0] and Y.shape[1] == X.shape[0]:
        Y = Y.T

    N, K = X.shape
    Ny, P = Y.shape

    if Ny != N:
        raise ValueError(f"Sample size mismatch: X={N}, Y={Ny}")
    if idx_hp.shape[0] != N:
        raise ValueError(f"idx_hp length should be {N}, got {idx_hp.shape[0]}")

    Xs, mu, std, zero = _standardize_hp(X, idx_hp)
    if np.any(zero) and verbose:
        zero_idx = np.where(zero)[0].tolist()
        print(f"[wscore] Constant covariate columns within HP set to zero: {zero_idx}")

    Z = np.zeros((N, P), dtype=float)
    R = np.zeros((N, P), dtype=float)
    Yhat = np.zeros((N, P), dtype=float)

    subj = np.arange(N)
    be_all = np.zeros((N, 1), dtype=int)
    be_tr = be_all[idx_hp]

    basis_mean = (
        BsplineBasisFunction(
            basis_column=list(range(K)),
            nknots=bspline_nknots,
            degree=bspline_degree,
        )
        if use_bspline
        else PolynomialBasisFunction(
            basis_column=list(range(K)),
            degree=poly_degree,
        )
    )

    try:
        pe = int(progress_every)
    except Exception:
        pe = 50
    if pe <= 0:
        pe = 50

    for j in range(P):
        y = Y[:, j].astype(float)

        if np.nanstd(y[idx_hp], ddof=1) < EPS:
            yhat = np.full(N, float(np.nanmean(y[idx_hp])))
            resid = y - yhat
            sd_hp = np.nanstd(resid[idx_hp], ddof=1)
            z = resid / sd_hp if sd_hp > EPS else np.zeros_like(resid)
            R[:, j], Z[:, j], Yhat[:, j] = resid, z, yhat
        else:
            train = NormData.from_ndarrays(
                name="train",
                X=Xs[idx_hp],
                Y=y[idx_hp][:, None],
                batch_effects=be_tr,
                subject_ids=subj[idx_hp],
            )
            test = NormData.from_ndarrays(
                name="test",
                X=Xs,
                Y=y[:, None],
                batch_effects=be_all,
                subject_ids=subj,
            )

            template = BLR(
                heteroskedastic=False,
                basis_function_mean=basis_mean,
                basis_function_var=None,
            )

            model = NormativeModel(
                template,
                savemodel=False,
                evaluate_model=False,
                saveresults=False,
                saveplots=False,
                inscaler="none",
                outscaler="none",
            )

            with _silence_stdout_stderr(silence_internal):
                pred = model.fit_predict(train, test)

            try:
                z = _get_var(pred, ["Z", "z", "zscore", "Zscores"])
            except KeyError:
                with _silence_stdout_stderr(silence_internal):
                    pred2 = model.predict(test)
                z = _get_var(pred2, ["Z", "z", "zscore", "Zscores"])

            try:
                yhat = _get_var(pred, ["yhat", "ymu", "Yhat", "Ymu"])
            except KeyError:
                with _silence_stdout_stderr(silence_internal):
                    pred2 = model.predict(test)
                yhat = _get_var(pred2, ["yhat", "ymu", "Yhat", "Ymu"])

            yhat = yhat.reshape(-1)
            z = z.reshape(-1)
            R[:, j], Z[:, j], Yhat[:, j] = (y - yhat), z, yhat

        if verbose and ((j + 1) % pe == 0 or j == P - 1):
            print(f"[wscore] Completed {j + 1}/{P}")

    return Z, R, Yhat
