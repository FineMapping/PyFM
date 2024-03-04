# TODO: Add citations to Caviarbf
#!/usr/bin/env python3
import argparse
import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt

# from utils import
from data import Data
from bf import calculate_BFs, calculate_scores
from configurations import CONFIGS_IMPL_MAP
from priors import calculate_priors
from selection import get_rhos, get_marginal_pip, print_marginal_pips, print_rhos


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-z",
        "--zfile",
        help="Marginal test statistics with 3 columns - rsid, standardized effect size, (optional) variance of variant",
        required=True,
    )
    parser.add_argument("-r", "--rfile", help="LD scores matrix", required=True)
    parser.add_argument("-o", "--outdir", required=True)
    parser.add_argument(
        "-c",
        "--max-causal",
        help="Maximum number of causal variants (inclusive)",
        default=3,
    )
    parser.add_argument("-n", "--sample-number", default=0)
    parser.add_argument("-t", "--prior-type", default=0)
    parser.add_argument("-a", "--prior-values", default="0.1 0.2 0.4 0.8 1.6")
    parser.add_argument(
        "-e",
        "--epsilon",
        help="the value added to the diagonal of the correlation matrix.",
        default=0,
    )
    parser.add_argument(
        "--approx-bf",
        help="Calculate approx BF (recommended for Binary trait) instead of exact BF (recommended for Quant. trait)",
        action="store_true",
    )
    parser.add_argument(
        "--configs-method",
        help=f"choose one causal configurations exploration method from {CONFIGS_IMPL_MAP.keys()}",
        default="AllConfigurations",
    )
    parser.add_argument(
        "--rho",
        help="the cutoff for rho-confidence set",
        default=0.95,
    )
    params = parser.parse_args()
    # TODO: Move the code below into a function that can be called if this project is imported as a python library
    n = int(params.sample_number)
    approx_bf = params.approx_bf

    # TODO: make sure zfile and rfile exist
    data = Data(params.zfile, params.rfile, n, approx_bf)

    # TODO: make sure outdir exists otherwise make folders and subfolders as needed
    outdir = params.outdir
    prior_type = bool(params.prior_type)
    pve_for_prior = int(params.prior_type) > 0 # switch arg for prior_type
    prior_values = np.array([float(x) for x in params.prior_values.split()])
    e = float(params.epsilon)

    if not pve_for_prior:
        prior_values = prior_values**2  # sigmaa^2
    else:
        prior_values = prior_values / (1 - prior_values)

    max_causal = int(params.max_causal)

    configs_method = params.configs_method

    rho = float(params.rho)
    return (
        data,
        n,
        pve_for_prior,
        prior_values,
        e,
        max_causal,
        configs_method,
        approx_bf,
        rho,
        outdir,
    )


if __name__ == "__main__":
    (
        data,
        n,
        pve_for_prior,
        prior_values,
        e,
        max_causal,
        configs_method,
        approx_bf,
        rho,
        outdir,
    ) = parse_args()

    # set scores = log(BF)
    config_scores, best_configs = calculate_BFs(
        data,
        n,
        pve_for_prior,
        prior_values,
        e,
        max_causal,
        configs_method,
        approx_bf,
        os.path.join(outdir, "BF.tsv"),
    )

    # calculate priors
    n_causal2log_prior = calculate_priors(data.m, max_causal)

    # set scores = log (BF x prior/max_BF)
    max_BF = max([config_scores[k][tuple(best_configs[k])] for k in config_scores])
    config_scores, total_score = calculate_scores(
        config_scores, n_causal2log_prior, max_BF
    )

    # set null_model_score
    null_model_score = 10 ** (
        n_causal2log_prior[0] - max_BF
    )  # BF = 1, log_BF = 0 for null model
    total_score += null_model_score
    print("totalscore = ", total_score)

    # get marginal PIP for each variant
    marginal_pip, ranking = get_marginal_pip(data.m, config_scores, total_score)

    print_marginal_pips(marginal_pip, ranking, data, os.path.join(outdir, "pips.tsv"))

    # rho-confidence set
    rho_scores = get_rhos(ranking, config_scores, max_causal)
    print_rhos(rho_scores, ranking, data, os.path.join(outdir, "rhos.tsv"))
