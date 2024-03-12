# TODO: Add citations to Caviarbf
#!/usr/bin/env python3
import argparse
import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import time

# from utils import
from data import Data
from bf import calculate_BFs, calculate_scores
import configurations
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
    parser.add_argument("-a", "--prior-values", default="0.1,0.2,0.4")
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
        "--SSS-iterations",
        help=f"Number of iterations for Shotgun Stochastic Search to run",
        default=100,
    )
    parser.add_argument(
        "--SSS-alpha1",
        help=f"""Temperature Parameter for SSS sampling stage 1. 
        When this parameter is high, the chosen model in each group 
        will be the highest scoring model (greedy search). 

        If this parameter is low, the chosen model in each group will
        be more random.""",
        default=1.5,
    )
    parser.add_argument(
        "--SSS-alpha2",
        help=f"""Temperature Parameter for SSS sampling stage 2. 
        When this parameter is high, the group chosen
        will be the highest scoring model (greedy search). 

        If this parameter is low, the group chosen will
        be more random.""",
        default=1.5,
    )
    parser.add_argument(
        "--Random-Seed",
        help=f"""Random Seed For Shotgun Stochastic Search. Set to 'None' For No Random Seeding""",
        default="None",
    )
    parser.add_argument(
        "-p",
        "--rho",
        default=1,
    )
    params = parser.parse_args()
    # TODO: Move the code below into a function that can be called if this project is imported as a python library
    n = int(params.sample_number)
    approx_bf = params.approx_bf

    # TODO: make sure zfile and rfile exist
    data = Data(params.zfile, params.rfile, n, approx_bf)

    outdir = params.outdir
    os.makedirs(outdir, exist_ok=True)
    prior_type = bool(params.prior_type)
    pve_for_prior = int(params.prior_type) > 0  # switch arg for prior_type
    prior_values = np.array([float(x) for x in params.prior_values.split(",")])
    e = float(params.epsilon)

    if not pve_for_prior:
        prior_values = prior_values**2  # sigmaa^2
    else:
        prior_values = prior_values / (1 - prior_values)

    max_causal = int(params.max_causal)

    configs_method = params.configs_method
    optimization_params = {
        "SSS_iterations": int(params.SSS_iterations),
        "SSS_alpha1": float(params.SSS_alpha1),
        "SSS_alpha2": float(params.SSS_alpha2),
        "Random_Seed": eval(params.Random_Seed),
    }

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
        optimization_params,
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
        optimization_params,
        outdir,
    ) = parse_args()
    tic = time.time()

    # set scores = log(BF)
    config_scores, max_BF = calculate_BFs(
        data,
        n,
        pve_for_prior,
        prior_values,
        e,
        max_causal,
        configs_method,
        approx_bf,
        optimization_params,
        outdir=outdir,
    )

    # calculate priors
    n_causal2log_prior = calculate_priors(data.m, max_causal)

    # set scores = log (BF x prior/max_BF)
    config_scores, total_score = calculate_scores(
        config_scores, n_causal2log_prior, max_BF
    )

    # if configs_method == 'SSSConfigurations':
    #     # Returning early for SSS since we don't have all the scores for each model
    #     exit(0)

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
    rho_scores = get_rhos(rho, ranking, config_scores, max_causal, total_score)
    print_rhos(rho_scores, ranking, data, os.path.join(outdir, "rhos.tsv"))
    print(f"PyFM ran successfully in {round(time.time() - tic,2)} seconds")
