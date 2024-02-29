# TODO: Add citations to Caviarbf
#!/usr/bin/env python3
import argparse
import pandas as pd
import numpy as np
import os

# from utils import
from Data import Data
from BF import outputBF
from Configurations import CONFIGS_IMPL_MAP


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
    params = parser.parse_args()
    # TODO: Move the code below into a function that can be called if this project is imported as a python library
    n = int(params.sample_number)
    approx_bf = params.approx_bf

    # TODO: make sure zfile and rfile exist
    data = Data(params.zfile, params.rfile, n, approx_bf)

    # TODO: make sure outdir exists otherwise make folders and subfolders as needed
    outdir = params.outdir
    prior_type = bool(params.prior_type)
    pve_for_prior = int(params.prior_type) > 0
    prior_values = np.array([float(x) for x in params.prior_values.split()])
    e = float(params.epsilon)

    if not pve_for_prior:
        prior_values = prior_values**2  # sigmaa^2
    else:
        prior_values = prior_values / (1 - prior_values)

    max_causal = int(params.max_causal)

    configs_method = params.configs_method
    return (
        data,
        n,
        pve_for_prior,
        prior_values,
        e,
        max_causal,
        configs_method,
        approx_bf,
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
        outdir,
    ) = parse_args()
    BF = outputBF(
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
