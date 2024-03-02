import numpy as np
from numpy import log10


def calculate_logprior(m, n_causal):
    p = 1.0 / m  # binomial param, prob of one snp being causal
    return log10(p) * n_causal + log10(1 - p) * (m - n_causal)


def calculate_priors(m, max_causal):
    n_causal2log_prior = dict()
    for n_causal in range(max_causal + 1):
        n_causal2log_prior[n_causal] = calculate_logprior(m, n_causal)
    return n_causal2log_prior
