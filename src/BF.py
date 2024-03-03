import numpy as np
from scipy.linalg import cholesky, cho_solve
from numpy import log, abs, exp
from configurations import ConfigurationsFactory
import time


def calculate_BFs(
    data,
    n,
    pve_for_prior,
    prior_values,
    e,
    max_causal,
    configs_method,
    approx_bf,
    filename,
):

    count = 0
    # TODO: consider putting u in data
    u = (data.pve.values if data.pve is not None else np.ones(data.m)) * (
        n if approx_bf else (n - 1)
    )
    config_BFs = dict()  # map n_causal to all visited config scores
    best_configs = dict()  # map n_causal to all visited config scores
    for n_causal in range(1, max_causal + 1):
        config_iter = ConfigurationsFactory(configs_method)(
            n_causal=n_causal,
            m=data.m,
            score_config=lambda t: config_BF(
                t, data, n, pve_for_prior, prior_values, e, n_causal, approx_bf, u
            ),
        )
        tic = time.time()
        config_iter.search(
            # lambda iter: print("config:", iter.config, "score:", iter.current_score)
        )
        print(
            f"Expored models with {n_causal} causal variant in {time.time() - tic} seconds"
        )
        print("n_causal", n_causal)
        print("best config", config_iter.best_config)
        print("best score", config_iter.best_score)
        config_BFs[n_causal] = config_iter.visited_configs
        best_configs[n_causal] = config_iter.best_config
    return config_BFs, best_configs


def config_BF(t, data, n, pve_for_prior, prior_values, e, n_causal, approx_bf, u):
    if pve_for_prior:
        if data.pve is not None:
            sumOfVarianceSelected = data.pve[t - 1].sum()
            u_sel = u[t - 1] / sumOfVarianceSelected
        else:
            u_sel = np.ones(n_causal) * u[0] / n_causal
    else:
        u_sel = u[t - 1]

    z_sel = data.z[t - 1]
    ld_sel = data.ld[t - 1][:, t - 1]
    logBFs = []  # a BF for each prior_value
    for prior_value in prior_values:
        if approx_bf:
            logBFs.append(approx_BF(z_sel, ld_sel, u_sel * prior_value, e))
        else:
            logBFs.append(exact_BF(z_sel, ld_sel, u_sel * prior_value, n, e))
    mean_logBF = (max(logBFs) + log(exp(np.array(logBFs) - max(logBFs)).mean())) / log(
        10
    )
    return mean_logBF


def exact_BF(z, ld, u, n, e):
    m = len(z)
    for i in range(m):
        ld[i, i] = 1 + (1 / u[i] if u[i] != 0 else 0) + e
    llt_ld = cholesky(ld)
    zTSigmaXz = z.T @ cho_solve((llt_ld, False), z)

    logDet = sum([log(u[i]) + log(abs(2 * llt_ld[i, i])) for i in range(m)])

    logBF = -0.5 * logDet - 0.5 * n * log(1 - zTSigmaXz)

    # print("u: ",u)
    # print("z: ",z)
    # print("ld: ",ld)
    # print("llt_ld: ",llt_ld)
    # print("zTSigmaXz: ",zTSigmaXz)
    # print("logDet: ",logDet)
    # print("logBF: ",logBF)

    return logBF


def approx_BF(z, ld, u, e):
    # TODO: implement approx_BF
    pass


def calculate_scores(config_scores, n_causal2log_prior, max_BF):
    # for each model, set scores = log (BF x prior/max_BF)
    total_score = 0
    for k in config_scores:
        for config in config_scores[k]:
            config = tuple(config)
            config_scores[k][config] += n_causal2log_prior[k] - max_BF
            total_score += 10 ** config_scores[k][config]
    return config_scores, total_score
