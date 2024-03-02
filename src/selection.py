from itertools import combinations
import numpy as np
import sys


def get_rhos(ranking, config_scores, max_causal):
    # TODO: use greedy search as described in Caviar and Caviarbf
    # Currently using naive impl but with up to top 5 variants (or max_causal if that's smaller)
    rho_scores = []
    subsets_score = 0
    subsets = []
    for k in range(1, min(max_causal + 1, 6)):
        # k_combinations = list(combinations(rho_set,k))
        # subsets.extend(k_combinations)
        new_subsets = [[ranking[k - 1]]]
        if len(subsets) > 0:
            new_subsets.extend([subset + [ranking[k - 1]] for subset in subsets])
        for new_subset in new_subsets:
            subsets_score += 10 ** config_scores[len(new_subset)][tuple(new_subset)]
        rho_scores.append(subsets_score)
        subsets.extend(new_subsets)
    return rho_scores


def get_marginal_pip(m, config_scores, total_score):
    # for each variant, set marginal prob = sum of model scores contating variant / sum of all model scores
    marginal_pip = np.zeros(m)
    for k in config_scores:
        for config in config_scores[k]:
            config = tuple(config)
            for i in config:
                marginal_pip[i - 1] += 10 ** config_scores[k][config] / total_score
        # rank variants by marginal pip
    ranking = (
        np.argsort(marginal_pip)[::-1] + 1
    )  # 1-based indices of variants in order of marginal prob
    return marginal_pip, ranking


def print_marginal_pips(marginal_pip, ranking, data, filename=None):
    e_pip = int(np.ceil(marginal_pip.sum()))
    f = open(filename, "w") if filename is not None else sys.stdout
    print("# Expected number of causal SNPs: ", marginal_pip.sum(), file=f)
    f.close()
    f = open(filename, "a") if filename is not None else sys.stdout
    print("# Ranked Marginal PIPs:", file=f)
    print("rank\tid\trsid\tpip", file=f)
    for rank in range(e_pip):
        i = ranking[rank]
        print(f"{rank+1}\t{i}\t{data.rsid[i-1]}\t{marginal_pip[i-1]}", file=f)
    if f is not sys.stdout:
        f.close()


def print_rhos(rhos, ranking, data, filename=None):
    f = open(filename, "w") if filename is not None else sys.stdout
    print("# Rho-confidence sets:", file=f)
    f.close()
    f = open(filename, "a") if filename is not None else sys.stdout
    print("step\tid\trsid\trho", file=f)
    for rank in range(len(rhos)):
        i = ranking[rank]
        print(f"{rank+1}\t{i}\t{data.rsid[i-1]}\t{rhos[rank]}", file=f)
    if f is not sys.stdout:
        f.close()
