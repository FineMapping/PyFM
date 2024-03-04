import numpy as np
from copy import deepcopy
import time
from collections import defaultdict

class Configurations:
    def __init__(self, **kwarg):
        self.max_causal = kwarg["max_causal"]
        self.m = kwarg["m"]
        self.visited_config_scores = dict() # stores scores of configs already visited
        self.score_config = kwarg[
            "score_config"
        ]  # method to calculate score for a config
        # child class should initialize:
        self.config = None
        self.current_score = None

    def next(self):
        """
        shift self.config to the next configuration or None if no next config
        """
        raise NotImplementedError

    def ended(self):
        return self.config is None

    def search(self, do=None):
        """
        keep exploring configurations till the end
        params:
        do - to be called on the config after every iteration
        """
        tic = time.time()
        prev_n_causal = self.n_causal
        while not self.ended():
            if self.n_causal > prev_n_causal:
                print(
                    f"Expored models with {self.n_causal} causal variant in {time.time() - tic} seconds"
                )
                print("n_causal", self.n_causal)
                print("best config", self.best_config)
                print("best score", self.best_score)
                tic = time.time()
                prev_n_causal = self.n_causal
            if do is not None:
                do(self)
            self.next()

    def get_scores_by_n_causal(self):
        scores_by_n_causal = defaultdict(dict)
        for config,score in self.visited_config_scores.items():
            scores_by_n_causal[len(config)][config] = score
        return scores_by_n_causal


class AllConfigurations(Configurations):
    """
    explore all unique configurations for given n_causal
    """

    def __init__(self, **kwarg):
        """
        initialize vars, including calculating score of the first config
        """
        super().__init__(**kwarg)
        self.n_causal = 1
        self.config = np.arange(
            1, self.n_causal + 1
        )  # list of indices of selected variants
        self.current_score = self.score_config(self.config)
        self.best_config = deepcopy(self.config)
        self.best_score = self.current_score
        self.visited_config_scores[tuple(self.config)] = self.current_score

    def next(self):
        """
        moves from config (..., p-1, n-k-1, n-k-2, ... n-2, n-1, n) to (..., p, n-k-1, n-k-2, ... n-2, n-1, n) till p <= n-k
        where n = self.n_causal, k = self.m
        """
        i = (
            self.n_causal
        )  # 1-based index at which config will change, i.e index of p-1 -> p
        while i > 0 and self.config[i - 1] - i == self.m - self.n_causal:
            i -= 1
        if i == 0:
            # no next config
            if self.n_causal < self.max_causal:
                self.n_causal += 1
                self.config = np.arange(1, self.n_causal + 1)
            else:
                self.config = None
        else:
            self.config[i - 1 :] = (
                self.config[i - 1] + 1 + np.arange(self.n_causal - i + 1)
            )
            self.current_score = self.score_config(self.config)
            self.visited_config_scores[tuple(self.config)] = self.current_score
            if self.current_score > self.best_score:
                self.best_config = deepcopy(self.config)
                self.best_score = self.current_score


# TODO: implement next for this class
class SSSConfigurations(Configurations):
    """
    Shotgun Stochastic Search, such as used by FINEMAP
    """

    def __init__(self, **kwarg):
        super().__init__(**kwarg)


#############################################################
#        Configurations Factory
#############################################################

CONFIGS_IMPL_MAP = {
    "AllConfigurations": AllConfigurations,
    "SSSConfigurations": SSSConfigurations,
}


def ConfigurationsFactory(configs_impl):
    if configs_impl not in CONFIGS_IMPL_MAP:
        raise ValueError(
            f"Must choose one configurations implementation from {CONFIGS_IMPL_MAP.keys()}"
        )
    else:
        return CONFIGS_IMPL_MAP[configs_impl]
