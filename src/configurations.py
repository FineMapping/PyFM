import numpy as np
from copy import deepcopy
import time
from collections import defaultdict
from tqdm import tqdm
import os
class Configurations:
    def __init__(self, **kwarg):
        self.max_causal = kwarg["max_causal"]
        self.m = kwarg["m"]
        self.visited_config_scores = dict() # stores scores of configs already visited
        self.score_config = kwarg[
            "score_config"
        ]  # method to calculate score for a config
        self.optimization_params = kwarg["optimization_params"]
        self.outdir = kwarg["outdir"]
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
            if self.n_causal > prev_n_causal: # if n_causal increments, log some stats
                print(
                    f"Expored models with {prev_n_causal} causal variant in {time.time() - tic} seconds"
                )
                print("n_causal", prev_n_causal)
                print("best config", self.best_config)
                print("best score", self.best_score, end='\n\n')
                tic = time.time()
                prev_n_causal = self.n_causal
            if do is not None:
                do(self)
            self.next()

        print(
            f"Expored models with {self.n_causal} causal variant in {time.time() - tic} seconds"
        )
        print("n_causal", self.n_causal)
        print("best config", self.best_config)
        print("best score", self.best_score, end='\n\n')


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


class HashableArray(np.ndarray):
    # https://stackoverflow.com/questions/26821169/python-how-to-wrap-the-np-ndarray-class
    def __new__(cls, *args, **kwargs):
        this = np.array(*args, **kwargs)
        this = np.asarray(this).view(cls)
        return this

    def __array_finalize__(self, obj):
        if obj.ndim == 4:
            self.myShape = (self.shape[0]*self.shape[1], self.shape[2]*self.shape[3])
        else:
            self.myShape = self.shape

    def hash(self):
        # https://stackoverflow.com/questions/16589791/most-efficient-property-to-hash-for-numpy-array
        hash_value = sum(map(lambda x: 1 << x,self.data))
        return hash_value

    def __hash__(self):
        return self.hash()

    def __eq__(self, x):
        return self.hash() == x.hash()


class SSSConfigurations(Configurations):
    """
    Shotgun Stochastic Search, such as used by FINEMAP
    """

    def __init__(self, **kwarg):
        super().__init__(**kwarg)
        current_causal = (self.max_causal//2) if self.max_causal > 1 else 1
        self.current_model = HashableArray(np.arange(1,current_causal+1))
        self.not_included = set(range(current_causal+1, self.m+1))
        
        self.best_config = self.current_model.copy()
        self.best_score = self.score_config(self.current_model)
        # Optimization Parameters
        self.num_iterations = self.optimization_params['SSS_iterations']
        self.alpha1 = self.optimization_params['SSS_alpha1']
        self.alpha2 = self.optimization_params['SSS_alpha2']

        np.random.seed(self.optimization_params['Random_Seed'])

    def get_negative_neighborhood(self):
        if not len(self.current_model) - 1:
            return []
        copy = self.current_model.copy()[:-1]
        yield copy
        for i in range(len(self.current_model)-1):
            copy[i] = self.current_model[-1]
            yield copy
            copy[i] = self.current_model[i]

    def get_same_neighborhood(self):
        copy = self.current_model.copy()
        if not len(self.current_model):
            return []
        else:
            for snp in self.not_included:
                for i in range(len(self.current_model)):
                    copy[i] = snp
                    yield copy
                    copy[i] = self.current_model[i]
                    
    def get_positive_neighborhood(self):
        if len(self.current_model) == self.max_causal:
            return []
        copy = HashableArray(np.concatenate([self.current_model, np.array([-1])]))
        for snp in self.not_included:
            copy[-1] = snp
            yield copy
    
    def sample_from_distribution(self, generator):
        neighbor_scores = np.array([self.score_config(x) for x in generator()])
        if not len(neighbor_scores):
            return None, None
        best_config_index = np.argmax(neighbor_scores)
        
        # Weighted by score
        weights = np.power(neighbor_scores, self.alpha1)
        
        random_model_index = np.random.choice(np.arange(len(neighbor_scores)), p=weights/weights.sum())        
        models = {i: x.copy() for i, x in enumerate(generator()) if i in (best_config_index,random_model_index)}
        if neighbor_scores.max() > self.best_score:
            self.best_score = neighbor_scores.max()
            self.best_config = models[best_config_index].copy()
        return models[random_model_index], neighbor_scores[random_model_index]
        
    def search(self):
        tic = time.time()
        for _ in tqdm(range(self.num_iterations)):
            neighbors = list(map(self.sample_from_distribution, 
                (self.get_negative_neighborhood, self.get_same_neighborhood, self.get_positive_neighborhood)
            ))
            weights = np.power(np.array([neighbor[1] for neighbor in neighbors if neighbor[0] is not None]), self.alpha2)
            neighbors = [neighbor[0] for neighbor in neighbors if neighbor[0] is not None]
            next_model_index = np.random.choice(np.arange(len(neighbors)), p=weights/weights.sum())
            
            self.current_model = neighbors[next_model_index]
            self.not_included = set(np.arange(1, self.m+1)) - set(self.current_model)

        print(
            f"Expored models up to {self.max_causal} causal variant in {time.time() - tic} seconds"
        )
        print("n_causal", len(self.best_config))
        print("best config", np.array(sorted(self.best_config)))
        print("best score", self.best_score)

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
