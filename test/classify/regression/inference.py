from typing import Dict, List, Set

import pymc as pm
import numpy as np
import arviz as az


def formulate_matrices(matches: Dict[str, List[int]], counts: Dict[str, int]):
    rfx: Set[int] = set()
    for lst in matches.values():
        rfx |= set(lst)
    idxmap = dict(zip(sorted(rfx), range(len(rfx))))
    
    X = np.zeros(shape=(len(matches), len(idxmap)))
    y = np.zeros(shape=(len(matches)))
    
    for x, (kmer, refs) in enumerate(matches.items()):
        for ref in refs:
            X[x, idxmap[ref]] = 1
        y[x] = counts[kmer]
            
    return X, y

def estimate_abundances(X, y):
    # X: kmer × genome matrix
    # y: observed k-mer counts

    with pm.Model() as model:

        # composition (sums to 1)
        theta = pm.Dirichlet("theta", a=np.ones(X.shape[1]))

        # sequencing depth
        lam = pm.Gamma("lam", alpha=2, beta=1)

        # global error
        epsilon = pm.Gamma("epsilon", alpha=1, beta=10)

        # dispersion
        phi = pm.Gamma("phi", alpha=2, beta=1)

        # expected counts
        mu = lam * pm.math.dot(X, theta) + epsilon

        # likelihood
        y_obs = pm.NegativeBinomial("y_obs", mu=mu, alpha=phi, observed=y)

        trace = pm.sample(1000, tune=1000)
        
        # parameter estimates
        summary_df = az.summary(trace, ci_prob=0.95)
        
        return summary_df
