import pandas as pd
import json
import numpy as np
import pyswarms as ps
from scipy.optimize import minimize
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from scipy.stats import nbinom, poisson
from .pseudotimeAPI import *
from .pseudotimeEstInfer import *

import warnings
warnings.filterwarnings("ignore")

def main(gene_index = 100, marginal="ZIP", iter=50, data_dir=None, save_dir=None):
    print("Loading data......")

    ## LOAD DATA
    data = pd.read_csv(data_dir)
    print("Loading finished!")

    ## TAKE NEEDED DATA
    t = data.iloc[:, 1]
    y1 = np.floor(data.iloc[:, gene_index])
    gene_name = data.columns[gene_index]

    ## ESTIMATION
    print("\nWe are estimating gene %d with marginal %s." % (gene_index, marginal))

    result = {}
    gcost, gbest = estimation(y1, t, marginal, iter)
    result['negative_log_likelihood'] = gcost

    if gcost > 1e2:
        print("\nBest negative log-likelihood: ", np.round(gcost, 2), "\n")
    else:
        print("\nAlgorithm fails to find reasonable estimation.\n")

    if marginal == "ZIP":
        result['mu'] = gbest[0]; result['k1'] = gbest[1]
        result['k2'] = gbest[2]; result['t0'] = gbest[3]; result['phi'] = "Nah"
        result['alpha'] = gbest[4]; result['beta'] = gbest[5]


        print("Best parameter estimation:\n",
              "mu , k1 , k2 , t0 , p:\n",
              np.round(gbest, 2), "\n")
    elif marginal == "ZINB":
        gbest[4] = np.maximum(np.floor(gbest[-2]), 1)
        result['mu'] = gbest[0]; result['k1'] = gbest[1]; result['k2'] = gbest[2]
        result['t0'] = gbest[3]; result['phi'] = gbest[4]; result['alpha'] = gbest[5]; result['beta'] = gbest[6]

        print("Best parameter estimation:\n",
              "mu , k1 , k2 , t0 , phi , p:\n",
              np.round(gbest, 2), "\n")
    elif marginal == "Poisson":
        result['mu'] = gbest[0]
        result['k1'] = gbest[1]
        result['k2'] = gbest[2]
        result['t0'] = gbest[3]
        result['phi'] = "Nah"
        result['p'] = "Nah"

        print("Best parameter estimation:\n",
              "mu , k1 , k2 , t0:\n",
              np.round(gbest[:-1], 2), "\n")
    else:
        gbest[-2] = np.maximum(np.floor(gbest[-2]), 1)
        result['mu'] = gbest[0]
        result['k1'] = gbest[1]
        result['k2'] = gbest[2]
        result['t0'] = gbest[3]
        result['phi'] = gbest[4]
        result['p'] = "Nah"

        print("Best parameter estimation:\n",
              "mu , k1 , k2 , t0 , phi:\n",
              np.round(gbest[:-1], 2), "\n")

    ## PLOTTING
    fig, ax = plt.subplots(figsize=(10, 8))
    log_data = np.log(y1 + 1)
    plt.scatter(t, log_data, s=10, c=log_data, cmap=plt.get_cmap('PRGn'))
    plt.ylim(np.min(log_data) - 1, np.max(log_data) + 1)

    plot_result(gbest, t, ['red', 'pink', 'orange', 'darkgreen'], marginal=marginal)

    plt.title("Marginal: " + marginal + ". Gene: " + str(gene_name) + '.', fontsize=15)
    #plt.show()

    plt.savefig(save_dir + str(gene_index) + marginal + ".png")

    ## FISHER INFORMATION
    fisher, var, t0_lower, t0_upper = inference(t, gbest, marginal)
    result['t0_lower'] = t0_lower
    result['t0_upper'] = t0_upper

    #print("Inverse Fisher information matrix of first 4 parameters or t0 alone:\n",
    #      var , "\n")

    ## CONFIDENCE INTERVAL
    print("The 95% confidence interval of the activation time t0:\n" +
          "t0 : (" + str(t0_lower) + ", " + str(t0_upper) + ")")

    result['k1_lower'] = result['k1_upper'] = result['k2_lower'] = \
        result['k2_upper'] = result['mu_lower'] = result['mu_upper'] =  "Nah"

    if np.ndim(var) > 1:
        k1_lower = np.round(gbest[1] - 1.96 * np.sqrt(var[1, 1]), 3)
        k1_upper = np.round(gbest[1] + 1.96 * np.sqrt(var[1, 1]), 3)
        k2_lower = np.round(gbest[2] - 1.96 * np.sqrt(var[2, 2]), 3)
        k2_upper = np.round(gbest[2] + 1.96 * np.sqrt(var[2, 2]), 3)
        mu_lower = np.round(gbest[0] - 1.96 * np.sqrt(var[3, 3]), 3)
        mu_upper = np.round(gbest[0] + 1.96 * np.sqrt(var[3, 3]), 3)
        print("\nThe 95% CIs for activation strength k1 and k2:\n" +
            " k1 : (" + str(k1_lower) + ", " + str(k1_upper) + ")\n",
            "k2 : (" + str(k2_lower) + ", " + str(k2_upper) + ")\n"
             )
        result['k1_lower'] = k1_lower; result['k1_upper'] = k1_upper
        result['k2_lower'] = k2_lower; result['k2_upper'] = k2_upper
        result['mu_lower'] = mu_lower; result['mu_upper'] = mu_upper
        result['Fisher'] = 'Non-singular'
    else:
        var = fisher
        var[0, 0] = 1 / (var[0, 0] + 1e-100)
        var[1, 1] = 1 / (var[1, 1] + 1e-100)
        var[2, 2] = 1 / (var[2, 2] + 1e-100)
        var[3, 3] = 1 / (var[3, 3] + 1e-100)
        k1_lower = np.round(gbest[1] - 1.96 * np.sqrt(var[1, 1]), 3)
        k1_upper = np.round(gbest[1] + 1.96 * np.sqrt(var[1, 1]), 3)
        k2_lower = np.round(gbest[2] - 1.96 * np.sqrt(var[2, 2]), 3)
        k2_upper = np.round(gbest[2] + 1.96 * np.sqrt(var[2, 2]), 3)
        mu_lower = np.round(gbest[0] - 1.96 * np.sqrt(var[3, 3]), 3)
        mu_upper = np.round(gbest[0] + 1.96 * np.sqrt(var[3, 3]), 3)
        print("\nThe 95% CIs for activation strength k1 and k2:\n" +
            " k1 : (" + str(k1_lower) + ", " + str(k1_upper) + ")\n",
            "k2 : (" + str(k2_lower) + ", " + str(k2_upper) + ")\n"
             )
        result['k1_lower'] = k1_lower; result['k1_upper'] = k1_upper
        result['k2_lower'] = k2_lower; result['k2_upper'] = k2_upper
        result['mu_lower'] = mu_lower; result['mu_upper'] = mu_upper
        result['Fisher'] = 'Singular'

    ## SAVE ESTIMATION RESULTS
    with open(save_dir + str(gene_index) + marginal + '.json', 'w') as fp:
        json.dump(result, fp)

    return result
