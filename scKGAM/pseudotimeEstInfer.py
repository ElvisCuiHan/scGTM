import numpy as np
from pseudotimeAPI import *
import pyswarms as ps

def estimation(y, t, marginal, iter=50):

    n = 30
    if marginal in ["Poisson", "ZIP"]:
        d = 6
        bounds = [np.array([np.log((np.min(y + 1))) , -np.inf, -np.inf, t.min(), -np.inf, -np.inf]),
                  np.array([np.log((np.max(y))) , np.inf, np.inf, t.max(), np.inf, np.inf])]
    elif marginal in ["NB", "ZINB"]:
        d = 7
        bounds = [np.array([np.log(np.min(y) + 1) , -np.inf, -np.inf, t.min() - 0.5, 1, 0, 0]),
                  np.array([np.log(np.max(y) + 1) , np.inf, np.inf, t.max() + 0.5, 100, 100, 100])]
    else:
        raise ValueError("Enter a valid marginal distribution: [NB, ZINB, Poisson, ZIP]!")

    np.random.seed(123)
    b = np.random.random((n, d))

    # Set-up hyperparameters and correct initial position
    options = {'c1': 1.2, 'c2': 0.3, 'w': 0.9}
    b[:, 0] = np.log(np.mean(y) + 1)
    b[:, -1] = 0.1
    b[:, 1] += 5
    b[:, 2] += 5
    if d == 7:
        b[:, -3] += 1

    # Call instance of PSO
    gmodel = ps.single.GlobalBestPSO(n_particles=n, dimensions=b.shape[1], options=options, bounds=bounds, init_pos=b)

    # Perform optimization
    gcost, gbest = gmodel.optimize(pso_obj_fct, iters=iter, y=y, t=t, marginal=marginal)

    return [gcost, gbest]

def inference(t, para, marginal):
    """
    This function outputs both the estimated Fisher information matrix
    and the estimated inverse Fisher information matrix for the first
    four parameters of the model: t0, k1, k2 and mu.

    The latter is also known as the covariance matrix in likelihood theory.
    :param t: Pseudotime that are used to estimate parameters.

    :param para: Estimated parameters given by "estimation" function.
    :param marginal: The marginal distribution of the model.
    :return:
        fisher: The 4x4 estimated Fisher information matrix.
        Cov: The 4x4 estimated covariance matrix of four parameters.
        lower: Lower bound of the 95% confidence interval.
        Upper: Upper bound of the 95% confidence interval.
    """

    fisher = Fisher_info(t, para, marginal)
    try:
        cov = np.linalg.inv(fisher)
        lower = np.round(para[3] - 1.96 * np.sqrt(cov[0, 0]), 3)
        upper = np.round(para[3] + 1.96 * np.sqrt(cov[0, 0]), 3)
    except:
        print("The Fisher information matrix is singular: estimated t0 is either less than 0 or greater than 1.\n",
              "Calculating the sub-Fisher information matrix of t0 instead.\n")
        cov = 1 / fisher[0, 0]
        lower = np.round(para[3] - 1.96 * np.sqrt(cov), 3)
        upper = np.round(para[3] + 1.96 * np.sqrt(cov), 3)
    return [fisher, cov, lower, upper]