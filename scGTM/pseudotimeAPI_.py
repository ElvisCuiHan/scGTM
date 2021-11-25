import numpy as np
from scipy.stats import nbinom, poisson
import matplotlib.pyplot as plt

## Shared parameters log-bell
def link(t, mu, k1, k2, t0):
    """

    :param t: pseudotime
    :param mu: peak expression value
    :param k1: activation strength or how quickly a gene is up regulated (increasing)
    :param k2: activation strength or how quickly a gene is down regulated (decreasing)
    :param t0: turning point
    :return:
    """
    part1 = mu * np.exp(- np.abs(k1) * (t - t0) ** 2)
    part2 = mu * np.exp(- np.abs(k2) * (t - t0) ** 2)

    return part1 * (t <= t0) + part2 * (t > t0)

## Poisson
def single_gene_log_likelihood_Poisson(y, t, mu, k1, k2, t0, x=0, flag=False):
    """

    :param y: observation
    :param t: pseudotime
    :param mu: peak expression value
    :param k1: activation strength or how quickly a gene is up regulated (increasing)
    :param k2: activation strength or how quickly a gene is down regulated (decreasing)
    :param t0: turning point
    :return:
    """
    bell = link(t, mu, k1, k2, t0)
    if flag:
        bell = - bell + x
    else:
        bell -= 0.01*x
    mut = np.maximum(np.exp(bell) , 0.1)
    return -np.log(poisson.pmf(y, mut) + 1e-300).sum()

## Zero-inflated Poisson
def single_gene_log_likelihood_ZIP(y, t, mu, k1, k2, t0, alpha, beta, x=0, flag=False):
    """

    :param y:
    :param t:
    :param mu:
    :param k1:
    :param k2:
    :param t0:
    :param p:
    :return:
    """
    bell = link(t, mu, k1, k2, t0)
    if flag:
        bell = - bell + x
    else:
        bell -= 0.01*x
    mut = np.maximum(np.exp(bell), 0.1)
    cache = poisson.pmf(y, mut) + 1e-300

    ## Zero-inflation
    p = 1 / (1 + np.exp(alpha * np.log(mut) + beta))

    return - np.log(cache * (1 - p) + p * (y == 0)).sum()

## Negative Binomial
def single_gene_log_likelihood_NB(y, t, mu, k1, k2, t0, phi, x=0, flag=False):
    """

    :param y:
    :param t:
    :param mu:
    :param k1:
    :param k2:
    :param t0:
    :param phi:
    :return:
    """
    bell = link(t, mu, k1, k2, t0)
    if flag:
        bell = - bell + x
    else:
        bell -= 0.01*x
    mut = np.maximum(np.exp(bell) , 0.1)
    phi = np.maximum(np.floor(phi), 1)
    p0 = mut / (mut + phi)
    cache = nbinom.pmf(y, phi, 1 - p0) + 1e-300
    return -np.log(cache).sum()

## Zero-inflated Negative Binomial
def single_gene_log_likelihood_ZINB(y, t, mu, k1, k2, t0, phi, alpha, beta, x=0, flag=False):
    """

    :param y:
    :param t:
    :param mu:
    :param k1:
    :param k2:
    :param t0:
    :param phi:
    :param p:
    :return:
    """
    bell = link(t, mu, k1, k2, t0)
    if flag:
        bell = - bell + x
    else:
        bell -= 0.01*x
    mut = np.maximum(np.exp(bell) , 0.1)
    phi = np.maximum(np.floor(phi), 1)
    p0 = phi / (mut + phi)
    cache = nbinom.pmf(y, phi, p0) + 1e-300

    ## Zero-inflation
    p = 1 / (1 + np.exp(alpha * np.log(mut) + beta))

    return - np.log(cache * (1 - p) + p * (y == 0)).sum()

def pso_obj_fct(b, **kwargs):
    """

    :param b: an Nxd particle matrix where N is the number of particles,
              and d is the dimension of each particle.
    :param kwargs:
    :return:
    """
    N, d = b.shape
    y, t, marginal = kwargs.values()

    cost = np.zeros(N)
    for i in range(N):
        if marginal == "ZINB":
            mu, k1, k2, t0, phi, alpha, beta, x = b[i, :]
            cost[i] = single_gene_log_likelihood_ZINB(y, t, mu, k1, k2, t0, phi, alpha, beta, x, flag=False)
        elif marginal == "NB":
            mu, k1, k2, t0, phi, alpha, beta, x = b[i, :]
            cost[i] = single_gene_log_likelihood_NB(y, t, mu, k1, k2, t0, phi, x, flag=False)
        elif marginal == "ZIP":
            mu, k1, k2, t0, alpha, beta, x = b[i, :]
            cost[i] = single_gene_log_likelihood_ZIP(y, t, mu, k1, k2, t0, alpha, beta, x, flag=False)
        else:
            mu, k1, k2, t0, alpha, beta, x = b[i, :]
            cost[i] = single_gene_log_likelihood_Poisson(y, t, mu, k1, k2, t0, x, flag=True)
    return cost

def pso_obj_fct_valley(b, **kwargs):
    """

    :param b: an Nxd particle matrix where N is the number of particles,
              and d is the dimension of each particle.
    :param kwargs:
    :return:
    """
    N, d = b.shape
    y, t, marginal = kwargs.values()

    cost = np.zeros(N)
    for i in range(N):
        if marginal == "ZINB":
            mu, k1, k2, t0, phi, alpha, beta, x = b[i, :]
            cost[i] = single_gene_log_likelihood_ZINB(y, t, mu, k1, k2, t0, phi, alpha, beta, x, flag=True)
        elif marginal == "NB":
            mu, k1, k2, t0, phi, alpha, beta, x = b[i, :]
            cost[i] = single_gene_log_likelihood_NB(y, t, mu, k1, k2, t0, phi, x, flag=True)
        elif marginal == "ZIP":
            mu, k1, k2, t0, alpha, beta, x = b[i, :]
            cost[i] = single_gene_log_likelihood_ZIP(y, t, mu, k1, k2, t0, alpha, beta, x, flag=True)
        else:
            mu, k1, k2, t0, alpha, beta, x = b[i, :]
            cost[i] = single_gene_log_likelihood_Poisson(y, t, mu, k1, k2, t0, x, flag=True)
    return cost

## Plot the results
def plot_result(para, t, color, marginal, flag, y1):
    mu_fit, k1_fit, k2_fit, t0_fit = para[:4]
    x_fit = para[-1]
    log_mut_fit = link(np.sort(t), mu_fit, k1_fit, k2_fit, t0_fit)

    if flag:
        log_mut_fit = - log_mut_fit + x_fit
    else:
        log_mut_fit -= x_fit

    p_fit = 1 / (1 + np.exp(para[-1] + para[-2] * np.exp(log_mut_fit)))

    ZIlog_mut_fit = np.maximum(log_mut_fit + np.log(1 - p_fit), 0)

    plt.plot(np.sort(t), log_mut_fit, linewidth=3, c=color[0], label="Fitted curve (mean function)")
    if marginal in ["ZIP", "ZINB"]:
        plt.plot(np.sort(t), ZIlog_mut_fit - 0.1, linewidth=1.2, c=color[1], label="With dropout (minus -log(1-p))")
    plt.xlabel("Pseudotime", fontsize=20)
    plt.ylabel("Expression log(count+1)", fontsize=20)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    if t0_fit <= 1 and t0_fit >= 0:
        plt.axvline(t0_fit, linewidth=4, color=color[2], linestyle="--",
                    ymin=0.02, ymax=0.96, label="Estimated t0 (turning point)")
    #plt.legend()

    if marginal in ["ZIP", "ZINB"]:
        plt.twinx()
        plt.ylabel("Dropout Rate", c=color[3], fontsize=14)
        plt.plot(np.sort(t), p_fit, linewidth=1.5, c=color[3], label="Dropout rate")

def Fisher_info(t, para, marginal):
    """
    This function calculates the Fisher information matrix of parameters of interest.
    Specifically, mu, k1, k2 and t0. Zero-inflation parameter p will also be included
    if the marginal is either ZIP or ZINB.
    :param t: Pseudotime, should be a vector.
    :param para: Estimated parameters, should be a vector.
    :param marginal: Marginal distribution, should be one of ["Poisson", "NB", "ZIP", "ZINB"].
    :return: a 4x4 Fisher information matrix.
    """
    if marginal == "ZIP" or marginal == "Poisson":
        mu_fit, k1_fit, k2_fit, t0_fit, alpha_fit, beta_fit, x_fit = para
    else:
        mu_fit, k1_fit, k2_fit, t0_fit, phi_fit, alpha_fit, beta_fit, x_fit = para
    log_mut_fit = link(t, mu_fit, k1_fit, k2_fit, t0_fit)

    mut = np.maximum(np.exp(log_mut_fit), 0.1)

    t0_deri = 2 * (k1_fit * (t <= t0_fit) + k2_fit * (t > t0_fit)) * (t - t0_fit) * log_mut_fit  # * (1 + 1/mut)
    k1_deri = (t - t0_fit) ** 2 * log_mut_fit * (t <= t0_fit)  # * (1 + 1/mut)
    k2_deri = (t - t0_fit) ** 2 * log_mut_fit * (t > t0_fit)  # * (1 + 1/mut)
    mu_deri = log_mut_fit / mu_fit  # * (1 + 1/mut)

    if marginal == "Poisson":
        t0_deri *= (1 + 1 / mut); k1_deri *= (1 + 1 / mut)
        k2_deri *= (1 + 1 / mut); mu_deri *= (1 + 1 / mut)
    elif marginal == "ZIP":
        t0_deri *= (1 + 1 / mut); k1_deri *= (1 + 1 / mut)
        k2_deri *= (1 + 1 / mut); mu_deri *= (1 + 1 / mut)
    else:
        t0_deri *= phi_fit * (mut + 1) / (mut + phi_fit)
        k1_deri *= phi_fit * (mut + 1) / (mut + phi_fit)
        k2_deri *= phi_fit * (mut + 1) / (mut + phi_fit)
        mu_deri *= phi_fit * (mut + 1) / (mut + phi_fit)

    cache = np.vstack([t0_deri, k1_deri, k2_deri, mu_deri])

    return (cache.dot(cache.T)) #/ len(t)