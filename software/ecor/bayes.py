import numpy as np 
import scipy.stats
import scipy.optimize

def steady_state_growth_rate_log_posterior(params, time, od, neg=True):
    """
    Computes and returns the log posterior for steady state growth. 

    Parameters
    ----------
    params : list
        The parameters that define the steady state growth rate. They should 
        be provided in the following order:

        lam : float, bounded [0, inf]
            The growth rate in units of inverse time.
        od_init : float, bounded [0, inf]
            The initial optical density in arbitrary units. 
        sigma : float, bouned [0, inf]
            The homoscedastic error.
    
    time : float or numpy array, bounded [0, inf]
        The time range over which the posterior should be evaluated.
    od : float or numpy array, bounded [0, inf]
        The observed optical density in arbitrary units.
    neg : bool, optional
        If True, the negative log posterior is returned.

    Returns 
    --------
    log_post : float
        The value of the log posterior for the provided parameters and data.

    Notes 
    ------
    The posterior distribution is composed of the product of likelihood function 
    and the prior distributions for each parameter. The log posterior is the 
    sum of the likelihood and prior distributions. Here, the likelihood function
    is chosen to be a normal distribution parameterized as:

    .. math::
        f(OD, t \\vert \lambda, OD_0, \sigma) \propto \\frac{1}{\sqrt{2\pi\sigma^2}\exp\left[-\\frac{\left(OD - OD^*\right)^2}{2\sigma^2}}\right]

    where :math:`OD^*` is the theoretically predicted optical density in
    steady-state exponential growth given the parameters :math:`\lambda` and 
    :math:`OD_0`, 

    .. math::
        OD^*(\lambda, OD_0\,\\vert\,t) = OD_0\exp\left[\lambda t\right].
    
    The prior distributions for all three parameters :math:`\theta` 
    (:math:`[\lambda, OD_0, \sigma] \in \theta` ):math:`\lambda`, :math:`OD_0`,
    and :math:`\sigma` are taken as gamma distributions parameterized as 

    .. math::
        g(\theta\,\\vert\,\alpha, \beta) = \\frac{1}{\Gamma(\alpha)}\\frac{(\beta\theta)^\alpha}{\theta}\exp\left[-\beta\theta\right].
    
    In this case, the parameters :math:`\alpha`, and :math:`\beta` are chosen to 
    be 2.20 and 5, respectively. This choice results in 95% of the density 
    lying between 0 and 1 for the parameter of interest. 
    """

    def _theoretical_od(time, lam, od_init):
        """
        Computes the theoretical optical density as a function of time 
        given the growth rate (lam) and the initial optical density (od_init)
        """
        return od_init * np.exp(lam * time)

    def _log_prior(lam, od_init, sigma):
        """
        Defines the log prior for the growth rate (lam), initial optical density
        (od_init), and the homoscedastic error (sigma) as a gamma distribution. 
        """
        lam_lp = scipy.stats.gamma.logpdf(lam, 2.20, loc=0, scale=0.2)
        sig_lp = scipy.stats.gamma.logpdf(sigma, 2.20, loc=0, scale=0.2)
        od_init_lp = scipy.stats.gamma.logpdf(od_init, 2.20, loc=0, scale=0.2)
        return lam_lp + sig_lp + od_init_lp

    def _log_likelihood(time, od, lam, od_init, sigma):
        """
        Defines the log likelihood for the growth rate. Likelihood function is 
        assumed to be a normal distribution with a mean defined by the 
        theoretical exponential growth and a homoscedastic error.
        """
        mu = _theoretical_od(time, lam, od_init)
        return np.sum(scipy.stats.norm.logpdf(od, mu, sigma))

    # Determine the prefactor to be returned
    if neg:
        prefactor = -1
    else:
        prefactor = 1  

    # unpack the parameters
    lam, od_init, sigma = params

    # Ensure that the bounds are obeyed 
    if (np.array(params) < 0).any():
        return prefactor * (-np.inf)

    # Compute the prior
    lp = _log_prior(lam, od_init, sigma)
    if lp == -np.inf:
        return lp

    # Compute the likelihood
    like = _log_likelihood(time, od, lam, od_init, sigma)
    log_post = prefactor * (lp + like)
    return log_post
    
def compute_percentiles(values, 
                        percs={95:[2.5, 97.5],
                               75: [12.5, 87.5],
                               50: [25, 75],
                               25: [37.5, 62.5],
                               5: [45, 55]}):
    """
    Computes the percentiles of a given array.

    Parameters 
    -----------
    values: 1d numpy array 
        The values of which the percentiles will be computed.
    percs: dict, optional
        The desired percentiles as a dictionary with the percentile name 
        and the upper and lower bounds. Default computes the 95th, 75th, 50th, 
        25th, and 5th percentiles of the data

    Returns
    ------- 
    perc_dict : dict
        A dictionary of the percentiles with the lower and upper bounds.
    """
    perc_dict = {}
    for k, v in percs.items():
        perc_dict[k] = np.percentile(values, v)
    return perc_dict
        
# def compute_growth_curve_fit(samples, 
#                              time_range, 
#                              params={'mu':'mu', 
#                                     'od_init':'od_init'},
#                              percs =                             ):
#     mu = samples[params['mu']].values
#     od_init = samples[params['od_init']].values
#     perc_dfs = []
#     for k, v in percs.items():     
#         cred_region = np.zeros((2, len(time_range)))
#         for i, t in enumerate(time_range):
#             _fit = od_init * np.exp(mu * t)
#             cred_region[:, i] = np.percentile(_fit, (v[0], v[1]))
#         _df = pd.DataFrame(cred_region.T, columns=['low', 'high'])
#         _df['time'] = time_range
#         _df['percentile'] = k
#         perc_dfs.append(_df)
#     df = pd.concat(perc_dfs, sort=False)
#     return df
