"""
Module for handy python functions used foe measuring bacterial growth rates from 
temporal measurements.
"""

import pandas as pd
from .bayes import *
import tqdm
import statsmodels.tools.numdiff as smnd
import matplotlib.pyplot as plt 
import seaborn as sns
import scipy.optimize
import scipy.stats



def infer_growth_rate(data, 
                      od_bounds=None, 
                      convert_time=True, 
                      groupby=None,
                      cols={'time':'clock_time', 'od':'od_600nm'}, 
                      return_opts=True,
                      print_params=True,
                      **kwargs):                  
    """
    Infers the maximal a posteriori (MAP) parameter set for the steady state growth 
    rate given measurements of optical density. This is performed via optimization 
    by minimization.

    Parameters
    ----------
    data : pandas DataFrame
        A tidy long-form pandas DataFrame with columns corresponding to the 
        measurement time and optical density.
    od_bounds : list of floats
        The lower and upper bounds of the optical density range to be considered.
        The default bounds assumed are [0.04, 0.41] inclusive.
    convert_time : bool
        If `True`, the provided time column needs to be converted to elapsed 
        time. In this case, the provided time column is assumed to be 
        the clock time of the measurement and is converted to minutes elapsed.
    groupby : list of str, optional
        The column names for the groupby operation to operate upon. For example,
        if there are multiple strains measured in the data set, a groupby of 
        `['strain']` will yield a growth rate estimate for each strain in 
        the data. A groupby of `['strain', 'replicate']` will return a growth 
        rate estimate for each strain and biological replicate.
    cols : dict, keys 'time', and 'od'
        The column names of the time and optical density measurements in the 
        DataFrame. 
    return_opts : bool
        If `True`, the approximated covariance matrix, optimal parameters, and 
        approximate hessian matrix for each grouping is returned as a dictionary.
    print_params : bool 
        If `True`, the estimated parameters will be printed to the screen when 
        the estimation is finished. 
    
    Returns
    -------
    data_df : pandas DataFrame
        A pandas DataFrame with the converted time measurements cropped to the 
        provided optical density bounds.
    param_df : pandas DataFrame
        A pandas DataFrame containing the parameters, values, and their 95% credible 
        intervals for each obejct in the provided groupby.
    opts : dict 
        If `return_opts = True`, the estimated covariance matrix, optimal parameters, 
        and approximate Hessian matrix is returned. 

    Notes
    -----
    This function infers the "maximal a posteriori" parameter set using a 
    Bayesian definition of probability. This function  calls the posterior 
    defined by `cremerlab.bayes.steady_state_log_posterior` which contains 
    more information about the construction of the statistical model.
    """

    # TODO: Include type checks

    if (groupby is not None) & (type(groupby) is not list):
        groupby = [groupby]

    # Unpack the time col
    time_col = cols['time'] 
    od_col = cols['od']

    # Determine the OD bounds to consider
    if od_bounds is not None:
        data = data[(data[od_col] >= od_bounds[0]) & (data[od_col] <= od_bounds[1])]

    faux_groupby = False
    if groupby is None:            
        faux_groupby = True
        data['__group_idx'] = 1   
        groupby=['__group_idx']
        iterator = data.groupby(groupby)    
    else:
        iterator = tqdm.tqdm(data.groupby(groupby), desc='Estimating parameters...')

    # Iterate through each grouping
    data_dfs = []
    param_dfs = []
    opts = {'groupby':groupby}
    iter = 0  # Iterator for opts
    output = """\n
============================================================
Parameter Estimate Summary
============================================================
\n
"""

    for g, d in iterator:
        
        # Convert time if necessary
        if convert_time:
            d[time_col] = pd.to_datetime(d[time_col])
            d.sort_values(by=time_col, inplace=True)
            d['elapsed_time_hr'] = d[time_col].values - d[time_col].values[0]
            d['elapsed_time_hr'] = (d['elapsed_time_hr'].astype('timedelta64[m]')
                                    )/60
            _time = d['elapsed_time_hr'].values
            _od = d[od_col].values
        else:
            _time  = d[time_col].values
            _od = d[od_col].values

        # Define the arguments and initial guesses of the parameters
        # lam_guess = np.mean(np.diff(np.log(_od)) / np.diff(_time))
        params = [1, _od.min(), 0.1]
        args = (_time, _od)
    
        # Compute the MAP
        res = scipy.optimize.minimize(steady_state_growth_rate_log_posterior,
                                      params, args=args, method="powell")
        # Get the optimal parameters
        popt = res.x

        # Compute the Hessian and covariance matrix
        hes = smnd.approx_hess(popt, steady_state_growth_rate_log_posterior, args=args)
        cov = np.linalg.inv(hes)

        # Extract the MAP parameters and CI
        lam_MAP, od_init_MAP, sigma_MAP = popt
        lam_CI =  1.96 * np.sqrt(cov[0, 0]) 
        od_init_CI =  1.96 * np.sqrt(cov[1, 1]) 
        sigma_CI =  1.96 * np.sqrt(cov[2, 2]) 

        if print_params:
            if faux_groupby == False: 
                header =  f"""Parameter Estimates for grouping {groupby}: {g}
------------------------------------------------------------
""" 
            else:
                header =  """Parameter Estimates
------------------------------------------------------------
"""
            output +=  header + f"""growth rate,  λ = {lam_MAP:0.2f} ± {lam_CI:0.3f} [per unit time]
initial OD, OD_0 = {od_init_MAP:0.2f} ± {lam_CI:0.3f} [a.u.]
homoscedastic error, σ = {sigma_MAP:0.2f} ± {sigma_CI:0.3f} [a.u.]
\n
"""
        # Assemble the data dataframe
        _data_df = pd.DataFrame([])
        if convert_time:
            _data_df['elapsed_time_hr'] = _time
        else:
            _data_df[time_col] = _time
        _data_df[od_col] = _od

        # Include other columns that were not carried through
        colnames = [k for k in d.keys() if k not in [time_col, od_col]]
        for c in colnames:
            _data_df[c] = d[c].values
        if '__group_idx' in _data_df.keys():
            _data_df.drop(columns=['__group_idx'], inplace=True)
        _data_df.rename(columns={'od':od_col})
    
        # Assemble the parameter dataframe
        _param_df = pd.DataFrame([])
        for title, MAP, CI in zip(['growth_rate', 'od_init', 'sigma'],
                                  [lam_MAP, od_init_MAP, sigma_MAP],
                                  [lam_CI, od_init_CI, sigma_CI]): 
            _param_df = _param_df.append({'parameter':title,
                                          'map_val': MAP,
                                          'cred_int': CI},
                                          ignore_index=True)

        # Add grouping identifier if provided
        if groupby is not None:
            # if type(g) is not list:
                # _g = [g]
            _g = g
            for title, value in zip(groupby, _g):
                _data_df[title] = value 
                _param_df[title] = value
        
        #  Append the dataframes to the storage lists
        param_dfs.append(_param_df)
        data_dfs.append(_data_df)
        opts[iter] = {'groupby': g, 'cov':cov, 'popt':popt, 'hessian':hes}  
        iter += 1

    # Concatenate the dataframes for return
    if len(data_dfs) == 1:
        data_df = data_dfs[0]
        param_df = param_dfs[0]
    else:
        data_df = pd.concat(data_dfs, sort=False)
        param_df = pd.concat(param_dfs, sort=False)
    
    if print_params:
        print(output)

    if return_opts:
        return_obj = [data_df, param_df, opts]
    else:
        return_obj = [data_df, param_df]

    return return_obj