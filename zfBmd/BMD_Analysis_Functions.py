#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Paritosh Pande
Pacific Northwest National Lab, Richland, WA
Original created on: April 2020
"""

# Load modules
import numpy as np
from scipy import stats
from statsmodels.base.model import GenericLikelihoodModel

import warnings
warnings.filterwarnings('ignore')

BMR = 0.1
BMR_50 = 0.5

########################
## LOGISTIC FUNCTIONS ##
########################

def logistic_fun(dose, params):
    alpha_ = params[0].astype('float')
    beta_ = params[1].astype('float')
    dose = dose.astype('float')
    prob_dose = 1/(1 + np.exp(-alpha_ - beta_*dose))
    return prob_dose

class Logistic(GenericLikelihoodModel):
    def __init__(self, endog, exog=None, **kwds):
        super(Logistic, self).__init__(endog, exog, **kwds)


    def nloglikeobs(self, params):
        alpha_ = params[0]
        beta_ = params[1]
        dose = self.endog[:,0].flatten()

        params = [alpha_, beta_]
        probs = logistic_fun(dose, params)
        num_affected = self.endog[:,1].flatten()
        num_total = self.endog[:,2].flatten()

        log_lhood = (num_affected * np.log(probs)) \
                  + ((num_total - num_affected) * (np.log(1 - probs)))
        return -log_lhood

    def fit(self, start_params = None, maxiter = 10000, maxfun = 5000, **kwds):
        if start_params is None:
            mu_0 = self.endog[:,0].flatten().mean()
            s_0  = np.sqrt(3)*np.std(self.endog[:,0].flatten())/np.pi
            alpha_0 = -mu_0/s_0
            beta_0 = 1/s_0
            start_params = np.array([alpha_0, beta_0])

        return super(Logistic, self).fit(start_params = start_params, maxiter = maxiter, maxfun = maxfun, method = 'lbfgs', bounds = [[None, None],[1e-5,None]], disp=0, **kwds)


class Logistic_BMD(GenericLikelihoodModel):
    def __init__(self, endog, exog=None, **kwds):
        super(Logistic_BMD, self).__init__(endog, exog=None, **kwds)

    def nloglikeobs(self, params):

        alpha_ = params[0]
        bmdl_ = params[1]

        p_0 = 1/(1 + np.exp(-alpha_))

        chi_ = (1 - p_0) * BMR + p_0
        xi_ = np.log((1 - chi_)/chi_)
        beta_reparam = -(alpha_ + xi_)/bmdl_
        dose = self.endog[:,0].flatten()
        num_affected = self.endog[:,1].flatten()
        num_total = self.endog[:,2].flatten()

        probs = logistic_fun(dose, [alpha_, beta_reparam])
        log_lhood = (num_affected * np.log(probs)) \
                  + ((num_total - num_affected) * (np.log(1 - probs)))
        return -log_lhood

    def profile_ll_fit(self, start_params = None, maxiter = 10000, maxfun = 5000, **kwds):
        return super(Logistic_BMD, self).fit(start_params = start_params, maxiter = maxiter, maxfun = maxfun, method = 'lbfgs', bounds = [[None,None],[start_params[1],start_params[1]]],  disp = 0, **kwds)


# Gamma function

def gamma_fun(dose, params):
    g_ = params[0].astype('float')
    alpha_ = params[1].astype('float')
    beta_ = params[2].astype('float')
    dose = dose.astype('float')
    prob_dose = g_ + (1 - g_) * stats.gamma.cdf(dose, a = alpha_, scale = 1/beta_)
    return prob_dose

class Gamma(GenericLikelihoodModel):
    def __init__(self, endog, exog=None, **kwds):
        super(Gamma, self).__init__(endog, exog, **kwds)

    def nloglikeobs(self, params):
        g_ = params[0].astype('float')
        alpha_ = params[1].astype('float')
        beta_ = params[2].astype('float')
        dose = self.endog[:,0].flatten()
        params = [g_, alpha_, beta_]
        probs = gamma_fun(dose, params)

        num_affected = self.endog[:,1].flatten()
        num_total = self.endog[:,2].flatten()

        log_lhood = (num_affected * np.log(probs)) \
                  + ((num_total - num_affected) * (np.log(1 - probs)))
        return -log_lhood

    def fit(self, start_params = None, maxiter = 10000, maxfun = 5000, **kwds):
        if start_params is None:
            g_0 = 0.1
            beta_0 = self.endog[:,0].flatten().mean()/self.endog[:,0].flatten().var()
            alpha_0 = self.endog[:,0].flatten().mean() * beta_0
            start_params = np.array([g_0, alpha_0, beta_0])

        return super(Gamma, self).fit(start_params = start_params, maxiter = maxiter, maxfun = maxfun, method = 'lbfgs', bounds = [[1e-5,0.99],[0.2, 18],[1e-5, None]],  disp = 0, **kwds)

class Gamma_BMD(GenericLikelihoodModel):
    def __init__(self, endog, exog=None, **kwds):
        super(Gamma_BMD, self).__init__(endog, exog=None, **kwds)

    def nloglikeobs(self, params):

        g_ = params[0].astype('float')
        alpha_ = params[1].astype('float')
        bmdl_ = params[2].astype('float')

        beta_reparam = stats.gamma.ppf(BMR, alpha_)/bmdl_

        dose = self.endog[:,0].flatten()
        num_affected = self.endog[:,1].flatten()
        num_total = self.endog[:,2].flatten()

        probs = gamma_fun(dose, [g_, alpha_, beta_reparam])

        log_lhood = (num_affected * np.log(probs)) \
                  + ((num_total - num_affected) * (np.log(1 - probs)))
        return -log_lhood

    def profile_ll_fit(self, start_params = None, maxiter = 10000, maxfun = 5000, **kwds):
        return super(Gamma_BMD, self).fit(start_params = start_params, maxiter = maxiter, maxfun = maxfun, method = 'lbfgs', bounds = [[1e-5,0.99],[0.2, 18],[start_params[2],start_params[2]]], disp = 0, **kwds)

# Weibull function

def weibull_fun(dose, params):
    g_ = params[0].astype('float')
    alpha_ = params[1].astype('float')
    beta_ = params[2].astype('float')
    dose = dose.astype('float')
    prob_dose = g_ + (1 - g_) * (1 - np.exp(-beta_ * (dose.astype('float') ** alpha_)))
    return prob_dose

class Weibull(GenericLikelihoodModel):
    def __init__(self, endog, exog=None, **kwds):
        super(Weibull, self).__init__(endog, exog, **kwds)

    def nloglikeobs(self, params):
        g_ = params[0].astype('float')
        alpha_ = params[1].astype('float')
        beta_ = params[2].astype('float')
        dose = self.endog[:,0].flatten()
        params = [g_, alpha_, beta_]
        probs = weibull_fun(dose, params)

        num_affected = self.endog[:,1].flatten()
        num_total = self.endog[:,2].flatten()

        log_lhood = (num_affected * np.log(probs)) \
                  + ((num_total - num_affected) * (np.log(1 - probs)))
        return -log_lhood

    def fit(self, start_params = None, maxiter = 10000, maxfun = 5000, **kwds):
        if start_params is None:
            g_0 = 0.1

            dose = self.endog[:,0].flatten()
            num_affected = self.endog[:,1].flatten()
            num_total = self.endog[:,2].flatten()
            frac_affected = num_affected/num_total

            X = np.append(np.ones([len(dose[1:]),1]),np.log(np.reshape(dose[1:],(len(dose[1:]),1))),1)
            Y = np.array(np.reshape(np.log(-np.log(1 - frac_affected[1:])),(len(dose[1:]),1)))

            betas = np.linalg.inv((X.T).dot(X)).dot(X.T).dot(Y)

            alpha_0 = betas[1]
            beta_0 = np.exp(betas[0])

            start_params = np.array([g_0, alpha_0, beta_0])

        return super(Weibull, self).fit(start_params = start_params, maxiter = maxiter, maxfun = maxfun, method = 'lbfgs', bounds = [[1e-5,0.99],[1e-5,None],[1e-9,None]], disp = 0, **kwds)

class Weibull_BMD(GenericLikelihoodModel):
    def __init__(self, endog, exog=None, **kwds):
        super(Weibull_BMD, self).__init__(endog, exog=None, **kwds)

    def nloglikeobs(self, params):
        g_ = params[0].astype('float')
        alpha_ = params[1].astype('float')
        bmdl_ = params[2].astype('float')

        beta_reparam = -np.log(1-BMR)/(bmdl_**alpha_)

        dose = self.endog[:,0].flatten()
        num_affected = self.endog[:,1].flatten()
        num_total = self.endog[:,2].flatten()

        probs = weibull_fun(dose, [g_, alpha_, beta_reparam])

        log_lhood = (num_affected * np.log(probs)) \
                  + ((num_total - num_affected) * (np.log(1 - probs)))
        return -log_lhood

    def profile_ll_fit(self, start_params = None, maxiter = 10000, maxfun = 5000, **kwds):
        return super(Weibull_BMD, self).fit(start_params = start_params, maxiter = maxiter, maxfun = maxfun, method = 'lbfgs', bounds = [[1e-5,0.99],[1e-5,None],[start_params[2],start_params[2]]], disp = 0,**kwds)

# Log-logistic function

def log_logistic_fun(dose, params):
    g_ = params[0].astype('float')
    alpha_ = params[1].astype('float')
    beta_ = params[2].astype('float')
    dose = dose.astype('float')
    dose_nonzero = dose.copy()
    dose_nonzero[dose_nonzero == 0] = 1e-9
    prob_dose = g_ + (1 - g_)/(1 + np.exp(-alpha_ - beta_*np.log(dose_nonzero.astype('float'))))
    return prob_dose

class Log_Logistic(GenericLikelihoodModel):
    def __init__(self, endog, exog=None, **kwds):
        super(Log_Logistic, self).__init__(endog, exog, **kwds)

    def nloglikeobs(self, params):
        g_ = params[0].astype('float')
        alpha_ = params[1].astype('float')
        beta_ = params[2].astype('float')
        dose = self.endog[:,0].flatten()

        params = [g_,alpha_, beta_]
        probs = log_logistic_fun(dose, params)
        num_affected = self.endog[:,1].flatten()
        num_total = self.endog[:,2].flatten()

        log_lhood = (num_affected * np.log(probs)) \
                  + ((num_total - num_affected) * (np.log(1 - probs)))
        return -log_lhood

    def fit(self, start_params = None, maxiter = 10000, maxfun = 5000, **kwds):
        if start_params is None:
            doses = self.endog[:,0].copy().flatten()
            nonzero_doses = doses[1:]
            g_0 = 0.1
            mu_0 = np.log(nonzero_doses).mean()
            s_0  = np.sqrt(3)*np.std(np.log(nonzero_doses))/np.pi
            alpha_0 = -mu_0/s_0
            beta_0 = 1/s_0
            start_params = np.array([g_0, alpha_0, beta_0])

        return super(Log_Logistic, self).fit(start_params = start_params, maxiter = maxiter, maxfun = maxfun, method = 'lbfgs', bounds = [[1e-5,0.99],[None, None],[None, None]],  disp = 0, **kwds)

class Log_Logistic_BMD(GenericLikelihoodModel):
    def __init__(self, endog, exog=None, **kwds):
        super(Log_Logistic_BMD, self).__init__(endog, exog=None, **kwds)

    def nloglikeobs(self, params):

        g_ = params[0].astype('float')
        beta_ = params[1].astype('float')
        bmdl_ = params[2].astype('float')

        alpha_reparam = np.log(BMR/(1-BMR)) - beta_*np.log(bmdl_)

        dose = self.endog[:,0].flatten()
        num_affected = self.endog[:,1].flatten()
        num_total = self.endog[:,2].flatten()

        probs = log_logistic_fun(dose, [g_, alpha_reparam, beta_])

        log_lhood = (num_affected * np.log(probs)) \
                  + ((num_total - num_affected) * (np.log(1 - probs)))
        return -log_lhood

    def profile_ll_fit(self, start_params = None, maxiter = 10000, maxfun = 5000, **kwds):
        return super(Log_Logistic_BMD, self).fit(start_params = start_params, maxiter = maxiter, maxfun = maxfun, method = 'lbfgs', bounds = [[1e-5,0.99],[start_params[1]/2,start_params[1]*2],[start_params[2],start_params[2]]], disp = 0,**kwds)

# Probit function

def probit_fun(dose, params):
    alpha_ = params[0].astype('float')
    beta_ = params[1].astype('float')
    dose = dose.astype('float')
    prob_dose = stats.norm.cdf((alpha_ + beta_ * dose), loc=0, scale=1)
    return prob_dose

class Probit(GenericLikelihoodModel):
    def __init__(self, endog, exog=None, **kwds):
        super(Probit, self).__init__(endog, exog, **kwds)

    def nloglikeobs(self, params):
        alpha_ = params[0].astype('float')
        beta_ = params[1].astype('float')
        dose = self.endog[:,0].flatten()

        probs = probit_fun(dose, [alpha_, beta_])

        num_affected = self.endog[:,1].flatten()
        num_total = self.endog[:,2].flatten()

        log_lhood = (num_affected * np.log(probs)) \
                  + ((num_total - num_affected) * (np.log(1 - probs)))
        return -log_lhood

    def fit(self, start_params = None, maxiter = 10000, maxfun = 5000, **kwds):
        if start_params is None:
            dose = self.endog[:,0].flatten()
            num_affected = self.endog[:,1].flatten()
            num_total = self.endog[:,2].flatten()
            alpha_0 = stats.norm.ppf(num_affected[0]/num_total[0])
            beta_0 = (stats.norm.ppf(num_affected[-1]/num_total[-1]) - alpha_0)/dose[-1]
            start_params = np.array([alpha_0, beta_0])

        return super(Probit, self).fit(start_params = start_params, maxiter = maxiter, maxfun = maxfun, method = 'lbfgs', bounds = [[None, None],[1e-5,None]],  disp = 0, **kwds)

class Probit_BMD(GenericLikelihoodModel):
    def __init__(self, endog, exog=None, **kwds):
        super(Probit_BMD, self).__init__(endog, exog=None, **kwds)

    def nloglikeobs(self, params):

        alpha_ = params[0].astype('float')
        bmdl_ = params[1].astype('float')

        p_0 = stats.norm.cdf(alpha_)
        xi_ = (1 - p_0) * BMR + p_0

        beta_reparam = (stats.norm.ppf(xi_) - alpha_)/bmdl_

        dose = self.endog[:,0].flatten()
        num_affected = self.endog[:,1].flatten()
        num_total = self.endog[:,2].flatten()

        probs = probit_fun(dose, [alpha_, beta_reparam])

        log_lhood = (num_affected * np.log(probs)) \
                  + ((num_total - num_affected) * (np.log(1 - probs)))
        return -log_lhood

    def profile_ll_fit(self, start_params = None, maxiter = 10000, maxfun = 5000, **kwds):
        return super(Probit_BMD, self).fit(start_params = start_params, maxiter = maxiter, maxfun = maxfun, method = 'lbfgs', bounds = [[None,None],[start_params[1],start_params[1]]], disp = 0,**kwds)

# Log-probit function

def log_probit_fun(dose, params):
    g_ = params[0].astype('float')
    alpha_ = params[1].astype('float')
    beta_ = params[2].astype('float')
    dose_nonzero = dose.copy().astype('float')
    dose_nonzero[dose_nonzero == 0] = 1e-9
    prob_dose = g_ + (1 - g_) * stats.norm.cdf((alpha_ + beta_ * np.log(dose_nonzero)), loc=0, scale=1)
    return prob_dose

class Log_Probit(GenericLikelihoodModel):
    def __init__(self, endog, exog=None, **kwds):
        super(Log_Probit, self).__init__(endog, exog, **kwds)

    def nloglikeobs(self, params):
        g_ = params[0]
        alpha_ = params[1]
        beta_ = params[2]
        dose = self.endog[:,0].flatten()

        probs = log_probit_fun(dose, [g_, alpha_, beta_])
        num_affected = self.endog[:,1].flatten()
        num_total = self.endog[:,2].flatten()

        log_lhood = (num_affected * np.log(probs)) \
                  + ((num_total - num_affected) * (np.log(1 - probs)))
        return -log_lhood

    def fit(self, start_params = None, maxiter = 10000, maxfun = 5000, **kwds):
        if start_params is None:
            g_0 = 0.1

            dose = self.endog[:,0].flatten()
            num_affected = self.endog[:,1].flatten()
            num_total = self.endog[:,2].flatten()
            frac_affected = num_affected/num_total
            X = np.array([[1, dose[1]],[1, dose[-1]]])
            Y = (np.array([stats.norm.ppf(1 - frac_affected[1]), \
                           stats.norm.ppf(1 - frac_affected[-1])])).T

            betas = np.linalg.inv((X.T).dot(X)).dot(X.T).dot(Y)

            alpha_0 = betas[0]
            beta_0 = max(1e-5,betas[1])

            dose = self.endog[:,0].flatten()
            start_params = np.array([g_0, alpha_0, beta_0])

        return super(Log_Probit, self).fit(start_params = start_params, maxiter = maxiter, maxfun = maxfun, method = 'lbfgs', bounds = [[1e-5,0.99],[None,None],[1e-5,None]],  disp = 0, **kwds)

class Log_Probit_BMD(GenericLikelihoodModel):
    def __init__(self, endog, exog=None, **kwds):
        super(Log_Probit_BMD, self).__init__(endog, exog=None, **kwds)

    def nloglikeobs(self, params):

        g_ = params[0].astype('float')
        alpha_ = params[1].astype('float')
        bmdl_ = params[2].astype('float')

        beta_reparam = (stats.norm.ppf(BMR) - alpha_)/np.log(bmdl_)

        dose = self.endog[:,0].flatten()
        num_affected = self.endog[:,1].flatten()
        num_total = self.endog[:,2].flatten()

        probs = log_probit_fun(dose, [g_, alpha_, beta_reparam])

        log_lhood = (num_affected * np.log(probs)) \
                  + ((num_total - num_affected) * (np.log(1 - probs)))
        return -log_lhood

    def profile_ll_fit(self, start_params = None, maxiter = 10000, maxfun = 5000, **kwds):
        return super(Log_Probit_BMD, self).fit(start_params = start_params, maxiter = maxiter, maxfun = maxfun, method = 'lbfgs', bounds = [[1e-9,0.99],[None,None],[start_params[2],start_params[2]]],  disp = 0, **kwds)

# Multistage(degree 2)

def multistage_2_fun(dose, params):
    g_ = params[0]
    beta1_ = params[1]
    beta2_ = params[2]
    prob_dose = g_ + (1 - g_) * (1 - np.exp(-(beta1_ * dose) \
                                            -(beta2_ * (dose ** 2))))
    return prob_dose

class Multistage_2(GenericLikelihoodModel):
    def __init__(self, endog, exog=None, **kwds):
        super(Multistage_2, self).__init__(endog, exog, **kwds)

    def nloglikeobs(self, params):
        g_ = params[0].astype('float')
        beta1_ = params[1].astype('float')
        beta2_ = params[2].astype('float')
        dose = self.endog[:,0].flatten().astype('float')

        probs = multistage_2_fun(dose, [g_, beta1_, beta2_])
        num_affected = self.endog[:,1].flatten()
        num_total = self.endog[:,2].flatten()

        log_lhood = (num_affected * np.log(probs)) \
                  + ((num_total - num_affected) * (np.log(1 - probs)))
        return -log_lhood

    def fit(self, start_params = None, maxiter = 10000, maxfun = 5000, **kwds):
        if start_params is None:
            g_0 = 0.05

            dose = self.endog[:,0].flatten()
            num_affected = self.endog[:,1].flatten()
            num_total = self.endog[:,2].flatten()
            frac_affected = num_affected/num_total

            X = np.append(np.reshape(dose[1:],(len(dose[1:]),1)),(np.reshape(dose[1:],(len(dose[1:]),1)))**2,1)
            Y = np.array(-np.log(1 - frac_affected[1:]))

            betas = np.linalg.inv((X.T).dot(X)).dot(X.T).dot(Y)

            beta1_0 = betas[0]
            beta2_0 = betas[1]

            start_params = np.array([g_0, beta1_0, beta2_0])

        return super(Multistage_2, self).fit(start_params = start_params, maxiter = maxiter, maxfun = maxfun, method = 'lbfgs', bounds = [[1e-9,0.99],[1e-9,None],[1e-9,None]],  disp = 0, **kwds)

class Multistage_2_BMD(GenericLikelihoodModel):
    def __init__(self, endog, exog=None, **kwds):
        super(Multistage_2_BMD, self).__init__(endog, exog=None, **kwds)

    def nloglikeobs(self, params):

        g_ = params[0].astype('float')
        beta1_ = params[1].astype('float')
        bmdl_ = params[2].astype('float')

        beta2_reparam = -(np.log(1-BMR) + beta1_*bmdl_)/(bmdl_**2)

        dose = self.endog[:,0].flatten()
        num_affected = self.endog[:,1].flatten()
        num_total = self.endog[:,2].flatten()

        probs = multistage_2_fun(dose, [g_, beta1_, beta2_reparam])

        log_lhood = (num_affected * np.log(probs)) \
                  + ((num_total - num_affected) * (np.log(1 - probs)))
        return -log_lhood

    def profile_ll_fit(self, start_params = None, maxiter = 10000, maxfun = 5000, **kwds):
        return super(Multistage_2_BMD, self).fit(start_params = start_params, maxiter = maxiter, maxfun = maxfun, method = 'lbfgs', bounds = [[1e-9,0.99],[1e-9,None],[start_params[2],start_params[2]]], disp = 0, **kwds)

# Quantal-linear model

def quantal_linear_fun(dose, params):
    g_ = params[0].astype('float')
    beta_ = params[1].astype('float')
    dose = dose.astype('float')
    prob_dose = g_ + (1 - g_) * (1 - np.exp(-(beta_ * dose)))
    return prob_dose

class Quantal_Linear(GenericLikelihoodModel):
    def __init__(self, endog, exog=None, **kwds):
        super(Quantal_Linear, self).__init__(endog, exog, **kwds)

    def nloglikeobs(self, params):
        g_ = params[0].astype('float')
        beta_ = params[1].astype('float')

        dose = self.endog[:,0].flatten()

        probs = quantal_linear_fun(dose, [g_, beta_])
        num_affected = self.endog[:,1].flatten()
        num_total = self.endog[:,2].flatten()

        log_lhood = (num_affected * np.log(probs)) \
                  + ((num_total - num_affected) * (np.log(1 - probs)))
        return -log_lhood

    def fit(self, start_params = None, maxiter = 10000, maxfun = 5000, **kwds):
        if start_params is None:
            g_0 = 0.1
            beta_0 = 1/((self.endog[:,0].flatten().mean()))/np.log(2)
            start_params = np.array([g_0, beta_0])

        return super(Quantal_Linear, self).fit(start_params = start_params, maxiter = maxiter, maxfun = maxfun, method = 'lbfgs', bounds = [[1e-5,0.99],[1e-5,None]],  disp = 0, **kwds)

class Quantal_Linear_BMD(GenericLikelihoodModel):
    def __init__(self, endog, exog=None, **kwds):
        super(Quantal_Linear_BMD, self).__init__(endog, exog=None, **kwds)

    def nloglikeobs(self, params):

        g_ = params[0].astype('float')
        bmdl_ = params[1].astype('float')
        beta_reparam = -np.log(1-BMR)/bmdl_

        dose = self.endog[:,0].flatten()
        num_affected = self.endog[:,1].flatten()
        num_total = self.endog[:,2].flatten()

        probs = quantal_linear_fun(dose, [g_, beta_reparam])

        log_lhood = (num_affected * np.log(probs)) \
                  + ((num_total - num_affected) * (np.log(1 - probs)))
        return -log_lhood

    def profile_ll_fit(self, start_params = None, maxiter = 10000, maxfun = 5000, **kwds):
        return super(Quantal_Linear_BMD, self).fit(start_params = start_params, maxiter = maxiter, maxfun = maxfun, method = 'lbfgs', bounds = [[1e-5,0.99],[start_params[1],start_params[1]]], disp = 0, **kwds)
