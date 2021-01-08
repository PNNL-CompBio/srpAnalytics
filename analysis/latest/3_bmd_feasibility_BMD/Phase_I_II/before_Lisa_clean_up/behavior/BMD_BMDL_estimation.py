#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Paritosh Pande
Pacific Northwest National Lab, Richland, WA
Original created on: April 2020
Modified: June 2020, Added BMD50 calculation
"""
import pandas as pd
from scipy import stats
import numpy as np
import BMD_Analysis_Functions as baf

BMR = "global"
BMR = 0.1

BMR_50 = "global"
BMR_50 = 0.5

# Maximum iteration limit for BMDL convergence
bmdl_max_iter = 1000
tol_thresh = 1e-4

def analyze_dose_response_data(test_dose_response, model_names=None, bmdl_analysis_flag=1):
    # Setup the dataframe for data reporting
    model_predictions = pd.DataFrame(columns = ['Model', 'Chi-squared', \
                                                'p-val', 'AIC', 'BMD10', \
                                                'BMDL10', 'BMD50', 'Scaled Residuals', \
                                                'Optimized Params', 'Conv Flag', 'BMDL Conv'])
    if model_names is None:
        model_names = ['logistic', 'gamma', 'weibull', 'log_logistic', \
                   'probit', 'log_probit', 'multistage_2', 'quantal_linear']
    
   
    # Fit different models
    for model_index, model_name in enumerate(model_names):
        
        print(model_name)
        # Set flags for model and bmdl_convergence convergence
        model_converge_flag = 1
        bmdl_converge_flag = 1
        bmdl_mid_val = np.nan
        
        # Setup model
        if(model_name == 'logistic'):
            model = baf.Logistic(test_dose_response[['dose','num_affected', 'total_num']].copy())
            res = model.fit()
            
            print('Model Convergence: ' + str(res.mle_retvals['converged']))
            
            # Check model convergence
            if(res.mle_retvals['converged'] is True):
                model_converge_flag = 0
                
                pred_vals = baf.logistic_fun(test_dose_response.dose, res.params)
    
                alpha_chi_square = 0.9
                bmdl_llv_thresh = res.llf - stats.chi2.ppf(alpha_chi_square, 1)/2
                
                # Estimate BMD 
                alpha_ = res.params[0]
                beta_ = res.params[1]
                
                bmd = np.log((1 + np.exp(-alpha_)*BMR)/(1-BMR))/beta_
                bmd_50 = np.log((1 + np.exp(-alpha_)*BMR_50)/(1-BMR))/beta_
                
                print('Fit Paramas:' + str([alpha_,beta_]))
                print('BMD 10: ' + str(bmd))
                print('BMD 50: ' + str(bmd_50))
            else:
                pred_vals = np.nan
                bmd = np.nan
                bmd_50 = np.nan
            
            if((bmdl_analysis_flag) and (res.mle_retvals['converged'] is True)):
                print('Estimating BMDL ...')
                # Using a bisection method for finding the BMDL value
                bmdl_val_lo = bmd/10
                bmdl_val_hi = bmd

                tol = 1
                bmdl_iter_count = 0
                while ((tol > tol_thresh) and (bmdl_iter_count < bmdl_max_iter)):
                    bmdl_iter_count+=1
                    # If maximum iterations are reached, set convergence flag to False and quit
                    if(bmdl_iter_count == bmdl_max_iter):
                        bmdl_converge_flag = 1
                        bmdl_mid_val = np.nan
                        break
                        
                    bmdl_mid_val = (bmdl_val_lo + bmdl_val_hi)/2
           
                    # Examine the sign of (llv at bmdl_mid_val) - bmdl_llv_thresh
                    # and update bmdl_lo and hi values appropriately
                    model_logistic_bmd = baf.Logistic_BMD(test_dose_response[['dose','num_affected', 'total_num']].copy())
                    logistic_bmr_fit = model_logistic_bmd.profile_ll_fit([alpha_, bmdl_mid_val])
                    
                    # Check model convergence
                    if(logistic_bmr_fit.mle_retvals['converged'] is True):
                        bmdl_converge_flag = 0
                    else:
                        bmdl_converge_flag = 1
                        break
                    
                    bmdl_llv = logistic_bmr_fit.llf

                    if((bmdl_llv - bmdl_llv_thresh)>0):
                        bmdl_val_hi = bmdl_mid_val
                    else:
                        bmdl_val_lo = bmdl_mid_val
                
                    tol = abs(bmdl_llv - bmdl_llv_thresh)
                    
                print('BMDL Convergence: ' + str(not(bool(bmdl_converge_flag))))
            
        elif(model_name == 'gamma'):
            model = baf.Gamma(test_dose_response[['dose','num_affected', 'total_num']].copy())
            res = model.fit()

            print('Model Convergence: ' + str(res.mle_retvals['converged']))
            
            # Check model convergence
            if(res.mle_retvals['converged'] is True):
                model_converge_flag = 0            
            
                pred_vals = baf.gamma_fun(test_dose_response.dose, res.params)
     
                alpha_chi_square = 0.9
                bmdl_llv_thresh = res.llf - stats.chi2.ppf(alpha_chi_square, 1)/2
            
                g_ = res.params[0]
                alpha_ = res.params[1]
                beta_ = res.params[2]
            
                # Estimate BMD
                bmd = stats.gamma.ppf(BMR, alpha_, scale = 1/beta_)
                bmd_50 = stats.gamma.ppf(BMR_50, alpha_, scale = 1/beta_)
                print('Fit Paramas:' + str([g_,alpha_,beta_]))
                print('BMD 10: ' + str(bmd))
                print('BMD 50: ' + str(bmd_50))
            else:
                pred_vals = np.nan
                bmd = np.nan
                bmd_50 = np.nan
            
            if((bmdl_analysis_flag) and (res.mle_retvals['converged'] is True)):
                print('Estimating BMDL ...')
                # Using a bisection method for finding the BMDL value
                bmdl_val_lo = bmd/10
                bmdl_val_hi = bmd

                tol = 1
                bmdl_iter_count = 0
                while ((tol > tol_thresh) and (bmdl_iter_count < bmdl_max_iter)):
                    bmdl_iter_count+=1
                    # If maximum iterations are reached, set convergence flag to False and quit
                    if(bmdl_iter_count == bmdl_max_iter):
                        bmdl_converge_flag = 1
                        bmdl_mid_val = np.nan
                        break                    
                
                    bmdl_mid_val = (bmdl_val_lo + bmdl_val_hi)/2
                               
                    # Examine the sign of (llv at bmdl_mid_val) - bmdl_llv_thresh
                    # and update bmdl_lo and hi values appropriately
                    model_gamma_bmd = baf.Gamma_BMD(test_dose_response[['dose','num_affected', 'total_num']].copy())
                    gamma_bmr_fit = model_gamma_bmd.profile_ll_fit([g_, alpha_, bmdl_mid_val])

                    # Check model convergence
                    if(gamma_bmr_fit.mle_retvals['converged'] is True):
                        bmdl_converge_flag = 0
                    else:
                        bmdl_converge_flag = 1
                        break
                    
                    bmdl_llv = gamma_bmr_fit.llf

                    if((bmdl_llv - bmdl_llv_thresh)>0):
                        bmdl_val_hi = bmdl_mid_val
                    else:
                        bmdl_val_lo = bmdl_mid_val
                
                    tol = abs(bmdl_llv - bmdl_llv_thresh)                    
    
                print('BMDL Convergence: ' + str(not(bool(bmdl_converge_flag))))
                           
        elif(model_name == 'weibull'):
            model = baf.Weibull(test_dose_response[['dose','num_affected', 'total_num']].copy())
            res = model.fit()

            print('Model Convergence: ' + str(res.mle_retvals['converged']))            
            # Check model convergence
            if(res.mle_retvals['converged'] is True):
                model_converge_flag = 0            
            
                pred_vals = baf.weibull_fun(test_dose_response.dose, res.params)
    
                alpha_chi_square = 0.9
                bmdl_llv_thresh = res.llf - stats.chi2.ppf(alpha_chi_square, 1)/2
                
                g_ = res.params[0]
                alpha_ = res.params[1]
                beta_ = res.params[2]
            
                # Estimate BMD
                bmd = (-np.log(1 - BMR)/beta_)**(1/alpha_)
                bmd_50 = (-np.log(1 - BMR_50)/beta_)**(1/alpha_)
                print('Fit Paramas:' + str([g_,alpha_,beta_]))
                print('BMD 10: ' + str(bmd))
                print('BMD 50: ' + str(bmd_50))
            else:
                pred_vals = np.nan
                bmd = np.nan
                bmd_50 = np.nan
                
            
            if((bmdl_analysis_flag) and (res.mle_retvals['converged'] is True)):
                print('Estimating BMDL ...')
                # Using a bisection method for finding the BMDL value
                bmdl_val_lo = bmd/10
                bmdl_val_hi = bmd

                tol = 1
                bmdl_iter_count = 0
                while ((tol > tol_thresh) and (bmdl_iter_count < bmdl_max_iter)):
                    bmdl_iter_count+=1
                    
                    # If maximum iterations are reached, set convergence flag to False and quit
                    if(bmdl_iter_count == bmdl_max_iter):
                        bmdl_converge_flag = 1
                        bmdl_mid_val = np.nan
                        break
                    
                    bmdl_mid_val = (bmdl_val_lo + bmdl_val_hi)/2
           
                    # Examine the sign of (llv at bmdl_mid_val) - bmdl_llv_thresh
                    # and update bmdl_lo and hi values appropriately
                    model_weibull_bmd = baf.Weibull_BMD(test_dose_response[['dose','num_affected', 'total_num']].copy())
                    weibull_bmr_fit = model_weibull_bmd.profile_ll_fit([g_, alpha_, bmdl_mid_val])

                    # Check model convergence
                    if(weibull_bmr_fit.mle_retvals['converged'] is True):
                        bmdl_converge_flag = 0
                    else:
                        bmdl_converge_flag = 1
                        break
                    
                    bmdl_llv = weibull_bmr_fit.llf

                    if((bmdl_llv - bmdl_llv_thresh)>0):
                        bmdl_val_hi = bmdl_mid_val
                    else:
                        bmdl_val_lo = bmdl_mid_val
                
                    tol = abs(bmdl_llv - bmdl_llv_thresh)
                print('BMDL Convergence: '+ str(not(bool(bmdl_converge_flag))))
                
        elif(model_name == 'log_logistic'):
            model = baf.Log_Logistic(test_dose_response[['dose','num_affected', 'total_num']].copy())  
            res = model.fit()

            print('Model Convergence: ' + str(res.mle_retvals['converged']))            
            # Check model convergence
            if(res.mle_retvals['converged'] is True):
                model_converge_flag = 0            
            
                pred_vals = baf.log_logistic_fun(test_dose_response.dose, res.params)
        
                alpha_chi_square = 0.9
                bmdl_llv_thresh = res.llf - stats.chi2.ppf(alpha_chi_square, 1)/2
                
                g_ = res.params[0]
                alpha_ = res.params[1]
                beta_ = res.params[2]
            
                # Estimate BMD
                bmd = np.exp((np.log(BMR/(1-BMR)) - alpha_)/beta_)
                bmd_50 = np.exp((np.log(BMR_50/(1-BMR_50)) - alpha_)/beta_)
                print('Fit Paramas:' + str([g_,alpha_,beta_]))
                print('BMD 10: ' + str(bmd))
                print('BMD 50: ' + str(bmd_50))
                
            else:
                pred_vals = np.nan
                bmd = np.nan
                bmd_50 = np.nan
            
            if((bmdl_analysis_flag) and (res.mle_retvals['converged'] is True)):
                print('Estimating BMDL ...')
                # Using a bisection method for finding the BMDL value
                bmdl_val_lo = bmd/10
                bmdl_val_hi = bmd

                tol = 1
                bmdl_iter_count = 0
                while ((tol > tol_thresh) and (bmdl_iter_count < bmdl_max_iter)):
                    bmdl_iter_count+=1
                    
                    # If maximum iterations are reached, set convergence flag to False and quit
                    if(bmdl_iter_count == bmdl_max_iter):
                        bmdl_converge_flag = 1
                        bmdl_mid_val = np.nan
                        break
                    
                    bmdl_mid_val = (bmdl_val_lo + bmdl_val_hi)/2

                    # Examine the sign of (llv at bmdl_mid_val) - bmdl_llv_thresh
                    # and update bmdl_lo and hi values appropriately
                    model_log_logistic_bmd = baf.Log_Logistic_BMD(test_dose_response[['dose','num_affected', 'total_num']].copy())
                    log_logistic_bmr_fit = model_log_logistic_bmd.profile_ll_fit([g_, beta_, bmdl_mid_val])

                    # Check model convergence
                    if(log_logistic_bmr_fit.mle_retvals['converged'] is True):
                        bmdl_converge_flag = 0
                    else: 
                        bmdl_converge_flag = 1
                        break
                    
                    bmdl_llv = log_logistic_bmr_fit.llf

                    if((bmdl_llv - bmdl_llv_thresh)>0):
                        bmdl_val_hi = bmdl_mid_val
                    else:
                        bmdl_val_lo = bmdl_mid_val

                    tol = abs(bmdl_llv - bmdl_llv_thresh)          
                print('BMDL Convergence: ' + str(not(bool(bmdl_converge_flag))))

        elif(model_name == 'probit'):
            model = baf.Probit(test_dose_response[['dose','num_affected', 'total_num']].copy())        
            res = model.fit()

            print('Model Convergence: ' + str(res.mle_retvals['converged']))            
            # Check model convergence
            if(res.mle_retvals['converged'] is True):
                model_converge_flag = 0            
            
                pred_vals = baf.probit_fun(test_dose_response.dose, res.params)
    
                alpha_chi_square = 0.9
                bmdl_llv_thresh = res.llf - stats.chi2.ppf(alpha_chi_square, 1)/2
                
                alpha_ = res.params[0]
                beta_ = res.params[1]
            
                # Estimate BMD
                p_0 = stats.norm.cdf(alpha_)
                p_BMD = p_0 + (1 - p_0)*BMR
                p_BMD_50 = p_0 + (1 - p_0)*BMR_50
                bmd = (stats.norm.ppf(p_BMD) - alpha_)/beta_
                bmd_50 = (stats.norm.ppf(p_BMD_50) - alpha_)/beta_
                print('Fit Paramas:' + str([alpha_,beta_]))
                print('BMD 10: ' + str(bmd))
                print('BMD 50: ' + str(bmd_50))
                
            else:
                pred_vals = np.nan
                bmd = np.nan
                bmd_50 = np.nan
            
            if((bmdl_analysis_flag) and (res.mle_retvals['converged'] is True)):
                print('Estimating BMDL ...')
                # Using a bisection method for finding the BMDL value
                bmdl_val_lo = bmd/10
                bmdl_val_hi = bmd

                tol = 1
                bmdl_iter_count = 0
                while ((tol > tol_thresh) and (bmdl_iter_count < bmdl_max_iter)):
                    bmdl_iter_count+=1

                    # If maximum iterations are reached, set convergence flag to False and quit
                    if(bmdl_iter_count == bmdl_max_iter):
                        bmdl_converge_flag = 1
                        bmdl_mid_val = np.nan
                        break
                    
                    bmdl_mid_val = (bmdl_val_lo + bmdl_val_hi)/2

                    # Examine the sign of (llv at bmdl_mid_val) - bmdl_llv_thresh
                    # and update bmdl_lo and hi values appropriately
                    model_probit_bmd = baf.Probit_BMD(test_dose_response[['dose','num_affected', 'total_num']].copy())
                    probit_bmr_fit = model_probit_bmd.profile_ll_fit([alpha_, bmdl_mid_val])

                    # Check model convergence
                    if(probit_bmr_fit.mle_retvals['converged'] is True):
                        bmdl_converge_flag = 0
                    else:
                        bmdl_converge_flag = 1
                        break
                    
                    bmdl_llv = probit_bmr_fit.llf

                    if((bmdl_llv - bmdl_llv_thresh)>0):
                        bmdl_val_hi = bmdl_mid_val
                    else:
                        bmdl_val_lo = bmdl_mid_val

                    tol = abs(bmdl_llv - bmdl_llv_thresh)
                print('BMDL Convergence: ' + str(not(bool(bmdl_converge_flag))))

        elif(model_name == 'log_probit'):
            model = baf.Log_Probit(test_dose_response[['dose','num_affected', 'total_num']].copy())        
            res = model.fit()

            print('Model Convergence: ' + str(res.mle_retvals['converged']))
            # Check model convergence
            if(res.mle_retvals['converged'] is True):
                model_converge_flag = 0            
            
                pred_vals = baf.log_probit_fun(test_dose_response.dose, res.params)
    
                alpha_chi_square = 0.9
                bmdl_llv_thresh = res.llf - stats.chi2.ppf(alpha_chi_square, 1)/2
                
                g_ = res.params[0]
                alpha_ = res.params[1]
                beta_ = res.params[2]
            
                # Estimate BMD
                bmd = np.exp((stats.norm.ppf(BMR) - alpha_)/beta_)
                bmd_50 = np.exp((stats.norm.ppf(BMR_50) - alpha_)/beta_)
                print('Fit Paramas:' + str([g_,alpha_,beta_]))
                print('BMD 10: ' + str(bmd))
                print('BMD 50: ' + str(bmd_50))
                
            else:
                pred_vals = np.nan
                bmd = np.nan
                bmd_50 = np.nan
                
            
            if((bmdl_analysis_flag) and (res.mle_retvals['converged'] is True)):
                print('Estimating BMDL ...')
                # Using a bisection method for finding the BMDL value
                bmdl_val_lo = bmd/10
                bmdl_val_hi = bmd

                tol = 1
                bmdl_iter_count = 0
                while ((tol > tol_thresh) and (bmdl_iter_count < bmdl_max_iter)):
                    bmdl_iter_count+=1

                    # If maximum iterations are reached, set convergence flag to False and quit
                    if(bmdl_iter_count == bmdl_max_iter):
                        bmdl_converge_flag = 1
                        bmdl_mid_val = np.nan
                        break

                    bmdl_mid_val = (bmdl_val_lo + bmdl_val_hi)/2

                    # Examine the sign of (llv at bmdl_mid_val) - bmdl_llv_thresh
                    # and update bmdl_lo and hi values appropriately
                    model_log_probit_bmd = baf.Log_Probit_BMD(test_dose_response[['dose','num_affected', 'total_num']].copy())
                    log_probit_bmr_fit = model_log_probit_bmd.profile_ll_fit([g_, alpha_, bmdl_mid_val])
                    
                    # Check model convergence
                    if(log_probit_bmr_fit.mle_retvals['converged'] is True):
                        bmdl_converge_flag = 0
                    else:
                        bmdl_converge_flag = 1
                        break
                    
                    bmdl_llv = log_probit_bmr_fit.llf

                    if((bmdl_llv - bmdl_llv_thresh)>0):
                        bmdl_val_hi = bmdl_mid_val
                    else:
                        bmdl_val_lo = bmdl_mid_val

                    tol = abs(bmdl_llv - bmdl_llv_thresh)
                print('BMDL Convergence: ' + str(not(bool(bmdl_converge_flag))))
            
        elif(model_name == 'multistage_2'):
            model = baf.Multistage_2(test_dose_response[['dose','num_affected', 'total_num']].copy())        
            res = model.fit()

            print('Model Convergence: ' + str(res.mle_retvals['converged']))            
            # Check model convergence
            if(res.mle_retvals['converged'] is True):
                model_converge_flag = 0            
            
                pred_vals = baf.multistage_2_fun(test_dose_response.dose, res.params)
                
                alpha_chi_square = 0.9
                bmdl_llv_thresh = res.llf - stats.chi2.ppf(alpha_chi_square, 1)/2
                
                g_ = res.params[0]
                beta1_ = res.params[1]
                beta2_ = res.params[2]
            
                # Estimate BMD
                bmd = (-beta1_ + np.sqrt((beta1_**2) \
                       - (4*beta2_*np.log(1 - BMR))))/(2*beta2_)
                bmd_50 = (-beta1_ + np.sqrt((beta1_**2) \
                       - (4*beta2_*np.log(1 - BMR_50))))/(2*beta2_)
                print('Fit Paramas:' + str([g_,beta1_,beta2_]))
                print('BMD 10: ' + str(bmd))
                print('BMD 50: ' + str(bmd_50))
                
            else:
                pred_vals = np.nan
                bmd = np.nan
                bmd_50 = np.nan
            
            if((bmdl_analysis_flag) and (res.mle_retvals['converged'] is True)):
                print('Estimating BMDL ...')
                # Using a bisection method for finding the BMDL value
                bmdl_val_lo = bmd/10
                bmdl_val_hi = bmd

                tol = 1
                bmdl_iter_count = 0
                while ((tol > tol_thresh) and (bmdl_iter_count < bmdl_max_iter)):
                    bmdl_iter_count+=1

                    # If maximum iterations are reached, set convergence flag to False and quit
                    if(bmdl_iter_count == bmdl_max_iter):
                        bmdl_converge_flag = 1
                        bmdl_mid_val = np.nan
                        break

                    bmdl_mid_val = (bmdl_val_lo + bmdl_val_hi)/2
                    print('bmdl_mid_val: ' + str(bmdl_mid_val))

                    # Examine the sign of (llv at bmdl_mid_val) - bmdl_llv_thresh
                    # and update bmdl_lo and hi values appropriately
                    model_multistage_2_bmd = baf.Multistage_2_BMD(test_dose_response[['dose','num_affected', 'total_num']].copy())
                    multistage_2_bmr_fit = model_multistage_2_bmd.profile_ll_fit([g_, beta1_, bmdl_mid_val])

                    # Check model convergence
                    if(multistage_2_bmr_fit.mle_retvals['converged'] is True):
                        bmdl_converge_flag = 0
                    else:
                        bmdl_converge_flag = 1
                        break
                    
                    bmdl_llv = multistage_2_bmr_fit.llf

                    if((bmdl_llv - bmdl_llv_thresh)>0):
                        bmdl_val_hi = bmdl_mid_val
                    else:
                        bmdl_val_lo = bmdl_mid_val

                    tol = abs(bmdl_llv - bmdl_llv_thresh)    
                print('BMDL Convergence: ' + str(not(bool(bmdl_converge_flag))))
                
        elif(model_name == 'quantal_linear'):
            model = baf.Quantal_Linear(test_dose_response[['dose','num_affected', 'total_num']].copy())        
            res = model.fit()
            
            print('Model Convergence: ' + str(res.mle_retvals['converged']))            
            # Check model convergence
            if(res.mle_retvals['converged'] is True):
                model_converge_flag = 0            
            
                pred_vals = baf.quantal_linear_fun(test_dose_response.dose, res.params)
                
                alpha_chi_square = 0.9
                bmdl_llv_thresh = res.llf - stats.chi2.ppf(alpha_chi_square, 1)/2
                
                g_ = res.params[0]
                beta_ = res.params[1]
            
                # Estimate BMD
                bmd = -np.log(1 - BMR)/beta_
                bmd_50 = -np.log(1 - BMR_50)/beta_
                print('Fit Paramas:' + str([g_,beta_]))
                print('BMD 10: ' + str(bmd))
                print('BMD 50: ' + str(bmd_50))
                
            else:
                pred_vals = np.nan
                bmd = np.nan
                bmd_50 = np.nan
            
            if((bmdl_analysis_flag) and (res.mle_retvals['converged'] is True)):
                print('Estimating BMDL ...')
                # Using a bisection method for finding the BMDL value
                bmdl_val_lo = bmd/10
                bmdl_val_hi = bmd

                tol = 1
                bmdl_iter_count = 0
                while ((tol > tol_thresh) and (bmdl_iter_count < bmdl_max_iter)):
                    bmdl_iter_count+=1

                    # If maximum iterations are reached, set convergence flag to False and quit
                    if(bmdl_iter_count == bmdl_max_iter):
                        bmdl_converge_flag = 1
                        bmdl_mid_val = np.nan
                        break

                    bmdl_mid_val = (bmdl_val_lo + bmdl_val_hi)/2

                    # Examine the sign of (llv at bmdl_mid_val) - bmdl_llv_thresh
                    # and update bmdl_lo and hi values appropriately
                    model_quantal_linear_bmd = baf.Quantal_Linear_BMD(test_dose_response[['dose','num_affected', 'total_num']].copy())
                    quantal_linear_bmr_fit = model_quantal_linear_bmd.profile_ll_fit([g_, bmdl_mid_val])

                    # Check model convergence
                    if(quantal_linear_bmr_fit.mle_retvals['converged'] is True):
                        bmdl_converge_flag = 0
                    else:
                        bmdl_converge_flag = 1
                        break
                    
                    bmdl_llv = quantal_linear_bmr_fit.llf

                    if((bmdl_llv - bmdl_llv_thresh)>0):
                        bmdl_val_hi = bmdl_mid_val
                    else:
                        bmdl_val_lo = bmdl_mid_val

                    tol = abs(bmdl_llv - bmdl_llv_thresh)        
                print('BMDL Convergence: ' + str(not(bool(bmdl_converge_flag))))
            
        # Estimate fit statistics
        
        if((model_converge_flag == 0) and (bmd > 0)):
            
            param_est = res.params
            num_params = len(param_est)
            AIC_val = -2*(res.llf) + (2 * num_params)
            exp_vals = test_dose_response.num_affected/test_dose_response.total_num
            num_vals = test_dose_response.total_num
            chi_square = 0
            scaled_residuals = np.zeros(len(exp_vals))
        
            for index in range(len(num_vals)):
                scaled_residuals[index] = (pred_vals[index] - exp_vals[index])/np.sqrt((pred_vals[index]*(1-pred_vals[index]))/num_vals[index])
                chi_square+=(num_vals[index]/(pred_vals[index]*(1-pred_vals[index])))*(exp_vals[index] - pred_vals[index])**2
        else:
            param_est = np.nan
            scaled_residuals = np.array([np.nan]*len(test_dose_response['dose']))
            chi_square = np.nan
            AIC_val = np.nan
            bmd = np.nan
            bmdl_mid_val = np.nan
            model_converge_flag = 1
                

        if(bmdl_converge_flag):
            bmdl_mid_val = np.nan
            print('BMDL: ' + 'NaN' + '\n')
        else:
            print('BMDL: ' + str(bmdl_mid_val) + '\n')
                
            
        # Populate data reporting dataframe
        model_predictions.loc[model_index] = [model_name, chi_square, stats.chi2.sf(chi_square, len(test_dose_response.dose) - len(res.params)), \
                                              AIC_val, bmd, bmdl_mid_val, bmd_50, scaled_residuals, param_est, model_converge_flag, bmdl_converge_flag]

        
    return model_predictions

def select_model(model_preds):
    '''Assumes that the min BMDL value yields just one model'''
    
    selected_model = 'None'
    no_unique_model_found_flag = 0
    model_select_flag = 0
    pval_thresh_flag = 0
    
    # 1: No model fit converged
    # 2: No fitted models satisfied p-val threshold
    # 3: Multiple models were found to have same AIC and BMD values but no valid BMDL values
    # 4: Multiple models found, user advised to look at the results of analysis    
    
    # Discard models for which fit convergence was not achieved
    model_preds_valid = model_preds[model_preds['Conv Flag'] == 0]
    
    if(model_preds_valid.empty is True):
        no_unique_model_found_flag = 1
        model_select_flag = 1
    else:
        #Discard models with p-value < 0.1
        model_preds_pvalue_cutoff = model_preds_valid[model_preds_valid['p-val'] > 0.1]
        
        if(model_preds_pvalue_cutoff.empty is True):
            #no_unique_model_found_flag = 1
            model_select_flag = 2
            pval_thresh_flag = 1
            model_preds_aic_cutoff = model_preds_valid[model_preds_valid['AIC'] == min(model_preds_valid['AIC'])]
        else:
            model_preds_aic_cutoff = model_preds_pvalue_cutoff[model_preds_pvalue_cutoff['AIC'] == min(model_preds_pvalue_cutoff['AIC'])]
            
        if(model_preds_aic_cutoff.shape[0] == 1):
            #Best unique model found
            selected_model = model_preds_aic_cutoff['Model']
        else:
            #Pick the model with minimum BMD
            selected_model_bmd_cutoff = model_preds_aic_cutoff[model_preds_aic_cutoff['Conv Flag'] == min(model_preds_aic_cutoff['Conv Flag'])]
            if(selected_model_bmd_cutoff.shape[0] == 1):
                # Best unique model found
                selected_model = selected_model_bmd_cutoff['Model']
            else:
                # Select model with min BMDL
                selected_model_valid_BMDL = selected_model_bmd_cutoff[selected_model_bmd_cutoff['BMDL Conv'] == 0]
                # If no models exist then the search terminates
                if(selected_model_valid_BMDL.empty is True):
                    no_unique_model_found_flag = 1
                    model_select_flag = 3
                else:
                    # Select the model with lowest BMDL value, if there are models with valid BMDLs
                    selected_model_bmdl = selected_model_valid_BMDL[selected_model_valid_BMDL['BMDL Conv'] == min(selected_model_valid_BMDL['BMDL Conv'])]
                    if(selected_model_bmdl.shape[0] == 1):
                        selected_model = selected_model_bmdl['Model']
                    else:
                        no_unique_model_found_flag = 1
                        model_select_flag = 4
                            
    selected_model_params = dict();  
    selected_model_params['model'] = selected_model
    selected_model_params['no_unique_model_found_flag'] = no_unique_model_found_flag
    selected_model_params['model_select_flag'] = model_select_flag
    selected_model_params['pval_thresh_flag'] = pval_thresh_flag
                           
    return selected_model_params