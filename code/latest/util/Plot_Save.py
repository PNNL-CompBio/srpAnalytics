#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Paritosh Pande
Pacific Northwest National Lab, Richland, WA
Original created on: May 2020
Modified: June 2020, Added AUC calculation and adding dose range to results
"""

import datetime
import numpy as np
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
from astropy import stats as astrostats
import BMD_Analysis_Functions as baf
import os

import warnings
warnings.filterwarnings('ignore')

from datetime import datetime as dt

today = dt.now()  
time_now_date = today.strftime('%Y_%m_%d')

report = True
#report = False

def save_results_poor_data_or_no_convergence(test_dose_response, qc_flag, chemical_id, end_point, selected_models = None):
    # Create the PdfPages object to which we will save the pages:
    # The with statement makes sure that the PdfPages object is closed properly at
    # the end of the block, even if an Exception occurs.
    #print(test_dose_response)
    #print(qc_flag)
    #print(chemical_id)
    #print(end_point)
    #print(selected_models)
    
    # Estimate AUC and min and mox doses
    if(not test_dose_response.empty):
        dose_response_auc = np.trapz(test_dose_response.num_affected/test_dose_response.total_num, x = test_dose_response.dose)
        dose_min = min(test_dose_response.dose)
        dose_max = max(test_dose_response.dose)
        dose_response_auc_norm = dose_response_auc/(dose_max - dose_min)
    else:
        dose_response_auc = np.nan
        dose_min = np.nan
        dose_max = np.nan 
        dose_response_auc_norm = np.nan
    
    if(not isinstance(chemical_id, str)):
        chemical_id = str(chemical_id)    
    filename = chemical_id + '_' + end_point + '.pdf'
    
    # Create dictionaries for various flags
    data_qc_flag_vals = {0 : 'Not enough dose groups for BMD analysis.' + '\n ' + 'BMD analysis not performed.',
                         1 : 'No trend detected in dose-response data.' + '\n' + 'BMD analysis not performed.',
                         2 : 'Dose-response data quality very good.',
                         3 : 'Dose-response data quality good.',
                         4 : 'Data resolution poor.' + '\n' + 'Caution advised.',
                         5 : 'Negative correlation detected in dose-response data.' + '\n' + 'Caution advised.' }
 
    # Filenames for csv files containing the results of analysis
    bmd_vals_file_name = 'bmd_vals_' + str(time_now_date) + '.csv'
    dose_response_vals_file_name = 'dose_response_vals_' + str(time_now_date) + '.csv'
    fit_vals_file_name = 'fit_vals_' + str(time_now_date) + '.csv'
    
    # Generate text for report
    if(selected_models is not None):
        text_for_report = 'Convergence not achieved for any dose-response model.'
    else:
        text_for_report = data_qc_flag_vals[qc_flag]
        
    with PdfPages(filename) as pdf:
        # Output data summary
        fig, ax = plt.subplots()
        # hide axes
        fig.patch.set_visible(False)

        ax.axis('off')        
        ax.axis('tight')
        #fig.text(0.1,0.7,' '.join(map(str, text_for_report)), transform=fig.transFigure, size=10, ha="left")
        fig.text(0.1,0.7,text_for_report, transform=fig.transFigure, size=10, ha="left")
        plt.title('Summary of Analysis')
        pdf.savefig()
        plt.close()
         
        # Plot dose-response data        
        CI_bounds = np.zeros([2, len(test_dose_response.dose)])
        # in save_results_poor_data_or_no_convergence fn
        for index in range(len(test_dose_response.dose)):
            print (f"test_dose_response.num_affected[index]:{test_dose_response.num_affected[index]}")
            print (f"test_dose_response.total_num[index]:{test_dose_response.total_num[index]}")
            CI = astrostats.binom_conf_interval(test_dose_response.num_affected[index], test_dose_response.total_num[index], confidence_level = 0.95)
            CI = np.abs(CI - test_dose_response.num_affected[index]/test_dose_response.total_num[index])
            CI_bounds[0, index] = CI[0]
            CI_bounds[1, index] = CI[1]
        fig, ax = plt.subplots()
        
        # Setting the values for all axes.
        custom_ylim = (0,1)
        plt.setp(ax, ylim=custom_ylim)
        
        ax.set_xscale("linear")
        ax.errorbar(test_dose_response.dose, test_dose_response.num_affected/test_dose_response.total_num, CI_bounds, marker ='s', mfc='red', fmt='.')
        
        ax.set_xlabel('Dose')
        ax.set_ylabel('Fractional Response')
        ax.set_title('Dose-response Data')
        pdf.savefig()  # saves the current figure into a pdf page
        plt.close()
        
        # Create dataframes to apprend to write to csv files
        bmd_vals = pd.DataFrame(columns = ['Chemical_ID', 'End_Point', 'Model', 'BMD10', 'BMDL', 'BMD50', 'AUC', 'Min_Dose', 'Max_Dose', 'AUC_Norm', 'DataQC_Flag', 'BMD_Analysis_Flag', 'BMD10_Flag', 'BMD50_Flag'])
        dose_response_vals = pd.DataFrame(columns = ['Chemical_ID', 'End_Point', 'Dose', 'Response', 'CI_Lo', 'CI_Hi'])
        fit_vals = pd.DataFrame(columns = ['Chemical_ID', 'End_Point', 'X_vals', 'Y_vals'])
        
        # Populate dataframes
        bmd_vals['Chemical_ID'] = [chemical_id]
        bmd_vals['End_Point'] = [end_point]
        bmd_vals['Model'] = np.nan
        bmd_vals['BMD10'] = np.nan
        bmd_vals['BMDL'] = np.nan
        bmd_vals['BMD50'] = np.nan
        bmd_vals['DataQC_Flag'] = qc_flag
        bmd_vals['AUC'] = dose_response_auc
        bmd_vals['Min_Dose'] = dose_min
        bmd_vals['Max_Dose'] = dose_max
        bmd_vals['AUC_Norm'] = dose_response_auc_norm
        bmd_vals['BMD_Analysis_Flag'] = np.nan
        bmd_vals['BMD10_Flag'] = np.nan
        bmd_vals['BMD50_Flag'] = np.nan
        
        assign_nan = False
        try: # 53_ANY24
            bogus = test_dose_response.dose[0]
            #print ("test_dose_response.dose[0]:"+str(test_dose_response.dose[0]))                
        except: # 1532_ANY24
            assign_nan = True
            # print ("test_dose_response.dose:"+str(test_dose_response.dose))
            # Series([], Name: dose, dtype: object)
        
        if (assign_nan):            
            dose_response_vals['Chemical_ID'] = [chemical_id]
            dose_response_vals['End_Point'] = [end_point]
            dose_response_vals['Dose'] = np.nan    
            dose_response_vals['Response'] = np.nan
            dose_response_vals['CI_Lo'] = np.nan
            dose_response_vals['CI_Hi'] = np.nan
        else:
            dose_response_vals['Chemical_ID'] = [chemical_id]*len(test_dose_response.dose)
            dose_response_vals['End_Point'] = [end_point]*len(test_dose_response.dose)
            dose_response_vals['Dose'] = test_dose_response.dose
            dose_response_vals['Response'] = test_dose_response.num_affected/test_dose_response.total_num
            dose_response_vals['CI_Lo'] = CI_bounds[0, :]
            dose_response_vals['CI_Hi'] = CI_bounds[1, :]
        
        fit_vals['Chemical_ID'] = [chemical_id]
        fit_vals['End_Point'] = [end_point]
        fit_vals['X_vals'] = np.nan
        fit_vals['Y_vals'] = np.nan
        
        if not os.path.isfile(bmd_vals_file_name):
            bmd_vals.to_csv(bmd_vals_file_name, header='column_names', index=False,na_rep='NULL')
        else: # else it exists so append without writing the header
            bmd_vals.to_csv(bmd_vals_file_name, mode='a', header=False, index=False,na_rep='NULL')
            
        if not os.path.isfile(dose_response_vals_file_name):
            dose_response_vals.to_csv(dose_response_vals_file_name, header='column_names', index=False,na_rep='NULL')
        else: # else it exists so append without writing the header
            dose_response_vals.to_csv(dose_response_vals_file_name, mode='a', header=False, index=False, na_rep='NULL')
            
        if not os.path.isfile(fit_vals_file_name):
            fit_vals.to_csv(fit_vals_file_name, header='column_names', index=False, na_rep='NULL')
        else: # else it exists so append without writing the header
            fit_vals.to_csv(fit_vals_file_name, mode='a', header=False, index=False, na_rep='NULL')
        
        # We can also set the file's metadata via the PdfPages object:
        d = pdf.infodict()
        d['Author'] = 'Paritosh Pande'
        d['CreationDate'] = datetime.datetime.today()
###### end of save_results_poor_data_or_no_convergence()


def save_results_good_data_unique_model(test_dose_response, qc_flag, model_preds, selected_models, chemical_id, end_point):
    # Create the PdfPages object to which we will save the pages:
    # The with statement makes sure that the PdfPages object is closed properly at
    # the end of the block, even if an Exception occurs.
    
    # Estimate AUC and min and mox doses
    if(not test_dose_response.empty):        
        dose_response_auc = np.trapz(test_dose_response.num_affected/test_dose_response.total_num, x = test_dose_response.dose)
        dose_min = min(test_dose_response.dose)
        dose_max = max(test_dose_response.dose)
        dose_response_auc_norm = dose_response_auc/(dose_max - dose_min)
    else:
        dose_response_auc = np.nan
        dose_min = np.nan
        dose_max = np.nan 
        dose_response_auc_norm = np.nan
    
    if(not isinstance(chemical_id, str)):
        chemical_id = str(chemical_id)
        
    filename = chemical_id + '_' + end_point + '.pdf'
    model_preds = model_preds.round(8)
    
    # Extract subset of results table 
    model_preds_basic_stats = model_preds[['Model', 'Chi-squared', 'p-val', 'AIC', 'BMD10', 'BMDL10']]
    residual_column_names = [('dose'+ str(i)) for i in range(len(test_dose_response['dose']))] 
    model_preds_residuals = pd.DataFrame(columns = ['Model'] + residual_column_names)   
    model_preds_residuals['Model'] = model_preds['Model']
    model_preds_residuals_matrix = np.empty((model_preds['Scaled Residuals'].shape[0],len(test_dose_response['dose'])))
    model_preds_residuals_matrix[:] = np.nan
    model_preds_residuals[residual_column_names] = model_preds_residuals_matrix
    

    for model_pred_index in range(model_preds['Scaled Residuals'].shape[0]):


        if(not any(np.isnan(model_preds['Scaled Residuals'][model_pred_index]))):
            if (report):
                print (f"model_preds_residuals:\n{model_preds_residuals}")
                print (f"model_pred_index:\n{model_pred_index}")
                
                print (f"np.matrix(model_preds['Scaled Residuals'][model_pred_index].tolist()).round(8):\n{np.matrix(model_preds['Scaled Residuals'][model_pred_index].tolist()).round(8)}")
                #[[-0.71275987 -0.04841195  1.2423122   0.22264199  0.0676003  -0.32941879  -1.42086953  1.03044301]]
            
            model_preds_residuals.iloc[model_pred_index,1:] = np.matrix(model_preds['Scaled Residuals'][model_pred_index].tolist()).round(8)

    # Create dictionaries for various flags
    data_qc_flag_vals = {0 : 'Not enough dose groups for BMD analysis.' + '\n ' + 'BMD analysis not performed.',
                         1 : 'No trend detected in dose-response data.' + '\n' + 'BMD analysis not performed.',
                         2 : 'Dose-response data quality very good.',
                         3 : 'Dose-response data quality good.',
                         4 : 'Data resolution poor. Caution advised.',
                         5 : 'Negative correlation detected in dose-response data.' + '\n' + 'Caution advised.' }

    bmd_analysis_flag_vals = {1 : 'Convergence not achieved for any dose-response model.',
                              2 : 'Model fit might be unreliable.' + '\n' + 'p-val for chi-squared statistic was < 0.1 for all converged models.',
                              3 : 'A unique model could not be determined.' + '\n' + 'Multiple models had the same AIC and BMD values but no valid BMDL values.',
                              4 : 'Multiple models found.' + '\n' + 'User advised to look at the results of analysis to choose the best model.'}

    txt_for_model_selection = 'Best model found:' + selected_models['model'].values 
    unique_model_flag_vals = {0 : txt_for_model_selection,
                              1 : 'Best model could not be determined'}
    
    bmd_analysis_flag = selected_models['model_select_flag']
    unique_model_flag = selected_models['no_unique_model_found_flag']
    
    # Filenames for csv files containing the results of analysis
    bmd_vals_file_name = 'bmd_vals_' + str(time_now_date) + '.csv'
    dose_response_vals_file_name = 'dose_response_vals_' + str(time_now_date) + '.csv'
    fit_vals_file_name = 'fit_vals_' + str(time_now_date) + '.csv'
    
    text_for_report = data_qc_flag_vals[qc_flag]
    
    # Generate text for report
    if((unique_model_flag == 0) and (bmd_analysis_flag != 2)):
        text_for_report = text_for_report + '\n' + unique_model_flag_vals[unique_model_flag]
    elif((unique_model_flag == 0) and (bmd_analysis_flag == 2)):
        text_for_report = text_for_report + '\n' + unique_model_flag_vals[unique_model_flag] + '\n' + bmd_analysis_flag_vals[bmd_analysis_flag]    
    else:
        # Specify reason for non-uniqueness
        text_for_report = text_for_report + '\n' + unique_model_flag_vals[unique_model_flag] + '\n' + bmd_analysis_flag_vals[bmd_analysis_flag]
      
    with PdfPages(filename) as pdf:
        # Output data summary
        fig, ax = plt.subplots()
        # hide axes
        fig.patch.set_visible(False)
        
        ax.axis('off')
        ax.axis('tight')
        fig.text(0.1,0.7,' '.join(map(str, text_for_report)), transform=fig.transFigure, size=10, ha="left")
        plt.title('Summary of Analysis')   
        pdf.savefig()
        plt.close()
        
        # Print Model Predictions
        fig, ax = plt.subplots()
        # hide axes
        fig.patch.set_visible(False)
        
        ax.axis('off')
        ax.axis('tight')
        ax.table(cellText=model_preds_basic_stats.values, colLabels=model_preds_basic_stats.columns, loc='center')
        
        plt.title('Model Predictions')
        fig.tight_layout()
        pdf.savefig()  # saves the current figure into a pdf page
        plt.close()

        # Print residuals for different models        
        fig, ax = plt.subplots()
        # hide axes
        fig.patch.set_visible(False)
        
        ax.axis('off')
        ax.axis('tight')
        ax.table(cellText=model_preds_residuals.values, colLabels=model_preds_residuals.columns, loc='center')
        plt.title('Scaled Residuals')
        fig.tight_layout()
        pdf.savefig()  # saves the current figure into a pdf page
        plt.close()
        
        # Extract data for best model found and save it for portal
        # and plot fit for selected model
        model_name = selected_models['model'].values
        optimized_params = model_preds.loc[model_preds['Model'] == model_name[0], 'Optimized Params'].values[0]
                
        CI_bounds = np.zeros([2, len(test_dose_response.dose)])
        
        for index in range(len(test_dose_response.dose)):
            CI = astrostats.binom_conf_interval(test_dose_response.num_affected[index], test_dose_response.total_num[index], confidence_level = 0.95)
            CI = np.abs(CI - test_dose_response.num_affected[index]/test_dose_response.total_num[index])
            CI_bounds[0, index] = CI[0]
            CI_bounds[1, index] = CI[1]

        fig, ax = plt.subplots()

        # Setting the values for all axes.
        custom_ylim = (0,1)
        plt.setp(ax, ylim=custom_ylim)

        ax.set_xscale("linear")
        ax.errorbar(test_dose_response.dose, test_dose_response.num_affected/test_dose_response.total_num, CI_bounds, marker ='s', mfc='red', fmt='.')
        
        ax.set_xlabel('Dose')
        ax.set_ylabel('Fractional Response')
        ax.set_title(' '.join(map(str, 'Dose-response with best fit model (' + model_name + ')')))
        
        int_steps = 10
        dose_x_vals = gen_uneven_spacing(test_dose_response.dose, int_steps)
        np.append(dose_x_vals, dose_x_vals[-1] + (dose_x_vals[-1] - dose_x_vals[-2])/int_steps)

        if(model_name != 'None'):
            if(model_name == 'logistic'):
                ax.plot(dose_x_vals, baf.logistic_fun(dose_x_vals, optimized_params),'b-')
                y_vals = baf.logistic_fun(dose_x_vals, optimized_params)
            elif(model_name == 'log_logistic'):
                ax.plot(dose_x_vals, baf.log_logistic_fun(dose_x_vals, optimized_params),'b-')
                y_vals = baf.log_logistic_fun(dose_x_vals, optimized_params)
            elif(model_name == 'gamma'):
                ax.plot(dose_x_vals, baf.gamma_fun(dose_x_vals, optimized_params),'b-')
                y_vals = baf.gamma_fun(dose_x_vals, optimized_params)
            elif(model_name == 'weibull'):
                ax.plot(dose_x_vals, baf.weibull_fun(dose_x_vals, optimized_params),'b-')
                y_vals = baf.weibull_fun(dose_x_vals, optimized_params)
            elif(model_name == 'probit'):
                ax.plot(dose_x_vals, baf.probit_fun(dose_x_vals, optimized_params),'b-')
                y_vals = baf.probit_fun(dose_x_vals, optimized_params)
            elif(model_name == 'log_probit'):
                ax.plot(dose_x_vals, baf.log_probit_fun(dose_x_vals, optimized_params),'b-')
                y_vals = baf.log_probit_fun(dose_x_vals, optimized_params)
            elif(model_name == 'multistage_2'):
                ax.plot(dose_x_vals, baf.multistage_2_fun(dose_x_vals, optimized_params),'b-') 
                y_vals = baf.multistage_2_fun(dose_x_vals, optimized_params)
            elif(model_name == 'quantal_linear'):
                ax.plot(dose_x_vals, baf.quantal_linear_fun(dose_x_vals, optimized_params),'b-') 
                y_vals = baf.quantal_linear_fun(dose_x_vals, optimized_params)
        
        pdf.savefig()  # saves the current figure into a pdf page
        plt.close()
        
        # Create dataframes to apprend to write to csv files
        bmd_vals = pd.DataFrame(columns = ['Chemical_ID', 'End_Point', 'Model', 'BMD10', 'BMDL', 'BMD50', 'AUC', 'Min_Dose', 'Max_Dose', 'AUC_Norm', 'DataQC_Flag', 'BMD_Analysis_Flag', 'BMD10_Flag', 'BMD50_Flag'])
        dose_response_vals = pd.DataFrame(columns = ['Chemical_ID', 'End_Point', 'Dose', 'Response', 'CI_Lo', 'CI_Hi'])
        #fit_vals = pd.DataFrame(columns = ['Chemical_ID', 'End_Point', 'X_vals', 'Y_vals', 'Y_vals_diff'])
        fit_vals = pd.DataFrame(columns = ['Chemical_ID', 'End_Point', 'X_vals', 'Y_vals'])
        
        # Populate dataframes
        bmd_vals['Chemical_ID'] = [chemical_id]
        bmd_vals['End_Point'] = [end_point]
        bmd_vals['Model'] = model_name
        bmd_vals['BMD10'] = model_preds.loc[model_preds['Model'] == model_name[0], 'BMD10'].values
        bmd_vals['BMDL'] = model_preds.loc[model_preds['Model'] == model_name[0], 'BMDL10'].values
        bmd_vals['BMD50'] = model_preds.loc[model_preds['Model'] == model_name[0], 'BMD50'].values
        bmd_vals['DataQC_Flag'] = qc_flag
        bmd_vals['AUC'] = dose_response_auc
        bmd_vals['Min_Dose'] = dose_min
        bmd_vals['Max_Dose'] = dose_max
        bmd_vals['AUC_Norm'] = dose_response_auc_norm
        bmd_vals['BMD_Analysis_Flag'] = bmd_analysis_flag

        if(model_preds.loc[model_preds['Model'] == model_name[0], 'BMD10'].values < test_dose_response.dose[1]):
            bmd_vals['BMD10_Flag'] = -1
        elif(model_preds.loc[model_preds['Model'] == model_name[0], 'BMD10'].values > test_dose_response.dose.iloc[-1]):
            bmd_vals['BMD10_Flag'] = 1
        else:
            bmd_vals['BMD10_Flag'] = 0
            
        if(model_preds.loc[model_preds['Model'] == model_name[0], 'BMD50'].values < test_dose_response.dose[1]):
            bmd_vals['BMD50_Flag'] = -1
        elif(model_preds.loc[model_preds['Model'] == model_name[0], 'BMD50'].values > test_dose_response.dose.iloc[-1]):
            bmd_vals['BMD50_Flag'] = 1
        else:
            bmd_vals['BMD50_Flag'] = 0
        
        dose_response_vals['Chemical_ID'] = [chemical_id]*len(test_dose_response.dose)
        dose_response_vals['End_Point'] = [end_point]*len(test_dose_response.dose)
        dose_response_vals['Dose'] = test_dose_response.dose
        dose_response_vals['Response'] = test_dose_response.num_affected/test_dose_response.total_num
        dose_response_vals['CI_Lo'] = CI_bounds[0, :]
        dose_response_vals['CI_Hi'] = CI_bounds[1, :]
        
        if (report):
            print(len(dose_x_vals))
            print(len(y_vals))
        
        fit_vals['Chemical_ID'] = [chemical_id]*len(dose_x_vals)
        fit_vals['End_Point'] = [end_point]*len(dose_x_vals)
        fit_vals['X_vals'] = dose_x_vals
        fit_vals['Y_vals'] = y_vals
        #fit_vals['Y_vals_diff'] = y_vals
        
        if not os.path.isfile(bmd_vals_file_name):
            bmd_vals.to_csv(bmd_vals_file_name, header='column_names', index=False, na_rep='NULL')
        else: # else it exists so append without writing the header
            bmd_vals.to_csv(bmd_vals_file_name, mode='a', header=False, index=False, na_rep='NULL')
        
        if not os.path.isfile(dose_response_vals_file_name):
            dose_response_vals.to_csv(dose_response_vals_file_name, header='column_names', index=False, na_rep='NULL')
        else: # else it exists so append without writing the header
            dose_response_vals.to_csv(dose_response_vals_file_name, mode='a', header=False, index=False, na_rep='NULL')
        
        if not os.path.isfile(fit_vals_file_name):
            fit_vals.to_csv(fit_vals_file_name, header='column_names', index=False, na_rep='NULL')
        else: # else it exists so append without writing the header
            fit_vals.to_csv(fit_vals_file_name, mode='a', header=False, index=False, na_rep='NULL')
        
        # We can also set the file's metadata via the PdfPages object:
        d = pdf.infodict()
        d['Author'] = 'Paritosh Pande'
        d['CreationDate'] = datetime.datetime.today()
        
def save_results_good_data_nounique_model(test_dose_response, qc_flag, model_preds, selected_models, chemical_id, end_point):
    # Create the PdfPages object to which we will save the pages:
    # The with statement makes sure that the PdfPages object is closed properly at
    # the end of the block, even if an Exception occurs.
    
    # Estimate AUC and min and mox doses
    if(not test_dose_response.empty):
        print ("test_dose_response:"+str(test_dose_response))
        
        dose_response_auc = np.trapz(test_dose_response.num_affected/test_dose_response.total_num, x = test_dose_response.dose)
        dose_min = min(test_dose_response.dose)
        dose_max = max(test_dose_response.dose)
        dose_response_auc_norm = dose_response_auc/(dose_max - dose_min)
    else:
        dose_response_auc = np.nan
        dose_min = np.nan
        dose_max = np.nan
        dose_response_auc_norm = np.nan
    
    if(not isinstance(chemical_id, str)):
        chemical_id = str(chemical_id)
        
    filename = chemical_id + '_' + end_point + '.pdf'
    model_preds = model_preds.round(8)
    
    # Extract subset of results table 
    model_preds_basic_stats = model_preds[['Model', 'Chi-squared', 'p-val', 'AIC', 'BMD10', 'BMDL10']]

    residual_column_names = [('dose'+ str(i)) for i in range(len(test_dose_response['dose']))] 
    model_preds_residuals = pd.DataFrame(columns = ['Model'] + residual_column_names)   
    model_preds_residuals['Model'] = model_preds['Model']
    model_preds_residuals_matrix = np.empty((model_preds['Scaled Residuals'].shape[0],len(test_dose_response['dose'])))
    model_preds_residuals_matrix[:] = np.nan
    
    model_preds_residuals[residual_column_names] = model_preds_residuals_matrix
    for model_pred_index in range(model_preds['Scaled Residuals'].shape[0]):
        if(not any(np.isnan(model_preds['Scaled Residuals'][model_pred_index]))):
            model_preds_residuals.iloc[model_pred_index,1:] = np.matrix(model_preds['Scaled Residuals'][model_pred_index].tolist()).round(8)

        
    # Create dictionaries for various flags
    data_qc_flag_vals = {0 : 'Not enough dose groups for BMD analysis.' + '\n '+ 'BMD analysis not performed.',
                         1 : 'No trend detected in dose-response data.' + '\n' + 'BMD analysis not performed.',
                         2 : 'Dose-response data quality very good.',
                         3 : 'Dose-response data quality good.',
                         4 : 'Data resolution poor.' + '\n' + 'Caution advised.',
                         5 : 'Negative correlation detected in dose-response data.' + '\n' + 'Caution advised.' }

    bmd_analysis_flag_vals = {1 : 'Convergence not achieved for any dose-response model.',
                              2 : 'Model fit might be unreliable.' + '\n' + 'p-val for chi-squared statistic was < 0.1 for all converged models.',
                              3 : 'A unique model could not be determined.' + '\n' + 'Multiple models had the same AIC and BMD values but no valid BMDL values.',
                              4 : 'Multiple models found.' + '\n' + 'User advised to look at the results of analysis to choose the best model.'}

    #txt_for_model_selection = selected_models['model'].values + ' determined to be the best model'
    unique_model_flag_vals = {0 : 'None',
                              1 : 'Best model could not be determined'}
    
    bmd_analysis_flag = selected_models['model_select_flag']
    unique_model_flag = selected_models['no_unique_model_found_flag']
    
    # Filenames for csv files containing the results of analysis
    bmd_vals_file_name = 'bmd_vals_' + str(time_now_date) + '.csv'
    dose_response_vals_file_name = 'dose_response_vals_' + str(time_now_date) + '.csv'
    fit_vals_file_name = 'fit_vals_' + str(time_now_date) + '.csv'
    
    # Generate text for report
    text_for_report = data_qc_flag_vals[qc_flag]
    
    # Specify reason for non-unique model
    text_for_report += '\n' + unique_model_flag_vals[unique_model_flag]
    text_for_report += '\n' + bmd_analysis_flag_vals[bmd_analysis_flag]

    with PdfPages(filename) as pdf:
        # Output data summary
        fig, ax = plt.subplots()
        # hide axes
        fig.patch.set_visible(False)
        
        ax.axis('off')
        ax.axis('tight')
        fig.text(0.1,0.7,' '.join(map(str, text_for_report)), transform=fig.transFigure, size=10, ha="left")
        plt.title('Summary of Analysis')
        pdf.savefig()
        plt.close()        
        
        # Print Model Predictions
        fig, ax = plt.subplots()
        # hide axes
        fig.patch.set_visible(False)
        
        ax.axis('off')
        ax.axis('tight')
        ax.table(cellText=model_preds_basic_stats.values, colLabels=model_preds_basic_stats.columns, loc='center')
        
        plt.title('Model Predictions')
        fig.tight_layout()
        pdf.savefig()  # saves the current figure into a pdf page
        plt.close()

        # Print residuals for different models        
        fig, ax = plt.subplots()
        # hide axes
        fig.patch.set_visible(False)
        ax.axis('off')
        ax.axis('tight')

        ax.table(cellText=model_preds_residuals.values, colLabels=model_preds_residuals.columns, loc='center')
        
        plt.title('Scaled Residuals')
        fig.tight_layout()
        pdf.savefig()  # saves the current figure into a pdf page
        plt.close()
                
        CI_bounds = np.zeros([2, len(test_dose_response.dose)])
        for index in range(len(test_dose_response.dose)):
            CI = astrostats.binom_conf_interval(test_dose_response.num_affected[index], test_dose_response.total_num[index], confidence_level = 0.95)
            CI = np.abs(CI - test_dose_response.num_affected[index]/test_dose_response.total_num[index])
            CI_bounds[0, index] = CI[0]
            CI_bounds[1, index] = CI[1]

        fig, ax = plt.subplots()
        ax.set_xscale("linear")
        ax.errorbar(test_dose_response.dose, test_dose_response.num_affected/test_dose_response.total_num, CI_bounds, marker ='s', mfc='red', fmt='.')
        
        ax.set_xlabel('Dose')
        ax.set_ylabel('Fractional Response')
        ax.set_title('Dose-response Data')
        pdf.savefig()  # saves the current figure into a pdf page
        plt.close()
        
        # Create dataframes to apprend to write to csv files
        bmd_vals = pd.DataFrame(columns = ['Chemical_ID', 'End_Point', 'Model', 'BMD10', 'BMDL', 'BMD50', 'AUC', 'Min_Dose', 'Max_Dose', 'AUC_Norm', 'DataQC_Flag', 'BMD_Analysis_Flag', 'BMD10_Flag', 'BMD50_Flag'])
        dose_response_vals = pd.DataFrame(columns = ['Chemical_ID', 'End_Point', 'Dose', 'Response', 'CI_Lo', 'CI_Hi'])
        fit_vals = pd.DataFrame(columns = ['Chemical_ID', 'End_Point', 'X_vals', 'Y_vals'])
        
        # Populate dataframes
        bmd_vals['Chemical_ID'] = [chemical_id]
        bmd_vals['End_Point'] = [end_point]
        bmd_vals['Model'] = np.nan
        bmd_vals['BMD10'] = np.nan
        bmd_vals['BMDL'] = np.nan
        bmd_vals['BMD50'] = np.nan
        bmd_vals['DataQC_Flag'] = qc_flag
        bmd_vals['AUC'] = dose_response_auc
        bmd_vals['Min_Dose'] = dose_min
        bmd_vals['Max_Dose'] = dose_max
        bmd_vals['AUC_Norm'] = dose_response_auc_norm
        bmd_vals['BMD_Analysis_Flag'] = bmd_analysis_flag
        bmd_vals['BMD10_Flag'] = np.nan
        bmd_vals['BMD50_Flag'] = np.nan
                
        dose_response_vals['Chemical_ID'] = [chemical_id]*len(test_dose_response.dose)
        dose_response_vals['End_Point'] = [end_point]*len(test_dose_response.dose)
        dose_response_vals['Dose'] = test_dose_response.dose
        dose_response_vals['Response'] = test_dose_response.num_affected/test_dose_response.total_num
        dose_response_vals['CI_Lo'] = CI_bounds[0, :]
        dose_response_vals['CI_Hi'] = CI_bounds[1, :]
        
        fit_vals['Chemical_ID'] = [chemical_id]
        fit_vals['End_Point'] = [end_point]
        fit_vals['X_vals'] = np.nan
        fit_vals['Y_vals'] = np.nan
        
        if not os.path.isfile(bmd_vals_file_name):
            bmd_vals.to_csv(bmd_vals_file_name, header='column_names', index=False, na_rep='NULL')
        else: # else it exists so append without writing the header
            bmd_vals.to_csv(bmd_vals_file_name, mode='a', header=False, index=False, na_rep='NULL')
            
        if not os.path.isfile(dose_response_vals_file_name):
            dose_response_vals.to_csv(dose_response_vals_file_name, header='column_names', index=False, na_rep='NULL')
        else: # else it exists so append without writing the header
            dose_response_vals.to_csv(dose_response_vals_file_name, mode='a', header=False, index=False, na_rep='NULL')
            
        if not os.path.isfile(fit_vals_file_name):
            fit_vals.to_csv(fit_vals_file_name, header='column_names', index=False, na_rep='NULL')
        else: # else it exists so append without writing the header
            fit_vals.to_csv(fit_vals_file_name, mode='a', header=False, index=False, na_rep='NULL')
        
        # We can also set the file's metadata via the PdfPages object:
        d = pdf.infodict()
        d['Author'] = 'Paritosh Pande'
        d['CreationDate'] = datetime.datetime.today()

def gen_uneven_spacing(doses, int_steps):
    dose_samples = list()
    for dose_index in range(len(doses) - 1):
        dose_samples.extend(np.linspace(doses[dose_index],doses[dose_index + 1], int_steps).tolist())
    return np.unique(dose_samples)   


#save_results_good_data_unique_model(test_dose_response, qc_flag, model_preds, selected_models, chemical_id, end_point)