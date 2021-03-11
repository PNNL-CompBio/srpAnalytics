#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Paritosh Pande
Pacific Northwest National Lab, Richland, WA
Original created on: May 2020
"""

import numpy as np
import pandas as pd
from scipy import stats


# Get dose-respone data

def gen_dose_response(delta_mov_auc_data, end_point):
    # Remove the plates for which number of abnormal control wells > 50% total
    # control wells
      
    unique_plate_IDs = np.unique(delta_mov_auc_data['Plate'])
    #print('Unique number of plates:', len(unique_plate_IDs))
    
    abnormal_response_wells = [];
    
    # Count number of abnormal negative control wells wells for each chemical
    for plate_ID in unique_plate_IDs:
        plate_data_subset = delta_mov_auc_data.loc[delta_mov_auc_data['Plate'] == plate_ID]
    
        # Extract data for negative control wells
        plate_data_subset_nc = plate_data_subset.loc[plate_data_subset['CONC'] == 0]
    
        # Display the number of negative control wells for sanity check
        #print('Number of negative control wells on plate:', plate_ID, 'are:', plate_data_subset_nc.shape[0])
    
        # Count and display the number of wells with negative response values
        number_of_abnormal_nc_wells = (plate_data_subset_nc.loc[plate_data_subset_nc[end_point] < 0]).shape[0]
        #print('Number of abnormal response negative control wells:', number_of_abnormal_nc_wells)   

        if(number_of_abnormal_nc_wells >= (plate_data_subset_nc.shape[0])/2):
            abnormal_response_wells.append(plate_ID)
            print('Wells with abnormal response for negative control', abnormal_response_wells)

    # Remove data for abnormal wells
    delta_mov_auc_normal = delta_mov_auc_data[~delta_mov_auc_data['Plate'].isin(abnormal_response_wells)]

    # Pre-define dataframe for dose-respose
    dose_response_info = pd.DataFrame(columns = ['Concentration', 'Response', 
                                                 'Hypo', 'Hyper', 'Number_of_Wells'])
    dose_response = pd.DataFrame(columns = ['dose', 'num_affect', 'num_embryos'])
                                                 
    # Get data for a given concentration
    delta_mov_auc_data_compound = delta_mov_auc_normal
    all_neg_control_vals = delta_mov_auc_data_compound.loc[delta_mov_auc_data_compound['CONC'] == 0.0]

    concentrations = np.unique(delta_mov_auc_data_compound['CONC'])
    for concentration in concentrations:
        if concentration != 0.0:
            delta_mov_auc_data_compound_concentration = delta_mov_auc_data_compound.loc[delta_mov_auc_data_compound['CONC'] == concentration]
            plate_ids = np.unique(delta_mov_auc_data_compound_concentration['Plate'])
            #print("Plate(s) for concentration:", concentration, "are:", plate_ids)
            #print(delta_mov_auc_data_compound_concentration.shape)
            count_response_wells = 0
            count_total_wells = 0
            count_hypo_wells = 0
            count_hyper_wells = 0
            for plate_id in plate_ids:
                neg_control_plate_specific_vals = all_neg_control_vals.loc[all_neg_control_vals['Plate'] == plate_id]
                neg_control_ref_vals = neg_control_plate_specific_vals[end_point]
                response_vals_plate = delta_mov_auc_data_compound_concentration.loc[delta_mov_auc_data_compound_concentration['Plate'] == plate_id]
                response_vals = response_vals_plate[end_point]

                response_vals_positive = response_vals[response_vals >= 0]
                # Identify hypo- and hyperactive responses
                Q1_neg_control_vals = neg_control_ref_vals.quantile(0.25)
                Q3_neg_control_vals = neg_control_ref_vals.quantile(0.75)
                IQR_neg_control_vals = Q3_neg_control_vals - Q1_neg_control_vals
                hyper_response_vals = response_vals_positive[(response_vals_positive < (Q1_neg_control_vals - 1.5 * IQR_neg_control_vals)) | (response_vals_positive > (Q3_neg_control_vals + 1.5 * IQR_neg_control_vals))]
                hypo_response_vals = response_vals[response_vals <= 0]
                
                #Compare box-plots for visual inspection
                #fig, ax = plt.subplots()
                #bp = ax.boxplot([neg_control_ref_vals , response_vals])
                #ax.set_xticklabels(['Neg. Ctrl', 'Response'])
                count_hypo_wells+= len(hypo_response_vals)
                count_hyper_wells+= len(hyper_response_vals)
                count_response_wells += (len(hypo_response_vals) + len(hyper_response_vals))
                count_total_wells += len(response_vals)
                response_value = count_response_wells/count_total_wells
                
            # Populate dose response dataframe
            # Obtain the total number of wells for a compound and concentration before any filtering is performed
            dose_response_info = dose_response_info.append({'Concentration': concentration,
                                                            'Response': response_value ,
                                                            'Hypo': count_hypo_wells,
                                                            'Hyper': count_hyper_wells,
                                                            'Number_of_Wells': count_total_wells},
                                                             ignore_index = True)
            
            dose_response = dose_response.append({'dose': concentration,
                                                  'num_affect': (count_hypo_wells + count_hyper_wells),
                                                  'num_embryos': count_total_wells},
                                                   ignore_index = True)
                        
        else:
            # For concentration = 0.0 or negative controls
            delta_mov_auc_data_compound_concentration = delta_mov_auc_data_compound.loc[delta_mov_auc_data_compound['CONC'] == concentration]
            plate_ids = np.unique(delta_mov_auc_data_compound_concentration['Plate'])
            count_response_wells = 0
            count_total_wells = 0
            count_hypo_wells = 0
            count_hyper_wells = 0
            #print("Plate(s) for concentration:", concentration, "are:", plate_ids)
            for plate_id in plate_ids:
                neg_control_plate_specific_vals = all_neg_control_vals.loc[all_neg_control_vals['Plate'] == plate_id]
                neg_control_ref_vals = neg_control_plate_specific_vals[end_point]
                
                count_response_wells += len(neg_control_ref_vals[neg_control_ref_vals < 0])
                count_total_wells += len(neg_control_ref_vals)
                response_value = count_response_wells/count_total_wells
            
            # Populate dose response dataframe
            # Obtain the total number of wells for a compound and concentration before any filtering is performed
            dose_response_info = dose_response_info.append({'Concentration': concentration,
                                                            'Response': response_value ,
                                                            'Hypo': count_response_wells, 
                                                            'Hyper': np.nan,
                                                            'Number_of_Wells': count_total_wells},
                                                             ignore_index = True)

            dose_response = dose_response.append({'dose': concentration,
                                                  'num_affect': count_response_wells,
                                                  'num_embryos': count_total_wells},
                                                   ignore_index = True)
            
    return dose_response             

def BMD_feasibility_analysis(dose_response):
    '''This function performs feasibility analysis
    for dose respone data. The value returned is a 
    flag indicating data quality as defined below:
    0: Not enough dose groups for BMD analysis. BMD analysis not performed
    1: Not enough dose groups with response > 10%. BMD Analysis not performed
    2: Good dose-response data
    3: Dose-response data quality poor. BMD analysis might be unreliable
    4: Data resolution poor. BMD analysis might be unreliable
    5: No trend detected in dose-response data. BMD analysis not performed'''
    if(dose_response.shape[0] < 3):
        qc_flag = 0
    else:
        frac_response = dose_response['num_affect']/dose_response['num_embryos']      
        if((frac_response.iloc[-1] + frac_response.iloc[-2])/2.0 < (frac_response.iloc[0] + frac_response.iloc[1])/2.0):
            qc_flag = 1
        else:
            [t_stat, p_value] = stats.ttest_1samp(np.diff(frac_response),0)
            if(p_value < 0.05):
                # Good data
                qc_flag = 2
            elif((p_value >= 0.05) & (p_value < 0.32)):
                # Satisfactory data
                qc_flag = 3
            elif((p_value >= 0.32) & (p_value < 0.62)):
                # Poor data, few levels
                qc_flag = 4
            else:
                # No trend
                qc_flag = 1
    # Perform final check based on correlation
    if(qc_flag!=0 and qc_flag!=1):
        data_corr = stats.pearsonr(np.log10(dose_response['dose']+1e-15), frac_response)
        if(data_corr[0] < 0):
            qc_flag = 5

    return qc_flag

# Reformat dose-response data to be compatible with BMD analysis

def reformat_dose_response(dose_response):
    test_dose_response = pd.DataFrame(columns = ['dose', 'num_affected', 'total_num'])
    test_dose_response['dose'] = dose_response['dose']
    test_dose_response['num_affected'] = dose_response['num_affect']
    test_dose_response['total_num'] = dose_response['num_embryos']
    #index = np.arange(0,len(test_dose_response.dose))
    #test_dose_response.reset_index()
    test_dose_response.reset_index(inplace = True, drop = True) 
    return test_dose_response

