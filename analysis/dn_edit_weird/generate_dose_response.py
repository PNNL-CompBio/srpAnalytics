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
def gen_dose_response(data_ep_cid, end_point):
    dose_response = pd.DataFrame(columns = ['dose', 'num_affect', 'frac_affect', 'num_embryos', 'tot_wells'])
    # Remove all wells for plates for which number of hits for negative controls > 50% wells
    for plate_id in np.unique(data_ep_cid['plate.id']):
        # Count number of wells corresponding to negative controls
        data_ep_cid_plate = data_ep_cid.loc[data_ep_cid['plate.id'] == plate_id]
        neg_ctrl_wells = data_ep_cid_plate.loc[data_ep_cid_plate['conc'] == 0]
        num_neg_ctrl_wells = neg_ctrl_wells.shape[0]
        num_neg_ctrl_hits = (neg_ctrl_wells[end_point]).sum(axis=0,skipna=True,min_count=1)
        if(num_neg_ctrl_hits > 0.5*num_neg_ctrl_wells):
            # Delete all wells corresponding to that plate
            data_ep_cid = data_ep_cid[data_ep_cid['plate.id'] != plate_id]
    
    for concentration_id in np.unique(data_ep_cid['conc']):
        data_ep_cid_concs = data_ep_cid.loc[(data_ep_cid['conc'] == concentration_id)]
        # Get total number of wells for a given concentration
        tot_wells = len(data_ep_cid_concs[end_point])
        num_nonnan_wells = sum(~np.isnan(data_ep_cid_concs[end_point]))
        num_affected = (data_ep_cid_concs[end_point]).sum(axis=0,skipna=True,min_count=1)
        if(num_nonnan_wells == 0):
            fraction_affected = np.nan
        else:
            fraction_affected = num_affected / num_nonnan_wells
        dose_response = dose_response.append({'dose': concentration_id, 'num_affect': num_affected , 'frac_affect': fraction_affected, 'num_embryos': num_nonnan_wells, 'tot_wells': tot_wells}, ignore_index = True)
        
    # Filter dose-responsee to delete those dose groups where number of wells are
    # less than 50% of total wells
    # First get rid of nan values
    dose_response = dose_response.dropna()
    delete_count = 0
    for dr_index in range(dose_response.shape[0]):
        dr_index_original = dr_index
        dr_index = dr_index - delete_count
        if((dose_response.iloc[dr_index].num_embryos) < (0.5*(dose_response.iloc[dr_index].tot_wells))):
            dose_response = dose_response[dose_response.index != dr_index_original]
            delete_count+=1
    
     
    return dose_response


# Get data QC code
def BMD_feasibility_analysis(dose_response):
    '''This function performs feasibility analysis
    for dose respone data. The value returned is a 
    flag indicating data quality as defined below:
    0: Not enough dose groups for BMD analysis. BMD analysis not performed
    1: No trend detected in dose-response data. BMD analysis not performed
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
            #elif((frac_response.iloc[-1] + frac_response.iloc[-2])/2.0 < (frac_response.iloc[0] + frac_response.iloc[1])/2.0):
                # No trend
            #    qc_flag = 1
            else:
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

