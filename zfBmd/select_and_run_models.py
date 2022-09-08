#!/usr/bin/env python
# coding: utf-8

import numpy as np
import pandas as pd
import scipy.stats as stats

import BMD_Analysis_Functions as baf

def model_fitting(dose_response, BMD_Flags):
    """
    Description
    
    Returns
    ----

    
    Parameters
    ----

    """

    ###########################################
    ## CALCULATE VALUES FOR LOW QUALITY DATA ##
    ###########################################

    low_quality = dose_response[dose_response["ids"].isin(BMD_Flags[BMD_Flags["flag"].isin([0,1])]["ids"])].groupby("ids")

    # Calculate low quality metrics and start new BMDS file 
    BMDS_LowQual = low_quality.apply(lambda df: np.trapz(df["frac.affected"], x = df["conc"])).reset_index().rename(columns = {0: "AUC"})
    BMDS_LowQual[["Model", "BMD10", "BMDL", "BMD50"]] = np.nan
    BMDS_LowQual["Min_Dose"] = round(low_quality[["ids", "conc"]].min("conc").reset_index()["conc"], 4)
    BMDS_LowQual["Max_Dose"] = round(low_quality[["ids", "conc"]].max("conc").reset_index()["conc"], 4)
    BMDS_LowQual["AUC_Norm"] = BMDS_LowQual["AUC"] / (BMDS_LowQual["Max_Dose"] - BMDS_LowQual["Min_Dose"])

    # Order columns
    BMDS_LowQual = BMDS_LowQual[["ids", "Model", "BMD10", "BMDL", "BMD50", "AUC", "Min_Dose", "Max_Dose", "AUC_Norm"]]

    ################
    ## RUN MODELS ##
    ################

    # Define a function to run models
    def select_and_run_models(ID):
        '''For each ID, run the 8 regression models and return the best one.'''
        
        Data = dose_response[dose_response["ids"] == ID]

        # Get the number of nonNA samples per dose
        NonNATotals = Data["num.nonna"].tolist() 

        # Regression model function
        def run_regression_model(modelfun, fittedfun):
            '''Fit the regression model and return the parameters, fitted_values, and the p_value'''

            # Run the model
            model = modelfun(Data[["conc", "num.affected", "num.nonna"]].astype('float').copy())

            # Get the model parameters
            model_params = model.fit().params

            # Get the model's fitted values
            model_fittedvals = fittedfun(Data["conc"], model_params)

            # Get the p_value
            model_pval = calc_p_value(model_fittedvals, model_params)

            # Get the AIC
            AIC = -2*model.fit().llf + (2 * len(model_params))

            # Return a list
            return([model, model_params, model_fittedvals, model_pval, AIC])

        # Calculate P-Value Function
        def calc_p_value(PredictedValues, Params):
            '''Return a p-value of model fit for each unique ID and Model dataframe pairing'''

            # Get the experimental values 
            ExperimentalValues = Data["frac.affected"].tolist()

            # Now, calculate the chi squared value
            ChiSquared = ((NonNATotals / (PredictedValues * (1 - PredictedValues))) * (ExperimentalValues - PredictedValues)**2).sum()

            # Calculate a p-value of fit 
            return(1 - stats.chi2.sf(ChiSquared, len(NonNATotals) - len(Params)))

        # Run regression models in a dictionary
        models = {

            ## Logistic ##
            "Logistic": run_regression_model(baf.Logistic, baf.logistic_fun),

            ## Gamma ## 
            "Gamma": run_regression_model(baf.Gamma, baf.gamma_fun),

            ## Weibull ##
            "Weibull": run_regression_model(baf.Weibull, baf.weibull_fun),

            ## Log-Logistic ##
            "Log Logistic": run_regression_model(baf.Log_Logistic, baf.log_logistic_fun),

            ## Probit ##
            "Probit": run_regression_model(baf.Probit, baf.probit_fun),

            ## Log-Probit ##
            "Log Probit": run_regression_model(baf.Log_Probit, baf.log_probit_fun),

            ## Multistage ##
            "Multistage": run_regression_model(baf.Multistage_2, baf.multistage_2_fun),

            ## Quantal Linear ##
            "Quantal Linear": run_regression_model(baf.Quantal_Linear, baf.quantal_linear_fun),

        }

        # Iterate through all p-values 
        p_values = {}
        for key in models.keys():
            p_values[key] = models[key][3]

        # Iterate through all AICs 
        aics = {}
        for key in models.keys():
            aics[key] = models[key][4]

        # Determine the best model
        BestModel = min(aics, key=lambda k: aics[k]) 

        # Return results 
        return(
            [
                p_values,
                models[BestModel],
                BestModel,
                aics
            ]
        )

    # Run the models 
    model_results = {}
    for id in BMD_Flags[BMD_Flags["flag"].isin([2,3,4,5])]["ids"]:
        model_results[id] = select_and_run_models(id)

    # Make p-value dataframe
    p_value_list = []
    for id in model_results.keys():
        theDict = model_results[id][0]
        theDict["ids"] = id
        p_value_list.append(theDict)

    p_value_df = pd.DataFrame(p_value_list)

    ########################################
    ## CALCULATE FLAGS FOR MODEL MATCHING ##
    ########################################

    # Get flags of whether p_values are all below 0.1
    Analysis_Flags = p_value_df
    Analysis_Flags.iloc[:,:8] = Analysis_Flags.iloc[:,:8] < 0.1
    Analysis_Flags = Analysis_Flags.groupby("ids").apply(lambda df: int(df.all(axis = 1))).reset_index().rename(columns = {0:"BMD_Analysis_Flag"})
    BMD_Flags = pd.merge(BMD_Flags, Analysis_Flags)

    #######################################
    ## PULL AKAIKE INFORMATION CRITERION ##
    #######################################

    # Make aic dataframe
    aic_list = []
    for id in model_results.keys():
        theDict = model_results[id][3]
        theDict["ids"] = id
        aic_list.append(theDict)

    aic_df = pd.DataFrame(aic_list)

    ##############################
    ## CALCULATE BENCHMARK DOSE ##
    ##############################

    def Calculate_BMD(Model, params, BenchmarkResponse = 0.1):
        '''Calculate a benchmark dose'''

        # For each model, extract the parameters and run the calculations
        if (Model == "Logistic"): 
            alpha_, beta_ = params
            return(np.log((1 + np.exp(-alpha_)*BenchmarkResponse)/(1-BenchmarkResponse))/beta_)
        elif (Model == "Gamma"):
            g_, alpha_, beta_ = params
            return(stats.gamma.ppf(BenchmarkResponse, alpha_, scale = 1/beta_))
        elif (Model == "Weibull"):
            g_, alpha_, beta_ = params
            return((-np.log(1 - BenchmarkResponse)/beta_)**(1/alpha_))
        elif (Model == "Log Logistic"):
            g_, alpha_, beta_ = params
            return(np.exp((np.log(BenchmarkResponse/(1-BenchmarkResponse)) - alpha_)/beta_))
        elif (Model == "Probit"):
            alpha_, beta_ = params
            p_0 = stats.norm.cdf(alpha_)
            p_BMD = p_0 + (1 - p_0)*BenchmarkResponse
            return((stats.norm.ppf(p_BMD) - alpha_)/beta_)
        elif (Model == "Log Probit"):
            g_, alpha_, beta_ = params
            return(np.exp((stats.norm.ppf(BenchmarkResponse) - alpha_)/beta_))
        elif (Model == "Multistage"):
            g_, beta_, beta2_ = params
            return((-beta_ + np.sqrt((beta_**2) - (4*beta2_*np.log(1 - BenchmarkResponse))))/(2*beta2_))
        elif (Model == "Quantal Linear"):
            g_, beta_ = params
            return(-np.log(1 - BenchmarkResponse)/beta_)
        else:
            print(Model, "was not recognized as an acceptable model choice.")

    def Calculate_BMDL(Model, FittedModelObj, Data, BMD10, params, MaxIterations = 100, ToleranceThreshold = 1e-4):
        '''Run BMDL Function Test'''

        # Reformat data
        Data = Data[["conc", "num.affected", "num.nonna"]].astype('float').copy()

        # Define an initial low and high threshold
        BMD_Low = BMD10/10
        BMD_High = BMD10

        # Start a counter and set tolerance to 1
        Iteration_Count = 0
        Tolerance = 1

        # Set a LLV Threhold
        BMDL_LLV_Thresh = FittedModelObj.fit().llf - stats.chi2.ppf(0.9, 1)/2

        # Start a while condition loop
        while ((Tolerance > ToleranceThreshold) and (Iteration_Count <= MaxIterations)):
            
            # Add to the iteration counters 
            Iteration_Count+=1

            # If maximum iterations are reached, set BMDL to NA and break
            if (Iteration_Count == MaxIterations):
                BMDL = np.nan
                break

            # BMDL should be the mid value between the low and high estimate
            BMDL = (BMD_Low + BMD_High)/2
            ModelObj = np.nan

            # Select the correct BMD model
            if (Model == "Logistic"):
                ModelObj = baf.Logistic_BMD(Data).profile_ll_fit([params[0], BMDL]) # Value is alpha
            elif (Model == "Gamma"):
                ModelObj = baf.Gamma_BMD(Data).profile_ll_fit([params[0], params[1], BMDL]) # Value is g and alpha
            elif (Model == "Weibull"):
                ModelObj = baf.Weibull_BMD(Data).profile_ll_fit([params[0], params[1], BMDL]) # Value is g and alpha
            elif (Model == "Log Logistic"):
                ModelObj = baf.Log_Logistic_BMD(Data).profile_ll_fit([params[0], params[2], BMDL]) # Value is g and beta
            elif (Model == "Probit"):
                ModelObj = baf.Probit_BMD(Data).profile_ll_fit([params[0], BMDL]) # Value is alpha
            elif (Model == "Log Probit"):
                ModelObj = baf.Log_Probit_BMD(Data).profile_ll_fit([params[0], params[1], BMDL]) # Value is g and alpha
            elif (Model == "Multistage"):
                ModelObj = baf.Multistage_2_BMD(Data).profile_ll_fit([params[0], params[1], BMDL]) # Value is g and beta 1
            elif (Model == "Quantal Linear"):
                ModelObj = baf.Quantal_Linear_BMD(Data).profile_ll_fit([params[0], BMDL]) # Value is g
            else:
                print(Model, "was not recognized as an acceptable model choice.")

            # Pull the llf 
            LLF = ModelObj.llf

            # If the calculated LLF is not within the threshold, set high to BMDL and run again 
            if((LLF - BMDL_LLV_Thresh) > 0):
                BMD_High = BMDL
            # Otherwise, set low to BMDL
            else:
                BMD_Low = BMDL

            Tolerance = abs(LLF - BMDL_LLV_Thresh)

        return(BMDL)


    # Build BMD table for fitted data 
    BMDS_Model = []

    for id in model_results.keys():

        # Get the model name 
        Model = model_results[id][2]

        # Get the fitted modeol object
        FittedModelObj = model_results[id][1][0]

        # Get the parameters 
        params = model_results[id][1][1]

        # Get the BMD10 value 
        BMD10 = Calculate_BMD(Model, params, 0.1)

        # Get the dose response data
        Data = dose_response[dose_response["ids"] == id]

        # Get the AUC, min, and max dose 
        AUC = np.trapz(Data["frac.affected"], x = Data["conc"])
        Min_Dose = round(min(Data["conc"]), 4)
        Max_Dose = round(max(Data["conc"]), 4)

        # Return results in a dictionary
        rowDict = {
            "ids": id,
            "Model": Model,
            "BMD10": BMD10, 
            "BMDL": Calculate_BMDL(Model, FittedModelObj, Data, BMD10, params),
            "BMD50": Calculate_BMD(Model, params, 0.5),
            "AUC": AUC,
            "Min_Dose": Min_Dose,
            "Max_Dose": Max_Dose,
            "AUC_Norm": AUC / (Max_Dose - Min_Dose)
        }
        BMDS_Model.append(rowDict)

    BMDS_Model = pd.DataFrame(BMDS_Model)
    return(BMDS_Model)


