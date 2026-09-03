# Benchmark Dose Curves

## Input Data

A **binary class** object was created using data in **long** format. The following column names were set:

|Parameter|Column Name|
|---------|-----------|
|Chemical|chemical_id|
|Plate|plate_id|
|Well|well|
|Concentration|concentration|
|Endpoint|endpoint|
|Value|value|

## Pre-Processing

#### **Combine & Make New Endpoints**

New endpoints were made using existing endpoints using 'or', which means that if there is any endpoints with a '1', this new endpoint will also have a '1', regardless of how many zeroes there are in the other endpoints. See a summary table of added endpoints below:

|New Endpoint Name|Combined Existing Endpoints|
|---|---|
|ANY24|MO24, DP24, SM24, NC24|
|ANY120|MORT, YSE_, AXIS, EYE_, SNOU, JAW_, OTIC, PE__, BRAI, SOMI, PFIN, CFIN, PIG_, CIRC, TRUN, SWIM, NC__, TR__, ANY24|
|TOT_MORT|MO24, MORT|
|ALL_BUT_MORT|DP24, SM24, NC24, YSE_, AXIS, EYE_, SNOU, JAW_, OTIC, PE__, BRAI, SOMI, PFIN, CFIN, PIG_, CIRC, TRUN, SWIM, NC__, TR__|
|BRN_|BRAI, OTIC, PFIN|
|CRAN|EYE_, SNOU, JAW_|
|EDEM|YSE_, PE__|
|LTRK|TRUN, CFIN|
|MUSC|CIRC, SWIM, SOMI|
|SKIN|PIG_|
|TCHR|TR__|

#### **Set Invalid Wells to NA**

In some cases, like when a sample fish dies, many affected endpoints need to be set to NA. Here, the 'Endpoint Name' column denotes the specific endpoint that sets this rule. In this example, it could be MORT for mortality. Then, the endpoint value needs to be set, which in this case would be a 1 to indicate sample fish that did die. All endpoints would then be set to NA except for cases where the endpoint should not be affected, which are referred to as 'Endpoint Exceptions.'

|Endpoint Name|Endpoint Value|Endpoint Exceptions|
|---|---|---|
|DNC_|1|None|
|MORT|1|DP24, MO24, SM24, MORT|
|MO24|1|MO24|

#### **Remove Invalid Endpoints**

The following endpoints were removed: DNC_, PIG_, TR__, AUC2, AUC3, AUC4, BRAI, CFIN, CIRC, DP24, EYE_, JAW_, MORT, MUSC, NC__, PE__, SKIN, SM24, TCHR, YSE_, OTIC, PFIN, PIG_, SNOU, SOMI, SWIM, TR__, TRUN

## Filtering

#### **Negative Control Filter**

Plates with unusually high responses in negative control samples were filtered. The response threshold was set to **50**. See a summary below:

|Response|Number of Plates|Filter|
|---|---|---|
|0.0|16594|Keep|
|3.125|56|Keep|
|3.2258|1|Keep|
|3.3333|1|Keep|
|3.8462|2|Keep|
|4.1667|249|Keep|
|5.0|1|Keep|
|6.25|3378|Keep|
|6.6667|128|Keep|
|6.8966|1|Keep|
|7.1429|958|Keep|
|7.6923|54|Keep|
|8.3333|1415|Keep|
|9.0909|43|Keep|
|9.375|47|Keep|
|10.0|22|Keep|
|11.1111|7|Keep|
|12.5|1615|Keep|
|13.3333|12|Keep|
|14.2857|245|Keep|
|15.3846|10|Keep|
|15.625|11|Keep|
|16.6667|382|Keep|
|17.5|1|Keep|
|18.1818|5|Keep|
|18.75|758|Keep|
|20.0|4|Keep|
|20.8333|37|Keep|
|21.4286|59|Keep|
|21.875|1|Keep|
|23.0769|1|Keep|
|25.0|474|Keep|
|26.6667|2|Keep|
|27.2727|1|Keep|
|28.5714|17|Keep|
|29.1667|7|Keep|
|30.7692|1|Keep|
|31.25|220|Keep|
|33.3333|28|Keep|
|35.7143|4|Keep|
|37.5|119|Keep|
|41.6667|10|Keep|
|42.8571|1|Keep|
|43.75|61|Keep|
|44.4444|1|Keep|
|45.8333|5|Keep|
|50.0|24|Remove|
|54.1667|2|Remove|
|56.25|10|Remove|
|58.3333|2|Remove|
|62.5|4|Remove|
|68.75|7|Remove|
|75.0|2|Remove|
|100.0|8|Remove|

And here is the plot:
![Filter Negative Control](./filter_negative_control.png)

#### **Minimum Concentration Filter**

Endpoints with too few concentration measurements (non-NA) to model are removed. The minimum was set to **3**. See a summary below:

|Number of Concentrations|Number of Endpoints|Filter|
|---|---|---|
|24|18|Keep|
|20|63|Keep|
|18|23|Keep|
|17|48|Keep|
|16|111|Keep|
|15|68|Keep|
|14|61|Keep|
|12|836|Keep|
|11|61|Keep|
|10|1770|Keep|
|9|441|Keep|
|8|53|Keep|
|7|689|Keep|
|5|5619|Keep|
|3|72|Keep|
|1|1404|Remove|

And here is the plot:
![Filter Minimum Concentration](./filter_minimum_concentration.png)

#### **Correlation Score Filter**

Endpoints with little to no positive correlation with dose are unexpected and should be removed. The correlation threshold was set to **0.2**. See a summary below:

|Correlation Score Bin|Number of Endpoints|
|---|---|
|-1.0|133.0|
|-0.8|368.0|
|-0.6|540.0|
|-0.4|910.0|
|-0.2|1086.0|
|0.0|2332.0|
|0.2|1287.0|
|0.4|1099.0|
|0.6|1147.0|
|0.8|1031.0|

And here is the plot:
![Filter Correlation Score](./filter_correlation_score.png)

## Model Fitting & Output Modules

#### **Filter Summary**

Overall, 11337 endpoint and chemical combinations were considered. 4564 were deemed eligible for modeling, and 6773 were not based on filtering selections explained in the previous section. Of the 4564 deemed eligible for modeling, 923 did not pass modeling checks.

#### **Model Fitting Selections**

The following model fitting parameters were selected.

|Parameter|Value|Parameter Description|
|---|---|---|
|Goodness of Fit Threshold|0.1|Minimum p-value for fitting a model. Default is 0.1|
|Akaike Information Criterion (AIC) Threshold|2|Any models with an AIC within this value are considered an equitable fit. Default is 2.
|Model Selection|lowest BMDL|Either return one model with the lowest BMDL, or combine equivalent fits|

#### **Model Quality Summary**

Below is a summary table of the number of endpoints with a high quality fit and those with poor fit, as defined by each label below.

| Modeled Flag                          |   Count |
|:--------------------------------------|--------:|
| Fail - correlation score filter       |    5369 |
| Pass                                  |    3641 |
| Fail - other filter                   |    1404 |
| Fail - GOF check                      |     916 |
| Fail - GOF check; Fail - other filter |       7 |

And here is a summary delineating the good and moderate fits,based off of the following properties.

| Flag | Number of Non-Control Concentrations | Spearman Correlation | Goodness of Fit | BMD50 | Model Convergence |
| -- | -- | -- | -- | -- | -- |
| Not Fit | < 3 | < 0.2 | < 0.1 | Not within concentration range | No Models Converged |
| Moderate | >= 3 | 0.2 - 0.7 | >= 0.1 | Not within concentration range | At least 1 model converged |
| Good | >= 5 | > 0.7 | >= 0.1 | Within concentration range | At least 1 model converged |

| DataQC Flag   |   Count |
|:--------------|--------:|
| Not Fit       |    7696 |
| Moderate      |    3083 |
| Good          |     558 |

#### **Output Modules**

Below, see a table of useful methods for extracting outputs from bmdrc.

|Method|Description|
|---|---|
|.bmds|Table of fitted benchmark dose values|
|.bmds_filtered|Table of filtered models not eligible for benchmark dose calculations|
|.output_res_benchmark_dose|Table of benchmark doses for all models, regardless of whether they were filtered or not|
|.p_value_df|Table of goodness of fit p-values for every eligible endpoint|
|.aic_df|Table of Akaike Information Criterion values for every eligible endpoint|
|.response_curve|Plot a benchmark dose curve for an endpoint|

