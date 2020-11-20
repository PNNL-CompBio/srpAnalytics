import pandas as pd
from astropy import stats as astrostats

dose_response_data = {'dose': [0.0, 1.0, 5.0], 'num_affected': [0.0, 3.0, 8.0], 'total_num': [10, 10, 10]}
test_dose_response = pd.DataFrame(dose_response_data, columns = ['dose', 'num_affected', 'total_num'])
print(test_dose_response)

for index in range(len(test_dose_response.dose)):
    CI = astrostats.binom_conf_interval(test_dose_response.num_affected[index], test_dose_response.total_num[index], confidence_level = 0.95, interval='jeffreys')
    print ("\n\nCI with jeffreys:" + str(CI))
