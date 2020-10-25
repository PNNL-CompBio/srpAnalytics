just make/enter any conda environment
for Doonam, it is conda activate tox where some dependencies are installed like statsmodels

then, simply python (v3) Main_Wrapper_BMD_Analysis.py



hard-coded part


line 26
complete_file_path = '/Users/kimd999/Dropbox/script/python/toxicology/DN_try/Phase_I_II.csv'


line ~80
#end_points = ['ANY24','ANY120','TOT_MORT','ANY_MORT','BRN_','CRAN','EDEM','LTKR','MUSC','SKIN','TCHR']
#end_points = ['AXIS','NC__','MO24','DP24','SM24','MORT']
end_points = ['ANY24']
#for chemical_id in np.unique(morphological_data['chemical.id']):
for chemical_id in [53, 54]:


