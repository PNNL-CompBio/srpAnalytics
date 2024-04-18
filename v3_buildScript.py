'''
Build script moved to python for better extendability and interoperability.

'''


import os
import subprocess
import pandas as pd



def collectFiles(data_dir='https://raw.githubusercontent.com/PNNL-CompBio/srpAnalytics/main/data',filename='srp_build_files.csv'):
    '''
    every time the build file is updated, this script will collect the files and return
    a dictionary of files to be fed into each module
    '''
    df = pd.read_csv(data_dir+'/'+filename)
    return df


def fitCurveFiles(morpho_behavior_tuples):
    '''
    get new curve fits, list of tuples of morpho/behavior pairs
    '''
    


def combineFiles(location_list,ftype):
    '''
    helper function to combine duplicates
    '''
    dflist=[]
    required_columns = {'bmd':['Chemical_ID','End_Point','Model','BMD10','BMD50',"Min_Dose","Max_Dose",\
                                "AUC_Norm","DataQC_Flag","BMD_Analysis_Flag"],#,"BMD10_Flag","BMD50_Flag{"),
                          'dose':['Chemical_ID',"End_Point","Dose","Response","CI_Lo","CI_Hi"],\
                          'fit':['Chemical_ID',"End_Point","X_vals","Y_vals"]}

    print('concatenating '+ftype)
    for loc in location_list.location:
        f = pd.read_csv(loc)[required_columns[ftype]]
        dflist.append(f)
    fulldf=pd.concat(dflist)
    fulldf = fulldf.drop_duplicates()

    
    return fulldf.drop_duplicates()

def main():
    df = collectFiles()

    ##first find the morphology and behavior pairs for chemical sources
    chemdf = df.loc[df.sample_type=='chemical']
    morph = chemdf.loc[chemdf.data_type=='morphology']
    beh = chemdf.loc[chemdf.data_type=='behavior']
    tupes =[]
    for n in morph.name:
        tupes.append([morph.loc[morph.name==n].location,beh.loc[beh.name==n].location])

    ##call bmdrc on all morphology/behavior pairs for sample sources
    newbmds,newfits,newdoses =[],[],[]
    ##get output in bmds, fits, and curves

    #add chemical BMDS, fits, curves to existing data
    chemfiles=[]
    sampfiles=[]
    for st in ['chemical','extract']:
        for dt in ['bmd','fit','dose']:
            fdf = combineFiles(df.loc[df.sample_type==st].loc[df.data_type==dt],dt)
            fname = 'tmp_'+st+'_'+dt+'.csv'
            fdf.to_csv(fname,index=False)
            if st=='chemical':
                chemfiles.append(fname)
            else:
                sampfiles.append(fname)

    ##now map sample information
    sid = list(df.loc[df.name=='sampId'].location)[0]
    cid = list(df.loc[df.name=='chemId'].location)[0]
    cclass = list(df.loc[df.name=='class1'].location)[0]
    emap = list(df.loc[df.name=='endpointMap'].location)[0]
    fses = ','.join(list(df.loc[df.data_type=='sample'].location))
    ctfile = list(df.loc[df.name=='compTox'].location)[0]
    desfile = list(df.loc[df.name=='chemdesc'].location)[0]
    smap = list(df.loc[df.name=='sampMap'].location)[0]
                
    ##call script with sample files
    cmd = "Rscript bmd2Samps_v3/buildV3database.R --sample --drcFiles="+','.join(sampfiles)+\
        ' --sampId='+smap+' --chemId='+cid+' --epMap='+emap+' --chemClass='+cclass+\
        ' --compToxFile='+ctfile+' --sampleFiles='+fses+' --chemDesc='+desfile+\
        ' --sampMap='+smap
    
    print(cmd)
    os.system(cmd)
    ##call script with chemical files
            
main()
