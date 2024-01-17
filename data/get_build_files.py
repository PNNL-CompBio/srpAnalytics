'''
goal of this file is to download raw data to local repository so that they can be used to build new db
'''

import os
import pandas as pd
import wget

def get_build_files(version='1'):
    '''
    This function gets files with raw data for a particular version. this need to be
    run through the zfbmd data. all morphology and behavior files will have those
    words in the title so they can be queried by the next step
    '''
    # read in file location
    tab = pd.read_csv('srp_build_files.csv')

    fdict=dict()
    
    for i,row in tab.iterrows():
        #print(i,row)
        fname = row['name']
        dtype = row['data_type']
        samp  = row['sample_type']
        loc = row['location']
        ver = row['version']
        if ver !=  version:
            next
        print(fname)
        new_fname = fname+'_'+dtype+'.csv'
        tmpfile=wget.download(loc)
        os.rename(tmpfile,new_fname)
        if dtype=='morphology':
            fdict[fname]={'morphology': {'location':new_fname,'version':ver,'sample':samp}}
        else:
            fdict[fname]={'behavior': {'location':new_fname,'version':ver,'sample':samp}}
    
    return fdict

def get_processed_files():
    '''
    this function gets previously processed data to be added to portal. these
    files will all have `bmd` `dose_response` or `fit` in the title to be used
    in the database build
    '''
    fdict = dict()
    return fdict()


def get_sample_files():
    '''
    this function downloads sample files as part of the database build
    '''
rawfiles = get_build_files()

##buidlv1
##run docker images for v1 build

##buildv2
#run docker images for v2 build

##buildv3
#run docker images for v3 build