'''
goal of this file is to download raw data to local repository so that they can be used to build new db
'''

import os
import pandas as pd
import wget
import argparse


def get_build_files(tmpdir,version=1.0):
    '''
    This function gets files with raw data for a particular version. this need to be
    run through the zfbmd data. all morphology and behavior files will have those
    words in the title so they can be queried by the next step
    '''
    # read in file location
    tab = pd.read_csv('srp_build_files.csv')

    ##create dictionary for each type of file
    fdict = {'morphology':[],'behavior':[],'bmd':[],'fit':[],'dose':[],'sample':[],'expression':[],'reference':[],'other':[]}
    
    for i,row in tab.iterrows():
        #print(i,row)
        fname = row['name']
        dtype = row['data_type']
        samp  = row['sample_type']
        loc = row['location']
        ver = row['version']
        if ver !=  version or pd.isna(loc):
            continue
        print('\n'+fname+' '+dtype)
        new_fname = tmpdir+'/'+fname+'_'+dtype+'.csv'
        try:
            tmpfile=wget.download(loc)
        except:
            print("Cannot get file: "+loc)
            continue
        os.rename(tmpfile,new_fname)
        fdict[dtype].append({'fname':fname,'location':new_fname,'version':ver,'sample':samp})

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



def main():
    parser=argparse.ArgumentParser('Initiate srp analytics build')
    parser.add_options('--tmpdir',dest='tmpdir',default='temp',help='location to store files')
    args=parser.parse_args()
    
    tmpdir=opts.tmpdir
    os.mkdir(tmpdir)
    rawfiles = get_build_files(tmpdir)

    ##run BMD calculation on morophology files
    morph_files=['/tmp/'+a['location'] for a in rawfiles['morphology']]
    for val in rawfiles.items():
        print (val)
    bmdpull='docker pull sgosline/srp-zfbmd'
    bmdcmd='docker run -v $PWD:/tmp sgosline/srp-zfbmd'
    os.system(bmdcmd+' --morpho '+' '.join(morph_files))
    
    ##run BMD calculations on behavior files


    ##run zfbmd on fitted files and sample files
    
    
