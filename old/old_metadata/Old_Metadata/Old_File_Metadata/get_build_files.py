'''
goal of this file is to download raw data to local repository so that they can be used to build new db
'''

import os
import pandas as pd
import wget
import argparse
import dropbox


def get_build_files(tmpdir,ftype='morphology',version=1.0):
    '''
    This function gets files with raw data for a particular version. this need to be
    run through the zfbmd data. all morphology and behavior files will have those
    words in the title so they can be queried by the next step
    '''
    # read in file location
    tab = pd.read_csv('srp_build_files.csv')

    ##create dictionary for each type of file
    fdict = {'morphology':[],'behavior':[],'bmd':[],'fit':[],'dose':[],'sample':[],'expression':[],'reference':[],'other':[]}
 #   dbx = dropbox.dropbox_client()

    for i,row in tab.iterrows():
        #print(i,row)
        fname = row['name']
        dtype = row['data_type']
        if dtype!=ftype:
            continue
        samp  = row['sample_type']
        loc = row['location']
        ver = row['version']
        if ver !=  version or pd.isna(loc):
            continue
        print('\n'+fname+' '+dtype)
        new_fname = tmpdir+'/'+fname+'_'+dtype+'.csv'
        try:
      #      if 'dropbox' in loc:
      #          dbx.files_download_to_file(new_fname, loc)
      #      else:
            tmpfile = wget.download(loc)
            os.rename(tmpfile,new_fname)
        except:
            print("Cannot get file: "+loc)
            continue
        fdict[dtype].append({'fname':fname,'location':new_fname,'version':ver,'sample':samp})

    return fdict



def main():
    parser=argparse.ArgumentParser('Initiate srp analytics build')
    parser.add_argument('--tmpdir',dest='tmpdir',default='temp',help='location to store files')
    parser.add_argument('--filetype',destp='ftype',defult='morphology',help='Type of file to collect. Can be one of: morphology, behavior, bmd, dose, fit, sample, expression')
    args=parser.parse_args()
    
    tmpdir=args.tmpdir
    if not os.path.exists(tmpdir):
        os.mkdir(tmpdir)
    rawfiles = get_build_files(tmpdir,opts.ftype)

    print(rawfiles)

main()
    
