#!/usr/bin/env python
# coding: utf-8

import os, sys, time
import argparse

import re
from ingest import pull_raw_data, test_connection

OUT_FOLDER = '/tmp'
IF_EXITS = 'replace' # options: "append", "replace", "fail"
DB = 'develop' # options: "develop", "production"

parser = argparse.ArgumentParser('Take the files and moved to a database')

parser.add_argument('--validate', dest='validate', \
                    help='If this tag is added, then we validate existing files',\
                    action='store_true', default=False)

parser.add_argument('--update-db', dest='update_db', action='store_true', \
                    help='Include --update-db if you want to update the database',\
                    default=False)


def main():
    """
    main method for command line
    """
    start_time = time.time()
    args = parser.parse_args()

    allfiles = ['/tmp/'+a for a in os.listdir('/tmp') if 'csv' in a]
    print(allfiles)
    if args.validate:
        print("Validating existing files for database ingest")
        ##get files
        for fval in allfiles:
            valid.verify(pd.read_csv(fval, quotechar='"', quoting=1), re.sub('.csv', '', os.path.basename(fval)))
        ##validate
    if args.update_db:
        print('Saving to {}...'.format(DB))
        pull_raw_data(folder=OUT_FOLDER, if_exists=IF_EXITS, database=DB)
        print('Finished saving to database.')
  #  else: # if not saving to database, check connection to DB is okay
  #      print("Testing connection to database...", end='')
  #      okay, error = test_connection(database=DB)
  #      if okay:
  #          print('Connection OK')
  #      else:
  #          print('Connection failed, {}'.format(error))
    end_time = time.time()
    time_took = str(round((end_time-start_time), 1)) + " seconds"
    print ("Done, it took:" + str(time_took))


if __name__ == "__main__":
    main()
