import pandas as pd
import numpy as np
import sqlalchemy as db
from jsonschema import validate
import os
import sys
from json import load

import logging
logging.basicConfig()
logging.getLogger("sqlalchemy.engine").setLevel(logging.ERROR)
logging.getLogger("sqlalchemy.pool").setLevel(logging.ERROR)

username=os.environ.get('DB_USER')
password=os.environ.get('DB_PASS')
db_dev_ip=os.environ.get('DB_DEV_IP')
db_prod_ip=os.environ.get('DB_PROD_IP')
db_port=os.environ.get('DB_PORT')
db_name=os.environ.get('DB_NAME')

schemas = {
    'chemdoseResponseVals': load(open('./schemas/chemdoseResponseVals.json')),
    'chemicalByExtractSample': load(open('./schemas/chemicalByExtractSample.json')),
    'chemSummaryStats': load(open('./schemas/chemSummaryStats.json')),
    'chemXYcoords': load(open('./schemas/chemXYcoords.json')),
    'envSampdoseResponseVals': load(open('./schemas/envSampdoseResponseVals.json')),
    'envSampSummaryStats': load(open('./schemas/envSampSummaryStats.json')),
    'envSampXYcoords': load(open('./schemas/envSampXYcoords.json')),
}

def test_connection(database):
    if database == "production":
        engine = db.create_engine('mssql+pyodbc://{}:{}@{}:{}/{}?driver=ODBC Driver 17 for SQL Server'.format(username, password, db_prod_ip, db_port, db_name), 
                        pool_pre_ping=True, echo=False, hide_parameters=True, connect_args={'connect_timeout': 100}, fast_executemany=True)
    else:
        engine = db.create_engine('mssql+pyodbc://{}:{}@{}:{}/{}?driver=ODBC Driver 17 for SQL Server'.format(username, password, db_dev_ip, db_port, db_name), 
                        pool_pre_ping=True, echo=False, hide_parameters=True, connect_args={'connect_timeout': 100}, fast_executemany=True)
    try:
       engine.connect()
       return True, ""
    except Exception as e:
        return False, e


def pull_raw_data(folder, if_exists, database):
    """ Looks for CSV files in the provided path and then calls the read_and_save function. This function calls all of the others, this is the management function.

    Parameters: 
        folder (string): Path to folder where CSVs are stored
        if_exists (string): Tells Pandas what to do if table already exists
            Reference: https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.to_sql.html 

    Returns:
        nothing
    """
    if database == "production":
        engine = db.create_engine('mssql+pyodbc://{}:{}@{}:{}/{}?driver=ODBC Driver 17 for SQL Server'.format(username, password, db_prod_ip, db_port, db_name), 
                        pool_pre_ping=True, echo=False, hide_parameters=True, connect_args={'connect_timeout': 100}, fast_executemany=True)
    else:
        engine = db.create_engine('mssql+pyodbc://{}:{}@{}:{}/{}?driver=ODBC Driver 17 for SQL Server'.format(username, password, db_dev_ip, db_port, db_name), 
                        pool_pre_ping=True, echo=False, hide_parameters=True, connect_args={'connect_timeout': 100}, fast_executemany=True)
    for file in os.listdir(folder):
        filename, file_extension = os.path.splitext(file)
        if (file_extension == ".csv"):
            print("Analyzing {}".format(filename))
            read_and_save(os.path.join(folder, file), filename, if_exists, engine)
            print("Finished analyzing {}".format(filename))


def read_and_save(csv_file, table_name, if_exists, engine):
    """ Takes a CSV file, reads (using Pandas), and then saves as a Microsoft SQL Server table

    Parameters: 
        csv_file (string): Path to CSV file
        table_name (string): The name of the table the data will be stored as. Defaults to the same name as the filename.
        if_exists (string): Tells Pandas what to do if table already exists
            Reference: https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.to_sql.html 

    Returns:
        nothing
    """
    print("\tReading csv...")
    df = pd.read_csv(csv_file, sep=',')
    print("\t\tFinished reading csv.")
    print("\t\tWriting to {}...".format(table_name))
    # If any infinite value is found, replace with a NULL value.
    df.replace([np.inf, -np.inf], np.nan, inplace=True)
    try:
        validate(instance=df.to_json(), schema=schemas['table_name'])
    except Exception as e:
        print("Invalid schema: {}".format(e))
        return
    # TODO uncomment once done testing schema
    #df.to_sql(table_name, engine, if_exists=if_exists, index=False)
    print("\tFinished writing to {}.".format(table_name))
    verify(df, table_name, engine)


def verify(csv_df, table_name, engine):
    """ Takes a Pandas DataFrame and a Microsoft SQL Server table name to compare the number of rows in each. This function only works if the to_sql function replaces the table (if it appends the number of rows will obviously be off)

    Parameters: 
        csv_df (Pandas Dataframe): dataframe of the original data from CSV
        table_name (string): name of table that you want to compare against

    Returns:
        nothing
    """
    print("\tVerifying saved data...")
    original_numRows = len(csv_df)
    saved_numRows = len(pull(table_name, engine))
    if original_numRows < saved_numRows:
        print('\t\tExtra rows saved, if this is unexpected please review')
    elif saved_numRows < original_numRows: 
        print('\t\tMissing rows, if this is unexpected please review')
    else:
        print('\t\tNo known issues with number of rows')
    print("\tFinished verifying.")


def pull(table_name, engine):
    """ Pulls all data from a table

    Parameters: 
        table_name (string): The name of the table the data will be stored as. Defaults to the same name as the filename.

    Returns:
        Pandas DataFrame containing all of the rows in that table
    """
    print("\t\t\tPulling from database table: {}".format(table_name))
    table_name = table_name
    connection = engine.connect()
    metadata = db.MetaData()
    table = db.Table(table_name, metadata, autoload=True, autoload_with=engine)
    query = db.select([table])
    result = connection.execute(query)
    print("\t\t\tReturning Pandas dataframe of results from database table.")
    return pd.DataFrame(result.fetchall())


if __name__ == "__main__":
    # Default Values
    if_exists = 'append'
    folder = './out'
    database = 'dev'

    print('Usage: python3 ingest.py data_folder if_exists=["append*", "replace", "fail"] database=["development*", "production"]')
    if len(sys.argv) >= 2: # includes flag for if_exists
        folder = sys.argv[1]
        print("Overwriting default folder value with: {}".format(folder))
    if len(sys.argv) >= 3:
        if_exists = sys.argv[2]
        print("Overwriting default if_exists value with: {}".format(if_exists))
    if len(sys.argv) == 4:
        database = sys.argv[3]
        print("Overwriting default database configuration value with: {}".format(database))

    pull_raw_data(folder=folder, if_exists=if_exists, database=database)