from json import load
from jsonschema import Draft3Validator
import pandas as pd

schemas = {
    'chemdoseResponseVals': load(open('./schemas/chemdoseResponseVals.json')),
    'chemicalsByExtractSample': load(open('./schemas/chemicalByExtractSample.json')),
    'chemSummaryStats': load(open('./schemas/chemSummaryStats.json')),
    'chemXYcoords': load(open('./schemas/chemXYcoords.json')),
    'envSampdoseResponseVals': load(open('./schemas/envSampdoseResponseVals.json')),
    'envSampSummaryStats': load(open('./schemas/envSampSummaryStats.json')),
    'envSampXYcoords': load(open('./schemas/envSampXYcoords.json')),
}

def verify(df, table_name):
    """ Takes a Pandas DataFrame and a Microsoft SQL Server table name to compare the number of rows in each. This function only works if the to_sql function replaces the table (if it appends the number of rows will obviously be off)

    Parameters: 
        csv_df (Pandas Dataframe): dataframe of the original data from CSV
        table_name (string): name of table that you want to compare against

    Returns:
        nothing
    """
    print("\tVerifying schema...")
    df = df.where(pd.notnull(df), None)
    v = Draft3Validator(schemas[table_name])
    errors = set()
    for row in df.to_dict(orient='records'):
        for error in sorted(v.iter_errors(row), key=str):
            errors.add(str(error))   

    if errors:
        print('Validation errors when running schema check on {}'.format(table_name))
        with open("/tmp/{}_validation_errors.txt".format(table_name), 'w+') as fp:
            for error in errors:
                fp.write("{}\n\n\n".format(error))
        return False
    return True