from json import load
from jsonschema import Draft3Validator
import pandas as pd
import argparse

workdir='/srpAnalytics'

schemas = {
    'chemdoseResponseVals': load(open(workdir+'/schemas/chemdoseResponseVals.json')),
    'chemicalsByExtractSample': load(open(workdir+'/schemas/chemicalByExtractSample.json')),
    'chemSummaryStats': load(open(workdir+'/schemas/chemSummaryStats.json')),
    'chemXYcoords': load(open(workdir+'/schemas/chemXYcoords.json')),
    'envSampdoseResponseVals': load(open(workdir+'/schemas/envSampdoseResponseVals.json')),
    'envSampSummaryStats': load(open(workdir+'/schemas/envSampSummaryStats.json')),
    'envSampXYcoords': load(open(workdir+'/schemas/envSampXYcoords.json')),
}

proc_schemas = {
    'envSample': load(open(workdir+'/schemas/envSampleIntake.json')),
    'zebrafish': load(open(workdir+'/schemas/zebrafishDataIntake.json'))
}

def verify(df, table_name):
    """ Takes a Pandas DataFrame and a Microsoft SQL Server table name to compare the number of rows in each. This function only works if the to_sql function replaces the table (if it appends the number of rows will obviously be off)

    Parameters:
        csv_df (Pandas Dataframe): dataframe of the original data from CSV
        table_name (string): name of table that you want to compare against

    Returns:
        nothing
    """
    print("Verifying schema...")
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

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Validate CSVs are in proper schema')
    parser.add_argument('csv_file', type=str, help='Path to CSV you want to validate')
    parser.add_argument('schema', type=str, help='The schema you want to validate against. Acceptable options are: {}'.format(list(schemas.keys())))

    args = parser.parse_args()

    if args.schema not in schemas.keys():
        print("{} is not a valid schema, please enter one of the following: {}".format(args.schema, list(schemas.keys())))
        exit(-1)

    if verify(pd.read_csv(args.csv_file, quotechar='"', quoting=1), args.schema):
        print('CSV file is in expected format.')
    else:
        print('Schema validation failed.')
