docker build . -t srp-analytics
docker run -v $PWD:/tmp srp-analytics to_be_processed/7_PAH_zf_morphology_data_2020NOV11_tall.csv,to_be_processed/7_PAH_zf_morphology_data_2021_Jan11.csv
