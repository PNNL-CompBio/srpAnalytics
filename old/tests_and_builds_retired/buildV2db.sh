#!/bin/bash

##make a folder if one does not exist
if test -f temp; then
    echo "temp directory exists"
else
    mkdir temp
fi

#################################################### phase I,II data
##here we have the data from doo nam/lisa, we need to copy it into the temp directory
#phase12 data
p12_morph='raw_files/zf_morphology_data_335_chemicals_2020DEC16_fixed.csv'
p12_lpr_1='344_zf_LPR_data_phase_1_2_2020JUNE25_updated_plate_id_for_TX_tall_fixed_merged.csv'
#this second file was just a duplicate with data from chemicals that only had 15 time points or so
#p12_lpr_2='344_zf_LPR_data_phase_1_2_2020JUNE25_updated_plate_id_for_TX_tall_fixed_merged_full_240_timepoints.csv'
orig_lpr='raw_files/p12_lpr.tar.gz'


cmd='tar -xvzf '$orig_lpr
echo $cmd
$cmd

cp $p12_morph temp/morph0.csv
cp $p12_lpr_1 temp/lpr0_1.csv
#cp $p12_lpr_2 temp/lpr0_2.csv

#################################################### phase III data
##phase 3 data is on dropbox, we need to pull it to temp directory
#https://www.dropbox.com/sh/zg0q6wl13a3uo99/AAA0cdAK_fJwkJqpvF_HH6DWa?dl=0/
p3_morph='curl https://www.dropbox.com/sh/zg0q6wl13a3uo99/AACFZprOKkbvydjfoDI3oZo-a/Tanguay%20Phase%203%20zf%20morphology%20data%20PNNL%202021MAR23.csv -O temp/morph1.csv'
p3_lpr='wget https://www.dropbox.com/sh/zg0q6wl13a3uo99/AADd1QBStMguW9qYgzH2eatJa/Tanguay%20Phase%203%20zf%20LPR%20data%20PNNL%202021MAR23.csv -O temp/lpr1.csv'

echo $p3_morph
$p3_morph
echo $p3_lpr
$p3_lpr

#################################################### PFAS  data
##PFAS data is also on dropbox, we need to pull it
##downloading files from Lisa's dropbox PFAS data
#epr_cmd="wget https://www.dropbox.com/sh/69ootcq7yyvvx2h/AABgnmHtboM4LevxK1yxPIK-a/zf%20EPA%20PFAS%20EPR_PNNL_05-28-2021.csv" -O temp/epr.csv
#lpr_cmd='wget https://www.dropbox.com/sh/69ootcq7yyvvx2h/AABgzjaRPteU1EZIhnW9zv2Ka/zf%20EPA%20PFAS%20LPR_PNNL_05-28-2021.csv -O temp/lpr2.csv'
#mor_cmd='wget https://www.dropbox.com/sh/69ootcq7yyvvx2h/AABxsOLgwlv7-_HTZ0xaAIlNa/zf%20EPA%20PFAS%20morphology_PNNL_05-28-2021.csv -O temp/morph2.csv'

lpr_cmd="wget https://www.dropbox.com/s/ha7jduok03j82mf/zf%20EPA%20PFAS%20LPR_PNNL_05-28-2021.csv?dl=0 -O temp/lpr2.csv"
echo $lpr_cmd
$lpr_cmd

mor_cmd="wget https://www.dropbox.com/s/jma3b9al3u8hcny/zf%20EPA%20PFAS%20morphology_PNNL_05-28-2021.csv?dl=0 -O temp/morph2.csv"
echo $mor_cmd
$mor_cmd

################################################## run pipeline
#we have 3 pairs of files to run

dpath='/tmp/' ##path to files in docker images

all_lpr=$dpath"temp/lpr0_1.csv "$dpath"temp/lpr1.csv "$dpath"temp/lpr2.csv"
all_morph=$dpath"temp/morph0.csv "$dpath"temp/morph1.csv "$dpath"temp/morph2.csv"

##first we run validation on each
#docker pull sgosline/srp-schemadb
#docker run -v $PWD:/tmp sgosline/schemadb $all_lpr

##then we get the gene data
gpull="docker pull sgosline/srp-exposome"
echo $gpull
$gpull

grun="docker run -v "$PWD":/tmp sgosline/srp-exposome"
echo $grun
$grun

##then we run morph
dpull="docker pull sgosline/srp-zfbmd"
echo $dpull
$dpull

drun="docker run -v "$PWD":/tmp sgosline/srp-zfbmd --output=/tmp --morpho "$all_morph
echo $drun
$drun

##now rename these files
cpcmdb='mv new_bmds.csv new_bmds1.csv'
echo $cpcmdb
$cpcmdb

cpcmdf='mv new_fits.csv new_fits1.csv'
echo $cpcmdf
$cpcmdf

cpcmdd='mv new_dose.csv new_dose1.csv'
echo $cpcmdd
$cpcmdd

##then we concatentate them and run lpr

drun="docker run -v "$PWD":/tmp sgosline/srp-zfbmd --output /tmp --morpho "$all_morph" --LPR "$all_lpr
echo $drun
$drun

echo "BMDs\n"
cc='wc -l new_bmds.csv'
$cc

echo "Dose response\n"
wc -l new_dose.csv
echo "New Fits"
wc -l new_fits.csv
catcmd() {
   cat $1 >> $2
}


catcmdf='catcmd new_fits1.csv new_fits.csv'
echo $catcmdf
$catcmdf

catcmdd='catcmd new_dose1.csv new_dose.csv'
echo $catcmdd
$catcmdd

catcmdb='catcmd new_bmds1.csv new_bmds.csv'
echo $catcmdb
$catcmdb

##then we use output to build database
dpull="docker pull sgosline/bmd2samps:latest" # -t sgosline/srp-bmd2samps:latest"
echo $dpull
$dpull

##now build the database files
drun="docker run -v"$PWD":/tmp sgosline/srp-bmd2samps:latest --chemicals=/tmp/new_bmds.csv,/tmp/new_fits.csv,/tmp/new_dose.csv"
echo $drun
$drun
##then validate again and add to db

trun='tar -cvzf  srpCompendiumV2.tar.gz sigGeneStats.csv chemicals.csv samples.csv sampleToChemicals.csv zebrafishSampBMDs.csv zebrafishChemBMDs.csv zebrafishSampXYCoords.csv zebrafishChemXYCoords.csv zebrafishChemDoseResponse.csv  zebrafishSampDoseResponse.csv'
echo $trun
$trun
