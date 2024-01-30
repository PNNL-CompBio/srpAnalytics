#!/bin/bash
# V3 of the database includes the zebrafish data 

##make a folder if one does not exist
if test -f temp; then
    echo "temp directory exists"
else
    mkdir temp
fi

##Catcmd adds two files together
catcmd() {
   cat $1 >> $2
}


##build docker file for data
cmd='docker build -f data/Dockerfile -t srp-data data/ --build-arg HTTPS_PROXY=$HTTPS_PROXY temp/'
##run docker file for data
echo $cmd
$cmd


##run docker file for data
cmd='docker run srp-data'
echo $cmd
$cmd


##all data should be loaded into temp file

##build docker file for zfbmd
cmd='docker build -f zfbmd/Dockerfile -t srp-zfbmd zfbmd/ --build-arg HTTPS_PROXY=$HTTPS_PROXY'
##add in data path so it knows file files to run

echo $cmd
$cmd

##run docker file for zfbmd

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


##get the zf expression data
gpull='docker pull sogsline/srp-zfexp'
echo $gpull
$gpull

grun='docker run -v '$PWD':/tmp sgosline/srp-zfexp'
echo $grun
$grun

##then we have to concatenate the two together
catg='catcmd srpDEGstats.csv sigGeneStats.csv'
echo $catg
$catg

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
