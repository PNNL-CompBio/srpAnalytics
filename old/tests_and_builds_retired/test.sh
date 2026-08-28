#!/bin/bash


catcmd() {
	   cat $1 | tee -a $2
   }


cat new_fits1.csv | tee -a new_fits.csv
#echo $catcmdf
#$catcmdf

#cat new_dose1.csv | tee -a new_dose.csv
#echo $catcmdd
#$catcmdd

#cat new_bmds1.csv >> new_bmds.csv
#echo $catcmdd
##$catcmdd
