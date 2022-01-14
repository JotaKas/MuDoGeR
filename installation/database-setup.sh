#!/bin/bash

# this scripts asks for the desired location of installation for the databases.
# given the user input, the config file will be edited and all databases via wget, curl, etc and zipping

database_location="$1"
config_file=MuDoGeR/installation/config #Path problem! The script does not know where this config is
touch "$config_file"

mkdir "$database_location"

source "$config_file"
############################################### PROKARYOTES ###############################################
### CheckM
echo 'installing checkm database ...'
mkdir -p "$database_location"/checkm
cd  "$database_location"/checkm
if [ ! -f selected_marker_sets.tsv ]; then
wget https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz
tar -xvf checkm_data_2015_01_16.tar.gz
rm -fr checkm_data_2015_01_16.tar.gz
# On newer versions of CheckM, you would run:
#checkm data setRoot /path/to/your/dir/$MY_CHECKM_FOLDER
CHECKM_DB="$database_location"/checkm #Fixed? we need to test
echo CHECKM_DB="$CHECKM_DB" >> "$config_path"
else echo "-> your CheckM database is ready"
fi


### GTDB-tk
mkdir -p  "$database_location"/"gtdbtk"
cd "$database_location"/"gtdbtk"
if [ ! -d release* ]; then
wget https://data.gtdb.ecogenomic.org/releases/latest/auxillary_files/gtdbtk_data.tar.gz
tar xvzf gtdbtk_data.tar.gz
rm -fr gtdbtk_data.tar.gz
#echo  GTDBTK_DATA_PATH="$database_location"/gtdbtk/gtdbtk_r95_data >> "$config_file" ##FIX HERE
echo GTDBTK_DATA_PATH="$(ls "$database_location"/gtdbtk/release*)" >>  "$config_file" # fixed? we need to test
else echo "-> your GTDBtk database is ready"
fi

############################################### VIRUSES ###############################################
### CheckV
mkdir -p  "$database_location"/checkv
cd "$database_location"/checkv
if [ ! -d checkv-db-v1.0 ]; then
wget https://portal.nersc.gov/CheckV/checkv-db-v1.0.tar.gz
tar -zxvf checkv-db-v1.0.tar.gz
rm -fr checkv-db-v1.0.tar.gz
#ADD CHECKV DATABASE PATH TO CONFIG FILE
CHECKVDB="$database_location"/checkv/checkv-db-v1.0
echo CHECKVDB="$CHECKVDB" >> "$config_file"
else echo "-> your CheckV database is ready"
fi

############################################### EUKARYOTES ###############################################




