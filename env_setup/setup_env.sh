#!/bin/bash

cdir=$(pwd)

CONDA=$(which conda)

if [[ -z "$CONDA" ]];
then
	echo "No conda was found please install before running"
	exit
fi


conda env create -f afenv.yml

CENV=$(conda info --envs | grep alphafold2_env)

if [[ -z "$CENV" ]];
then
	echo "there was an error installing the env please try again, if it continues to fail please reach out to george@acpharis.com"
	exit
fi 



conda activate alphafold2_env
mkdir databases

./download_alphafold_params.sh $cdir/databases

cd $CONDA_PREFIX/lib/python3.11/site-packages
git clone https://github.com/gj85369/alphafold_acpharis.git 
cp -r alphafold_acpharis/alphafold .
cd $cdir


RUNIT=$(ls ../runit)

if [[ -z "$RUNIT" ]]
then
	echo "the runit file was not found, please run this script in the dir it was sent it, or you can manually change the runit file"
	exit 
fi 

sed -i "s@PATHTOENVDIR@$CONDA_PREFIX@g" ../runit
sed -i "s@PATHTOALPHAFOLDDIR@$CONDA_PREFIX/lib/python3.11/site-packages/alphafold@g" ../runit
sed -i "s@PATHTOPARAMSDIR@$cdir/databases@g" ../runit


