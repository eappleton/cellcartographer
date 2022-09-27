#This script makes folder structure for running tba from fastq to TF list

#Create data folders needed for running the script
mkdir -p ../resources/data/bed
mkdir -p ../resources/data/fastq
mkdir -p ../resources/data/HOMER
mkdir -p ../resources/data/IDR

#Create a python virtual environment and get the right software versions
python3 -m venv tbaVENV
source tbaVENV/bin/activate
pip install argparse
pip install mutiprocessing
pip install pandas==0.20.3
pip install numpy==1.13.3
pip install scipy==0.19.0
pip install sklearn==0.19.1
pip install biopython==1.7
deactivate
