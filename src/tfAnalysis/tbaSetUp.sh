#This script makes folder structure for running tba from fastq to TF list

#Create data folders needed for running the script

mkdir -p ../resources/data/bed
mkdir -p ../resources/data/fastq
mkdir -p ../resources/data/fastq/index
mkdir -p ../resources/data/HOMER
mkdir -p ../resources/data/IDR

#Create a python virtual environment and get the right software versions
python3 -m venv tbaVENV

source tbaVENV/bin/activate

pip install argparse
pip install multiprocess
pip install numpy==1.19.5
pip install pandas==0.25.3
pip install scipy==1.5.4
pip install sklearn==0.0
pip install biopython==1.79
pip install bio==1.0.2
pip install certifi==2020.12.5
pip install chardet==4.0.0
pip install cycler==0.10.0
pip install distance==0.1.3
pip install idna==2.10
pip install joblib==1.0.0
pip install kiwisolver==1.3.1
pip install matplotlib==2.2.2
pip install matplotlib_venn==0.11.6
pip install Pillow==8.0.1
pip install pyparsing==2.4.7
pip install python_dateutil==2.8.1
pip install pytz==2020.4
pip install requests==2.25.1
pip install scikit_learn==0.23.2
pip install seaborn==0.11.0
pip install setuptools==28.8.0
pip install six==1.15.0
pip install threadpoolctl==2.1.0
pip install urllib3==1.26.2

deactivate
