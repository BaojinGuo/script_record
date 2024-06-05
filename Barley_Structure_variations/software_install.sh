###install sniffles2.3.3 by conda
conda create -n sniffle -c bioconda -c conda-forge sniffles
###install Svisionpro
## Get the source code
git clone https://github.com/songbowang125/SVision-pro.git
cd SVision-pro

## Create a conda environment for SVision-pro
conda env create -n svisionpro -f ./environment.yml 

## Install from source
conda activate svisionpro
python setup.py install

Note: Please make sure you have installed the same version of dependencies as in ./environment.yml. 
      We recommand that you create a new conda env and install by the command lines above. 
