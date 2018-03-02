PYTHON="anaconda3-5.0.0"
pyenv install ${PYTHON}
pyenv local ${PYTHON}
pip install snakemake
pyenv rehash
