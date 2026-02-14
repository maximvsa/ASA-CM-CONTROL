@echo off
conda env create -f environment.yml
conda activate asa-cm-control-env
idaes get-extensions --extra petsc
python -c "import idaes; print('IDAES OK', idaes.__version__)"
