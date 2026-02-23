@echo off
conda env create -f environment.yaml
conda run -n asa-cm-control-env idaes get-extensions --extra petsc
conda run -n asa-cm-control-env python -c "import sys, idaes; print(sys.version); print('IDAES OK', idaes.__version__)"
