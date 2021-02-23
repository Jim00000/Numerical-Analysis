@echo off
@where /Q jupyter
if %ERRORLEVEL% NEQ 0 (
    @echo on
    @echo jupyter.exe wasn't found. Is jupyter installed? Download jupyter with $ python -m pip install jupyter jupyterlab
    pause
) else (
    jupyter-lab
)
