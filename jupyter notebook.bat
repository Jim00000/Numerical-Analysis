@echo off
@where /Q jupyter
if %ERRORLEVEL% NEQ 0 (
    @echo on
    @echo jupyter.exe wasn't found. Is jupyter installed? Download jupyter with $ python -m pip install --user jupyter
    pause
) else (
    jupyter notebook
)
