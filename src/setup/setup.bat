@echo off
REM Define the Conda environment name
SET ENV_NAME=nbody

REM init conda if it was not done before
CALL conda init

REM Create a new Conda environment with Python
conda create -y -n %ENV_NAME% python=3.10 pip

REM Deactivate first, because sometimes it bugs (for me, at least)
CALL conda deactivate

REM Activate the environment
CALL conda activate %ENV_NAME%

REM Install dependencies with pip inside the Conda environment
pip install -r requirements.txt

REM Confirm success
echo Conda environment "%ENV_NAME%" set up successfully with dependencies.