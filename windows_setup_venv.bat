@echo off

:: Create a Python virtual environment
python -m venv venv

:: Activate the virtual environment
call venv\Scripts\activate.bat

:: Update pip
python -m pip install --upgrade pip

:: Install requirements from a file
pip install -r requirements.txt

echo "Virtual environment has been set up."