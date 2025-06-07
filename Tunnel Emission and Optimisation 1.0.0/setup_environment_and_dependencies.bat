@echo off
REM Delete old virtual environment if it exists
rd /s /q .venv

REM Create a new virtual environment
python -m venv .venv

REM Activate the virtual environment
call .venv\Scripts\activate

REM Upgrade pip and install required packages
pip install --upgrade pip
pip install PyQt6==6.8.1 numpy sympy matplotlib pandas scipy pyomo pymoo scikit-learn alphashape descartes psutil

REM Export installed packages to requirements.txt
pip freeze > requirements.txt

echo Virtual environment setup complete.
pause
