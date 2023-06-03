#!/bin/bash

# Create a Python virtual environment
python3 -m venv venv

# Activate the virtual environment
source venv/bin/activate

# Update pip
pip install --upgrade pip

# Install requirements from a file
pip install -r requirements.txt

# Installing Aptfile
sudo apt-get install $(grep -vE "^\s*#" Aptfile  | tr "\n" " ") -y

echo "Virtual environment has been set up."
