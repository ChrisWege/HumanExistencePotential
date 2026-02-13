#!/bin/bash

VENV_DIR="venv"

#Updating package list and installing system dependencies
sudo apt update
sudo apt install -y python3 python3-venv 

#Creating Python virtual environment
python3 -m venv "$VENV_DIR"

#Activating virtual environment
source "$VENV_DIR/bin/activate"

#Upgrading pip
python -m pip install --upgrade pip

#Installing Python packages from requirements file
python -m pip install -r pyvenv_list.txt

#Installing basemap for plotting
python -m pip install basemap

#Deactivating virtual environment
deactivate

echo "Setup completed successfully."
