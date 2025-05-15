#!/usr/bin/env bash
# install.sh — auto-install Python requirement

set -e

# 1. check Python
if ! command -v python3 &> /dev/null; then
  echo "Error: python3 not installed"
  exit 1
fi

# 2. check pip
if ! python3 -m pip --version &> /dev/null; then
  echo "pip not intalled,try installing pip..."
  curl https://bootstrap.pypa.io/get-pip.py -o get-pip.py
  python3 get-pip.py
  rm get-pip.py
fi

# 3. install requirements
if [ -f requirements.txt ]; then
  echo "Installing requirements..."
  python3 -m pip install --user -r requirements.txt
else
  echo "No requirements.txt, skipp installing requirements."
fi

