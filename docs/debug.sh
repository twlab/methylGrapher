#!/bin/zsh


rm -rf build
make html
cd build/html


open http://localhost:8000
python -m http.server 8000



# https://cdn.jsdelivr.net/gh/twlab/methylGrapher/