#!/bin/zsh



# Push to Github



# Tag the release



# New build for PyPi
python3 setup.py sdist
twine upload dist/*


# Build the docker image
cd Docker
bash build.sh





