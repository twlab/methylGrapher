#!/bin/zsh



# Push to Github



# Tag the release



# New build for PyPi

# TODO: You need manually update the version in 2 places
# main.py and utility.py
python3 setup.py sdist
twine upload dist/*


# Build the docker image
# TODO: Update the version in build.sh
cd Docker
bash build.sh





