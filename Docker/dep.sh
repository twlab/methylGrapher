#!/bin/bash


# Trim Galore install?
# FastQC install?

pip install methylGrapher



tar -xvzf FastUniq-1.1.tar.gz -C ./
cd ./FastUniq/source/
make
cp fastuniq /bin/


cd /bin
wget https://github.com/vgteam/vg/releases/download/v1.55.0/vg


chmod +x vg
chmod +x fastuniq

