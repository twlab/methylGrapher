#!/bin/bash


# Trim Galore install?
# FastQC install?

pip install methylGrapher



tar -xvzf FastUniq-1.1.tar.gz -C ./
cd ./FastUniq/source/
make
cp fastuniq /bin/
mv GraphAligner /bin/


cd /bin
wget https://github.com/vgteam/vg/releases/download/v1.64.1/vg
mv vg vg_v1.64.1
wget https://github.com/vgteam/vg/releases/download/v1.61.0/vg


chmod +x vg
chmod +x vg_v1.64.1
chmod +x fastuniq

