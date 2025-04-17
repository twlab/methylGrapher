

Installation
====================================


Prerequisites
------------------------------
Adapter Trimming: `Trim Galore! <https://github.com/FelixKrueger/TrimGalore>`_

Deduplication: `FastUniq <https://sourceforge.net/projects/fastuniq/files/>`_

Alignment: `VG Giraffe <https://github.com/vgteam/vg>`_

Python3: Tested on 3.9, 3.10, 3.11 and 3.12.
I would recommend 3.12 as it is the latest version and has the best performance.


VG: Tested from 1.49 to 1.63.1

1.61.0 is recommended at the moment: from 1.49 to 1.61 all seems fine, however the short read mapping speed is around 3 times slower than the previous versions.

.. warning::
    You must ensure that both the VG toolkit and Python 3 are installed and added to your PATH.


Install via PIP
----------------------------------------

.. code-block:: shell

    pip install methylGrapher


Install via CONDA
----------------------------------------

.. code-block:: shell

    conda create -n methylGrapher python=3.12
    conda activate methylGrapher
    # Not recommended, but if you want to use conda to install vg, you can do so. Download VG from their release page instead.
    conda install bioconda::vg
    pip install methylGrapher


Use Docker
----------------------------------------
`Docker Hub <https://hub.docker.com/repository/docker/wenjin27/methylgrapher/general>`_

Source code
----------------------------------------
`The source code <https://github.com/twlab/MethylGrapher>`_

