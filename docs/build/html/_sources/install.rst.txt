

Installation
====================================


Prerequisites
------------------------------
Adapter Trimming: `Trim Galore! <https://github.com/FelixKrueger/TrimGalore>`_

Deduplication: `FastUniq <https://sourceforge.net/projects/fastuniq/files/>`_

Alignment: `VG Giraffe <https://github.com/vgteam/vg>`_

Python3: Tested on 3.9, 3.10 and 3.12.
I would recommend 3.12 as it is the latest version and has the best performance.

.. warning::
    You must ensure that both the VG toolkit and Python 3 are installed and added to your PATH.


Install via PIP
----------------------------------------

.. code-block:: shell

    pip install methylGrapher


Install via CONDA
----------------------------------------

.. code-block:: shell

    conda create -n methylGrapher python=3.9
    conda activate methylGrapher
    conda install bioconda::vg
    pip install methylGrapher


Use Docker
----------------------------------------
`Docker Hub <https://hub.docker.com/repository/docker/wenjin27/methylgrapher/general>`_

Source code
----------------------------------------
`The source code <https://github.com/twlab/MethylGrapher>`_

