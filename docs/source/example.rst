

Toy Example
====================================

This is a toy example to show how to use the methylGrapher package to analyze the methylation data.
Also you can use it to verify the installation of the package.

All the script and data are in the `example` folder below.
https://github.com/twlab/methylGrapher/blob/main/example/example_bs




.. code-block:: shell

    git clone git@github.com:twlab/methylGrapher.git
    cd methylGrapher/example/example_bs/
    bash pipeline.sh


The make_example.py script, which is incorporated within the pipeline.sh workflow, will generate the genome graph in GFA format along with the simulated FASTQ files in the example_bs folder.


The final output will be methylGrapher/example/graph.methyl


