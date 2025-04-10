

Process
====================


Genome Indexing
--------------------

Description
~~~~~~~~~~~~~~~~~~~~~~
This process involves transforming the genome graph from GFA format into two fully converted genome graphs: one depleted of C bases and another depleted of G bases. Additionally, if desired, you may include a spike-in genome in FASTA format to estimate the conversion rate in a single step further.

.. warning::
    Please use the same VG version for both genome indexing and alignment. Use different versions may cause VG to crash.


Example Usage
~~~~~~~~~~~~~~~~~~~~~~
.. code-block:: shell

    methylGrapher PrepareGenome -gfa YourGFAFilePath -lp LambdaPhageFastaFilePath -prefix OutputPrefix -t ThreadsUsed

Required arguments
~~~~~~~~~~~~~~~~~~~~~~
-gfa <gfa_file_path>
-prefix <output_prefix>

Optional arguments
~~~~~~~~~~~~~~~~~~~~~~
-lp <lambda_phage_genome_path>
-t <threads_used> (default: 1)


|
|
|
|
|


Deduplication
--------------------
Refer to FastUniq documentation for more details.

|
|
|
|
|




Main
--------------------

Description
~~~~~~~~~~~~~~~~~~~~~~

The core operation of methylGrapher comprises four essential steps. Firstly, fully converted FastQs are generated. Subsequently, these converted FastQs are aligned against methylGrapher indexes. Next, the resulting alignment files, formatted in GAF, undergo sorting to facilitate subsequent methylation extraction. Lastly, the methylation extraction process is executed.

Example Usage
~~~~~~~~~~~~~~~~~~~~~~
.. code-block:: shell

    methylGrapher Main -t ThreadsUsed -work_dir TheWorkingDir -index_prefix PrepareGenomeOutputPrefix


Required arguments
~~~~~~~~~~~~~~~~~~~~~~
-work_dir <working_directory>

-index_prefix <index_prefix>

-fq1 <fastq1_file_path>



Optional arguments
~~~~~~~~~~~~~~~~~~~~~~
-fq2 <fastq2_file_path>

-directional (default: Y)

-minimum_identity <minimum_identity> (default: 20)

-minimum_mapq <minimum_mapq> (default: 0)

-discard_multimapped (default: Y)

-batch_size <batch_size> (default: 4096)

-t <threads_used> (default: 1)



|
|
|
|
|












|
|
|
|
|




Align
--------------------

Description
~~~~~~~~~~~~~~~~~~~~~~

VG Giraffe alignment, please provide work directory and index prefix.

Example Usage
~~~~~~~~~~~~~~~~~~~~~~
.. code-block:: shell

    methylGrapher Align -t ThreadsUsed -work_dir TheWorkingDir -index_prefix PrepareGenomeOutputPrefix


Required arguments
~~~~~~~~~~~~~~~~~~~~~~
-work_dir <working_directory>

-index_prefix <index_prefix>




Optional arguments
~~~~~~~~~~~~~~~~~~~~~~

-t <threads_used> (default: 1)




|
|
|
|
|




MethylCall
--------------------

Description
~~~~~~~~~~~~~~~~~~~~~~

Methylation call from vg giraffe alignment result.

Example Usage
~~~~~~~~~~~~~~~~~~~~~~
.. code-block:: shell

    methylGrapher MethylCall  -t ThreadsUsed -work_dir TheWorkingDir -index_prefix PrepareGenomeOutputPrefix


Required arguments
~~~~~~~~~~~~~~~~~~~~~~
-work_dir <working_directory>

-index_prefix <index_prefix>




Optional arguments
~~~~~~~~~~~~~~~~~~~~~~

-minimum_identity <minimum_identity> (default: 20)

-minimum_mapq <minimum_mapq> (default: 0)

-discard_multimapped (default: Y)

-batch_size <batch_size> (default: 4096)

-t <threads_used> (default: 1)








|
|
|
|
|






ConversionRate
--------------------

Description
~~~~~~~~~~~~~~~~~~~~~~

Estimate the conversion rate from the spike-in genome.

Example Usage
~~~~~~~~~~~~~~~~~~~~~~
.. code-block:: shell

    methylGrapher ConversionRate -work_dir TheWorkingDir -index_prefix PrepareGenomeOutputPrefix


Required arguments
~~~~~~~~~~~~~~~~~~~~~~
-work_dir <working_directory>

-index_prefix <index_prefix>





|
|
|
|
|



MergeCpG
--------------------

Description
~~~~~~~~~~~~~~~~~~~~~~

Merge cytosine methylation call (graph.methyl) into CpG methylation call. During graph indexing, all CpG locations are identified and stored in a separate TSV file.
The graph CpG locations are stored under {index_prefix}cpg.tsv, with CpG id and both cytosine location on graph coordinate.
MergeCpG function will merge the cytosine methylation call (graph.methyl) into CpG methylation call using graph CpG id.

Example Usage
~~~~~~~~~~~~~~~~~~~~~~
.. code-block:: shell

    methylGrapher MergeCpG -work_dir TheWorkingDir -index_prefix PrepareGenomeOutputPrefix


Required arguments
~~~~~~~~~~~~~~~~~~~~~~
-work_dir <working_directory>

-index_prefix <index_prefix>





|
|
|
|
|







