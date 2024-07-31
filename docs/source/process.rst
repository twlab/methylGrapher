

Process
====================


Genome Indexing
--------------------

Description
~~~~~~~~~~~~~~~~~~~~~~
This process involves transforming the genome graph from GFA format into two fully converted genome graphs: one depleted of C bases and another depleted of G bases. Additionally, if desired, you may include a spike-in genome in FASTA format to estimate the conversion rate in a single step further.

.. note::
    It's important to note that in the C-to-T genome graph, if both C and T segments are positioned identically—meaning they share the same parent and child segments, and each has only one parent and one child—the T segments are removed. Additionally, any associated links and paths are redirected to the corresponding C segment. This principle also applies in the G-to-A genome graph.

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








PrepareLibrary
--------------------

Description
~~~~~~~~~~~~~~~~~~~~~~

It converts your input (FASTQs) into fully G->A and C->T converted FASTQ file.
For single-end reads, just provide FASTQ file path to -fq1 argument.

Example Usage
~~~~~~~~~~~~~~~~~~~~~~
.. code-block:: shell

    methylGrapher PrepareLibrary -work_dir TheWorkingDir -fq1 Fastq1FilePath -fq2 Fastq2FilePath


Required arguments
~~~~~~~~~~~~~~~~~~~~~~
-work_dir <working_directory>

-fq1 <fastq1_file_path>



Optional arguments
~~~~~~~~~~~~~~~~~~~~~~
-fq2 <fastq2_file_path>

-directional (default: Y)






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

