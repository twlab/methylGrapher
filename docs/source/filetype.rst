

File Types
====================================


GFA: Graphical Fragment Assembly format
------------------------------------------------------------
For detail explanation, refer to https://gfa-spec.github.io/GFA-spec/GFA1.html


GAF: Graphical mApping Format
------------------------------------------------------------
For detail explanation, refer to https://github.com/lh3/gfatools/blob/master/doc/rGFA.md


methyl: Graph methylation file
------------------------------------------------------------
The methylGrapher extraction output provides a list of each cytosine mapped to the genome graph, specified by its coordinates (segment ID, 0-based position on the segment, and strand relative to the segment).
For each cytosine, the output includes its context, the number of read-pairs supporting whether the cytosine is methylated or unmethylated, and the coverage, which is the sum of methylated and unmethylated read-pairs.
The methylation percentage is calculated as the ratio of methylated read-pairs to the total coverage.


Reference:
------------------------------
1. `The design and construction of reference pangenome graphs with minigraph <https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02168-z>`_
2. `Glossary of vg-related File Types <https://github.com/vgteam/vg/wiki/File-Types>`_





