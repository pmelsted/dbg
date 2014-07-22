dbg
===

Naive de Bruijn graph implementation in python. This is the supporting code for the blog post [http://pmelsted.wordpress.com/2013/11/23/naive-python-implementation-of-a-de-bruijn-graph/](http://pmelsted.wordpress.com/2013/11/23/naive-python-implementation-of-a-de-bruijn-graph/).

A sample fastq file is provided and the code should be run with

`python dbg.py 31 read_1.fq ...`

where 31 is the k-mer size and read_1.fq ... the input FASTQ files.

The output is a [GFA file](http://lh3.github.io/2014/07/19/a-proposal-of-the-grapical-fragment-assembly-format/)

Any GFA file can be converted into a dot file using

`python gfa2dot.py [file]`

which will read the file (or stdin if no file is given) and produce a dot file.
