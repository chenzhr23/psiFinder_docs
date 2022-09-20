Ψ-sites SNPs and SNVs Searching
==================================

.. contents::
    :local:

psiFinder ``Ψ-sites SNPs and SNVs Searching`` was developed for searching Ψ-sites around SNPs and SNVs.

.. image:: /images/SNPSNVsearch.png

Input
------
Users should choose to upload files (i.e. rtsSeeker result) in bed format, SNPs/SNVs reference annotation in bed format, and genome fasta file, to search SNPs/SNVs around Ψ-sites.

Predict Ψ-sites miRNA target
---------------------------------------------
Once click ``START``, psiFindeer will run ``SNPSNVsearch.sh``. `bedtools <https://buildmedia.readthedocs.org/media/pdf/bedtools/latest/bedtools.pdf>`_ will firstly extend the sequence of the input sites to 21 nt length long (with ``_win21.bed`` suffix), then run ``bedtools intersect`` based on the extended sequence file and the input SNPs/SNVs annotation file.

``SNPSNVsearch.sh``

.. code:: bash

    #! /bin/bash

    if [[ $# -ne 3 ]]
    then
        echo 'Usage: ./'$0 ' bedfile snpsnvbed genome '
        exit 1
    fi

    bedfile=$1
    snpsnvbed=$2
    genome=$3


    awk 'FS=OFS="\t" {print $1,$2-10,$3+10,""$1"_"$2"_"$3"_"$6"",$5,$6}' ${bedfile} >${bedfile%.bed}_win21.bed

    echo -e "searching $snpsnvbed"
    bedtools intersect -a ${bedfile%.bed}_win21.bed -b $snpsnvbed -wa -wb >${bedfile%.bed}_SN_result.bed

    echo -e "Finished: SNPSNVsearch done!\n"
    echo -e "SNPSNVsearch result in $(dirname ${bedfile})"


Output
--------

Result with ``_SN_result.bed`` suffix is the final SNPs/SNVs searching result.

.. code:: bash

    $ cd /the/directory/of/out_file_dir
    $ tree -L 1
    .
    ├── Day0_common_rep1.bed
    ├── Day0_common_rep1_SN_result.bed
    └── Day0_common_rep1_win21.bed

    0 directories, 3 files

.. note:: All user input will be recorded in a plain text file with suffix ``_SNPSNVsearch_config.txt`` in psiFinder/config and help users to easily reload the previous config (by simply clicking ``CONFIG`` button).
