Ψ-sites miRNA Target Prediction
==================================

.. contents::
    :local:

psiFinder ``Ψ-sites miRNA Target Prediction`` was developed for Ψ-sites miRNA target prediction (e.g. Ψ-modified mRNA regulated by miRNA).

.. image:: /images/MIscan.png

Input
------
Users should choose to upload files (i.e. rtsSeeker result) in bed format, miRNA sequence fasta file, and genome fasta file, to predict miRNA-(Ψ)target interaction.

miranda core algorithm parameters
*****************************************
``score/energy/scaling/strict`` options of miranda correspond to ``Set score threshold``, ``Set energy threshold``, ``Set scaling parameter``, and ``Demand strict 5' seed pairing``.

Predict Ψ-sites miRNA target
---------------------------------------------
Once click ``START``, psiFindeer will run ``MIscan.sh``. `bedtools <https://buildmedia.readthedocs.org/media/pdf/bedtools/latest/bedtools.pdf>`_ will firstly extend the sequence of the input sites to 21 nt length long (with ``_win21.bed`` suffix), then `miranda <https://cbio.mskcc.org/miRNA2003/miranda.html>`_ will be used to predict the miRNA target based on the input miRNA sequence file and the extended sequence file.

``MIscan.sh``

.. code:: bash

    #! /bin/bash

    if [[ $# -ne 7 ]]
    then
        echo 'Usage: ./'$0 ' bedfile miRNAseq genome score energy scaling strict'
        exit 1
    fi

    bedfile=$1
    miRNAseq=$2 #grep "Homo sapiens" mature.fa -A 1 | grep -v "\--" > only_human_mirna.fa or download directly from miBase
    genome=$3
    score=$4
    energy=$5
    scaling=$6
    strict_bool=$7
    strict=''

    if [[ "${strict_bool}" =~ "true"  ]] ;
    then
        "$strict" = "-strict"
    fi


    awk 'FS=OFS="\t" {print $1,$2-10,$3+10,""$1"_"$2"_"$3"_"$6"",$5,$6}' ${bedfile} >${bedfile%.bed}_win21.bed
    bedtools getfasta -fi $genome -bed ${bedfile%.bed}_win21.bed -name+ -s>${bedfile%.bed}_win21.fa

    echo "finding miRNA target"
    miranda $miRNAseq ${bedfile%.bed}_win21.fa -sc $score -en $energy -scale $scaling $strict > ${bedfile%.bed}_miRNAtarget.out


    grep 'Seq1,Seq2,Tot Score,Tot Energy,Max Score,Max Energy,Strand,Len1,Len2,Positions' -A 1 ${bedfile%.bed}_miRNAtarget.out >${bedfile%.bed}_miRNAtarget_tmp.txt
    sed -i -e 's/>//' -e 's/--//' -e 's/,/\t/g' ${bedfile%.bed}_miRNAtarget_tmp.txt
    sort -r -u ${bedfile%.bed}_miRNAtarget_tmp.txt >${bedfile%.bed}_miRNAtarget.txt

    echo -e "Finished: MIscan done!\n"
    echo -e "MIscan result in $(dirname ${bedfile})"


Output
--------

Result with ``_miRNAtarget.txt`` suffix is the final miRNA target prediction result.

.. code:: bash

    $ cd /the/directory/of/out_file_dir
    $ tree -L 1
    .
    ├── Day0_common_rep1.bed
    ├── Day0_common_rep1_miRNAtarget.out
    ├── Day0_common_rep1_miRNAtarget_tmp.txt
    ├── Day0_common_rep1_miRNAtarget.txt
    ├── Day0_common_rep1_win21.bed
    └── Day0_common_rep1_win21.fa

    0 directories, 6 files



.. note:: All user input will be recorded in a plain text file with suffix ``_MIscan_config.txt`` in psiFinder/config and help users to easily reload the previous config (by simply clicking ``CONFIG`` button).
