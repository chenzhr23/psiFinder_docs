Quality Control
=================

.. role:: red

.. contents::
    :local:

psiFinder ``Quality Control`` utilized `cutadapt <https://cutadapt.readthedocs.io/en/stable/guide.html>`_ to control quality for input sequencing reads.
Users should choose to upload input files suitable for `cutadapt <https://cutadapt.readthedocs.io/en/stable/guide.html>`_ software. Once ``START``, psiFinder will run ``1_QC_SE.sh`` or ``1_QC_PE.sh`` to gain processed fastq files.


.. tip:: **single-end**: Upload sequencing reads in ``fastq|fq|fastq.gz|fq.gz`` format from :red:`single-end` reverse transcription stop data.

.. image:: /images/QC_SE.png


.. tip:: **paired-end**: Upload sequencing reads in ``fastq|fq|fastq.gz|fq.gz`` format from :red:`paired-end` reverse transcription stop data.

.. image:: /images/QC_PE.png

Input
************************************
Users can choose to upload files with one of the format suffixes ``fastq|fq|fastq.gz|fq.gz`` for both CMC-input and CMC-treated groups. The ``Quality Control`` QT widget accept the argument of command-line adapter trimming tools `cutadapt <https://cutadapt.readthedocs.io/en/stable/guide.html>`_.

Perform quality control
************************************
Simple clicking of ``START`` button after uploading the data and setting the adapter sequence, uploaded data will be automatically controlled for customized output.

Output
*************************
Results in fastq fromat (with ``.trimmed.fastq`` suffix) will be output to the target directory.

.. note:: All user input will be recorded in a plain text file with suffix ``_cutadapt_SE_config.txt`` in psiFinder/config and help users to easily reload the previous config (by simply clicking ``CONFIG`` button).

