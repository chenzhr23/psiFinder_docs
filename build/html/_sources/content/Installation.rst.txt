Installation
=============

``Remember`` -- psiFinder are available for **Unix-based** system!

.. contents::
    :local:

Download psiFinder
--------------------
*psiFinder is free and open to all users.*

If you're here for the first time, download psiFinder first, and decompress it into a directory whatever you like.

`Download psiFinder v.0.1 for Linux x64 <https://mega.nz/fm/5CsUVbKC>`_

.. note:: For users who want to remote operation a linux system in windows, the SSH client must be capable of X11-Forwarding, such as `MobaXterm <https://mobaxterm.mobatek.net/>`_ et al.

Download genome sequence and annotation file
------------------------------------------------------------------

Users should download the genome file and annotation file for the species of interest.

.. note:: Current version of psiFinder mainly support genome data derived from human, more support for other species please contact us for new version development!

`Download hg38 genome <fa> <https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/>`_

`Download hg38 annotation <gtf> <https://www.gencodegenes.org/human/>`_

Download example data
-----------------------

psiFinder provide some example data (human) for the user to test. If want to do so, you should download the example data.

`Download example data <https://mega.nz/fm/wacTDQQK>`_


Install conda
---------------------------------
psiFinder needs to use ``conda`` to configure the environment, please install conda first:

* *Linux*

.. code:: bash

    wget -c https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
    chmod 777 Miniconda3-latest-Linux-x86_64.sh
    bash Miniconda3-latest-Linux-x86_64.sh
    source ~/.bashrc

Starting the software
---------------------------------

* *Linux*

Once the psiFinder mian program and the Genome annotation are downloaded and unzipped properly, go to the directory of psiFinder:

.. code:: bash

    $ gunzip psiFinder.zip
    $ cd /the/directory/of/unzipped_psiFinder
    $ tree -L 1
    .
    ├── config
    ├── lib
    ├── log
    ├── psiFinder
    ├── script
    └── snakemake

    5 directories, 2 files

The operation of psiFinder needs to grant executable permission to the script of the software. Please execute the following code before running the software:

.. code:: bash

    $ cd /the/directory/of/unzipped_psiFinder
    $ chmod 777 *

Then, upon executed the ``./psiFinder``, the following window appears:

.. code:: bash

	$ ./psiFinder

.. image:: /images/psiFinder_window.png

* *Windows WSL*

For Windows operate system that equipped with WSL (`Windows Subsystem for Linux <https://learn.microsoft.com/en-us/windows/wsl/install-manual>`_), users should firstly install `Xming <https://sourceforge.net/projects/xming/>`_ software in your personal computer. Each time before running ``./psiFinder``, launch Xming and run ``export DISPLAY=:0`` (to ensure and render display device working) in WSL terminal.

.. image:: /images/WSL.png

In some cases, when we run ``bedAnnotator`` (``Ψ-sites Annotation``) of psiFinder in WSL, the following error may occurred:

.. code:: bash

    $ ./bedAnnotator
    $ ./bedAnnotator: error while loading shared libraries: libssl.so.10: cannot open shared object file: No such file or directory

    $ ldd bedAnnotator
    linux-vdso.so.1 (0x00007fffcca6f000)
    libssl.so.10 => not found
    libcrypto.so.10 => not found
    libpng12.so.0 => not found
    libz.so.1 => /lib/x86_64-linux-gnu/libz.so.1 (0x00007f41b60c0000)
    libstdc++.so.6 => /lib/x86_64-linux-gnu/libstdc++.so.6 (0x00007f41b5ed0000)
    libm.so.6 => /lib/x86_64-linux-gnu/libm.so.6 (0x00007f41b5d81000)
    libgcc_s.so.1 => /lib/x86_64-linux-gnu/libgcc_s.so.1 (0x00007f41b5d50000)
    libpthread.so.0 => /lib/x86_64-linux-gnu/libpthread.so.0 (0x00007f41b5d2d000)
    libc.so.6 => /lib/x86_64-linux-gnu/libc.so.6 (0x00007f41b5b30000)
    /lib64/ld-linux-x86-64.so.2 (0x00007f41b60f3000)

This error occurred due to the loss of the ``libssl.so.10``, ``libcrypto.so.10``, and ``libpng12.so.0``, To avoid this error, we can copy these shared library to the ``/lib`` directory (Therefore, make sure the ``bedAnnotator`` can be run in your terminal, e.g. run ``./bedAnnotator`` to see if help information is showed).

.. code:: bash

    $ cd /the/directory/of/unzipped_psiFinder
    $ sudo cp script/lib/libssl.so.10 /lib
    $ sudo cp script/lib/libcrypto.so.10 /lib
    $ sudo cp script/lib/libpng12.so.0 /lib

Configuration
---------------------------
Before uploading data, you should first configure the operating environment of the software, this step can be achieved through clicking the option ``global->configuration`` at the upper of the main window.

.. image:: /images/configuration.png

Once click ``CHECK``, psiFindeer will run ``0_configuration.sh``.

``0_configuration.sh``

.. code:: bash

    #!/bin/bash

    Help() {
    echo
    "0_software_test.sh example:
     bash 0_software_test.sh
    "
    }


    usage() {                                      # Function: Print a help message.
      echo "Usage: $0 [ -h help]" 1>&2
    }
    exit_abnormal() {                              # Function: Exit with error.
      usage
      exit 1
    }


    #psiFinder requires pre-installation of mupdf, wkhtmltopdf, conda, perl, R (optparse, dplyr, tidyr, tidyverse, ggrepel, stringr, ggplot2, RColorBrewer, caTools, cowplot, pROC, e1071,
    # neuralnet, venndiagram, reshape2, limma, edgeR, gprofiler2, ggrepel, pathview, bedr, scales, motifStack), snakemake, seqtk, cutadapt, STAR, samtools, rtsSeeker, bedAnnotator, ACAscan,
    # miranda, homer, gtftogenepred, bedtools, minion, please click 'CHECK' button to check if it is already installed; if not installed, required softwares will be installed or set into environment
    # variables by conda.

    while getopts ":h:" options; do

      case "${options}" in
        h|:)
          usage
          Help
          exit 0
          ;;
        \?) # incorrect option
          echo "Error: -${OPTARG} Invalid option"
          exit_abnormal
          ;;
      esac
    done

    shift $(($OPTIND - 1))
    date  ## echo the date at start
    echo "Starting: software checking..."

    #install conda
    if [ -z  `which conda` ]
    then
    echo "conda is not installed, or not in environment variables. Configuration check will add it to environment variables..."
    eval "cd pipline_script/"
    conda_setup="bash conda_setup.sh"
    eval $conda_setup
    wait
    eval "cd ../"
    else
    echo "conda(miniconda3) is installed."
    fi

    #install perl
    if [ `which perl | grep -c "miniconda3/bin/perl"` == 0 ]
    then
    echo "perl is not installed, or not in environment variables. Configuration check will add it to environment variables..."
    conda_perl="conda install --yes -c anaconda perl"
    eval $conda_perl
    else
    echo "perl(miniconda3) is installed."
    fi

    #install R
    if [ `which R | grep -c "miniconda3/bin/R"` == 0 ]
    then
    echo "R is not installed, or not in environment variables. Configuration check will add it to environment variables..."
    conda_R="conda install --yes -c r r-base"
    eval $conda_R
    wait
    else
    echo "R(miniconda3) is installed."
    fi

    #install optparse R package
    echo "installing optparse..."
    conda install --yes -c bioconda r-optparse
    echo -e "\n"

    #install openxlsx R package
    echo "installing openxlsx..."
    conda install --yes -c bioconda r-openxlsx
    echo -e "\n"

    #install dplyr R package
    echo "installing dplyr..."
    conda install --yes -c r r-dplyr
    echo -e "\n"

    #install tidyr R package
    echo "installing tidyr..."
    conda install --yes -c r r-tidyr
    echo -e "\n"

    #install tidyverse R package
    echo "installing tidyverse..."
    conda install --yes -c r r-tidyverse
    echo -e "\n"

    #install stringr R package
    echo "installing stringr..."
    conda install --yes -c r r-stringr
    echo -e "\n"

    #install ggplot2 R package
    echo "installing ggplot2..."
    conda install --yes -c r r-ggplot2
    echo -e "\n"

    #install ggpol R package
    echo "installing ggpol..."
    conda install --yes -c conda-forge r-ggpol
    echo -e "\n"

    #install ggpubr R package
    echo "installing ggpubr..."
    conda install --yes -c conda-forge r-ggpubr
    echo -e "\n"

    #install RColorBrewer R package
    echo "installing RColorBrewer..."
    conda install --yes -c r r-rcolorbrewer
    echo -e "\n"

    #install cowplot R package
    echo "installing cowplot..."
    conda install --yes -c conda-forge r-cowplot
    echo -e "\n"

    #install gridextra R package
    echo "installing gridextra..."
    conda install --yes -c r r-gridextra
    echo -e "\n"

    #install pROC R package
    echo "installing pROC..."
    conda install --yes -c r r-proc
    echo -e "\n"

    #install e1071 R package
    echo "installing e1071..."
    conda install --yes -c r r-e1071
    echo -e "\n"

    #install neuralnet R package
    echo "installing neuralnet..."
    conda install --yes -c conda-forge r-neuralnet
    echo -e "\n"

    #install caTools R package
    echo "installing caTools..."
    conda install --yes -c r r-catools
    echo -e "\n"

    #install factoextra R package
    echo "installing factoextra..."
    conda install --yes -c conda-forge r-factoextra
    echo -e "\n"

    #install venndiagram R package
    echo "installing venndiagram..."
    conda install --yes -c bioconda r-venndiagram
    echo -e "\n"

    #install reshape2 R package
    echo "installing reshape2..."
    conda install --yes -c r r-reshape2
    echo -e "\n"

    #install limma R package
    echo "installing limma..."
    conda install --yes -c bioconda bioconductor-limma
    echo -e "\n"

    #install edger R package
    echo "installing edger..."
    conda install --yes -c bioconda bioconductor-edger
    echo -e "\n"

    #install gprofiler2 R package
    echo "installing gprofiler2..."
    conda install --yes -c conda-forge r-gprofiler2
    echo -e "\n"

    #install ggrepel R package
    echo "installing ggrepel..."
    conda install --yes -c bioconda r-ggrepel
    echo -e "\n"

    #install pathview R package
    echo "installing pathview..."
    conda install --yes -c bioconda bioconductor-pathview
    echo -e "\n"

    #install bedr R package
    echo "installing bedr..."
    conda install --yes -c bioconda r-bedr
    echo -e "\n"

    #install scales R package
    echo "installing scales..."
    conda install --yes -c r r-scales
    echo -e "\n"

    #install motifstack R package
    echo "installing motifstack..."
    conda install --yes -c bioconda bioconductor-motifstack
    echo -e "\n"

    #install snakemake
    if [ -z  `which snakemake` ]
    then
    echo "snakemake is not installed, or not in environment variables. Configuration check will add it to environment variables..."
    conda_snakemake="conda install --yes -c bioconda snakemake"
    eval $conda_snakemake
    else
    echo "snakemake is installed."
    fi

    #install mupdf-x11
    if [ -z  `which  mupdf-x11` ]
    then
    echo "mupdf is not installed, or not in environment variables. Configuration check will add it to environment variables..."
    conda_mupdf="conda install --yes -c conda-forge mupdf"
    eval $conda_mupdf
    else
    echo "mupdf is installed."
    fi

    #install wkhtmltopdf
    if [ -z  `which  wkhtmltopdf` ]
    then
    echo "wkhtmltopdf is not installed, or not in environment variables. Configuration check will add it to environment variables..."
    conda_wkhtmltopdf="conda install --yes -c conda-forge wkhtmltopdf"
    eval $conda_wkhtmltopdf
    else
    echo "wkhtmltopdf is installed."
    fi

    #install fastx_collapser
    # if [ -z  `which  fastx_collapser` ]
    # then
    # echo "fastx_toolkit is not installed, or not in environment variables. Configuration check will add it to environment variables..."
    # conda_fastx_toolkit="conda install --yes -c bioconda fastx_toolkit"
    # eval $conda_fastx_toolkit
    # else
    # echo "fastx_toolkit is installed."
    # fi

    #install cutadapt
    if [ -z  `which cutadapt` ]
    then
    echo "cutadapt is not installed, or not in environment variables. Configuration check will add it to environment variables..."
    conda_cutadapt="conda install --yes -c bioconda cutadapt"
    eval $conda_cutadapt
    else
    echo "cutadapt is installed."
    fi

    #install seqtk
    if [ -z  `which seqtk` ]
    then
    echo "seqtk is not installed, or not in environment variables. Configuration check will add it to environment variables..."
    conda_seqtk="conda install --yes -c bioconda seqtk"
    eval $conda_seqtk
    else
    echo "seqtk is installed."
    fi

    #install STAR
    if [ -z  `which STAR` ]
    then
    echo "STAR is not installed, or not in environment variables. Configuration check will add it to environment variables..."
    conda_STAR="conda install --yes -c bioconda star"
    eval $conda_STAR
    else
    echo "STAR is installed."
    fi

    #install samtools
    if [ -z  `which samtools` ]
    then
    echo "samtools is not installed, or not in environment variables. Configuration check will add it to environment variables..."
    conda_samtools="conda install --yes -c bioconda samtools"
    eval $conda_samtools
    else
    echo "samtools is installed."
    fi

    #install gtfToGenePred
    if [ -z  `which gtfToGenePred` ]
    then
    echo "gtfToGenePred is not installed, or not in environment variables. Configuration check will add it to environment variables..."
    conda_gtfToGenePred="conda install --yes -c bioconda ucsc-gtftogenepred"
    eval $conda_gtfToGenePred
    else
    echo "gtfToGenePred is installed."
    fi

    #install bedtools
    if [ -z  `which bedtools` ]
    then
    echo "bedtools is not installed, or not in environment variables. Configuration check will add it to environment variables..."
    conda_bedtools="conda install --yes -c bioconda bedtools"
    eval $conda_bedtools
    else
    echo "bedtools is installed."
    fi

    #locate rtsSeeker
    if [ `test -e /script/rtsSeeker` ]
    then
    echo "rtsSeeker is not exist, please check if script/rtsSeeker exist in your current directory (directory with psiFinder) or re-download the .zip package and unzip it..."
    else
    echo "rtsSeeker exist (Note: if you copy the psiFinder to other directory, make sure you copy the rtsSeeker/ directory with it as well)."
    fi

    #locate bedAnnotator
    if [ `test -e /script/bedAnnotator` ]
    then
    echo "bedAnnotator is not exist, please check if script/bedAnnotator exist in your current directory (directory with psiFinder) or re-download the .zip package and unzip it..."
    else
    echo "bedAnnotator exist (Note: if you copy the psiFinder to other directory, make sure you copy the annotation/ directory with it as well)."
    fi

    #locate ACAscan
    if [ `test -e /script/ACAscan` ]
    then
    echo "ACAscan is not exist, please check if annotation/ACAscan exist in your current directory (directory with psiFinder) or re-download the .zip package and unzip it..."
    else
    echo "ACAscan exist (Note: if you copy the psiFinder to other directory, make sure you copy the annotation/ directory with it as well)."
    fi

    #install viennarna
    echo "installing viennarna..."
    conda install --yes -c bioconda viennarna
    echo -e "\n"

    #install miranda
    if [ -z  `which miranda` ]
    then
    echo "miranda is not installed, or not in environment variables. Configuration check will add it to environment variables..."
    conda_miranda="conda install --yes -c bioconda miranda"
    eval $conda_miranda
    else
    echo "miranda is installed."
    fi

    #install homer
    if [ -z  `which homer` ]
    then
    echo "homer is not installed, or not in environment variables. Configuration check will add it to environment variables..."
    conda_homer="conda install --yes -c bioconda homer"
    eval $conda_homer
    else
    echo "homer is installed."
    fi

    echo "Finished: All software checking is done!"
    exit 0  # Exit normally.

.. note:: This step will automatically install some software through `conda <https://docs.conda.io/en/latest/>`_\ .

.. image:: /images/WSL_conda.png

.. tip:: For support or questions please make a post on `Biostars <http://biostars.org>`__. For feature requests or bug reports please open an issue on `github <https://github.com/worsteggs/psiFinder_readthedocs/issues>`__.
