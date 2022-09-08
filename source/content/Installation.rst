Installation
=============

``Remember`` -- psiFinder are available for **Linux os** as well as for
**Mac os**!

.. contents::
    :local:

Download psiFinder
--------------------
*psiFinder is free and open to all users.*

If you're here for the first time, download psiFinder first, and decompress it into a directory whatever you like.

`Download psiFinder v.0.1 for Linux x64 <https://mega.nz/file/wDdwyCZY#KasVu7WPJfKLDpSh_nnGfrBk5ho14QWnToQDHgraqaU>`_

`Download psiFinder v.0.1 for Mac OSX <https://mega.nz/file/ROsBkYiY#IFZ56zYR-3j7dCuz-34UF3r-LU7GZx-TkHdURTJ-5zI>`_

Windows support is currently a work-in-progress. Stay tuned and visit back for a Windows compatible executable.

.. note:: For users who want to remote operation a linux system in windows, the SSH client must be capable of X11-Forwarding, such as `MobaXterm <https://mobaxterm.mobatek.net/>`_ et al.

Download Genome annotation file
---------------------------------

Users should download the genome annotation file for the species of interest, and decompress these files into the same directory of psiFinder.

.. There are two types of genome annotation file, one is the indexes with spike-in pre-built by us. The other is HISAT2 official indexes.

* **pre-built indexes**

`Download hg38 <https://mega.nz/file/E4lH3IRD#EcmYwJPJEi9T706fyUgejvA2AFBaejdjCA9ACP8uKiE>`_

`Download mm10 <https://mega.nz/file/Ess2gYAZ#PyfipRKvbL17OA6w6kh60rPDP-bNm3hyqbszxiS4et8>`_

.. **HISAT2 indexes**


.. `Download H. sapiens GRCh38 genome_snp <https://genome-idx.s3.amazonaws.com/hisat/grch38_snp.tar.gz>`_

.. `Download M. musculus GRCm38 genome_snp <https://cloud.biohpc.swmed.edu/index.php/s/grcm38_snp/download>`_

.. .. note:: This step will automatically install some software through `conda <https://docs.conda.io/en/latest/>`_\ .

Download example data
-----------------------

psiFinder provide some example data (mouse) for the user to test. If want to do so, you should download the example data, and decompress these files into the same directory of psiFinder.

`Download example data <https://mega.nz/file/Axk3wa6S#eEXffvKWwFO3Zz07PaGr3elSSFdRF92swOeIJwUn05c>`_


Install conda
---------------------------------
psiFinder needs to use conda to configure the environment, please install conda first:

* *Linux*

.. code:: bash

    wget -c https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
    chmod 777 Miniconda3-latest-Linux-x86_64.sh
    bash Miniconda3-latest-Linux-x86_64.sh
    source ~/.bashrc

* *Mac OS*

.. code:: bash

    curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh
    chmod 777 Miniconda3-latest-MacOSX-x86_64.sh
    sh Miniconda3-latest-MacOSX-x86_64.sh
    source ~/.bash_profile


Starting the software
---------------------------------

* *Linux*

Once the psiFinder mian program and the Genome annotation are downloaded and unzipped properly, go to the directory of psiFinder:

.. code:: bash

    $ cd /the/directory/of/unzipped_psiFinder

    # see all files, don't run.
    $ tree -L 1
    .
    ├── AppRun -> psiFinder
    ├── CPAT-3.0.0
    ├── genome
    ├── default.desktop
    ├── psiFinder
    ├── psiFinder.sh
    ├── lib
    ├── libexec
    ├── libXss.so.1
    ├── pack.sh
    ├── plugins
    ├── qss
    ├── qt.conf
    ├── resources
    ├── snakemake
    └── translations

    9 directories, 7 files


The operation of psiFinder needs to grant executable permission to the script of the software. Please execute the following code before running the software:

.. code:: bash

 $ chmod 777 ./snakemake/script/*
 $ chmod 777 ./psiFinder

Then, upon executed the ./psiFinder, the following window appears:

.. code:: bash

	$ ./psiFinder.sh

.. image:: /images/introduction2.png

* *Mac OS*

For Mac users, after decompressing the file, open the ``Command Line`` file, cd to the ``psiFinder.app directory``, and input ``./psiFinder.app/Content/Macos/psiFinder`` at Command file to run.

.. note:: For Mac, you should place the downloaded genome file and/or the example data in the same directory of ``psiFinder`` flie mentioned above. (like ``psiFinder/genome`` and/or ``psiFinder/data``)

Configuration
---------------------------
Before uploading data, you should first configure the operating environment of the software, this step can be achieved through clicking the *Configuration* button at the bottom of the main window.

.. image:: /images/Configuration2.png

.. note:: This step will automatically install some software through `conda <https://docs.conda.io/en/latest/>`_\ .

.. tip:: For support or questions please make a post on `Biostars <http://biostars.org>`__. For feature requests or bug reports please open an issue on `github <https://github.com/worsteggs/psiFinder_readthedocs/issues>`__.
