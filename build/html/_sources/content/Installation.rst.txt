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

.. image:: /images/psiFinder_window.png

* *Mac OS*

For Mac users, after decompressing the file, open the ``Command Line`` file, cd to the ``psiFinder.app directory``, and input ``./psiFinder.app/Content/Macos/psiFinder`` at Command file to run.

.. note:: For Mac, you should place the downloaded genome file and/or the example data in the same directory of ``psiFinder`` flie mentioned above. (like ``psiFinder/genome`` and/or ``psiFinder/data``)

Configuration
---------------------------
Before uploading data, you should first configure the operating environment of the software, this step can be achieved through clicking the option ``global->configuration`` at the upper of the main window.

.. image:: /images/configuration.png

Once click ``CHECK``, psiFindeer will run ``0_configuration.sh``.

.. note:: This step will automatically install some software through `conda <https://docs.conda.io/en/latest/>`_\ .

.. tip:: For support or questions please make a post on `Biostars <http://biostars.org>`__. For feature requests or bug reports please open an issue on `github <https://github.com/worsteggs/psiFinder_readthedocs/issues>`__.
