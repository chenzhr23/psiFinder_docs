snakemake pipline
===========================



.. role:: red

.. contents::
    :local:

For a quick start of Ψ-sites data analyses, we developed a snakemake pipline and designed corresponding QT widget for easy running.

.. image:: /images/quickstart_snakemake.png

Input
------

Select Genome
**************
Current version of psiFinder snakemake pipline mainly support :red:`Group: Mammal->Genome: Homo_sapiens->Assembly: hg38`, more support for other species will coming soon.


Data upload
**************
Users should ``browse`` and locate `fastq|fasta fastq.gz/fasta.gz fq/fa fq.gz/fa.gz` file for both CMC-input and CMC-treated samples.

-  If ``Barcode removal`` is checked, it will enable 5'/3' barcode setting, remove ``LENGTH`` (`cutadapt <https://cutadapt.readthedocs.io/en/stable/guide.html>`_) bases from reads.
-  If ``Adapter removal`` is checked, it will enable 5'/3' adapter setting, sequence of an adapter ligated to the 5'/3' end will be cutadapt.

Basic analyses
**************
Current version of psiFinder snakemake pipline provides four basic analysis module of Ψ-sites: (1) identification; (2) annotation; (3) distribution; (4) snoRNA-guided target.

Mode options
**************
Different mode options are provided for selection on SVM/ANN/User-defined approaches.

Filter threshold
***********************
If User-defined ``Mode options`` is checked, then ``Filter thresholds`` is enable to customize the threshold of the six key metrics.


snakemake rules
----------------
Once click ``START``, psiFindeer will firstly run ``cp_unzip.sh`` to unzip and redirect the input sequencing reads files. Then, ``psiFinder.sh`` will be secondly run to apply snakemake rules selected by ``quickstart_config.yml``.

``cp_unzip.sh``

.. code:: bash

    #! /bin/bash

    if [[ $# -ne 5 ]]
    then
        echo 'Usage: ./'$0 ' input_sample input_path treat_sample treat_path log_file'
        exit 1
    fi

    input_sample=$1
    input_path=$2
    treat_sample=$3
    treat_path=$4
    log_file=$5


    echo -e "clear psiFinder.sh history"
    echo -e "Your quickstart job is starting, please wait..." > $log_file

    #.fastq.gz
    if [[ "$input_sample" == *.fastq.gz ]] ;
    then
        echo -e "$(tput setaf 2)input: unzipping and redirecting fastq file...\n$(tput sgr 0)" >> $log_file
        gunzip -c $input_sample  > $input_path &
    fi

    if [[ "$treat_sample" == *.fastq.gz ]] ;
    then
        echo -e "$(tput setaf 2)treat: unzipping and redirecting fastq file...\n$(tput sgr 0)" >> $log_file
        gunzip -c $treat_sample  > $treat_path &
    fi

    #.fasta.gz
    if [[ "$input_sample" == *.fasta.gz ]] ;
    then
        echo -e "$(tput setaf 2)input: unzipping and redirecting fasta file...\n$(tput sgr 0)" >> $log_file
        gunzip -c $input_sample  > $input_path
        seqtk seq -F '#' $input_path > ${input_path%.fasta}.fastq
        rm $input_path
    fi

    if [[ "$treat_sample" == *.fasta.gz ]] ;
    then
        echo -e "$(tput setaf 2)treat: unzipping and redirecting fasta file...\n$(tput sgr 0)" >> $log_file
        gunzip -c $treat_sample  > $treat_path
        seqtk seq -F '#' $treat_path > ${treat_path%.fasta}.fastq
        rm $treat_path
    fi

    #.fq.gz
    if [[ "$input_sample" == *.fq.gz ]] ;
    then
        echo -e "$(tput setaf 2)input: unzipping and redirecting fastq file...\n$(tput sgr 0)" >> $log_file
        gunzip -c $input_sample  > $input_path &
    fi

    if [[ "$treat_sample" == *.fq.gz ]] ;
    then
        echo -e "$(tput setaf 2)treat: unzipping and redirecting fastq file...\n$(tput sgr 0)" >> $log_file
        gunzip -c $treat_sample  > $treat_path &
    fi

    #.fa.gz
    if [[ "$input_sample" == *.fa.gz ]] ;
    then
        echo -e "$(tput setaf 2)input: unzipping and redirecting fasta file...\n$(tput sgr 0)" >> $log_file
        gunzip -c $input_sample  > $input_path
        seqtk seq -F '#' $input_path > ${input_path%.fasta}.fastq
        rm $input_path
    fi

    if [[ "$treat_sample" == *.fa.gz ]] ;
    then
        echo -e "$(tput setaf 2)treat: unzipping and redirecting fasta file...\n$(tput sgr 0)" >> $log_file
        gunzip -c $treat_sample  > $treat_path
        seqtk seq -F '#' $treat_path > ${treat_path%.fasta}.fastq
        rm $treat_path
    fi

    #.fastq
    if [[ "$input_sample" == *.fastq ]] ;
    then
        echo -e "$(tput setaf 2)input: redirecting fastq file...\n$(tput sgr 0)" >> $log_file
        cp -R $input_sample $input_path &
    fi

    if [[ "$treat_sample" == *.fastq ]] ;
    then
        echo -e "$(tput setaf 2)treat: redirecting fastq file...\n$(tput sgr 0)" >> $log_file
        cp -R $treat_sample $treat_path &
    fi

    #.fasta
    if [[ "$input_sample" == *.fasta ]] ;
    then
        echo -e "$(tput setaf 2)input: redirecting fasta file...\n$(tput sgr 0)" >> $log_file
        cp -R $input_sample $input_path
        seqtk seq -F '#' $input_path > ${input_path%.fasta}.fastq
        rm $input_path
    fi

    if [[ "$treat_sample" == *.fasta ]] ;
    then
        echo -e "$(tput setaf 2)treat: redirecting fasta file...\n$(tput sgr 0)" >> $log_file
        cp -R $treat_sample $treat_path
        seqtk seq -F '#' $treat_path > ${treat_path%.fasta}.fastq
        rm $treat_path
    fi

    #.fq
    if [[ "$input_sample" == *.fq ]] ;
    then
        echo -e "$(tput setaf 2)input: redirecting fastq file...\n$(tput sgr 0)" >> $log_file
        cp -R $input_sample $input_path &
    fi

    if [[ "$treat_sample" == *.fq ]] ;
    then
        echo -e "$(tput setaf 2)treat: redirecting fastq file...\n$(tput sgr 0)" >> $log_file
        cp -R $treat_sample $treat_path &
    fi

    #.fa
    if [[ "$input_sample" == *.fa ]] ;
    then
        echo -e "$(tput setaf 2)input: redirecting fasta file...\n$(tput sgr 0)" >> $log_file
        cp -R $input_sample $input_path
        seqtk seq -F '#' $input_path > ${input_path%.fasta}.fastq
        rm $input_path
    fi

    if [[ "$treat_sample" == *.fa ]] ;
    then
        echo -e "$(tput setaf 2)treat: redirecting fasta file...\n$(tput sgr 0)" >> $log_file
        cp -R $treat_sample $treat_path
        seqtk seq -F '#' $treat_path > ${treat_path%.fasta}.fastq
        rm $treat_path
    fi

    wait


    echo -e "Starting snakemake workflow..."
    #snakemake
    #fastq
    if [[ "$(basename ${input_path})" == *.fastq ]] && [[ "$(basename ${treat_path})" == *.fastq ]];
    then
        echo "check info: Input is fastq files"  >> $log_file
        nohup bash snakemake/psiFinder.sh -i $(basename ${input_path%.fastq}) -t $(basename ${treat_path%.fastq}) -c quickstart_config.yml >> $log_file 2>&1 &
        wait
    else
        echo "check info: Input is not fastq files" >> $log_file
    fi

    #fasta
    if [[ "$(basename ${input_path})" == *.fasta ]] && [[ "$(basename ${treat_path})" == *.fasta ]] ;
    then
        echo "check info: Input is fasta files" >> $log_file
        nohup bash snakemake/psiFinder.sh -i $(basename ${input_path%.fasta}) -t $(basename ${treat_path%.fasta}) -c quickstart_config.yml >> $log_file 2>&1 &
        wait
    else
        echo "check info: Input is not fasta files" >> $log_file
    fi

    echo -e "Finished snakemake workflow!" >> $log_file


``psiFinder.sh``

.. code:: bash

    #!/bin/bash

    usage() {                                      # Function: Print a help message.
      echo "Usage: $0 [ -h help] [ -i input ] [ -t treat ] [ -c quickstart_config.yml ]" 1>&2
    }
    exit_abnormal() {                              # Function: Exit with error.
      usage
      exit 1
    }


    while getopts ":h:i:t:c:" options; do

      case "${options}" in
        h|:)
          usage
          Help
          exit 0
          ;;
        i)
          input=${OPTARG}
          if ! [[ -n $input ]] ; then
            echo "You didn't set the input sample"
          fi
          ;;
        t)
          treat=${OPTARG}
          if ! [[ -n $treat ]] ; then
            echo "You didn't set the treat sample"
          fi
          ;;
        c)
          config=${OPTARG}
          if ! [[ -n $config ]] ; then
            echo "You didn't set the config yml file"
          fi
          ;;
        \?) # incorrect option
          echo "Error: -${OPTARG} Invalid option"
          exit_abnormal
          ;;
      esac
    done

    shift $(($OPTIND - 1))

    quickstart_config=$(cat $config)
    quickstart_array=($(echo $quickstart_config | tr ":" "\n"))

    echo -e "==Starting a quickstart job...==:\n"

    #cutadapt-se
    echo -e "ad_remove_input:${quickstart_array[41]}"
    if [[ "${quickstart_array[41]}" =~ "true"  ]] ;
    then
        echo -e "Getting input cutadapt fastq...\n"
        snakemake -s snakemake/Snakefile --cores 8 snakemake/cutadapt/input/"${input}".fastq
    else
        echo -e "Moving input cutadapt fastq...\n"
        mv snakemake/reads/input/"${input}".fastq snakemake/cutadapt/input
    fi

    echo -e "ad_remove_treat:${quickstart_array[43]}"
    if [[ "${quickstart_array[43]}" =~ "true"  ]] ;
    then
        echo -e "Getting treat cutadapt fastq...\n"
        snakemake -s snakemake/Snakefile --cores 8 snakemake/cutadapt/treat/"${treat}".fastq
    else
        echo -e "Moving treat cutadapt fastq...\n"
        mv snakemake/reads/treat/"${treat}".fastq snakemake/cutadapt/treat
    fi

    #STAR-index
    temp=$(echo -e "${quickstart_array[5]}" | sed "s/\"//g"  )
    DIRECTORY="snakemake/genome/${temp}"
    if [ -d "$DIRECTORY" ]; then
      echo "$DIRECTORY is not empty, STAR index have been built"

    else
      echo -e "No STAR index in $DIRECTORY, generating STAR index...\n"
      snakemake -s snakemake/Snakefile --cores 8 snakemake/genome/"${temp}"
    fi

    #rtsSeeker result
    echo -e "sites_identification:${quickstart_array[31]}"
    if [[ "${quickstart_array[31]}" =~ "true"  ]] ;
    then
        echo -e "Getting rtsSeeker result for Ψ sites identification...\n"
        snakemake -s snakemake/Snakefile --cores 8 snakemake/output/sites_identification/"${input}"_versus_"${treat}"/"${input}"_versus_"${treat}".bed
    fi

    #rtsSeeker ann
    echo -e "ann:${quickstart_array[49]}"
    if [[ "${quickstart_array[49]}" =~ "true"  ]] ;
    then
        echo -e "Getting ann result for Ψ sites identification...\n"
        snakemake -s snakemake/Snakefile --cores 8 snakemake/output/ann/"${input}"_versus_"${treat}"/"${input}"_versus_"${treat}"_ann_psi_prediction.bed
    fi

    #rtsSeeker svm
    echo -e "svm:${quickstart_array[45]}"
    if [[ "${quickstart_array[45]}" =~ "true"  ]] ;
    then
        echo -e "Getting svm result for Ψ sites identification...\n"
        snakemake -s snakemake/Snakefile --cores 8 snakemake/output/svm/"${input}"_versus_"${treat}"/"${input}"_versus_"${treat}"_svm_psi_prediction.bed
    fi

    #rtsSeeker user-defined
    echo -e "user-defined:${quickstart_array[51]}"
    if [[ "${quickstart_array[51]}" =~ "true"  ]] ;
    then
        echo -e "Getting user-defined result for Ψ sites identification...\n"
        snakemake -s snakemake/Snakefile --cores 8 snakemake/output/user_defined/"${input}"_versus_"${treat}"/"${input}"_versus_"${treat}"_user_defined_psi_prediction.bed
    fi

    #bedAnnotator
    echo -e "sites_annotation:${quickstart_array[33]}"
    if [[ "${quickstart_array[33]}" =~ "true"  ]] ;
    then
        echo -e "Getting bedAnnotator result for Ψ sites annotation...\n"
        snakemake -s snakemake/Snakefile --cores 8 snakemake/output/sites_annotation/"${input}"_versus_"${treat}"/"${input}"_versus_"${treat}"_add_seq_group_uniq.bed
    fi

    #metagene
    echo -e "metagene:${quickstart_array[47]}"
    if [[ "${quickstart_array[47]}" =~ "true"  ]] ;
    then
        echo -e "Getting metagene result for Ψ sites identification...\n"
        snakemake -s snakemake/Snakefile --cores 8 snakemake/output/meta_gene/"${input}"_versus_"${treat}"/"${input}"_versus_"${treat}"_metagene_pseudoU.pdf
    fi

    #ACAscan
    echo -e "sites_target_prediction:${quickstart_array[35]}"
    if [[ "${quickstart_array[35]}" =~ "true"  ]] ;
    then
        echo -e "Getting ACAscan result for Ψ sites target prediction...\n"
        snakemake -s snakemake/Snakefile --cores 8 snakemake/output/sites_target_prediction/"${input}"_versus_"${treat}"/"${input}"_versus_"${treat}"_fa.out
    fi

    #emit done signal
    echo -e "==All done!=="


``quickstart_config.yml``

.. code:: bash

    group: "Mammal"
    genome: "Homo_sapiens"
    assembly: "hg38"
    se_input: "/public/home/chenzr/PSI_Seq_brainCell/A1-A12-totalRNA-result/A1.cutadapt.extendedFrags.collapse.cutBarcodes.fa.gz"
    se_treat: "/public/home/chenzr/PSI_Seq_brainCell/A1-A12-totalRNA-result/A2.cutadapt.extendedFrags.collapse.cutBarcodes.fa.gz"
    se_input_prefix: "A1"
    se_treat_prefix: "A2"
    se_br5_input: "0"
    se_br3_input: "0"
    se_ad5_input: ""
    se_ad3_input: ""
    se_br5_treat: "0"
    se_br3_treat: "0"
    se_ad5_treat: ""
    se_ad3_treat: ""
    sites_identification: "true"
    sites_annotation: "true"
    sites_target_prediction: "true"
    br_remove_input: "false"
    br_remove_treat: "false"
    ad_remove_input: "false"
    ad_remove_treat: "false"
    svm: "true"
    metagene: "false"
    ann: "false"
    user_defined: "false"
    treatpre_Fold_thres: "0"
    preFold_FC_thres: "0"
    treataft_Fold_thres: "0"
    aftFold_FC_thres: "0"
    treatstoprate_thres: "0"

Output
--------
snakemake result are output to specific directory where psiFinder package is unziped.

.. code:: bash

    $ cd /the/directory/of/psiFinder/snakemake/
    $ tree -L 1
    .
    ├── cp_unzip.sh
    ├── genome
    ├── logs
    ├── output
    ├── psiFinder.sh
    ├── reads
    ├── script
    ├── Snakefile
    ├── star
    └── cutadapt

    7 directories, 3 files

    $ cd /the/directory/of/psiFinder/snakemake/output
    $ tree -L 1
    .
    ├── ann
    ├── meta_gene
    ├── sites_annotation
    ├── sites_identification
    ├── sites_target_prediction
    ├── svm
    └── user_defined

    7 directories, 0 files

.. note:: All user input will be recorded in a plain text file with suffix ``_quickstart_config.txt`` in psiFinder/config and help users to easily reload the previous config (by simply clicking ``CONFIG`` button).
