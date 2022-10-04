Transcripts Alignment
====================================

.. contents::
    :local:

psiFinder ``Transcripts Alignment`` utilized `STAR <https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/STARmanual.pdf>`_ to map sequencing reads to reference genome location and consist of two parts:

-  ``Generating genome indexes``
-  ``Running mapping jobs``

.. image:: /images/Transcripts_Index_Alignment.png

Once click ``BUILD``, then psiFinder will run ``2.0_STAR_index.sh`` to generate genome indexes.

.. code:: bash

    #! /bin/bash

    Help() {
    echo 'Usage: '$0 ' -a threads -b genome -c gtf -o outputpath'
    }


    usage() {                                      # Function: Print a help message.
      echo "Usage: $0 [ -h help] [ -a threads ] [ -b genome ] [ -c gtf ] [ -o outputpath ]" 1>&2
    }
    exit_abnormal() {                              # Function: Exit with error.
      usage
      exit 1
    }


    while getopts ":h:a:b:c:o:" options; do

      case "${options}" in
        h|:)
          usage
          Help
          exit 0
          ;;
        a)
          threads=${OPTARG}
          if ! [[ -n $threads ]] ; then
            echo "You didn't input the threads"
          fi
          ;;
        b)
          genome=${OPTARG}
          if ! [[ -n $genome ]] ; then
            echo "You didn't input the genome"
          fi
          ;;
        c)
          gtf=${OPTARG}
          if ! [[ -n $gtf ]] ; then
            echo "You didn't input the gtf"
          fi
          ;;
        o)
          outputpath=${OPTARG}
          if ! [[ -n $outputpath ]] ; then
            echo "You didn't input the outputpath"
          fi
          ;;
        \?) # incorrect option
          echo "Error: -${OPTARG} Invalid option"
          exit_abnormal
          ;;
      esac
    done

    shift $(($OPTIND - 1))
    date  ## echo the date at start


    echo star indexing start
    echo STAR --runThreadN $threads --runMode genomeGenerate --genomeDir $outputpath --genomeFastaFiles $genome --sjdbGTFfile $gtf --sjdbOverhang 99
        STAR \
        --runThreadN $threads \
        --runMode genomeGenerate \
        --genomeDir $outputpath\
        --genomeFastaFiles $genome \
        --sjdbGTFfile $gtf \
        --sjdbOverhang 99
    echo star indexing end


Once click ``START``, then psiFinder will run ``2.1_STAR_align.sh`` to run mapping jobs.

Input
*************************
Firstly, upload downloaded genome sequence and annotation file to build genome indexes. Next, locate the directory where your indexes were generated, and ``browse`` the syntropic read (the sequencing reads with the same strand direction as the genome) file for both CMC-input and CMC-treat samples. Finally, output file will be generated at the target out directory.

Perform transcripts mapping
*********************************

Built-in STAR argument for ``Running mapping jobs``:

.. code:: bash

    STAR \
	--genomeDir $starIndex \
	--readFilesIn $input \
	--runThreadN $threads \
	--genomeLoad NoSharedMemory \
	--alignEndsType EndToEnd \
	--outSAMtype BAM SortedByCoordinate \
	--outFileNamePrefix ${input%.fastq} \
	--outStd Log \
	--limitBAMsortRAM 60000000000 \
	--outFilterType BySJout \
	--outFilterMultimapScoreRange 0 \
	--outFilterMultimapNmax 30 \
	--outFilterMismatchNmax 15 \
	--outFilterMismatchNoverLmax 0.1 \
	--outFilterScoreMin 0 \
	--outFilterScoreMinOverLread 0 \
	--outFilterMatchNmin 15 \
	--outFilterMatchNminOverLread 0.8 \
	--alignIntronMin 20 \
	--alignIntronMax 1000000 \
	--alignMatesGapMax 1000000 \
	--seedSearchStartLmax 15 \
	--seedSearchStartLmaxOverLread 1 \
	--seedSearchLmax 0 \
	--seedMultimapNmax 20000 \
	--seedPerReadNmax 1000 \
	--seedPerWindowNmax 100 \
	--seedNoneLociPerWindow 20 \
	--alignSJoverhangMin 8 \
	--alignSJDBoverhangMin 1 \
	--quantMode TranscriptomeSAM GeneCounts \
	--outSAMmode Full \
	--outSAMattributes All \
	--outSAMunmapped None \
	--outSAMorder Paired \
	--outSAMprimaryFlag AllBestScore \
	--outSAMreadID Standard \
	--outReadsUnmapped Fastx \
	--alignEndsProtrude 150 ConcordantPair \
	--limitOutSJcollapsed 5000000

Output
*************************
Results in bam format (with ``Aligned.sortedByCoord.out.bam`` suffix) will be output to the target directory.

.. note:: All user input will be recorded in a plain text file with suffixes ``_STAR_build_config.txt`` and ``_STAR_align_config.txt`` in psiFinder/config and help users to easily reload the previous config (by simply clicking ``CONFIG`` button).
