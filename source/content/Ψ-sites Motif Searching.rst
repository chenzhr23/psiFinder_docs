Ψ-sites Motif Searching
=========================

.. contents::
    :local:

psiFinder ``Ψ-sites Motif Searching`` utilize `HOMER <http://homer.ucsd.edu/homer/motif/>`_ to generate motif searching result for Ψ-sites.

.. image:: /images/Motif_Searching.png


Input
---------------------------------------------

Users should choose to upload files (i.e. rtsSeeker result) to ``findMoitf`` QT widget in bed format and genome fasta file.


Search Ψ-sites Motif
---------------------------------------------

Once click ``START``, psiFindeer will run ``findMotif.sh`` and ``findMotifsGenome.r``.

``findMotif.sh``

.. code:: bash

    #! /bin/bash

    if [[ $# -ne 2 ]]
    then
        echo 'Usage: ./'$0 ' bedfile genome'
        exit 1
    fi

    bedfile=$1
    genome=$2


    awk 'FS=OFS="\t" {print $1,$2-10,$3+10,""$1"_"$2"_"$3"_"$6"",$5,$6}' ${bedfile} >${bedfile%.bed}_win21.bed
    echo "clearing previous Motif result"
    rm -rf ${bedfile%.bed}_motif

    echo "finding Motif for current task"
    findMotifsGenome.pl ${bedfile%.bed}_win21.bed $genome ${bedfile%.bed}_motif -size given -rna 2>${bedfile%.bed}_motif.log
    Rscript script/findMotifsGenome.r -f ${bedfile%.bed}_motif/homerResults/
    wkhtmltopdf -s A2 ${bedfile%.bed}_motif/homerResults.html ${bedfile%.bed}_motif/homerResults.pdf &> /dev/null
    mupdf-x11 ${bedfile%.bed}_motif/homerResults.pdf &> /dev/null

    echo -e "Finished: findMotif done!\n"
    echo -e "findMotif result in $(dirname ${bedfile})"


``findMotifsGenome.r``

.. code:: R

    #!/usr/bin/env Rscript

    suppressMessages(    if (!require('motifStack')) install.packages('motifStack');
    suppressMessages(library("motifStack"))

    option_list = list(
      make_option(c("-f", "--filepath"), type="character", default=NULL,
                  help="filepath [file]", metavar="character")
    );
    opt_parser = OptionParser(option_list=option_list);
    opt = parse_args(opt_parser);


    if (is.null(opt$filepath)){
      print_help(opt_parser);
      stop("Please provide -f filepath", call.=FALSE);
    }

    filepath = opt$filepath
    outFile_prefix = opt$outfile_prefix

    print(filepath)

    files <- c(list.files(path = filepath, pattern = "motif[0-9]*.motif$"))
    motiffile<-list()
    # invisible(lapply(1:10,function(i){motiffile[[i]]<<-read.table(paste(filepath,"motif",i,".motif",sep=""),header=F,skip=1)}))
    invisible(lapply(seq_along(files),function(i){motiffile[[i]]<<-read.table(paste(filepath,"motif",i,".motif",sep=""),header=F,skip=1)}))
    names(motiffile)<-paste("Ranked_Motif",seq_along(files),sep="")

    invisible(lapply(seq_along(files),function(j){
            names(motiffile[[j]]) <- c("A","C","G","U")
            A <- motiffile[[j]][,1]
            C <- motiffile[[j]][,2]
            G <- motiffile[[j]][,3]
            U <- motiffile[[j]][,4]
            data <- rbind(A,C,G,U)
            pcm <- data[,1:ncol(data)]
            rownames(pcm) <- c("A","C","G","U")
            motif <- new("pcm", mat=as.matrix(pcm), name=names(motiffile[j]))
            pdf.options(reset = TRUE, onefile = FALSE)
            pdf(paste(filepath,names(motiffile[j]),".pdf",sep=""),height=2,width=5)
            plot(motif)
            invisible(dev.off() )
    }))


Output
--------

Information
************

.. code:: bash

    $ cd /the/directory/of/out_file_dir
    $ tree -L 1
    .
    ├── Day0_Day4_common.bed
    ├── Day0_Day4_common_motif
    ├── Day0_Day4_common_motif.log
    └── Day0_Day4_common_win21.bed

    1 directory, 3 files

Diagram
********
File ``homerResults.pdf`` (in a directory which named by the prefix of the input rtsSeeker bed file) is a motif enrichment report of Ψ-sites motif on input Ψ-sites. Corresponding file with pattern ``Ranked_Motif*.pdf`` in ``homerResults/`` save all the motif diagram.

.. image:: /images/homerResults.png

.. image:: /images/homerResults_output.png

.. image:: /images/homerResults_demo.png

.. note:: All user input will be recorded in a plain text file with suffix ``_findMotif_config.txt`` in psiFinder/config and help users to easily reload the previous config (by simply clicking ``CONFIG`` button).
