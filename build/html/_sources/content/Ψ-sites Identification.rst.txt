Ψ-sites Identification
=======================

.. role:: red

.. contents::
    :local:

psiFinder ``Ψ-sites Identification`` was developed with a core C program--``rtsSeeker`` (as the core of psiFinder), for seeking the reverse transcriptase stop sites (RTS) and output the overall RTS information.

.. image:: /images/Identification.png

To build rRNA datasets for overall Ψ identification, we first download files from three databases: Modomics, snOPY, and snoRNABase. Since human 5.8s rRNA (chr21:8212571-8212727), 18s rRNA (chr21:8209631-8211499), and 28s rRNA (chr21:8213888-8218941) have multiple copies in the genome, we extracted one of the copies for subsequent analyses. We calculate and gain the relative position of Ψ sites on each rRNA. Centering on these sites and extending 10 nt of upstream and downstream, we then extract 21 nt sequences in bulk and map them to human genome (-v 0 -a -m 1000 --best --strata) using bowtie. A total of 96 Ψ sites were obtained in the above copy and used as a positive subset of the dataset to perform ROC evaluation or SVM modeling. Meanwhile, using bedtools (v2.30.0) intersect tools (-v -s), the sites that were identified as not Ψ sites were extracted from the above copy and used as a negative subset of the dataset.

.. image:: /images/identification_rtsSeeker_workflow.svg

Input
---------

Users should choose to upload files in bam format for both CMC-input and CMC-treated groups. Upload reference genome file in fasta format and set the output path of directory. Input the file name for the output information (i.e. file name with ``.bed`` suffix).

.. note:: STAR input/treat alignment Uniquely mapped reads % should be greater than 30%, or the rtsSeeker may fail to calculate the valid information.

To set the output path of directory, users can right click in the pop up file chooser (click ``browse``) of ``Out file dir<path>``

.. image:: /images/newfolder.png

Optional
---------
Multiple option are provided for rtsSeeker argument control.

customize rtsSeeker argument
************************************
if your want to customize the rtsSeeker argument, do not check :red:`get all psi information` or check :red:`plot ROC` option. Once ``START``, psiFinder will run ``3.0_rtsSeeker.sh``.

``3.0_rtsSeeker.sh``

.. code:: bash

    #!/bin/bash
    Help() {
    echo
    "rtsSeeker_wrap.sh example:
     bash rtsSeeker_wrap.sh -c rtsSeeker command -g genome
    "
    }


    usage() {                                      # Function: Print a help message.
      echo "Usage: $0 [ -h help] [ -c cmd ] [ -g genome]" 1>&2
    }
    exit_abnormal() {                              # Function: Exit with error.
      usage
      exit 1
    }


    while getopts ":h:c:g:" options; do

      case "${options}" in
        h|:)
          usage
          Help
          exit 0
          ;;
        c)
          cmd=${OPTARG}
          if ! [[ -n $cmd ]] ; then
            echo "You didn't input the cmd"
          fi
          ;;
        g)
          genome=${OPTARG}
          if ! [[ -n $genome ]] ; then
            echo "You didn't input the genome"
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
    echo "Starting: one rtsSeeker job is starting..."

    echo "Normal rtsSeeker mode: rtsSeeker start..." > rtsSeeker_wrap.log

    if [ ! -f "${genome}.fai"  ]
    then
    echo "generate fai file" >> rtsSeeker_wrap.log
    samtools faidx $genome
    else
    echo "fai file exist" >> rtsSeeker_wrap.log
    fi

    eval $cmd >> rtsSeeker_wrap.log 2>&1 &
    wait
    echo 'Genome-wise pseudouridylated sites identification is done by rtsSeeker!' >> rtsSeeker_wrap.log

    echo "rtsSeeker_wrap.sh running is done!"
    echo -e "\n $(tput setaf 3)Succeed: one rtsSeeker job is finished!$(tput sgr 0)"
    exit 0  # Exit normally.

get all psi information
*************************
If this option is checked, ``Detail`` argument input will be blocked and psiFinder will run ``get_all.sh``.

``get_all.sh``

.. code:: bash

    #! /bin/bash

    if [[ $# -ne 4 ]]
    then
        echo 'Usage: ./'$0 ' genome.fa treatment.bam input.bam out.bed'
        exit 1
    fi

    genome=$1
    treatment=$2
    input=$3
    out=$4
    outFileName=${out%.bed}

    echo "rtsSeeker get All mode: rtsSeeker start..."
    if [ ! -f "${genome}.fai"  ]
    then
    echo "generate faifile"
    samtools faidx $genome
    else
    echo "faifile exist"
    fi

    $(dirname "$0")/rtsSeeker --fa $genome --fai ${genome}.fai --treat $treatment --input $input -p 1.5 -t 0 -r 0 -M 1 --gene $(dirname "$0")/hg38.gencode.v30.tRNA.refseqNcRNA.geneAnno.bed12 -f 0 -m 0 -s -w 20 -o ${outFileName}_all.bed 2>${outFileName}_all_rtsSeeker.log
    awk 'FS=OFS="\t" {if($1~/^chr[0-9|a-z|A-Z]*$/ && $10=="T" && $14>5){print $0}}' ${outFileName}_all.bed > ${outFileName}_filt.bed
    awk 'FS=OFS="\t" {$25=sprintf("%.2f",$25);$27=sprintf("%.2f",$27);$28=sprintf("%.2f",$28);$30=sprintf("%.2f",$30);$31=sprintf("%.2f",$31);$33=sprintf("%.2f",$33)}1' ${outFileName}_filt.bed > tmp && mv tmp ${outFileName}_filt.bed
    echo -e "All rtsSeeker result in ${outFileName}_all.bed"


plot ROC
*************************
If this option is checked, ``Detail`` argument input will be blocked and psiFinder will run ``3.1_rtsSeeker_roc.sh`` and ``roc.r``.

``3.1_rtsSeeker_roc.sh``

.. code:: bash

    #! /bin/bash

    if [[ $# -ne 4 ]]
    then
        echo 'Usage: '$0 ' genome.fa treatment.bam input.bam out.bed'
        exit 1
    fi

    genome=$1
    treatment=$2
    input=$3
    out=$4
    outFileName=${out%.bed}

    echo "rtsSeeker ROC mode: rtsSeeker start..."
    if [ ! -f "${genome}.fai"  ]
    then
    echo "generate faifile"
    samtools faidx $genome
    else
    echo "faifile exist"
    fi

    # ROC
    echo "$(dirname "$0")/rtsSeeker --fa $genome --fai ${genome}.fai --treat $treatment --input $input -p 1.5 -t 5 -r 0.05 -M 1 --gene $(dirname "$0")/hg38.gencode.v30.tRNA.refseqNcRNA.geneAnno.bed12 -f 1 -m 0 -s -n -w 20 -o ${outFileName}_roc_all.bed 2>${outFileName}_roc_all_rtsSeeker.log" > ${outFileName}_rtsSeeker_cmd.log
    $(dirname "$0")/rtsSeeker --fa $genome --fai ${genome}.fai --treat $treatment --input $input -p 1.5 -t 5 -r 0.05 -M 1 --gene $(dirname "$0")/hg38.gencode.v30.tRNA.refseqNcRNA.geneAnno.bed12 -f 1 -m 0 -s -n -w 20 -o ${outFileName}_roc_all.bed 2>${outFileName}_roc_all_rtsSeeker.log # no -c remove duplicate change model to 1 and add gene annotation add -p 1.5
    cat ${outFileName}_roc_all.bed > $out

    awk 'FS=OFS="\t" {if($1~/^chr[0-9|a-z|A-Z]*$/ && $10=="T" && $14>10){print $0}}' ${outFileName}_roc_all.bed > ${outFileName}_roc_filt.bed
    bedtools intersect -a ${outFileName}_roc_filt.bed -b $(dirname "$0")/hg38_human_chr21_rRNA_known_pseudoU_SingleSites.bed -s > ${outFileName}_knowpse.bed
    bedtools intersect -a ${outFileName}_roc_filt.bed -b $(dirname "$0")/rrna_chr21.bed -s > ${outFileName}_rrna.bed
    bedtools intersect -a ${outFileName}_rrna.bed -b $(dirname "$0")/hg38_human_chr21_rRNA_known_pseudoU_SingleSites.bed -s -v > ${outFileName}_notpse.bed
    cat ${outFileName}_knowpse.bed |awk 'FS=OFS="\t" {print $0,"1"}' > ${outFileName}_knowpse.txt
    cat ${outFileName}_notpse.bed |awk 'FS=OFS="\t" {print $0,"0"}' > ${outFileName}_notpse.txt
    echo "getting roc_plot.txt"
    cat ${outFileName}_knowpse.txt ${outFileName}_notpse.txt > ${outFileName}_roc_plot.txt

    echo "generate ROC plot..."
    nohup Rscript $(dirname "$0")/roc.r -f ${outFileName}_roc_plot.txt -t ${outFileName}_roc_filt.bed -r $(dirname "$0")/hg38_human_chr21_rRNA_known_pseudoU_SingleSites.bed -s $(dirname "$0")/hg38.psiU.SingleSites.bed -o ${outFileName} > ${outFileName}_roc_bestthres.log 2>&1 &
    wait
    cat ${outFileName}_roc_bestthres.log

    mupdf-x11 ${outFileName}_six_variables_rRNA_violinplot.pdf &> /dev/null
    mupdf-x11 ${outFileName}_roc_summary.pdf &> /dev/null
    mupdf-x11 ${outFileName}_roc_best_evaluation.pdf &> /dev/null

    echo "rtsSeeker ROC end"
    echo -e "rtsSeeker ROC result in $(dirname ${outFileName}_roc_psi_prediction.bed)"

``roc.r``

.. code:: R

    suppressMessages(library("optparse"))
    suppressMessages(library("pROC"))
    suppressMessages(library("ggplot2"))
    suppressMessages(library("ggpubr"))
    suppressMessages(library("cowplot"))
    suppressMessages(library("dplyr"))
    suppressMessages(library("gridExtra"))
    suppressMessages(library("reshape2"))
    suppressMessages(library("stringr"))
    suppressMessages(library("RColorBrewer"))

    option_list = list(
      make_option(c("-f", "--rocfile"), type="character", default=NULL,
                  help="ROC of single sites [file]", metavar="character"),
      make_option(c("-r", "--rRNAfile"), type="character", default=NULL,
                  help="hg38_human_chr21_rRNA_known_pseudoU_SingleSites.bed [file]", metavar="character"),
      make_option(c("-s", "--rRNAfile2"), type="character", default=NULL,
                  help="hg38.psiU.SingleSites.bed [file]", metavar="character"),
      make_option(c("-t", "--filtfile"), type="character", default=NULL,
                  help="filt file [file]", metavar="character"),
      make_option(c("-o", "--outfile_prefix"), type="character", default=NULL,
                  help="output file name [default= %default]", metavar="character")
    );
    opt_parser = OptionParser(option_list=option_list);
    opt = parse_args(opt_parser);

    if (is.null(opt$rocfile)|| is.null(opt$filtfile) || is.null(opt$rRNAfile) || is.null(opt$rRNAfile2) || is.null(opt$outfile_prefix) ){
      print_help(opt_parser);
      stop("Please provide -f rocfile, -r rRNAfile, -s rRNAfile2, -t filtfile  and -o outfile_prefix option", call.=FALSE);
    }

    ROCfile = opt$rocfile
    rRNAfile = opt$rRNAfile
    rRNAfile2 = opt$rRNAfile2
    filtfile = opt$filtfile
    outFile_prefix = opt$outfile_prefix

    print(ROCfile)
    print(rRNAfile)
    print(rRNAfile2)
    print(filtfile)
    print(outFile_prefix)

    # roc_plot.txt
    ROC_data<-read.table(ROCfile,head=F)
    colnames(ROC_data)<-c("chrom",#1
        "chromStart",#2
        "chromEnd",#3
        "name",#4
        "foldChange",#5
        "strand",#6
        "geneName",#7
        "geneStart",#8
        "geneEnd",#9
        "base",#10
        "treatPval",#11
        "ctrlPval",#12
        "minusPval",#13
        "treatStopNum",#14
        "treatStopRPM",#15
        "treatPreStopNum",#16
        "treatAfterStopNum",#17
        "treatReadthroughNum",#18
        "ctrlStopNum",#19
        "ctrlStopRPM",#20
        "ctrlPreStopNum",#21
        "ctrlAfterStopNum",#22
        "ctrlReadthroughNum",#23
        "stopRpmFC",#24
        "treatPreRpmFold",#25
        "ctrlPreRpmFold",#26
        "preRpmFoldRatio",#27
        "treatAfterRpmFold",#28
        "ctrlAfterRpmFold",#29
        "afterRpmFoldRatio",#30
        "treatStopRatio",#31
        "ctrlStopRatio",#32
        "stopRatioFC",#33
        "treatStopMeanNum",#34
        "treatStopMeanFold",#35
        "ctrlStopMeanNum",#36
        "ctrlStopMeanFold",#37
        "treatStopMeanFoldRatio",#38
        "extendSeq",#39
        "class")#40
    ROC_data$base<-"T"
    ROC_data_sel<-ROC_data %>% select(treatPreRpmFold,treatAfterRpmFold,treatStopMeanFold,treatStopRatio,preRpmFoldRatio, afterRpmFoldRatio, stopRatioFC,treatStopMeanFoldRatio,class)


    #rRNA-psi-non-psi visualization
    ROC_data_melt<-melt(ROC_data[,c(11:38,40)],id.vars = "class")
    ROC_data_melt$class<-str_replace(as.character(ROC_data_melt$class),"0","non-psi")
    ROC_data_melt$class<-str_replace(as.character(ROC_data_melt$class),"1","psi")

    data_summary <- function(x) {
       m <- mean(x)
       ymin <- m-sd(x)
       ymax <- m+sd(x)
       return(c(y=m,ymin=ymin,ymax=ymax))
    }

    my_comparisons <- list( c("non-psi", "psi") )

    #treatPreRpmFold
    treatPreRpmFold<-ROC_data_melt%>%filter(str_detect(.$variable,"^treatPreRpmFold$"))
    treatPreRpmFold_bp <- ggplot(treatPreRpmFold, aes(x=class, y=log2(value), fill=class)) +
      stat_boxplot(geom = "errorbar",
                   width = 0.15) +
      geom_violin(trim=FALSE,alpha=0.8)+
      stat_summary(fun.data=data_summary,geom="crossbar", width=0.2,color="black")+
    labs(x="group", y = "log2(treatPreRpmFold)")
    treatPreRpmFold_bp<-treatPreRpmFold_bp + scale_fill_manual(values=c(brewer.pal(9,"Set1")[2],brewer.pal(9,"Set1")[1])) + theme_classic()+
    theme(text=element_text(size=8),legend.position = "none") + font("xy.text", size = 8)+
    stat_compare_means(method = "t.test",comparisons = my_comparisons,label.y = max(log2(treatPreRpmFold$value))+abs(quantile(log2(treatPreRpmFold$value))[1]))+
    scale_y_continuous(breaks = scales::pretty_breaks(n = 5))+
    scale_y_continuous(expand = expansion(mult = c(0, 0.1)))+
    scale_fill_manual(values=c("#129a92","#1e95d4"))

    #treatAfterRpmFold
    treatAfterRpmFold<-ROC_data_melt%>%filter(str_detect(.$variable,"^treatAfterRpmFold$"))
    treatAfterRpmFold_bp <- ggplot(treatAfterRpmFold, aes(x=class, y=log2(value), fill=class)) +
      stat_boxplot(geom = "errorbar",
                   width = 0.15) +
      geom_violin(trim=FALSE,alpha=0.8)+
      stat_summary(fun.data=data_summary,geom="crossbar", width=0.2,color="black")+
    labs(x="group", y = "log2(treatAfterRpmFold)")
    treatAfterRpmFold_bp<-treatAfterRpmFold_bp + scale_fill_manual(values=c(brewer.pal(9,"Set1")[2],brewer.pal(9,"Set1")[1])) + theme_classic()+
    theme(text=element_text(size=8), legend.position = "none") + font("xy.text", size = 8)+
    stat_compare_means(method = "t.test",comparisons = my_comparisons,label.y = max(log2(treatAfterRpmFold$value))+abs(quantile(log2(treatAfterRpmFold$value))[1]))+
    scale_y_continuous(breaks = scales::pretty_breaks(n = 5))+
    scale_y_continuous(expand = expansion(mult = c(0, 0.1)))+
    scale_fill_manual(values=c("#129a92","#1e95d4"))

    #preRpmFoldRatio
    preRpmFoldRatio<-ROC_data_melt%>%filter(str_detect(.$variable,"^preRpmFoldRatio$"))
    preRpmFoldRatio_bp <- ggplot(preRpmFoldRatio, aes(x=class, y=log2(value), fill=class)) +
      stat_boxplot(geom = "errorbar",
                   width = 0.15) +
      geom_violin(trim=FALSE,alpha=0.8)+
      stat_summary(fun.data=data_summary,geom="crossbar", width=0.2,color="black")+
    labs(x="group", y = "log2(preRpmFoldRatio)")
    preRpmFoldRatio_bp<-preRpmFoldRatio_bp + scale_fill_manual(values=c(brewer.pal(9,"Set1")[2],brewer.pal(9,"Set1")[1])) + theme_classic()+
    theme(text=element_text(size=8), legend.position = "none") + font("xy.text", size = 8)+
    stat_compare_means(method = "t.test",comparisons = my_comparisons,label.y = max(log2(preRpmFoldRatio$value))+abs(quantile(log2(preRpmFoldRatio$value))[1]))+
    scale_y_continuous(breaks = scales::pretty_breaks(n = 5))+
    scale_y_continuous(expand = expansion(mult = c(0, 0.1)))+
    scale_fill_manual(values=c("#129a92","#1e95d4"))

    #afterRpmFoldRatio
    afterRpmFoldRatio<-ROC_data_melt%>%filter(str_detect(.$variable,"^afterRpmFoldRatio$"))
    afterRpmFoldRatio_bp <- ggplot(afterRpmFoldRatio, aes(x=class, y=log2(value), fill=class)) +
      stat_boxplot(geom = "errorbar",
                   width = 0.15) +
      geom_violin(trim=FALSE,alpha=0.8)+
      stat_summary(fun.data=data_summary,geom="crossbar", width=0.2,color="black")+
      labs(x="group", y = "log2(afterRpmFoldRatio)")
    afterRpmFoldRatio_bp<-afterRpmFoldRatio_bp + scale_fill_manual(values=c(brewer.pal(9,"Set1")[2],brewer.pal(9,"Set1")[1])) + theme_classic()+
    theme(text=element_text(size=8),legend.position = "none") + font("xy.text", size = 8)+
    stat_compare_means(method = "t.test",comparisons = my_comparisons,label.y = max(log2(afterRpmFoldRatio$value))+abs(quantile(log2(afterRpmFoldRatio$value))[1]))+
    scale_y_continuous(breaks = scales::pretty_breaks(n = 5))+
    scale_y_continuous(expand = expansion(mult = c(0, 0.1)))+
    scale_fill_manual(values=c("#129a92","#1e95d4"))

    #treatStopRatio
    treatStopRatio<-ROC_data_melt%>%filter(str_detect(.$variable,"^treatStopRatio$"))
    treatStopRatio<-treatStopRatio[treatStopRatio$value!=0,]
    treatStopRatio_bp <- ggplot(treatStopRatio, aes(x=class, y=log2(value), fill=class)) +
      stat_boxplot(geom = "errorbar",
                   width = 0.15) +
      geom_violin(trim=FALSE,alpha=0.8)+
      stat_summary(fun.data=data_summary,geom="crossbar", width=0.2,color="black")+
      labs(x="group", y = "log2(treatStopRatio)")
    treatStopRatio_bp<-treatStopRatio_bp + scale_fill_manual(values=c(brewer.pal(9,"Set1")[2],brewer.pal(9,"Set1")[1])) + theme_classic()+
    theme(text=element_text(size=8),legend.position = "none") + font("xy.text", size = 8)+
    stat_compare_means(method = "t.test",comparisons = my_comparisons,label.y = max(log2(treatStopRatio$value))+abs(quantile(log2(treatStopRatio$value))[1]))+
    scale_y_continuous(breaks = scales::pretty_breaks(n = 5))+
    scale_y_continuous(expand = expansion(mult = c(0, 0.1)))+
    scale_fill_manual(values=c("#129a92","#1e95d4"))

    stopRatioFC<-ROC_data_melt%>%filter(str_detect(.$variable,"^stopRatioFC$"))
    stopRatioFC_bp <- ggplot(stopRatioFC, aes(x=class, y=log2(value), fill=class)) +
      stat_boxplot(geom = "errorbar",
                   width = 0.15) +
      geom_violin(trim=FALSE,alpha=0.8)+
      stat_summary(fun.data=data_summary,geom="crossbar", width=0.2,color="black")+
      labs(x="group", y = "log2(stopRatioFC)")
    stopRatioFC_bp<-stopRatioFC_bp + scale_fill_manual(values=c(brewer.pal(9,"Set1")[2],brewer.pal(9,"Set1")[1])) + theme_classic()+
    theme(text=element_text(size=8),legend.position = "none") + font("xy.text", size = 8)+
    stat_compare_means(method = "t.test",comparisons = my_comparisons,label.y = max(log2(stopRatioFC$value))+abs(quantile(log2(stopRatioFC$value))[1]))+
    scale_y_continuous(breaks = scales::pretty_breaks(n = 5))+
    scale_y_continuous(expand = expansion(mult = c(0, 0.1)))+
    scale_fill_manual(values=c("#129a92","#1e95d4"))


    pdf(paste(outFile_prefix,"_six_variables_rRNA_violinplot.pdf",sep=""),width=5,height=6)
    plot_grid(
      treatPreRpmFold_bp,
      preRpmFoldRatio_bp,
      treatAfterRpmFold_bp,
      afterRpmFoldRatio_bp,
      treatStopRatio_bp,
      stopRatioFC_bp,
      align = "hv",
      labels = c('A','B','C','D','E','F'),ncol=2,nrow=3)
    invisible(dev.off())


    cat("=====================treatPreRpmFold/controlPreRpmFold=========================\n")
    #treatPreRpmFold controlPreRpmFold
    pre_input=roc(ROC_data$class,ROC_data$ctrlPreRpmFold,smooth=FALSE, print.auc=TRUE, col="#e41a1c", plot = FALSE, print.thres="best", print.thres.best.method="youden", levels = c(0,1), direction='<',auc=T, ci=T)
    pdf(paste(outFile_prefix, '_roc_treatPreRpmFold.pdf', sep=""))
    pre_CMC=roc(ROC_data$class,ROC_data$treatPreRpmFold,smooth=FALSE, print.auc=TRUE, col="#e41a1c", plot = TRUE, print.thres="best", print.thres.best.method="youden", levels = c(0,1), direction='<',auc=T, ci=T)
    invisible(dev.off())
    pre_input_label = paste("controlPreRpmFold AUC:",sprintf("%.3f",pre_input$auc))
    pre_CMC_label = paste("treatPreRpmFold AUC:",sprintf("%.3f",pre_CMC$auc))
    preList=list(CMC=pre_CMC, Input=pre_input)
    names(preList) <- c(pre_CMC_label,pre_input_label)
    pre_plot<-ggroc(preList,size=0.8)
    pre_plot<-pre_plot+geom_abline(intercept = 1,slope = 1,col="black",size=0.6,linetype=2)+theme_classic()+theme(legend.position="top",legend.title=element_blank(),panel.background=element_rect(fill="white",color="black"))
    pre_plot<-pre_plot+ guides(colour = guide_legend(nrow = 2))
    pre_plot<-pre_plot+ annotate("text", x = .5, y = .5,
                              label = paste("Method: ",roc.test(pre_input, pre_CMC)$method,"\np.value: ",signif(roc.test(pre_input, pre_CMC)$p.value,3),sep=""),
                              size = 3)
    print(paste(outFile_prefix,"_pre_plot.pdf",sep=""))
    pdf(paste(outFile_prefix,"_pre_plot.pdf",sep=""),width =4 ,height = 3)
    pre_plot
    invisible(dev.off())
    print(paste(outFile_prefix,"_pre_plot.png",sep=""))
    png(paste(outFile_prefix,"_pre_plot.png",sep=""))
    pre_plot
    invisible(dev.off())

    roc.test(pre_input, pre_CMC)

    cat("=====================treatAfterRpmFold/controlAfterRpmFold=========================\n")
    #treatAfterRpmFold controlAfterRpmFold
    aft_input=roc(ROC_data$class,ROC_data$ctrlAfterRpmFold,smooth=FALSE, print.auc=TRUE, col="#e41a1c", plot = FALSE, print.thres="best", print.thres.best.method="youden", levels = c(0,1), direction='<',auc=T, ci=T)
    pdf(paste(outFile_prefix, '_roc_treatAfterRpmFold.pdf', sep=""))
    aft_CMC=roc(ROC_data$class,ROC_data$treatAfterRpmFold,smooth=FALSE, print.auc=TRUE, col="#e41a1c", plot = TRUE, print.thres="best", print.thres.best.method="youden", levels = c(0,1), direction='<',auc=T, ci=T)
    invisible(dev.off())
    aft_input_label = paste("controlAfterRpmFold AUC:",sprintf("%.3f",aft_input$auc))
    aft_CMC_label = paste("treatAfterRpmFold AUC:",sprintf("%.3f",aft_CMC$auc))
    aftList=list(CMC=aft_CMC, Input=aft_input)
    names(aftList) <- c(aft_CMC_label,aft_input_label)
    aft_plot<-ggroc(aftList,size=0.8)
    aft_plot<-aft_plot+geom_abline(intercept = 1,slope = 1,col="black",size=0.6,linetype=2)+theme_classic()+theme(legend.position="top",legend.title=element_blank(),panel.background=element_rect(fill="white",color="black"))
    aft_plot<-aft_plot+ guides(colour = guide_legend(nrow = 2))
    aft_plot<-aft_plot+ annotate("text", x = .5, y = .5,
                              label = paste("Method: ",roc.test(aft_input, aft_CMC)$method,"\np.value: ",signif(roc.test(aft_input, aft_CMC)$p.value,3),sep=""),
                              size = 3)
    print(paste(outFile_prefix,"_aft_plot.pdf",sep=""))
    pdf(paste(outFile_prefix,"_aft_plot.pdf",sep=""),width =4 ,height = 3)
    aft_plot
    invisible(dev.off())
    print(paste(outFile_prefix,"_aft_plot.png",sep=""))
    png(paste(outFile_prefix,"_aft_plot.png",sep=""))
    aft_plot
    invisible(dev.off())

    roc.test(aft_input, aft_CMC)


    cat("=====================treatStopMeanFold/controlStopMeanFold=========================\n")
    #treatStopMeanFold controlStopMeanFold
    controlStopMeanFold_input=roc(ROC_data$class,ROC_data$ctrlStopMeanFold,smooth=FALSE, print.auc=TRUE, col="#e41a1c", plot = FALSE, print.thres="best", print.thres.best.method="youden", levels = c(0,1), direction='<',auc=T, ci=T)
    pdf(paste(outFile_prefix, '_roc_treatStopMeanFold.pdf', sep=""))
    treatStopMeanFold_CMC=roc(ROC_data$class,ROC_data$treatStopMeanFold,smooth=FALSE, print.auc=TRUE, col="#e41a1c", plot = TRUE, print.thres="best", print.thres.best.method="youden", levels = c(0,1), direction='<',auc=T, ci=T)
    invisible(dev.off())
    controlStopMeanFold_label = paste("controlStopMeanFold AUC:",sprintf("%.3f",controlStopMeanFold_input$auc))
    treatStopMeanFold_label = paste("treatStopMeanFold AUC:",sprintf("%.3f",treatStopMeanFold_CMC$auc))
    StopMeanFoldList=list(CMC=treatStopMeanFold_CMC, Input=controlStopMeanFold_input)
    names(StopMeanFoldList) <- c(treatStopMeanFold_label,controlStopMeanFold_label)
    StopMeanFold_plot<-ggroc(StopMeanFoldList,size=0.8)
    StopMeanFold_plot<-StopMeanFold_plot+geom_abline(intercept = 1,slope = 1,col="black",size=0.6,linetype=2)+theme_classic()+theme(legend.position="top",legend.title=element_blank(),panel.background=element_rect(fill="white",color="black"))
    StopMeanFold_plot<-StopMeanFold_plot+ guides(colour = guide_legend(nrow = 2))
    StopMeanFold_plot<-StopMeanFold_plot+ annotate("text", x = .5, y = .5,
                              label = paste("Method: ",roc.test(controlStopMeanFold_input, treatStopMeanFold_CMC)$method,"\np.value: ",signif(roc.test(controlStopMeanFold_input, treatStopMeanFold_CMC)$p.value,3),sep=""),
                              size = 3)
    print(paste(outFile_prefix,"_StopMeanFold_plot.pdf",sep=""))
    pdf(paste(outFile_prefix,"_StopMeanFold_plot.pdf",sep=""),width =4 ,height = 3)
    StopMeanFold_plot
    invisible(dev.off())

    print(paste(outFile_prefix,"_StopMeanFold_plot.png",sep=""))
    png(paste(outFile_prefix,"_StopMeanFold_plot.png",sep=""))
    StopMeanFold_plot
    invisible(dev.off())

    roc.test(controlStopMeanFold_input, treatStopMeanFold_CMC)


    cat("=====================treatStopRatio/controlStopRatio=========================\n")
    #treatStopRatio controlStopRatio
    stopratio_input=roc(ROC_data$class,ROC_data$ctrlStopRatio,smooth=FALSE, print.auc=TRUE, col="#e41a1c", plot = FALSE, print.thres="best", print.thres.best.method="youden", levels = c(0,1), direction='<',auc=T, ci=T)
    pdf(paste(outFile_prefix, '_roc_treatStopRatio.pdf', sep=""))
    stopratio_CMC=roc(ROC_data$class,ROC_data$treatStopRatio,smooth=FALSE, print.auc=TRUE, col="#e41a1c", plot = TRUE, print.thres="best", print.thres.best.method="youden", levels = c(0,1), direction='<',auc=T, ci=T)
    invisible(dev.off())
    stopratio_input_label = paste("controlStopRatio AUC:",sprintf("%.3f",stopratio_input$auc))
    stopratio_CMC_label = paste("treatStopRatio AUC:",sprintf("%.3f",stopratio_CMC$auc))
    stopratioList=list(CMC=stopratio_CMC, Input=stopratio_input)
    names(stopratioList) <- c(stopratio_CMC_label,stopratio_input_label)
    stopratio_plot<-ggroc(stopratioList,size=0.8)
    stopratio_plot<-stopratio_plot+geom_abline(intercept = 1,slope = 1,col="black",size=0.6,linetype=2)+theme_classic()+theme(legend.position="top",legend.title=element_blank(),panel.background=element_rect(fill="white",color="black"))
    stopratio_plot<-stopratio_plot+ guides(colour = guide_legend(nrow = 2))
    stopratio_plot<-stopratio_plot+ annotate("text", x = .5, y = .5,
                              label = paste("Method: ",roc.test(stopratio_input, stopratio_CMC)$method,"\np.value: ",signif(roc.test(stopratio_input, stopratio_CMC)$p.value,3),sep=""),
                              size = 3)
    print(paste(outFile_prefix,"_stopratio_plot.pdf",sep=""))
    pdf(paste(outFile_prefix,"_stopratio_plot.pdf",sep=""),width =4 ,height = 3)
    stopratio_plot
    invisible(dev.off())
    print(paste(outFile_prefix,"_stopratio_plot.png",sep=""))
    png(paste(outFile_prefix,"_stopratio_plot.png",sep=""))
    stopratio_plot
    invisible(dev.off())

    roc.test(stopratio_input, stopratio_CMC)


    cat("=====================Ratio: preRpmFoldRatio, afterRpmFoldRatio, stopRatioFC,treatStopMeanFoldRatio=========================\n")

    mycol = brewer.pal(8, "Set2")[c(1:4)]

    #ratio: preRpmFoldRatio, afterRpmFoldRatio, stopRatioFC,treatStopMeanFoldRatio
    pdf(paste(outFile_prefix, '_roc_preRpmFoldRatio.pdf', sep=""))
    preRpmFoldRatio=roc(ROC_data$class,ROC_data$preRpmFoldRatio,smooth=FALSE, print.auc=TRUE, col="#e41a1c", plot = TRUE, print.thres="best", print.thres.best.method="youden", levels = c(0,1), direction='<',auc=T, ci=T)
    invisible(dev.off())
    pdf(paste(outFile_prefix, '_roc_afterRpmFoldRatio.pdf', sep=""))
    afterRpmFoldRatio=roc(ROC_data$class,ROC_data$afterRpmFoldRatio,smooth=FALSE, print.auc=TRUE, col="#e41a1c", plot = TRUE, print.thres="best", print.thres.best.method="youden", levels = c(0,1), direction='<',auc=T, ci=T)
    invisible(dev.off())
    pdf(paste(outFile_prefix, '_roc_stopRatioFC.pdf', sep=""))
    stopRatioFC=roc(ROC_data$class,ROC_data$stopRatioFC,smooth=FALSE, print.auc=TRUE, col="#e41a1c", plot = TRUE, print.thres="best", print.thres.best.method="youden", levels = c(0,1), direction='<',auc=T, ci=T)
    invisible(dev.off())
    pdf(paste(outFile_prefix, '_roc_treatStopMeanFoldRatio.pdf', sep=""))
    treatStopMeanFoldRatio=roc(ROC_data$class,ROC_data$treatStopMeanFoldRatio,smooth=FALSE, print.auc=TRUE, col="#e41a1c", plot = TRUE, print.thres="best", print.thres.best.method="youden", levels = c(0,1), direction='<',auc=T, ci=T)
    invisible(dev.off())
    preRpmFoldRatio_label = paste("preRpmFoldRatio AUC:",sprintf("%.3f",preRpmFoldRatio$auc))
    afterRpmFoldRatio_label = paste("afterRpmFoldRatio AUC:",sprintf("%.3f",afterRpmFoldRatio$auc))
    stopRatioFC_label = paste("stopRatioFC AUC:",sprintf("%.3f",stopRatioFC$auc))
    treatStopMeanFoldRatio_label = paste("treatStopMeanFoldRatio AUC:",sprintf("%.3f",treatStopMeanFoldRatio$auc))
    ratioList=list(preRpmFoldRatio = preRpmFoldRatio, afterRpmFoldRatio = afterRpmFoldRatio,  stopRatioFC = stopRatioFC, treatStopMeanFoldRatio = treatStopMeanFoldRatio)
    names(ratioList) <- c(preRpmFoldRatio_label, afterRpmFoldRatio_label, stopRatioFC_label, treatStopMeanFoldRatio_label)
    ratio_plot<-ggroc(ratioList,size=0.8)
    ratio_plot<-ratio_plot+geom_abline(intercept = 1,slope = 1,col="black",size=0.6,linetype=2)+theme_classic()+theme(legend.position="top",legend.title=element_blank(),panel.background=element_rect(fill="white",color="black"))
    ratio_plot<-ratio_plot+ guides(colour = guide_legend(nrow = 2, ncol = 2))+scale_color_manual(values = mycol)
    print(paste(outFile_prefix,"_ratio_plot.pdf",sep=""))
    pdf(paste(outFile_prefix,"_ratio_plot.pdf",sep=""),width =9 ,height = 9)
    ratio_plot
    invisible(dev.off())
    print(paste(outFile_prefix,"_ratio_plot.png",sep=""))
    png(paste(outFile_prefix,"_ratio_plot.png",sep=""))
    ratio_plot
    invisible(dev.off())


    cat("=====================six variables: treatPreRpmFold, preRpmFoldRatio, treatAfterRpmFold, afterRpmFoldRatio, stopRatioFC,treatStopRatio=========================\n")

    mycol = brewer.pal(8, "Dark2")[c(1:6)]

    #ratio: preRpmFoldRatio, afterRpmFoldRatio, stopRatioFC,treatStopRatio
    pdf(paste(outFile_prefix, '_roc_preRpmFoldRatio.pdf', sep=""))
    preRpmFoldRatio=roc(ROC_data$class,ROC_data$preRpmFoldRatio,smooth=FALSE, print.auc=TRUE, col="#e41a1c", plot = TRUE, print.thres="best", print.thres.best.method="youden", levels = c(0,1), direction='<',auc=T, ci=T)
    invisible(dev.off())
    pdf(paste(outFile_prefix, '_roc_afterRpmFoldRatio.pdf', sep=""))
    afterRpmFoldRatio=roc(ROC_data$class,ROC_data$afterRpmFoldRatio,smooth=FALSE, print.auc=TRUE, col="#e41a1c", plot = TRUE, print.thres="best", print.thres.best.method="youden", levels = c(0,1), direction='<',auc=T, ci=T)
    invisible(dev.off())
    pdf(paste(outFile_prefix, '_roc_stopRatioFC.pdf', sep=""))
    stopRatioFC=roc(ROC_data$class,ROC_data$stopRatioFC,smooth=FALSE, print.auc=TRUE, col="#e41a1c", plot = TRUE, print.thres="best", print.thres.best.method="youden", levels = c(0,1), direction='<',auc=T, ci=T)
    invisible(dev.off())
    pdf(paste(outFile_prefix, '_roc_treatStopRatio.pdf', sep=""))
    treatStopRatio=roc(ROC_data$class,ROC_data$treatStopRatio,smooth=FALSE, print.auc=TRUE, col="#e41a1c", plot = TRUE, print.thres="best", print.thres.best.method="youden", levels = c(0,1), direction='<',auc=T, ci=T)
    invisible(dev.off())
    preRpmFoldRatio_label = paste("preRpmFoldRatio AUC:",sprintf("%.3f",preRpmFoldRatio$auc))
    afterRpmFoldRatio_label = paste("afterRpmFoldRatio AUC:",sprintf("%.3f",afterRpmFoldRatio$auc))
    stopRatioFC_label = paste("stopRatioFC AUC:",sprintf("%.3f",stopRatioFC$auc))
    treatStopRatio_label = paste("treatStopRatio AUC:",sprintf("%.3f",treatStopRatio$auc))
    ratioList=list(treatPreRpmFold = pre_CMC, preRpmFoldRatio = preRpmFoldRatio, treatAfterRpmFold= aft_CMC, afterRpmFoldRatio = afterRpmFoldRatio,  stopRatioFC = stopRatioFC, treatStopRatio = treatStopRatio)
    names(ratioList) <- c(pre_CMC_label, preRpmFoldRatio_label, aft_CMC_label, afterRpmFoldRatio_label, stopRatioFC_label, treatStopRatio_label)
    six_variables_ratio_plot<-ggroc(ratioList,size=0.8)
    six_variables_ratio_plot<-six_variables_ratio_plot+geom_abline(intercept = 1,slope = 1,col="black",size=0.6,linetype=2)+
    theme_classic()+
    theme(legend.position="top",
      legend.key.width= unit(0.3, 'cm'),
      legend.text = element_text(size=7),
      legend.title=element_blank(),
      panel.background=element_rect(fill="white",color="black"))
    six_variables_ratio_plot<-six_variables_ratio_plot+ guides(colour = guide_legend(nrow = 2, ncol = 3))+scale_color_manual(values = mycol)
    print(paste(outFile_prefix,"_six_variables_ratio_plot.pdf",sep=""))
    pdf(paste(outFile_prefix,"_six_variables_plot.pdf",sep=""),width =6 ,height = 6)
    six_variables_ratio_plot
    invisible(dev.off())

    pdf(paste(outFile_prefix,"_roc_summary.pdf",sep=""),width=16,height=8)
    plot_grid(pre_plot, aft_plot, StopMeanFold_plot, stopratio_plot, ratio_plot,six_variables_ratio_plot, align = "hv",labels = c('A', 'B','C','D','E','F'))
    invisible(dev.off())

    invisible(lapply(seq_along(ROC_data_sel[,-length(ROC_data_sel)]),function(x){
        assign(colnames(ROC_data_sel)[x],roc(ROC_data$class,ROC_data_sel[,x],smooth=FALSE, print.auc=TRUE, col="#e41a1c", plot = FALSE, print.thres="best", print.thres.best.method="youden", levels = c(0,1), direction='<',auc=T, ci=T),envir = .GlobalEnv)
        assign(paste(colnames(ROC_data_sel)[x],"_thres",sep=""),coords(get(colnames(ROC_data_sel)[x]), "best",transpose = TRUE)[1],envir = .GlobalEnv)
    }))

    cat("=====================pre_CMC=========================\n")
    # treatPreRpmFold_thres<-coords(pre_CMC, "best",transpose = TRUE)[1]
    coords(pre_CMC, "best", ret = "all", transpose = TRUE)

    cat("=====================aft_CMC=========================\n")
    # treatAfterRpmFold_thres<-coords(aft_CMC, "best",transpose = TRUE)[1]
    coords(aft_CMC, "best", ret = "all", transpose = TRUE)

    cat("=====================treatStopMeanFold_CMC=========================\n")
    # treatStopMeanFold_thres<-coords(treatStopMeanFold_CMC, "best",transpose = TRUE)[1]
    coords(treatStopMeanFold_CMC, "best", ret = "all", transpose = TRUE)

    cat("=====================preRpmFoldRatio=========================\n")
    # preRpmFoldRatio_thres<-coords(preRpmFoldRatio, "best",transpose = TRUE)[1]
    coords(preRpmFoldRatio, "best", ret = "all", transpose = TRUE)

    cat("=====================afterRpmFoldRatio=========================\n")
    # afterRpmFoldRatio_thres<-coords(afterRpmFoldRatio, "best",transpose = TRUE)[1]
    coords(afterRpmFoldRatio, "best", ret = "all", transpose = TRUE)

    cat("=====================treatStopRatio=========================\n")
    # treatStopRatio_thres<-coords(treatStopRatio, "best",transpose = TRUE)[1]
    coords(treatStopRatio, "best", ret = "all", transpose = TRUE)

    cat("=====================stopRatioFC=========================\n")
    # stopRatioFC_thres<-coords(stopRatioFC, "best",transpose = TRUE)[1]
    coords(stopRatioFC, "best", ret = "all", transpose = TRUE)

    cat("=====================treatStopMeanFoldRatio=========================\n")
    # treatStopMeanFoldRatio_thres<-coords(treatStopMeanFoldRatio, "best",transpose = TRUE)[1]
    coords(treatStopMeanFoldRatio, "best", ret = "all", transpose = TRUE)


    thres<-cbind(treatPreRpmFold_thres,treatAfterRpmFold_thres,treatStopMeanFold_thres,treatStopRatio_thres,preRpmFoldRatio_thres, afterRpmFoldRatio_thres, stopRatioFC_thres,treatStopMeanFoldRatio_thres)
    thres
    write.table(thres,paste(outFile_prefix, '_roc_bestthres.txt', sep=""),sep="\t",row.names=FALSE, col.names=FALSE,quote=F)
    write.table(thres,paste(outFile_prefix, '_roc_bestthres_colname.txt', sep=""),sep="\t",row.names=FALSE, col.names=TRUE,quote=F)

    # Values of all the confusion matrix terms were calculated at the optimal threshold
    invisible(lapply(seq_along(ROC_data_sel[,-length(ROC_data_sel)]),function(x){
        assign(paste(colnames(ROC_data_sel)[x],"_TP",sep=""),dim(ROC_data_sel[as.numeric(ROC_data_sel$class)==1 & ROC_data_sel[,x] > get(paste(colnames(ROC_data_sel)[x],"_thres",sep="")), ])[1],envir = .GlobalEnv)#True Positives (TP)
        assign(paste(colnames(ROC_data_sel)[x],"_FP",sep=""),dim(ROC_data_sel[as.numeric(ROC_data_sel$class)==0 & ROC_data_sel[,x] > get(paste(colnames(ROC_data_sel)[x],"_thres",sep="")), ])[1],envir = .GlobalEnv)#False Positives (FP)
        assign(paste(colnames(ROC_data_sel)[x],"_TN",sep=""),dim(ROC_data_sel[as.numeric(ROC_data_sel$class)==0 & ROC_data_sel[,x] <= get(paste(colnames(ROC_data_sel)[x],"_thres",sep="")), ])[1],envir = .GlobalEnv)#True Negatives (TN)
        assign(paste(colnames(ROC_data_sel)[x],"_FN",sep=""),dim(ROC_data_sel[as.numeric(ROC_data_sel$class)==1 & ROC_data_sel[,x] <= get(paste(colnames(ROC_data_sel)[x],"_thres",sep="")), ])[1],envir = .GlobalEnv)#False Negatives (FN)
        assign(paste(colnames(ROC_data_sel)[x],"_TPR",sep=""),round(get(paste(colnames(ROC_data_sel)[x],"_TP",sep=""))/(get(paste(colnames(ROC_data_sel)[x],"_TP",sep=""))+get(paste(colnames(ROC_data_sel)[x],"_FN",sep=""))),3),envir = .GlobalEnv)#sensitivity (true positive rate, TPR)
        assign(paste(colnames(ROC_data_sel)[x],"_TNR",sep=""),round(get(paste(colnames(ROC_data_sel)[x],"_TN",sep=""))/(get(paste(colnames(ROC_data_sel)[x],"_TN",sep=""))+get(paste(colnames(ROC_data_sel)[x],"_FP",sep=""))),3),envir = .GlobalEnv)#specifity (selectivity or true negative rate, TNR)
        assign(paste(colnames(ROC_data_sel)[x],"_FPR",sep=""),round(1-get(paste(colnames(ROC_data_sel)[x],"_TNR",sep="")),3),envir = .GlobalEnv)#False Positive Rate (FPR) (1 - specificit = FP/​N = FP/(TN + FP), FPR)
        assign(paste(colnames(ROC_data_sel)[x],"_FNR",sep=""),round(1-get(paste(colnames(ROC_data_sel)[x],"_TPR",sep="")),3),envir = .GlobalEnv)#False Negative Rate, FNR)
        assign(paste(colnames(ROC_data_sel)[x],"_Prec",sep=""),round(get(paste(colnames(ROC_data_sel)[x],"_TP",sep=""))/(get(paste(colnames(ROC_data_sel)[x],"_TP",sep=""))+get(paste(colnames(ROC_data_sel)[x],"_FP",sep=""))),3),envir = .GlobalEnv)#Precision
        assign(paste(colnames(ROC_data_sel)[x],"_Recall",sep=""),round(get(paste(colnames(ROC_data_sel)[x],"_TP",sep=""))/(get(paste(colnames(ROC_data_sel)[x],"_TP",sep=""))+get(paste(colnames(ROC_data_sel)[x],"_FN",sep=""))),3),envir = .GlobalEnv)#Recall
        assign(paste(colnames(ROC_data_sel)[x],"_ACC",sep=""),round((get(paste(colnames(ROC_data_sel)[x],"_TP",sep=""))+get(paste(colnames(ROC_data_sel)[x],"_TN",sep="")))/(get(paste(colnames(ROC_data_sel)[x],"_TP",sep=""))+get(paste(colnames(ROC_data_sel)[x],"_TN",sep=""))+get(paste(colnames(ROC_data_sel)[x],"_FP",sep=""))+get(paste(colnames(ROC_data_sel)[x],"_FN",sep=""))),3),envir = .GlobalEnv)#accuracy
        assign(paste(colnames(ROC_data_sel)[x],"_F1_score",sep=""),round((2*get(paste(colnames(ROC_data_sel)[x],"_Recall",sep=""))*get(paste(colnames(ROC_data_sel)[x],"_Prec",sep="")))/(get(paste(colnames(ROC_data_sel)[x],"_Recall",sep=""))+get(paste(colnames(ROC_data_sel)[x],"_Prec",sep=""))),3),envir = .GlobalEnv)#F1_score
        assign(paste(colnames(ROC_data_sel)[x],"_eval",sep=""),cbind(get(paste(colnames(ROC_data_sel)[x],"_TP",sep="")),get(paste(colnames(ROC_data_sel)[x],"_FP",sep="")),get(paste(colnames(ROC_data_sel)[x],"_TN",sep="")),get(paste(colnames(ROC_data_sel)[x],"_FN",sep="")),get(paste(colnames(ROC_data_sel)[x],"_TPR",sep="")),get(paste(colnames(ROC_data_sel)[x],"_TNR",sep="")),get(paste(colnames(ROC_data_sel)[x],"_FPR",sep="")),get(paste(colnames(ROC_data_sel)[x],"_FNR",sep="")),get(paste(colnames(ROC_data_sel)[x],"_Prec",sep="")),get(paste(colnames(ROC_data_sel)[x],"_Recall",sep="")),get(paste(colnames(ROC_data_sel)[x],"_ACC",sep="")),get(paste(colnames(ROC_data_sel)[x],"_F1_score",sep=""))),envir = .GlobalEnv)
        tmp<-as.data.frame(get(paste(colnames(ROC_data_sel)[x],"_eval",sep="")))
        names(tmp)<-c(
          paste(colnames(ROC_data_sel)[x],"_TP",sep=""),
          paste(colnames(ROC_data_sel)[x],"_FP",sep=""),
          paste(colnames(ROC_data_sel)[x],"_TN",sep=""),
          paste(colnames(ROC_data_sel)[x],"_FN",sep=""),
          paste(colnames(ROC_data_sel)[x],"_TPR",sep=""),
          paste(colnames(ROC_data_sel)[x],"_TNR",sep=""),
          paste(colnames(ROC_data_sel)[x],"_FPR",sep=""),
          paste(colnames(ROC_data_sel)[x],"_FNR",sep=""),
          paste(colnames(ROC_data_sel)[x],"_Prec",sep=""),
          paste(colnames(ROC_data_sel)[x],"_Recall",sep=""),
          paste(colnames(ROC_data_sel)[x],"_ACC",sep=""),
          paste(colnames(ROC_data_sel)[x],"_F1_score",sep="")
          )
        assign(paste(colnames(ROC_data_sel)[x],"_eval",sep=""), tmp, envir = .GlobalEnv)
        }))


    confusion_matrix_and_indicators<-as.data.frame(rbind(
        t(as.data.frame(treatPreRpmFold_eval)),
        t(as.data.frame(treatAfterRpmFold_eval)),
        t(as.data.frame(treatStopMeanFold_eval)),
        t(as.data.frame(preRpmFoldRatio_eval)),
        t(as.data.frame(afterRpmFoldRatio_eval)),
        t(as.data.frame(treatStopRatio_eval)),
        t(as.data.frame(stopRatioFC_eval)),
        t(as.data.frame(treatStopMeanFoldRatio_eval))))
    colnames(confusion_matrix_and_indicators)<-"value_or_percentage"
    confusion_matrix_and_indicators_arrange<-arrange(confusion_matrix_and_indicators,desc(value_or_percentage))
    write.table(confusion_matrix_and_indicators,paste(outFile_prefix, '_roc_confusion_matrix_and_indicators.txt', sep=""),sep="\t",row.names=TRUE, col.names=TRUE,quote=F)
    write.table(confusion_matrix_and_indicators_arrange,paste(outFile_prefix, '_roc_confusion_matrix_and_indicators_arrange.txt', sep=""),sep="\t",row.names=TRUE, col.names=TRUE,quote=F)


    best_eval<-rownames(confusion_matrix_and_indicators_arrange)[str_detect(rownames(confusion_matrix_and_indicators_arrange),"_F1_score")][1]
    assign("best_eval_t_df",as.data.frame(t(as.data.frame(get(str_replace(best_eval,"_F1_score","_eval"))))))
    colnames(best_eval_t_df)<-"value_or_percentage"

    #show afterRpmFoldRatio_eval_t_df as pdf table
    tt3 <- ttheme_minimal(
      core=list(bg_params = list(fill = blues9[1:4], col=NA),
                fg_params=list(fontface=3)),
      colhead=list(fg_params=list(col="navyblue", fontface=4L)),
      rowhead=list(fg_params=list(col="orange", fontface=3L)))

    pdf(paste(outFile_prefix, '_roc_best_evaluation.pdf', sep=""), width = 7, height = 7) # Open a new pdf file
    grid.arrange(
      tableGrob(best_eval_t_df, theme=tt3),
      nrow=1)
    invisible(dev.off()) # Close the file

    #read hg38_human_chr21_rRNA_known_pseudoU_SingleSites.bed evidence
    evidence<-read.table(rRNAfile,head=F)#"hg38_human_chr21_rRNA_known_pseudoU_SingleSites.bed"
    colnames(evidence)<-c("chrom","chromStart","chromEnd","rRNA_anno","score","strand")
    evidence$rRNA_uniq_id<-paste(evidence$chrom,evidence$chromStart,evidence$chromEnd,evidence$strand,sep="_")

    #read hg38.psiU.SingleSites.bed
    evidence2<-read.table(rRNAfile2,head=F)#"hg38.psiU.SingleSites.bed"
    colnames(evidence2)<-c("chrom","chromStart","chromEnd","rRNA_anno","score","strand")
    evidence2$rRNA_uniq_id<-paste(evidence2$chrom,evidence2$chromStart,evidence2$chromEnd,evidence2$strand,sep="_")

    #filt by stopRatioFC_thres
    to_filt<-read.table(filtfile,head=F)
    colnames(to_filt)<-c("chrom",#1
        "chromStart",#2
        "chromEnd",#3
        "name",#4
        "foldChange",#5
        "strand",#6
        "geneName",#7
        "geneStart",#8
        "geneEnd",#9
        "base",#10
        "treatPval",#11
        "ctrlPval",#12
        "minusPval",#13
        "treatStopNum",#14
        "treatStopRPM",#15
        "treatPreStopNum",#16
        "treatAfterStopNum",#17
        "treatReadthroughNum",#18
        "ctrlStopNum",#19
        "ctrlStopRPM",#20
        "ctrlPreStopNum",#21
        "ctrlAfterStopNum",#22
        "ctrlReadthroughNum",#23
        "stopRpmFC",#24
        "treatPreRpmFold",#25
        "ctrlPreRpmFold",#26
        "preRpmFoldRatio",#27
        "treatAfterRpmFold",#28
        "ctrlAfterRpmFold",#29
        "afterRpmFoldRatio",#30
        "treatStopRatio",#31
        "ctrlStopRatio",#32
        "stopRatioFC",#33
        "treatStopMeanNum",#34
        "treatStopMeanFold",#35
        "ctrlStopMeanNum",#36
        "ctrlStopMeanFold",#37
        "treatStopMeanFoldRatio",#38
        "extendSeq")#39
    to_filt$base<-"T"

    to_filt$roc_class<-as.numeric(to_filt[,str_replace(best_eval,"_F1_score","")]>get(str_replace(best_eval,"_F1_score","_thres")))
    to_filt$roc_class<-str_replace(as.character(to_filt$roc_class),"0","non-psi")
    to_filt$roc_class<-str_replace(as.character(to_filt$roc_class),"1","psi")
    table(to_filt$roc_class)
    write.table(to_filt,paste(outFile_prefix, '_roc_total_prediction.txt', sep=""),sep="\t",row.names=FALSE, col.names=TRUE,quote=F)
    final_pred<-to_filt[to_filt[,str_replace(best_eval,"_F1_score","")]>get(str_replace(best_eval,"_F1_score","_thres")),]
    write.table(final_pred,paste(outFile_prefix, '_roc_psi_prediction.txt', sep=""),sep="\t",row.names=FALSE, col.names=TRUE,quote=F)
    write.table(final_pred,paste(outFile_prefix, '_roc_psi_prediction.bed', sep=""),sep="\t",row.names=FALSE, col.names=FALSE,quote=F)

    #known_data miss/hit hg38_human_chr21_rRNA_known_pseudoU_SingleSites.bed
    ROC_data$rRNA_uniq_id<-paste(ROC_data$chrom,ROC_data$chromStart,ROC_data$chromEnd,ROC_data$strand,sep="_")
    ROC_data_evidence<-ROC_data %>% left_join(evidence,by=c("rRNA_uniq_id"="rRNA_uniq_id")) %>% filter(!is.na(rRNA_anno))
    write.csv(ROC_data_evidence,paste(outFile_prefix,"_hg38_human_chr21_rRNA_known_pseudoU_SingleSites_rtsSeeker_hit.csv",sep=""))
    ROC_data_no_evidence<-evidence %>% left_join(ROC_data,by=c("rRNA_uniq_id"="rRNA_uniq_id")) %>% filter(is.na(class))
    write.csv(ROC_data_no_evidence,paste(outFile_prefix,"_hg38_human_chr21_rRNA_known_pseudoU_SingleSites_rtsSeeker_miss.csv",sep=""))
    recovery<-paste(round(length(unique(ROC_data_evidence$rRNA_uniq_id))/length(unique(evidence$rRNA_uniq_id))*100,2),"%",sep="")
    cat("rtsSeeker recover (hg38_human_chr21_rRNA_known_pseudoU_SingleSites.bed)",recovery,"rRNA psi sites in all known chrom21\n")

    #known_data miss/hit hg38.psiU.SingleSites.bed
    ROC_data$rRNA_uniq_id<-paste(ROC_data$chrom,ROC_data$chromStart,ROC_data$chromEnd,ROC_data$strand,sep="_")
    ROC_data_evidence2<-ROC_data %>% left_join(evidence2,by=c("rRNA_uniq_id"="rRNA_uniq_id")) %>% filter(!is.na(rRNA_anno))
    write.csv(ROC_data_evidence2,paste(outFile_prefix,"_hg38.psiU.SingleSites.bed_rtsSeeker_hit.csv",sep=""))
    ROC_data_no_evidence2<-evidence2 %>% left_join(ROC_data,by=c("rRNA_uniq_id"="rRNA_uniq_id")) %>% filter(is.na(class))
    write.csv(ROC_data_no_evidence2,paste(outFile_prefix,"_hg38.psiU.SingleSites.bed_rtsSeeker_miss.csv",sep=""))
    recovery<-paste(round(length(unique(ROC_data_evidence2$rRNA_uniq_id))/length(unique(evidence2$rRNA_uniq_id))*100,2),"%",sep="")
    cat("rtsSeeker recover (hg38.psiU.SingleSites.bed.bed)",recovery,"rRNA psi sites in all known chrom21\n")

    #final_pred miss/hit
    final_pred$rRNA_uniq_id<-paste(final_pred$chrom,final_pred$chromStart,final_pred$chromEnd,final_pred$strand,sep="_")
    final_pred_evidence<-final_pred %>% left_join(evidence,by=c("rRNA_uniq_id"="rRNA_uniq_id")) %>% filter(!is.na(rRNA_anno))
    write.csv(final_pred_evidence,paste(outFile_prefix,"_hg38_human_chr21_rRNA_known_pseudoU_SingleSites_rtsSeeker_roc_stopRatioFC_thres_hit.csv",sep=""))
    final_pred_no_evidence<-evidence %>% left_join(final_pred,by=c("rRNA_uniq_id"="rRNA_uniq_id")) %>% filter(is.na(base))
    write.csv(final_pred_no_evidence,paste(outFile_prefix,"_hg38_human_chr21_rRNA_known_pseudoU_SingleSites_rtsSeeker_roc_stopRatioFC_thres_miss.csv",sep=""))
    recovery<-paste(round(length(unique(final_pred_evidence$rRNA_uniq_id))/length(unique(evidence$rRNA_uniq_id))*100,2),"%",sep="")
    cat("rtsSeeker+ROC recover (hg38_human_chr21_rRNA_known_pseudoU_SingleSites.bed)",recovery,"rRNA psi sites in all known chrom21\n")

search psi sites
********************
Correspond to rtsSeeker ``-U/--psi`` argument.

skip read with soft clip
*************************
Default is false, correspond to rtsSeeker ``-S/--soft`` argument.

normalized reads to the locus number
********************************************
Correspond to rtsSeeker ``-n/--norm`` argument.

keep duplication
*************************
Default is false, correspond to rtsSeeker ``-c/--collapser`` argument.


Detail
-------
Detailed argume
nt for ``rtsSeeker`` output information control.

.. code:: bash

    Usage:  rtsSeeker [options] --fa <genome seq> --fai <fai file> --gene <bed12 gene file> --treat <treat alignments> --input <input alignments>
    rtsSeeker: for seeking the reverse transcriptase stop sites
    [options]
    -v/--verbose                   : verbose information
    -V/--version                   : rtsSeeker version
    -h/--help                      : help informations
    -U/--psi                       : search psi sites
    -s/--strand                    : strand-specific[default=false]
    -S/--soft                      : skip read with soft-clip[default=false]
    -n/--norm                      : normalized reads to the locus number
    -p/--pval                      : p-value[default=1.5]
    -c/--collapser                 : keep duplication, deault is false
    -M/--model <int>               : model[0 for genome, 1 for gene, 2 for both, default=1]
    --treat <string>               : treatment file<BAM format>
    --input <string>               : input file<BAM format>
    --gene <string>                : gene file <BED12 format>
    -o/--outfile <string>          : output file
    -t/--min-tag <double>          : minimum tag number for each psi, default>=5.0 read
    -r/--rpm <double>              : minimum rpm value for each psi, default>=0.05
    -f/--fold <int>                : minimum fold-change[default>=1.0]
    -m/--mfold <int>               : minimum fold-change for stop-tag/mean-tag
    -w/--window <int>              : window size around the rts position[default=50]
    -l/--min-len <int>             : minimum length of reads, default=15
    -L/--max-len <int>             : maximum length of reads, default=1000


Identify Ψ-sites using different approaches
---------------------------------------------

(1) SVM: Support Vector Machine
****************************************
If ``SVM`` QT wideget is clicked and popouped:

-  press ``BUILD`` buttor will run ``svm.sh`` and ``svm_build_totalRNA.r`` to generate file with ``_SVM_model.RData`` suffix, which is the support vector model built by SVM algorithm.
-  press ``PREDICT`` buttor will run ``svm_predict_totalRNA.sh`` and ``svm_predict_totalRNA.r`` to generate file with ``_svm_psi_prediction.bed`` suffix, which is the candidate Ψ-sites predicted by SVM algorithm.

.. image:: /images/SVM.png

``svm.sh``

.. code:: bash

    #! /bin/bash

    if [[ $# -ne 4 ]]
    then
        echo 'Usage: '$0 ' input.bam treatment.bam genome.fa out.bed'
        exit 1
    fi

    input=$1
    treatment=$2
    genome=$3
    out=$4
    outFileName=${out%.bed}

    echo "rtsSeeker SVM mode: rtsSeeker start..."
    if [ ! -f "${genome}.fai"  ]
    then
    echo "generate faifile"
    samtools faidx $genome
    else
    echo "faifile exist"
    fi

    echo "$(dirname "$0")/rtsSeeker --fa $genome --fai ${genome}.fai --treat $treatment --input $input -p 1.5 -t 5 -r 0.05 -M 1 --gene $(dirname "$0")/hg38.gencode.v30.tRNA.refseqNcRNA.geneAnno.bed12 -f 1 -m 0 -s -w 20 -o ${outFileName}_svm_all.bed 2>${outFileName}_svm_all_rtsSeeker.log" > ${outFileName}_svm_all_rtsSeeker.cmd # no -c remove duplicate
    $(dirname "$0")/rtsSeeker --fa $genome --fai ${genome}.fai --treat $treatment --input $input -p 1.5 -t 5 -r 0.05 -M 1 --gene $(dirname "$0")/hg38.gencode.v30.tRNA.refseqNcRNA.geneAnno.bed12 -f 1 -m 0 -s -w 20 -o ${outFileName}_svm_all.bed 2>${outFileName}_svm_all_rtsSeeker.log # no -c remove duplicate
    cat ${outFileName}_svm_all.bed > $out

    awk 'FS=OFS="\t" {if($1~/^chr[0-9|a-z|A-Z]*$/ && $10=="T" && $14>10){print $0}}' ${outFileName}_svm_all.bed > ${outFileName}_svm_filt_totalRNA.bed
    bedtools intersect -a ${outFileName}_svm_filt_totalRNA.bed -b $(dirname "$0")/hg38_human_chr21_rRNA_known_pseudoU_SingleSites.bed -s > ${outFileName}_knowpse.bed
    bedtools intersect -a ${outFileName}_svm_filt_totalRNA.bed -b $(dirname "$0")/rrna_chr21.bed -s > ${outFileName}_rrna.bed
    bedtools intersect -a ${outFileName}_rrna.bed -b $(dirname "$0")/hg38_human_chr21_rRNA_known_pseudoU_SingleSites.bed -s -v > ${outFileName}_notpse.bed
    cat ${outFileName}_knowpse.bed |awk 'FS=OFS="\t" {print $0,"1"}' > ${outFileName}_knowpse.txt
    cat ${outFileName}_notpse.bed |awk 'FS=OFS="\t" {print $0,"0"}' > ${outFileName}_notpse.txt
    echo "getting roc_plot.txt"
    cat ${outFileName}_knowpse.txt ${outFileName}_notpse.txt > ${outFileName}_svm_plot.txt

    nohup Rscript $(dirname "$0")/svm_build_totalRNA.r -f ${outFileName}_svm_plot.txt -k ${outFileName}_svm_filt_totalRNA.bed -r $(dirname "$0")/hg38_human_chr21_rRNA_known_pseudoU_SingleSites.bed -s $(dirname "$0")/hg38.psiU.SingleSites.bed  -o ${outFileName} > ${outFileName}_svm_evaluation_totalRNA.log 2>&1 &
    wait

    cat ${outFileName}_svm_evaluation_totalRNA.log
    mupdf-x11 ${outFileName}_six_variables_rRNA_violinplot.pdf &> /dev/null
    mupdf-x11 ${outFileName}_SVM_roc_test_data_plot.pdf &> /dev/null
    mupdf-x11 ${outFileName}_svm_evaluation.pdf &> /dev/null

    echo "total RNA: SVM program is done!"
    echo -e "total RNA: SVM result in $(dirname ${outFileName}_svm_psi_prediction.bed)"

``svm_build_totalRNA.r``

.. code:: R

    #!/usr/bin/env Rscript

    options(warn=-1)
    suppressMessages(library("optparse"))
    suppressMessages(library("e1071"))
    suppressMessages(library("dplyr"))
    suppressMessages(library("stringr"))
    suppressMessages(library("gridExtra"))
    suppressMessages(library("cowplot"))
    suppressMessages(library("pROC"))
    suppressMessages(library("ggplot2"))
    suppressMessages(library("ggpubr"))
    suppressMessages(library("reshape2"))
    suppressMessages(library("RColorBrewer"))
    suppressMessages(library("openxlsx"))
    suppressMessages(library("caTools"))
    suppressMessages(library("factoextra"))
    suppressMessages(library("ggpol"))


    options(warn=-1)

    option_list = list(
      make_option(c("-f", "--svmfile"), type="character", default=NULL,
                  help="ROC of single sites [file]", metavar="character"),
      make_option(c("-k", "--filtfile"), type="character", default=NULL,
                  help="filted file of single sites [file]", metavar="character"),
      make_option(c("-r", "--rRNAfile"), type="character", default=NULL,
                  help="rRNA of single sites [file]", metavar="character"),
      make_option(c("-s", "--rRNAfile2"), type="character", default=NULL,
                  help="hg38.psiU.SingleSites.bed [file]", metavar="character"),
      make_option(c("-o", "--outfile_prefix"), type="character", default=NULL,
                  help="output file name [default= %default]", metavar="character")
    );
    opt_parser = OptionParser(option_list=option_list);
    opt = parse_args(opt_parser);

    if (is.null(opt$svmfile)|| is.null(opt$filtfile) || is.null(opt$rRNAfile) || is.null(opt$rRNAfile2) || is.null(opt$outfile_prefix) ){
      print_help(opt_parser);
      stop("Please provide -f svmfile, -k filtfile, -r rRNAfile, -s rRNAfile2, and -o outfile_prefix option", call.=FALSE);
    }

    svmfile = opt$svmfile
    filtfile = opt$filtfile
    rRNAfile = opt$rRNAfile
    rRNAfile2 = opt$rRNAfile2
    outFile_prefix = opt$outfile_prefix

    print(svmfile)
    print(filtfile)
    print(rRNAfile)
    print(rRNAfile2)
    print(outFile_prefix)


    #load data
    cat("\n\n=====================Load data to perform classification model (factor as input)=====================\n")
    SVM_data<-read.table(svmfile)
    colnames(SVM_data)<-c("chrom",#1
      "chromStart",#2
      "chromEnd",#3
      "name",#4
      "foldChange",#5
      "strand",#6
      "geneName",#7
      "geneStart",#8
      "geneEnd",#9
      "base",#10
      "treatPval",#11
      "ctrlPval",#12
      "minusPval",#13
      "treatStopNum",#14
      "treatStopRPM",#15
      "treatPreStopNum",#16
      "treatAfterStopNum",#17
      "treatReadthroughNum",#18
      "ctrlStopNum",#19
      "ctrlStopRPM",#20
      "ctrlPreStopNum",#21
      "ctrlAfterStopNum",#22
      "ctrlReadthroughNum",#23
      "stopRpmFC",#24
      "treatPreRpmFold",#25 *
      "ctrlPreRpmFold",#26
      "preRpmFoldRatio",#27 *
      "treatAfterRpmFold",#28 *
      "ctrlAfterRpmFold",#29
      "afterRpmFoldRatio",#30 *
      "treatStopRatio",#31 *
      "ctrlStopRatio",#32
      "stopRatioFC",#33 *
      "treatStopMeanNum",#34
      "treatStopMeanFold",#35
      "ctrlStopMeanNum",#36
      "ctrlStopMeanFold",#37
      "treatStopMeanFoldRatio",#38
      "extendSeq",#39
      "class")#40
    SVM_data$base<-"T"

    #rRNA-psi-non-psi visualization
    SVM_data_sel<-SVM_data %>% select(treatPreRpmFold,preRpmFoldRatio,treatAfterRpmFold,afterRpmFoldRatio,treatStopRatio,stopRatioFC,class)
    SVM_data_melt<-melt(SVM_data_sel,id.vars = "class")
    SVM_data_melt$class<-str_replace(as.character(SVM_data_melt$class),"0","non-psi")
    SVM_data_melt$class<-str_replace(as.character(SVM_data_melt$class),"1","psi")

    data_summary <- function(x) {
       m <- mean(x)
       ymin <- m-sd(x)
       ymax <- m+sd(x)
       return(c(y=m,ymin=ymin,ymax=ymax))
    }

    my_comparisons <- list( c("non-psi", "psi") )

    #treatPreRpmFold
    treatPreRpmFold<-SVM_data_melt%>%filter(str_detect(.$variable,"^treatPreRpmFold$"))
    treatPreRpmFold_bp <- ggplot(treatPreRpmFold, aes(x=class, y=log2(value), fill=class)) +
      stat_boxplot(geom = "errorbar",
                   width = 0.15) +
      geom_violin(trim=FALSE,alpha=0.8)+
      stat_summary(fun.data=data_summary,geom="crossbar", width=0.2,color="black")+
    labs(x="group", y = "log2(treatPreRpmFold)")
    treatPreRpmFold_bp<-treatPreRpmFold_bp + scale_fill_manual(values=c(brewer.pal(9,"Set1")[2],brewer.pal(9,"Set1")[1])) + theme_classic()+
    theme(text=element_text(size=8),legend.position = "none") + font("xy.text", size = 8)+
    stat_compare_means(method = "t.test",comparisons = my_comparisons,label.y = max(log2(treatPreRpmFold$value))+abs(quantile(log2(treatPreRpmFold$value))[1]))+
    scale_y_continuous(breaks = scales::pretty_breaks(n = 5))+
    scale_y_continuous(expand = expansion(mult = c(0, 0.1)))+
    scale_fill_manual(values=c("#129a92","#1e95d4"))

    #treatAfterRpmFold
    treatAfterRpmFold<-SVM_data_melt%>%filter(str_detect(.$variable,"^treatAfterRpmFold$"))
    treatAfterRpmFold_bp <- ggplot(treatAfterRpmFold, aes(x=class, y=log2(value), fill=class)) +
      stat_boxplot(geom = "errorbar",
                   width = 0.15) +
      geom_violin(trim=FALSE,alpha=0.8)+
      stat_summary(fun.data=data_summary,geom="crossbar", width=0.2,color="black")+
    labs(x="group", y = "log2(treatAfterRpmFold)")
    treatAfterRpmFold_bp<-treatAfterRpmFold_bp + scale_fill_manual(values=c(brewer.pal(9,"Set1")[2],brewer.pal(9,"Set1")[1])) + theme_classic()+
    theme(text=element_text(size=8), legend.position = "none") + font("xy.text", size = 8)+
    stat_compare_means(method = "t.test",comparisons = my_comparisons,label.y = max(log2(treatAfterRpmFold$value))+abs(quantile(log2(treatAfterRpmFold$value))[1]))+
    scale_y_continuous(breaks = scales::pretty_breaks(n = 5))+
    scale_y_continuous(expand = expansion(mult = c(0, 0.1)))+
    scale_fill_manual(values=c("#129a92","#1e95d4"))

    #preRpmFoldRatio
    preRpmFoldRatio<-SVM_data_melt%>%filter(str_detect(.$variable,"^preRpmFoldRatio$"))
    preRpmFoldRatio_bp <- ggplot(preRpmFoldRatio, aes(x=class, y=log2(value), fill=class)) +
      stat_boxplot(geom = "errorbar",
                   width = 0.15) +
      geom_violin(trim=FALSE,alpha=0.8)+
      stat_summary(fun.data=data_summary,geom="crossbar", width=0.2,color="black")+
    labs(x="group", y = "log2(preRpmFoldRatio)")
    preRpmFoldRatio_bp<-preRpmFoldRatio_bp + scale_fill_manual(values=c(brewer.pal(9,"Set1")[2],brewer.pal(9,"Set1")[1])) + theme_classic()+
    # theme(text=element_text(size=8), axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1),legend.position = "none") + font("xy.text", size = 8)+
    theme(text=element_text(size=8), legend.position = "none") + font("xy.text", size = 8)+
    stat_compare_means(method = "t.test",comparisons = my_comparisons,label.y = max(log2(preRpmFoldRatio$value))+abs(quantile(log2(preRpmFoldRatio$value))[1]))+
    scale_y_continuous(breaks = scales::pretty_breaks(n = 5))+
    scale_y_continuous(expand = expansion(mult = c(0, 0.1)))+
    scale_fill_manual(values=c("#129a92","#1e95d4"))

    #afterRpmFoldRatio
    afterRpmFoldRatio<-SVM_data_melt%>%filter(str_detect(.$variable,"^afterRpmFoldRatio$"))
    afterRpmFoldRatio_bp <- ggplot(afterRpmFoldRatio, aes(x=class, y=log2(value), fill=class)) +
      stat_boxplot(geom = "errorbar",
                   width = 0.15) +
      geom_violin(trim=FALSE,alpha=0.8)+
      stat_summary(fun.data=data_summary,geom="crossbar", width=0.2,color="black")+
      labs(x="group", y = "log2(afterRpmFoldRatio)")
    afterRpmFoldRatio_bp<-afterRpmFoldRatio_bp + scale_fill_manual(values=c(brewer.pal(9,"Set1")[2],brewer.pal(9,"Set1")[1])) + theme_classic()+
    theme(text=element_text(size=8),legend.position = "none") + font("xy.text", size = 8)+
    stat_compare_means(method = "t.test",comparisons = my_comparisons,label.y = max(log2(afterRpmFoldRatio$value))+abs(quantile(log2(afterRpmFoldRatio$value))[1]))+
    scale_y_continuous(breaks = scales::pretty_breaks(n = 5))+
    scale_y_continuous(expand = expansion(mult = c(0, 0.1)))+
    scale_fill_manual(values=c("#129a92","#1e95d4"))

    #treatStopRatio
    treatStopRatio<-SVM_data_melt%>%filter(str_detect(.$variable,"^treatStopRatio$"))
    treatStopRatio<-treatStopRatio[treatStopRatio$value!=0,]
    treatStopRatio_bp <- ggplot(treatStopRatio, aes(x=class, y=log2(value), fill=class)) +
      stat_boxplot(geom = "errorbar",
                   width = 0.15) +
      geom_violin(trim=FALSE,alpha=0.8)+
      stat_summary(fun.data=data_summary,geom="crossbar", width=0.2,color="black")+
      labs(x="group", y = "log2(treatStopRatio)")
    treatStopRatio_bp<-treatStopRatio_bp + scale_fill_manual(values=c(brewer.pal(9,"Set1")[2],brewer.pal(9,"Set1")[1])) + theme_classic()+
    theme(text=element_text(size=8),legend.position = "none") + font("xy.text", size = 8)+
    stat_compare_means(method = "t.test",comparisons = my_comparisons,label.y = max(log2(treatStopRatio$value))+abs(quantile(log2(treatStopRatio$value))[1]))+
    scale_y_continuous(breaks = scales::pretty_breaks(n = 5))+
    scale_y_continuous(expand = expansion(mult = c(0, 0.1)))+
    scale_fill_manual(values=c("#129a92","#1e95d4"))

    stopRatioFC<-SVM_data_melt%>%filter(str_detect(.$variable,"^stopRatioFC$"))
    stopRatioFC_bp <- ggplot(stopRatioFC, aes(x=class, y=log2(value), fill=class)) +
      stat_boxplot(geom = "errorbar",
                   width = 0.15) +
      geom_violin(trim=FALSE,alpha=0.8)+
      stat_summary(fun.data=data_summary,geom="crossbar", width=0.2,color="black")+
      labs(x="group", y = "log2(stopRatioFC)")
    stopRatioFC_bp<-stopRatioFC_bp + scale_fill_manual(values=c(brewer.pal(9,"Set1")[2],brewer.pal(9,"Set1")[1])) + theme_classic()+
    theme(text=element_text(size=8),legend.position = "none") + font("xy.text", size = 8)+
    stat_compare_means(method = "t.test",comparisons = my_comparisons,label.y = max(log2(stopRatioFC$value))+abs(quantile(log2(stopRatioFC$value))[1]))+
    scale_y_continuous(breaks = scales::pretty_breaks(n = 5))+
    scale_y_continuous(expand = expansion(mult = c(0, 0.1)))+
    scale_fill_manual(values=c("#129a92","#1e95d4"))


    pdf(paste(outFile_prefix,"_six_variables_rRNA_violinplot.pdf",sep=""),width=5,height=6)
    plot_grid(
      treatPreRpmFold_bp,
      preRpmFoldRatio_bp,
      treatAfterRpmFold_bp,
      afterRpmFoldRatio_bp,
      treatStopRatio_bp,
      stopRatioFC_bp,
      align = "hv",
      labels = c('A','B','C','D','E','F'),ncol=2,nrow=3)
    invisible(dev.off())


    SVM_data_sel<-SVM_data %>% select(treatPreRpmFold,preRpmFoldRatio,treatAfterRpmFold,afterRpmFoldRatio,treatStopRatio,stopRatioFC,class)
    rownames(SVM_data_sel)<-SVM_data$name
    cpm.of.diff.gene <- t(SVM_data_sel %>% select(-class))
    pca <- prcomp(cpm.of.diff.gene, scale = TRUE)
    fviz_eig(pca,addlabels = T)
    group.list <- c("treatPreRpmFold","preRpmFoldRatio","treatAfterRpmFold","afterRpmFoldRatio","treatStopRatio","stopRatioFC")
    group.list <- as.factor(group.list)
    fviz_pca_ind(pca, geom.ind = c("point"), col.ind = group.list,
                 palette = c(brewer.pal(9,"Set1")[1],brewer.pal(9,"Set1")[2],brewer.pal(9,"Set1")[3],brewer.pal(9,"Set1")[4],brewer.pal(9,"Set1")[5],brewer.pal(9,"Set1")[7]),
                 legend.title = "Groups")
    ggsave(paste(outFile_prefix,"_pca.pdf",sep=""))

    SVM_data_sel$class<-str_replace_all(SVM_data_sel$class, "0", "non-psi")
    SVM_data_sel$class<-str_replace_all(SVM_data_sel$class, "1", "psi")
    #Encoding the target feature as factor
    SVM_data_sel$class<-as.factor(SVM_data_sel$class)
    cat("total classification: ")
    table(SVM_data_sel$class)

    # Splitting the dataset into the Training set and Test set
    set.seed(123)
    split = sample.split(SVM_data_sel$class, SplitRatio = 0.7)
    training_set = subset(SVM_data_sel, split == TRUE)
    test_set = subset(SVM_data_sel, split == FALSE)
    training_set_mean<-apply(training_set %>% select(-class),2,mean)# training_set_sd<-apply(training_set %>% select(-class),2,sd)
    training_set_sd<-c(2,2,4,2,0.05,10)
    training_set[-7] = scale(training_set[-7],training_set_mean,training_set_sd)# test_set[-7] = scale(test_set[-7],center=training_set_mean,scale=training_set_sd)
    training_set_origin = subset(SVM_data, split == TRUE)
    test_set_origin = subset(SVM_data, split == FALSE)
    cat("\n","training set classification using sample.split: ")
    table(training_set$class)
    summary(training_set)#mean equals to mymodel[["x.scale"]][["scaled:center"]]; sd equals to mymodel[["x.scale"]][["scaled:scale"]]
    cat("\n","test set classification using sample.split: ")
    table(test_set$class)
    summary(test_set)

    #evaluate contributation of each variables
    cat("\n\n=====================Evaluate contributation for each variables=====================\n")
    fit2 <- svm(class ~ ., data = SVM_data_sel)
    w <- t(fit2$coefs) %*% fit2$SV                 # weight vectors
    w <- apply(w, 2, function(v){sqrt(sum(v^2))})  # weight
    w <- sort(w, decreasing = T)
    print(w)


    cat("\n\n=====================Get best model using tune (for modeling)=====================\n")

    set.seed(100)
    tmodel = tune(svm,
                  class~.,
                  data=training_set,
                  type="C-classification",
                  kernel="radial",
                  ranges=list(cost=10^(-1:2), gamma=c(.5,1,2)),probability=TRUE,scale=FALSE
    )

    pdf(paste(outFile_prefix, '_tmodel_plot.pdf', sep=""))
    plot(tmodel)
    invisible(dev.off())

    summary(tmodel)
    mymodel <- tmodel$best.model
    mymodel$scaled<-as.logical(rep("TRUE",6))
    mymodel[["x.scale"]][["scaled:center"]]<-training_set_mean
    attr(mymodel[["x.scale"]][["scaled:center"]],"names")<-colnames(SVM_data_sel)[1:6]
    mymodel[["x.scale"]][["scaled:scale"]]<-training_set_sd
    attr(mymodel[["x.scale"]][["scaled:scale"]],"names")<-colnames(SVM_data_sel)[1:6]

    summary(mymodel)
    str(mymodel)
    save(mymodel, file = paste(outFile_prefix, '_SVM_model.RData', sep=""))#"my-svm.RData"
    saveRDS(mymodel, file = paste(outFile_prefix, '_SVM_model.rds', sep=""))#"my-svm.rds"
    write.svm(mymodel,scale.file = paste(outFile_prefix, '_SVM_model.scale', sep=""), yscale.file = paste(outFile_prefix, '_SVM_model.yscale', sep=""))
    pred <- predict(mymodel,test_set,probability=TRUE, decision.values=TRUE)

    #Visualize the prediction effect of the model
    plot_confusion_matrix <- ggplot() +
    geom_confmat(aes(x = test_set$class, y = pred),
                            normalize = TRUE, text.perc = TRUE) +
      labs(x = "Reference", y = "Prediction") +
      scale_fill_gradient2(low = "darkblue", high = "#ec2f2f") +
      theme_bw() +
      theme(plot.margin = unit(c(6, 5, 6, 5), "cm"))

    pdf(paste(outFile_prefix, '_svm_best_test_confusion_matrix.pdf', sep = ""))
    plot_confusion_matrix
    invisible(dev.off())

    #get attr
    pred.decision.values<-as.vector(attr(pred, "decision.values"))
    pred.probabilities<-attr(pred, "probabilities")
    pred.probabilities<-as.data.frame(pred.probabilities)

    pdf(paste(outFile_prefix, '_SVM_roc_test_data_plot.pdf', sep=""))
    svm_roc_test<-roc(main="Test Data ROC",test_set$class,pred.probabilities$psi,smooth=FALSE, print.auc=TRUE, col="#e41a1c", plot = TRUE, print.thres="best", print.thres.best.method="youden", levels = c("non-psi","psi"), direction='<',auc=T, ci=T)
    invisible(dev.off())

    #add attr to SVM_test_data
    SVM_test_data<-test_set_origin
    SVM_test_data$pred.decision.values<-pred.decision.values
    coords(svm_roc_test, "best", ret = "all", transpose = TRUE)
    SVM_test_data$svm_test_prob_class<-ifelse(pred.probabilities$psi>coords(svm_roc_test, "best", ret = "all", transpose = TRUE)[1],"psi","non-psi")
    SVM_test_data<-cbind(SVM_test_data,pred.probabilities,pred)
    write.xlsx(SVM_test_data,paste(outFile_prefix, '_SVM_test_data.xlsx', sep=""),overwrite = TRUE)

    #output model info
    cat("\n\n=====================tune best model confusion matrix (for modeling)=====================\n")
    tab <- table(Predicted = pred,Actual = test_set$class)
    tab
    cat("tune best model error rate (for modeling): ",1-sum(diag(tab))/sum(tab),"\n")
    cat("tune best model correct rate (for modeling): ",sum(diag(tab))/sum(tab),"\n")

    #calculate evaluation indicators
    cat("\n\n=====================Calculate evaluation indicators=====================\n")
    confusion_matrix<-as.data.frame(tab)
    SVM_TP <- confusion_matrix[confusion_matrix$Predicted=="psi"&confusion_matrix$Actual=="psi",]$Freq#True Positives (TP)
    SVM_FP <- confusion_matrix[confusion_matrix$Predicted=="psi"&confusion_matrix$Actual=="non-psi",]$Freq#False Positives (FP)
    SVM_TN <- confusion_matrix[confusion_matrix$Predicted=="non-psi"&confusion_matrix$Actual=="non-psi",]$Freq#True Negatives (TN)
    SVM_FN <- confusion_matrix[confusion_matrix$Predicted=="non-psi"&confusion_matrix$Actual=="psi",]$Freq#False Negatives (FN)
    SVM_TPR <- SVM_TP / (SVM_TP + SVM_FN)#sensitivity (true positive rate, TPR)
    SVM_TNR <- SVM_TN / (SVM_TN + SVM_FP)#specifity (selectivity or true negative rate, TNR)
    SVM_FPR <- 1-SVM_TNR#False Positive Rate (FPR) (1 - specificit = FP/​N = FP/(TN + FP), FPR)
    SVM_FNR <- 1-SVM_TPR#False Negative Rate, FNR)
    SVM_Prec <- SVM_TP / (SVM_TP + SVM_FP)#Precision
    SVM_Recall <- SVM_TP / (SVM_TP + SVM_FN)#Recall
    SVM_ACC <- (SVM_TP + SVM_TN) / (SVM_TP + SVM_TN + SVM_FP + SVM_FN)#accuracy
    SVM_F1_score <- (2*SVM_Recall*SVM_Prec) / (SVM_Recall + SVM_Prec)#F1_score
    eval<-cbind(SVM_TP,SVM_FP,SVM_TN,SVM_FN,SVM_TPR,SVM_TNR,SVM_FPR,SVM_FNR,SVM_Prec,SVM_Recall,SVM_ACC,SVM_F1_score)
    eval<-round(eval,3)
    eval
    write.table(eval,paste(outFile_prefix, '_SVM_eval.txt', sep=""),sep="\t",row.names=FALSE, col.names=TRUE,quote=F)

    #show svm evaluation as pdf table
    svm_eval_t_df<-as.data.frame(t(as.data.frame(eval)))
    colnames(svm_eval_t_df)<-"value_or_percentage"
    tt3 <- ttheme_minimal(
      core=list(bg_params = list(fill = blues9[1:4], col=NA),
                fg_params=list(fontface=3)),
      colhead=list(fg_params=list(col="navyblue", fontface=4L)),
      rowhead=list(fg_params=list(col="orange", fontface=3L)))

    pdf(paste(outFile_prefix, '_svm_evaluation.pdf', sep=""), width = 7, height = 7) # Open a new pdf file
    grid.arrange(
      tableGrob(svm_eval_t_df, theme=tt3),
      nrow=1)
    invisible(dev.off()) # Close the file

    #filt by svm best model
    cat("\n\n=====================Filt by svm best model=====================\n")
    cat("below is your input data ready to be predicted...\n")
    # get prediction
    to_pred<-read.table(filtfile)
    colnames(to_pred)<-c("chrom",#1
      "chromStart",#2
      "chromEnd",#3
      "name",#4
      "foldChange",#5
      "strand",#6
      "geneName",#7
      "geneStart",#8
      "geneEnd",#9
      "base",#10
      "treatPval",#11
      "ctrlPval",#12
      "minusPval",#13
      "treatStopNum",#14
      "treatStopRPM",#15
      "treatPreStopNum",#16
      "treatAfterStopNum",#17
      "treatReadthroughNum",#18
      "ctrlStopNum",#19
      "ctrlStopRPM",#20
      "ctrlPreStopNum",#21
      "ctrlAfterStopNum",#22
      "ctrlReadthroughNum",#23
      "stopRpmFC",#24
      "treatPreRpmFold",#25
      "ctrlPreRpmFold",#26
      "preRpmFoldRatio",#27
      "treatAfterRpmFold",#28
      "ctrlAfterRpmFold",#29
      "afterRpmFoldRatio",#30
      "treatStopRatio",#31
      "ctrlStopRatio",#32
      "stopRatioFC",#33
      "treatStopMeanNum",#34
      "treatStopMeanFold",#35
      "ctrlStopMeanNum",#36
      "ctrlStopMeanFold",#37
      "treatStopMeanFoldRatio",#38
      "extendSeq")#39
    to_pred$base<-"T"
    to_pred_var<-to_pred %>% select(treatPreRpmFold,preRpmFoldRatio,treatAfterRpmFold,afterRpmFoldRatio,treatStopRatio,stopRatioFC)
    str(to_pred_var)

    #get prediction
    pred <- predict(mymodel,to_pred_var,decision.values=TRUE,probability=TRUE)

    #get attr
    pred.decision.values<-as.vector(attr(pred, "decision.values"))
    pred.probabilities<-attr(pred, "probabilities")
    pred.probabilities<-as.data.frame(pred.probabilities)

    #add attr to SVM_pred_data
    SVM_pred_data<-to_pred
    SVM_pred_data$pred.decision.values<-pred.decision.values
    SVM_pred_data$svm_pred_desc_class<-ifelse(SVM_pred_data$pred.decision.values>0,"psi","non-psi")
    SVM_pred_data<-cbind(SVM_pred_data,pred.probabilities,pred)
    SVM_pred_data$svm_pred_prob_class<-ifelse(SVM_pred_data$psi>=0.5,"psi","non-psi")
    table(SVM_pred_data$svm_pred_prob_class)
    write.xlsx(SVM_pred_data,paste(outFile_prefix, '_svm_total_prediction.xlsx', sep=""),overwrite = TRUE)
    write.table(SVM_pred_data,paste(outFile_prefix, '_svm_total_prediction.txt', sep=""),sep="\t",row.names=FALSE, col.names=TRUE,quote=F)
    write.table(SVM_pred_data,paste(outFile_prefix, '_svm_total_prediction.bed', sep=""),sep="\t",row.names=FALSE, col.names=FALSE,quote=F)

    #read hg38_human_chr21_rRNA_known_pseudoU_SingleSites.bed evidence
    evidence<-read.table(rRNAfile,head=F)#"hg38_human_chr21_rRNA_known_pseudoU_SingleSites.bed"
    colnames(evidence)<-c("chrom","chromStart","chromEnd","rRNA_anno","score","strand")
    evidence$rRNA_uniq_id<-paste(evidence$chrom,evidence$chromStart,evidence$chromEnd,evidence$strand,sep="_")

    #read hg38.psiU.SingleSites.bed
    evidence2<-read.table(rRNAfile2,head=F)#"hg38.psiU.SingleSites.bed"
    colnames(evidence2)<-c("chrom","chromStart","chromEnd","rRNA_anno","score","strand")
    evidence2$rRNA_uniq_id<-paste(evidence2$chrom,evidence2$chromStart,evidence2$chromEnd,evidence2$strand,sep="_")

    #known_data miss/hit hg38_human_chr21_rRNA_known_pseudoU_SingleSites
    SVM_data$rRNA_uniq_id<-paste(SVM_data$chrom,SVM_data$chromStart,SVM_data$chromEnd,SVM_data$strand,sep="_")
    SVM_data_evidence<-SVM_data %>% left_join(evidence,by=c("rRNA_uniq_id"="rRNA_uniq_id")) %>% filter(!is.na(rRNA_anno))
    write.csv(SVM_data_evidence,paste(outFile_prefix,"_hg38_human_chr21_rRNA_known_pseudoU_SingleSites_rtsSeeker_hit.csv",sep=""))
    SVM_data_no_evidence<-evidence %>% left_join(SVM_data,by=c("rRNA_uniq_id"="rRNA_uniq_id")) %>% filter(is.na(class))
    write.csv(SVM_data_no_evidence,paste(outFile_prefix,"_hg38_human_chr21_rRNA_known_pseudoU_SingleSites_rtsSeeker_miss.csv",sep=""))
    recovery<-paste(round(length(unique(SVM_data_evidence$rRNA_uniq_id))/length(unique(evidence$rRNA_uniq_id))*100,2),"%",sep="")
    cat("rtsSeeker recover (hg38_human_chr21_rRNA_known_pseudoU_SingleSites.bed)",recovery,"rRNA psi sites in all known chrom21\n")

    #known_data miss/hit hg38.psiU.SingleSites.bed
    SVM_data$rRNA_uniq_id<-paste(SVM_data$chrom,SVM_data$chromStart,SVM_data$chromEnd,SVM_data$strand,sep="_")
    SVM_data_evidence2<-SVM_data %>% left_join(evidence2,by=c("rRNA_uniq_id"="rRNA_uniq_id")) %>% filter(!is.na(rRNA_anno))
    write.csv(SVM_data_evidence2,paste(outFile_prefix,"_hg38.psiU.SingleSites.bed_rtsSeeker_hit.csv",sep=""))
    SVM_data_no_evidence2<-evidence2 %>% left_join(SVM_data,by=c("rRNA_uniq_id"="rRNA_uniq_id")) %>% filter(is.na(class))
    write.csv(SVM_data_no_evidence2,paste(outFile_prefix,"_hg38.psiU.SingleSites.bed_rtsSeeker_miss.csv",sep=""))
    recovery<-paste(round(length(unique(SVM_data_evidence2$rRNA_uniq_id))/length(unique(evidence2$rRNA_uniq_id))*100,2),"%",sep="")
    cat("rtsSeeker recover (hg38.psiU.SingleSites.bed.bed)",recovery,"rRNA psi sites in all known chrom21\n")

    #final_pred miss/hit
    final_pred<-SVM_pred_data[SVM_pred_data$svm_pred_prob_class=="psi",]
    final_pred<-final_pred[final_pred$foldChange>2,]
    final_pred<-final_pred %>% arrange(desc(pred.decision.values))
    write.table(final_pred,paste(outFile_prefix, '_svm_psi_prediction.txt', sep=""),sep="\t",row.names=FALSE, col.names=TRUE,quote=F)
    write.table(final_pred,paste(outFile_prefix, '_svm_psi_prediction.bed', sep=""),sep="\t",row.names=FALSE, col.names=FALSE,quote=F)

    final_pred$rRNA_uniq_id<-paste(final_pred$chrom,final_pred$chromStart,final_pred$chromEnd,final_pred$strand,sep="_")
    final_pred_evidence<-final_pred %>% left_join(evidence,by=c("rRNA_uniq_id"="rRNA_uniq_id")) %>% filter(!is.na(rRNA_anno))
    write.csv(final_pred_evidence,paste(outFile_prefix,"_hg38_human_chr21_rRNA_known_pseudoU_SingleSites_rtsSeeker_svm_hit.csv",sep=""))
    final_pred_no_evidence<-evidence %>% left_join(final_pred,by=c("rRNA_uniq_id"="rRNA_uniq_id")) %>% filter(is.na(base))
    write.csv(final_pred_no_evidence,paste(outFile_prefix,"_hg38_human_chr21_rRNA_known_pseudoU_SingleSites_rtsSeeker_svm_miss.csv",sep=""))
    recovery<-paste(round(length(unique(final_pred_evidence$rRNA_uniq_id))/length(unique(evidence$rRNA_uniq_id))*100,2),"%",sep="")
    cat("rtsSeeker+SVM recover (hg38_human_chr21_rRNA_known_pseudoU_SingleSites.bed)",recovery,"rRNA psi sites in all known chrom21\n")

``svm_predict_totalRNA.sh``

.. code:: bash


    #! /bin/bash

    if [[ $# -ne 2 ]]
    then
        echo 'Usage: '$0 ' input SVM_model'
        exit 1
    fi

    input=$1
    SVM_model=$2

    awk 'FS=OFS="\t" {if($1~/^chr[0-9|a-z|A-Z]*$/ && $10=="T" && $14>10){print $0}}' $input > ${input%.bed}_svm_filt.bed
    nohup Rscript $(dirname "$0")/svm_predict_totalRNA.r -f ${input%.bed}_svm_filt.bed -k $SVM_model -r $(dirname "$0")/hg38_human_chr21_rRNA_known_pseudoU_SingleSites.bed -s $(dirname "$0")/hg38.psiU.SingleSites.bed  -o ${input%.bed} > ${input%.bed}.log
    wait
    cat ${input%.bed}.log

    echo "total RNA: SVM predict program is done!"
    echo -e "total RNA: SVM predict result in $(dirname ${input})"


``svm_predict_totalRNA.r``

.. code:: R

    #!/usr/bin/env Rscript

    options(warn=-1)
    suppressMessages(library("optparse"))
    suppressMessages(library("dplyr"))
    suppressMessages(library("stringr"))
    suppressMessages(library("reshape2"))
    suppressMessages(library("openxlsx"))
    suppressMessages(library("e1071"))
    suppressMessages(library("cowplot"))
    options(warn=-1)

    option_list = list(
      make_option(c("-f", "--filtfile"), type="character", default=NULL,
                  help="filted file of single sites [file]", metavar="character"),
      make_option(c("-k", "--svmmodelfile"), type="character", default=NULL,
                  help="svmmodelfile [file]", metavar="character"),
      make_option(c("-r", "--rRNAfile"), type="character", default=NULL,
                  help="rRNA of single sites [file]", metavar="character"),
      make_option(c("-s", "--rRNAfile2"), type="character", default=NULL,
                  help="hg38.psiU.SingleSites.bed [file]", metavar="character"),
      make_option(c("-o", "--outfile_prefix"), type="character", default=NULL,
                  help="output file name [default= %default]", metavar="character")
    );
    opt_parser = OptionParser(option_list=option_list);
    opt = parse_args(opt_parser);

    if ( is.null(opt$filtfile) || is.null(opt$svmmodelfile)|| is.null(opt$rRNAfile) || is.null(opt$rRNAfile2) || is.null(opt$outfile_prefix) ){
      print_help(opt_parser);
      stop("Please provide -f svmmodelfile, -k filtfile, -r rRNAfile, -s rRNAfile2, and -o outfile_prefix option", call.=FALSE);
    }

    filtfile = opt$filtfile
    svmmodelfile = opt$svmmodelfile
    rRNAfile = opt$rRNAfile
    rRNAfile2 = opt$rRNAfile2
    outFile_prefix = opt$outfile_prefix

    print(filtfile)
    print(svmmodelfile)
    print(rRNAfile)
    print(rRNAfile2)
    print(outFile_prefix)



    # mymodel<-readRDS(file = svmmodelfile)#"my-svm.RData"
    load(file = svmmodelfile)#"my-svm.RData"

    #filt by svm best model
    cat("\n\n=====================Filt by svm best model=====================\n")
    cat("below is your input data ready to be predicted...\n")
    # get prediction
    to_pred<-read.table(filtfile)
    colnames(to_pred)<-c("chrom",#1
      "chromStart",#2
      "chromEnd",#3
      "name",#4
      "foldChange",#5
      "strand",#6
      "geneName",#7
      "geneStart",#8
      "geneEnd",#9
      "base",#10
      "treatPval",#11
      "ctrlPval",#12
      "minusPval",#13
      "treatStopNum",#14
      "treatStopRPM",#15
      "treatPreStopNum",#16
      "treatAfterStopNum",#17
      "treatReadthroughNum",#18
      "ctrlStopNum",#19
      "ctrlStopRPM",#20
      "ctrlPreStopNum",#21
      "ctrlAfterStopNum",#22
      "ctrlReadthroughNum",#23
      "stopRpmFC",#24
      "treatPreRpmFold",#25
      "ctrlPreRpmFold",#26
      "preRpmFoldRatio",#27
      "treatAfterRpmFold",#28
      "ctrlAfterRpmFold",#29
      "afterRpmFoldRatio",#30
      "treatStopRatio",#31
      "ctrlStopRatio",#32
      "stopRatioFC",#33
      "treatStopMeanNum",#34
      "treatStopMeanFold",#35
      "ctrlStopMeanNum",#36
      "ctrlStopMeanFold",#37
      "treatStopMeanFoldRatio",#38
      "extendSeq")#39
    to_pred$base<-"T"
    to_pred_var<-to_pred %>% select(treatPreRpmFold, preRpmFoldRatio, treatAfterRpmFold, afterRpmFoldRatio, treatStopRatio, stopRatioFC)
    str(to_pred_var)

    pred <- predict(mymodel,to_pred_var,decision.values=TRUE,probability=TRUE)

    #get attr
    pred.decision.values<-as.vector(attr(pred, "decision.values"))
    pred.probabilities<-attr(pred, "probabilities")
    pred.probabilities<-as.data.frame(pred.probabilities)

    #add attr to SVM_pred_data
    SVM_pred_data<-to_pred
    SVM_pred_data$pred.decision.values<-pred.decision.values
    SVM_pred_data$svm_pred_desc_class<-ifelse(SVM_pred_data$pred.decision.values>0,"psi","non-psi")
    SVM_pred_data<-cbind(SVM_pred_data,pred.probabilities,pred)
    SVM_pred_data$svm_pred_prob_class<-ifelse(SVM_pred_data$psi>=0.5,"psi","non-psi")
    print("Psi predicted by SVM model:")
    table(SVM_pred_data$svm_pred_prob_class)
    print("Psi predicted by SVM model and fold change (stopRatioFC) >2:")
    SVM_pred_data<-SVM_pred_data[SVM_pred_data$foldChange>2,]
    write.xlsx(SVM_pred_data,paste(outFile_prefix, '_svm_total_prediction.xlsx', sep=""),overwrite = TRUE)
    write.table(SVM_pred_data,paste(outFile_prefix, '_svm_total_prediction.txt', sep=""),sep="\t",row.names=FALSE, col.names=TRUE,quote=F)
    write.table(SVM_pred_data,paste(outFile_prefix, '_svm_total_prediction.bed', sep=""),sep="\t",row.names=FALSE, col.names=FALSE,quote=F)

    #read hg38_human_chr21_rRNA_known_pseudoU_SingleSites.bed evidence
    evidence<-read.table(rRNAfile,head=F)#"hg38_human_chr21_rRNA_known_pseudoU_SingleSites.bed"
    colnames(evidence)<-c("chrom","chromStart","chromEnd","rRNA_anno","score","strand")
    evidence$rRNA_uniq_id<-paste(evidence$chrom,evidence$chromStart,evidence$chromEnd,evidence$strand,sep="_")

    #read hg38.psiU.SingleSites.bed
    evidence2<-read.table(rRNAfile2,head=F)#"hg38.psiU.SingleSites.bed"
    colnames(evidence2)<-c("chrom","chromStart","chromEnd","rRNA_anno","score","strand")
    evidence2$rRNA_uniq_id<-paste(evidence2$chrom,evidence2$chromStart,evidence2$chromEnd,evidence2$strand,sep="_")


    final_pred<-SVM_pred_data[SVM_pred_data$svm_pred_prob_class=="psi",]
    final_pred<-final_pred %>% arrange(desc(pred.decision.values))
    write.table(final_pred,paste(outFile_prefix, '_svm_psi_prediction.txt', sep=""),sep="\t",row.names=FALSE, col.names=TRUE,quote=F)
    write.table(final_pred,paste(outFile_prefix, '_svm_psi_prediction.bed', sep=""),sep="\t",row.names=FALSE, col.names=FALSE,quote=F)

    final_pred$rRNA_uniq_id<-paste(final_pred$chrom,final_pred$chromStart,final_pred$chromEnd,final_pred$strand,sep="_")
    final_pred_evidence<-final_pred %>% left_join(evidence,by=c("rRNA_uniq_id"="rRNA_uniq_id")) %>% filter(!is.na(rRNA_anno))
    write.csv(final_pred_evidence,paste(outFile_prefix,"_hg38_human_chr21_rRNA_known_pseudoU_SingleSites_rtsSeeker_svm_hit.csv",sep=""))
    final_pred_no_evidence<-evidence %>% left_join(final_pred,by=c("rRNA_uniq_id"="rRNA_uniq_id")) %>% filter(is.na(base))
    write.csv(final_pred_no_evidence,paste(outFile_prefix,"_hg38_human_chr21_rRNA_known_pseudoU_SingleSites_rtsSeeker_svm_miss.csv",sep=""))
    recovery<-paste(round(length(unique(final_pred_evidence$rRNA_uniq_id))/length(unique(evidence$rRNA_uniq_id))*100,2),"%",sep="")
    cat("rtsSeeker+SVM recover (hg38_human_chr21_rRNA_known_pseudoU_SingleSites.bed)",recovery,"rRNA psi sites in all known chrom21\n")

.. note:: All user input will be recorded in a plain text file with suffixes ``_svm_totalRNA_build.txt`` and ``_svm_totalRNA_predict.txt`` in psiFinder/config and help users to easily reload the previous config (by simply clicking ``CONFIG`` button).

(2) ANN: Artificial Neural Network
*****************************************
If ``ANN`` QT wideget is clicked and popouped:

-  press ``Build`` buttor will run ``ann.sh`` and ``ann_build_totalRNA.r`` to generate file with ``_ANN_model.RData`` suffix, which is the artifical neuralnet model (resilient backpropagation with weight backtracking, two layers, each 12 neuron).
-  press ``PREDICT`` buttor will run ``ann_predict_totalRNA.sh`` and ``ann_predict_totalRNA.r`` to generate file with ``_ann_psi_prediction.bed`` suffix, which is the candidate Ψ-sites predicted by ANN algorithm.

.. image:: /images/ANN.png

``ann.sh``

.. code:: bash

    #! /bin/bash

    if [[ $# -ne 4 ]]
    then
        echo 'Usage: '$0 ' input.bam treatment.bam genome.fa out.bed'
        exit 1
    fi

    input=$1
    treatment=$2
    genome=$3
    out=$4
    outFileName=${out%.bed}

    echo "rtsSeeker ANN mode: rtsSeeker start..."
    if [ ! -f "${genome}.fai"  ]
    then
    echo "generate faifile"
    samtools faidx $genome
    else
    echo "faifile exist"
    fi

    echo "$(dirname "$0")/rtsSeeker --fa $genome --fai ${genome}.fai --treat $treatment --input $input -p 1.5 -t 5 -r 0.05 -M 1 --gene $(dirname "$0")/hg38.gencode.v30.tRNA.refseqNcRNA.geneAnno.bed12 -f 1 -m 0 -s -w 20 -o ${outFileName}_ann_all.bed 2>${outFileName}_ann_all_rtsSeeker.log" > ${outFileName}_ann_all_rtsSeeker.cmd # no -c remove duplicate
    $(dirname "$0")/rtsSeeker --fa $genome --fai ${genome}.fai --treat $treatment --input $input -p 1.5 -t 5 -r 0.05 -M 1 --gene $(dirname "$0")/hg38.gencode.v30.tRNA.refseqNcRNA.geneAnno.bed12 -f 1 -m 0 -s -w 20 -o ${outFileName}_ann_all.bed 2>${outFileName}_ann_all_rtsSeeker.log # no -c remove duplicate
    cat ${outFileName}_ann_all.bed > $out

    awk 'FS=OFS="\t" {if($1~/^chr[0-9|a-z|A-Z]*$/ && $10=="T" && $14>10){print $0}}' ${outFileName}_ann_all.bed > ${outFileName}_ann_filt_totalRNA.bed
    bedtools intersect -a ${outFileName}_ann_filt_totalRNA.bed -b $(dirname "$0")/hg38_human_chr21_rRNA_known_pseudoU_SingleSites.bed -s > ${outFileName}_knowpse.bed
    bedtools intersect -a ${outFileName}_ann_filt_totalRNA.bed -b $(dirname "$0")/rrna_chr21.bed -s > ${outFileName}_rrna.bed
    bedtools intersect -a ${outFileName}_rrna.bed -b $(dirname "$0")/hg38_human_chr21_rRNA_known_pseudoU_SingleSites.bed -s -v > ${outFileName}_notpse.bed
    cat ${outFileName}_knowpse.bed |awk 'FS=OFS="\t" {print $0,"1"}' > ${outFileName}_knowpse.txt
    cat ${outFileName}_notpse.bed |awk 'FS=OFS="\t" {print $0,"0"}' > ${outFileName}_notpse.txt
    echo "getting roc_plot.txt"
    cat ${outFileName}_knowpse.txt ${outFileName}_notpse.txt > ${outFileName}_ann_plot.txt

    nohup Rscript $(dirname "$0")/ann_build_totalRNA.r -f ${outFileName}_ann_plot.txt -k ${outFileName}_ann_filt_totalRNA.bed -r $(dirname "$0")/hg38_human_chr21_rRNA_known_pseudoU_SingleSites.bed -s $(dirname "$0")/hg38.psiU.SingleSites.bed  -o ${outFileName} > ${outFileName}_ann_evaluation_totalRNA.log 2>&1 &
    wait

    cat ${outFileName}_ann_evaluation_totalRNA.log
    mupdf-x11 ${outFileName}_six_variables_rRNA_violinplot.pdf &> /dev/null
    mupdf-x11 ${outFileName}_ann_roc_test_data_plot.pdf &> /dev/null
    mupdf-x11 ${outFileName}_ann_evaluation.pdf &> /dev/null

    echo "total RNA: ann program is done!"
    echo -e "total RNA: ann result in $(dirname ${outFileName}_ann_psi_prediction.bed)"

``ann_build_totalRNA.r``

.. code:: R

    #!/usr/bin/env Rscript

    options(warn = -1)
    suppressMessages(library("optparse"))
    suppressMessages(library("devtools"))
    suppressMessages(library("caTools"))
    suppressMessages(library("neuralnet"))
    suppressMessages(library("NeuralNetTools"))
    suppressMessages(library("dplyr"))
    suppressMessages(library("stringr"))
    suppressMessages(library("gridExtra"))
    suppressMessages(library("cowplot"))
    suppressMessages(library("pROC"))
    suppressMessages(library("ggplot2"))
    suppressMessages(library("ggpol"))
    suppressMessages(library("ggpubr"))
    suppressMessages(library("RColorBrewer"))
    suppressMessages(library("openxlsx"))
    suppressMessages(library("reshape2"))
    suppressMessages(library("factoextra"))

    options(warn=-1)

    option_list = list(
      make_option(c("-f", "--annfile"), type="character", default=NULL,
                  help="ROC of single sites [file]", metavar="character"),
      make_option(c("-k", "--filtfile"), type="character", default=NULL,
                  help="filted file of single sites [file]", metavar="character"),
      make_option(c("-r", "--rRNAfile"), type="character", default=NULL,
                  help="rRNA of single sites [file]", metavar="character"),
      make_option(c("-s", "--rRNAfile2"), type="character", default=NULL,
                  help="hg38.psiU.SingleSites.bed [file]", metavar="character"),
      make_option(c("-o", "--outfile_prefix"), type="character", default=NULL,
                  help="output file name [default= %default]", metavar="character")
    );
    opt_parser = OptionParser(option_list=option_list);
    opt = parse_args(opt_parser);

    if (is.null(opt$annfile)|| is.null(opt$filtfile) || is.null(opt$rRNAfile) || is.null(opt$rRNAfile2) || is.null(opt$outfile_prefix) ){
      print_help(opt_parser);
      stop("Please provide -f annfile, -k filtfile, -r rRNAfile, -s rRNAfile2, and -o outfile_prefix option", call.=FALSE);
    }

    annfile = opt$annfile
    filtfile = opt$filtfile
    rRNAfile = opt$rRNAfile
    rRNAfile2 = opt$rRNAfile2
    outFile_prefix = opt$outfile_prefix

    print(annfile)
    print(filtfile)
    print(rRNAfile)
    print(rRNAfile2)
    print(outFile_prefix)


    #load data
    cat("\n\n=====================Load data to perform classification model (factor as input)=====================\n")
    ann_data<-read.table(annfile)
    colnames(ann_data)<-c("chrom",#1
      "chromStart",#2
      "chromEnd",#3
      "name",#4
      "foldChange",#5
      "strand",#6
      "geneName",#7
      "geneStart",#8
      "geneEnd",#9
      "base",#10
      "treatPval",#11
      "ctrlPval",#12
      "minusPval",#13
      "treatStopNum",#14
      "treatStopRPM",#15
      "treatPreStopNum",#16
      "treatAfterStopNum",#17
      "treatReadthroughNum",#18
      "ctrlStopNum",#19
      "ctrlStopRPM",#20
      "ctrlPreStopNum",#21
      "ctrlAfterStopNum",#22
      "ctrlReadthroughNum",#23
      "stopRpmFC",#24
      "treatPreRpmFold",#25 *
      "ctrlPreRpmFold",#26
      "preRpmFoldRatio",#27 *
      "treatAfterRpmFold",#28 *
      "ctrlAfterRpmFold",#29
      "afterRpmFoldRatio",#30 *
      "treatStopRatio",#31 *
      "ctrlStopRatio",#32
      "stopRatioFC",#33 *
      "treatStopMeanNum",#34
      "treatStopMeanFold",#35
      "ctrlStopMeanNum",#36
      "ctrlStopMeanFold",#37
      "treatStopMeanFoldRatio",#38
      "extendSeq",#39
      "class")#40
    ann_data$base<-"T"

    #rRNA-psi-non_psi visualization
    ann_data_sel<-ann_data %>% select(treatPreRpmFold,preRpmFoldRatio,treatAfterRpmFold,afterRpmFoldRatio,treatStopRatio,stopRatioFC,class)
    ann_data_melt<-melt(ann_data_sel,id.vars = "class")
    ann_data_melt$class<-str_replace(as.character(ann_data_melt$class),"0","non_psi")
    ann_data_melt$class<-str_replace(as.character(ann_data_melt$class),"1","psi")

    data_summary <- function(x) {
       m <- mean(x)
       ymin <- m-sd(x)
       ymax <- m+sd(x)
       return(c(y=m,ymin=ymin,ymax=ymax))
    }

    my_comparisons <- list( c("non_psi", "psi") )

    #treatPreRpmFold
    treatPreRpmFold<-ann_data_melt%>%filter(str_detect(.$variable,"^treatPreRpmFold$"))
    treatPreRpmFold_bp <- ggplot(treatPreRpmFold, aes(x=class, y=log2(value), fill=class)) +
      stat_boxplot(geom = "errorbar",
                   width = 0.15) +
      geom_violin(trim=FALSE,alpha=0.8)+
      stat_summary(fun.data=data_summary,geom="crossbar", width=0.2,color="black")+
    labs(x="group", y = "log2(treatPreRpmFold)")
    treatPreRpmFold_bp<-treatPreRpmFold_bp + scale_fill_manual(values=c(brewer.pal(9,"Set1")[2],brewer.pal(9,"Set1")[1])) + theme_classic()+
    theme(text=element_text(size=8),legend.position = "none") + font("xy.text", size = 8)+
    stat_compare_means(method = "t.test",comparisons = my_comparisons,label.y = max(log2(treatPreRpmFold$value))+abs(quantile(log2(treatPreRpmFold$value))[1]))+
    scale_y_continuous(breaks = scales::pretty_breaks(n = 5))+
    scale_y_continuous(expand = expansion(mult = c(0, 0.1)))+
    scale_fill_manual(values=c("#129a92","#1e95d4"))

    #treatAfterRpmFold
    treatAfterRpmFold<-ann_data_melt%>%filter(str_detect(.$variable,"^treatAfterRpmFold$"))
    treatAfterRpmFold_bp <- ggplot(treatAfterRpmFold, aes(x=class, y=log2(value), fill=class)) +
      stat_boxplot(geom = "errorbar",
                   width = 0.15) +
      geom_violin(trim=FALSE,alpha=0.8)+
      stat_summary(fun.data=data_summary,geom="crossbar", width=0.2,color="black")+
    labs(x="group", y = "log2(treatAfterRpmFold)")
    treatAfterRpmFold_bp<-treatAfterRpmFold_bp + scale_fill_manual(values=c(brewer.pal(9,"Set1")[2],brewer.pal(9,"Set1")[1])) + theme_classic()+
    theme(text=element_text(size=8), legend.position = "none") + font("xy.text", size = 8)+
    stat_compare_means(method = "t.test",comparisons = my_comparisons,label.y = max(log2(treatAfterRpmFold$value))+abs(quantile(log2(treatAfterRpmFold$value))[1]))+
    scale_y_continuous(breaks = scales::pretty_breaks(n = 5))+
    scale_y_continuous(expand = expansion(mult = c(0, 0.1)))+
    scale_fill_manual(values=c("#129a92","#1e95d4"))

    #preRpmFoldRatio
    preRpmFoldRatio<-ann_data_melt%>%filter(str_detect(.$variable,"^preRpmFoldRatio$"))
    preRpmFoldRatio_bp <- ggplot(preRpmFoldRatio, aes(x=class, y=log2(value), fill=class)) +
      stat_boxplot(geom = "errorbar",
                   width = 0.15) +
      geom_violin(trim=FALSE,alpha=0.8)+
      stat_summary(fun.data=data_summary,geom="crossbar", width=0.2,color="black")+
    labs(x="group", y = "log2(preRpmFoldRatio)")
    preRpmFoldRatio_bp<-preRpmFoldRatio_bp + scale_fill_manual(values=c(brewer.pal(9,"Set1")[2],brewer.pal(9,"Set1")[1])) + theme_classic()+
    theme(text=element_text(size=8), legend.position = "none") + font("xy.text", size = 8)+
    stat_compare_means(method = "t.test",comparisons = my_comparisons,label.y = max(log2(preRpmFoldRatio$value))+abs(quantile(log2(preRpmFoldRatio$value))[1]))+
    scale_y_continuous(breaks = scales::pretty_breaks(n = 5))+
    scale_y_continuous(expand = expansion(mult = c(0, 0.1)))+
    scale_fill_manual(values=c("#129a92","#1e95d4"))

    #afterRpmFoldRatio
    afterRpmFoldRatio<-ann_data_melt%>%filter(str_detect(.$variable,"^afterRpmFoldRatio$"))
    afterRpmFoldRatio_bp <- ggplot(afterRpmFoldRatio, aes(x=class, y=log2(value), fill=class)) +
      stat_boxplot(geom = "errorbar",
                   width = 0.15) +
      geom_violin(trim=FALSE,alpha=0.8)+
      stat_summary(fun.data=data_summary,geom="crossbar", width=0.2,color="black")+
      labs(x="group", y = "log2(afterRpmFoldRatio)")
    afterRpmFoldRatio_bp<-afterRpmFoldRatio_bp + scale_fill_manual(values=c(brewer.pal(9,"Set1")[2],brewer.pal(9,"Set1")[1])) + theme_classic()+
    theme(text=element_text(size=8),legend.position = "none") + font("xy.text", size = 8)+
    stat_compare_means(method = "t.test",comparisons = my_comparisons,label.y = max(log2(afterRpmFoldRatio$value))+abs(quantile(log2(afterRpmFoldRatio$value))[1]))+
    scale_y_continuous(breaks = scales::pretty_breaks(n = 5))+
    scale_y_continuous(expand = expansion(mult = c(0, 0.1)))+
    scale_fill_manual(values=c("#129a92","#1e95d4"))

    #treatStopRatio
    treatStopRatio<-ann_data_melt%>%filter(str_detect(.$variable,"^treatStopRatio$"))
    treatStopRatio<-treatStopRatio[treatStopRatio$value!=0,]
    treatStopRatio_bp <- ggplot(treatStopRatio, aes(x=class, y=log2(value), fill=class)) +
      stat_boxplot(geom = "errorbar",
                   width = 0.15) +
      geom_violin(trim=FALSE,alpha=0.8)+
      stat_summary(fun.data=data_summary,geom="crossbar", width=0.2,color="black")+
      labs(x="group", y = "log2(treatStopRatio)")
    treatStopRatio_bp<-treatStopRatio_bp + scale_fill_manual(values=c(brewer.pal(9,"Set1")[2],brewer.pal(9,"Set1")[1])) + theme_classic()+
    theme(text=element_text(size=8),legend.position = "none") + font("xy.text", size = 8)+
    stat_compare_means(method = "t.test",comparisons = my_comparisons,label.y = max(log2(treatStopRatio$value))+abs(quantile(log2(treatStopRatio$value))[1]))+
    scale_y_continuous(breaks = scales::pretty_breaks(n = 5))+
    scale_y_continuous(expand = expansion(mult = c(0, 0.1)))+
    scale_fill_manual(values=c("#129a92","#1e95d4"))

    stopRatioFC<-ann_data_melt%>%filter(str_detect(.$variable,"^stopRatioFC$"))
    stopRatioFC_bp <- ggplot(stopRatioFC, aes(x=class, y=log2(value), fill=class)) +
      stat_boxplot(geom = "errorbar",
                   width = 0.15) +
      geom_violin(trim=FALSE,alpha=0.8)+
      stat_summary(fun.data=data_summary,geom="crossbar", width=0.2,color="black")+
      labs(x="group", y = "log2(stopRatioFC)")
    stopRatioFC_bp<-stopRatioFC_bp + scale_fill_manual(values=c(brewer.pal(9,"Set1")[2],brewer.pal(9,"Set1")[1])) + theme_classic()+
    theme(text=element_text(size=8),legend.position = "none") + font("xy.text", size = 8)+
    stat_compare_means(method = "t.test",comparisons = my_comparisons,label.y = max(log2(stopRatioFC$value))+abs(quantile(log2(stopRatioFC$value))[1]))+
    scale_y_continuous(breaks = scales::pretty_breaks(n = 5))+
    scale_y_continuous(expand = expansion(mult = c(0, 0.1)))+
    scale_fill_manual(values=c("#129a92","#1e95d4"))

    pdf(paste(outFile_prefix,"_six_variables_rRNA_violinplot.pdf",sep=""),width=5,height=6)
    plot_grid(
      treatPreRpmFold_bp,
      preRpmFoldRatio_bp,
      treatAfterRpmFold_bp,
      afterRpmFoldRatio_bp,
      treatStopRatio_bp,
      stopRatioFC_bp,
      align = "hv",
      labels = c('A','B','C','D','E','F'),ncol=2,nrow=3)
    invisible(dev.off())

    ann_data_sel<-ann_data %>% select(treatPreRpmFold,preRpmFoldRatio,treatAfterRpmFold,afterRpmFoldRatio,treatStopRatio,stopRatioFC,class)
    rownames(ann_data_sel)<-ann_data$name
    cpm.of.diff.gene <- t(ann_data_sel %>% select(-class))
    pca <- prcomp(cpm.of.diff.gene, scale = TRUE)
    fviz_eig(pca,addlabels = T)
    group.list <- c("treatPreRpmFold","preRpmFoldRatio","treatAfterRpmFold","afterRpmFoldRatio","treatStopRatio","stopRatioFC")
    group.list <- as.factor(group.list)
    fviz_pca_ind(pca, geom.ind = c("point"), col.ind = group.list,
                 palette = c(brewer.pal(9,"Set1")[1],brewer.pal(9,"Set1")[2],brewer.pal(9,"Set1")[3],brewer.pal(9,"Set1")[4],brewer.pal(9,"Set1")[5],brewer.pal(9,"Set1")[7]),
                 legend.title = "Groups")
    ggsave(paste(outFile_prefix,"_pca.pdf",sep=""))

    ann_data_sel$class<-str_replace_all(ann_data_sel$class, "0", "non_psi")
    ann_data_sel$class<-str_replace_all(ann_data_sel$class, "1", "psi")
    #Encoding the target feature as factor
    ann_data_sel$class<-as.factor(ann_data_sel$class)
    cat("total classification: ")
    table(ann_data_sel$class)

    # Splitting the dataset into the Training set and Test set
    set.seed(123)
    ann_data_sel$psi <- ann_data_sel$class == "psi"
    ann_data_sel$non_psi <- ann_data_sel$class == "non_psi"
    split = sample.split(ann_data_sel$class, SplitRatio = 0.7)
    training_set = subset(ann_data_sel, split == TRUE)
    test_set = subset(ann_data_sel, split == FALSE)

    # feature scaling
    training_set_mean <- apply(training_set %>% select(-class, -psi, -non_psi), 2, mean)
    training_set_sd <- apply(training_set %>% select(-class, -psi, -non_psi), 2, sd)
    training_set[, c(-7, -8, -9)] = scale(training_set[, c(-7, -8, -9)], center = training_set_mean, scale = training_set_sd)
    test_set[, c(-7, -8, -9)] = scale(test_set[, c(-7, -8, -9)], center = training_set_mean, scale = training_set_sd)

    training_set_origin = subset(ann_data, split == TRUE)
    test_set_origin = subset(ann_data, split == FALSE)
    cat("\n", "training set classification using sample.split: ")
    table(training_set$class)
    summary(training_set)
    cat("\n", "test set classification using sample.split: ")
    table(test_set$class)
    summary(test_set)

    #Network Aplication
    ann_data_sel.net <- neuralnet(psi + non_psi ~
                          treatPreRpmFold + preRpmFoldRatio + treatAfterRpmFold + afterRpmFoldRatio + treatStopRatio + stopRatioFC,
                          data = training_set, hidden = c(12, 12), rep = 5, err.fct = "ce",
                          linear.output = FALSE, lifesign = "minimal", stepmax = 1000000,
                          threshold = 0.001)

    #Get weights for a neural network in an organized list by extracting values from a neural network object
    neuralweights(ann_data_sel.net)

    #import the function from Github
    source_url('https://gist.githubusercontent.com/fawda123/7471137/raw/466c1474d0a505ff044412703516c34f1a4684a5/nnet_plot_update.r')
    #plot each model
    pdf(paste(outFile_prefix, '_ann_nnet_train_data_plot.pdf', sep = ""),width=16,height=8)
    plot.nnet(ann_data_sel.net)
    invisible(dev.off())

    #Visualize the network structure of the fully connected model
    pdf(paste(outFile_prefix, '_ann_best_train_data_plot.pdf', sep = ""))
    plot(ann_data_sel.net, rep = "best")
    invisible(dev.off())

    #Visualize the prediction effect of the model
    mlppre <- predict(ann_data_sel.net, test_set)
    colnames(mlppre)<-c("psi","non_psi")
    mlpprelab <- apply(mlppre, 1, which.max)
    mlpprelab<-sub("1","psi",as.character(mlpprelab))
    mlpprelab<-sub("2","non_psi",as.character(mlpprelab))
    plot_confusion_matrix <- ggplot() +
    geom_confmat(aes(x = test_set$class, y = mlpprelab),
                            normalize = TRUE, text.perc = TRUE) +
      labs(x = "Reference", y = "Prediction") +
      scale_fill_gradient2(low = "darkblue", high = "#ec2f2f") +
      theme_bw() +
      theme(plot.margin = unit(c(6, 5, 6, 5), "cm"))

    pdf(paste(outFile_prefix, '_ann_best_test_confusion_matrix.pdf', sep = ""))
    plot_confusion_matrix
    invisible(dev.off())

    # Predicting Result
    ann_data_sel.prediction <- neuralnet::compute(ann_data_sel.net, test_set[-7:-9])
    idx <- apply(ann_data_sel.prediction$net.result, 1, which.max)
    net.result <- as.data.frame(ann_data_sel.prediction$net.result)
    colnames(net.result) <- c('psi', 'non_psi')
    pred <- c('psi', 'non_psi')[idx]
    table(pred)

    pdf(paste(outFile_prefix, '_ann_roc_test_data_plot.pdf', sep=""))
    ann_roc_test <- roc(main = "Test Data ROC", test_set$class, net.result$psi, smooth = FALSE, print.auc = TRUE, col = "#e41a1c", plot = TRUE, print.thres = "best", print.thres.best.method = "youden", levels = c("non_psi", "psi"), direction = '<', auc = T, ci = T)
    invisible(dev.off())
    cat("best roc threshold: ")
    coords(ann_roc_test, "best", ret = "all", transpose = TRUE)[1]

    summary(ann_data_sel.net)
    str(ann_data_sel.net)
    save(ann_data_sel.net, ann_roc_test, training_set_mean, training_set_sd, file = paste(outFile_prefix, '_ANN_model.RData', sep = "")) #"my-nn.RData"
    # saveRDS(mymodel, file = paste(outFile_prefix, '_ann_model.rds', sep=""))#"my-nn.rds"

    #add attr to ann_test_data
    ann_test_data <- test_set_origin
    # ann_test_data$ann_test_prob_class <- ifelse(net.result$psi >= coords(ann_roc_test, "best", ret = "all", transpose = TRUE)[1], "psi", "non_psi")
    ann_test_data <- cbind(ann_test_data, net.result, pred)
    write.xlsx(ann_test_data, paste(outFile_prefix, '_ann_test_data.xlsx', sep = ""), overwrite = TRUE)

    #output model info
    cat("\n\n=====================tune best model confusion matrix (for modeling)=====================\n")
    tab <- table(Predicted = pred, Actual = test_set$class)
    tab
    cat("tune best model error rate (for modeling): ", 1 - sum(diag(tab)) / sum(tab), "\n")
    cat("tune best model correct rate (for modeling): ", sum(diag(tab)) / sum(tab), "\n")

    #calculate evaluation indicators
    cat("\n\n=====================Calculate evaluation indicators=====================\n")
    confusion_matrix <- as.data.frame(tab)
    ann_TP <- confusion_matrix[confusion_matrix$Predicted == "psi" & confusion_matrix$Actual == "psi",]$Freq #True Positives (TP)
    ann_FP <- confusion_matrix[confusion_matrix$Predicted == "psi" & confusion_matrix$Actual == "non_psi",]$Freq #False Positives (FP)
    ann_TN <- confusion_matrix[confusion_matrix$Predicted == "non_psi" & confusion_matrix$Actual == "non_psi",]$Freq #True Negatives (TN)
    ann_FN <- confusion_matrix[confusion_matrix$Predicted == "non_psi" & confusion_matrix$Actual == "psi",]$Freq #False Negatives (FN)
    ann_TPR <- ann_TP / (ann_TP + ann_FN) #sensitivity (true positive rate, TPR)
    ann_TNR <- ann_TN / (ann_TN + ann_FP) #specifity (selectivity or true negative rate, TNR)
    ann_FPR <- 1 - ann_TNR #False Positive Rate (FPR) (1 - specificit = FP/​N = FP/(TN + FP), FPR)
    ann_FNR <- 1 - ann_TPR #False Negative Rate, FNR)
    ann_Prec <- ann_TP / (ann_TP + ann_FP) #Precision
    ann_Recall <- ann_TP / (ann_TP + ann_FN) #Recall
    ann_ACC <- (ann_TP + ann_TN) / (ann_TP + ann_TN + ann_FP + ann_FN) #accuracy
    ann_F1_score <- (2 * ann_Recall * ann_Prec) / (ann_Recall + ann_Prec) #F1_score
    eval <- cbind(ann_TP, ann_FP, ann_TN, ann_FN, ann_TPR, ann_TNR, ann_FPR, ann_FNR, ann_Prec, ann_Recall, ann_ACC, ann_F1_score)
    eval <- round(eval, 3)
    eval
    write.table(eval, paste(outFile_prefix, '_ann_eval.txt', sep = ""), sep = "\t", row.names = FALSE, col.names = TRUE, quote = F)

    #show nn evaluation as pdf table
    ann_eval_t_df <- as.data.frame(t(as.data.frame(eval)))
    colnames(ann_eval_t_df) <- "value_or_percentage"
    tt3 <- ttheme_minimal(
      core = list(bg_params = list(fill = blues9[1:4], col = NA),
                fg_params = list(fontface = 3)),
      colhead = list(fg_params = list(col = "navyblue", fontface = 4L)),
      rowhead = list(fg_params = list(col = "orange", fontface = 3L)))

    pdf(paste(outFile_prefix, '_ann_evaluation.pdf', sep = ""), width = 7, height = 7) # Open a new pdf file
    grid.arrange(
      tableGrob(ann_eval_t_df, theme = tt3),
      nrow = 1)
    invisible(dev.off()) # Close the file


    #filt by nn best model
    cat("\n\n=====================Filt by nn best model=====================\n")
    cat("below is your input data ready to be predicted...\n")
    # get prediction
    to_pred <- read.table(filtfile)
    colnames(to_pred)<-c("chrom",#1
      "chromStart",#2
      "chromEnd",#3
      "name",#4
      "foldChange",#5
      "strand",#6
      "geneName",#7
      "geneStart",#8
      "geneEnd",#9
      "base",#10
      "treatPval",#11
      "ctrlPval",#12
      "minusPval",#13
      "treatStopNum",#14
      "treatStopRPM",#15
      "treatPreStopNum",#16
      "treatAfterStopNum",#17
      "treatReadthroughNum",#18
      "ctrlStopNum",#19
      "ctrlStopRPM",#20
      "ctrlPreStopNum",#21
      "ctrlAfterStopNum",#22
      "ctrlReadthroughNum",#23
      "stopRpmFC",#24
      "treatPreRpmFold",#25
      "ctrlPreRpmFold",#26
      "preRpmFoldRatio",#27
      "treatAfterRpmFold",#28
      "ctrlAfterRpmFold",#29
      "afterRpmFoldRatio",#30
      "treatStopRatio",#31
      "ctrlStopRatio",#32
      "stopRatioFC",#33
      "treatStopMeanNum",#34
      "treatStopMeanFold",#35
      "ctrlStopMeanNum",#36
      "ctrlStopMeanFold",#37
      "treatStopMeanFoldRatio",#38
      "extendSeq")#39
    # to_pred$ann_class<-"1"
    to_pred$base<-"T"

    to_pred_var <- to_pred %>% select(treatPreRpmFold,preRpmFoldRatio,treatAfterRpmFold,afterRpmFoldRatio,treatStopRatio,stopRatioFC)
    to_pred_var = scale(to_pred_var, center = training_set_mean, scale = training_set_sd)
    str(to_pred_var)

    # Predicting Result
    ann_data_sel.prediction <- neuralnet::compute(ann_data_sel.net, to_pred_var)
    idx <- apply(ann_data_sel.prediction$net.result, 1, which.max)
    net.result <- as.data.frame(ann_data_sel.prediction$net.result)
    colnames(net.result) <- c('psi', 'non_psi')
    pred <- c('psi', 'non_psi')[idx]
    table(pred)


    #add attr to ann_pred_data
    ann_pred_data <- to_pred
    ann_pred_data <- cbind(ann_pred_data, net.result, pred)
    write.xlsx(ann_pred_data, paste(outFile_prefix, '_ann_total_prediction.xlsx', sep = ""), overwrite = TRUE)

    #read hg38_human_chr21_rRNA_known_pseudoU_SingleSites.bed evidence
    evidence<-read.table(rRNAfile,head=F)#"hg38_human_chr21_rRNA_known_pseudoU_SingleSites.bed"
    colnames(evidence)<-c("chrom","chromStart","chromEnd","rRNA_anno","score","strand")
    evidence$rRNA_uniq_id<-paste(evidence$chrom,evidence$chromStart,evidence$chromEnd,evidence$strand,sep="_")

    #read hg38.psiU.SingleSites.bed
    evidence2<-read.table(rRNAfile2,head=F)#"hg38.psiU.SingleSites.bed"
    colnames(evidence2)<-c("chrom","chromStart","chromEnd","rRNA_anno","score","strand")
    evidence2$rRNA_uniq_id<-paste(evidence2$chrom,evidence2$chromStart,evidence2$chromEnd,evidence2$strand,sep="_")

    #known_data miss/hit hg38_human_chr21_rRNA_known_pseudoU_SingleSites
    ann_data$rRNA_uniq_id<-paste(ann_data$chrom,ann_data$chromStart,ann_data$chromEnd,ann_data$strand,sep="_")
    ann_data_evidence<-ann_data %>% left_join(evidence,by=c("rRNA_uniq_id"="rRNA_uniq_id")) %>% filter(!is.na(rRNA_anno))
    write.csv(ann_data_evidence,paste(outFile_prefix,"_hg38_human_chr21_rRNA_known_pseudoU_SingleSites_psiFinder_hit.csv",sep=""))
    ann_data_no_evidence<-evidence %>% left_join(ann_data,by=c("rRNA_uniq_id"="rRNA_uniq_id")) %>% filter(is.na(class))
    write.csv(ann_data_no_evidence,paste(outFile_prefix,"_hg38_human_chr21_rRNA_known_pseudoU_SingleSites_psiFinder_miss.csv",sep=""))
    recovery<-paste(round(length(unique(ann_data_evidence$rRNA_uniq_id))/length(unique(evidence$rRNA_uniq_id))*100,2),"%",sep="")
    cat("psiFinder recover (hg38_human_chr21_rRNA_known_pseudoU_SingleSites.bed)",recovery,"rRNA psi sites in all known chrom21\n")

    #known_data miss/hit hg38.psiU.SingleSites.bed
    ann_data$rRNA_uniq_id<-paste(ann_data$chrom,ann_data$chromStart,ann_data$chromEnd,ann_data$strand,sep="_")
    ann_data_evidence2<-ann_data %>% left_join(evidence2,by=c("rRNA_uniq_id"="rRNA_uniq_id")) %>% filter(!is.na(rRNA_anno))
    write.csv(ann_data_evidence2,paste(outFile_prefix,"_hg38.psiU.SingleSites.bed_psiFinder_hit.csv",sep=""))
    ann_data_no_evidence2<-evidence2 %>% left_join(ann_data,by=c("rRNA_uniq_id"="rRNA_uniq_id")) %>% filter(is.na(class))
    write.csv(ann_data_no_evidence2,paste(outFile_prefix,"_hg38.psiU.SingleSites.bed_psiFinder_miss.csv",sep=""))
    recovery<-paste(round(length(unique(ann_data_evidence2$rRNA_uniq_id))/length(unique(evidence2$rRNA_uniq_id))*100,2),"%",sep="")
    cat("psiFinder recover (hg38.psiU.SingleSites.bed.bed)",recovery,"rRNA psi sites in all known chrom21\n")

    #final_pred miss/hit
    final_pred<-ann_pred_data[ann_pred_data$pred=="psi",]
    final_pred<-final_pred[final_pred$foldChange>2,]
    write.table(final_pred,paste(outFile_prefix, '_ann_psi_prediction.txt', sep=""),sep="\t",row.names=FALSE, col.names=TRUE,quote=F)
    write.table(final_pred,paste(outFile_prefix, '_ann_psi_prediction.bed', sep=""),sep="\t",row.names=FALSE, col.names=FALSE,quote=F)

    final_pred$rRNA_uniq_id<-paste(final_pred$chrom,final_pred$chromStart,final_pred$chromEnd,final_pred$strand,sep="_")
    final_pred_evidence<-final_pred %>% left_join(evidence,by=c("rRNA_uniq_id"="rRNA_uniq_id")) %>% filter(!is.na(rRNA_anno))
    write.csv(final_pred_evidence,paste(outFile_prefix,"_hg38_human_chr21_rRNA_known_pseudoU_SingleSites_psiFinder_ann_hit.csv",sep=""))
    final_pred_no_evidence<-evidence %>% left_join(final_pred,by=c("rRNA_uniq_id"="rRNA_uniq_id")) %>% filter(is.na(base))
    write.csv(final_pred_no_evidence,paste(outFile_prefix,"_hg38_human_chr21_rRNA_known_pseudoU_SingleSites_psiFinder_ann_miss.csv",sep=""))
    recovery<-paste(round(length(unique(final_pred_evidence$rRNA_uniq_id))/length(unique(evidence$rRNA_uniq_id))*100,2),"%",sep="")
    cat("psiFinder+ANN recover (hg38_human_chr21_rRNA_known_pseudoU_SingleSites.bed)",recovery,"rRNA psi sites in all known chrom21\n")

``ann_predict_totalRNA.sh``

.. code:: bash

    #! /bin/bash

    if [[ $# -ne 2 ]]
    then
        echo 'Usage: '$0 ' input ANN_model'
        exit 1
    fi

    input=$1
    ANN_model=$2

    awk 'FS=OFS="\t" {if($1~/^chr[0-9|a-z|A-Z]*$/ && $10=="T" && $14>10){print $0}}' $input > ${input%.bed}_ann_filt.bed
    nohup Rscript $(dirname "$0")/ann_predict_totalRNA.r -f ${input%.bed}_ann_filt.bed -k $ANN_model -r $(dirname "$0")/hg38_human_chr21_rRNA_known_pseudoU_SingleSites.bed -o ${input%.bed} > ${input%.bed}.log
    wait
    cat ${input%.bed}.log

    echo "total RNA: ANN predict program is done!"
    echo -e "total RNA: ANN predict result in $(dirname ${input})"


``ann_predict_totalRNA.r``

.. code:: R

    #!/usr/bin/env Rscript

    options(warn=-1)
    suppressMessages(library("optparse"))
    suppressMessages(library("dplyr"))
    suppressMessages(library("stringr"))
    suppressMessages(library("reshape2"))
    suppressMessages(library("openxlsx"))
    suppressMessages(library("neuralnet"))
    options(warn=-1)

    option_list = list(
      make_option(c("-f", "--filtfile"), type="character", default=NULL,
                  help="filted file of single sites [file]", metavar="character"),
      make_option(c("-k", "--anpsiodelfile"), type="character", default=NULL,
                  help="anpsiodelfile [file]", metavar="character"),
      make_option(c("-r", "--rRNAfile"), type="character", default=NULL,
                  help="rRNA of single sites [file]", metavar="character"),
      make_option(c("-o", "--outfile_prefix"), type="character", default=NULL,
                  help="output file name [default= %default]", metavar="character")
    );
    opt_parser = OptionParser(option_list=option_list);
    opt = parse_args(opt_parser);

    if ( is.null(opt$filtfile) || is.null(opt$anpsiodelfile)|| is.null(opt$rRNAfile) || is.null(opt$outfile_prefix) ){
      print_help(opt_parser);
      stop("Please provide -f anpsiodelfile, -k filtfile, -r rRNAfile, and -o outfile_prefix option", call.=FALSE);
    }

    filtfile = opt$filtfile
    anpsiodelfile = opt$anpsiodelfile
    rRNAfile = opt$rRNAfile
    outFile_prefix = opt$outfile_prefix

    print(filtfile)
    print(anpsiodelfile)
    print(rRNAfile)
    print(outFile_prefix)

    load(file = anpsiodelfile)#"my-ann.RData"

    #filt by ann best model
    cat("\n\n=====================Filt by ann best model=====================\n")
    cat("below is your input data ready to be predicted...\n")
    # get prediction
    to_pred<-read.table(filtfile)
    colnames(to_pred)<-c("chrom",#1
      "chromStart",#2
      "chromEnd",#3
      "name",#4
      "foldChange",#5
      "strand",#6
      "geneName",#7
      "geneStart",#8
      "geneEnd",#9
      "base",#10
      "treatPval",#11
      "ctrlPval",#12
      "minusPval",#13
      "treatStopNum",#14
      "treatStopRPM",#15
      "treatPreStopNum",#16
      "treatAfterStopNum",#17
      "treatReadthroughNum",#18
      "ctrlStopNum",#19
      "ctrlStopRPM",#20
      "ctrlPreStopNum",#21
      "ctrlAfterStopNum",#22
      "ctrlReadthroughNum",#23
      "stopRpmFC",#24
      "treatPreRpmFold",#25
      "ctrlPreRpmFold",#26
      "preRpmFoldRatio",#27
      "treatAfterRpmFold",#28
      "ctrlAfterRpmFold",#29
      "afterRpmFoldRatio",#30
      "treatStopRatio",#31
      "ctrlStopRatio",#32
      "stopRatioFC",#33
      "treatStopMeanNum",#34
      "treatStopMeanFold",#35
      "ctrlStopMeanNum",#36
      "ctrlStopMeanFold",#37
      "treatStopMeanFoldRatio",#38
      "extendSeq")#39
    to_pred$base<-"T"
    to_pred_var<-to_pred %>% select(treatPreRpmFold, preRpmFoldRatio, treatAfterRpmFold, afterRpmFoldRatio, treatStopRatio, stopRatioFC)
    to_pred_var = scale(to_pred_var, center = training_set_mean, scale = training_set_sd)
    str(to_pred_var)

    # Predicting Result
    ann_data_sel.prediction <- neuralnet::compute(ann_data_sel.net, to_pred_var)
    idx <- apply(ann_data_sel.prediction$net.result, 1, which.max)
    net.result <- as.data.frame(ann_data_sel.prediction$net.result)
    colnames(net.result) <- c('psi', 'non_psi')
    pred <- c('psi', 'non_psi')[idx]
    table(pred)

    #add attr to ann_pred_data
    ann_pred_data <- to_pred
    ann_pred_data <- cbind(ann_pred_data, net.result, pred)
    write.xlsx(ann_pred_data, paste(outFile_prefix, '_ann_total_prediction.xlsx', sep = ""), overwrite = TRUE)

    final_pred <- ann_pred_data[ann_pred_data$pred == "psi",]
    print("psi predicted by nn model and fold change (foldChange) >2:")
    final_pred <- final_pred[final_pred$foldChange > 2,]
    print("Final psi prediction:")
    table(final_pred$pred)
    write.xlsx(final_pred, paste(outFile_prefix, '_ann_psi_prediction.xlsx', sep = ""), overwrite = TRUE)
    write.table(final_pred,paste(outFile_prefix, '_ann_psi_prediction.txt', sep=""),sep="\t",row.names=FALSE, col.names=TRUE,quote=F)
    write.table(final_pred,paste(outFile_prefix, '_ann_psi_prediction.bed', sep=""),sep="\t",row.names=FALSE, col.names=FALSE,quote=F)

    #known evidence
    evidence <- read.table(rRNAfile, head = F)
    colnames(evidence) <- c("chrom", "chromStart", "chromEnd", "rRNA_anno", "score", "strand")
    evidence$rRNA_uniq_id <- paste(evidence$chrom, evidence$chromStart, evidence$chromEnd, evidence$strand, sep = "_")
    final_pred$rRNA_uniq_id <- paste(final_pred$chrom, final_pred$chromStart, final_pred$chromEnd, final_pred$strand, sep = "_")
    final_pred_evidence <- final_pred %>% left_join(evidence, by = c("rRNA_uniq_id" = "rRNA_uniq_id")) %>% filter(!is.na(rRNA_anno))
    write.csv(final_pred_evidence, paste(outFile_prefix, "_rtsSeeker_ann_hit.csv", sep = ""))
    final_pred_no_evidence <- evidence %>% left_join(final_pred, by = c("rRNA_uniq_id" = "rRNA_uniq_id")) %>% filter(is.na(base))
    write.csv(final_pred_no_evidence, paste(outFile_prefix, "_rtsSeeker_ann_miss.csv", sep = ""))
    recovery <- paste(round(length(unique(final_pred_evidence$rRNA_uniq_id)) / length(unique(evidence$rRNA_uniq_id)) * 100, 2), "%", sep = "")
    cat("rtsSeeker+ANN recover (hg38_human_chr21_rRNA_known_pseudoU_SingleSites.bed)", recovery, "rRNA psi sites in all known chrom21\n")

.. note:: All user input will be recorded in a plain text file with suffixes ``_ann_totalRNA_build.txt`` and ``_ann_totalRNA_predict.txt`` in psiFinder/config and help users to easily reload the previous config (by simply clicking ``CONFIG`` button).

(3) ROC: Receiver Operating Characteristic
*************************************************
If ``plot ROC`` option is checked,  ``3.1_rtsSeeker_roc.sh`` will run ``roc.r`` to generate file with ``_roc_psi_prediction.bed`` suffix, which is the candidate Ψ-sites gained by ROC best threshold (determined by one of the six variables with the highest F1 score).


(4) User-defined: Customize key thresholds
********************************************
By user-defined thresholds with stringent statistical controls, psiFinder allows users to filter and preserve highly reliable true positive Ψ sites (expertise of threshold determination is needed).

Once click ``FILT``, ``User-defined`` QT widget will run ``user_defined.sh`` and ``user_defined_evaluation.r``

.. image:: /images/User-defined.png

``user_defined.sh``

.. code:: bash

    #! /bin/bash

    if [[ $# -ne 7 ]]
    then
        echo 'Usage: '$0 ' input.bed treatPreRpmFold preRpmFoldRatio treatAfterRpmFold afterRpmFoldRatio treatStopRatio stopRatioFC'
        exit 1
    fi

    out=$1
    treatPreRpmFold=$2
    preRpmFoldRatio=$3
    treatAfterRpmFold=$4
    afterRpmFoldRatio=$5
    treatStopRatio=$6
    stopRatioFC=$7
    outFileName=${out%.bed}

    echo "rtsSeeker User-defined mode: rtsSeeker start..."

    awk 'FS=OFS="\t" {if($1~/^chr[0-9|a-z|A-Z]*$/ && $10=="T" && $14>10){print $0}}' $out > ${outFileName}_user_defined_filt.bed
    echo -e "treatPreRpmFold\tpreRpmFoldRatio\ttreatAfterRpmFold\tafterRpmFoldRatio\ttreatStopRatio\tstopRatioFC\n$treatPreRpmFold\t$preRpmFoldRatio\t$treatAfterRpmFold\t$afterRpmFoldRatio\t$treatStopRatio\t$stopRatioFC" >${outFileName}_user_defined_thres_colname.txt
    echo -e "$treatPreRpmFold\t$preRpmFoldRatio\t$treatAfterRpmFold\t$afterRpmFoldRatio\t$treatStopRatio\t$stopRatioFC" >${outFileName}_user_defined_thres.txt

    cat ${outFileName}_user_defined_filt.bed |awk -v sample=$(basename ${outFileName}) -v treatpre=`(cut -f 1 ${outFileName}_user_defined_thres.txt)` -v preFoldFC=`(cut -f 2 ${outFileName}_user_defined_thres.txt)` -v treataft=`(cut -f 3 ${outFileName}_user_defined_thres.txt)` -v aftFoldFC=`(cut -f 4 ${outFileName}_user_defined_thres.txt)` -v treatstoprate=`(cut -f 5 ${outFileName}_user_defined_thres.txt)` -v stoprateFC=`(cut -f 6 ${outFileName}_user_defined_thres.txt)` 'FS=OFS="\t" {if( $25>=treatpre && $27>=preFoldFC && $28>=treataft && $30>=aftFoldFC && $31>= treatstoprate && $33>=stoprateFC ) {print $0,sample}}' > ${outFileName}_user_defined_psi_prediction.bed
    cat ${outFileName}_user_defined_filt.bed |awk -v treatPreRpmFold=`(cut -f 1 ${outFileName}_user_defined_thres.txt)` -v preRpmFoldRatio=`(cut -f 2 ${outFileName}_user_defined_thres.txt)` -v treatAfterRpmFold=`(cut -f 3 ${outFileName}_user_defined_thres.txt)` -v afterRpmFoldRatio=`(cut -f 4 ${outFileName}_user_defined_thres.txt)` -v treatStopRatio=`(cut -f 5 ${outFileName}_user_defined_thres.txt)` -v stopRatioFC=`(cut -f 6 ${outFileName}_user_defined_thres.txt)` 'FS=OFS="\t" {if( $25>=treatPreRpmFold && $27>=preRpmFoldRatio && $28>=treatAfterRpmFold && $30>=afterRpmFoldRatio && $31>= treatStopRatio && $33>=stopRatioFC) {print $0,"1"} else{ print $0,"0"} }' > ${outFileName}_user_defined_prediction_total.bed

    bedtools intersect -a ${outFileName}_user_defined_prediction_total.bed -b $(dirname "$0")/hg38_human_chr21_rRNA_known_pseudoU_SingleSites.bed -s >${outFileName}_knowpse.bed
    bedtools intersect -a ${outFileName}_user_defined_prediction_total.bed -b $(dirname "$0")/rrna_chr21.bed -s >${outFileName}_rrna.bed

    bedtools intersect -a ${outFileName}_rrna.bed -b $(dirname "$0")/hg38_human_chr21_rRNA_known_pseudoU_SingleSites.bed -s -v >${outFileName}_notpse.bed
    cat ${outFileName}_knowpse.bed |awk 'FS=OFS="\t" {print $0,"1"}' > ${outFileName}_knowpse.txt
    cat ${outFileName}_notpse.bed |awk 'FS=OFS="\t" {print $0,"0"}' > ${outFileName}_notpse.txt

    echo "getting roc_plot.txt"
    cat ${outFileName}_knowpse.txt ${outFileName}_notpse.txt > ${outFileName}_roc_plot.txt

    echo "generate ROC plot..."
    echo "user-defined ROC ploting preparating"
    nohup Rscript $(dirname "$0")/user_defined_evaluation.r -f ${outFileName}_roc_plot.txt -r $(dirname "$0")/hg38_human_chr21_rRNA_known_pseudoU_SingleSites.bed -s $(dirname "$0")/hg38.psiU.SingleSites.bed -t ${outFileName}_user_defined_prediction_total.bed -o ${outFileName} > ${outFileName}_roc_bestthres.log 2>&1 &
    wait

    cat ${outFileName}_roc_bestthres.log
    mupdf-x11 ${outFileName}_six_variables_rRNA_violinplot.pdf &> /dev/null
    mupdf-x11 ${outFileName}_user_defined_evaluation.pdf &> /dev/null

    echo "User-defined program is done!"
    echo -e "User-defined: User-defined result in $(dirname ${outFileName}_user_defined_psi_prediction.bed)"

``user_defined_evaluation.r``

.. code:: R

    #!/usr/bin/env Rscript

    suppressMessages(library("optparse"))
    suppressMessages(library("ggplot2"))
    suppressMessages(library("ggpubr"))
    suppressMessages(library("cowplot"))
    suppressMessages(library("dplyr"))
    suppressMessages(library("gridExtra"))
    suppressMessages(library("reshape2"))
    suppressMessages(library("stringr"))
    suppressMessages(library("RColorBrewer"))

    option_list = list(
      make_option(c("-f", "--rocfile"), type="character", default=NULL,
                  help="ROC of single sites [file]", metavar="character"),
      make_option(c("-r", "--rRNAfile"), type="character", default=NULL,
                  help="hg38_human_chr21_rRNA_known_pseudoU_SingleSites.bed [file]", metavar="character"),
      make_option(c("-s", "--rRNAfile2"), type="character", default=NULL,
                  help="hg38.psiU.SingleSites.bed [file]", metavar="character"),
      make_option(c("-t", "--filtfile"), type="character", default=NULL,
                  help="filt file [file]", metavar="character"),
      make_option(c("-o", "--outfile_prefix"), type="character", default=NULL,
                  help="output file name [default= %default]", metavar="character")
    );
    opt_parser = OptionParser(option_list=option_list);
    opt = parse_args(opt_parser);

    if (is.null(opt$rocfile)|| is.null(opt$filtfile) || is.null(opt$rRNAfile) || is.null(opt$rRNAfile2) || is.null(opt$outfile_prefix) ){
      print_help(opt_parser);
      stop("Please provide -f rocfile, -r rRNAfile, -s rRNAfile2, -t filtfile  and -o outfile_prefix option", call.=FALSE);
    }

    ROCfile = opt$rocfile
    rRNAfile = opt$rRNAfile
    rRNAfile2 = opt$rRNAfile2
    filtfile = opt$filtfile
    outFile_prefix = opt$outfile_prefix

    print(ROCfile)
    print(rRNAfile)
    print(rRNAfile2)
    print(filtfile)
    print(outFile_prefix)

    # roc_plot.txt
    ROC_data<-read.table(ROCfile,head=F)
    colnames(ROC_data)<-c("chrom",#1
      "chromStart",#2
      "chromEnd",#3
      "name",#4
      "foldChange",#5
      "strand",#6
      "geneName",#7
      "geneStart",#8
      "geneEnd",#9
      "base",#10
      "treatPval",#11
      "ctrlPval",#12
      "minusPval",#13
      "treatStopNum",#14
      "treatStopRPM",#15
      "treatPreStopNum",#16
      "treatAfterStopNum",#17
      "treatReadthroughNum",#18
      "ctrlStopNum",#19
      "ctrlStopRPM",#20
      "ctrlPreStopNum",#21
      "ctrlAfterStopNum",#22
      "ctrlReadthroughNum",#23
      "stopRpmFC",#24
      "treatPreRpmFold",#25
      "ctrlPreRpmFold",#26
      "preRpmFoldRatio",#27
      "treatAfterRpmFold",#28
      "ctrlAfterRpmFold",#29
      "afterRpmFoldRatio",#30
      "treatStopRatio",#31
      "ctrlStopRatio",#32
      "stopRatioFC",#33
      "treatStopMeanNum",#34
      "treatStopMeanFold",#35
      "ctrlStopMeanNum",#36
      "ctrlStopMeanFold",#37
      "treatStopMeanFoldRatio",#38
      "extendSeq",#39
      "pred_class",#40
      "real_class")#41
    ROC_data$base<-"T"
    print("real class:\n")
    table(ROC_data$real_class)
    print("prediction class:\n")
    table(ROC_data$pred_class)

    tab <- table(Predicted = ROC_data$pred_class,Actual = ROC_data$real_class)
    tab

    confusionMatrix(factor(ROC_data$pred_class,levels=c("0","1")),factor(ROC_data$real_class,levels=c("0","1")))
    ROC_data_sel<-ROC_data %>% select(treatPreRpmFold,treatAfterRpmFold,treatStopMeanFold,treatStopRatio,preRpmFoldRatio, afterRpmFoldRatio, stopRatioFC,treatStopMeanFoldRatio,real_class)
    print("total summary (rRNA psi and rRNA non-psi)")
    summary(ROC_data_sel)
    print("psi summary (all real rRNA psi)")
    real_rRNA_psi<-ROC_data_sel %>% filter(real_class=="1")
    summary(real_rRNA_psi)
    pdf(paste(outFile_prefix,"_real_rRNA_psi_datadensity.pdf",sep=""))
    datadensity(real_rRNA_psi, lwd = 1,group=cut2(real_rRNA_psi$treatStopRatio,g=2))#cut tretRtsRatio into 2 color group
    dev.off()

    #rRNA-psi-non-psi visualization
    ROC_data_melt<-melt(ROC_data[,c(11:38,41)],id.vars = "real_class")
    ROC_data_melt$real_class<-str_replace(as.character(ROC_data_melt$real_class),"0","non-psi")
    ROC_data_melt$real_class<-str_replace(as.character(ROC_data_melt$real_class),"1","psi")

    data_summary <- function(x) {
       m <- mean(x)
       ymin <- m-sd(x)
       ymax <- m+sd(x)
       return(c(y=m,ymin=ymin,ymax=ymax))
    }

    my_comparisons <- list( c("non-psi", "psi") )

    #treatPreRpmFold
    treatPreRpmFold<-ROC_data_melt%>%filter(str_detect(.$variable,"^treatPreRpmFold$"))
    treatPreRpmFold_bp <- ggplot(treatPreRpmFold, aes(x=real_class, y=log2(value), fill=real_class)) +
      stat_boxplot(geom = "errorbar",
                   width = 0.15) +
      geom_violin(trim=FALSE,alpha=0.8)+
      stat_summary(fun.data=data_summary,geom="crossbar", width=0.2,color="black")+
    labs(x="group", y = "log2(treatPreRpmFold)")
    treatPreRpmFold_bp<-treatPreRpmFold_bp + scale_fill_manual(values=c(brewer.pal(9,"Set1")[2],brewer.pal(9,"Set1")[1])) + theme_classic()+
    theme(text=element_text(size=8),legend.position = "none") + font("xy.text", size = 8)+
    stat_compare_means(method = "t.test",comparisons = my_comparisons,label.y = max(log2(treatPreRpmFold$value))+abs(quantile(log2(treatPreRpmFold$value))[1]))+
    scale_y_continuous(breaks = scales::pretty_breaks(n = 5))+
    scale_y_continuous(expand = expansion(mult = c(0, 0.1)))+
    scale_fill_manual(values=c("#129a92","#1e95d4"))

    #treatAfterRpmFold
    treatAfterRpmFold<-ROC_data_melt%>%filter(str_detect(.$variable,"^treatAfterRpmFold$"))
    treatAfterRpmFold_bp <- ggplot(treatAfterRpmFold, aes(x=real_class, y=log2(value), fill=real_class)) +
      stat_boxplot(geom = "errorbar",
                   width = 0.15) +
      geom_violin(trim=FALSE,alpha=0.8)+
      stat_summary(fun.data=data_summary,geom="crossbar", width=0.2,color="black")+
    labs(x="group", y = "log2(treatAfterRpmFold)")
    treatAfterRpmFold_bp<-treatAfterRpmFold_bp + scale_fill_manual(values=c(brewer.pal(9,"Set1")[2],brewer.pal(9,"Set1")[1])) + theme_classic()+
    theme(text=element_text(size=8), legend.position = "none") + font("xy.text", size = 8)+
    stat_compare_means(method = "t.test",comparisons = my_comparisons,label.y = max(log2(treatAfterRpmFold$value))+abs(quantile(log2(treatAfterRpmFold$value))[1]))+
    scale_y_continuous(breaks = scales::pretty_breaks(n = 5))+
    scale_y_continuous(expand = expansion(mult = c(0, 0.1)))+
    scale_fill_manual(values=c("#129a92","#1e95d4"))

    #preRpmFoldRatio
    preRpmFoldRatio<-ROC_data_melt%>%filter(str_detect(.$variable,"^preRpmFoldRatio$"))
    preRpmFoldRatio_bp <- ggplot(preRpmFoldRatio, aes(x=real_class, y=log2(value), fill=real_class)) +
      stat_boxplot(geom = "errorbar",
                   width = 0.15) +
      geom_violin(trim=FALSE,alpha=0.8)+
      stat_summary(fun.data=data_summary,geom="crossbar", width=0.2,color="black")+
    labs(x="group", y = "log2(preRpmFoldRatio)")
    preRpmFoldRatio_bp<-preRpmFoldRatio_bp + scale_fill_manual(values=c(brewer.pal(9,"Set1")[2],brewer.pal(9,"Set1")[1])) + theme_classic()+
    theme(text=element_text(size=8), legend.position = "none") + font("xy.text", size = 8)+
    stat_compare_means(method = "t.test",comparisons = my_comparisons,label.y = max(log2(preRpmFoldRatio$value))+abs(quantile(log2(preRpmFoldRatio$value))[1]))+
    scale_y_continuous(breaks = scales::pretty_breaks(n = 5))+
    scale_y_continuous(expand = expansion(mult = c(0, 0.1)))+
    scale_fill_manual(values=c("#129a92","#1e95d4"))

    #afterRpmFoldRatio
    afterRpmFoldRatio<-ROC_data_melt%>%filter(str_detect(.$variable,"^afterRpmFoldRatio$"))
    afterRpmFoldRatio_bp <- ggplot(afterRpmFoldRatio, aes(x=real_class, y=log2(value), fill=real_class)) +
      stat_boxplot(geom = "errorbar",
                   width = 0.15) +
      geom_violin(trim=FALSE,alpha=0.8)+
      stat_summary(fun.data=data_summary,geom="crossbar", width=0.2,color="black")+
      labs(x="group", y = "log2(afterRpmFoldRatio)")
    afterRpmFoldRatio_bp<-afterRpmFoldRatio_bp + scale_fill_manual(values=c(brewer.pal(9,"Set1")[2],brewer.pal(9,"Set1")[1])) + theme_classic()+
    theme(text=element_text(size=8),legend.position = "none") + font("xy.text", size = 8)+
    stat_compare_means(method = "t.test",comparisons = my_comparisons,label.y = max(log2(afterRpmFoldRatio$value))+abs(quantile(log2(afterRpmFoldRatio$value))[1]))+
    scale_y_continuous(breaks = scales::pretty_breaks(n = 5))+
    scale_y_continuous(expand = expansion(mult = c(0, 0.1)))+
    scale_fill_manual(values=c("#129a92","#1e95d4"))

    #treatStopRatio
    treatStopRatio<-ROC_data_melt%>%filter(str_detect(.$variable,"^treatStopRatio$"))
    treatStopRatio<-treatStopRatio[treatStopRatio$value!=0,]
    treatStopRatio_bp <- ggplot(treatStopRatio, aes(x=real_class, y=log2(value), fill=real_class)) +
      stat_boxplot(geom = "errorbar",
                   width = 0.15) +
      geom_violin(trim=FALSE,alpha=0.8)+
      stat_summary(fun.data=data_summary,geom="crossbar", width=0.2,color="black")+
      labs(x="group", y = "log2(treatStopRatio)")
    treatStopRatio_bp<-treatStopRatio_bp + scale_fill_manual(values=c(brewer.pal(9,"Set1")[2],brewer.pal(9,"Set1")[1])) + theme_classic()+
    theme(text=element_text(size=8),legend.position = "none") + font("xy.text", size = 8)+
    stat_compare_means(method = "t.test",comparisons = my_comparisons,label.y = max(log2(treatStopRatio$value))+abs(quantile(log2(treatStopRatio$value))[1]))+
    scale_y_continuous(breaks = scales::pretty_breaks(n = 5))+
    scale_y_continuous(expand = expansion(mult = c(0, 0.1)))+
    scale_fill_manual(values=c("#129a92","#1e95d4"))

    stopRatioFC<-ROC_data_melt%>%filter(str_detect(.$variable,"^stopRatioFC$"))
    stopRatioFC_bp <- ggplot(stopRatioFC, aes(x=real_class, y=log2(value), fill=real_class)) +
      stat_boxplot(geom = "errorbar",
                   width = 0.15) +
      geom_violin(trim=FALSE,alpha=0.8)+
      stat_summary(fun.data=data_summary,geom="crossbar", width=0.2,color="black")+
      labs(x="group", y = "log2(stopRatioFC)")
    stopRatioFC_bp<-stopRatioFC_bp + scale_fill_manual(values=c(brewer.pal(9,"Set1")[2],brewer.pal(9,"Set1")[1])) + theme_classic()+
    theme(text=element_text(size=8),legend.position = "none") + font("xy.text", size = 8)+
    stat_compare_means(method = "t.test",comparisons = my_comparisons,label.y = max(log2(stopRatioFC$value))+abs(quantile(log2(stopRatioFC$value))[1]))+
    scale_y_continuous(breaks = scales::pretty_breaks(n = 5))+
    scale_y_continuous(expand = expansion(mult = c(0, 0.1)))+
    scale_fill_manual(values=c("#129a92","#1e95d4"))

    pdf(paste(outFile_prefix,"_six_variables_rRNA_violinplot.pdf",sep=""),width=5,height=6)
    plot_grid(
      treatPreRpmFold_bp,
      preRpmFoldRatio_bp,
      treatAfterRpmFold_bp,
      afterRpmFoldRatio_bp,
      treatStopRatio_bp,
      stopRatioFC_bp,
      align = "hv",
      labels = c('A','B','C','D','E','F'),ncol=2,nrow=3)
    invisible(dev.off())

    #calculate evaluation indicators
    cat("\n\n=====================Calculate evaluation indicators=====================\n")
    confusion_matrix<-as.data.frame(tab)
    confusion_matrix$Predicted<-str_replace(confusion_matrix$Predicted,"1","psi")
    confusion_matrix$Predicted<-str_replace(confusion_matrix$Predicted,"0","non-psi")
    confusion_matrix$Actual<-str_replace(confusion_matrix$Actual,"1","psi")
    confusion_matrix$Actual<-str_replace(confusion_matrix$Actual,"0","non-psi")
    ud_TP <- confusion_matrix[confusion_matrix$Predicted=="psi"&confusion_matrix$Actual=="psi",]$Freq#True Positives (TP)
    ud_FP <- confusion_matrix[confusion_matrix$Predicted=="psi"&confusion_matrix$Actual=="non-psi",]$Freq#False Positives (FP)
    ud_TN <- confusion_matrix[confusion_matrix$Predicted=="non-psi"&confusion_matrix$Actual=="non-psi",]$Freq#True Negatives (TN)
    ud_FN <- confusion_matrix[confusion_matrix$Predicted=="non-psi"&confusion_matrix$Actual=="psi",]$Freq#False Negatives (FN)
    ud_TPR <- ud_TP / (ud_TP + ud_FN)#sensitivity (true positive rate, TPR)
    ud_TNR <- ud_TN / (ud_TN + ud_FP)#specifity (selectivity or true negative rate, TNR)
    ud_FPR <- 1-ud_TNR#False Positive Rate (FPR) (1 - specificit = FP/​N = FP/(TN + FP), FPR)
    ud_FNR <- 1-ud_TPR#False Negative Rate, FNR)
    ud_Prec <- ud_TP / (ud_TP + ud_FP)#Precision
    ud_Recall <- ud_TP / (ud_TP + ud_FN)#Recall
    ud_ACC <- (ud_TP + ud_TN) / (ud_TP + ud_TN + ud_FP + ud_FN)#accuracy
    ud_F1_score <- (2*ud_Recall*ud_Prec) / (ud_Recall + ud_Prec)#F1_score
    eval<-cbind(ud_TP,ud_FP,ud_TN,ud_FN,ud_TPR,ud_TNR,ud_FPR,ud_FNR,ud_Prec,ud_Recall,ud_ACC,ud_F1_score)
    eval<-round(eval,3)
    eval
    write.table(eval,paste(outFile_prefix, '_ud_eval.txt', sep=""),sep="\t",row.names=FALSE, col.names=TRUE,quote=F)

    #show ud evaluation as pdf table
    ud_eval_t_df<-as.data.frame(t(as.data.frame(eval)))
    colnames(ud_eval_t_df)<-"value_or_percentage"
    tt3 <- ttheme_minimal(
      core=list(bg_params = list(fill = blues9[1:4], col=NA),
                fg_params=list(fontface=3)),
      colhead=list(fg_params=list(col="navyblue", fontface=4L)),
      rowhead=list(fg_params=list(col="orange", fontface=3L)))

    pdf(paste(outFile_prefix, '_user_defined_evaluation.pdf', sep=""), width = 7, height = 7) # Open a new pdf file
    grid.arrange(
      tableGrob(ud_eval_t_df, theme=tt3),
      nrow=1)
    invisible(dev.off()) # Close the file

    #filt by best F1 score threshold
    to_filt<-read.table(filtfile,head=F)
    colnames(to_filt)<-c(
      "chrom",#1
      "chromStart",#2
      "chromEnd",#3
      "name",#4
      "foldChange",#5
      "strand",#6
      "geneName",#7
      "geneStart",#8
      "geneEnd",#9
      "base",#10
      "treatPval",#11
      "ctrlPval",#12
      "minusPval",#13
      "treatStopNum",#14
      "treatStopRPM",#15
      "treatPreStopNum",#16
      "treatAfterStopNum",#17
      "treatReadthroughNum",#18
      "ctrlStopNum",#19
      "ctrlStopRPM",#20
      "ctrlPreStopNum",#21
      "ctrlAfterStopNum",#22
      "ctrlReadthroughNum",#23
      "stopRpmFC",#24
      "treatPreRpmFold",#25
      "ctrlPreRpmFold",#26
      "preRpmFoldRatio",#27
      "treatAfterRpmFold",#28
      "ctrlAfterRpmFold",#29
      "afterRpmFoldRatio",#30
      "treatStopRatio",#31
      "ctrlStopRatio",#32
      "stopRatioFC",#33
      "treatStopMeanNum",#34
      "treatStopMeanFold",#35
      "ctrlStopMeanNum",#36
      "ctrlStopMeanFold",#37
      "treatStopMeanFoldRatio",#38
      "extendSeq",#39
      "pred_class")#40
    to_filt$base<-"T"
    to_filt$pred_class<-str_replace(as.character(to_filt$pred_class),"0","non-psi")
    to_filt$pred_class<-str_replace(as.character(to_filt$pred_class),"1","psi")

    table(to_filt$pred_class)
    write.table(to_filt,paste(outFile_prefix, '_user_defined_total_prediction.txt', sep=""),sep="\t",row.names=FALSE, col.names=TRUE,quote=F)
    final_pred<-to_filt %>% filter(pred_class=="psi")
    write.table(final_pred,paste(outFile_prefix, '_ud_psi_prediction.txt', sep=""),sep="\t",row.names=FALSE, col.names=TRUE,quote=F)
    write.table(final_pred,paste(outFile_prefix, '_ud_psi_prediction.bed', sep=""),sep="\t",row.names=FALSE, col.names=FALSE,quote=F)

    #read hg38_human_chr21_rRNA_known_pseudoU_SingleSites.bed evidence
    evidence<-read.table(rRNAfile,head=F)
    colnames(evidence)<-c("chrom","chromStart","chromEnd","rRNA_anno","score","strand")
    evidence$rRNA_uniq_id<-paste(evidence$chrom,evidence$chromStart,evidence$chromEnd,evidence$strand,sep="_")
    ROC_data$rRNA_uniq_id<-paste(ROC_data$chrom,ROC_data$chromStart,ROC_data$chromEnd,ROC_data$strand,sep="_")
    ROC_data_evidence<-ROC_data %>% left_join(evidence,by=c("rRNA_uniq_id"="rRNA_uniq_id")) %>% filter(!is.na(rRNA_anno))
    write.csv(ROC_data_evidence,paste(outFile_prefix,"_hg38_human_chr21_rRNA_known_pseudoU_SingleSites_rtsSeeker_hit.csv",sep=""))
    ROC_data_no_evidence<-evidence %>% left_join(ROC_data,by=c("rRNA_uniq_id"="rRNA_uniq_id")) %>% filter(is.na(class))
    write.csv(ROC_data_no_evidence,paste(outFile_prefix,"_hg38_human_chr21_rRNA_known_pseudoU_SingleSites_rtsSeeker_miss.csv",sep=""))
    recovery<-paste(round(length(unique(ROC_data_evidence$rRNA_uniq_id))/length(unique(evidence$rRNA_uniq_id))*100,2),"%",sep="")
    cat("rtsSeeker recover (hg38_human_chr21_rRNA_known_pseudoU_SingleSites.bed)",recovery,"rRNA psi sites in all known chrom21\n")

    #final_pred miss/hit
    final_pred$rRNA_uniq_id<-paste(final_pred$chrom,final_pred$chromStart,final_pred$chromEnd,final_pred$strand,sep="_")
    final_pred_evidence<-final_pred %>% left_join(evidence,by=c("rRNA_uniq_id"="rRNA_uniq_id")) %>% filter(!is.na(rRNA_anno))
    write.csv(final_pred_evidence,paste(outFile_prefix,"_hg38_human_chr21_rRNA_known_pseudoU_SingleSites_rtsSeeker_ud_thres_hit.csv",sep=""))
    final_pred_no_evidence<-evidence %>% left_join(final_pred,by=c("rRNA_uniq_id"="rRNA_uniq_id")) %>% filter(is.na(base))
    write.csv(final_pred_no_evidence,paste(outFile_prefix,"_hg38_human_chr21_rRNA_known_pseudoU_SingleSites_rtsSeeker_ud_thres_miss.csv",sep=""))
    recovery<-paste(round(length(unique(final_pred_evidence$rRNA_uniq_id))/length(unique(evidence$rRNA_uniq_id))*100,2),"%",sep="")
    cat("rtsSeeker+ud recover (hg38_human_chr21_rRNA_known_pseudoU_SingleSites.bed)",recovery,"rRNA psi sites in all known chrom21\n")

Schematic overview of the six key metrics:

.. image:: /images/principle.png

.. note:: All user input will be recorded in a plain text file with suffixes ``_user_defined_predict.txt`` in psiFinder/config and help users to easily reload the previous config (by simply clicking ``CONFIG`` button).


Output
--------

rtsSeeker result
********************

rtsSeeker result will be well-organized as a bed file with 36 ``(-M 0)`` or 39 ``(-M 1)`` column.

.. code:: bash

    #run with -M 0, rtsSeeker will output below file information
    ##column description
    # 1	#chrom	chr3
    # 2	chromStart	132621249
    # 3	chromEnd	132621250
    # 4	name	rtsSeeker-1
    # 5	foldChange	1
    # 6	strand	+
    # 7	base	G
    # 8	treatPval	-5.01125
    # 9	ctrlPval	-3.06332
    # 10	minusPval	-1.94793
    # 11	treatStopNum	22
    # 12	treatStopRPM	2.62337
    # 13	treatPreStopNum	1
    # 14	treatAfterStopNum	21
    # 15	treatReadthroughNum	22
    # 16	ctrlStopNum	14.33333
    # 17	ctrlStopRPM	2.21748
    # 18	ctrlPreStopNum	1
    # 19	ctrlAfterStopNum	13.6
    # 20	ctrlReadthroughNum	14.33333
    # 21	stopRpmFC	1.18304
    # 22	treatPreRpmFold	22
    # 23	ctrlPreRpmFold	14.33333
    # 24	preRpmFoldRatio	1.53488
    # 25	treatAfterRpmFold	1.04762
    # 26	ctrlAfterRpmFold	1.05392
    # 27	afterRpmFoldRatio	0.99402
    # 28	treatStopRatio	1
    # 29	ctrlStopRatio	1
    # 30	stopRatioFC	1
    # 31	treatStopMeanNum	14
    # 32	treatStopMeanFold	1.57143
    # 33	ctrlStopMeanNum	8.97778
    # 34	ctrlStopMeanFold	1.59653
    # 35	treatStopMeanFoldRatio	0.98427
    # 36	extendSeq	TCTTAGAAAGAGGAGYTTGCCTCCTTAGCGC

    #run with -M 1 --gene $(dirname "$0")/hg38.gencode.v30.tRNA.refseqNcRNA.geneAnno.bed12, rtsSeeker will output below file information
    ##column description
    # 1	#chrom	chr14
    # 2	chromStart	49586636
    # 3	chromEnd	49586637
    # 4	name	rtsSeeker-1
    # 5	foldChange	0.04179
    # 6	strand	+
    # 7	geneName	ENST00000618786.1|RN7SL1-201|ENSG00000276168.1|RN7SL1|misc_RNA
    # 8	geneStart	57
    # 9	geneEnd	58
    # 10	base	G
    # 11	treatPval	0
    # 12	ctrlPval	-2.67E-09
    # 13	minusPval	2.67E-09
    # 14	treatStopNum	2.91439
    # 15	treatStopRPM	0.07151
    # 16	treatPreStopNum	8.32692
    # 17	treatAfterStopNum	6.92437
    # 18	treatReadthroughNum	10574.08742
    # 19	ctrlStopNum	20.62316
    # 20	ctrlStopRPM	0.52566
    # 21	ctrlPreStopNum	17.01255
    # 22	ctrlAfterStopNum	36.85115
    # 23	ctrlReadthroughNum	3126.84178
    # 24	stopRpmFC	0.13604
    # 25	treatPreRpmFold	0.35
    # 26	ctrlPreRpmFold	1.21223
    # 27	preRpmFoldRatio	0.28872
    # 28	treatAfterRpmFold	0.42089
    # 29	ctrlAfterRpmFold	0.55963
    # 30	afterRpmFoldRatio	0.75208
    # 31	treatStopRatio	0.00028
    # 32	ctrlStopRatio	0.0066
    # 33	stopRatioFC	0.04179
    # 34	treatStopMeanNum	186.30904
    # 35	treatStopMeanFold	0.01564
    # 36	ctrlStopMeanNum	117.26259
    # 37	ctrlStopMeanFold	0.17587
    # 38	treatStopMeanFoldRatio	0.08894
    # 39	extendSeq	AGGCTGAGGCTGGAGYATCGCTTGAGTCCAG


ROC evaluation
********************
The ``plot ROC`` function return information and plot for ROC evaluation. File with ``_roc_psi_prediction.bed`` suffix is the filter result by the best ROC threshold (with best F1 score). Running Information will be recorded in a log file, e.g. ``A1_A2_roc_bestthres.log``.

.. code:: bash

    $ cd /the/directory/of/out_file_dir
    $ tree -L 1
    .
    ├── A1_A2_aft_plot.pdf
    ├── A1_A2_aft_plot.png
    ├── A1_A2.bed
    ├── A1_A2_hg38_human_chr21_rRNA_known_pseudoU_SingleSites_rtsSeeker_hit.csv
    ├── A1_A2_hg38_human_chr21_rRNA_known_pseudoU_SingleSites_rtsSeeker_miss.csv
    ├── A1_A2_hg38_human_chr21_rRNA_known_pseudoU_SingleSites_rtsSeeker_roc_stopRatioFC_thres_hit.csv
    ├── A1_A2_hg38_human_chr21_rRNA_known_pseudoU_SingleSites_rtsSeeker_roc_stopRatioFC_thres_miss.csv
    ├── A1_A2_hg38.psiU.SingleSites.bed_rtsSeeker_hit.csv
    ├── A1_A2_hg38.psiU.SingleSites.bed_rtsSeeker_miss.csv
    ├── A1_A2_knowpse.bed
    ├── A1_A2_knowpse.txt
    ├── A1_A2_notpse.bed
    ├── A1_A2_notpse.txt
    ├── A1_A2_pre_plot.pdf
    ├── A1_A2_pre_plot.png
    ├── A1_A2_ratio_plot.pdf
    ├── A1_A2_ratio_plot.png
    ├── A1_A2_roc_afterRpmFoldRatio.pdf
    ├── A1_A2_roc_all.bed
    ├── A1_A2_roc_best_evaluation.pdf
    ├── A1_A2_roc_bestthres_colname.txt
    ├── A1_A2_roc_bestthres.log
    ├── A1_A2_roc_bestthres.txt
    ├── A1_A2_roc_confusion_matrix_and_indicators_arrange.txt
    ├── A1_A2_roc_confusion_matrix_and_indicators.txt
    ├── A1_A2_roc_filt.bed
    ├── A1_A2_roc_plot.txt
    ├── A1_A2_roc_preRpmFoldRatio.pdf
    ├── A1_A2_roc_psi_prediction.bed
    ├── A1_A2_roc_psi_prediction.txt
    ├── A1_A2_roc_stopRatioFC.pdf
    ├── A1_A2_roc_summary.pdf
    ├── A1_A2_roc_total_prediction.txt
    ├── A1_A2_roc_treatAfterRpmFold.pdf
    ├── A1_A2_roc_treatPreRpmFold.pdf
    ├── A1_A2_roc_treatStopMeanFold.pdf
    ├── A1_A2_roc_treatStopMeanFoldRatio.pdf
    ├── A1_A2_roc_treatStopRatio.pdf
    ├── A1_A2_rrna.bed
    ├── A1_A2_rtsSeeker_cmd.log
    ├── A1_A2_six_variables_plot.pdf
    ├── A1_A2_six_variables_rRNA_violinplot.pdf
    ├── A1_A2_StopMeanFold_plot.pdf
    ├── A1_A2_StopMeanFold_plot.png
    ├── A1_A2_stopratio_plot.pdf
    └── A1_A2_stopratio_plot.png

    0 directories, 46 files


.. note:: ``A1_A2_roc_bestthres_colname.txt`` record the best ROC thresholds of the six key metrics, users can apply them to ``User-defined`` QT widget.

.. image:: /images/ROC_violinplot.png

.. image:: /images/ROC_summary.png

.. image:: /images/ROC_best_threshold.png

Model Building
********************

Different identification approaches bring different sensitivity and specificity. For ``SVM/ANN/ROC/User-defined``, ROC evaluation running information will be recorded in a log file, e.g. ``A1_A2_svm_evaluation_totalRNA.log`` and corresponding ROC plot will popup if your input is validated by programs.

A ``SVM`` QT widget model building result example:

.. code:: bash

    $ cd /the/directory/of/out_file_dir
    $ tree -L 1
    .
    ├── A1_A2.bed
    ├── A1_A2_hg38_human_chr21_rRNA_known_pseudoU_SingleSites_rtsSeeker_hit.csv
    ├── A1_A2_hg38_human_chr21_rRNA_known_pseudoU_SingleSites_rtsSeeker_miss.csv
    ├── A1_A2_hg38_human_chr21_rRNA_known_pseudoU_SingleSites_rtsSeeker_svm_hit.csv
    ├── A1_A2_hg38_human_chr21_rRNA_known_pseudoU_SingleSites_rtsSeeker_svm_miss.csv
    ├── A1_A2_hg38.psiU.SingleSites.bed_rtsSeeker_hit.csv
    ├── A1_A2_hg38.psiU.SingleSites.bed_rtsSeeker_miss.csv
    ├── A1_A2_knowpse.bed
    ├── A1_A2_knowpse.txt
    ├── A1_A2_notpse.bed
    ├── A1_A2_notpse.txt
    ├── A1_A2_pca.pdf
    ├── A1_A2_rrna.bed
    ├── A1_A2_six_variables_rRNA_violinplot.pdf
    ├── A1_A2_svm_all.bed
    ├── A1_A2_svm_all_rtsSeeker.cmd
    ├── A1_A2_SVM_eval.txt
    ├── A1_A2_svm_evaluation.pdf
    ├── A1_A2_svm_evaluation_totalRNA.log
    ├── A1_A2_svm_filt_totalRNA.bed
    ├── A1_A2_SVM_model.RData
    ├── A1_A2_SVM_model.rds
    ├── A1_A2_SVM_model.scale
    ├── A1_A2_svm_plot.txt
    ├── A1_A2_svm_psi_prediction.bed
    ├── A1_A2_svm_psi_prediction.txt
    ├── A1_A2_SVM_roc_test_data_plot.pdf
    ├── A1_A2_SVM_test_data.xlsx
    ├── A1_A2_svm_total_prediction.bed
    ├── A1_A2_svm_total_prediction.txt
    └── A1_A2_svm_total_prediction.xlsx

    0 directory, 31 files



.. image:: /images/SVM_roc.png

.. image:: /images/SVM_evaluation.png

Model Prediction
********************

After ``Model Building``, users can use ``PREDICT`` in ``SVM/ANN`` QT widget to predict Ψ-sites by a built ``SVM/ANN model``.

A ``SVM`` QT widget model prediction result example:

.. code:: bash

    $ cd /the/directory/of/out_file_dir
    $ tree -L 1
    .
    ├── A1_A2_svm_filt_totalRNA.bed
    ├── A1_A2_svm_filt_totalRNA_hg38_human_chr21_rRNA_known_pseudoU_SingleSites_rtsSeeker_svm_hit.csv
    ├── A1_A2_svm_filt_totalRNA_hg38_human_chr21_rRNA_known_pseudoU_SingleSites_rtsSeeker_svm_miss.csv
    ├── A1_A2_svm_filt_totalRNA.log
    ├── A1_A2_svm_filt_totalRNA_svm_filt.bed
    ├── A1_A2_svm_filt_totalRNA_svm_psi_prediction.bed
    ├── A1_A2_svm_filt_totalRNA_svm_psi_prediction.txt
    ├── A1_A2_svm_filt_totalRNA_svm_total_prediction.bed
    ├── A1_A2_svm_filt_totalRNA_svm_total_prediction.txt
    └── A1_A2_svm_filt_totalRNA_svm_total_prediction.xlsx

    0 directory, 10 files
