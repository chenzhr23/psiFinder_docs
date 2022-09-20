Ψ-sites Differential Modification
===========================================

.. contents::
    :local:

psiFinder ``Ψ-sites Differential Modification`` utilize `limma <https://www.bioconductor.org/packages/devel/bioc/vignettes/limma/inst/doc/usersguide.pdf>`_ to generate result for Ψ-sites differential modification based on ``stopRatioFC`` of input rtsSeeker result files.

.. image:: /images/Differential_Modification.png

Input
------
Users should choose to upload files (i.e. rtsSeeker result with ``_all.bed`` suffix) to ``pseUdiff`` QT widget in bed format, , to construct Ψ-sites differential modification comparison.

.. note:: **Group1**: normally the control group, input required at least two samples; **Group2**: normally the experiment group , input required at least two samples.

Compare Ψ-sites modification abundance
---------------------------------------------
Once click ``START``, psiFindeer will run ``pseUdiff.sh``, ``pseUdiff.r``, ``limma_bed_annotation.sh``, and ``limma_remove_redundancy.r``.

``pseUdiff.sh``

.. code:: bash

    #! /bin/bash

    if [[ $# -ne 9 ]]
    then
        echo 'Usage: '$0 ' group1_name group2_name group1filelist group2filelist uplog2FC downlog2FC pval mincount outdir'
        exit 1
    fi

    group1_name=$1
    group2_name=$2
    group1filelist=$3
    group2filelist=$4
    uplog2FC=$5
    downlog2FC=$6
    pval=$7
    mincount=$8
    outdir=$9


    Rscript $(dirname "$0")/pseudiff.r -d $group1_name -e $group2_name -f $group1filelist -g $group2filelist -l $mincount -o $outdir > "${outdir}/${group2_name}-${group1_name}"_limma.log
    cat "${outdir}/${group2_name}-${group1_name}"_limma.log
    bash $(dirname "$0")/limma_bed_annotation.sh $uplog2FC $downlog2FC $pval "${group2_name}-${group1_name}" "${outdir}/${group2_name}-${group1_name}"_treat_control_overall_filt.bed $(dirname "$0")/hg38.genecode.v30.tRNA.snoRNA.miRNA.rmsk.exonFeatures.bed6 $group1filelist $group2filelist > "${outdir}/${group2_name}-${group1_name}"_limma_bed_annotation.log
    cat "${outdir}/${group2_name}-${group1_name}"_limma_bed_annotation.log

``pseUdiff.r``

.. code:: R

    #!/usr/bin/env Rscript
    options(warn=-1)
    suppressMessages(library("optparse"))
    suppressMessages(library("openxlsx"))
    suppressMessages(library("limma"))
    suppressMessages(library("edgeR"))
    suppressMessages(library("dplyr"))
    suppressMessages(library("tidyverse"))

    options(warn=-1)

    option_list = list(
      make_option(c("-d", "--group1_name"), type="character", default=NULL,
                  help="file name of group1 for limma [file]", metavar="character"),
      make_option(c("-e", "--group2_name"), type="character", default=NULL,
                  help="file name of group1 for limma [file]", metavar="character"),
      make_option(c("-f", "--group1filelist"), type="character", default=NULL,
                  help="file list of group1 for limma [file]", metavar="character"),
      make_option(c("-g", "--group2filelist"), type="character", default=NULL,
                  help="file list of group2 for limma [file]", metavar="character"),
      make_option(c("-l", "--mincount"), type="numeric", default=NULL,
                  help="mincount for limma [file]", metavar="numeric"),
      make_option(c("-o", "--outdir"), type="character", default=NULL,
                  help="output file name [default= %default]", metavar="character")
    );
    opt_parser = OptionParser(option_list=option_list);
    opt = parse_args(opt_parser);

    if (is.null(opt$group1_name) || is.null(opt$group2_name) || is.null(opt$group1filelist) || is.null(opt$group2filelist)  || is.null(opt$mincount) || is.null(opt$outdir) ){
      print_help(opt_parser);
      stop("Please provide -d group1_name, -e group2_name, -f group1filelist, -g group2filelist, -l mincount, and -o outdir option", call.=FALSE);
    }

    group1_name = opt$group1_name
    group2_name = opt$group2_name
    group1filelist = opt$group1filelist
    group1filelist<-unlist(strsplit(group1filelist,","))
    group2filelist = opt$group2filelist
    group2filelist<-unlist(strsplit(group2filelist,","))
    mincount = opt$mincount
    outdir = opt$outdir

    print(paste("name for group1:",group1_name,sep=""))
    print(paste("name for group2:",group2_name,sep=""))
    print(paste("control group (group1) list:",group1filelist,sep=""))
    print(paste("experiment group (group2) list:",group2filelist,sep=""))
    print(paste("down-regulated log2(fold change):",downlog2FC,sep=""))
    print(paste("minimum count for each group (replicates):",mincount,sep=""))
    print(paste("output directory",outdir,sep=""))
    print("limma program start")

    ####group1####
    if(!is.null(group1filelist)){
      print("getting group 1 foldChange")
      datalist = lapply(group1filelist, function(x)read.delim(x, header=T))
      names(datalist)<-group1filelist
      if(grepl("*.bed",basename(group1filelist))){
        print("input files are bed format")
        invisible(lapply(seq_along(datalist),function(x){
        colnames(datalist[[x]])[1]<<-"chrom"
      }))
      }else{
        print("input files are txt format")
      }
      invisible(lapply(seq_along(datalist),function(x){
        datalist[[x]]$base<<-str_replace_all(datalist[[x]]$base,"TRUE","T")
      }))
      foldChange_datalist<-list()
      invisible(lapply(seq_along(datalist),function(x){
        foldChange_datalist[[x]]<<-datalist[[x]] %>% filter(base=="T") %>% dplyr::select(chrom,chromStart,chromEnd,name,foldChange,treatStopNum,treatPreRpmFold,preRpmFoldRatio,treatAfterRpmFold,afterRpmFoldRatio,treatStopRatio,stopRatioFC,strand,extendSeq)
      }))
      #give names to datalist
      names(foldChange_datalist)<-group1filelist
      foldChange_datalist_df<-foldChange_datalist %>% reduce(full_join, by=c("extendSeq"))
      colnames(foldChange_datalist_df)

      foldChange_datalist_df<-foldChange_datalist_df %>% distinct(extendSeq,.keep_all = TRUE)
      col_name<-colnames(foldChange_datalist_df)

      chrom_col_name<-col_name[which(str_detect(col_name,"^chrom$|chrom[.]"))]
      chrom_col_name_list<-as.list(foldChange_datalist_df[,chrom_col_name])

      chromStart_col_name<-col_name[which(str_detect(col_name,"^chromStart$|chromStart[.]"))]
      chromStart_col_name_list<-as.list(foldChange_datalist_df[,chromStart_col_name])

      chromEnd_col_name<-col_name[which(str_detect(col_name,"^chromEnd$|chromEnd[.]"))]
      chromEnd_col_name_list<-as.list(foldChange_datalist_df[,chromEnd_col_name])

      strand_col_name<-col_name[which(str_detect(col_name,"^strand$|strand[.]"))]
      strand_col_name_list<-as.list(foldChange_datalist_df[,strand_col_name])

      foldChange_datalist_df<-foldChange_datalist_df %>% mutate(chrom.x =do.call(coalesce,chrom_col_name_list)) %>% mutate(chromStart.x =do.call(coalesce,chromStart_col_name_list)) %>% mutate(chromEnd.x =do.call(coalesce,chromEnd_col_name_list)) %>% mutate(strand.x =do.call(coalesce,strand_col_name_list))

      foldChange_datalist_df_chrloc<-foldChange_datalist_df %>% dplyr::select(chrom.x,chromStart.x,chromEnd.x,strand.x)
      foldChange_datalist_df_foldChange<-foldChange_datalist_df %>% dplyr::select(contains("foldChange"))
      foldChange_datalist_df_treatStopNum<-foldChange_datalist_df %>% dplyr::select(contains("treatStopNum"))
      foldChange_datalist_df_treatPreRpmFold<-foldChange_datalist_df %>% dplyr::select(contains("treatPreRpmFold"))
      foldChange_datalist_df_preRpmFoldRatio<-foldChange_datalist_df %>% dplyr::select(contains("preRpmFoldRatio"))
      foldChange_datalist_df_treatAfterRpmFold<-foldChange_datalist_df %>% dplyr::select(contains("treatAfterRpmFold"))
      foldChange_datalist_df_afterRpmFoldRatio<-foldChange_datalist_df %>% dplyr::select(contains("afterRpmFoldRatio"))
      foldChange_datalist_df_treatStopRatio<-foldChange_datalist_df %>% dplyr::select(contains("treatStopRatio"))
      foldChange_datalist_df_stopRatioFC<-foldChange_datalist_df %>% dplyr::select(contains("stopRatioFC"))
      foldChange_datalist_df_extendSeq<-foldChange_datalist_df %>% dplyr::select(contains("extendSeq"))
      foldChange_datalist_df_foldChange<-cbind(foldChange_datalist_df_chrloc,foldChange_datalist_df_extendSeq,foldChange_datalist_df_foldChange,foldChange_datalist_df_treatStopNum,foldChange_datalist_df_treatPreRpmFold,foldChange_datalist_df_preRpmFoldRatio,foldChange_datalist_df_treatAfterRpmFold,foldChange_datalist_df_afterRpmFoldRatio,foldChange_datalist_df_treatStopRatio,foldChange_datalist_df_stopRatioFC)

      colnames(foldChange_datalist_df_foldChange)[c(6:length(colnames(foldChange_datalist_df_foldChange)))]<-paste(gsub(".[xy]","",colnames(foldChange_datalist_df_foldChange)[c(6:length(colnames(foldChange_datalist_df_foldChange)))]),rep(paste("rep",1:length(group1filelist),sep=""),1,each=1),sep="_")
      foldChange_datalist_df_foldChange[is.na(foldChange_datalist_df_foldChange)] = 0
      group1_overall<-foldChange_datalist_df_foldChange
      colnames(group1_overall)[1:4]<-c("chrom","chromStart","chromEnd","strand")

    }

    ####group2####
    if(!is.null(group2filelist)){
      print("getting group 2 foldChange")
      datalist = lapply(group2filelist, function(x)read.delim(x, header=T))
      names(datalist)<-group2filelist
      if(grepl("*.bed",basename(group1filelist))){
        print("input files are bed format")
        invisible(lapply(seq_along(datalist),function(x){
        colnames(datalist[[x]])[1]<<-"chrom"
      }))
      }else{
        print("input files are txt format")
      }
      invisible(lapply(seq_along(datalist),function(x){
        datalist[[x]]$base<<-str_replace_all(datalist[[x]]$base,"TRUE","T")
      }))
      foldChange_datalist<-list()
      invisible(lapply(seq_along(datalist),function(x){
        foldChange_datalist[[x]]<<-datalist[[x]] %>% filter(base=="T") %>% dplyr::select(chrom,chromStart,chromEnd,name,foldChange,treatStopNum,treatPreRpmFold,preRpmFoldRatio,treatAfterRpmFold,afterRpmFoldRatio,treatStopRatio,stopRatioFC,strand,extendSeq)
      }))
      #give names to datalist
      names(foldChange_datalist)<-group2filelist
      foldChange_datalist_df<-foldChange_datalist %>% reduce(full_join, by=c("extendSeq"))
      colnames(foldChange_datalist_df)

      foldChange_datalist_df<-foldChange_datalist_df %>% distinct(extendSeq,.keep_all = TRUE)
      col_name<-colnames(foldChange_datalist_df)

      chrom_col_name<-col_name[which(str_detect(col_name,"^chrom$|chrom[.]"))]
      chrom_col_name_list<-as.list(foldChange_datalist_df[,chrom_col_name])

      chromStart_col_name<-col_name[which(str_detect(col_name,"^chromStart$|chromStart[.]"))]
      chromStart_col_name_list<-as.list(foldChange_datalist_df[,chromStart_col_name])

      chromEnd_col_name<-col_name[which(str_detect(col_name,"^chromEnd$|chromEnd[.]"))]
      chromEnd_col_name_list<-as.list(foldChange_datalist_df[,chromEnd_col_name])

      strand_col_name<-col_name[which(str_detect(col_name,"^strand$|strand[.]"))]
      strand_col_name_list<-as.list(foldChange_datalist_df[,strand_col_name])

      foldChange_datalist_df<-foldChange_datalist_df %>% mutate(chrom.x =do.call(coalesce,chrom_col_name_list)) %>% mutate(chromStart.x =do.call(coalesce,chromStart_col_name_list)) %>% mutate(chromEnd.x =do.call(coalesce,chromEnd_col_name_list)) %>% mutate(strand.x =do.call(coalesce,strand_col_name_list))

      foldChange_datalist_df_chrloc<-foldChange_datalist_df %>% dplyr::select(chrom.x,chromStart.x,chromEnd.x,strand.x)
      foldChange_datalist_df_foldChange<-foldChange_datalist_df %>% dplyr::select(contains("foldChange"))
      foldChange_datalist_df_treatStopNum<-foldChange_datalist_df %>% dplyr::select(contains("treatStopNum"))
      foldChange_datalist_df_treatPreRpmFold<-foldChange_datalist_df %>% dplyr::select(contains("treatPreRpmFold"))
      foldChange_datalist_df_preRpmFoldRatio<-foldChange_datalist_df %>% dplyr::select(contains("preRpmFoldRatio"))
      foldChange_datalist_df_treatAfterRpmFold<-foldChange_datalist_df %>% dplyr::select(contains("treatAfterRpmFold"))
      foldChange_datalist_df_afterRpmFoldRatio<-foldChange_datalist_df %>% dplyr::select(contains("afterRpmFoldRatio"))
      foldChange_datalist_df_treatStopRatio<-foldChange_datalist_df %>% dplyr::select(contains("treatStopRatio"))
      foldChange_datalist_df_stopRatioFC<-foldChange_datalist_df %>% dplyr::select(contains("stopRatioFC"))
      foldChange_datalist_df_extendSeq<-foldChange_datalist_df %>% dplyr::select(contains("extendSeq"))
      foldChange_datalist_df_foldChange<-cbind(foldChange_datalist_df_chrloc,foldChange_datalist_df_extendSeq,foldChange_datalist_df_foldChange,foldChange_datalist_df_treatStopNum,foldChange_datalist_df_treatPreRpmFold,foldChange_datalist_df_preRpmFoldRatio,foldChange_datalist_df_treatAfterRpmFold,foldChange_datalist_df_afterRpmFoldRatio,foldChange_datalist_df_treatStopRatio,foldChange_datalist_df_stopRatioFC)

      colnames(foldChange_datalist_df_foldChange)[c(6:length(colnames(foldChange_datalist_df_foldChange)))]<-paste(gsub(".[xy]","",colnames(foldChange_datalist_df_foldChange)[c(6:length(colnames(foldChange_datalist_df_foldChange)))]),rep(paste("rep",1:length(group2filelist),sep=""),1,each=1),sep="_")
      foldChange_datalist_df_foldChange[is.na(foldChange_datalist_df_foldChange)] = 0
      group2_overall<-foldChange_datalist_df_foldChange
      colnames(group2_overall)[1:4]<-c("chrom","chromStart","chromEnd","strand")

    }

    if(!is.null(group1_overall) & !is.null(group2_overall)){
      print("calculating contrast significance of group2/group1 using limma voom")
      group_list<-list(group1_overall,group2_overall)
      names(group_list)<-c("group1","group2")
      foldChange_group_list_df<-group_list %>% reduce(full_join, by=c("extendSeq"))
      colnames(foldChange_group_list_df)

      foldChange_group_list_df<-foldChange_group_list_df %>% distinct(extendSeq,.keep_all = TRUE)
      foldChange_group_list_df<-foldChange_group_list_df %>% mutate(chrom.x = coalesce(chrom.x,chrom.y)) %>% mutate(chromStart.x = coalesce(chromStart.x,chromStart.y)) %>% mutate(chromEnd.x = coalesce(chromEnd.x,chromEnd.y)) %>% mutate(strand.x = coalesce(strand.x,strand.y))

      foldChange_group_list_df_chrloc<-foldChange_group_list_df %>% dplyr::select(chrom.x,chromStart.x,chromEnd.x,strand.x)
      foldChange_group_list_df_foldChange<-foldChange_group_list_df %>% dplyr::select(contains("foldChange"))

      #treatStopNum
      foldChange_group_list_df_treatStopNum<-foldChange_group_list_df %>% dplyr::select(contains("treatStopNum"))
      foldChange_group_list_df_treatStopNum[is.na(foldChange_group_list_df_treatStopNum)] = 0
      group1_treatStopNumBaseMean<-rowMeans(foldChange_group_list_df_treatStopNum[,1:length(group1filelist)])
      group2_treatStopNumBaseMean<-rowMeans(foldChange_group_list_df_treatStopNum[,(length(group1filelist)+1):(length(group1filelist)+length(group2filelist))])

      #treatPreRpmFold
      foldChange_group_list_df_treatPreRpmFold<-foldChange_group_list_df %>% dplyr::select(contains("treatPreRpmFold"))
      foldChange_group_list_df_treatPreRpmFold[is.na(foldChange_group_list_df_treatPreRpmFold)] = 0
      group1_treatPreRpmFoldBaseMean<-rowMeans(foldChange_group_list_df_treatPreRpmFold[,1:length(group1filelist)])
      group2_treatPreRpmFoldBaseMean<-rowMeans(foldChange_group_list_df_treatPreRpmFold[,(length(group1filelist)+1):(length(group1filelist)+length(group2filelist))])

      #preRpmFoldRatio
      foldChange_group_list_df_preRpmFoldRatio<-foldChange_group_list_df %>% dplyr::select(contains("preRpmFoldRatio"))
      foldChange_group_list_df_preRpmFoldRatio[is.na(foldChange_group_list_df_preRpmFoldRatio)] = 0
      group1_preRpmFoldRatioBaseMean<-rowMeans(foldChange_group_list_df_preRpmFoldRatio[,1:length(group1filelist)])
      group2_preRpmFoldRatioBaseMean<-rowMeans(foldChange_group_list_df_preRpmFoldRatio[,(length(group1filelist)+1):(length(group1filelist)+length(group2filelist))])

      #treatAfterRpmFold
      foldChange_group_list_df_treatAfterRpmFold<-foldChange_group_list_df %>% dplyr::select(contains("treatAfterRpmFold"))
      foldChange_group_list_df_treatAfterRpmFold[is.na(foldChange_group_list_df_treatAfterRpmFold)] = 0
      group1_treatAfterRpmFoldBaseMean<-rowMeans(foldChange_group_list_df_treatAfterRpmFold[,1:length(group1filelist)])
      group2_treatAfterRpmFoldBaseMean<-rowMeans(foldChange_group_list_df_treatAfterRpmFold[,(length(group1filelist)+1):(length(group1filelist)+length(group2filelist))])

      #afterRpmFoldRatio
      foldChange_group_list_df_afterRpmFoldRatio<-foldChange_group_list_df %>% dplyr::select(contains("afterRpmFoldRatio"))
      foldChange_group_list_df_afterRpmFoldRatio[is.na(foldChange_group_list_df_afterRpmFoldRatio)] = 0
      group1_afterRpmFoldRatioBaseMean<-rowMeans(foldChange_group_list_df_afterRpmFoldRatio[,1:length(group1filelist)])
      group2_afterRpmFoldRatioBaseMean<-rowMeans(foldChange_group_list_df_afterRpmFoldRatio[,(length(group1filelist)+1):(length(group1filelist)+length(group2filelist))])

      #treatStopRatio
      foldChange_group_list_df_treatStopRatio<-foldChange_group_list_df %>% dplyr::select(contains("treatStopRatio"))
      foldChange_group_list_df_treatStopRatio[is.na(foldChange_group_list_df_treatStopRatio)] = 0
      group1_treatStopRatioBaseMean<-rowMeans(foldChange_group_list_df_treatStopRatio[,1:length(group1filelist)])
      group2_treatStopRatioBaseMean<-rowMeans(foldChange_group_list_df_treatStopRatio[,(length(group1filelist)+1):(length(group1filelist)+length(group2filelist))])

      #stopRatioFC
      foldChange_group_list_df_stopRatioFC<-foldChange_group_list_df %>% dplyr::select(contains("stopRatioFC"))
      foldChange_group_list_df_stopRatioFC[is.na(foldChange_group_list_df_stopRatioFC)] = 0
      group1_stopRatioFCBaseMean<-rowMeans(foldChange_group_list_df_stopRatioFC[,1:length(group1filelist)])
      group2_stopRatioFCBaseMean<-rowMeans(foldChange_group_list_df_stopRatioFC[,(length(group1filelist)+1):(length(group1filelist)+length(group2filelist))])

      foldChange_group_list_df_extendSeq<-foldChange_group_list_df %>% dplyr::select(contains("extendSeq"))
      foldChange_group_list_df_foldChange<-cbind(foldChange_group_list_df_chrloc,foldChange_group_list_df_extendSeq,foldChange_group_list_df_foldChange)

      colnames(foldChange_group_list_df_foldChange)[c(6:length(colnames(foldChange_group_list_df_foldChange)))]<-paste(gsub(".[xy]","",colnames(foldChange_group_list_df_foldChange)[c(6:length(colnames(foldChange_group_list_df_foldChange)))]),c(rep(paste("group",1,sep=""),times=length(group1filelist)),rep(paste("group",2,sep=""),times=length(group2filelist))),sep="_")
      foldChange_group_list_df_foldChange[is.na(foldChange_group_list_df_foldChange)] = 0

      colnames(foldChange_group_list_df_foldChange)[1:4]<-c("chrom","chromStart","chromEnd","strand")

      #construct sample design
      design <- model.matrix(~ 0+factor(c(rep(1,length(group1filelist)),rep(2,length(group2filelist)))))
      colnames(design) <- c("group1", "group2")
      foldChange<-foldChange_group_list_df_foldChange[,6:length(foldChange_group_list_df_foldChange)]
      rownames(design) <- colnames(foldChange)

      foldChange<-foldChange+0.001
      log2foldChange <- log2(foldChange)
      fit <- lmFit(log2foldChange, design)
      cont.matrix <- makeContrasts(group2-group1, levels=design)
      fit2 <- contrasts.fit(fit, cont.matrix)
      fit2 <- eBayes(fit2,trend=TRUE,robust=TRUE)

      ####overall results####
      overall=topTable(fit2, n=Inf, sort.by="none")
      overall<-cbind(foldChange_group_list_df_foldChange,overall)
      overall<-cbind(overall,group1_treatStopNumBaseMean,group2_treatStopNumBaseMean,group1_treatPreRpmFoldBaseMean,group2_treatPreRpmFoldBaseMean,group1_preRpmFoldRatioBaseMean,group2_preRpmFoldRatioBaseMean,group1_treatAfterRpmFoldBaseMean,group2_treatAfterRpmFoldBaseMean,group1_afterRpmFoldRatioBaseMean,group2_afterRpmFoldRatioBaseMean,group1_treatStopRatioBaseMean,group2_treatStopRatioBaseMean,group1_stopRatioFCBaseMean,group2_stopRatioFCBaseMean)
      write.xlsx(overall,paste(outdir,"/",group2_name,"-",group1_name,"_treat_control_overall.xlsx",sep=""),overwrite =T)

      overall_filt<-overall %>%
        rowwise() %>%
        mutate(
          aux_condition1 = sum(c_across(cols = ends_with("_group1")) > mincount),#group1 aux_condition1
          aux_condition2 = sum(c_across(cols = ends_with("_group1")) == 0),#group1 aux_condition2
          aux_condition3 = sum(c_across(cols = ends_with("_group2")) > mincount),#group2 aux_condition3
          aux_condition4 = sum(c_across(cols = ends_with("_group2")) == 0)#group2 aux_condition4
        ) %>%
        ungroup() %>%
        filter((aux_condition1 == length(group1filelist) & aux_condition3 == length(group2filelist))|
               (aux_condition1 == length(group1filelist) & aux_condition4 == length(group2filelist))|
               # (aux_condition3 == length(group1filelist) & aux_condition1 == length(group2filelist))|
               (aux_condition3 == length(group1filelist) & aux_condition2 == length(group2filelist))) %>% dplyr::select(-starts_with("aux_"))
      write.xlsx(overall_filt,paste(outdir,"/",group2_name,"-",group1_name,"_treat_control_overall_filt.xlsx",sep=""),overwrite =T)

      overall_filt$name<-paste("group2-group1-",1:length(overall_filt$chrom),sep="")
      overall_filt_bed<-overall_filt[,c(1:3,length(overall_filt),5+length(group1filelist)+length(group2filelist)+1,4:5,6:(6+length(group1filelist)+length(group2filelist)-1),(5+length(group1filelist)+length(group2filelist)+1+1):(length(overall_filt)-1))]
      write.table(overall_filt_bed,paste(outdir,"/",group2_name,"-",group1_name,"_treat_control_overall_filt.bed",sep=""),row.names = F,quote = F,col.names=F,sep="\t")
    }

    print("limma program end")

``limma_bed_annotation.sh``

.. code:: bash

    #! /bin/bash

    if [[ $# -ne 8 ]]
    then
        echo 'Usage: ./'$0 ' uplog2FC downlog2FC pval group_name bedfile bed6file group1filelist group2filelist'
        exit 1
    fi

    uplog2FC=$1
    downlog2FC=$2
    pval=$3
    group_name=$4
    bedfile=$5
    bed6file=$6
    group1filelist=$7
    group2filelist=$8

    echo "generating ${bedfile%.bed}_anno.bed"
    cut -f 1-6 ${bedfile} > ${bedfile}6
    bedAnnotator_cmd1="$(dirname "$0")/bedAnnotator -s 1 --anno $bed6file --bed ${bedfile}6 -o ${bedfile%.bed}_anno.bed"
    bedAnnotator_cmd2="$(dirname "$0")/bedAnnotator -s 1 --anno $bed6file --bed $bedfile -o ${bedfile%.bed}_anno_append.bed"
    eval $bedAnnotator_cmd1 &
    eval $bedAnnotator_cmd2 &
    wait

    #add sequence and remove redundancy
    Rscript $(dirname "$0")/limma_remove_redundancy.r -l $uplog2FC -m $downlog2FC -n $pval -q $group_name -i $group1filelist -j $group2filelist -f ${bedfile%.bed}_anno.bed -g $bedfile -s $(dirname "$0")/hg38.psiU.SingleSites.bed -e $(dirname "$0")/human.hg38.Pseudo.result.col29.xlsx -o ${bedfile%.bed}
    echo -e "Finished: bedAnnotator done\n"

    echo -e "bedAnnotator result in $(dirname ${bedfile%.bed})"

    mupdf-x11 ${bedfile%.bed}_DGE_volcano_plot.pdf &> /dev/null
    mupdf-x11 ${bedfile%.bed}_gene_biotype_piechart.pdf &> /dev/null
    mupdf-x11 ${bedfile%.bed}_plot_table.pdf &> /dev/null

``limma_remove_redundancy.r``

.. code:: R

    #!/usr/bin/env Rscript
    options(warn=-1)
    suppressMessages(library("optparse"))
    suppressMessages(library("dplyr"))
    suppressMessages(library("stringr"))
    suppressMessages(library("ggrepel"))
    suppressMessages(library("openxlsx"))
    suppressMessages(library("ggplot2"))
    suppressMessages(library("RColorBrewer"))
    suppressMessages(library("tidyr"))
    suppressMessages(library("gridExtra"))


    option_list = list(
      make_option(c("-l", "--uplog2FC"), type="numeric", default=NULL,
                  help="uplog2FC", metavar="numeric"),
      make_option(c("-m", "--downlog2FC"), type="numeric", default=NULL,
                  help="downlog2FC", metavar="numeric"),
      make_option(c("-n", "--pval"), type="numeric", default=NULL,
                  help="pval", metavar="numeric"),
      make_option(c("-q", "--group_name"), type="character", default=NULL,
                  help="group_name", metavar="character"),
      make_option(c("-i", "--group1filelist"), type="character", default=NULL,
                  help="group1filelist", metavar="character"),
      make_option(c("-j", "--group2filelist"), type="character", default=NULL,
                  help="group2filelist", metavar="character"),
      make_option(c("-f", "--infile1"), type="character", default=NULL,
                  help="input file 1 [file]", metavar="character"),
      make_option(c("-g", "--infile2"), type="character", default=NULL,
                  help="input file 2 [file]", metavar="character"),
      make_option(c("-s", "--infile3"), type="character", default=NULL,
                  help="input file 3 [file]", metavar="character"),
      make_option(c("-e", "--infile4"), type="character", default=NULL,
                  help="input file 4 [file]", metavar="character"),
      make_option(c("-o", "--outfile"), type="character", default=NULL,
                  help="output file name [default= %default]", metavar="character")
    );
    opt_parser = OptionParser(option_list=option_list);
    opt = parse_args(opt_parser);

    if (is.null(opt$uplog2FC) || is.null(opt$downlog2FC) || is.null(opt$pval) ||is.null(opt$group_name) || is.null(opt$group1filelist) || is.null(opt$group2filelist) ||is.null(opt$infile1) || is.null(opt$infile2) ||is.null(opt$infile3) ||is.null(opt$infile4) ||is.null(opt$outfile) ){
      print_help(opt_parser);
      stop("Please provide -l uplog2FC, -m downlog2FC, -n pval, -g group_name, -i group1filelist, -j group2filelist, -f infile1 -g infile2 and -s infile3, -e infile4, and -o outfile option", call.=FALSE);
    }

    uplog2FC = opt$uplog2FC
    downlog2FC = opt$downlog2FC
    pval = opt$pval
    group_name = opt$group_name
    group1filelist = opt$group1filelist
    group1filelist<-unlist(strsplit(group1filelist,","))
    group2filelist = opt$group2filelist
    group2filelist<-unlist(strsplit(group2filelist,","))
    anno.biotype.bed.file = opt$infile1
    input.bed.file = opt$infile2
    hg38.psiU.SingleSites.bed.file = opt$infile3
    human.hg38.Pseudo.result.col29.xlsx.file = opt$infile4
    outFile = opt$outfile

    print(uplog2FC)
    print(downlog2FC)
    print(pval)
    print(group_name)
    print(group1filelist)
    print(group2filelist)
    print(anno.biotype.bed.file)
    print(input.bed.file)
    print(hg38.psiU.SingleSites.bed.file)
    print(human.hg38.Pseudo.result.col29.xlsx.file)
    print(outFile)

    anno.biotype.bed<-read.delim(anno.biotype.bed.file,sep="\t",head=F)
    anno.biotype.bed<-separate(data = anno.biotype.bed, col = V7, into = c("Transcript_id", "Transcript_name","Gene_id","Gene_name","type","Region"), sep = "\\|")

    #pie chart for all snoRNA-guided RNAs
    mRNA<-data.frame(gene_biotype="mRNA",feature_type=c("protein_coding"))
    pseudogene<-data.frame(gene_biotype="pseudogene",feature_type=c("rRNA_pseudogene", "transcribed_unprocessed_pseudogene", "translated_processed_pseudogene", "unprocessed_pseudogene", "processed_pseudogene", "transcribed_processed_pseudogene", "unitary_pseudogene", "pseudogene", "polymorphic_pseudogene", "transcribed_unitary_pseudogene", "TR_V_pseudogene", "TR_J_pseudogene", "IG_V_pseudogene", "IG_C_pseudogene", "IG_D_pseudogene", "IG_pseudogene", "IG_J_pseudogene"))
    lncRNA<-data.frame(gene_biotype="lncRNA",feature_type=c("lncRNA","processed_transcript","lincRNA","non_coding","3prime_overlapping_ncRNA","3prime_overlapping_ncrna","sense_intronic","antisense","sense_overlapping","known_ncrna","macro_lncRNA","bidirectional_promoter_lncRNA","retained_intron","TEC"))
    sncRNA<-data.frame(gene_biotype="sncRNA",feature_type=c("snRNA","snoRNA","misc_RNA","miscRNA","miRNA","ribozyme","sRNA","scRNA","scaRNA","srpRNA","tRNA-Deu","tRNA-RTE","piRNA","siRNA"))
    rRNA<-data.frame(gene_biotype="rRNA",feature_type=c("rRNA","Mt_rRNA"))
    tRNA<-data.frame(gene_biotype="tRNA",feature_type=c("tRNA","Mt_tRNA","vaultRNA"))
    IG_gene<-data.frame(gene_biotype="IG_gene",feature_type=c("IG_C_gene","IG_D_gene","IG_J_gene","IG_LV_gene","IG_V_gene"))
    TR_gene<-data.frame(gene_biotype="TR_gene",feature_type=c("TR_C_gene","TR_J_gene","TR_V_gene","TR_D_gene"))
    repeatMasker<-data.frame(gene_biotype="repeatMasker",feature_type=c("5S-Deu-L2","Alu","centr","CR1","DNA","DNA?","ERV1","ERV1?","ERVK","ERVL","ERVL?","ERVL-MaLR","Gypsy","Gypsy?","hAT","hAT?","hAT-Ac","hAT-Blackjack","hAT-Charlie","hAT-Tag1","hAT-Tip100","hAT-Tip100?","Helitron","Helitron?","L1","L2","Low_complexity","LTR","LTR?","MIR","MULE-MuDR","nonsense_mediated_decay","non_stop_decay","Penelope","PiggyBac","PiggyBac?","RNA","RTE-BovB","RTE-X","Satellite","Simple_repeat","SVA","TcMar?","TcMar-Mariner","TcMar-Tc2","TcMar-Tigger","telo","Unknown","acro","Crypton","Dong-R4","I-Jockey","Kolobok","L1-Tx1","Merlin","MULE-MuDR?","PIF-Harbinger","SINE?","TcMar","TcMar-Pogo","TcMar-Tc1"))
    intergenic<-data.frame(gene_biotype="intergenic",feature_type=c("intergenic"))
    circRNA<-data.frame(gene_biotype="circRNA",feature_type=c("circRNA"))
    category<-rbind(mRNA,pseudogene,lncRNA,sncRNA,rRNA,tRNA,IG_gene,TR_gene,repeatMasker,intergenic,circRNA)

    anno.biotype.bed<-anno.biotype.bed %>% left_join(category,by=c("type"="feature_type"))
    anno.biotype.bed<-anno.biotype.bed %>% mutate(type = ifelse(as.character(Transcript_id) == "intergenic", "intergenic", as.character(type)))
    input.bed<-read.delim(input.bed.file,sep="\t",header=F)
    probe<-which(input.bed$V4 %in% anno.biotype.bed$V4)
    input.bed<-input.bed[probe,]
    add_seq <- anno.biotype.bed %>% left_join(input.bed,by=c("V4"="V4"))


    extendSeq_split<-as.data.frame(str_split_fixed(add_seq[,"V7"], "Y", 2))
    colnames(extendSeq_split)<-c("extendSeq_bef","extendSeq_aft")
    add_seq<-cbind(add_seq,extendSeq_split)
    seq_group<-data.frame(uniq_seq=unique(add_seq$extendSeq_aft))
    seq_group$seq_id<-1:length(seq_group$uniq_seq)
    add_seq <- add_seq %>% left_join(seq_group,by=c("extendSeq_aft"="uniq_seq"))

    priority<-c("rRNA","tRNA","sncRNA","mRNA","lncRNA","circRNA","IG_gene","TR_gene","pseudogene","repeatMasker","intergenic")
    priority_table<-data.frame(type=priority,priority_rank=c(seq(length(priority))))
    add_seq <- add_seq %>% left_join(priority_table,by=c("gene_biotype"="type"))
    add_seq$group_id<-paste("group",add_seq$seq_id,add_seq$priority_rank,sep="-")
    write.table(add_seq,paste(outFile,"_add_seq_group.bed",sep=""),row.names=F,col.names=F,quote=F,sep="\t")
    write.xlsx(add_seq,paste(outFile,"_add_seq_group.xlsx",sep=""), overwrite = TRUE)


    add_seq_group_list<-split(add_seq,add_seq$seq_id)

    invisible(lapply(seq_along(add_seq_group_list),function(x){
        add_seq_group_list[[x]]<<-head(arrange(add_seq_group_list[[x]],priority_rank),1)
        }))

    add_seq_group_uniq<-do.call(rbind.data.frame, add_seq_group_list)
    dim(add_seq_group_uniq)


    colnames(add_seq_group_uniq)<-c("chrom",#1
    "chromStart",#2
    "chromEnd",#3
    "name",#4
    "logFC",#5
    "strand",#6
    "Transcript_id",#7
    "Transcript_name",#8
    "Gene_id",#9
    "Gene_name",#10
    "type",#11
    "Region",#12
    "gene_feature",#13
    "tolBaseNum",#14
    "tqueryCov",#15
    "tsampCov",#16
    "tupDist",#17
    "tdownDist",#18
    "gene_biotype",#19
    "chrom.y",#20
    "chromStart.y",#21
    "chromEnd.y",#22
    "logFC.y",#23
    "strand.y",#24
    "extendSeq",#25
    paste("group1_rep",1:length(group1filelist),sep=""),
    paste("group2_rep",1:length(group2filelist),sep=""),
    "AveExpr",
    "t",
    "P.Value",
    "adj.P.Val",
    "B",
    "group1_treatStopNumBaseMean",
    "group2_treatStopNumBaseMean",
    "group1_treatPreRpmFoldBaseMean",
    "group2_treatPreRpmFoldBaseMean",
    "group1_preRpmFoldRatioBaseMean",
    "group2_preRpmFoldRatioBaseMean",
    "group1_treatAfterRpmFoldBaseMean",
    "group2_treatAfterRpmFoldBaseMean",
    "group1_afterRpmFoldRatioBaseMean",
    "group2_afterRpmFoldRatioBaseMean",
    "group1_treatStopRatioBaseMean",
    "group2_treatStopRatioBaseMean",
    "group1_stopRatioFCBaseMean",
    "group2_stopRatioFCBaseMean",
    "extendSeq_bef",
    "extendSeq_aft",
    "seq_id",
    "priority_rank",
    "group_id")
    add_seq_group_uniq$uniq_id<-paste(add_seq_group_uniq$chrom,add_seq_group_uniq$chromStart,add_seq_group_uniq$chromEnd,add_seq_group_uniq$strand,sep="_")

    #volcano plot visualization
    add_seq_group_uniq<-add_seq_group_uniq %>% arrange(desc(logFC))

    add_seq_group_uniq$status<-as.factor(
      ifelse((add_seq_group_uniq$logFC > uplog2FC & add_seq_group_uniq$P.Value < pval),"Up_regulated",
             ifelse((add_seq_group_uniq$logFC < downlog2FC & add_seq_group_uniq$P.Value < pval),"Down_regulated",
                    ifelse(add_seq_group_uniq$P.Value > pval,"Not_significant","Not_expected"))))
    table(add_seq_group_uniq$status)
    # levels(add_seq_group_uniq$status)<-c("Down_regulated","Up_regulated","Not_significant","Not_expected")
    mycolors<-c(Down_regulated=brewer.pal(6,"PuBu")[6],Up_regulated=brewer.pal(3,"Set1")[1],Not_significant=brewer.pal(5,"Greys")[2],Not_expected=brewer.pal(5,"YlGn")[2])#blue/red/grey/green
    # barplot(c(2,5,7,9), col=mycolors)

    unique_gene_biotype<-unique(add_seq_group_uniq$gene_biotype)

    top3_Transcript<-c(head(which(add_seq_group_uniq$status=="Up_regulated"),3),tail(which(add_seq_group_uniq$status=="Down_regulated"),3))
    # add_seq_group_uniq$Gene_name[setdiff(seq_along(add_seq_group_uniq$Gene_name),top3_gene)]<-""
    add_seq_group_uniq$Transcript_name_top<-add_seq_group_uniq$Transcript_name
    add_seq_group_uniq$Transcript_name_top[setdiff(seq_along(add_seq_group_uniq$Transcript_name_top),top3_Transcript)]<-""

    if(dim(add_seq_group_uniq)[1] != 0){
      # Set it globally:
      options(ggrepel.max.overlaps = Inf)
      p_volcano <- ggplot(add_seq_group_uniq,
                  aes(x=logFC,
                      y=-log10(P.Value),
                      colour=str_wrap(status,width = 2),
                      label =Transcript_name_top) ) +
        geom_jitter(size=2,alpha=0.6)
      p_volcano <- p_volcano + scale_color_manual(values = mycolors)
      p_volcano <- p_volcano + geom_vline(xintercept=c(downlog2FC,uplog2FC), linetype="longdash", size=0.2)
      p_volcano <- p_volcano + geom_hline(yintercept=c(-log10(pval)), linetype="longdash", size=0.2)
      p_volcano<-p_volcano+
        # labs(title = "Group2 - Group1")+
        labs(title = group_name)+
        xlab("log2FC")+
        theme_bw()+
        theme(legend.position = "top",
              plot.title=element_text(hjust=0.5),
              legend.title =element_blank(),
              legend.key.width = unit(1.5, "cm"),
              legend.key.size = unit(0, "cm"),
              legend.spacing.x = unit(0,"cm"),
              legend.text.align = 0.7
              # plot.margin=unit(c(1,1,1,1),units="cm")
              )+
        scale_x_continuous(limits = c(-max(abs(min(add_seq_group_uniq$logFC)),max(add_seq_group_uniq$logFC)), max(abs(min(add_seq_group_uniq$logFC)),max(add_seq_group_uniq$logFC))))+
        guides(color=guide_legend(nrow=2,ncol=2,byrow=TRUE))

      p_volcano_text<-p_volcano+geom_text_repel(color="black",
                                family="sans",
                                force = 15)


      pdf(paste(outFile,"_DGE_volcano_plot.pdf",sep=""),width=8,height=8)
      print(p_volcano_text)
      dev.off()
    }

    write.xlsx(add_seq_group_uniq,paste(outFile,"_add_seq_group_uniq.xlsx",sep=""), overwrite = TRUE)
    add_seq_group_uniq_filt<-add_seq_group_uniq[add_seq_group_uniq$status=="Up_regulated" | add_seq_group_uniq$status=="Down_regulated",]

    if(dim(add_seq_group_uniq_filt)[1] != 0){

      print("pie plot for gene_biotype")
      pseudoU_anno_genetype_num<-arrange(as.data.frame(table(add_seq_group_uniq_filt$gene_biotype)),desc(Freq))
      write.table(pseudoU_anno_genetype_num,paste(outFile,"_anno_gene_biotype_num.sort",sep=""),row.names=F,col.names=F,quote=F,sep="\t")
      percent<-round(100*pseudoU_anno_genetype_num$Freq/sum(pseudoU_anno_genetype_num$Freq),2)
      percent <-paste('(',percent, "%", ", ", pseudoU_anno_genetype_num$Freq,')', sep = "")
      pseudoU_anno_genetype_num$Var1 <- paste(pseudoU_anno_genetype_num$Var1, percent, sep = '')
      pseudoU_anno_genetype_num$Var1 <- factor(pseudoU_anno_genetype_num$Var1, levels = pseudoU_anno_genetype_num$Var1)
      mycol = c(brewer.pal(12, "Set3"),brewer.pal(12, "Paired"),brewer.pal(11, "Spectral"),brewer.pal(9, "Set1"),brewer.pal(8, "Set2"),brewer.pal(8, "Pastel1"),brewer.pal(8, "Pastel2"),brewer.pal(8, "Dark2"),brewer.pal(8, "Accent"))
      pie_plot<-ggplot(data = pseudoU_anno_genetype_num, mapping = aes(x = 'Content', y = Freq, fill = Var1)) +
          geom_bar(stat = 'identity', position = 'stack', width = 1) + coord_polar(theta = 'y') + labs(x = '', y = '', title = '') +
          theme(axis.text = element_blank()) +
            theme(axis.ticks = element_blank()) +
            theme(panel.grid = element_blank(),
            panel.background=element_blank(),
            axis.text = element_blank(),
            axis.title = element_blank(),
            axis.ticks = element_blank()) +
            guides(fill = guide_legend(title = "gene biotype"))+
            scale_fill_manual(values = mycol)

      pdf(paste(outFile,"_gene_biotype_piechart.pdf",sep=""))
      print(pie_plot)
      invisible(dev.off())

      pseudoU_anno_genetype_num<-arrange(as.data.frame(table(add_seq_group_uniq_filt$gene_feature)),desc(Freq))
      percent<-round(100*pseudoU_anno_genetype_num$Freq/sum(pseudoU_anno_genetype_num$Freq),2)
      percent <-paste('(',percent, "%", ", ", pseudoU_anno_genetype_num$Freq,')', sep = "")
      pseudoU_anno_genetype_num$Var1 <- paste(pseudoU_anno_genetype_num$Var1, percent, sep = '')
      pseudoU_anno_genetype_num$Var1 <- factor(pseudoU_anno_genetype_num$Var1, levels = pseudoU_anno_genetype_num$Var1)
      mycol = c(brewer.pal(12, "Set3"),brewer.pal(12, "Paired"),brewer.pal(11, "Spectral"),brewer.pal(9, "Set1"),brewer.pal(8, "Set2"),brewer.pal(8, "Pastel1"),brewer.pal(8, "Pastel2"),brewer.pal(8, "Dark2"),brewer.pal(8, "Accent"))
      pie_plot<-ggplot(data = pseudoU_anno_genetype_num, mapping = aes(x = 'Content', y = Freq, fill = Var1)) +
          geom_bar(stat = 'identity', position = 'stack', width = 1) + coord_polar(theta = 'y') + labs(x = '', y = '', title = '') +
          theme(axis.text = element_blank()) +
            theme(axis.ticks = element_blank()) +
            theme(panel.grid = element_blank(),
            panel.background=element_blank(),
            axis.text = element_blank(),
            axis.title = element_blank(),
            axis.ticks = element_blank()) +
            guides(fill = guide_legend(title = "gene feature"))+
            scale_fill_manual(values = mycol)

      pdf(paste(outFile,"_gene_feature_piechart.pdf",sep=""))
      print(pie_plot)
      invisible(dev.off())

      write.table(add_seq_group_uniq_filt,paste(outFile,"_add_seq_group_uniq_filt.bed",sep=""),row.names=F,col.names=F,quote=F,sep="\t")
      write.xlsx(add_seq_group_uniq_filt,paste(outFile,"_add_seq_group_uniq_filt.xlsx",sep=""), overwrite = TRUE)
      write.xlsx(add_seq_group_uniq_filt,paste(outFile,"_anno_group_redundance.xlsx",sep=""), overwrite = TRUE)

      add_seq_group_uniq_filt_high_confid<-add_seq_group_uniq_filt %>% filter((group1_treatPreRpmFoldBaseMean >4 | group2_treatPreRpmFoldBaseMean >4) & (group1_treatAfterRpmFoldBaseMean >4 | group2_treatAfterRpmFoldBaseMean >4) & (group1_stopRatioFCBaseMean >4 | group2_stopRatioFCBaseMean >4))
      write.xlsx(add_seq_group_uniq_filt_high_confid,paste(outFile,"_anno_group_redundance_high_confid.xlsx",sep=""), overwrite = TRUE)


      plot_table<-as.data.frame(table(add_seq_group_uniq_filt$gene_feature))
      colnames(plot_table)<-c("Gene Feature","Num")
      to_plot<-arrange(plot_table,desc(Num))
      tt3 <- ttheme_minimal(
      core=list(bg_params = list(fill = blues9[1:4], col=NA),
                fg_params=list(fontface=3)),
      colhead=list(fg_params=list(col="navyblue", fontface=4L)),
      rowhead=list(fg_params=list(col="orange", fontface=3L)))

      pdf(paste(outFile, '_plot_table.pdf', sep=""), width = 7, height = 7) # Open a new pdf file
      grid.arrange(
      tableGrob(to_plot, theme=tt3),
      nrow=1)
      invisible(dev.off()) # Close the file

      add_seq_group_uniq_filt_rRNA<-add_seq_group_uniq_filt[add_seq_group_uniq_filt$gene_biotype%in%"rRNA",]
      if(dim(add_seq_group_uniq_filt_rRNA)[1] != 0){
        print("pie plot for RMBase evidence")
        #read human.hg38.Pseudo.result.col29.xlsx
        cat("Detecting novel psi (human.hg38.Pseudo.result.col29.xlsx: all known pseudouridylation sites download from RMBase)\n")
        human.hg38.Pseudo.result.col29.xlsx<-read.xlsx(human.hg38.Pseudo.result.col29.xlsx.file)#"human.hg38.Pseudo.result.col29.xlsx"
        human.hg38.Pseudo.result.col29.xlsx.tab<-as.data.frame(table(add_seq_group_uniq_filt$uniq_id%in%human.hg38.Pseudo.result.col29.xlsx$uniq_id))
        colnames(human.hg38.Pseudo.result.col29.xlsx.tab)<-c("known","hit")
        print(human.hg38.Pseudo.result.col29.xlsx.tab, row.names = F)
        novel_psi<-add_seq_group_uniq_filt[add_seq_group_uniq_filt$uniq_id%in%human.hg38.Pseudo.result.col29.xlsx$uniq_id=="FALSE",]
        common_psi<-add_seq_group_uniq_filt[add_seq_group_uniq_filt$uniq_id%in%human.hg38.Pseudo.result.col29.xlsx$uniq_id=="TRUE",]
        common_psi<-common_psi %>% left_join(human.hg38.Pseudo.result.col29.xlsx,by=c("uniq_id"="uniq_id"))
        write.xlsx(novel_psi,paste(outFile,"human.hg38.Pseudo.result.col29_novel.xlsx",sep=""), overwrite = TRUE)
        write.xlsx(common_psi,paste(outFile,"human.hg38.Pseudo.result.col29_common.xlsx",sep=""), overwrite = TRUE)

        known_num<-arrange(as.data.frame(table(add_seq_group_uniq_filt$uniq_id%in%human.hg38.Pseudo.result.col29.xlsx$uniq_id)),desc(Freq))
        known_num$Var1<-str_replace(known_num$Var1,"TRUE","Known")
        known_num$Var1<-str_replace(known_num$Var1,"FALSE","Novel")
        percent<-round(100*known_num$Freq/sum(known_num$Freq),2)
        percent <-paste('(',percent, "%", ", ", known_num$Freq,')', sep = "")
        known_num$Var1 <- paste(known_num$Var1, percent, sep = '')
        known_num$Var1 <- factor(known_num$Var1, levels = known_num$Var1)
        mycol = c(brewer.pal(12, "Set3"),brewer.pal(12, "Paired"),brewer.pal(11, "Spectral"))
        pie_plot<-ggplot(data = known_num, mapping = aes(x = 'Content', y = Freq, fill = Var1)) +
            geom_bar(stat = 'identity', position = 'stack', width = 1) + coord_polar(theta = 'y') + labs(x = '', y = '', title = '') +
            theme(axis.text = element_blank()) +
              theme(axis.ticks = element_blank()) +
              theme(panel.grid = element_blank(),
              panel.background=element_blank(),
              axis.text = element_blank(),
              axis.title = element_blank(),
              axis.ticks = element_blank()) +
              guides(fill = guide_legend(title = "Total Known/Novel psi"))+
              scale_fill_manual(values = mycol)

        pdf(paste(outFile,"_known_novel_psi_piechart.pdf",sep=""))
        print(pie_plot)
        invisible(dev.off())


        print("pie plot for known rRNA evidence")
        #read hg38.psiU.SingleSites.bed
        add_seq_group_uniq_filt_rRNA<-add_seq_group_uniq_filt [add_seq_group_uniq_filt$gene_biotype=="rRNA" ,]
        cat("Detecting novel rRNA psi (hg38.psiU.SingleSites.bed: all known pseudouridylation sites in rRNA)\n")
        hg38.psiU.SingleSites.bed<-read.table(hg38.psiU.SingleSites.bed.file,head=F)#"hg38.psiU.SingleSites.bed"
        colnames(hg38.psiU.SingleSites.bed)<-c("chrom","chromStart","chromEnd","rRNA_anno","score","strand")
        hg38.psiU.SingleSites.bed$uniq_id<-paste(hg38.psiU.SingleSites.bed$chrom,hg38.psiU.SingleSites.bed$chromStart,hg38.psiU.SingleSites.bed$chromEnd,hg38.psiU.SingleSites.bed$strand,sep="_")
        hg38.psiU.SingleSites.bed.tab<-as.data.frame(table(add_seq_group_uniq_filt_rRNA$uniq_id%in%hg38.psiU.SingleSites.bed$uniq_id))
        colnames(hg38.psiU.SingleSites.bed.tab)<-c("known","hit")
        print(hg38.psiU.SingleSites.bed.tab, row.names = F)
        novel_rRNA_psi<-add_seq_group_uniq_filt_rRNA[add_seq_group_uniq_filt_rRNA$uniq_id%in%hg38.psiU.SingleSites.bed$uniq_id=="FALSE",]
        common_rRNA_psi<-add_seq_group_uniq_filt_rRNA[add_seq_group_uniq_filt_rRNA$uniq_id%in%hg38.psiU.SingleSites.bed$uniq_id=="TRUE",]
        common_rRNA_psi<-common_rRNA_psi %>% left_join(hg38.psiU.SingleSites.bed,by=c("uniq_id"="uniq_id"))
        write.xlsx(novel_rRNA_psi,paste(outFile,"hg38.psiU.SingleSites_novel.xlsx",sep=""), overwrite = TRUE)
        write.xlsx(common_rRNA_psi,paste(outFile,"hg38.psiU.SingleSites_common.xlsx",sep=""), overwrite = TRUE)

        known_rRNA_num<-arrange(as.data.frame(table(add_seq_group_uniq_filt_rRNA$uniq_id%in%hg38.psiU.SingleSites.bed$uniq_id)),desc(Freq))
        known_rRNA_num$Var1<-str_replace(known_rRNA_num$Var1,"TRUE","Known")
        known_rRNA_num$Var1<-str_replace(known_rRNA_num$Var1,"FALSE","Novel")
        percent<-round(100*known_rRNA_num$Freq/sum(known_rRNA_num$Freq),2)
        percent <-paste('(',percent, "%", ", ", known_rRNA_num$Freq,')', sep = "")
        known_rRNA_num$Var1 <- paste(known_rRNA_num$Var1, percent, sep = '')
        known_rRNA_num$Var1 <- factor(known_rRNA_num$Var1, levels = known_rRNA_num$Var1)
        mycol = c(brewer.pal(12, "Set3"),brewer.pal(12, "Paired"),brewer.pal(11, "Spectral"))
        pie_plot<-ggplot(data = known_rRNA_num, mapping = aes(x = 'Content', y = Freq, fill = Var1)) +
            geom_bar(stat = 'identity', position = 'stack', width = 1) + coord_polar(theta = 'y') + labs(x = '', y = '', title = '') +
            theme(axis.text = element_blank()) +
              theme(axis.ticks = element_blank()) +
              theme(panel.grid = element_blank(),
              panel.background=element_blank(),
              axis.text = element_blank(),
              axis.title = element_blank(),
              axis.ticks = element_blank()) +
              guides(fill = guide_legend(title = "rRNA Known/Novel psi"))+
              scale_fill_manual(values = mycol)

        pdf(paste(outFile,"_known_novel_rRNA_psi_piechart.pdf",sep=""))
        print(pie_plot)
        invisible(dev.off())
        }else{
          print("No known/novel rRNA detected!")
        }

    }else{
      print("No differentially expressed psi-sites detected!")
    }


Optional Parameters
----------------------

up-log2FC
********************
Users can customize the up-regulated fold change threshold by setting ``up-log2FC``, for example, if ``up-log2FC`` is set as ``1.00``, then result with log2(fold change) > 1.00 is retained.

down-log2FC
**********************
Users can customize the down-regulated fold change threshold by setting ``down-log2FC``, for example, if ``down-log2FC`` is set as ``-1.00``, then result with log2(fold change) < -1.00 is retained.


P.Value
********************
Users can customize the statistical test threshold by setting ``P.Value``, for example, if ``P.Value`` is set as ``0.05``, then result with p value < 0.05 is retained.

min-fold
**********************
Users can customize the minimum stopRatioFC threshold of each group by setting ``min-fold``, for example, if ``min-fold`` is set as ``2.00``, then result with minimum ``stopRatioFC > 2.00`` (each sample) is retained.

Output
--------

Information
************

Result with ``_group_redundance_high_confid.xlsx`` suffix is the final differential modification result, which is a filter result of ``_group_redundance.xlsx`` (filter by ``Optional Parameters``).

.. code:: bash

    $ cd /the/directory/of/out_file_dir
    $ tree -L 1
    .
    ├── B11_B12_svm_all.bed
    ├── B1_B2_svm_all.bed
    ├── B3_B4_svm_all.bed
    ├── B5_B6_svm_all.bed
    ├── B7_B8_svm_all.bed
    ├── B9_B10_svm_all.bed
    ├── polyaRNA_Day4-polyaRNA_Day0_limma_bed_annotation.log
    ├── polyaRNA_Day4-polyaRNA_Day0_limma.log
    ├── polyaRNA_Day4-polyaRNA_Day0_treat_control_overall_filt_add_seq_group.bed
    ├── polyaRNA_Day4-polyaRNA_Day0_treat_control_overall_filt_add_seq_group_uniq_filt.bed
    ├── polyaRNA_Day4-polyaRNA_Day0_treat_control_overall_filt_add_seq_group_uniq_filt.xlsx
    ├── polyaRNA_Day4-polyaRNA_Day0_treat_control_overall_filt_add_seq_group_uniq.xlsx
    ├── polyaRNA_Day4-polyaRNA_Day0_treat_control_overall_filt_add_seq_group.xlsx
    ├── polyaRNA_Day4-polyaRNA_Day0_treat_control_overall_filt_anno_append.bed
    ├── polyaRNA_Day4-polyaRNA_Day0_treat_control_overall_filt_anno.bed
    ├── polyaRNA_Day4-polyaRNA_Day0_treat_control_overall_filt_anno_gene_biotype_num.sort
    ├── polyaRNA_Day4-polyaRNA_Day0_treat_control_overall_filt_anno_group_redundance_high_confid.xlsx
    ├── polyaRNA_Day4-polyaRNA_Day0_treat_control_overall_filt_anno_group_redundance.xlsx
    ├── polyaRNA_Day4-polyaRNA_Day0_treat_control_overall_filt.bed
    ├── polyaRNA_Day4-polyaRNA_Day0_treat_control_overall_filt.bed6
    ├── polyaRNA_Day4-polyaRNA_Day0_treat_control_overall_filt_DGE_volcano_plot.pdf
    ├── polyaRNA_Day4-polyaRNA_Day0_treat_control_overall_filt_gene_biotype_piechart.pdf
    ├── polyaRNA_Day4-polyaRNA_Day0_treat_control_overall_filt_gene_feature_piechart.pdf
    ├── polyaRNA_Day4-polyaRNA_Day0_treat_control_overall_filthg38.psiU.SingleSites_common.xlsx
    ├── polyaRNA_Day4-polyaRNA_Day0_treat_control_overall_filthg38.psiU.SingleSites_novel.xlsx
    ├── polyaRNA_Day4-polyaRNA_Day0_treat_control_overall_filthuman.hg38.Pseudo.result.col29_common.xlsx
    ├── polyaRNA_Day4-polyaRNA_Day0_treat_control_overall_filthuman.hg38.Pseudo.result.col29_novel.xlsx
    ├── polyaRNA_Day4-polyaRNA_Day0_treat_control_overall_filt_known_novel_psi_piechart.pdf
    ├── polyaRNA_Day4-polyaRNA_Day0_treat_control_overall_filt_known_novel_rRNA_psi_piechart.pdf
    ├── polyaRNA_Day4-polyaRNA_Day0_treat_control_overall_filt_plot_table.pdf
    ├── polyaRNA_Day4-polyaRNA_Day0_treat_control_overall_filt.xlsx
    └── polyaRNA_Day4-polyaRNA_Day0_treat_control_overall.xlsx

    0 directories, 32 files

Diagram
********
File with suffix ``_volcano_plot.pdf`` is a volcano plot of Ψ-sites differential modification on input Ψ-sites file (rtsSeeker result).

.. image:: /images/Differential_Modification_volcano_plot.png


.. note:: All user input will be recorded in a plain text file with suffix ``_pseUdiff_config.txt`` in psiFinder/config and help users to easily reload the previous config (by simply clicking ``CONFIG`` button).
