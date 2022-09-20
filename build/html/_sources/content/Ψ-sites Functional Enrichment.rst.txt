Ψ-sites Functional Enrichment
===========================================

.. contents::
    :local:

psiFinder ``Ψ-sites Functional Enrichment`` utilize `gprofiler2 <https://cran.r-project.org/web/packages/gprofiler2/vignettes/gprofiler2.html>`_ and `pathview <https://bioconductor.org/packages/release/bioc/vignettes/pathview/inst/doc/pathview.pdf>`_ to generate Gene ontology (GO) result for Ψ-sites functional enrichment based on ``stopRatioFC`` of input rtsSeeker result files.

.. image:: /images/Functional_Enrichment_1.png

If ``ADD`` button is clicked, then a new input text edit field will be added. Users can input multiple gene list and each gene list result will be outputed to the same target directory, for example:

.. code:: bash

    $ cd /the/directory/of/out_file_dir
    $ tree -L 1
    .
    ├── Day0_com_gprofiler.csv
    ├── Day0_com_gprofiler.pdf
    ├── Day0_com_top10_REAC_bubble.pdf
    ├── Day0_com_top12GO_CC_BP_MF_bubble.pdf
    ├── Day0_com.txt
    ├── Day4_com_Day0_com.pdf
    ├── Day4_com_gprofiler.csv
    ├── Day4_com_gprofiler.pdf
    ├── Day4_com_top10_REAC_bubble.pdf
    ├── Day4_com_top12GO_CC_BP_MF_bubble.pdf
    └── Day4_com.txt

    0 directories, 11 files

.. image:: /images/Functional_Enrichment_2.png

``pseUfun`` also support KEGG pathway enrichment, if column of gene value is added, then KEGG enrichment result will be outputed as well.

.. code:: bash

    $ cd /the/directory/of/out_file_dir
    $ tree -L 1
    .
    ├── Day4_com_gprofiler.csv
    ├── Day4_com_gprofiler.pdf
    ├── Day4_com_pathview.pdf
    ├── Day4_com.pdf
    ├── Day4_com_top10_REAC_bubble.pdf
    ├── Day4_com_top12GO_CC_BP_MF_bubble.pdf
    ├── Day4_com_topKEGG.csv
    ├── Day4_com.txt
    ├── hsa03010.png
    ├── hsa03010.Ribosome.png
    ├── hsa03010.xml
    ├── hsa05012.Parkinson_disease.png
    ├── hsa05012.png
    ├── hsa05012.xml
    ├── hsa05171.Coronavirus_disease_-_COVID-19.png
    ├── hsa05171.png
    └── hsa05171.xml

    0 directories, 17 files

.. image:: /images/Functional_Enrichment_3.png

Input
------
Users should input gene IDs (proteins, transcripts, microarray IDs, etc), SNP IDs, chromosomal intervals or term IDs id to ``pseUfun`` QT widget, to construct Ψ-sites functional enrichment.

Construct Ψ-sites functional enrichment
---------------------------------------------
Once click ``START``, psiFindeer will run ``pseUfun.sh`` and ``pseUfun.r``.

``pseUfun.sh``

.. code:: bash

    #!/bin/bash
    Help() {
    echo
    "pseUfun_wrap.sh example:
     bash pseUfun_wrap.sh -c cmd -o outfile_path
    "
    }


    usage() {                                      # Function: Print a help message.
      echo "Usage: $0 [ -h help] [ -c cmd ] [ -l filelist_join2 ] [ -o outfile_path ]" 1>&2
    }
    exit_abnormal() {                              # Function: Exit with error.
      usage
      exit 1
    }


    while getopts ":h:c:l:o:" options; do

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
        l)
          filelist_join2=${OPTARG}
          if ! [[ -n $filelist_join2 ]] ; then
            echo "You didn't input the filelist_join2"
          fi
          ;;
        o)
          outfile_path=${OPTARG}
          if ! [[ -n $outfile_path ]] ; then
            echo "You didn't input the outfile_path"
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
    echo "Starting: one pseUfun diagram job is starting..."

    eval $cmd
    echo 'pseUfun diagram plotting is done!'

    if [ ! -f "${outfile_path}/${filelist_join2}.pdf"  ]
    then
    echo -e "${outfile_path}/${filelist_join2}.pdf not exist"
    else
    echo -e "browse ${outfile_path}/${filelist_join2}.pdf"
    mupdf-x11 ${outfile_path}/${filelist_join2}.pdf &> /dev/null
    fi

    if [ ! -f "${outfile_path}/${filelist_join2}.txt"  ]
    then
    echo -e "${outfile_path}/${filelist_join2}.pdf not exist"
    else
    echo -e "browse pdf in ${outfile_path}/${filelist_join2}.txt"
    txt_var=`awk 'BEGIN{ORS="";}{ for (i=NF; i>0; i--) print $i," "; }' ${outfile_path}/${filelist_join2}.txt`
    # echo $txt_var
    mutool merge -o ${outfile_path}/${filelist_join2}_pathview.pdf $txt_var &> /dev/null
    mupdf-x11 ${outfile_path}/${filelist_join2}_pathview.pdf &> /dev/null
    fi

    exit 0  # Exit normally.

``pseUfun.r``

.. code:: R

    suppressMessages(library("gprofiler2"))
    suppressMessages(library("gridExtra"))
    suppressMessages(library("RColorBrewer"))
    suppressMessages(library("ggplot2"))
    suppressMessages(library("ggrepel"))
    suppressMessages(library("dplyr"))
    suppressMessages(library("optparse"))
    suppressMessages(library("stringr"))
    suppressMessages(library("pathview"))

    option_list = list(
      make_option(c("-f","--filelist"),default="geneset1.txt",
                   help="comma separated list of files (default %default)"),
      make_option(c("-s","--species"),default="hsapiens",
                   help="Species choice (default %default)")
    );
    opt_parser = OptionParser(option_list=option_list);
    opt = parse_args(opt_parser);

    if (is.null(opt$filelist) || is.null(opt$species) ){
      print_help(opt_parser);
      stop("Please provide -f group1_name, -s species option", call.=FALSE);
    }

    filelist <- as.list(unlist(strsplit(opt$filelist, ",")))


    ####pathview2####
    pathview2<-function (gene.data = NULL, cpd.data = NULL, pathway.id, species = "hsa",
                         kegg.dir = ".", cpd.idtype = "kegg", gene.idtype = "entrez",
                         gene.annotpkg = NULL, min.nnodes = 3, kegg.native = TRUE,
                         map.null = TRUE, expand.node = FALSE, split.group = FALSE,
                         map.symbol = TRUE, map.cpdname = TRUE, node.sum = "sum",
                         discrete = list(gene = FALSE, cpd = FALSE), limit = list(gene = 1,
                                                                                  cpd = 1), bins = list(gene = 10, cpd = 10), both.dirs = list(gene = T,
                                                                                                                                               cpd = T), trans.fun = list(gene = NULL, cpd = NULL),
                         low = list(gene = "green", cpd = "blue"), mid = list(gene = "gray",
                                                                              cpd = "gray"), high = list(gene = "red", cpd = "yellow"),
                         na.col = "transparent", ...)
    {
      dtypes = !is.null(gene.data) + !is.null(cpd.data)
      cond0 = dtypes == 1 & is.numeric(limit) & length(limit) >
        1
      if (cond0) {
        if (limit[1] != limit[2] & is.null(names(limit)))
          limit = list(gene = limit[1:2], cpd = limit[1:2])
      }
      if (is.null(trans.fun))
        trans.fun = list(gene = NULL, cpd = NULL)
      arg.len2 = c("discrete", "limit", "bins", "both.dirs", "trans.fun",
                   "low", "mid", "high")
      for (arg in arg.len2) {
        obj1 = eval(as.name(arg))
        if (length(obj1) == 1)
          obj1 = rep(obj1, 2)
        if (length(obj1) > 2)
          obj1 = obj1[1:2]
        obj1 = as.list(obj1)
        ns = names(obj1)
        if (length(ns) == 0 | !all(c("gene", "cpd") %in% ns))
          names(obj1) = c("gene", "cpd")
        assign(arg, obj1)
      }
      if (is.character(gene.data)) {
        gd.names = gene.data
        gene.data = rep(1, length(gene.data))
        names(gene.data) = gd.names
        both.dirs$gene = FALSE
        ng = length(gene.data)
        nsamp.g = 1
      }
      else if (!is.null(gene.data)) {
        if (length(dim(gene.data)) == 2) {
          gd.names = rownames(gene.data)
          ng = nrow(gene.data)
          nsamp.g = 2
        }
        else if (is.numeric(gene.data) & is.null(dim(gene.data))) {
          gd.names = names(gene.data)
          ng = length(gene.data)
          nsamp.g = 1
        }
        else stop("wrong gene.data format!")
      }
      else if (is.null(cpd.data)) {
        stop("gene.data and cpd.data are both NULL!")
      }
      gene.idtype = toupper(gene.idtype)
      data(bods)
      if (species != "ko") {
        species.data = kegg.species.code(species, na.rm = T,
                                         code.only = FALSE)
      }
      else {
        species.data = c(kegg.code = "ko", entrez.gnodes = "0",
                         kegg.geneid = "K01488", ncbi.geneid = NA, ncbi.proteinid = NA,
                         uniprot = NA)
        gene.idtype = "KEGG"
        msg.fmt = "Only KEGG ortholog gene ID is supported, make sure it looks like \"%s\"!"
        msg = sprintf(msg.fmt, species.data["kegg.geneid"])
        message("Note: ", msg)
      }
      if (length(dim(species.data)) == 2) {
        message("Note: ", "More than two valide species!")
        species.data = species.data[1, ]
      }
      species = species.data["kegg.code"]
      entrez.gnodes = species.data["entrez.gnodes"] == 1
      if (is.na(species.data["ncbi.geneid"])) {
        if (!is.na(species.data["kegg.geneid"])) {
          msg.fmt = "Mapping via KEGG gene ID (not Entrez) is supported for this species,\nit looks like \"%s\"!"
          msg = sprintf(msg.fmt, species.data["kegg.geneid"])
          message("Note: ", msg)
        }
        else {
          stop("This species is not annotated in KEGG!")
        }
      }
      if (is.null(gene.annotpkg))
        gene.annotpkg = bods[match(species, bods[, 3]), 1]
      if (length(grep("ENTREZ|KEGG|NCBIPROT|UNIPROT", gene.idtype)) <
          1 & !is.null(gene.data)) {
        if (is.na(gene.annotpkg))
          stop("No proper gene annotation package available!")
        if (!gene.idtype %in% gene.idtype.bods[[species]])
          stop("Wrong input gene ID type!")
        gene.idmap = id2eg(gd.names, category = gene.idtype,
                           pkg.name = gene.annotpkg, unique.map = F)
        gene.data = mol.sum(gene.data, gene.idmap)
        gene.idtype = "ENTREZ"
      }
      if (gene.idtype != "KEGG" & !entrez.gnodes & !is.null(gene.data)) {
        id.type = gene.idtype
        if (id.type == "ENTREZ")
          id.type = "ENTREZID"
        kid.map = names(species.data)[-c(1:2)]
        kid.types = names(kid.map) = c("KEGG", "ENTREZID", "NCBIPROT",
                                       "UNIPROT")
        kid.map2 = gsub("[.]", "-", kid.map)
        kid.map2["UNIPROT"] = "up"
        if (is.na(kid.map[id.type]))
          stop("Wrong input gene ID type for the species!")
        message("Info: Getting gene ID data from KEGG...")
        gene.idmap = keggConv(kid.map2[id.type], species)
        message("Info: Done with data retrieval!")
        kegg.ids = gsub(paste(species, ":", sep = ""), "", names(gene.idmap))
        in.ids = gsub(paste0(kid.map2[id.type], ":"), "", gene.idmap)
        gene.idmap = cbind(in.ids, kegg.ids)
        gene.data = mol.sum(gene.data, gene.idmap)
        gene.idtype = "KEGG"
      }
      if (is.character(cpd.data)) {
        cpdd.names = cpd.data
        cpd.data = rep(1, length(cpd.data))
        names(cpd.data) = cpdd.names
        both.dirs$cpd = FALSE
        ncpd = length(cpd.data)
      }
      else if (!is.null(cpd.data)) {
        if (length(dim(cpd.data)) == 2) {
          cpdd.names = rownames(cpd.data)
          ncpd = nrow(cpd.data)
        }
        else if (is.numeric(cpd.data) & is.null(dim(cpd.data))) {
          cpdd.names = names(cpd.data)
          ncpd = length(cpd.data)
        }
        else stop("wrong cpd.data format!")
      }
      if (length(grep("kegg", cpd.idtype)) < 1 & !is.null(cpd.data)) {
        data(rn.list)
        cpd.types = c(names(rn.list), "name")
        cpd.types = tolower(cpd.types)
        cpd.types = cpd.types[-grep("kegg", cpd.types)]
        if (!tolower(cpd.idtype) %in% cpd.types)
          stop("Wrong input cpd ID type!")
        cpd.idmap = cpd2kegg(cpdd.names, in.type = cpd.idtype)
        cpd.data = mol.sum(cpd.data, cpd.idmap)
      }
      warn.fmt = "Parsing %s file failed, please check the file!"
      if (length(grep(species, pathway.id)) > 0) {
        pathway.name = pathway.id
        pathway.id = gsub(species, "", pathway.id)
      }
      else pathway.name = paste(species, pathway.id, sep = "")
      kfiles = list.files(path = kegg.dir, pattern = "[.]xml|[.]png")
      npath = length(pathway.id)
      out.list = list()
      tfiles.xml = paste(pathway.name, "xml", sep = ".")
      tfiles.png = paste(pathway.name, "png", sep = ".")
      if (kegg.native)
        ttype = c("xml", "png")
      else ttype = "xml"
      xml.file <- paste(kegg.dir, "/", tfiles.xml, sep = "")
      for (i in 1:npath) {
        if (kegg.native)
          tfiles = c(tfiles.xml[i], tfiles.png[i])
        else tfiles = tfiles.xml[i]
        if (!all(tfiles %in% kfiles)) {
          dstatus = download.kegg(pathway.id = pathway.id[i],
                                  species = species, kegg.dir = kegg.dir, file.type = ttype)
          if (dstatus == "failed") {
            warn.fmt = "Failed to download KEGG xml/png files, %s skipped!"
            warn.msg = sprintf(warn.fmt, pathway.name[i])
            message("Warning: ", warn.msg)
            return(invisible(0))
          }
        }
        if (kegg.native) {
          node.data = try(node.info(xml.file[i]), silent = T)
          if (class(node.data) == "try-error") {
            warn.msg = sprintf(warn.fmt, xml.file[i])
            message("Warning: ", warn.msg)
            return(invisible(0))
          }
          node.type = c("gene", "enzyme", "compound", "ortholog")
          sel.idx = node.data$type %in% node.type
          nna.idx = !is.na(node.data$x + node.data$y + node.data$width +
                             node.data$height)
          sel.idx = sel.idx & nna.idx
          if (sum(sel.idx) < min.nnodes) {
            warn.fmt = "Number of mappable nodes is below %d, %s skipped!"
            warn.msg = sprintf(warn.fmt, min.nnodes, pathway.name[i])
            message("Warning: ", warn.msg)
            return(invisible(0))
          }
          node.data = lapply(node.data, "[", sel.idx)
        }
        else {
          gR1 = try(parseKGML2Graph2(xml.file[i], genes = F,
                                     expand = expand.node, split.group = split.group),
                    silent = T)
          node.data = try(node.info(gR1), silent = T)
          if (class(node.data) == "try-error") {
            warn.msg = sprintf(warn.fmt, xml.file[i])
            message("Warning: ", warn.msg)
            return(invisible(0))
          }
        }
        if (species == "ko")
          gene.node.type = "ortholog"
        else gene.node.type = "gene"
        if ((!is.null(gene.data) | map.null) & sum(node.data$type ==
                                                   gene.node.type) > 1) {
          plot.data.gene = node.map(gene.data, node.data, node.types = gene.node.type,
                                    node.sum = node.sum, entrez.gnodes = entrez.gnodes)
          plot.data.gene<-plot.data.gene[rowSums(plot.data.gene[,c("x","y","width","height")])!=4,]
          kng = plot.data.gene$kegg.names
          kng.char = gsub("[0-9]", "", unlist(kng))
          if (any(kng.char > ""))
            entrez.gnodes = FALSE
          if (map.symbol & species != "ko" & entrez.gnodes) {
            if (is.na(gene.annotpkg)) {
              warn.fmt = "No annotation package for the species %s, gene symbols not mapped!"
              warn.msg = sprintf(warn.fmt, species)
              message("Warning: ", warn.msg)
            }
            else {
              plot.data.gene$labels = eg2id(as.character(plot.data.gene$kegg.names),
                                            category = "SYMBOL", pkg.name = gene.annotpkg)[,
                                                                                           2]
              mapped.gnodes = rownames(plot.data.gene)
              node.data$labels[mapped.gnodes] = plot.data.gene$labels
            }
          }
          cols.ts.gene = node.color(plot.data.gene, limit$gene,
                                    bins$gene, both.dirs = both.dirs$gene, trans.fun = trans.fun$gene,
                                    discrete = discrete$gene, low = low$gene, mid = mid$gene,
                                    high = high$gene, na.col = na.col)
        }
        else plot.data.gene = cols.ts.gene = NULL
        if ((!is.null(cpd.data) | map.null) & sum(node.data$type ==
                                                  "compound") > 1) {
          plot.data.cpd = node.map(cpd.data, node.data, node.types = "compound",
                                   node.sum = node.sum)
          if (map.cpdname & !kegg.native) {
            plot.data.cpd$labels = cpdkegg2name(plot.data.cpd$labels)[,
                                                                      2]
            mapped.cnodes = rownames(plot.data.cpd)
            node.data$labels[mapped.cnodes] = plot.data.cpd$labels
          }
          cols.ts.cpd = node.color(plot.data.cpd, limit$cpd,
                                   bins$cpd, both.dirs = both.dirs$cpd, trans.fun = trans.fun$cpd,
                                   discrete = discrete$cpd, low = low$cpd, mid = mid$cpd,
                                   high = high$cpd, na.col = na.col)
        }
        else plot.data.cpd = cols.ts.cpd = NULL
        if (kegg.native) {
          pv.pars = keggview.native(plot.data.gene = plot.data.gene,
                                    cols.ts.gene = cols.ts.gene, plot.data.cpd = plot.data.cpd,
                                    cols.ts.cpd = cols.ts.cpd, node.data = node.data,
                                    pathway.name = pathway.name[i], kegg.dir = kegg.dir,
                                    limit = limit, bins = bins, both.dirs = both.dirs,
                                    discrete = discrete, low = low, mid = mid, high = high,
                                    na.col = na.col, ...)
        }
        else {
          pv.pars = keggview.graph(plot.data.gene = plot.data.gene,
                                   cols.ts.gene = cols.ts.gene, plot.data.cpd = plot.data.cpd,
                                   cols.ts.cpd = cols.ts.cpd, node.data = node.data,
                                   path.graph = gR1, pathway.name = pathway.name[i],
                                   map.cpdname = map.cpdname, split.group = split.group,
                                   limit = limit, bins = bins, both.dirs = both.dirs,
                                   discrete = discrete, low = low, mid = mid, high = high,
                                   na.col = na.col, ...)
        }
        plot.data.gene = cbind(plot.data.gene, cols.ts.gene)
        if (!is.null(plot.data.gene)) {
          cnames = colnames(plot.data.gene)[-(1:8)]
          nsamp = length(cnames)/2
          if (nsamp > 1) {
            cnames[(nsamp + 1):(2 * nsamp)] = paste(cnames[(nsamp +
                                                              1):(2 * nsamp)], "col", sep = ".")
          }
          else cnames[2] = "mol.col"
          colnames(plot.data.gene)[-(1:8)] = cnames
        }
        plot.data.cpd = cbind(plot.data.cpd, cols.ts.cpd)
        if (!is.null(plot.data.cpd)) {
          cnames = colnames(plot.data.cpd)[-(1:8)]
          nsamp = length(cnames)/2
          if (nsamp > 1) {
            cnames[(nsamp + 1):(2 * nsamp)] = paste(cnames[(nsamp +
                                                              1):(2 * nsamp)], "col", sep = ".")
          }
          else cnames[2] = "mol.col"
          colnames(plot.data.cpd)[-(1:8)] = cnames
        }
        out.list[[i]] = list(plot.data.gene = plot.data.gene,
                             plot.data.cpd = plot.data.cpd)
      }
      if (npath == 1)
        out.list = out.list[[1]]
      else names(out.list) = pathway.name
      return(invisible(out.list))
    }


    ####publish_gostplot_repel####
    publish_gostplot_repel<-function (p, highlight_terms = NULL, filename = NULL, width = NA,
                                      height = NA)
    {
      if (!("ggplot" %in% class(p))) {
        warning("Highlighting terms in a Manhattan plot is available for a ggplot object only.\nPlease set 'interactive = F' in the gostplot() function and try again.")
        return(NULL)
      }
      term_id <- logpval <- term_size_scaled <- id <- query <- p_value <- NULL
      if (!is.null(highlight_terms)) {
        if (is.data.frame(highlight_terms)) {
          message("The input 'highlight_terms' is a data.frame and therefore the column 'term_id' will be used for detection.")
          if ("term_id" %in% colnames(highlight_terms)) {
            highlight_terms <- highlight_terms$term_id
          }
          else {
            stop("No column named 'term_id'.")
          }
        }
        df <- p$data
        subdf <- base::subset(df, term_id %in% highlight_terms)
        if (nrow(subdf) == 0) {
          message("None of the term IDs in the 'highlight_terms' was found from the results.")
          return(p)
        }
        highlight_terms <- unique(highlight_terms)
        subdf$id <- match(subdf$term_id, highlight_terms)
        p <- p + ggplot2::geom_point(data = subdf, ggplot2::aes(x = order,
                                                                y = logpval, size = term_size_scaled), pch = 21,
                                     colour = "black")
        p <- p  + ggrepel::geom_text_repel(box.padding = 0.5, max.overlaps = Inf,data = subdf,arrow = arrow(length = unit(0.008, "npc"), type = "closed", ends = "first"),force = 40,
                                                                                    size = 4, colour = "black", fontface = "bold", ggplot2::aes(label = as.character(id)),
                                                                                    hjust = -1.2, vjust = -0.05)
        pseudo_gostres <- list(result = data.frame(subdf), meta = list(query_metadata = list(queries = sapply(unique(subdf$query),
                                                                                                              function(x) NULL))))
        tb <- publish_gosttable(pseudo_gostres, highlight_terms = highlight_terms,
                                use_colors = TRUE, show_columns = c("source", "term_name",
                                                                    "term_size","intersection_size"), filename = NULL, ggplot = FALSE)
        h <- grid::unit.c(grid::unit(1, "null"), sum(tb$heights) +
                            grid::unit(3, "mm"))
        w <- grid::unit.c(grid::unit(1, "null"))
        tg <- gridExtra::grid.arrange(p, tb, ncol = 1, heights = h,
                                      widths = w, newpage = TRUE)
        p <- ggplot2::ggplot() + ggplot2::annotation_custom(tg) +
          ggplot2::geom_blank() + ggplot2::theme_void()
      }
      if (is.null(filename)) {
        return(p)
      }
      else {
        imgtype <- strsplit(basename(filename), split = "\\.")[[1]][-1]
        if (length(imgtype) == 0) {
          filename = paste0(filename, ".pdf")
        }
        if (tolower(imgtype) %in% c("png", "pdf", "jpeg", "tiff",
                                    "bmp")) {
          if (is.na(width)) {
            width = max(grDevices::dev.size()[1], 8)
          }
          if (is.na(height)) {
            height = max(grDevices::dev.size()[2], 6)
          }
          ggplot2::ggsave(filename = filename, plot = p, width = width,
                          height = height, limitsize = F)
          message("The image is saved to ", filename)
          return(p)
        }
        else {
          stop("The given file format is not supported.\nPlease use one of the following extensions: .png, .pdf, .jpeg, .tiff, .bmp")
        }
      }
    }

    ####merge.png.pdf####
    merge.png.pdf <- function(pdfFile, pngFiles,imagename, deletePngFiles=FALSE) {

      #### Package Install ####
      pngPackageExists <- require ("png")
      suppressMessages(library("grid"))

      if ( !pngPackageExists ) {
        install.packages ("png")
        library ("png")
      }

      pdf(pdfFile)
      n <- length(pngFiles)
      plot.new()
      for( i in 1:n) {
        pngFile <- pngFiles[i]
        pngRaster <- readPNG(pngFile)
        grid.raster(pngRaster, width=unit(0.5, "npc"),height= unit(0.5, "npc"))
        text(x=0.25,y=0, labels=imagename)
        if (i < n) plot.new()
      }
      dev.off()

      if (deletePngFiles) {
        unlink(pngFiles)
      }
    }

    f <- file.path(args, filelist)
    query <- lapply(f, read.table)
    # lapply(query, names)
    names(query) <- gsub(".*/(.*)\\..*", "\\1", f)
    lapply(seq_along(query), function(i)names(query[[i]]$V1)<-names(query[i]))
    pp<-list()
    pdfFiles <- c()

    print("start for loop for graphs saving...")
    for (i in 1:length(query)){
      # names(query[[i]][,1])<-names(query[i])
      tmp_list <- list(genes=query[[i]][,1])
      names(tmp_list)<-names(query[i])
      gostres <- gost(query = tmp_list,
                      organism = species, ordered_query = FALSE,
                      multi_query = FALSE, significant = TRUE, exclude_iea = FALSE,
                      measure_underrepresentation = FALSE, evcodes = TRUE,
                      user_threshold = 0.05, correction_method = "g_SCS",
                      domain_scope = "annotated", custom_bg = NULL,
                      numeric_ns = "", sources = NULL)
      if(!is.null(gostres)){
        highlight_terms<- gostres$result %>% filter(source == "REAC") %>% select(term_id) %>% head(10)
        p <- gostplot(gostres, capped = TRUE, interactive = FALSE)

        pdf(paste(args,"/",names(query[i]),"_gprofiler.pdf",sep=""),width = 12, height = 8)
        pp[[i]] <- publish_gostplot_repel(p, highlight_terms =highlight_terms$term_id ,width = NA, height = NA, filename = NULL )
        dev.off()

        gostres$result$parents<-unlist(lapply(seq_along(gostres$result$parents),function(i){paste(gostres$result$parents[[i]],collapse="/")}))
        write.csv(as.data.frame(gostres$result),paste(args,"/",names(query[i]),"_gprofiler.csv",sep=""),row.names = F)

        ####add a GO bubble plot####
        to_bubble<-as.data.frame(gostres$result)
        to_bubble<-to_bubble %>% filter(source=="REAC") %>% select(term_name, intersection_size,p_value) %>% head(10)
        to_bubble$p_value<-log10(to_bubble$p_value)
        to_bubble$p_value<-(-to_bubble$p_value)
        colnames(to_bubble)[3]<-"-log10_p_val"
        mycol=brewer.pal(9, "Set1")
        p_bubble<-ggplot(to_bubble, aes(x=`-log10_p_val`, y=reorder(term_name,`-log10_p_val`), size = intersection_size,colour =`-log10_p_val`)) +
          ylab("term name")+
          xlab("")+
          geom_point()+
          theme_bw()+
          theme(plot.margin = unit(c(8,4,8,4), "cm"))+
          scale_color_gradient(low = mycol[2],high = mycol[1])

        p_bubble
        ggsave(paste(args,"/",names(query[i]),"_top10_REAC_bubble.pdf",sep=""),width = 12, height = 10)
        # dev.off()

        to_bubble<-as.data.frame(gostres$result)
        to_bubble_list<-split(to_bubble, to_bubble$source)
        lapply(seq_along(to_bubble_list), function(x){
          to_bubble_list[[x]]<<-head(to_bubble_list[[x]],4)
        })
        to_bubble<-rbind(to_bubble_list$`GO:CC`,to_bubble_list$`GO:BP`,to_bubble_list$`GO:MF`)
        to_bubble<-to_bubble %>% select(term_name, intersection_size,p_value,source)
        to_bubble$p_value<-log10(to_bubble$p_value)
        to_bubble$p_value<-(-to_bubble$p_value)
        colnames(to_bubble)[3]<-"-log10_p_val"
        mycol=brewer.pal(9, "Set1")
        p_bubble<-ggplot(to_bubble, aes(x=`-log10_p_val`, y=reorder(term_name,`-log10_p_val`), size = intersection_size,colour =`-log10_p_val`,shape=source)) +
          ylab("term name")+
          xlab("")+
          geom_point()+
          theme_bw()+
          theme(
            plot.margin = unit(c(8,5,8,5), "cm")
          )+
          scale_color_gradient(low = mycol[2],
                               high = mycol[1])

        p_bubble
        ggsave(paste(args,"/",names(query[i]),"_top12GO_CC_BP_MF_bubble.pdf",sep=""),width = 12, height = 10)

        }else{
        print(paste("No results to show for ",names(query[i]),". reasons: 1)wrong input gene set, check your input gene set 2)No functional enrichment result for the input genes in deed",sep=""))
      }

      ## KEGG for target genes##
      top_KEGG_id<-gostres$result %>% filter(source == "KEGG") %>% select(term_id) %>% head(10)

      if(!is.null(top_KEGG_id) & length(query[[i]])==2){
        if(is.numeric(query[[i]][,2])){
          top_KEGG_id$term_id<-gsub('KEGG:','hsa',top_KEGG_id$term_id)
          top_KEGG<-gostres$result %>% filter(source == "KEGG") %>% select(c("term_id","term_name","intersection"))%>% head(10)
          write.csv(as.data.frame(top_KEGG),paste(args,"/",names(query[i]),"_topKEGG.csv",sep=""),row.names = F)
          top_KEGG$term_id<-top_KEGG_id$term_id

          gene_list<-as.list(top_KEGG$intersection)
          names(gene_list)<-top_KEGG$term_id

          toppathway<-top_KEGG_id$term_id

          pngFiles <- c()
          lapply(seq_along(toppathway),function(k){

            eg<-id2eg(ids=unlist(str_split(gene_list[[toppathway[k]]],",")), category='SYMBOL', org='Hs')[,"ENTREZID"]
            # FC<-rnorm(length(eg),-0.5,1)
            FC<-query[[i]][,2]
            names(FC)<-eg
            pathview2(gene.data = FC,
                      kegg.dir = args,
                      pathway = toppathway[k],
                      species = "hsa",
                      out.suffix = gsub(" ","_",top_KEGG$term_name[k]),
                      kegg.native = TRUE)
            flist <- list.files("./", paste("hsa.*.",gsub(" ","_",top_KEGG$term_name[k]),sep=""), full.names = TRUE)
            file.copy(flist,args)
            file.remove(flist)


            pngFile <- paste(paste0(args,"/",toppathway[k]), gsub(" ","_",top_KEGG$term_name[k]),"png",sep=".")
            pngFiles[k] <<- pngFile

          })
          kegg_pathname <- paste(args,"/",names(query[i]),"_pathview",".pdf",sep="")
          pdfFiles[i]<-kegg_pathname
          merge.png.pdf(pdfFile = kegg_pathname, pngFiles = pngFiles,imagename = names(query[i]), deletePngFiles = F)
        }else{
          print("column 2 of input genes is not numeric! Please set them as numeric, such as the stop rate fold change or others.")
        }

      }



    }


    pp<- pp[lapply(pp,length)>0]
    ggsave(
      filename = paste(args,"/",paste(names(query),collapse = "_"),".pdf",sep=""),
      plot = marrangeGrob(pp, nrow=1, ncol=1),
      width = 12, height = 8
    )

    if(!is.null(pdfFiles)){
      write.table(pdfFiles,paste(args,"/",paste(names(query),collapse = "_"),".txt",sep=""),row.names=F,col.names=F,quote=F)
    }


Output
--------

Information
************

Result with ``_gprofiler.csv`` suffix is the final functional enrichment result.

.. code:: bash

    $ cd /the/directory/of/out_file_dir
    $ tree -L 1
    .
    ├── Day4_com_gprofiler.csv
    ├── Day4_com_gprofiler.pdf
    ├── Day4_com.pdf
    ├── Day4_com_top10_REAC_bubble.pdf
    ├── Day4_com_top12GO_CC_BP_MF_bubble.pdf
    └── Day4_com.txt

    0 directories, 6 files

Diagram
************
File with suffix ``_gprofiler.pdf`` is a gprofiler graphical summary of Ψ-sites functional enrichment on input gene id.

.. image:: /images/Functional_Enrichment_gprofiler.png

File with suffix ``_top10_REAC_bubble.pdf`` is a bubble plot of Ψ-sites functional enrichment on input gene id.

.. image:: /images/Functional_Enrichment_REAC_bubble.png

File with suffix ``_top12GO_CC_BP_MF_bubble.pdf`` is a bubble plot of Ψ-sites functional enrichment on input gene id.

.. image:: /images/Functional_Enrichment_GO_bubble.png

File with suffix ``_pathview.pdf`` is a KEGG pathway summary plot of Ψ-sites functional enrichment on input gene id.

.. image:: /images/Functional_Enrichment_pathview.png

.. note:: All user input will be recorded in a plain text file with suffix ``_pseUfun_config.txt`` in psiFinder/config and help users to easily reload the previous config (by simply clicking ``CONFIG`` button).
