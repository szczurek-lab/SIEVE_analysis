# Title     : Real data visualization.
# Objective : Provide helper functions.
# Created by: senbaikang
# Created on: 20.03.21

repo <- "https://stat.ethz.ch/CRAN/"

if (!"BiocManager" %in% installed.packages())
    install.packages("BiocManager")

if(!"ggplot2" %in% installed.packages())
  install.packages("ggplot2", repos = repo)
library(ggplot2, quietly = TRUE)

if(!"ggrepel" %in% installed.packages())
  install.packages("ggrepel", repos = repo)
library(ggrepel, quietly = TRUE)

if(!"ggtree" %in% installed.packages())
  BiocManager::install("ggtree")
library(ggtree, quietly = TRUE)

if(!"ggnewscale" %in% installed.packages())
  BiocManager::install("ggnewscale")
library(ggnewscale, quietly = TRUE)

if(!"VariantAnnotation" %in% installed.packages())
  BiocManager::install("VariantAnnotation")
library(VariantAnnotation, quietly = TRUE)

if(!"ape" %in% installed.packages())
  install.packages("ape", repos = repo)
library(ape, quietly = TRUE)

if(!"treeio" %in% installed.packages())
  install.packages("treeio", repos = repo)
library(treeio, quietly = TRUE)

if(!"tibble" %in% installed.packages())
  install.packages("tibble", repos = repo)
library(tibble, quietly = TRUE)

if(!"dplyr" %in% installed.packages())
  install.packages("dplyr", repos = repo)
library(dplyr, quietly = TRUE)

if(!"tidyr" %in% installed.packages())
  install.packages("tidyr", repos = repo)
library(tidyr, quietly = TRUE)

if(!"pbapply" %in% installed.packages())
  install.packages("pbapply", repos = repo)
library(pbapply, quietly = TRUE)

if(!"RColorBrewer" %in% installed.packages())
  install.packages("RColorBrewer", repos = repo)
library(RColorBrewer, quietly = TRUE)

if(!"wesanderson" %in% installed.packages())
  install.packages("wesanderson", repos = repo)
library(wesanderson, quietly = TRUE)

if(!"gridExtra" %in% installed.packages())
  install.packages("gridExtra", repos = repo)
library(gridExtra, quietly = TRUE)

if(!"stringr" %in% installed.packages())
  install.packages("stringr", repos = repo)
library(stringr, quietly = TRUE)

if(!"fastcluster" %in% installed.packages())
  BiocManager::install("fastcluster")
library(fastcluster, quietly = TRUE)

if(!"ComplexHeatmap" %in% installed.packages())
  BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap, quietly = TRUE)
library(circlize, quietly = TRUE)

if(!"gridtext" %in% installed.packages())
  BiocManager::install("gridtext")
library(gridtext, quietly = TRUE)

if(!"magick" %in% installed.packages())
  install.packages("magick", repos = repo)
library(magick, quietly = TRUE)

if(!"tiff" %in% installed.packages())
  install.packages("tiff", repos = repo)
library(tiff, quietly = TRUE)


get.ggtree <- function(
  tree,
  layout = "rectangular",
  branch.alpha = 0.65,
  branch.size = 1.0 / .pt,
  plot.branch.length = TRUE,
  scale.longest.branch = 1
) {
  if (plot.branch.length) {
    branch.length <- "branch.length"

    if ("length" %in% colnames(tree@data)) {
      tree@data <- tree@data %>%
          mutate(ori_length = length)
    }

    if (scale.longest.branch > 1) {
      branch.lengths <- sort(tree@phylo$edge.length, decreasing = TRUE)
      message(paste0("Threshold of branch length: ", branch.lengths[scale.longest.branch]))

      if ("length" %in% colnames(tree@data)) {
        tree@data <- tree@data %>%
          mutate(length = if_else(as.numeric(length) > branch.lengths[scale.longest.branch], as.character(branch.lengths[scale.longest.branch]), length))
      }

      tree@phylo$edge.length[tree@phylo$edge.length > branch.lengths[scale.longest.branch]] <- branch.lengths[scale.longest.branch]
    }
  } else {
    branch.length <- "none"
  }

  tree.plot <- ggtree(
    tree,
    layout = layout,
    branch.length = branch.length,
    aes(
      size = branch,
      alpha = branch.alpha
    )
  ) +
    scale_size_continuous(
      range = c(branch.size, branch.size)
    ) +
    guides(
      size = "none",
      alpha = "none"
    )

  tree.plot
}


normalize.support <- function(vals) {
  sapply(
    vals,
    function(x) {
      ifelse(grepl("\\d+", x), as.character(as.numeric(x) / 100), "")
    },
    USE.NAMES = FALSE
  )
}


get.tip.color <- function(i = 1) {
  if (i == 1) {
    return(
      c("#FFA319FF", "#350E20FF", "#767676FF", "#C16622FF")
    )
  } else if (i == 2) {
    return(
      c("#F39B7FFF", "#7E6148FF", "#B09C85FF", "#8491B4FF")
    )
  } else if (i == 3) {
    return(
      c("#440154FF", "#472D7BFF", "#3A538BFF", "#20928CFF", "#62CB5FF")
    )
  }
}

plot.nodes <- function(
  tree.plot,
  anno.tip.text = TRUE,
  anno.tip.label = TRUE,
  tip.data = NULL,
  tip.point.size = 2,
  tip.label.align = TRUE,
  tip.label.linetype = "dotted",
  tip.label.linesize = 0.3,
  tip.label.offset = 0.0000008,
  tip.label.size = 6,
  tip.label.legend.title = "Cell origin",
  tip.label.color = character(),
  x.limit = c(0, NA),
  legend.spacing = unit(0, "mm"),
  legend.pos = c(0, 1),
  legend.just = c("left", "top"),
  legend.direction = "vertical",
  legend.key.size = unit(2, "pt"),
  legend.text.size = element_text(size = 6),
  lengend.title.size = element_text(size = 7),
  plot.margin = unit(c(0, 0, 0, 0), "mm"),
  anno.hpd = FALSE,
  hpd.color = "red",
  hpd.size = 2.5 / .pt,
  hpd.alpha = 1.0,
  anno.branch.length = FALSE,
  branch.length.digit = 3L,
  branch.length.text.size = 5L / .pt,
  branch.length.color = "lightblue",
  branch.length.alpha = 1.0,
  label.repel.arrow.length = unit(1, "mm"),
  label.repel.force = 1,
  label.repel.force_pull = 1,
  label.repel.min.segment.length = 0,
  label.repel.max.overlaps = 10,
  label.repel.label.padding = unit(0.1, "lines"),
  label.repel.box.padding = unit(0.1, "lines"),
  label.repel.label.r = unit(0.1, "lines"),
  label.repel.label.size = unit(0.15, "mm"),
  anno.support = TRUE,
  support.var = "posterior",
  support.threshold = 0.5,
  internal.node.text.size = 5,
  internal.node.label.padding = unit(0.1, "lines"),
  internal.node.label.r = unit(0.1, "lines"),
  internal.node.label.size = unit(0.15, "mm"),
  internal.node.label.alpha = 1.0,
  internal.node.label.fill.color = "lightblue",
  seed = 1L,
  out.file = NULL,
  out.format = "pdf",
  out.width = 180,
  out.height = 80,
  out.dpi = 300,
  out.limitsize = FALSE
) {
  if (!is.null(tip.data)) {
    ret.plot <- tree.plot %<+% tip.data
  } else {
    ret.plot <- tree.plot
  }

  if (anno.tip.label) {
    if (!is.null(tip.data)) {
      ret.plot <- ret.plot +
        geom_tippoint(
          inherit.aes = FALSE,
          aes(
            x = x,
            y = y,
            color = cell_groups
          ),
          size = tip.point.size,
          alpha = 1.0
        )
    } else {
      ret.plot <- ret.plot +
        geom_tippoint(
          inherit.aes = FALSE,
          aes(
            x = x,
            y = y
          ),
          size = tip.point.size,
          alpha = 1.0
        )
    }
  }

  if (anno.tip.text) {
    if (!is.null(tip.data)) {
      ret.plot <- ret.plot +
        geom_tiplab(
          inherit.aes = FALSE,
          aes(
            x = x,
            y = y,
            color = cell_groups
          ),
          align = tip.label.align,
          linetype = tip.label.linetype,
          linesize = tip.label.linesize,
          offset = tip.label.offset,
          size = tip.label.size,
          alpha = 1.0
        )
    } else {
      ret.plot <- ret.plot +
        geom_tiplab(
          inherit.aes = FALSE,
          aes(
            x = x,
            y = y
          ),
          align = tip.label.align,
          linetype = tip.label.linetype,
          linesize = tip.label.linesize,
          offset = tip.label.offset,
          size = tip.label.size,
          alpha = 1.0
        )
    }
  }

  if (anno.tip.label | anno.tip.text) {
    ret.plot <- ret.plot +
      scale_color_manual(
        name = tip.label.legend.title,
        values = tip.label.color
      )
  }

  ret.plot <- ret.plot +
    xlim(x.limit) +
    theme(
      legend.spacing = legend.spacing,
      legend.position = legend.pos,
      legend.justification = legend.just,
      legend.direction = legend.direction,
      legend.key.size = legend.key.size,
      legend.text = legend.text.size,
      title = lengend.title.size,
      plot.margin = plot.margin
    )

  if (anno.hpd) {
    ret.plot <- ret.plot +
      geom_range(
        "length_0.95_HPD",
        color = hpd.color,
        size = hpd.size,
        alpha = hpd.alpha
      )
  }

  if (anno.branch.length) {
    ret.plot <- ret.plot +
      geom_label_repel(
        inherit.aes = FALSE,
        show.legend = FALSE,
        data = ret.plot$data[which(as.numeric(ret.plot$data$ori_length) > 0),],
        mapping = aes(
          x = branch,
          y = y,
          label = signif(as.numeric(ori_length), branch.length.digit)
        ),
        size = branch.length.text.size,
        fill = branch.length.color,
        alpha = branch.length.alpha,
        segment.size = 1 / .pt,
        # arrow = arrow(
        #   length = label.repel.arrow.length,
        #   type = "open"
        # ),
        force = label.repel.force,
        force_pull = label.repel.force_pull,
        max.overlaps = label.repel.max.overlaps,
        box.padding = label.repel.box.padding,
        label.padding = label.repel.label.padding,
        label.r = label.repel.label.r,
        label.size = label.repel.label.size,
        min.segment.length = label.repel.min.segment.length,
        seed = seed
      )
  }

  if (anno.support) {
    ret.plot <- ret.plot +
      geom_label2(
        aes(
          x = x,
          y = y,
          label = format(round(as.numeric(get(support.var)), 2), nsmall = 2),
          subset = as.numeric(get(support.var)) > support.threshold
        ),
        size = internal.node.text.size,
        label.padding = internal.node.label.padding,
        label.r = internal.node.label.r,
        label.size = internal.node.label.size,
        alpha = internal.node.label.alpha,
        fill = internal.node.label.fill.color,
        inherit.aes = FALSE
      )
  }

  if (!is.null(out.file)) {
    ggsave(
      out.file,
      plot = ret.plot,
      device = out.format,
      width = out.width,
      height = out.height,
      units = "mm",
      dpi = out.dpi,
      limitsize = out.limitsize
    )
  }

  ret.plot
}


get.tip.info <- function(
  tree,
  pat2names = NULL
) {
  if (is.null(pat2names)) return(NULL)

  tip.names <- get.tree(tree)$tip.label

  tip.groups <- rep("", length(tip.names))
  for (pat2name in pat2names) {
    tip.groups[grep(pat2name[1], tip.names)] <- pat2name[2]
  }

  return(
    data.frame(
      cell_names = tip.names,
      cell_groups = factor(tip.groups)
    )
  )
}


get.ordered.tip.names <- function(tree) {
  tree.tmp <- subset(fortify(tree), isTip)
  return(tree.tmp$label[order(tree.tmp$y, decreasing = TRUE)])
}


.combine.duplicate.genes <- function(
  x,
  sep = " \U00D7 ",
  criterion = tibble(
    num = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200),
    group = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20)
  ),
  collapse = "  "
) {
  if (all(is.na(x))) return(NA)

  .x <- tibble(ori = sort(x)) %>%
    group_by(ori) %>%
    mutate(id = 1:n()) %>%
    mutate(str = ifelse(max(id) > 1, paste(max(id), ori, sep = sep), ori)) %>%
    ungroup() %>%
    select(ori, str) %>%
    unique()

  .num.group <- criterion$group[(function(i, j){
    for (.i in seq_len(length(i))) {
      if (j <= i[.i]) return(.i)
    }
  })(criterion$num, length(.x$str))]

  tryCatch(
    {
      .gene.groups <- tibble(str = .x$str) %>%
        mutate(id = 1:n()) %>%
        mutate(group = (id - 1) %/% .num.group) %>%
        split(f = .$group)
    },
    error = function(cond) {
      message(cond)
    },
    warning = function(cond) {
      message(cond)
    }
  )

  ret <- sapply(.gene.groups, function(i) {
    paste(i$str, collapse = collapse)
  })
  ret
}


#` mode: 1 - extract all genes; 2 - extract genes following ISA; 3 - extract genes following FSA
.extract.genes <- function(
  tree.plot,
  mark.complex.genes = TRUE,
  complex.genes = list(
    data = NULL,
    et = "event_type_fsa",
    ge = "gene_fsa"
  ),
  .names = list(
    et = "event_type",
    eet = "expanded_event_type",
    ge = "gene",
    cg = "collapsed_gene"
  ),
  dup.gene.sep = " \U00D7 ",
  criterion = tibble(
    num = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100),
    group = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
  ),
  collapse = "  "
) {
  df <- separate_rows(tree.plot$data, !!as.name(.names$et), !!as.name(.names$ge), sep = "; ")

  if (mark.complex.genes & !is.null(complex.genes$data)) {
    df <- df %>%
      left_join(
        complex.genes$data %>% mutate(is_parallel = TRUE),
        by = c(setNames(nm = .names$et, complex.genes$et), setNames(nm = .names$ge, complex.genes$ge))
      ) %>%
      mutate(
        !!as.name(.names$et) := case_when(
          is.null(is_parallel) | !is_parallel ~ !!as.name(.names$et),
          is_parallel ~ paste0("P", !!as.name(.names$et)),
          TRUE ~ !!as.name(.names$et)
        )
      )
  }

  df <- df %>%
    group_by(node, !!as.name(.names$et)) %>%
    mutate(
      !!as.name(.names$cg) := paste(
        .combine.duplicate.genes(
          !!as.name(.names$ge),
          sep = dup.gene.sep,
          collapse = collapse,
          criterion = criterion
        ),
        collapse = "\n"
      )
    ) %>%
    select(parent, node, label, !!as.name(.names$et), isTip, x, y, branch, angle, !!as.name(.names$cg)) %>%
    unique() %>%
    drop_na(!!as.name(.names$et)) %>%
    mutate(!!as.name(.names$eet) := case_when(
      grepl("^SM$", !!as.name(.names$et), ignore.case = TRUE) ~ "Single mutation (0/0 -> 0/1)",
      grepl("^PSM$", !!as.name(.names$et), ignore.case = TRUE) ~ "Parallel single mutation (0/0 -> 0/1)",
      grepl("^HoSDM$", !!as.name(.names$et), ignore.case = TRUE) ~ "Homozygous simultaneous double mutation (0/0 -> 1/1)",
      grepl("^PHoSDM$", !!as.name(.names$et), ignore.case = TRUE) ~ "Parallel homozygous simultaneous double mutation (0/0 -> 1/1)",
      grepl("^HeSDM$", !!as.name(.names$et), ignore.case = TRUE) ~ "Heterozygous simultaneous double mutation (0/0 -> 1/1')",
      grepl("^PHeSDM$", !!as.name(.names$et), ignore.case = TRUE) ~ "Parallel heterozygous simultaneous double mutation (0/0 -> 1/1')",
      grepl("^SB$", !!as.name(.names$et), ignore.case = TRUE) ~ "Single back mutation (0/1 -> 0/0)",
      grepl("^PSB$", !!as.name(.names$et), ignore.case = TRUE) ~ "Parallel single back mutation (0/1 -> 0/0)",
      grepl("^DB$", !!as.name(.names$et), ignore.case = TRUE) ~ "Double back mutation (1/1 -> 0/0)",
      grepl("^PDB$", !!as.name(.names$et), ignore.case = TRUE) ~ "Parallel double back mutation (1/1 -> 0/0)",
      grepl("^HoSMA$", !!as.name(.names$et), ignore.case = TRUE) ~ "Homozygous single mutation addition (0/1 -> 1/1)",
      grepl("^PHoSMA$", !!as.name(.names$et), ignore.case = TRUE) ~ "Parallel homozygous single mutation addition (0/1 -> 1/1)",
      grepl("^HeSMA$", !!as.name(.names$et), ignore.case = TRUE) ~ "Heterozygous single mutation addition (0/1 -> 1/1')",
      grepl("^PHeSMA$", !!as.name(.names$et), ignore.case = TRUE) ~ "Parallel heterozygous single mutation addition (0/1 -> 1/1')",
      grepl("^HoSSM$", !!as.name(.names$et), ignore.case = TRUE) ~ "Homozygous substitute single mutation (1/1' -> 1/1)",
      grepl("^PHoSSM$", !!as.name(.names$et), ignore.case = TRUE) ~ "Parallel homozygous substitute single mutation (1/1' -> 1/1)",
      grepl("^HeSSM$", !!as.name(.names$et), ignore.case = TRUE) ~ "Heterozygous substitute single mutation (1/1 -> 1/1')",
      grepl("^PHeSSM$", !!as.name(.names$et), ignore.case = TRUE) ~ "Parallel heterozygous substitute single mutation (1/1 -> 1/1')",
      grepl("^D$", !!as.name(.names$et), ignore.case = TRUE) ~ "Deletion",
      grepl("^PD$", !!as.name(.names$et), ignore.case = TRUE) ~ "Parallel deletion",
      grepl("^I$", !!as.name(.names$et), ignore.case = TRUE) ~ "Insertion",
      grepl("^PI$", !!as.name(.names$et), ignore.case = TRUE) ~ "Parallel insertion",
      TRUE ~ !!as.name(.names$et)
    ))
  df
}


get.color.for.genes <- function(name = "Spectral") {
  .tmp <- brewer.pal(n = 9, name = name)
  names(.tmp) <- c(
    "SM",
    "HoSDM",
    "HeSDM",
    "SB",
    "DB",
    "HoSMA",
    "HeSMA",
    "HoSSM",
    "HeSSM"
  )
  .tmp
}


get.color.for.genes.2 <- function() {
  .tmp <- c(brewer.pal(n = 3, name = "Paired")[2], brewer.pal(n = 8, name = "Dark2"), brewer.pal(n = 12, name = "Paired")[c(1, 3:12)])
  names(.tmp) <- c(
    "Single mutation (0/0 -> 0/1)",
    "Homozygous substitute single mutation (1/1' -> 1/1)",
    "Homozygous simultaneous double mutation (0/0 -> 1/1)",
    "Heterozygous simultaneous double mutation (0/0 -> 1/1')",
    "Parallel single mutation (0/0 -> 0/1)",
    "Single back mutation (0/1 -> 0/0)",
    "Homozygous single mutation addition (0/1 -> 1/1)",
    "Heterozygous single mutation addition (0/1 -> 1/1')",
    "Double back mutation (1/1 -> 0/0)",
    "Heterozygous substitute single mutation (1/1 -> 1/1')",
    "Parallel homozygous simultaneous double mutation (0/0 -> 1/1)",
    "Parallel heterozygous simultaneous double mutation (0/0 -> 1/1')",
    "Parallel single back mutation (0/1 -> 0/0)",
    "Parallel double back mutation (1/1 -> 0/0)",
    "Parallel homozygous single mutation addition (0/1 -> 1/1)",
    "Parallel heterozygous single mutation addition (0/1 -> 1/1')",
    "Parallel homozygous substitute single mutation (1/1' -> 1/1)",
    "Parallel heterozygous substitute single mutation (1/1 -> 1/1')",
    "Parallel deletion",
    "Parallel insertion"
  )
  .tmp
}


.get.complex.genes <- function(
  data,
  strict.mode = TRUE,
  et = "event_type_fsa",
  ge = "gene_fsa"
) {
  if (is.null(data[[as.name(et)]])) return(NULL)

  df <- separate_rows(data, !!as.name(et), !!as.name(ge), sep = "; ") %>%
    select(node, !!as.name(et), !!as.name(ge)) %>%
    unique() %>%
    drop_na()

  if (strict.mode) {
    df <- df %>%
      group_by(!!as.name(et), !!as.name(ge)) %>%
      mutate(cnt = 1:n()) %>%
      mutate(include = cnt > 1) %>%
      ungroup() %>%
      filter(include == TRUE) %>%
      select(-node, -cnt, -include) %>%
      unique()
  } else {
    df <- df %>%
      select(!!as.name(et), !!as.name(ge)) %>%
      unique()
  }
  df %>% select(!!as.name(et), !!as.name(ge))
}


#` mode: 1 - extract all genes; 2 - extract genes following ISA; 3 - extract genes following FSA
plot.genes <- function(
  tree.plot,
  mode = 1L,
  mark.complex.genes = list(
    general = TRUE, # mark genes or not
    strict = TRUE # only mark genes with parallel evolution or not
  ),
  dup.gene.sep = " \U00D7 ",
  criterion = tibble(
    num = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100),
    group = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
  ),
  collapse = "  ",
  text.alpha = 1.0,
  text.size = 5 / .pt,
  arrow.length = unit(1, "mm"),
  arrow.type = "open",
  text.repel.force = 1,
  text.repel.force_pull = 1,
  text.repel.min.segment.length = 0,
  text.repel.max.overlaps = 10,
  text.repel.label.padding = unit(0.1, "lines"),
  text.repel.box.padding = unit(0.2, "lines"),
  text.repel.point.padding = unit(0, "lines"),
  label.repel.label.r = unit(0.1, "lines"),
  label.repel.label.size = unit(0.15, "mm"),
  legend.key.size = unit(2, "pt"),
  legend.title = "Evolutionary event",
  repel.type = "label", # or "text"
  text.color.type = "fill", # or "color" / "colour"
  direction = "both",
  seed = 1L,
  out.file = NULL,
  out.format = "pdf",
  out.width = 180,
  out.height = 150,
  out.dpi = 300,
  out.limitsize = FALSE
) {
  stopifnot(mode == 1L || mode == 2L || mode == 3L)
  stopifnot(text.color.type == "fill" || text.color.type == "color" || text.color.type == "colour")
  stopifnot(repel.type == "text" || repel.type == "label")

  if (mode == 1L) {
    .names <- list(
      et = "event_type",
      eet = "expanded_event_type",
      ge = "gene",
      cg = "collapsed_gene"
    )
  } else if (mode == 2L) {
    .names <- list(
      et = "event_type_isa",
      eet = "expanded_event_type_isa",
      ge = "gene_isa",
      cg = "collapsed_gene_isa"
    )
  } else if (mode == 3L) {
    .names <- list(
      et = "event_type_fsa",
      eet = "expanded_event_type_fsa",
      ge = "gene_fsa",
      cg = "collapsed_gene_fsa"
    )
  }

  mark.complex.genes$general <- mode == 1L & mark.complex.genes$general

  if (mark.complex.genes$general) {
    complex.genes <- .get.complex.genes(
      tree.plot$data,
      strict.mode = mark.complex.genes$strict,
      et = "event_type_fsa",
      ge = "gene_fsa"
    )
  } else {
    complex.genes <- NULL
  }

  if (!is.null(complex.genes)) {
    message("Genes with parallel evolution:")
    print.data.frame(complex.genes)
  }

  data <- .extract.genes(
    tree.plot = tree.plot,
    mark.complex.genes$general,
    complex.genes = list(
      data = complex.genes,
      et = "event_type_fsa",
      ge = "gene_fsa"
    ),
    .names = .names,
    dup.gene.sep = dup.gene.sep,
    criterion = criterion,
    collapse = collapse
  )

  if (text.color.type == "fill") {
    repel.aes <- aes(
        x = branch,
        y = y,
        label = !!as.name(.names$cg),
        fill = as.factor(!!as.name(.names$eet))
    )
    repel.scale.func <- scale_fill_manual
    ret.plot <- tree.plot
  } else if (text.color.type == "color" || text.color.type == "colour") {
    repel.aes <- aes(
        x = branch,
        y = y,
        label = !!as.name(.names$cg),
        color = as.factor(!!as.name(.names$eet))
    )
    repel.scale.func <- scale_color_manual
    ret.plot <- tree.plot + new_scale_color()
  }

  if (repel.type == "text") {
    ret.plot <- ret.plot +
      geom_text_repel(
        inherit.aes = FALSE,
        data = data,
        mapping = repel.aes,
        alpha = text.alpha,
        size = text.size,
        family = "sans",
        segment.size = 1 / .pt,
        # arrow = arrow(
        #   length = arrow.length,
        #   type = arrow.type
        # ),
        force = text.repel.force,
        force_pull = text.repel.force_pull,
        min.segment.length = text.repel.min.segment.length,
        max.overlaps = text.repel.max.overlaps,
        label.padding = text.repel.label.padding,
        box.padding = text.repel.box.padding,
        point.padding = text.repel.point.padding,
        direction = direction,
        seed = seed
      )
  } else if (repel.type == "label") {
    ret.plot <- ret.plot +
      geom_label_repel(
        inherit.aes = FALSE,
        data = data,
        mapping = repel.aes,
        alpha = text.alpha,
        size = text.size,
        family = "sans",
        segment.size = 1 / .pt,
        # arrow = arrow(
        #   length = arrow.length,
        #   type = arrow.type
        # ),
        force = text.repel.force,
        force_pull = text.repel.force_pull,
        min.segment.length = text.repel.min.segment.length,
        max.overlaps = text.repel.max.overlaps,
        label.padding = text.repel.label.padding,
        box.padding = text.repel.box.padding,
        point.padding = text.repel.point.padding,
        label.r = label.repel.label.r,
        label.size = label.repel.label.size,
        direction = direction,
        seed = seed
      )
  }

  ret.plot <- ret.plot +
    repel.scale.func(
      name = legend.title,
      values = get.color.for.genes.2()[sort(levels(as.factor(data[[.names$eet]])))],
      labels = sort(levels(as.factor(data[[.names$eet]])))
    )

  if (!is.null(out.file)) {
    ggsave(
      out.file,
      plot = ret.plot,
      device = out.format,
      width = out.width,
      height = out.height,
      units = "mm",
      dpi = out.dpi,
      limitsize = out.limitsize
    )
  }

  ret.plot
}


# loadSamples <- function(
#   file,
#   num.cells,
#   num.lines.to.read = 100L,
#   burn.in = 0.1,
#   gt.category = c("0/0", "0/1", "1/1", "1/2"),
#   ado.category = c("0", "1"),
#   cores = 4L,
#   local.copy = NULL
# ) {
#   stopifnot(!is.null(file) && file.exists(file))
#
#   if (burn.in > 1.0 || burn.in < 0.0) burn.in <- 0.1
#
#   .gt.ado.sample.file <- file(file, "r")
#   tryCatch(
#     {
#       f.lines <- readLines(.gt.ado.sample.file, n = 1L)
#       sites <- tail(
#         str_split(
#           string = f.lines,
#           pattern = "\t",
#           simplify = TRUE
#         )[1,],
#         n = -1L
#       )
#
#       num.sites <- length(sites)
#       cell.names <- character()
#       reach.site.map <- FALSE
#       site.map <- character()
#
#       samples <- list()
#
#       while (TRUE) {
#         f.lines <- readLines(.gt.ado.sample.file, n = num.lines.to.read)
#         if (length(f.lines) == 0L) break
#
#         is.comment <- startsWith(x = f.lines, prefix = "#")
#         f.lines.sample <- f.lines[!is.comment]
#         f.lines.comment <- f.lines[is.comment]
#
#         if (length(f.lines.sample) > 0) {
#           samples <- c(samples, pbapply(
#             str_split(
#               string = f.lines.sample,
#               pattern = "\t",
#               simplify = TRUE
#             ),
#             MARGIN = 1L,
#             function(x) {
#               parsed.f.lines <- str_split(
#                 string = tail(
#                   x,
#                   n = -1L
#                 ),
#                 pattern = ";",
#                 simplify = TRUE
#               )
#
#               ret <- list(
#                 gt = integer(length = num.sites * num.cells * length(gt.category)),
#                 ado = integer(length = num.sites * num.cells * length(ado.category))
#               )
#
#               for (i in seq_len(dim(parsed.f.lines)[1])) {
#                 for (j in seq_len(dim(parsed.f.lines)[2])) {
#                   .data <- str_split(
#                     string = parsed.f.lines[i, j],
#                     pattern = ",",
#                     simplify = TRUE
#                   )
#
#                   ret$gt[(i - 1) * num.cells * length(gt.category) + (j - 1) * length(gt.category) + match(.data[1], gt.category)] %+=% 1L
#                   ret$ado[(i - 1) * num.cells * length(ado.category) + (j - 1) * length(ado.category) + match(.data[2], ado.category)] %+=% 1L
#                 }
#               }
#
#               return(ret)
#             },
#             cl = cores
#           ))
#         }
#
#         if (length(f.lines.comment) > 0) {
#           if (reach.site.map)
#             site.map <- c(site.map, .get.site.map(f.lines.comment))
#           else {
#             for (x in seq_len(length(f.lines.comment))) {
#               if (reach.site.map) {
#                 site.map <- c(site.map, .get.site.map(f.lines.comment[x:length(f.lines.comment)]))
#                 break
#               } else if (startsWith(f.lines.comment[x], "#Cells: ")) {
#                 cell.names <- str_split(
#                   string = sub(
#                     pattern = "\\s",
#                     replacement = "",
#                     x = sub(
#                       pattern = "^#Cells: ",
#                       replacement = "",
#                       x = f.lines.comment[x]
#                     )
#                   ),
#                   pattern = ",",
#                   simplify = TRUE
#                 )[1, ]
#               } else if (startsWith(f.lines.comment[x], "#Sites map "))
#                 reach.site.map <- TRUE
#             }
#           }
#         }
#       }
#
#       samples <- tail(
#         samples,
#         n = as.integer(-ceiling(length(samples) * burn.in))
#       )
#
#       ret <- list(
#         gt = array(
#           data = rowSums(sapply(samples, "[[", 1L)),
#           dim = c(num.sites, num.cells, length(gt.category)),
#           dimnames = list(
#             site.map[sites],
#             cell.names,
#             gt.category
#           )
#         ),
#         ado = array(
#           data = rowSums(sapply(samples, "[[", 2L)),
#           dim = c(num.sites, num.cells, length(ado.category)),
#           dimnames = list(
#             site.map[sites],
#             cell.names,
#             ado.category
#           )
#         )
#       )
#     },
#     finally = close(.gt.ado.sample.file)
#   )
#
#   if (!is.null(local.copy))
#     saveRDS(
#       object = ret,
#       file = local.copy
#     )
#
#   return(ret)
# }
#
#
# loadSamples2 <- function(
#   file,
#   num.cells,
#   num.lines.to.read = 100L,
#   burn.in = 0.1,
#   per.sample.lapply = function(l, ...) pblapply(l, cl = 8L, ...),
#   per.site.lapply = function(l, ...) lapply(l, ...),
#   local.copy = NULL
# ) {
#   stopifnot(!is.null(file) && file.exists(file))
#
#   .start <- Sys.time()
#
#   if (burn.in > 1.0 || burn.in < 0.0) burn.in <- 0.1
#
#   .gt.ado.sample.file <- file(file, "r")
#   tryCatch(
#     {
#       f.lines <- readLines(.gt.ado.sample.file, n = 1L)
#
#       psedu.site.names <- tail(
#         str_split(
#           string = f.lines,
#           pattern = "\t",
#           simplify = TRUE
#         )[1L,],
#         n = -1L
#       )
#       num.sites <- length(psedu.site.names)
#       psedu.cell.names <- paste0("cell", 1L:num.cells)
#       reach.site.map <- FALSE
#       site.map <- list()
#
#       gt.ado.sample <- list()
#
#       iter.sample <- 0L
#       iter.site.map <- 0L
#       while (TRUE) {
#         f.lines <- readLines(.gt.ado.sample.file, n = num.lines.to.read)
#         if (length(f.lines) == 0L) break
#
#         is.comment <- startsWith(x = f.lines, prefix = "#")
#         f.lines.sample <- f.lines[!is.comment]
#         f.lines.comment <- f.lines[is.comment]
#
#         if (length(f.lines.sample) > 0L) {
#           iter.sample %+=% 1L
#
#           sample.site <- str_split(
#             string = f.lines.sample,
#             pattern = "\t",
#             simplify = TRUE
#           )
#
#           stopifnot(dim(sample.site)[2L] == num.sites + 1L)
#
#           gt.ado.sample[[iter.sample]] <- do.call(
#             rbind,
#             per.sample.lapply(
#               seq_len(dim(sample.site)[1L]),
#               function(x) {
#                 site.cell <- str_split(
#                   string = sample.site[x, 2L:(num.sites + 1L)],
#                   pattern = ";",
#                   simplify = TRUE
#                 )
#
#                 stopifnot(dim(site.cell)[2L] == num.cells)
#
#                 do.call(
#                   rbind,
#                   per.site.lapply(
#                     seq_len(num.sites),
#                     function(y) {
#                       as.data.frame(t(sapply(
#                         seq_len(num.cells),
#                         function(z) {
#                           gt.ado <- str_split(
#                             string = site.cell[y, z],
#                             pattern = ",",
#                             simplify = TRUE
#                           )[1L,]
#
#                           ret <- c(sample.site[x, 1L], psedu.site.names[y], psedu.cell.names[z], gt.ado[1L], gt.ado[2L])
#                           names(ret) <- c("sample", "psedu_site", "psedu_cell", "gt", "ado")
#                           return(ret)
#                         }
#                       )))
#                     }
#                   )
#                 )
#               }
#             )
#           )
#         }
#
#         if (length(f.lines.comment) > 0L) {
#           if (reach.site.map) {
#             iter.site.map %+=% 1L
#             site.map[[iter.site.map]] <- .get.site.map2(f.lines.comment)
#           } else {
#             for (x in seq_len(length(f.lines.comment))) {
#               if (reach.site.map) {
#                 iter.site.map %+=% 1L
#                 site.map[[iter.site.map]] <- .get.site.map2(f.lines.comment[x:length(f.lines.comment)])
#                 break
#               } else if (startsWith(f.lines.comment[x], "#Cells: ")) {
#                 cell.names <- tibble(
#                   psedu_cell = psedu.cell.names,
#                   cell = str_split(
#                     string = sub(
#                       pattern = "\\s",
#                       replacement = "",
#                       x = sub(
#                         pattern = "^#Cells: ",
#                         replacement = "",
#                         x = f.lines.comment[x]
#                       )
#                     ),
#                     pattern = ",",
#                     simplify = TRUE
#                   )[1, ]
#                 )
#               } else if (startsWith(f.lines.comment[x], "#Sites map ")) {
#                 reach.site.map <- TRUE
#               }
#             }
#           }
#         }
#       }
#
#       gt.ado.sample <- do.call(rbind, gt.ado.sample)
#       gt.ado.sample$sample <- as.integer(gt.ado.sample$sample)
#       gt.ado.sample$ado <- as.integer(gt.ado.sample$ado)
#       site.map <- do.call(rbind, site.map)
#
#       samples <- sort(unique(gt.ado.sample$sample))
#       gt.ado.sample <- gt.ado.sample %>%
#         filter(sample > samples[as.integer(ceiling(length(samples) * burn.in))]) %>%
#         left_join(cell.names, by = "psedu_cell") %>%
#         left_join(site.map, by = "psedu_site") %>%
#         select(-psedu_cell, -psedu_site) %>%
#         rename(site = collapsed_site)
#
#       gt.ado.names <- c("gt", "ado")
#       ret <- pbsapply(
#         gt.ado.names,
#         function(x) {
#           excluded <- gt.ado.names[x != gt.ado.names]
#           gt.ado.sample %>%
#             select(-!!as.name(excluded)) %>%
#             add_count(cell, site, !!as.name(x), name = "count") %>%
#             select(-sample) %>%
#             unique() %>%
#             group_by(cell, site) %>%
#             mutate(freq = count / sum(count)) %>%
#             ungroup() %>%
#             arrange(chr, pos, cell, !!as.name(x))
#         },
#         simplify = FALSE,
#         USE.NAMES = TRUE,
#         cl = 2L
#       )
#     },
#     finally = close(.gt.ado.sample.file)
#   )
#
#   if (!is.null(local.copy))
#     saveRDS(
#       object = ret,
#       file = local.copy
#     )
#
#   .end <- Sys.time()
#   message(paste0("Finished in ", (.end - .start), " seconds."))
#
#   return(ret)
# }
#
#
# loadSamples3 <- function(
#   file,
#   num.cells,
#   burn.in = 0.1,
#   local.copy = NULL
# ) {
#   stopifnot(!is.null(file) && file.exists(file))
#
#   .start <- Sys.time()
#
#   if (burn.in > 1.0 || burn.in < 0.0) burn.in <- 0.1
#
#   psedu.cell.names <- paste0("cell", 1L:num.cells)
#
#   gt.ado.sample <- read.table(
#     file = file,
#     header = TRUE
#   ) %>%
#     column_to_rownames(var = "Sample")
#
#   gt.ado.sample <- tibble(
#       data = as.vector(as.matrix(gt.ado.sample)),
#       sample = as.integer(rep(rownames(gt.ado.sample), dim(gt.ado.sample)[2L])),
#       psedu_site = as.vector(sapply(colnames(gt.ado.sample), function(x) rep(x, dim(gt.ado.sample)[1L]))),
#       psedu_cell = rep(paste(psedu.cell.names, collapse = ";"), dim(gt.ado.sample)[1L] * dim(gt.ado.sample)[2L])
#     ) %>%
#     separate_rows(data, psedu_cell, sep = ";") %>%
#     separate(col = data, into = c("gt", "ado"), sep = ",", convert = TRUE)
#
#   .gt.ado.sample.file <- file(file, "r")
#   tryCatch(
#     {
#       readLines(.gt.ado.sample.file, n = 1L)
#
#       reach.site.map <- FALSE
#       site.map <- list()
#
#       iter.site.map <- 0L
#       while (TRUE) {
#         f.lines <- readLines(.gt.ado.sample.file, n = 1000L)
#         if (length(f.lines) == 0L) break
#
#         is.comment <- startsWith(x = f.lines, prefix = "#")
#         f.lines.comment <- f.lines[is.comment]
#
#         if (length(f.lines.comment) > 0L) {
#           if (reach.site.map) {
#             iter.site.map %+=% 1L
#             site.map[[iter.site.map]] <- .get.site.map2(f.lines.comment)
#           } else {
#             for (x in seq_len(length(f.lines.comment))) {
#               if (reach.site.map) {
#                 iter.site.map %+=% 1L
#                 site.map[[iter.site.map]] <- .get.site.map2(f.lines.comment[x:length(f.lines.comment)])
#                 break
#               } else if (startsWith(f.lines.comment[x], "#Cells: ")) {
#                 cell.names <- tibble(
#                   psedu_cell = psedu.cell.names,
#                   cell = str_split(
#                     string = sub(
#                       pattern = "\\s",
#                       replacement = "",
#                       x = sub(
#                         pattern = "^#Cells: ",
#                         replacement = "",
#                         x = f.lines.comment[x]
#                       )
#                     ),
#                     pattern = ",",
#                     simplify = TRUE
#                   )[1, ]
#                 )
#               } else if (startsWith(f.lines.comment[x], "#Sites map ")) {
#                 reach.site.map <- TRUE
#               }
#             }
#           }
#         }
#       }
#
#       site.map <- do.call(rbind, site.map)
#
#       samples <- sort(unique(gt.ado.sample$sample))
#       gt.ado.sample <- gt.ado.sample %>%
#         filter(sample > samples[as.integer(ceiling(length(samples) * burn.in))]) %>%
#         left_join(cell.names, by = "psedu_cell") %>%
#         left_join(site.map, by = "psedu_site") %>%
#         select(-psedu_cell, -psedu_site) %>%
#         rename(site = collapsed_site)
#
#       gt.ado.names <- c("gt", "ado")
#       ret <- pbsapply(
#         gt.ado.names,
#         function(x) {
#           excluded <- gt.ado.names[x != gt.ado.names]
#           gt.ado.sample %>%
#             select(-!!as.name(excluded)) %>%
#             add_count(cell, site, !!as.name(x), name = "count") %>%
#             select(-sample) %>%
#             unique() %>%
#             group_by(cell, site) %>%
#             mutate(freq = count / sum(count)) %>%
#             ungroup() %>%
#             arrange(chr, pos, cell, !!as.name(x))
#         },
#         simplify = FALSE,
#         USE.NAMES = TRUE,
#         cl = 2L
#       )
#     },
#     finally = close(.gt.ado.sample.file)
#   )
#
#   if (!is.null(local.copy))
#     saveRDS(
#       object = ret,
#       file = local.copy
#     )
#
#   .end <- Sys.time()
#   message(paste0("Finished in ", (.end - .start), " seconds."))
#
#   return(ret)
# }
#
#
# .get.site.map <- function(x) {
#   ret <- str_split(
#     string = sub(
#       pattern = "^#",
#       replacement = "",
#       x = x
#     ),
#     pattern = "->",
#     simplify = TRUE
#   )
#
#   stopifnot(dim(ret)[2L] == 2L)
#
#   colnames(ret) <- c("psedu", "real")
#   ret <- as.data.frame(ret) %>%
#     mutate(p_real = .parse.real.site(real))
#
#   site.names <- ret$p_real
#   names(site.names) <- ret$psedu
#   return(site.names)
# }
#
#
# .get.site.map2 <- function(x) {
#   ret <- str_split(
#     string = sub(
#       pattern = "^#",
#       replacement = "",
#       x = x
#     ),
#     pattern = "->",
#     simplify = TRUE
#   )
#
#   stopifnot(dim(ret)[2L] == 2L)
#
#   colnames(ret) <- c("psedu_site", "ori_site")
#   ret <- as.data.frame(ret) %>%
#     left_join(.parse.real.site2(.$ori_site), by = "ori_site") %>%
#     select(-ori_site)
#   return(ret)
# }
#
#
# .parse.real.site <- function(x) {
#   ret <- str_split(
#     string = x,
#     pattern = ",",
#     simplify = TRUE
#   )
#
#   stopifnot(dim(ret)[2L] == 4L)
#
#   colnames(ret) <- c("chr", "pos", "ref", "alt")
#   ret <- as.data.frame(ret) %>%
#     mutate(collapsed = paste0(chr, ":", pos, "_", ref, "/", alt))
#   return(ret$collapsed)
# }
#
#
# .parse.real.site2 <- function(x) {
#   ret <- str_split(
#     string = x,
#     pattern = ",",
#     simplify = TRUE
#   )
#
#   stopifnot(dim(ret)[2L] == 4L)
#
#   colnames(ret) <- c("chr", "pos", "ref", "alt")
#   ret <- as.data.frame(ret) %>%
#     mutate(collapsed_site = paste0(chr, ":", pos, "_", ref, "/", alt)) %>%
#     mutate(ori_site = x)
#
#   ret$pos <- as.integer(ret$pos)
#
#   return(ret)
# }


loadVcf <- function(
  file,
  genome = "hg19",
  vcf.geno = c("GT", "GQ", "ADO", "ADOQ", "DP"),
  gt = "GT",
  rename.cells = NULL,
  included.cells = NULL,
  included.sites = NULL,
  destination.file = NULL,
  collapse.deletions = TRUE,
  cores = 4L,
  local.copy = NULL
) {
  stopifnot((!is.null(included.sites) && !is.null(destination.file)) || (is.null(included.sites) && is.null(destination.file)))

  # Filter sites if necessary.
  if (!is.null(destination.file) && !file.exists(destination.file)) {
    message("Filtering VCF file...")

    filterVcf(
      file = TabixFile(file),
      genome = genome,
      destination = destination.file,
      verbose = TRUE,
      prefilters = FilterRules(
        list(
          sites = function(x) {
            pbsapply(
              sub("\t", ":", sub("\t$", "", str_extract(x, ".*?\t\\d+?\t"))),
              FUN = function(y) sum(grepl(paste0("^", y, "($|\\D+)"), included.sites)) > 0,
              cl = cores
            )
          }
        )
      ),
      param = ScanVcfParam(
        geno = vcf.geno
      )
    )
    message("VCF file filtered.")
  }

  # Load VCF file.
  message("Loading VCF file...")
  vcf <- readVcf(
    file = ifelse(!is.null(destination.file) && file.exists(destination.file), destination.file, file),
    genome = genome,
    param = ScanVcfParam(
      geno = vcf.geno
    )
  )
  message("VCF file loaded.")

  # Rename cells if necessary.
  if (!is.null(rename.cells)) {
    colnames(vcf) <- vapply(
      colnames(vcf),
      FUN = function(x) {
        for (y in rename.cells) {
          x <- sub(y[1], y[2], x)
        }
        x
      },
      FUN.VALUE = character(1L),
      USE.NAMES = FALSE
    )
  }

  # Some cells should be filtered out if necessary.
  if (!is.null(included.cells)) {
    kept.cells <- colnames(vcf) %in% included.cells

    vcf@assays@data@listData <- lapply(
      seq_len(length(vcf.geno)),
      function(x) geno(vcf)[[x]][, kept.cells]
    )
    names(vcf@assays@data@listData) <- vcf.geno

    vcf@colData@rownames <- vcf@colData@rownames[kept.cells]
    vcf@colData@listData$Samples <- as.integer(seq_len(length(vcf@colData@rownames)))
    vcf@colData@nrows <- as.integer(length(vcf@colData@rownames))
  }

  # Convert alphabet genotypes to ternary counterparts.
  message("Converting genotypes from alphabetic to ternary...")
  vcf@assays@data@listData[["GTT"]] <- pbapply(
    geno(vcf)[[gt]],
    MARGIN = 2L,
    function(x) {
      sapply(
        x,
        convert.alpha2ternary,
        collapse.deletions = collapse.deletions,
        USE.NAMES = FALSE
      )
    },
    cl = cores
  )
  message("Genotype converted.")

  # Save VCF file if cell names are renamed or filtered.
  if (!is.null(destination.file)) {
    destination.file.2 <- file.path(
      dirname(destination.file),
      sub("(.*)\\.(.*)$", "\\1_cells.\\2", basename(destination.file))
    )

    if (!file.exists(destination.file.2) && (!is.null(rename.cells) || !is.null(included.cells))) {
      message("Saving VCF file...")
      writeVcf(
        vcf,
        destination.file.2
      )
      message("VCF file saved.")
    }
  }

  if (!is.null(local.copy))
    saveRDS(
      object = vcf,
      file = local.copy
    )

  return(vcf)
}


.get.chr <- function(x) {
  return(as.character(str_split(
    string = x,
    pattern = ":",
    simplify = TRUE
  )[, 1L]))
}


.get.pos <- function(x) {
  return(as.integer(str_split(
    string = str_split(
      string = x,
      pattern = ":",
      simplify = TRUE
    )[, 2L],
    pattern = "_",
    simplify = TRUE
  )[, 1L]))
}


`%+=%` <- function(e1, e2) eval.parent(substitute(e1 <- e1 + e2))


get.ternary <- function(data) {
  return(geno(data)[["GTT"]])
}


get.ado.state <- function(data) {
  return(geno(data)[["ADO"]])
}


convert2sifit <- function(
  vcf,
  prefix = '.',
  binary = TRUE
) {
  dir.create(prefix, showWarnings = FALSE)

  file.cell.names <- file.path(prefix, "cell_names")
  file.data <- file.path(prefix, "data")

  # Get the list of ternary genotypes.
  genotypes <- get.ternary(vcf)

  # Get the list of site names.
  sites <- rownames(genotypes)

  # Get the list of cell names.
  cells <- colnames(genotypes)

  writeLines(
    text = cells,
    con = file.cell.names,
    sep = " "
  )

  con <- file(file.data, open = "w")
  lapply(
    seq_along(sites),
    function(x) {
      writeLines(
        text = c(
          sites[x],
          .convert2sifit.per.site(genotypes[x,], binary = binary)
        ),
        con = con,
        sep = " "
      )
      writeLines(
        text = "\n",
        con = con,
        sep = ""
      )
    }
  )
  close(con)
}


.convert2sifit.per.site <- function(
  vals,
  binary = TRUE
) {
  vals[vals > 2] <- 2
  vals[vals == -3] <- 3
  vals[vals == -2] <- 2
  vals[vals == -1] <- 0

  if (binary) {
    vals[vals == 2] <- 1
  }

  return(vals)
}


loadGT <- function(
  file,
  prefix_cells = "^",
  suffix_cells = "$",
  prefix_filter_cells = "",
  filter_cells = NULL,
  suffix_filter_cells = "",
  filter_positions = NULL,
  cores = 4L
) {
  ret <- list()

  print("Loading file...")
  genotypes <- readGT(file)
  print("File loaded.")

  # rename rows
  pos.names <- rownames(genotypes)
  rownames(genotypes) <- vapply(
    pos.names,
    FUN = function(x) {
      tmp <- strsplit(x, "_")
      stopifnot(length(tmp) == 1)
      return(tmp[[1]][1])
    },
    FUN.VALUE = character(1),
    USE.NAMES = FALSE
  )

  # rename columns
  cell.names <- colnames(genotypes)
  colnames(genotypes) <- vapply(
    cell.names,
    FUN = function (x) {
      y <- str_match(x, paste0(prefix_cells, "(.*?)", suffix_cells))
      stopifnot(length(y) == 2)
      return(y[2])
    },
    FUN.VALUE = character(1),
    USE.NAMES = FALSE
  )

  # filter cells if required
  cell.indices <- 1:dim(genotypes)[2]
  if (!is.null(filter_cells)) {
    cell.names.ori <- colnames(genotypes)

    cell.indices <- vapply(
      filter_cells,
      FUN = function(x) {
        grep(paste0(prefix_filter_cells, "(", x, ")", suffix_filter_cells), cell.names.ori, value = FALSE)
      },
      FUN.VALUE = integer(1),
      USE.NAMES = FALSE
    )
  }

  # filter positions if required
  pos <- list()
  pos$index <- 1:dim(genotypes)[1]
  pos$missing.pos <- character()
  if (!is.null(filter_positions)) {
    pos.names.ori <- rownames(genotypes)

    pos.indices <- lapply(
      filter_positions,
      FUN = function(x) {
        tmp <- match(x, pos.names.ori)

        return(list(tmp, if_else(is.na(tmp), x, NA_character_)))
      }
    )

    pos$index <- sapply(pos.indices, "[[", 1, USE.NAMES = FALSE)
    pos$missing.pos <- sapply(pos.indices, "[[", 2, USE.NAMES = FALSE)
  }

  # save data
  ret$gt.alpha <- genotypes[pos$index[!is.na(pos$index)], cell.indices]
  if (!is.null(filter_cells))
    colnames(ret$gt.alpha) <- filter_cells

  ret$missing.pos <- pos$missing.pos[!is.na(pos$missing.pos)]

  # convert alphabet genotypes to ternary counterparts
  print("Converting genotypes from alphabetic to ternary...")
  plan(multicore, workers = cores)
  ternary <- future_apply(
    ret$gt.alpha,
    MARGIN = 2L,
    function(x) {
      sapply(
        x,
        convert.alpha2ternary,
        USE.NAMES = FALSE
      )
    }
  )
  names(ternary) <- names(ret$gt.alpha)
  print("Genotype converted.")

  ret$gt.ternary <- ternary

  return(ret)
}


#` Convert alphabet genotypes to ternary counterparts.
#` -3: missing genotype
#` -2: 1/. when collapse.deletions = FALSE
#` -1: 0/. when collapse.deletions = FALSE
#` 0: homozygous reference (0/0 or 0/. when collapse.deletions = TRUE)
#` 1: heterozygous mutation (0/1)
#` 2: homozygous mutation with homozygous alternative nucletides (1/1 or 1/. when collapse.deletions = TRUE)
#` 3: homozygous mutation with heterozygous alternative nucleotides (1/1')
convert.alpha2ternary <- function(
  gt,
  collapse.deletions = TRUE
) {
  splitted <- sort(
    strsplit(
      gt,
      if_else(
        grepl("/", gt, fixed = TRUE),
        "/",
        "|"),
      fixed = TRUE
    )[[1L]]
  )

  stopifnot(length(splitted) == 2L)

  no.digit <- splitted[grepl("\\D", splitted)]

  if (length(no.digit) == 2L)
    return(-3L)
  else if (length(no.digit) == 1L) {
    digit <- splitted[grepl("\\d", splitted)]

    if (as.integer(digit) == 0L)
      return(if_else(collapse.deletions, 0L, -1L))
    else
      return(if_else(collapse.deletions, 2L, -2L))
  } else {
    splitted <- as.integer(splitted)

    zero <- splitted[splitted == 0L]

    if (length(zero) == 2L)
      return(0L)
    else if (length(zero) == 1L)
      return(1L)
    else {
      if (length(levels(as.factor(splitted))) == 1L)
        return(2L)
      else
        return(
          if_else(
            sum(splitted) <= 2L,
            2L,
            3L
          )
        )
    }
  }

  return(NA)
}


filter.invariant.sites <- function(gt, threshold = 1) {
  rownames(gt)[apply(
    gt,
    MARGIN = 1L,
    FUN = function(x) {
      y <- x[x != 0]

      return(if_else(length(y) <= threshold, TRUE, FALSE))
    }
  )]
}


get.ordered.intersactions <- function(ordered.ref, target) {
  stopifnot(!is.null(ordered.ref) && !is.null(target))

  ref.df <- tibble(ori = ordered.ref) %>%
    separate(
      ori,
      c("chr", "pos", "ref", "alt"),
      ":|_|/",
      convert = TRUE
    ) %>%
    select(-ref, -alt) %>%
    mutate(ref_id = 1:n())

  target.df <- tibble(ori = target) %>%
    separate(
      ori,
      c("chr", "pos", "ref", "alt"),
      ":|_|/",
      remove = FALSE,
      convert = TRUE
    ) %>%
    mutate(target_id = 1:n()) %>%
    select(-ref, -alt)

  ret <- ref.df %>%
    inner_join(
      target.df,
      by = c("chr", "pos")
    ) %>%
    arrange(ref_id) %>%
    select(ori, target_id) %>%
    rename(pos = ori, id = target_id)

  ret
}


get.gt.hm.color <- function(i) {
  if (i == 1)
    colors <- c(
      wes_palette("Moonrise3")[1],
      "white",
      wes_palette("Moonrise3")[2],
      wes_palette("Moonrise3")[5],
      wes_palette("Moonrise3")[3]
    )
  else if (i == 2)
    colors <- c(
      wes_palette("Darjeeling1")[2], #("Cavalcanti1")[1],
      "white",
      wes_palette("GrandBudapest2")[4],
      wes_palette("GrandBudapest1")[2],
      wes_palette("GrandBudapest1")[3]
    )
  else if (i == 3)
    colors <- c(
          "#B09C85FF",
          "white",
          "#3C5488FF",
          "#E64B35FF",
          "#4DBBD5FF"
        )

  names(colors) <- c("-3", "0", "1", "2", "3")
  colors
}


get.legend <- function(
  data,
  is.ado = FALSE
) {
  if (is.ado) {
    category <-  c("No ADO", "Having ADO")
    names(category) <- c("0", "1")
  } else {
    category <- c("Missing", "0/0", "0/1", "1/1", "1/1'")
    names(category) <- c("-3", "0", "1", "2", "3")
  }

  vals <- as.character(sort(unique(as.vector(as.matrix(data)))))
  return(category[vals])
}


get.genotype.stats <- function(data, is.sieve) {
  genotypes <- get.ternary(data)

  prop <- sapply(
    c(-3, 0, 1, 2, 3),
    function(x) sum(genotypes == x)
  ) / prod(dim(genotypes))

  names(prop) <- c("Missing", "0/0", "0/1", "1/1", "1/1'")

  if (is.sieve) {
    prop["Missing"] <- NA
  } else {
    prop["1/1'"] <- NA
  }

  prop
}


get.formatted.num.as.char <- function(val) {
  .val <- str_split(
    string = as.character(val),
    pattern = "\\.",
    simplify = TRUE
  )

  int <- sapply(
    seq_len(dim(.val)[1]),
    function(x) {
      str_split(
        string = .val[x, 1],
        pattern = ""
      )
    }
  )

  ret <- sapply(
    seq_len(length(int)),
    function(x) {
      dot.num <- ceiling(length(int[[x]]) / 3) - 1L

      .ret <- character()
      i <- length(int[[x]])
      j <- 1L
      while (dot.num > 0L || i > 0L) {
        .ret <- c(.ret, int[[x]][i])
        i <- i - 1L

        if (j == 3) {
          j <- 1L

          if (dot.num > 0L) {
            .ret <- c(.ret, ",")
            dot.num <- dot.num - 1L
          }
        } else {
          j <- j + 1L
        }
      }

      paste0(rev(.ret), collapse = "")
    }
  )

  if (dim(.val)[2L] > 1L) {
    ret <- tibble(
      int = ret,
      dec = .val[, 2L]
    ) %>%
      mutate(
        full = if_else(
          dec != "",
          paste(int, dec, sep = "."),
          int
        )
      )

    return(ret$full)
  } else {
    return(ret)
  }
}

