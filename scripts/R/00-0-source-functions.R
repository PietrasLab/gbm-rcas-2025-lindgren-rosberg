.palette_discrete_bp <- function(discrete.vec){
  require(viridis)
  # For barplots
  # x is the full vector - unique entries from these are taken
  # to generate a palette
  x <- discrete.vec
  colourCount = length(unique(x))
  my.pal <- rev(viridis::rocket(colourCount+1)[-1])
  names(my.pal) <- unique(x)
  return(my.pal)
}


#' stacked barplot (identity)
#' @description wrapper for plotting data using ggplot geom_barplot
#' @family ggplot
#' @param x a character vector of values
#' @param random if shuffle/randomize colors
#' @return named vector with color palette of lengh unique x
.colorPalDiscrete <- function(x, random=FALSE){
  # for scatterplots
  # x is the full vector - unique entries from these are taken
  # to generate a palette
  colourCount = length(unique(x))
  if(colourCount < 3){
    my.pal <- (viridis::rocket(colourCount+1))
  }else if(colourCount < 9){
    my.pal <- rev(viridis::rocket(colourCount+1))
    # }else if (colourCount < 9 ){
    #   my.pal <- RColorBrewer::brewer.pal(9, "Set1")[-9  ]
  }else{
    getPalette = grDevices::colorRampPalette(RColorBrewer::brewer.pal(9, "Set1")[-9])
    # getPalette = grDevices::colorRampPalette(
    #   wesanderson::wes_palette(name = "Darjeeling1", n = 5)
    #   )
    my.pal <- getPalette(colourCount)
  }
  if(random){
    my.pal <- sample(my.pal, size = colourCount, replace = FALSE)
  }
  names(my.pal) <- unique(x)
  return(my.pal)
}


.barplot_stacked <- function(
    plot_df,
  my.pal = "",
  group.var = "",
  color.var = "",
  wrap.var = "",
  scaled.y=FALSE,
  scaled.y.wihtin.color = FALSE,
  swap.xy=FALSE,
  return.df.only = FALSE
){
  if(swap.xy){
    tmp <- group.var
    group.var <- color.var
    color.var <- tmp
    rm(tmp)
  }

  if( wrap.var %in% c("","--none--")){
    # if no wrap var
    ylab <- "cell count"
    df <- plot_df[,c(group.var, color.var)]
    colnames(df) <- c("group.var","color.var")
    df_sum <- df %>%
      dplyr::group_by(group.var) %>%
      summarize(group_sum = n())
    df_sum_color <- df %>%
      group_by(color.var) %>%
      summarize(color_var_sum = n())
    df <- df %>%
      dplyr::group_by(group.var, color.var) %>%
      summarize(y = n()) %>%
      dplyr::left_join(df_sum) %>%
      mutate(y_scaled = y/group_sum) %>%
      dplyr::left_join(df_sum_color) %>%
      mutate(y_scaled_within_color = y/color_var_sum)
    df_scaled_all <- df %>%
      dplyr::group_by(group.var) %>%
      summarize(scaled_color_sum = sum(y_scaled_within_color))
    df <- df %>%
      dplyr::left_join(df_scaled_all) %>%
      mutate(y_scaled_all = y_scaled_within_color/scaled_color_sum)

  }else{
    # include the wrap var
    ylab <- "cell fraction relative group"
    df <- plot_df[,c(group.var, color.var, wrap.var)]
    colnames(df) <- c("group.var","color.var","wrap.var")
    df_sum <- df %>%
      dplyr::group_by(group.var) %>%
      summarize(group_sum = n())
    df_sum_color <- df %>%
      group_by(color.var) %>%
      summarize(color_var_sum = n())
    df <- df %>%
      dplyr::group_by(group.var, color.var, wrap.var) %>%
      summarize(y = n()) %>%
      dplyr::left_join(df_sum) %>%
      mutate(y_scaled = y/group_sum) %>%
      dplyr::left_join(df_sum_color) %>%
      mutate(y_scaled_within_color = y/color_var_sum)
  }
  if (scaled.y && !scaled.y.wihtin.color){
    df <- df %>%
      mutate(y = y_scaled)
  } else if (!scaled.y && scaled.y.wihtin.color){
    df <- df %>%
      mutate(y = y_scaled_within_color)
  } else if (scaled.y && scaled.y.wihtin.color) {
    df <- df %>%
      mutate( y = y_scaled_all)
  }


  if(return.df.only){
    return(df)
  }

  if(my.pal[1] == ""){
    my.pal <- .palette_discrete_bp(df$color.var)
  }
  g <- ggplot(
    df,
    aes(x = group.var,
      y = y,
      fill = color.var)) +
    geom_bar(position = "stack", stat = "identity", color="black")

  if(! wrap.var %in% c("","--none--")){
    g <- g +
      facet_wrap(~ wrap.var, scales=("free"), drop = TRUE)
  }

  g <- g +
    scale_fill_manual(color.var, values = my.pal) +
    theme(axis.text.x=element_text(angle = 45, hjust = 1)) +
    labs(x = group.var, y = ylab)


  return(g)
}



.violin <- function(
    plot_df,
  group.var = "",
  score.var = "",
  wrap.var = "",
  scaled_y=FALSE,
  log_y=FALSE,
  my.pal=""
){
  require(ggplot2)

  # ylab <- "score"
  if( wrap.var %in% c("","--none--")){
    # cat("no wrap label defined")0
    # if no wrap var
    df <- plot_df[,c(group.var, score.var)]
    colnames(df) <- c("group.var","score.var")
  }else{
    # include the wrap var
    # cat("found wrap label")
    df <- plot_df[,c(group.var, score.var, wrap.var)]
    colnames(df) <- c("group.var","score.var","wrap.var")
  }

  if(my.pal[1] == ""){
    my.pal <- .palette_discrete_bp(df$group.var)
  }

  g <- ggplot(
    df,
    aes(x = group.var,
      y = score.var,
      fill = group.var
    )) +
    geom_violin(
      alpha=0.5, width=1, trim = TRUE,
      scale = "width", adjust = 0.5, color="black") +
    geom_boxplot(
      width=0.2, outlier.colour="grey25",
      notch = FALSE, notchwidth = .75, alpha = 0.75,
      colour = "grey10", outlier.size=0.3) +
    scale_fill_manual(values = my.pal, aesthetics = ) +
    guides(fill="none")

  if(! wrap.var %in% c("","--none--")){
    # cat("adding wrap")
    g <- g +
      facet_wrap(~ wrap.var, scales=("free"), drop = TRUE)
  }

  g <- g +
    scale_fill_manual(score.var, values = my.pal) +
    theme(axis.text.x=element_text(angle = 45, hjust = 1)) +
    labs(x = group.var, y = score.var)

  if(log_y){
    g <- g +
      scale_y_continuous(trans='log10')
  }

  return(g)
}


.barplot_y <- function(
    plot_df,
  group.var = "",
  y.var = "",
  wrap.var = ""
){

  if( wrap.var %in% c("","--none--")){
    # if no wrap var
    # ylab <- "cell count"
    df <- plot_df[,c(group.var, y.var)]
    colnames(df) <- c("group.var","y.var")
  }else{
    # include the wrap var
    # ylab <- "cell fraction relative group"
    df <- plot_df[,c(group.var, y.var, wrap.var)]
    colnames(df) <- c("group.var","y.var","wrap.var")
  }

  my.pal <- .palette_discrete_bp(df$group.var)

  g <- ggplot(
    df,
    aes(x = group.var,
      y = y.var,
      fill = group.var)) +
    geom_bar(stat = "identity", color="black")

  if(! wrap.var %in% c("","--none--")){
    g <- g +
      facet_wrap(~ wrap.var, scales=("free"), drop = TRUE)
  }

  g <- g +
    scale_fill_manual(y.var, values = my.pal) +
    theme(axis.text.x=element_text(angle = 45, hjust = 1)) +
    labs(x = group.var, y = y.var)


  return(g)
}

setClass(
  Class = "seuratFilterObject",
  slots = list(
    filter.name="character",
    cell.filters = "list",
    features.keep = "character",
    features.remove = "character",
    min.cells.expressed = "numeric",
    assert.dim = "numeric"
  ))
.CreateSeuratFilter <- function(
    filter.name = vector(mode="character", length = 0L),
  min.cells.expressed = vector(mode="numeric", length = 0L),
  cell.filters=list(),
  features.keep=vector(mode="character", length = 0L),
  features.remove=vector(mode="character", length = 0L),
  assert.dim=vector(mode="numeric", length = 0L)
){
  new(
    Class = "seuratFilterObject",
    filter.name = filter.name,
    min.cells.expressed = min.cells.expressed,
    cell.filters = cell.filters,
    features.keep = features.keep,
    features.remove = features.remove,
    assert.dim=assert.dim
  )
}


.seuratFilterFoo <- function(
    seurat.object,
  filter.object = NULL,
  filter.cfg
){
  # new.env()
  cli::cli_h1(paste0(
    cli::symbol$star,
    " SeuratFilterFoo ",
    cli::symbol$star
  ))
  # cli::cli_h2(paste0(
  #   " Validate input "
  # ))


  #
  # Check input
  #
  cat(paste0(cli::symbol$arrow_right, " Validating input args"))
  cl <- TRUE
  start_time <- Sys.time()

  # check if seruat object
  if (!(class(seurat.object) == "Seurat")) {
    if (cl) {
      cat("\n")
      cl <- FALSE
    }
    cli::cli_abort(
      "{.var seurat.object} needs to be an object of class SeuratObject.")
  }

  # define input mode (if filter.object from file or direct input)
  # if (exists(x = "filter.object")) {
  if (!is.null(filter.object)) {
    read.filter.config = FALSE

    if ( class(filter.object) != "seuratFilterObject") {
      if (cl) {
        cat("\n")
        cl <- FALSE
      }
      cli::cli_abort(
        "{.var filter.object} must be a .CreateSeuratFilter
        supplied as `filter.object `or defined in the 'filter.cfg' file:")
    }
  } else {
    read.filter.config = TRUE
  }

  # check if config is file (if seurat.object not provided directly)
  if ( read.filter.config && !(file.exists(filter.cfg)) ) {
    if (cl) {
      cat("\n")
      cl <- FALSE
    }
    cli::cli_abort(
      " if 'filter.object' is not provided {.var filter.cfg}
      must be a config file sourcable by R in which the 'filter.object'
      is defined through .CreateSeuratFilter")
  }

  # end section
  end_time <- Sys.time()
  if (cl) {
    # ANSI code for clearing the line and reset cursor to the beginning
    # of the line
    cat("\33[2K\r")
  }
  cli::cli_alert_success(paste0(
    " Input args ok ",
    cli::col_blue(
      "(in ", prettyunits::pretty_dt(end_time - start_time), ")"
    )
  ))

  #
  # source the filter cfg
  #
  if (read.filter.config){
    cat(paste0(cli::symbol$arrow_right, " Sourcing filter cfg file "))
    cl <- TRUE
    start_time <- Sys.time()

    source(filter.cfg, local = TRUE)

    if ( class(filter.object) != "seuratFilterObject" ) {
      if (cl) {
        cat("\n")
        cl <- FALSE
      }
      cli::cli_abort(
        "{.var filter.object} must be defined in the config file:
      a .seuratFilterObject and named `filter.object`")
    }

    end_time <- Sys.time()
    if (cl) {
      # ANSI code for clearing the line and reset cursor to the beginning
      # of the line
      cat("\33[2K\r")
    }
    cli::cli_alert_success(paste0(
      " Sourced filter cfg file ",
      cli::col_blue(
        "(in ", prettyunits::pretty_dt(end_time - start_time), ")"
      )
    ))
  }

  #
  # check filter.object
  #
  cat(paste0(
    cli::symbol$arrow_right,
    " Checking supplied filter.object "
  ))
  cl <- TRUE
  start_time <- Sys.time()

  ## check assert.dim.slot

  if( length(filter.object@assert.dim) ){
    do_assert_dims = TRUE
    if( !is.numeric(filter.object@assert.dim) &&
        length(filter.object@assert.dim) != 2 ){
      if (cl) {
        cat("\n")
        cl <- FALSE
      }
      cli::cli_abort(
        "{.var filter.object@assert.dim} must be numeric vector
        of length 2 representing the number of genes and cells,
        respectively, after filtering.")
    }
  } else {
    do_assert_dims = FALSE
  }

  # check if min.cells.expressed
  if( length(filter.object@min.cells.expressed) > 0 ){
    do_filter_min_expressed = TRUE
    if( !is.numeric(filter.object@min.cells.expressed) &&
        length(filter.object@min.cells.expressed) != 1 ){
      if (cl) {
        cat("\n")
        cl <- FALSE
      }
      cli::cli_abort(
        "{.var filter.object@min.cells.expressed} must be numeric vector
        of length 1 representing the minimum number of cells a feature must be expressed in .")
    }
  } else {
    do_filter_min_expressed = FALSE
  }

  ## Check cell.filters slot
  if ( length(filter.object@cell.filters) >0 ){
    do_cell_filter = TRUE
    cell.filters <- filter.object@cell.filters
    if ( !all(unlist(lapply(cell.filters, length))==1) ){
      if (cl) {
        cat("\n")
        cl <- FALSE
      }
      cli::cli_abort(
        "{.var filter.object@cell.filters} must be a list of
        named character vectors each defining a filter operation
        on a metadata.colum")
    }
  } else {
    do_cell_filter = FALSE
  }

  # wrap up this section
  end_time <- Sys.time()
  if (cl) {
    # ANSI code for clearing the line and reset cursor to the beginning
    # of the line
    cat("\33[2K\r")
  }
  cli::cli_alert_success(paste0(
    " Sourced file filter.object ok. ",
    cli::col_blue(
      "(in ", prettyunits::pretty_dt(end_time - start_time), ")"
    )
  ))

  #
  # Info Input object
  #
  # cat("\n")
  # cli::cli_h2(paste0(
  #   " Input SeuratObject & Filters"
  # ))
  cli::cli_alert_info("filter.name: {.var {filter.object@filter.name}}")
  cat(paste0(
    cli::col_yellow("-"),
    cli::col_yellow(cli::symbol$pointer),
    cli::col_white("  A Seurat object with "),
    cli::col_yellow(ncol(seurat.object)," cells"),
    cli::col_white(" and "),
    cli::col_yellow(nrow(seurat.object)," genes")
  ))
  cat("\n")

  ## if. do_filter_min_expressed
  if (do_filter_min_expressed) {
    cat(paste0(
      cli::col_yellow(cli::symbol$info),
      cli::col_white("  `min.cells.expressed` slot: "),
      cli::col_yellow((filter.object@min.cells.expressed))
    ))
    cat("\n")
  } else {
    cat(paste0(
      cli::col_yellow(cli::symbol$warning),
      cli::col_white("  No 'min.cells.expressed' supplied: skipping filter")
    ))
    cat("\n")
  }


  ## if do_cell_filter
  if (do_cell_filter) {
    cat(paste0(
      cli::col_yellow(cli::symbol$info),
      cli::col_white("  `cell.filters` slot: "),
      cli::col_yellow(length(filter.object@cell.filters)),
      cli::col_white(" filter(s) supplied. ")
    ))
    cat("\n")
  } else {
    cat(paste0(
      cli::col_yellow(cli::symbol$warning),
      cli::col_white("  No cell filter supplied: skipping filter")
    ))
    cat("\n")
  }

  ## Check features.keep and features.remove
  if (length(filter.object@features.keep) > 0){
    do_filter_features_keep = TRUE
    features.keep <- filter.object@features.keep
    cat(paste0(
      cli::col_green(cli::symbol$info),
      cli::col_white("  `features.keep` slot:"),
      cli::col_yellow(length((features.keep))),
      cli::col_white(" supplied. "),
      cli::col_yellow(length(sort(unique(features.keep)))),
      cli::col_white(" unique. ")
    ))
    cat("\n")
    features.keep <- base::intersect(
      sort(unique(features.keep)),
      rownames(seurat.object))
    if( length(features.keep)==0 ){
      do_filter_features_keep = FALSE
      features.keep <- c()
      cat(paste0(
        cli::col_yellow(cli::symbol$warning),
        cli::col_white("  None of the supplied features are present in your seurat.object. skipping filter")
      ))
      cat("\n")
    }
  } else {
    do_filter_features_keep = FALSE
    features.keep <- c()
    cat(paste0(
      cli::col_yellow(cli::symbol$warning),
      cli::col_white("  `features.keep` not supplied: Will not filter out any genes.")
    ))
    cat("\n")
  }

  ## Check features.remove
  if (length(filter.object@features.remove) > 0){
    do_filter_features_remove = TRUE
    features.remove <- filter.object@features.remove
    cat(paste0(
      cli::col_green(cli::symbol$info),
      cli::col_white("  `features.remove` slot: "),
      cli::col_yellow(length((features.remove))),
      cli::col_white(" supplied. "),
      cli::col_yellow(length(sort(unique(features.remove)))),
      cli::col_white(" unique. ")
    ))
    cat("\n")
    features.remove <- base::intersect(
      sort(unique(features.remove)),
      rownames(seurat.object))

    if( length(features.remove)==0 ){
      do_filter_features_remove = FALSE
      features.remove <- c()
      cat(paste0(
        cli::col_yellow(cli::symbol$warning),
        cli::col_white("  None of the supplied features are present
          in your seurat.object. skipping filter")
      ))
      cat("\n")
    }
  } else {
    do_filter_features_remove = FALSE
    features.remove <- c()
    cat(paste0(
      cli::col_yellow(cli::symbol$warning),
      cli::col_white("  `features.remove` not supplied. Will not remove any genes.")
    ))
    cat("\n")
  }

  #
  #  Cell Metadata Filters Seurat::WhichCells
  #
  # cat("\n")
  if ( do_cell_filter ){
    # cli::cli_h2(paste0(
    #   " Cell Metadata Filters "
    # ))

    ## loop all filters
    args.list <- filter.object@cell.filters
    my.df <- seurat.object@meta.data
    cat(paste0(
      cli::col_yellow(""),
      cli::col_green(cli::symbol$info),
      cli::col_white("  Looping filters - apply each on meta.data table ")
    ))
    cat("\n")
    for (i in 1:length(args.list)) {
      filter.arg <- args.list[i]
      filter.name <- names(args.list)[i]

      cat(paste0(
        cli::col_yellow(""),
        cli::col_green(cli::symbol$circle),
        cli::col_white("  Filter "),
        cli::col_yellow(filter.name),
        cli::col_white(" - "),
        cli::col_green(filter.arg)
        #cli::col_white(": ")
      ))
      cat("\n")
      cells.before <- nrow(my.df)
      y <- paste0(
        "my.df %>% ",
        "dplyr::filter(",filter.arg,")"
      )
      my.df <- eval(parse(text = y))
      cat(paste0(
        cli::col_yellow(""),
        cli::col_yellow(cli::symbol$pointer),
        cli::col_white("  Metadata - rows removed: "),
        cli::col_yellow(cells.before - nrow(my.df)),
        cli::col_white(". kept: "),
        cli::col_yellow(nrow(my.df))
        #cli::col_white(": ")
      ))
      cat("\n")
    } # end loop all whichCells

    # Apply cell filter based on metadata rows left after filter loop
    cat(paste0(
      cli::col_yellow(cli::symbol$circle),
      cli::col_white("  Filtering seurat.object using filtered metadata table")
      #cli::col_yellow("'cell.filters' on seurat.object")
    ))
    cl <- TRUE
    start_time <- Sys.time()


    selected.c <- rownames(my.df)
    seurat.object <- subset(
      seurat.object,
      cells = selected.c
    )
    cat("\33[2K\r")
    cat(paste0(
      cli::col_green(cli::symbol$tick)
    ))
    cat("\n")

    end_time <- Sys.time()
    cli::cli_alert_success(paste0(
      " Done looping all 'cell.filters' ",
      cli::col_blue(
        "(in ", prettyunits::pretty_dt(end_time - start_time), ")"
      )
    ))

  }

  #
  #  filter_min_expressed
  #
  cat("\n")
  if (do_filter_min_expressed) {
    # cli::cli_h2(paste0(
    #   " FeatureFilter: Min cells expressed "
    # ))

    cl <- TRUE
    start_time <- Sys.time()

    selected.f <- rownames(seurat.object)[
      Matrix::rowSums(Seurat::GetAssayData(seurat.object, layer = 'counts'))  >
        filter.object@min.cells.expressed ]
    # str(selected.f)

    cat(paste0(
      cli::col_yellow(cli::symbol$circle),
      cli::col_white("  Filter "),
      cli::col_yellow("'min.cells.expressed' "),
      cli::col_white("- keeping: "),
      cli::col_yellow(length(selected.f)),
      cli::col_white(" features. Removing: "),
      cli::col_yellow(nrow(seurat.object) - length(selected.f))
    ))
    cat("\n")

    cat(paste0(
      cli::symbol$arrow_right,
      "  Subsetting seurat object "
    ))
    seurat.object <- subset(
      seurat.object,
      features = selected.f
    )

    cat("\33[2K\r")
    cat(paste0(
      cli::col_green(cli::symbol$tick)
    ))
    cat("\n")

    end_time <- Sys.time()
    cli::cli_alert_success(paste0(
      " Done looping all 'cell.filters' ",
      cli::col_blue(
        "(in ", prettyunits::pretty_dt(end_time - start_time), ")"
      )
    ))
  }

  #
  # features.keep filter
  #
  # update the vec
  features.keep <- base::intersect(
    sort(unique(features.keep)),
    rownames(seurat.object))

  if ( do_filter_features_keep &&
      length(features.keep) > 0
  ){
    cli::cli_h2(paste0(
      " Feature filter - features.keep"
    ))

    # update the vec to filtered object
    cl <- TRUE
    start_time <- Sys.time()

    cat(paste0(
      cli::col_yellow(cli::symbol$circle),
      cli::col_white("  Filter "),
      cli::col_yellow("'features.keep' "),
      cli::col_white("- keeping: "),
      cli::col_yellow(length(features.keep)),
      cli::col_white(" features")
    ))

    seurat.object <- subset(
      seurat.object,
      features = features.keep,
      invert=FALSE
    )

    cat("\33[2K\r")
    cat(paste0(
      cli::col_green(cli::symbol$tick)
    ))
    cat("\n")

    end_time <- Sys.time()
    cli::cli_alert_success(paste0(
      " done 'features.keep' filter ",
      cli::col_blue(
        "(in ", prettyunits::pretty_dt(end_time - start_time), ")"
      )
    ))
  } else {
    cat(paste0(
      cli::col_yellow(cli::symbol$warning),
      cli::col_white("  No features left to filter - skipping 'features.keep'")
    ))
    cat("\n")
  } # end if filter keep

  #
  # features.remove filter
  #

  # update the remove vec
  features.remove <- base::intersect(
    sort(unique(features.remove)),
    rownames(seurat.object))

  if ( do_filter_features_remove &&
      length(features.remove) > 0
  ){
    cli::cli_h2(paste0(
      " Gene List filters - features.remove"
    ))

    cl <- TRUE
    start_time <- Sys.time()

    cat(paste0(
      cli::col_yellow(cli::symbol$circle),
      cli::col_white("  Filter "),
      cli::col_yellow("'features.remove' "),
      cli::col_white("- removing: "),
      cli::col_yellow(length(features.remove)),
      cli::col_white(" features")
    ))

    u <- match(features.remove, rownames(seurat.object))
    seurat.object <- seurat.object[-u,]

    cat("\33[2K\r")
    cat(paste0(
      cli::col_green(cli::symbol$tick)
    ))
    # cat("\n")
    #
    # cat(paste0(
    #   cli::col_yellow("-"),
    #   cli::col_yellow(cli::symbol$pointer),
    #   cli::col_white("  A Seurat object with "),
    #   cli::col_yellow(ncol(seurat.object)," cells"),
    #   cli::col_white(" and "),
    #   cli::col_yellow(nrow(seurat.object)," genes")
    # ))

    cat("\n")
    end_time <- Sys.time()
    cli::cli_alert_success(paste0(
      " done features.remove filter ",
      cli::col_blue(
        "(in ", prettyunits::pretty_dt(end_time - start_time), ")"
      )
    ))

  } else {
    cat(paste0(
      cli::col_yellow(cli::symbol$warning),
      cli::col_white("  No features left to filter - skipping 'features.remove'")
    ))
    cat("\n")
  } # end if filter remove



  ## IF ASSERT DIMENSIONS
  if (do_assert_dims) {
    cli::cli_h2(paste0(
      " Assert dimensions of filtered Object "
    ))

    cat(paste0(cli::symbol$arrow_right, " Check dims ..."))
    cl <- TRUE
    start_time <- Sys.time()

    # check if seruat object
    if ( !(identical(
      as.numeric(dim(seurat.object)),
      filter.object@assert.dim ))) {
      if (cl) {
        cat("\n")
        cl <- FALSE
      }

      cat(paste0(
        cli::col_green("-"),
        cli::col_green(cli::symbol$pointer),
        cli::col_white("     "),
        cli::col_yellow(nrow(seurat.object)," genes"),
        cli::col_white(" and "),
        cli::col_yellow(ncol(seurat.object)," cells")
      ))
      cat("\n")
      cli::cli_abort(paste0(
        cli::col_white(" Dimensions supplied in 'assert.dim' slot does not match "),
        cli::col_white(" demessions of filtered object. You supplied: "),
        cli::col_yellow(filter.object@assert.dim[1]," genes"),
        cli::col_white(" and "),
        cli::col_yellow(filter.object@assert.dim[2]," cells")
      ))
    }

    cat("\33[2K\r")
    cat(paste0(
      cli::col_green(cli::symbol$tick)
    ))
    cat("\n")
    end_time <- Sys.time()
    cli::cli_alert_success(paste0(
      " Done Assert Dimensions of object ",
      cli::col_blue(
        "(in ", prettyunits::pretty_dt(end_time - start_time), ")"
      )
    ))

  }

  ## RETURN

  cat(paste0(
    cli::col_green("---"),
    cli::col_green(cli::symbol$pointer),
    cli::col_white("  "),
    cli::col_yellow(nrow(seurat.object)," genes"),
    cli::col_white(" and "),
    cli::col_yellow(ncol(seurat.object)," cells")
  ))
  cat("\n")
  cat("\n")
  cli::cli_h1(paste0(
    #cli::symbol$star,
    " done "
  ))
  return(seurat.object)
}


.seuratFactorizeMdata <-
  function(
      seurat.object
  ){
    mdata <- seurat.object@meta.data
    for(mcol in colnames(mdata)){ # mcol = "nCount_RNA_bin"
      if (mcol %in% names(.color_pal)) {
        cli::cli_alert(". Setting factor levels for {.var {mcol}}")
        if (! all(
          unique(mdata[[mcol]])[!is.na(unique(mdata[[mcol]]))] %in%
            names(.color_pal[[mcol]])
        )){
          cli::cli_alert_warning("Not all levels of {.var {mcol}} in .color_pal, skipping!")
          next(mcol)
        }
        seurat.object[[mcol]] <- factor(
          mdata[[mcol]],
          levels = names(.color_pal[[mcol]]))
        # table(mdata[[mcol]])
        seurat.object[[mcol]] <- droplevels(seurat.object[[mcol]])
      }}

    return( seurat.object )
  }


.FactorizeMdata <- function(
    mdata
){
  for(mcol in colnames(mdata)){ # mcol = "nCount_RNA_bin"
    if (mcol %in% names(.color_pal)) {
      cli::cli_alert(". Setting factor levels for {.var {mcol}}")
      if (! all(
        unique(mdata[[mcol]])[!is.na(unique(mdata[[mcol]]))] %in%
          names(.color_pal[[mcol]])
      )){
        cli::cli_alert_warning("Not all levels of {.var {mcol}} in .color_pal, skipping!")
        next(mcol)
      }
      mdata[[mcol]] <- factor(
        mdata[[mcol]],
        levels = names(.color_pal[[mcol]]))
      # table(mdata[[mcol]])
      mdata[[mcol]] <- droplevels(mdata[[mcol]])
    }}
  return( mdata )
}

.seuratFeaturePlotHexbin <-
  function(
      seurat.object,
    reduction.name = "",
    feature,
    n.bins = 250
  ){
    library(Seurat)
    library(ggplot2)
    library(dplyr)
    library(cli)
    library(glue)

    if(reduction.name == ""){
      reduction.name <-  names(seurat.object@reductions)[1]
      cli::cli_alert(glue("No reduction name provided, using first available: ", {reduction.name}))
    }

    my.gene <- feature
    umap_df <- Seurat::FetchData(
      object = seurat.object,
      #vars = c(paste0(gsub("_","",reduction.name), c("_1", "_2"), append = ""), my.gene)
      vars = c(names(seurat.object@reductions[[reduction.name]]), my.gene)
    ) %>%
      dplyr::rename(Expression = all_of(my.gene))   # Rename gene column for clarity
    colnames(umap_df)[1:2] <- c("UMAP_1","UMAP_2")

    # View the first few rows
    # head(umap_df)
    # Extract axis labels dynamically from the Seurat reduction
    axis_labels <- names(seurat.object@reductions[[reduction.name]])

    # Make sure the reduction has exactly two components for x and y
    if (length(axis_labels) >= 2) {
      x_label <- axis_labels[1]
      y_label <- axis_labels[2]
    } else {
      stop("Reduction does not contain at least two components.")
    }
    cli::cli_alert(glue("Using {reduction.name} for UMAP coordinates"))
    cli::cli_alert(glue("Using {my.gene} for expression values"))
    cli::cli_alert(glue("Using {n.bins} bins for hexbin plot"))
    cli::cli_alert(glue("Using {x_label} for x-axis label"))
    cli::cli_alert(glue("Using {y_label} for y-axis label"))

    p <- ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2, z = Expression)) +
      stat_summary_hex(
        fun = median,  # Aggregate expression within hexagons
        bins = n.bins  # Adjust bin count for resolution
      ) +
      scale_fill_viridis_c(option = "H", name = my.gene) +  # Color scale
      theme_minimal() +
      labs(
        title = paste("Hexbin Feature Plot:", my.gene),
        x = x_label,  # Set x-axis label dynamically
        y = y_label   # Set y-axis label dynamically
      ) +
      theme(
        panel.grid = element_blank(),  # Remove grid lines
        axis.line = element_line(color = "black"),  # Add x and y axis lines
        axis.ticks = element_line(color = "black"),  # Ensure ticks are visible
        axis.title = element_text(face = "bold")  # **Make sure labels are shown**
      )
    return(p)
  }



# Wrapper that preserves Seurat’s *original* color scale
.spatialFeaturePlot <- function(
  object,
  gene.i,
  limits = NULL,
  palette = NULL,
  bg.color = "#F5F3EE",
  title = NULL,
  ...) {
  #if (is.null(palette)) {palette = rev(RColorBrewer::brewer.pal(11, "Spectral"))}
  if (is.null(palette)) {
    palette=c("#E4E3DB", "#BFE8DA", "#FEE08B", "#F7B36B", "#EA8969", "#E07BA3", "#9E3D8D")

  }
  args <- list(
    object         = object,
    features       = gene.i,
    pt.size.factor = 2.0,
    image.alpha    = 0
    # rasterize         = TRUE
  )
  args <- utils::modifyList(args, list(...))

  p <- do.call(SpatialFeaturePlot, args) +
    theme_void() +
    theme(
      panel.background = element_rect(fill = bg.color, color = NA),
      plot.background  = element_rect(fill = bg.color, color = NA),
      legend.background= element_rect(fill = bg.color, color = NA),
      legend.position  = "bottom",
      legend.title     = element_text(size = 10, color = "black"),
      legend.text      = element_text(size = 8,  color = "black"),
      plot.title       = element_text(size = 14, face = "bold", hjust = 0.5, color = "black")
    )
    if(is.null(title)) {
      p <- p + ggtitle(glue::glue("{gene.i} expression"))
    }else{
      p <- p + ggtitle(title)
      }


  if (!is.null(limits)) {
    p <- p +
      scale_fill_gradientn(
        colours = palette,
        limits  = limits,
        oob     = scales::squish,
        breaks  = pretty(limits, n = 4)
      )
  }

  return(p)
}


# Discrete (categorical) spatial plot with the same look/feel as .spatialFeaturePlot
.spatialLabelPlot <- function(
  object,
  group.name,                     # metadata column or vector
  include.labels = NULL,
  color.vec = NULL,               # named vector (names = levels)
  color.other = "#E8E9DD",
  bg.color = "#F5F3EE",
  na.color = "#E8E9DD",
  legend.ncol = 3,
  title = NULL,
  ...
) {

  # --- get annotation column or vector
  annot <- if (is.character(group.name) && length(group.name) == 1 &&
      group.name %in% colnames(object@meta.data)) {
    object@meta.data[[group.name]]
  } else if (length(group.name) == ncol(object@meta.data)) {
    group.name
  } else stop("`group.name` must be a metadata column name or vector of length ncol(object).")

  if (!is.factor(annot)) annot <- factor(annot)

  # --- filter to included labels and recode others
  if (!is.null(include.labels)) {
    include.labels <- setdiff(include.labels, NA_character_)
    annot <- forcats::fct_other(annot, keep = include.labels, other_level = "other")
    annot <- forcats::fct_relevel(annot, c(include.labels, "other"))
  }

  # --- handle NA as its own level
  if (anyNA(annot)) {
    annot <- forcats::fct_na_value_to_level(annot, level = "NA")
    annot <- forcats::fct_relevel(annot, "NA", after = Inf)
    color.vec <- c(color.vec, "NA" = na.color)
  }
  my.pal <- color.vec[levels(annot)]


  # Attach and plot
  object$annot <- annot
  p <- SpatialDimPlot(
    object         = object,
    group.by       = "annot",
    cols           = my.pal,   # now guaranteed non-empty and ordered
    image.alpha    = 0,
    pt.size.factor = 2.0,
    combine        = TRUE,
    ...
  ) +
    theme_void() +
    theme(
      panel.background = element_rect(fill = bg.color, color = NA),
      plot.background  = element_rect(fill = bg.color, color = NA),
      legend.background= element_rect(fill = bg.color, color = NA),
      legend.position  = "bottom",
      legend.title     = element_text(size = 10),
      legend.text      = element_text(size = 8),
      plot.title       = element_text(size = 14, face = "bold", hjust = 0.5)
    ) +
    guides(fill = guide_legend(ncol = legend.ncol, override.aes = list(size = 4))) +
    ggtitle(if (!is.null(title)) title else group.name)
  p
}


.plot_annotation_grid <- function(
    meta_data,
  plot_annotations,
  color_palettes,
  reverse = FALSE,
  sort_levels = NULL,
  keep_input_order = FALSE            # <-- NEW ARGUMENT
) {
  library(dplyr)
  library(ggplot2)
  library(patchwork)

  # Set order of annotations
  sorting_order <- if (reverse) rev(plot_annotations) else plot_annotations

  # Set ordering of first annotation if sort_levels is provided
  if (!is.null(sort_levels)) {
    meta_data[[sorting_order[1]]] <- factor(meta_data[[sorting_order[1]]], levels = sort_levels)
  } else {
    meta_data[[sorting_order[1]]] <- factor(meta_data[[sorting_order[1]]])
  }

  # ---- KEY ADDITION: preserve input ordering if requested ----
  if (keep_input_order) {
    meta_data$cell_index <- seq_len(nrow(meta_data))     # keep original row order
  } else {
    meta_data <- meta_data %>%
      arrange(across(all_of(sorting_order))) %>%
      mutate(cell_index = row_number())
  }

  # Stack annotations
  plot_df <- purrr::map_dfr(seq_along(plot_annotations), function(i) {
    annot <- plot_annotations[i]
    data.frame(
      x = meta_data$cell_index,
      y = i,
      label = meta_data[[annot]],
      annotation = annot,
      stringsAsFactors = FALSE
    )
  })

  # Generate each bar
  plot_list <- purrr::map(plot_annotations, function(annot) {
    plot_df_sub <- dplyr::filter(plot_df, annotation == annot)
    colors <- color_palettes[[annot]]

    # Check palette mapping
    lab_levels <- sort(unique(as.character(plot_df_sub$label)))
    miss <- base::setdiff(lab_levels, names(colors))
    if (length(miss)) message("Palette missing for ", annot, ": ", paste(miss, collapse=", "))

    ggplot(plot_df_sub, aes(x = x, y = 1L, fill = label)) +
      geom_tile(width = 1, height = 1, color = NA) +
      scale_x_continuous(expand = c(0, 0)) +
      scale_y_continuous(limits = c(0.5, 1.5), expand = c(0, 0)) +
      scale_fill_manual(values = colors, na.value = "white", drop = FALSE) +
      theme_void() +
      theme(legend.position = "top", plot.margin = margin(2, 2, 2, 2)) +
      labs(title = annot)
  })

  wrap_plots(plot_list, ncol = 1)
}

.create_zoom_alluvial_plot <- function(
    data,
  annot,
  color_palette = NULL,
  project.from.mid = FALSE,         # Option to project from the middle. Produces 2 plots.
  project.from.equal.sizes = FALSE, # Option to make start groups equal in size.
  scale.all.ratio = 5,              # Default scaling ratio when no filter is provided.
  filter_column = NULL,             # Supplied filter will override the scale-all mode.
  filter_value = NULL,
  title = "Zoomed Alluvial Plot",
  project_downward = FALSE          # Option to project downward.
) {
  library(ggalluvial)
  library(ggplot2)
  library(dplyr)
  library(cowplot)
  library(cli)
  library(glue)

  # Check if the annotation column exists in the data
  missing_annotations <- base::setdiff(annot, colnames(data))
  if (length(missing_annotations) > 0) {
    stop(glue("The following annotations are missing from the data: {paste(missing_annotations, collapse = ', ')}"))
  }

  # Calculate the total frequency of the input data
  input_freq <- data %>%
    group_by(across(all_of(annot))) %>%
    summarise(Freq = n(), .groups = "drop") %>%
    pull(Freq) %>%
    sum()
  cli::cli_alert_info(glue("Total cells in input data: {input_freq}"))

  # Filter the data if filter_column and filter_value are provided
  if (!is.null(filter_column) && !is.null(filter_value)) {
    data <- data %>% dplyr::filter(!!sym(filter_column) %in% filter_value)
    scale.mode <- "scale-filtered"
    cli::cli_h2(glue("Scale mode: {scale.mode}"))
    cli::cli_h2(glue("Filtering data based on {filter_column} with value: {filter_value}"))
  } else {
    scale.mode <- "scale-all"
    cli::cli_h2(glue("Scale mode: {scale.mode}"))
    cli::cli_alert(glue("Scaling all input samples by: {scale.all.ratio}"))
  }
  # Summarize the data for the zoomed group
  grouped_data <- data %>%
    group_by(across(all_of(annot))) %>%
    summarise(Freq = n(), .groups = "drop") %>%
    mutate(x = "no scale")



  # Calculate scaling factor
  grouped_freq <- sum(grouped_data$Freq)
  cli::cli_alert_info(glue("Total cells in filtered data: {grouped_freq}"))
  zoomed_scale_factor <- if_else(scale.mode == "scale-filtered", input_freq / grouped_freq, scale.all.ratio)
  cli::cli_alert_info(glue("Zoomed group size: {grouped_freq}"))
  cli::cli_alert_info(glue("Scaling factor: {zoomed_scale_factor}"))
  # Main plot if not projecting from the middle
  if (!project.from.mid) {

    zoomed_data <- grouped_data %>%
      mutate(Freq = Freq * zoomed_scale_factor, x = "scaled")

    plot_data <- grouped_data %>% bind_rows(zoomed_data)
    n_groups <- nrow(grouped_data)

    if (project.from.equal.sizes) {
      plot_data <- plot_data %>%
        mutate(Freq = if_else(x == "no scale", grouped_freq / n_groups, Freq))
    }

    # Reverse factor levels if projecting downward
    plot_data[[annot]] <- if (project_downward) factor(plot_data[[annot]], levels = rev(levels(plot_data[[annot]]))) else factor(plot_data[[annot]])

    # Create the alluvial plot
    g <- ggplot(plot_data, aes(x = x, alluvium = !!sym(annot), y = Freq)) +
      geom_flow(aes(fill = !!sym(annot)), width = 0.1, alpha = 0.85) +
      (if (!is.null(color_palette)) scale_fill_manual(values = color_palette, guide = guide_legend(reverse = project_downward)) else NULL) +
      theme_minimal(base_size = 14) +
      theme(panel.grid = element_blank(), panel.background = element_blank(), axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank()) +
      (if (project_downward) scale_y_reverse() else NULL) +
      labs(title = title, x = NULL, y = NULL)

  } else {
    # Project from the middle (two plots)
    total_samples <- sum(grouped_data$Freq)
    half_samples <- total_samples / 2

    grouped_data <- grouped_data %>%
      arrange(!!sym(annot)) %>%
      mutate(cumsum_freq = cumsum(Freq))

    # is plot data categories uneven?
    is.uneven <- half_samples %% 2 != 0

    # Split into upper and lower halves
    upper_half <- grouped_data %>%
      filter(cumsum_freq <= half_samples | (cumsum_freq - Freq) < half_samples) %>%
      mutate(Freq = if_else(cumsum_freq > half_samples, half_samples - (cumsum_freq - Freq), Freq)) %>%
      select(-cumsum_freq)

    zoomed_upper <- upper_half %>%
      mutate(Freq = Freq * zoomed_scale_factor, x = "scaled")

    plot_upper <- upper_half %>% bind_rows(zoomed_upper)

    lower_half <- grouped_data %>%
      dplyr::filter(!(!!sym(annot) %in% upper_half[[annot]]) | cumsum_freq > half_samples)
    if(is.uneven) lower_half$Freq[1] <- lower_half$Freq[1] /2

    lower_half[[annot]] <- factor(lower_half[[annot]], levels = rev(levels(grouped_data[[annot]])))

    zoomed_lower <- lower_half %>%
      mutate(Freq = Freq * zoomed_scale_factor, x = "scaled")

    plot_lower <- lower_half %>% bind_rows(zoomed_lower)

    # Equal size adjustment for both halves
    if (project.from.equal.sizes) {
      plot_upper <- plot_upper %>%
        mutate(Freq = if_else(x == "no scale", half_samples / nrow(plot_upper), Freq))
      plot_lower <- plot_lower %>%
        mutate(Freq = if_else(x == "no scale", half_samples / nrow(plot_lower), Freq))

      # Handle overlapping groups between upper and lower halves
      split.groups <- intersect(plot_upper[[annot]], plot_lower[[annot]])

      if (length(split.groups) > 0) {
        plot_upper <- plot_upper %>%
          mutate(Freq = if_else(!!sym(annot) %in% split.groups & x == "no scale", Freq / 2, Freq))

        plot_lower <- plot_lower %>%
          mutate(Freq = if_else(!!sym(annot) %in% split.groups & x == "no scale", Freq / 2, Freq))
      }
    }

    # Upper plot
    g1 <- ggplot(plot_upper, aes(x = x, alluvium = !!sym(annot), y = Freq)) +
      geom_flow(aes(fill = !!sym(annot)), width = 0.1, alpha = 0.85) +
      (if (!is.null(color_palette)) scale_fill_manual(values = color_palette, guide = guide_legend(reverse = FALSE)) else NULL) +
      theme_void() +
      theme(plot.margin = margin(0, 0, 0, 0))

    # Lower plot (reversed)
    g2 <- ggplot(plot_lower, aes(x = x, alluvium = !!sym(annot), y = Freq)) +
      geom_flow(aes(fill = !!sym(annot)), width = 0.1, alpha = 0.85) +
      (if (!is.null(color_palette)) scale_fill_manual(values = color_palette, guide = guide_legend(reverse = TRUE)) else NULL) +
      scale_y_reverse() +
      theme_void() +
      theme(plot.margin = margin(0, 0, 0, 0))

    # Combine the two plots
    g <- cowplot::plot_grid(g1, g2, nrow = 2, align = "v", rel_heights = c(0.5, 0.5))
  }

  return(g)  # Return the final plot
}

.piechart_plot <- function(
    plot_df,
  group.var = "",
  color.var = "",
  my.pal = "",
  scaled.y = FALSE,
  scaled.y.within.color = FALSE,
  return.df.only = FALSE
) {
  library(ggplot2)
  library(dplyr)

  # Summarize data
  df <- plot_df %>%
    dplyr::select(all_of(c(group.var, color.var))) %>%
    setNames(c("group.var", "color.var"))

  df_group_sum <- df %>%
    dplyr::group_by(group.var) %>%
    summarize(group_sum = n(), .groups = "drop")

  df_color_sum <- df %>%
    dplyr::group_by(color.var) %>%
    summarize(color_var_sum = n(), .groups = "drop")

  df <- df %>%
    dplyr::group_by(group.var, color.var) %>%
    summarize(y = n(), .groups = "drop") %>%
    dplyr::left_join(df_group_sum, by = "group.var") %>%
    mutate(y_scaled = y / group_sum) %>%
    dplyr::left_join(df_color_sum, by = "color.var") %>%
    mutate(y_scaled_within_color = y / color_var_sum)

  # Normalize across all clusters
  df_scaled_all <- df %>%
    dplyr::group_by(group.var) %>%
    summarize(scaled_color_sum = sum(y_scaled_within_color), .groups = "drop")

  df <- df %>%
    dplyr::left_join(df_scaled_all, by = "group.var") %>%
    mutate(y_scaled_all = y_scaled_within_color / scaled_color_sum)

  # Apply scaling logic
  if (scaled.y && !scaled.y.within.color) {
    df <- df %>% mutate(y = y_scaled)
  } else if (!scaled.y && scaled.y.within.color) {
    df <- df %>% mutate(y = y_scaled_within_color)
  } else if (scaled.y && scaled.y.within.color) {
    df <- df %>% mutate(y = y_scaled_all)
  }

  if (return.df.only) {
    return(df)
  }

  # Define color palette
  if (my.pal[1] == "") {
    my.pal <- .palette_discrete_bp(df$color.var)
  }

  # Create the pie chart
  g <- ggplot(df, aes(x = "", y = y, fill = color.var)) +
    geom_bar(width = 1, stat = "identity", color = "black") +
    coord_polar(theta = "y") +  # Convert to pie chart
    scale_fill_manual(values = my.pal) +
    theme_void() +  # Remove background and grid
    labs(fill = color.var) +
    facet_wrap(~ group.var, scales = "free")  # One pie chart per group

  return(g)
}


.create_alluvial_plot <- function(
  data,
  annotations,
  fill_by,
  title = "Alluvial Plot",
  color_palette = NULL,
  filter_column = NULL,
  filter_value = NULL
) {
  library(ggalluvial)
  library(ggplot2)
  library(cli)
  library(dplyr)
  library(glue)

  # Check if all annotations are present in the data
  missing_annotations <- base::setdiff(annotations, colnames(data))
  if (length(missing_annotations) > 0) {
    stop(cli::cli_alert_danger(glue("The following annotations are missing from the data: {paste(missing_annotations, collapse = ', ')}")))
  }

  # Filter the data if filter_column and filter_value are provided
  if (!is.null(filter_column) && !is.null(filter_value)) {
    data <- data %>% dplyr::filter(!!sym(filter_column) == filter_value)
    cli::cli_alert_info(glue("Filtering data based on {filter_column} = {filter_value}"))
  }

  # Summarize the data on supplied annotations
  grouped_data <- data %>%
    dplyr::group_by(across(all_of(annotations))) %>%
    dplyr::summarise(Freq = n(), .groups = "drop")

  # Pretty Print Strata Information with cli
  cli::cli_h1(glue("Strata Summary for {title}"))
  total_cells <- sum(grouped_data$Freq)
  cli::cli_alert_info(glue("Total cells: {total_cells}"))

  for (annotation in annotations) {
    cli::cli_h2(glue("Annotation: {annotation}"))

    strata_summary <- grouped_data %>%
      group_by(!!sym(annotation)) %>%
      summarise(Total = sum(Freq), .groups = "drop") %>%
      arrange(desc(Total))

    # Pretty print each strata with colors
    for (i in seq_len(nrow(strata_summary))) {
      cli::cli_alert_success(glue("{strata_summary[[annotation]][i]}: {strata_summary$Total[i]} cells"))
    }
  }

  # Dynamically build the ggplot aesthetics for the axes
  aes_axes <- setNames(as.list(syms(annotations)), paste0("axis", seq_along(annotations)))

  # Create the plot
  p <- ggplot(grouped_data, aes(y = Freq, !!!aes_axes)) +
    geom_alluvium(aes(fill = !!sym(fill_by)), alpha = 1) +  # Full opacity for alluvium
    geom_stratum(width = 1/12, fill = "gray80", color = "black") +
    geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3, vjust = -0.5) +
    scale_x_discrete(limits = annotations, expand = c(0.05, 0.05)) +  # Use annotation names for x-axis
    (if (!is.null(color_palette)) scale_fill_manual(values = color_palette) else NULL) +
    theme_minimal(base_size = 14) +
    theme(
      panel.grid = element_blank(),  # Remove gridlines
      panel.background = element_blank(),  # Remove background
      axis.line = element_blank(),  # Remove axis lines
      axis.ticks = element_blank(),  # Remove axis ticks
      axis.text = element_blank(),  # Remove axis text
      axis.title.y = element_blank(),  # Remove y-axis title
      legend.position = "right"  # Keep the legend on the right
    ) +
    ggtitle(title) +
    labs(y = NULL, x = "Annotations")  # Keep x-axis title but remove y-axis

  return(p)
}



.cellchat_plots <- function(
  cellchat, # cellchat object
  cellchat.object.name = "",
  annotation,
  analysis.id,
  main.dir,
  cell_types_to_plot = c("Neoplastic", "Astrocyte reactive")
) {
  # Load required packages
  library(ggplot2)
  library(ggrepel)
  library(glue)
  library(CellChat)
  cli::cli_h1("CellChat Plot Wrapper")
  analysis.dir <- file.path(main.dir, analysis.id)
  plot.sub.dir <- file.path(main.dir, analysis.id, glue("Plots_", gsub("CellChat_object_|.rds","",cellchat.object.name)))
  cli::cli_alert_info("Creating plot directory {.file {plot.sub.dir}}")
  dir.create(plot.sub.dir, showWarnings = FALSE, recursive = TRUE)
  table(cellchat@meta$ident, cellchat@meta$sample_type)
  # Save cellchat metadata
  csv.name <- file.path(analysis.dir, glue("cellchat_idents_table_{cellchat.object.name}.csv"))
  write.csv(table(cellchat@meta$ident, cellchat@meta$sample_type), file = csv.name, row.names = T)


  # Determine group sizes and included idents
  groupSize <- as.numeric(table(cellchat@idents))
  include.idents <- names(table(cellchat@idents))

  # Define color palette
  color.use <- .color_pal[[annotation]][include.idents]

  # Plot interaction circles
  cli::cli_alert_info("Plotting interaction circles")
  pdf(file.path(plot.sub.dir, "InteractionCircle.pdf"), width = 10, height = 10)
  par(mfrow = c(1,1), xpd=TRUE)
  netVisual_circle(cellchat@net$count, sources.use = include.idents, targets.use = include.idents,
    vertex.weight = groupSize, weight.scale = TRUE, label.edge= FALSE,
    title.name = "Number of interactions")
  netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = TRUE,
    label.edge= FALSE, title.name = "Interaction weights/strength")
  dev.off()

  # Plot individual interaction weights
  cli::cli_alert_info("Plotting individual interaction weights")
  mat <- cellchat@net$weight
  pdf(file.path(plot.sub.dir, "InteractionCircles_individual_Weight.pdf"), width = 20, height = 20)
  par(mfrow = c(4,4), xpd=TRUE)
  for (i in 1:nrow(mat)) {
    mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
    mat2[i, ] <- mat[i, ]
    netVisual_circle(mat2, color.use = color.use, vertex.weight = groupSize,
      weight.scale = TRUE, arrow.size = 0.75, edge.weight.max = max(mat),
      title.name = rownames(mat)[i])
  }
  dev.off()

  # Define function to plot interactions for a specific cell type
  plot_interactions <- function(cell_type, role, cellchat, plot.sub.dir, color.use) {
    if (!(cell_type %in% rownames(cellchat@net$weight))) return(NULL)

    if (role == "receiver") {
      weights <- cellchat@net$weight[, cell_type]
      counts <- cellchat@net$count[, cell_type]
      xlab <- "Sender Cell Type"
      title_suffix <- "as_Receiver"
    } else if (role == "sender") {
      weights <- cellchat@net$weight[cell_type, ]
      counts <- cellchat@net$count[cell_type, ]
      xlab <- "Receiver Cell Type"
      title_suffix <- "as_Sender"
    } else {
      stop("Role must be 'sender' or 'receiver'")
    }
    cli::cli_alert_info("Plotting interactions for {cell_type} as {role}")
    df <- data.frame(
      CellType = names(weights),
      Weight = as.numeric(weights),
      Count = as.numeric(counts),
      Size = as.numeric(table(cellchat@idents)[names(weights)]),
      Color = color.use[names(weights)]
    )

    pdf(file.path(plot.sub.dir, glue("Interactions_{cell_type}_{title_suffix}.pdf")), width = 6, height = 5)

    g <- ggplot(df, aes(x = reorder(CellType, -Weight), y = Weight)) +
      geom_bar(stat = "identity", fill = "steelblue") +
      labs(title = glue("Interactions {title_suffix} ({cell_type}, Weight)"), x = xlab, y = "Interaction Weight") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    print(g)

    g <- ggplot(df, aes(x = reorder(CellType, -Count), y = Count)) +
      geom_bar(stat = "identity", fill = "steelblue") +
      labs(title = glue("Interactions {title_suffix} ({cell_type}, Count)"), x = xlab, y = "Interaction Count") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    print(g)

    g <- ggplot(df, aes(x = Count, y = Weight, label = CellType)) +
      geom_point(aes(color = CellType, size = Size), show.legend = FALSE) +
      scale_color_manual(values = color.use) +
      scale_size_continuous(range = c(2, 8)) +
      ggrepel::geom_text_repel(size = 3) +
      labs(
        title = glue("Interaction Count vs Weight {title_suffix} ({cell_type})"),
        x = "Number of Interactions (Count)",
        y = "Cumulative Interaction Strength (Weight)"
      ) +
      theme_minimal()
    print(g)

    dev.off()
  }

  # Plot interactions for specific cell types
  cli::cli_alert("Plotting interactions for specific cell types")
  for (cell_type in cell_types_to_plot) {
    plot_interactions(cell_type, "receiver", cellchat=cellchat, color.use = color.use, plot.sub.dir=plot.sub.dir)
    plot_interactions(cell_type, "sender", cellchat=cellchat, color.use = color.use, plot.sub.dir=plot.sub.dir)
  }

  # RankNet plots
  cli::cli_alert("Plotting RankNet")
  pdf(file.path(plot.sub.dir, "CellChat_Networks_Barplots.pdf"), width = 6, height = 14)
  p <- rankNet(cellchat, mode = "single", slot.name = "netP")
  print(p + labs(title = "Total interactions. All Receiving Cells (from all sending)"))

  for (cell_type in cell_types_to_plot) {
    if (cell_type %in% rownames(cellchat@net$weight)) {
      p_source <- rankNet(cellchat, mode = "single", slot.name = "netP", sources.use = cell_type)
      print(p_source + labs(title = glue("{cell_type} as source")))
      p_target <- rankNet(cellchat, mode = "single", slot.name = "netP", targets.use = cell_type)
      print(p_target + labs(title = glue("{cell_type} as target")))
    }
  }
  dev.off()

  # Signaling role analysis
  cli::cli_alert_info("Computing signaling role analysis")
  cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
  pdf(file.path(plot.sub.dir, "CellChat_netAnalysis_signalingRole_scatter.pdf"), width = 6, height = 6)
  p <- netAnalysis_signalingRole_scatter(cellchat, color.use = color.use,
    title = "CellChat_netAnalysis_signalingRole")
  plot(p)
  dev.off()

  cli::cli_alert_info("Computing signaling role heatmaps")
  pdf(file.path(plot.sub.dir, "CellChat_netAnalysis_signalingRole_heatmaps.pdf"), width = 8, height = 14)
  p <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing", color.use = color.use,
    title = "Outgoing interactions", width = 10, height = 28)
  plot(p)
  p <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming", color.use = color.use,
    title = "Incoming interactions", width = 10, height = 28)
  plot(p)
  dev.off()

  # Communication patterns
  cli::cli_alert_info("Computing communication patterns ... SelectK")
  pdf(file.path(plot.sub.dir, "CellChat_SelectK.pdf"), width = 12, height = 8)
  p <- selectK(cellchat, pattern = "outgoing", title = "CellChat_SelectK_outgoing")
  print(p)
  p <- selectK(cellchat, pattern = "incoming", title = "CellChat_SelectK_incoming")
  print(p)
  dev.off()

  nPatterns.out <- 6
  nPatterns.in <- 5
  pdf(file.path(plot.sub.dir, "CellChat_netAnalysis_identifyCommunicationPatterns_Outgoing_Incoming.pdf"),
    width = 12, height = 20)
  cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns.out,
    width = 6, height = 28)
  cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = nPatterns.in,
    width = 6, height = 28)
  dev.off()

  # River and dot plots
  cli::cli_alert_info("Computing river and dot plots")
  pdf(file.path(plot.sub.dir, "CellChat_netAnalysis_river.pdf"), width = 14, height = 16)
  p <- netAnalysis_river(cellchat, pattern = "outgoing")
  plot(p + labs(title = "Outgoing communication patterns"))
  p <- netAnalysis_river(cellchat, pattern = "incoming")
  plot(p + labs(title = "Incoming communication patterns"))
  p <- netAnalysis_dot(cellchat, pattern = "outgoing")
  plot(p + labs(title = "Outgoing communication patterns"))
  p <- netAnalysis_dot(cellchat, pattern = "incoming")
  plot(p + labs(title = "Incoming communication patterns"))
  dev.off()

  # Embedding plots
  cli::cli_alert_info("Computing embedding plots")
  cellchat <- computeNetSimilarity(cellchat, type = "functional")
  cellchat <- netEmbedding(cellchat, type = "functional")
  cellchat <- netClustering(cellchat, type = "functional", do.parallel = FALSE)
  cellchat <- computeNetSimilarity(cellchat, type = "structural")
  cellchat <- netEmbedding(cellchat, type = "structural")
  cellchat <- netClustering(cellchat, type = "structural", do.parallel = FALSE)

  cli::cli_alert_info("Plotting embedding and clustering")
  pdf(file.path(plot.sub.dir, "CellChat_netClustering.pdf"), width = 8, height = 8)
  p <- netVisual_embedding(cellchat, type = "functional", label.size = 3.5, title = "Functional embedding")
  plot(p)
  p <- netVisual_embedding(cellchat, type = "structural", label.size = 3.5, title = "Structural embedding")
  plot(p)
  dev.off()

  return(cellchat)
}



.liana.wrapper <- function(
  analysis.id,
  analysis.main.dir,
  seurat.object,
  liana.annot,
  n_varFeatures = 5000,
  liana.downsample.n = 1000,
  liana.min.cells = 15
  ){
  gc()

  require(Seurat)
  require(dplyr)
  require(ggplot2)
  require(Matrix)
  require(ggalluvial)
  require(cli)
  require(conflicted)
  require(glue)
  require(BiocGenerics)
  require(ComplexHeatmap)
  require(msigdbr)
  require(SCPA)
  require(liana)

  conflicted::conflicts_prefer(dplyr::filter)
  conflicted::conflicts_prefer(base::setdiff)
  set.seed(169)


  cli::cli_h1(" ===== ")
  cli::cli_h1("Starting LIANA FOO {.var {analysis.id}}")
  cli::cli_h1(" ===== ")

  # parse names
  # liana.object.name <- deparse(substitute(seurat.object))
  # cli::cli_alert("dataset: {.var {seurat.object}}")

  cli::cli_h2("Checking input")
  # results.dir must be a valid dir
  if(!dir.exists(analysis.main.dir)){
    stop("analysis.main.dir does not exist")
  }


  # stop if ref.ident is not a column in the metadata of ref.object
  if(!all(liana.annot %in% colnames(seurat.object@meta.data))){
    stop("ref.idents must be a column in the metadata of ref.object")
  }


  #  Prep folders and names
  cli::cli_h2("Prep folders and names")

  # start by creating folder structure
  out.dir <- file.path(analysis.main.dir, analysis.id)
  dir.create(out.dir, recursive = TRUE, showWarnings = FALSE)

  # Print all params to file
  params.file <- file.path(out.dir, "params.txt")
  input.params <- list(

    #liana.object.name = liana.object.name,
    seurat.object = glue("ngenes: {dim(seurat.object)[1]}, ncells: {dim(seurat.object)[2]}"),
    liana.annot = liana.annot,
    liana.downsample.n = liana.downsample.n,
    liana.min.cells = liana.min.cells
  )

  # Function to pretty print parameters
  pretty_print_params <- function(params, file_path) {
    # Open a connection to the file
    file_conn <- file(file_path, "w")

    # Write the parameters to the file
    for (param_name in names(params)) {
      param_value <- params[[param_name]]
      if (is.null(param_value)) {
        param_value <- "NULL"
      } else if (is.vector(param_value)) {
        param_value <- paste(param_value, collapse = ", ")
      }
      cat(sprintf("%s: %s\n", param_name, param_value), file = file_conn)
    }

    # Close the file connection
    close(file_conn)
  }

  # Call the function to pretty print parameters to the file
  pretty_print_params(input.params, params.file)

  # Verify the file contents
  cat(readLines(params.file), sep = "\n")

  #  Downsample data --
  cli_alert_info("Downsampling to max {.var {liana.downsample.n}} cells per subcluster")
  Idents(seurat.object) <- liana.annot
  data.ds <- subset(seurat.object, downsample = liana.downsample.n)
  table(data.ds@meta.data[[liana.annot]])

  # - Drop groups with less than liana.min.cells cells --
  data.ds@meta.data %>%
    dplyr::group_by(!!rlang::sym(liana.annot)) %>%
    # calculate n
    dplyr::summarise(n = n()) %>%
    # what groups have more than liana.min.cells n
    dplyr::filter(n > liana.min.cells) %>%
    dplyr::pull(!!rlang::sym(liana.annot)) -> groups.to.keep
  cli::cli_alert_warning("Dropping groups with less than {.var {liana.min.cells}} cells")
  cli::cli_alert_info(glue("Keeping groups: {paste0(levels(data.ds@meta.data[[liana.annot]])[groups.to.keep], collapse = ', ')}"))
  cli::cli_alert_info(glue("Dropping groups: {paste0(setdiff(levels(data.ds@meta.data[[liana.annot]]), groups.to.keep), collapse = ', ')}"))
  data.ds <- subset(data.ds, subset = !!rlang::sym(liana.annot) %in% groups.to.keep)
  data.ds@meta.data <- droplevels(data.ds@meta.data)

  # ---- Normalize ---
  cli::cli_h2("Normalize ref data")
  data.ds <- data.ds %>%
    FindVariableFeatures(nfeatures = n_varFeatures) %>%
    NormalizeData()

  # - Run LIANA --
  cli::cli_alert("Running LIANA `liana_wrap` on object with {ncol(data.ds)} cells")
  liana_res <- liana::liana_wrap(

    data.ds,
    # method = "all",
    resource = "MouseConsensus",
    min_cells = liana.min.cells,

    # added 24 apr
    method.params = list(
      cellphonedb = list(
        permutation.params = list(
          nperms = 1000,       # Increased permutations
          parallelize = TRUE,  # Use parallel computation
          workers = 4          # Number of cores
        )
      )
    ),
    # expr_prop = 0.1, # General expression filter for all methods, never implemented
    return_all = FALSE  # TRUE
    # supp_columns = c("ligand.expr", "receptor.expr", "ligand.stat", "receptor.stat",
    #   "ligand.pval", "receptor.pval", "ligand.FDR", "receptor.FDR")

  )

  liana_aggr <- liana_res %>%
    liana_aggregate()

  cli_alert_success("LIANA completed successfully")

  liana.res.file <- file.path(out.dir, "liana_res.rds")
  cli_alert_info("Saving results: {.file {liana.res.file}}")
  saveRDS(liana_res, file = liana.res.file)

  liana.aggr.file <- file.path(out.dir, "liana_aggr.rds")
  cli_alert_info("Saving results: {.file {liana.aggr.file}}")
  saveRDS(liana_aggr, file = liana.aggr.file)

  cli_alert_success("LIANA results saved successfully")
  cli_alert_success("All done!")
}



.lianaDotplotCombined <- function(
    liana_aggr_results_list,
  plot.object,
  annotation,
  source_groups,
  target_groups,
  plot.n.interactions = 25,
  analysis.name,
  alternative.plot = TRUE,
  highlight_differential = TRUE,
  expression_palette = rev(viridis::plasma(100))  # customizable
) {
  # Load required libraries
  library(dplyr)
  library(ggplot2)
  library(glue)
  library(patchwork)
  library(cli)
  library(viridis)

  # Define analysis ID
  analysis.id <- glue("{analysis.name}_UnionTop{plot.n.interactions}")

  # Initialize variables
  sample_types <- names(liana_aggr_results_list)
  mean_expr_list <- list()

  # Add interaction.complex and interaction_label early
  liana_aggr_results_list <- lapply(liana_aggr_results_list, function(x) {
    x %>% mutate(
      interaction.complex = paste(ligand.complex, receptor.complex, sep = "_"),
      interaction_label = paste(source, target, interaction.complex, sep = "__")
    )
  })

  # Prepare interaction label selection
  if (highlight_differential) {
    combined_interactions <- liana_aggr_results_list %>%
      purrr::imap_dfr(~ .x %>%
          dplyr::filter(source %in% source_groups, target %in% target_groups) %>%
          mutate(sample_type = .y))

    all_labels <- unique(combined_interactions$interaction_label)

    complete_data <- expand.grid(interaction_label = all_labels, sample_type = sample_types) %>%
      left_join(combined_interactions, by = c("interaction_label", "sample_type")) %>%
      mutate(
        aggregate_rank = ifelse(is.na(aggregate_rank), 1, aggregate_rank),
        source = ifelse(is.na(source), sub("__.*", "", interaction_label), source),
        target = ifelse(is.na(target), sub(".*__(.*?)__.*", "\\1", interaction_label), target),
        interaction.complex = ifelse(is.na(interaction.complex), sub(".*__", "", interaction_label), interaction.complex),
        ligand.complex = ifelse(is.na(ligand.complex), NA, ligand.complex),
        receptor.complex = ifelse(is.na(receptor.complex), NA, receptor.complex),
        neg_log10_rank = -log10(aggregate_rank)
      )

    significant_labels <- complete_data %>%
      dplyr::group_by(interaction_label) %>%
      dplyr::filter(any(aggregate_rank < 0.05)) %>%
      pull(interaction_label) %>%
      unique()

    diff_data <- complete_data %>%
      dplyr::filter(interaction_label %in% significant_labels) %>%
      group_by(interaction_label) %>%
      summarize(
        max_diff = max(neg_log10_rank) - min(neg_log10_rank),
        n_samples = n_distinct(sample_type)
      ) %>%
      dplyr::filter(n_samples >= 2) %>%
      arrange(desc(max_diff)) %>%
      dplyr::slice_head(n = plot.n.interactions)

    liana_selected <- complete_data %>%
      filter(interaction_label %in% diff_data$interaction_label)

    interaction_labels <- unique(liana_selected$interaction_label)

  } else {
    extract_top_interactions <- function(data, source_groups, target_groups, top_n = 50) {
      data %>%
        dplyr::filter(source %in% source_groups, target %in% target_groups) %>%
        group_by(source) %>%
        arrange(aggregate_rank) %>%
        dplyr::slice_head(n = top_n) %>%
        ungroup()
    }

    interactions_list <- lapply(names(liana_aggr_results_list), function(name) {
      extract_top_interactions(liana_aggr_results_list[[name]], source_groups, target_groups, plot.n.interactions)
    })
    names(interactions_list) <- names(liana_aggr_results_list)

    interactions_list <- lapply(interactions_list, function(x) {
      x %>% mutate(interaction_label = paste(source, target, interaction.complex, sep = "__"))
    })
    interaction_labels <- unique(unlist(lapply(interactions_list, function(x) x$interaction_label)))
  }

  # Combine selected interactions
  filtered_list <- lapply(names(liana_aggr_results_list), function(name) {
    liana_aggr_results_list[[name]] %>%
      dplyr::filter(source %in% source_groups, target %in% target_groups) %>%
      dplyr::filter(interaction_label %in% interaction_labels) %>%
      mutate(sample_type = name) %>%
      group_by(interaction_label) %>%
      dplyr::slice(1) %>%
      ungroup()
  })

  liana_selected <- bind_rows(filtered_list) %>%
    mutate(interaction_label = factor(interaction_label, levels = unique(interaction_label))) %>%
    arrange(aggregate_rank)

  # Compute cell type sizes
  source_sizes <- plot.object@meta.data %>%
    dplyr::filter(!!sym(annotation) %in% source_groups) %>%
    group_by(!!sym(annotation), sample_type) %>%
    summarize(n_cells_source = n(), .groups = "drop") %>%
    dplyr::rename(source = !!sym(annotation))

  target_sizes <- plot.object@meta.data %>%
    dplyr::filter(!!sym(annotation) %in% target_groups) %>%
    group_by(!!sym(annotation), sample_type) %>%
    summarize(n_cells_target = n(), .groups = "drop") %>%
    dplyr::rename(target = !!sym(annotation))

  # Merge interaction data with cell sizes
  liana_plot <- liana_selected %>%
    left_join(source_sizes, by = c("source", "sample_type")) %>%
    left_join(target_sizes, by = c("target", "sample_type")) %>%
    group_by(interaction_label) %>%
    arrange(desc(n_cells_target), .by_group = TRUE) %>%
    ungroup() %>%
    mutate(
      neg_log10_rank = -log10(aggregate_rank),
      source_order = match(source, source_groups),
      target_order = match(target, target_groups)
    ) %>%
    group_by(interaction_label) %>%
    mutate(max_neg_log10_rank = max(neg_log10_rank, na.rm = TRUE)) %>%
    ungroup() %>%
    arrange(source_order, target_order, desc(max_neg_log10_rank), desc(n_cells_target)) %>%
    mutate(interaction_label = factor(interaction_label, levels = unique(interaction_label)))

  # Global size range
  size_range <- range(c(liana_plot$n_cells_source, liana_plot$n_cells_target), na.rm = TRUE)

  # Main dotplot
  p_main <- ggplot(liana_plot, aes(x = neg_log10_rank, y = interaction_label, color = sample_type, size = n_cells_target)) +
    geom_point(alpha = 0.85) +
    scale_color_manual(values = .color_pal[["sample_type"]]) +
    scale_size(range = c(2, 10), limits = size_range, name = "Target group size (cells)") +
    scale_y_discrete(limits = rev) +
    theme_minimal(base_size = 13) +
    labs(
      x = "-log10(Aggregate Rank)",
      y = "Source_Target_Interaction",
      size = "Group size (cells)",
      color = "Sample Type"
    ) +
    theme(
      panel.grid = element_line(color = "grey90"),
      axis.line = element_line(color = "black"),
      axis.ticks = element_line(color = "black"),
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      legend.position = "right",
      legend.key = element_rect(fill = "white", color = "black"),
      plot.margin = margin(t = 10, r = 30, b = 10, l = 30)
    ) +
    ggtitle("Top ranked interactions: ligand → receptor") +
    geom_vline(xintercept = -log10(0.05), linetype = "dashed", color = "black") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black")

  # Expression calculation
  ligand_genes <- unique(unlist(strsplit(unique(liana_plot$ligand.complex), "_")))
  receptor_genes <- unique(unlist(strsplit(unique(liana_plot$receptor.complex), "_")))
  all_genes <- unique(c(ligand_genes, receptor_genes))

  for (st in sample_types) {
    cells_st <- colnames(plot.object)[plot.object$sample_type == st]
    seurat_st <- subset(plot.object, cells = cells_st)

    for (src in source_groups) {
      cells_src <- colnames(seurat_st)[seurat_st[[annotation]] == src]
      if (length(cells_src) > 0) {
        mean_expr_list[[paste(st, src, "source", sep = "_")]] <- colMeans(FetchData(seurat_st, vars = all_genes, cells = cells_src, layer = "data"))
      }
    }

    for (tgt in target_groups) {
      cells_tgt <- colnames(seurat_st)[seurat_st[[annotation]] == tgt]
      if (length(cells_tgt) > 0) {
        mean_expr_list[[paste(st, tgt, "target", sep = "_")]] <- colMeans(FetchData(seurat_st, vars = all_genes, cells = cells_tgt, layer = "data"))
      }
    }
  }

  # Expression range
  all_expr_values <- unlist(lapply(mean_expr_list, function(x) unname(unlist(x))))
  expr_range <- range(all_expr_values, na.rm = TRUE)

  # Return
  return(list(
    combined_plot = p_main,
    liana_plot = liana_plot,
    mean_expr_list = mean_expr_list,
    expr_range = expr_range
  ))
}


.liana_per_sample_plot <- function(
    liana_plot,
  mean_expr_list,
  expr_range,
  sample_type_input,  # <--- argument must match name used in call
  source_groups,
  target_groups,
  expression_palette = rev(viridis::plasma(100))
) {
  # Helper to get expression values
  get_expr_value <- function(gene_complex, celltype, sample_type, type) {
    key <- paste(sample_type, celltype, type, sep = "_")
    genes <- strsplit(gene_complex, "_")[[1]]
    if (!key %in% names(mean_expr_list)) return(list(expr = NA_real_, gene = NA_character_))
    expr_vector <- mean_expr_list[[key]]
    expr_values <- sapply(genes, function(g) expr_vector[[g]])
    expr_values <- expr_values[!is.na(expr_values)]
    if (length(expr_values) == 0) return(list(expr = NA_real_, gene = NA_character_))
    max_expr <- max(expr_values)
    max_gene <- names(expr_values)[which.max(expr_values)]
    return(list(expr = max_expr, gene = max_gene))
  }

  # Threshold
  significance_threshold <- -log10(0.05)

  # Prepare df for one sample
  df <- liana_plot %>%
    filter(sample_type == sample_type_input) %>%
    mutate(
      x_source = -5,
      x_target = 5,
      y_pos = as.numeric(interaction_label),
      is_significant = neg_log10_rank > significance_threshold
    ) %>%
    rowwise() %>%
    mutate(
      source_expr_info = list(get_expr_value(ligand.complex, source, sample_type_input, "source")),
      target_expr_info = list(get_expr_value(receptor.complex, target, sample_type_input, "target")),
      source_expr = source_expr_info$expr,
      target_expr = target_expr_info$expr,
      receptor_gene_label = target_expr_info$gene,
      fill_source = ifelse(is_significant, source_expr, NA_real_),
      fill_target = ifelse(is_significant, target_expr, NA_real_)
    ) %>%
    ungroup()

  label_offset <- 1.5

  # Plot
  p <- ggplot(df) +
    geom_point(aes(x = x_source, y = y_pos, size = n_cells_source, fill = fill_source), shape = 21, na.rm = TRUE) +
    geom_point(aes(x = x_target, y = y_pos, size = n_cells_target, fill = fill_target), shape = 21, na.rm = TRUE) +
    geom_segment(
      data = df %>% filter(is_significant),
      aes(x = x_source + 3, xend = x_target - 3, y = y_pos, yend = y_pos),
      arrow = arrow(length = unit(0.2, "cm")), color = "black"
    ) +
    geom_text(aes(x = x_source - label_offset, y = y_pos, label = ligand.complex), hjust = 1.2, size = 3) +
    geom_text(aes(x = x_target + label_offset, y = y_pos, label = receptor_gene_label), hjust = 0, size = 3) +
    scale_size(range = c(2, 10), limits = range(c(df$n_cells_source, df$n_cells_target), na.rm = TRUE), name = "Group size (cells)") +
    scale_fill_gradientn(
      colours = expression_palette,
      name = "Expression",
      limits = expr_range,
      na.value = "transparent"
    ) +
    scale_x_continuous(limits = c(-10, 10)) +
    scale_y_reverse() +
    theme_void() +
    ggtitle(sample_type_input) +
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.position = "right",
      legend.key = element_rect(fill = "white", color = "black"),
      plot.margin = margin(t = 10, r = 30, b = 10, l = 40)
    )

  return(p)
}

.liana_generate_all_per_sample_plots <- function(
    sample_types,
  liana_plot,
  mean_expr_list,
  expr_range,
  source_groups,
  target_groups,
  expression_palette = rev(viridis::plasma(100))
) {
  per_sample_plots <- lapply(sample_types, function(st) {
    .liana_per_sample_plot(
      liana_plot = liana_plot,
      mean_expr_list = mean_expr_list,
      expr_range = expr_range,
      sample_type_input = st,  # pass it under new name
      source_groups = source_groups,
      target_groups = target_groups,
      expression_palette = expression_palette
    )
  })

  # Name the list for convenience
  names(per_sample_plots) <- sample_types

  # Replace NULL with empty plot (for consistency, though per_sample_plot shouldn't return NULL)
  per_sample_plots <- lapply(per_sample_plots, function(p) {
    if (is.null(p)) ggplot() + theme_void() else p
  })

  return(per_sample_plots)
}

# ---- .color_pal ----
.color_pal <- list(
  # 25 clusters SNN_p15_res05_qF03dbl10
  "ClusterNames_v1_qF03_snn05" = c(
    "01-Tumor" = "#F4E409",
    "02-Tumor" = "#EB442C",
    "05-Tumor-G1SG2" = "#A63C06",
    "08-Tumor-TR" = "#FFB915",
    "14-Tumor" = "#CA7B09",
    "06-Tumor-G1SG2" = "#8C0000",
    "04-Tumor-OPC-TH" = "#FFC8C4",
    "18-Tumor-OPC-TH" = "#FF97C4",
    "17-Astro-T" = "#F05E16",
    "16-Astrocyte-H" = "#F00000",
    "09-Astrocyte-H" = "#E12C2C",
    "10-Oligodendrocyte" = "#3F1D0B",
    "21-Mural" = "#0B6E4F",
    "11-Endothelial" = "#08A045",
    "15-Neuron" = "#6BBF59",
    "07-Macrophage-TR" = "#03045e",
    "03-Microglia-H" = "#0077b6",
    "13-Microglia-T" = "#00b4d8",
    "12-Microglia" = "#003F66",
    "22-Monocytes" = "#00A5CF",
    "25-Granulocyte" = "#caf0f8",
    "23-Macrophage" = "#002A45",
    "20-NK-Tcell" = "#AF9AB2",
    "19-Choroid-H" = "#21D375",
    "24-Ependymal" = "#D1C3D3"
  ),

  "ClusterNames_v2_qF03ds_snn05" = c(
    "01-Tumor" = "#F4E409",
    "02-Tumor" = "#EB442C",
    "06-Tumor-G1SG2" = "#A63C06",
    "05-Tumor-TR" = "#FFB915",
    "14-Tumor-OPC-TH" = "#FFC8C4",
    "04-Tumor-OPC-TH" = "#FF97C4",
    "12-Astro-T" = "#F05E16",
    "13-Astrocyte-H" = "#F00000",
    "10-Astrocyte-H" = "#E12C2C",
    "09-Oligodendrocyte" = "#3F1D0B",
    "20-Fibroblasts" = "#CA7B09",
    "15-Mural" = "#8C0000",
    "22-Pericytes" = "#8C0000",
    "11-Endothelial" = "#08A045",
    "17-Neuron" = "#6BBF59",
    "19-Neuron" = "#0B6E4F",
    "08-Macrophage-TR" = "#03045e",
    "03-Microglia-H" = "#0077b6",
    "07-Microglia" = "#003F66",
    "18-Monocytes" = "#00A5CF",
    "18-BloodCells" = "#00A5CF",
    "24-Granulocyte" = "#caf0f8",
    "21-Macrophage" = "#002A45",
    "15-NK-Tcell" = "#AF9AB2",
    "16-Choroid-H" = "#21D375",
    "23-Ependymal" = "#D1C3D3"
  ),


  m1v2_Tx_snn05_Pds10k_C1G1 = c(
    "07-Astro-H" = "#8CAB00",
    "14-Astro-H" = "#00BFC4",

    "09-Olig" = "#24B700",

    "18-Neuron-H" = "#42A0FF",

    "17-Choroid" = "#00ACFC",
    "25-Ependymal-H" = "#FF65AC",

    "02-Microglia-H" = "#EE8043",
    "04-Microglia-TA" = "#D19300",
    "08-Macrophage-TA" = "#68B100",
    "23-Macrophage-H" = "#F962DD",

    "19-NKT" = "#8B93FF",
    "20-Leukocyte-TA" = "#B684FF",
    "26-Leukocyte-TA" = "#FF6D8E",

    "22-Fibroblast" = "#EB69F0",
    "24-Pericyte-H" = "#FF61C6",
    "13-Endothelial" = "#00C1AB",

    "15-OPC-TH" = "#00BBDA",
    "16-Astro-TA" = "#00B5ED",
    "21-Neuron-TH" = "#D575FE",
    "06-Tumor-OPC" = "#A8A400",

    "01-Tumor" = "#F8766D",
    "11-Tumor" = "#00BE70",
    "12-Tumor" = "#00C090",
    "05-Tumor-G1SG2" = "#BE9C00",
    "03-Tumor-G2M" = "#E18A00",
    "10-Tumor-TR" = "#00BB49"
  ),



  ## ---- Healthy Tx ----
  "m1v2_Tx_snn05_Pds10k_Healthy_C1G1" = c(
    "Healthy" = "#F8D3BAFF",
    "Primary" = "#EC4B3EFF",
    "Recurrent" = "#611F53FF",
    "Tumor" = "#EC4B3EFF",

    "02-Astro-H" = "#F4E409",
    "06-Astro-H" = "#FFB915",
    "15-Astro-H" = "#CA7B09",
    "Astro-H" = "#CA7B09",
    "04-Olig" = "#7700FF",
    "04-Olig-H" = "#7700FF",
    "11-Olig-H" = "#03045e",
    "Olig-H" = "#03045e",
    "09-OPC-H" = "#FF00C3",
    "14-OPC-H" = "#FF97C4",
    "10-Neuron-H" = "#FFC8C4",
    "10-Neuron-gaba-H" = "#FFC8C4",
    "16-Neuron-H" = "#caf0f8",
    "16-Neuron-bdnf-sst-H" = "#caf0f8",
    "19-Neuron-H" = "#FF9898FF",
    "08-Choroid-H" = "#855C75FF",
    "17-Ependymal-H"  = "#8C785DFF",
    "01-Microglia-H" = "#B4BF3AFF",
    "03-Microglia-H" = "#5E9432FF",
    "07-Microglia-H"= "#88AB38FF",
    "14-Microglia-H" = "#21D375",
    "Microglia-H" = "#21D375",
    "12-Macrophage-H" = "#225F2FFF",
    "24-Macrophage-H" = "#237F3FFF",
    "24-Macro-dend-H" = "#237F3FFF",
    "23-Tcell-H" = "#05AFF2FF",
    "23-NKT-H" = "#05AFF2FF",
    "25-BCell-H" =  "#00A5CF",
    "05-Endo-H" = "#800000FF",
    "13-Pericyte-H"  = "#E12C2C",
    "18-Mural-H" = "#F05E16",
    "21-Fibro-H" = "#D6BB3BFF",
    "22-Fibro-H"  = "#A63C06",
    "21-Stromal-H" = "#D6BB3BFF",
    "22-Stromal-H"  = "#A63C06",
    "other" = "gray"
  ),

  ## ---- Sample id ----
  "sample_id" = c(
    # BRAIN_ group: shades of green
    "BRAIN_01" = "#B4BF3A",  # Light green
    "BRAIN_02" = "#A2C523",  # Slightly darker light green
    "BRAIN_03" = "#90C12F",  # Medium-light green
    "BRAIN_04" = "#8FBC8F", # Dark sea green
    "BRAIN_04a" = "#8FBC8F", # Dark sea green
    "BRAIN_04b" = "#3CB371", # Medium sea green
    "BRAIN_05" = "#2E8B57",  # Sea green
    "BRAIN_06" = "#006400",  # Dark green

    # TP_ group: variants of blue, pink, and purple
    "TP_01" = "#05DBF2",     # Light blue
    "TP_01a" = "#00E5FF",    # Light blue (same as TP_01)
    "TP_01b" = "#00C4E4",    # Light blue (same as TP_01)
    "TP_02" = "#87CEEB",     # Sky blue
    "TP_03" = "#76B5E5",     # Slightly darker sky blue
    "TP_04" = "#6495ED",     # Cornflower blue
    "TP_05" = "#4169E1",     # Royal blue
    "TP_06" = "#6A5ACD",     # Slate blue
    "TP_07" = "#7B68EE",     # Medium slate blue
    "TP_08" = "#8A2BE2",     # Blue violet

    # TR_ group: variants of yellow, red, and brown
    "TR_01" = "#FFD700",     # Gold
    "TR_02" = "#FF8C00",     # Dark orange
    "TR_03" = "#FF5722",     # Orange red
    "TR_04" = "#FF0000",     # Bright red
    "TR_05" = "#DC143C",     # Crimson (deep red)
    "TR_06" = "#8B4513",     # Saddle brown
    "TR_07" = "#A52A2A"     # Brown

    # Default for other groups
    # "other" = "gray"         # Gray for other
  ),



  ## ---- orig_ident ----
  orig_ident = c(
    # BRAIN_ group: shades of green
    "BRAIN_01" = "#B4BF3A",  # Light green
    "BRAIN_02" = "#A2C523",  # Slightly darker light green
    "BRAIN_03" = "#90C12F",  # Medium-light green
    "BRAIN_04" = "#8FBC8F", # Dark sea green
    "BRAIN_04a" = "#8FBC8F", # Dark sea green
    "BRAIN_04b" = "#3CB371", # Medium sea green
    "BRAIN_05" = "#2E8B57",  # Sea green
    "BRAIN_06" = "#006400",  # Dark green

    # TP_ group: variants of blue, pink, and purple
    "TP_01" = "#05DBF2",     # Light blue
    "TP_01a" = "#00E5FF",    # Light blue (same as TP_01)
    "TP_01b" = "#00C4E4",    # Light blue (same as TP_01)
    "TP_02" = "#87CEEB",     # Sky blue
    "TP_03" = "#76B5E5",     # Slightly darker sky blue
    "TP_04" = "#6495ED",     # Cornflower blue
    "TP_05" = "#4169E1",     # Royal blue
    "TP_06" = "#6A5ACD",     # Slate blue
    "TP_07" = "#7B68EE",     # Medium slate blue
    "TP_08" = "#8A2BE2",     # Blue violet

    # TR_ group: variants of yellow, red, and brown
    "TR_01" = "#FFD700",     # Gold
    "TR_02" = "#FF8C00",     # Dark orange
    "TR_03" = "#FF5722",     # Orange red
    "TR_04" = "#FF0000",     # Bright red
    "TR_05" = "#DC143C",     # Crimson (deep red)
    "TR_06" = "#8B4513",     # Saddle brown
    "TR_07" = "#A52A2A"     # Brown

    # Default for other groups
    # "other" = "gray"         # Gray for other
  ),

  ## ---- SNNN SNN_clusters_all -----
  SNN_clusters_all = c(
    # --- non-neoplastic (16 shades of sage; light → mid) ---
    "02" = "#EAF2E9",
    "04" = "#E1ECDD",
    "07" = "#D6E6D2",
    "08" = "#CCE0C6",
    "09" = "#C2DAC0",
    "13" = "#B9D4B6",
    "14" = "#B0CEAD",
    "17" = "#A6C8A4",
    "18" = "#9EC39C",
    "19" = "#97BE96",
    "20" = "#90B990",
    "22" = "#88B489",
    "23" = "#81AF84",
    "24" = "#7AAA7E",
    "25" = "#73A478",
    "26" = "#6C9F72",

    # --- mixed (4 shades of ochre; light → deeper) ---
    "06" = "#F2DFBF",
    "15" = "#E9C898",
    "16" = "#D8A372",
    "21" = "#C3895B",

    # --- neoplastic (6 shades of plum; light → deeper) ---
    "01" = "#E8DCE7",
    "03" = "#D6BDD1",
    "05" = "#C39EBB",
    "10" = "#B27FA5",
    "11" = "#A26790",
    "12" = "#955C80"
  ),
#
#   SNN_clusters_all = c(
#     # --- Glial / Neural lineage (soft oranges → warm pinks → muted purples) ---
#     "02" = "#F5B971",  # Astro_TE
#     "07" = "#F28C5C",  # Astro_NT
#     "15" = "#E87A7A",  # Astro_TA
#     "24" = "#F1C6A5",  # Ependymal
#     "11" = "#C97EBA",  # Interneuron
#     "18" = "#C4A3E0",  # Mixed neurons
#
#     # --- OPC / Oligo (cooler purples and blue-greys) ---
#     "03" = "#9FA8DA",  # Oligodendrocyte
#     "14" = "#7986CB",  # OPC
#     "12" = "#6E85B7",  # OPC_TA
#     "20" = "#5C6BC0",  # COP_TA
#
#     # --- Immune – lymphoid / myeloid (greens → brownish neutrals) ---
#     "13" = "#A5D6A7",  # NKT
#     "25" = "#81C784",  # B/pDC
#     "01" = "#B0BEC5",  # Microglia
#     "09" = "#90A4AE",  # Microglia
#     "08" = "#78909C",  # Microglia_TAM
#     "06" = "#8BC34A",  # TAM (macrophage phenotype)
#     "16" = "#689F38",  # BMD-TAM
#     "04" = "#558B2F",  # Macrophage_BMD_TAM
#     "17" = "#AED581",  # Dendritic
#     "19" = "#9CCC65",  # BAM
#     "26" = "#BDBDBD",  # Neutrophil (muted neutral)
#
#     # --- Vascular / Stromal (warm beiges → olive browns) ---
#     "10" = "#FFF176",  # Choroid
#     "05" = "#D7CCC8",  # Endothelial
#     "23" = "#BCAAA4",  # Endothelial (tumor-assoc)
#     "21" = "#A1887F",  # Mural
#     "22" = "#8D6E63",  # Fibroblast
#     "27" = "#AFAFAF",  # Fibroblast_meningeal (neutral grey)
#
#     # --- Neoplastic (rose to muted red-pink) ---
#     "18" = "#E0B1B8",  # (if reused)
#     "12" = "#E57373",  # Neoplastic / OPC-like
#     "15" = "#EC9BAE"   # Neoplastic (blended hue)
#   ),

  ## ---- SNN clusters broad ----
  # SNN_Clusters_broad = c(
  #   "Non-neoplastic Clusters" = "#FAEBDDFF",
  #   "Mixed Clusters" = "#F4875EFF",
  #   "Neoplastic Clusters" = "#E91E63"
  # ),

  SNN_Clusters_broad = c(
    "Non-neoplastic Clusters" = "#8FB98B",  # soft beige/sand
    "Mixed Clusters"          = "#D8A372",  # muted sage green
    "Neoplastic Clusters"     = "#955C80"   # dusty plum/purple-brown
  ),

  ## ---- SNN_clusters_nonNeoplastic----
  SNN_clusters_nonNeoplastic = c(
    # Immune – myeloid & microglia
    "01" = "#ED8A68",  # Microglia
    "04" = "#CC1212",  # Macrophage (BMD-TAM)
    "06" = "#CC1212",  # Macrophage (BMD-TAM, macrophage phenotype)
    "08" = "#ED8A68",  # Microglia-like TAM (microglia phenotype)
    "09" = "#ED8A68",  # Microglia
    "16" = "#CC1212",  # Macrophage (BMD-TAM)
    "17" = "#FF69B4",  # Dendritic (cDC)
    "19" = "#CC1212",  # BAM → macrophage family
    "26" = "#800080",  # Neutrophil

    # Glial lineage
    "02" = "#0B3800",  # Astrocyte TE  (Level_3AC: Astrocyte TE)
    "07" = "#82AC7C",  # Astrocyte NT  (Level_3AC: Astrocyte NT)
    "15" = "#B8D53D",  # Astrocyte TA  → mapped to Astrocyte R (reactive)

    "24" = "#FFEA70",  # Ependymal

    # Neuronal
    "11" = "#3CB371",  # Inhibitory interneuron → Neural
    "18" = "#3CB371",  # Neuron mixed → Neural

    # OPC / Oligo
    "03" = "#0065A2",  # Oligodendrocyte → OPC/COP/OLG
    "12" = "#0065A2",  # OPC_TA → OPC/COP/OLG
    "14" = "#0065A2",  # OPC → OPC/COP/OLG
    "20" = "#0065A2",  # COP_TA → OPC/COP/OLG

    # Immune – lymphoid & pDCs
    "13" = "#87CEEB",  # NK/T/B → NKTB
    "25" = "#87CEEB",  # B + pDC mix → NKTB family color (split across NKTB/Dendritic)

    # Vascular & stromal
    "05" = "#60330D",  # Endothelial (healthy-dominated)
    "10" = "#FFD700",  # Choroid
    "21" = "#DF8008",  # Mural (Pericyte/SMC)
    "22" = "#E4D938",  # Fibroblast
    "23" = "#60330D",  # Endothelial (tumor-assoc., mixed; kept as Endothelial)
    "27" = "#E4D938"   # Fibroblast (meningeal)
  ),


  ## ---- project id ----
  "project_id" = c(
    "2023_069_1" = "#FAEBDDFF",
    "2023_069_2" = "#F6BB97FF",
    "2023_141_1" = "#F4875EFF",
    "2023_141_2" = "#EC4B3EFF",
    "2023_141_3" = "#961C5BFF"
  ),
  "sample_type" = c(
    "Healthy" = "#9FC49A",
    "Primary" = "#6BAED6",
    "Recurrent" = "#D9A55A"
  ),
  "disease_state" = c(
    "Healthy" = "#F6A178FF",
    "Tumor" = "#611F53FF"
  ),
  "matrisome" = c(
    "ECM-affiliated Protein" = "#EEBA0B",
    "ECM Glycoproteins" = "#CA7B09",
    "ECM Regulators" = "#A63C06",
    "Non-matrisome" = "#F6A178FF",
    "Proteoglycans" = "#F00000",
    "ECM Regulators" = "#0B6E4F",
    "Secreted Factors" = "#0077b6",
    "Collagens" = "#00b4d8"
  ),
  "flunarizine_Condition" = c(
    "nonIR_DMSO" = "#F4E409",
    "IR_DMSO" = "#FFB915",
    "nonIR_Flunarizine" = "#A63C06",
    "IR_Flunarizine" = "#EB442C",
    "nonIR_Flunarizine_r1-r3" = "#FFC8C4",
    "nonIR_Flunarizine_r4-r5" = "#FF97C4",
    "IR_Flunarizine_r1-r3" = "#EB442C",
    "IR_Flunarizine_r4-r5" = "#A63C06"

  ),
  ## ---- RCAS ----
  "hPDGFB_bin" =  c(
    "count-00" = "#FAEBDDFF",
    "count-01" = "#F4875EFF",
    "count-02-05" = "#CB1B4FFF",
    "count-06+" = "#611F53FF"
  ),
  "RFP_1_bin" =  c(
    "count-00" = "#FAEBDDFF",
    "count-01" = "#F4875EFF",
    "count-02-05" = "#CB1B4FFF",
    "count-06+" = "#611F53FF"
  ),
  "RFP_2_bin" =  c(
    "count-00" = "#FAEBDDFF",
    "count-01" = "#F4875EFF",
    "count-02-05" = "#CB1B4FFF",
    "count-06+" = "#611F53FF"
  ),

  "rcas_strict" =  c(
    "rcas_neg" = "#FAEBDDFF",
    "rcas_strict" = "#99300CFF"
  ),
  "rcas_both" =  c(
    "rcas_neg" = "#FAF7D2FF",
    "hpdgfb_pos" = "#BE852CFF",
    "rfp_pos" = "#E1CA89FF",
    "rcas_pos" = "#99300CFF"
  ),
  "rcas_both_stricter" =  c(
    "rcas_neg" = "#FAEBDDFF",
    "hpdgfb_pos" = "#BE852CFF",
    "rfp_pos" = "#E1CA89FF",
    "rcas_pos" = "#99300CFF"
  ),
  "rcas_both_3bin" =  c(
    "rcas_neg" = "#FAF7D2FF",
    "rcas_ambiguous" = "#E1CA89FF",
    #"rcas_ambiguous" = "#BE852CFF",
    "rcas_pos" = "#99300CFF"
  ),
  "rcas_core_signature_3bin" =  c(
    "rcas_neg" = "#FAEBDDFF",
    "rcas_ambiguous" = "#99300CFF",
    "rcas_pos" = "#F00000"
  ),
  "rcas_pdgfb_nha_signature_3bin" =  c(
    "pdgfb_nha_neg" = "#FAEBDDFF",
    "pdgfb_nha_ambiguous" = "#99300CFF",
    "pdgfb_nha_pos" = "#F00000"
  ),

  "rcas_rfp_12_signature_3bin" =  c(
    "rfp_12_neg" = "#FAEBDDFF",
    "rfp_12_ambiguous" = "#99300CFF",
    "rfp_12_pos" = "#F00000"
  ),

  "rcas_module_bin" =  c(
    "rcas_module_neg" = "#FAEBDDFF",
    "rcas_module_amb" = "#D6BB3B",
    "rcas_module_pos" = "#F00000"
  ),
  "rcas_counts_bin" =  c(
    "rcas_neg"= "#FAEBDDFF",
    "rcas_ambiguous" = "#EAD98E",
    "rcas_pos" ="#F00000"
  ),
  "rcas_call" =  c(
    "rcas_neg"= "#E4D6C8",
    "rcas_ambiguous" = "#E1CA66",
    "rcas_low" = "#E1CA66", # same as rcas amb
    "rcas_pos" ="#F00000",
    "rcas_high" ="#F00000" # same as rcas pos
  ),
  "SNN_Clusters_broad" =  c(
    "Non-neoplastic"= "#FAEBDDFF",
    "Mixed" = "#80100CFF",
    "Neoplastic" ="#F00000"
  ),

  ## ---- subclusters ----
  "subcluster" =  c(
    "01" = "#F6A178FF",
    "02" = "#F26A48FF",
    "03" = "#CB1B4FFF",
    "04" = "#611F53FF",
    "other" = "#FAEBDDFF"

  ),

  "subcluster_sd_sub" =  c(
    "01_sd.subset" = "#F6BB97FF",
    "02_sd.subset" = "#F26A48FF",
    "03_sd.subset" = "#CB1B4FFF",
    "04_sd.subset" = "#611F53FF",
    "05_sd.subset" = "#EEBA0B",
    "06_sd.subset" = "#CA7B09",
    "07_sd.subset" = "#00b4d8",
    "08_sd.subset" = "#0077b6",
    "09_sd.subset" = "#03045e",
    "10_sd.subset" = "#0B6E4F",
    "11_sd.subset" = "#08A045",
    "12_sd.subset" = "#6BBF59",
    "13_sd.subset" = "#F4E409",
    "14_sd.subset" = "#FFB915",
    "15_sd.subset" = "#caf0f8",
    "other" = "#FAEBDDFF"
  ),
  "snn025_sd.subset" =  c(
    "01_sd.subset" = "#F6BB97FF",
    "02_sd.subset" = "#F26A48FF",
    "03_sd.subset" = "#CB1B4FFF",
    "04_sd.subset" = "#611F53FF",
    "05_sd.subset" = "#EEBA0B",
    "06_sd.subset" = "#CA7B09",
    "07_sd.subset" = "#00b4d8",
    "08_sd.subset" = "#0077b6",
    "09_sd.subset" = "#03045e",
    "10_sd.subset" = "#0B6E4F",
    "11_sd.subset" = "#08A045",
    "12_sd.subset" = "#6BBF59",
    "13_sd.subset" = "#F4E409",
    "14_sd.subset" = "#FFB915",
    "15_sd.subset" = "#caf0f8",
    "other" = "#FAEBDDFF"
  ),

  ### ---- Level_1 ----
  Level_1 = c(
    "Non-Neoplastic" = "#BBD8B8",  # same as before
    "Neoplastic"     = "#F4B7C7",
    "Ambiguous"      = "gray",
    "other"          = "gray",
    "NA"             = "gray"
  ),
  ### ---- Level_2 ----
  Level_2 = c(
    "Neural"         = "#82AC7C",   # from Glial-Neuronal
    "Immune"         = "#CC1212",   # use Myeloid red as Immune umbrella
    "Vascular/Stromal" = "#ED8A68", # from Vascular
    "Neoplastic"     = "#F4B7C7",
    "Ambiguous"      = "gray",
    "other"          = "gray"
  ),
  ### ---- Level_3 ----
  Level_3 = c(
    # Neural
    "Astrocyte"      = "#82AC7C",
    "OPC/COP/OLG"      = "#0065A2",
    "Neural"         = "#3CB371",
    "Ependymal"      = "#FFEA70",

    # Immune
    "Microglia"      = "#ED8A68",
    "Macrophage"     = "#CC1212",
    "Dendritic"      = "#FF69B4",
    "Neutrophil"    = "#800080",
    "NKTB"           = "#87CEEB",

    # Vascular/Stromal
    "Choroid"        = "#FFD700",
    "Endothelial"    = "#60330D",
    "Mural"          = "#DF8008",
    "Fibroblast"     = "#E4D938",

    # Neoplastic & Ambiguous
    "Neoplastic"     = "#F4B7C7",
    "Ambiguous"      = "gray",
    "other"          = "gray"
  ),

  ### ---- Level_3AC ----
  Level_3AC = c(
    # Neural
    "Astrocyte TE"       = "#0B3800",
    "Astrocyte NT"       = "#82AC7C",
    "Astrocyte TE/NT"       = "#82AC7C",
    "Astrocyte R" = "#B8D53D",
    "OPC/COP/OLG"          = "#0065A2",
    "Neural"             = "#3CB371",
    "Ependymal"          = "#FFEA70",

    # Immune
    "Microglia"          = "#ED8A68",
    "Macrophage"         = "#CC1212",
    "Dendritic"          = "#FF69B4",
    "Neutrophil"        = "#800080",
    "NKTB"               = "#87CEEB",

    # Vascular/Stromal
    "Choroid"            = "#FFD700",
    "Endothelial"        = "#60330D",
    "Mural"              = "#DF8008",
    "Fibroblast"         = "#E4D938",

    # Neoplastic & Ambiguous
    "Neoplastic"         = "#F4B7C7",
    "Ambiguous"          = "gray",
    "other"              = "gray"
  ),
  ### ---- Level_3ACM ----
  Level_3ACM = c(
    # Neural
    "Astrocyte TE/NT"       = "#82AC7C",
    "Astrocyte R" = "#B8D53D",
    "OPC/COP/OLG"          = "#0065A2",
    "Neural"             = "#3CB371",
    "Choroid"            = "#FFD700",
    "Ependymal"          = "#FFEA70",

    # Immune
    "Microglia"          = "#ED8A68",
    "Macrophage"         = "#CC1212",
    "Dendritic"          = "#FF69B4",
    "Neutrophil"        = "#800080",
    "NKTB"               = "#87CEEB",

    # Vascular/Stromal
    "Fibroblast"         = "#E4D938",
    "Mural"              = "#DF8008",
    "Endothelial"        = "#60330D",

    # Neoplastic & Ambiguous
    "Neoplastic"         = "#F4B7C7",
    "Ambiguous"          = "gray",
    "other"              = "gray"
  ),
  ### ---- Level_4 ----

  Level_4 = c(

    # Non-neoplastic
    "Astrocyte TE"       = "#0B3800",
    "Astrocyte NT"       = "#82AC7C",
    "Astrocyte TE/NT"       = "#82AC7C",
    "Astrocyte R" = "#B8D53D",
    "OPC/COP/OLG"         = "#7D3C98",   # PURPLE - kept clearly distinct from the neoplastic OPC/COP/Bulk blues
    "Neural"              = "#3CB371",
    "Ependymal"           = "#FFEA70",

    # Immune
    "Microglia"           = "#ED8A68",
    "Macrophage"          = "#CC1212",
    "Dendritic"           = "#FF69B4",
    "Neutrophil"          = "#800080",
    "NKTB"                = "#87CEEB",

    # Vascular / Stromal
    "Choroid"             = "#FFD700",
    "Endothelial"         = "#60330D",
    "Mural"               = "#DF8008",
    "Fibroblast"          = "#E4D938",

    # Neoplastic
    "Neopl-Bulk"        = "#00BFFF",   # Tumor I
    "Neopl-CC-I"        = "#F8BBD0",   # Tumor II
    "Neopl-CC-II"          = "#F48FB1",   # CC-I
    "Neopl-CC-III"          = "#D5006D",   # CC-II
    #"Neopl-Bulk-B"        = "#0091B0",   # Tumor II
    "Neopl-OPC"           = "#0065A2",   # Tumor OPC A
    "Neopl-COP"           = "#00B0BA",   # Tumor OPC B
    "Neopl-NC"           = "#50CF57",   # Tumor NC-sig
    "Neopl-ACR"           = "#F86814",
    "Neopl-ECM"           = "#F5A623",
    "Neopl-RNA-low"       = "#AED6F1",

    # Other
    "Ambiguous"       = "#B3B3B3",
    "other"           = "gray"

  ),
  Level_4ACM = c(

    # Non-neoplastic
    "Astrocyte TE/NT"       = "#82AC7C",
    "Astrocyte R" = "#B8D53D",
    "OPC/COP/OLG"         = "#7D3C98",   # PURPLE - kept clearly distinct from the neoplastic OPC/COP/Bulk blues
    # The 14-type collapse names as first-class keys (spatial maps): tuned so
    # Myeloid is NOT close to Neopl-ACR (orange); leaves the finer Microglia/Dendritic keys untouched
    "Myeloid"             = "#D24E8F",   # darker pink, kept distinct from Microglia salmon and Neopl-ACR orange
    "Immune-other"        = "#F7B6D2",   # light pink
    "Neural"              = "#3CB371",
    "Ependymal"           = "#FFEA70",

    # Immune
    "Microglia"           = "#ED8A68",
    "Macrophage"          = "#CC1212",
    "Dendritic"           = "#FF69B4",
    "Neutrophil"          = "#800080",
    "NKTB"                = "#87CEEB",

    # Vascular / Stromal
    "Choroid"             = "#FFD700",
    "Endothelial"         = "#60330D",
    "Mural"               = "#DF8008",
    "Fibroblast"          = "#E4D938",

    # Neoplastic
    "Neopl-Bulk"        = "#00BFFF",   # Tumor I
    "Neopl-CC-I"        = "#F8BBD0",   # Tumor II
    "Neopl-CC-II"          = "#F48FB1",   # CC-I
    "Neopl-CC-III"          = "#D5006D",   # CC-II
    "Neopl-OPC"           = "#0065A2",   # Tumor OPC A
    "Neopl-COP"           = "#00B0BA",   # Tumor OPC B
    "Neopl-NC"           = "#50CF57",   # Tumor NC-sig
    "Neopl-ACR"           = "#F86814",
    "Neopl-ECM"           = "#F5A623",
    "Neopl-RNA-low"       = "#AED6F1",

    # Other
    "Ambiguous"       = "#B3B3B3",
    "other"           = "gray"

  ),
  TransferLabels = c(

    # Non-neoplastic
    "Astrocyte TE"       = "#0B3800",
    "Astrocyte NT"       = "#82AC7C",
    "Astrocyte TE/NT"       = "#82AC7C",
    "Astrocyte R"         = "#B8D53D",
    "OPC/COP/OLG"         = "#7D3C98",   # PURPLE - kept clearly distinct from the neoplastic OPC/COP/Bulk blues
    "Neural"              = "#9DD8E2",
    "Ependymal"           = "#FFF4B3",

    # Immune
    "Microglia"           = "#ED8A68",
    "Macrophage"          = "#CC1212",
    "Myeloid"             = "#CC1212",
    "Dendritic"           = "#FF69B4",
    "Neutrophil"          = "#800080",
    "NKTB"                = "#87CEEB",

    # Vascular / Stromal
    "Choroid"             = "#FFD700",
    "Endothelial"         = "#60330D",
    "Mural"               = "#DF8008",
    "Endothelial/Mural"   = "#60330D",
    "Fibroblast"          = "#C4B48A",

    # Neoplastic
    "Neopl-Bulk"        = "#00BFFF",   # Tumor I
    "Neopl-CC"        = "#D5006D",   # Tumor II
    "Neopl-OPC-COP"           = "#0065A2",   # Tumor OPC A
    "Neopl-NC"           = "#7AD9B1",   # Tumor NC-sig
    "Neopl-ACR"           = "#F86814",
    "Neopl-ECM"           = "#F5A623",
    "Neopl-RNA-low"       = "#AED6F1",

    # Other
    "Ambiguous"       = "#B3B3B3",
    "other"           = "#E8E9DD"

  ),


  ### ---- allen_celltype ----
  allen_celltype  = c(
    # Astrocytes
    "318 Astro-NT NN"        = "#0B3800",  # dark green
    "319 Astro-TE NN"        = "#82AC7C",  # light green
    "320 Astro-OLF NN"       = "#B8D53D",  # lime/astro reactive

    # Ependymal lineages
    "321 Astroependymal NN"  = "#FFEA70",  # ependymal yellow
    "322 Tanycyte NN"        = "#FFD700",  # golden
    "323 Ependymal NN"       = "#FFEA70",
    "324 Hypendymal NN"      = "#FFE699",  # pale yellow

    # Choroid
    "325 CHOR NN"            = "#FFD700",  # choroid plexus gold

    # OPC / Oligodendrocyte lineage
    "326 OPC NN"             = "#1BCFD5",  # cyan
    "327 Oligo NN"           = "#9013FE",  # purple
    "328 OEC NN"             = "#BD10E0",  # magenta (olfactory ensheathing)

    # Other non-neuronal niche
    "329 ABC NN"             = "#E4D938",  # fibroblast CAF yellow
    "330 VLMC NN"            = "#A52A2A",  # fibroblast meningeal brown
    "331 Peri NN"            = "#CCBAAC",  # pericyte beige
    "332 SMC NN"             = "#D4AC0D",  # vascular smooth muscle

    # Endothelial / vascular
    "333 Endo NN"            = "#60330D",  # endothelial brown
    "334 Microglia NN"       = "#ED8A68",  # microglia orange
    "335 BAM NN"             = "#800B0B",  # macrophage BAM red

    # Myeloid / lymphoid
    "336 Monocytes NN"       = "#CC1212",  # red
    "337 DC NN"              = "#FF69B4",  # pink
    "338 Lymphoid NN"        = "#6EB4D1",  # light blue

    # Neurons
    "Dopa"                   = "#FFD700",  # golden dopamine
    "GABA"                   = "#174BDA",  # blue (matches OPC B)
    "Glut"                   = "#3CB371"   # medium green
  ),

  ## ---- GB Map annotations ----
  annotation_level_1 = c(
    "Neoplastic" = "#5C3317",
    "Non-neoplastic" = "#3CB371"
  ),

  annotation_level_2 = c(
    # neoplastic
    "Differentiated-like" = "#F5DEB3",
    "Stem-like" = "#FF4500",

    "Glial-Neuronal" = "#FFD700",
    "Myeloid" = "#6B8E23",
    "Lymphoid" = "#87CEEB",
    "Vascular" = "#800000"
  ),

  annotation_level_3 = c(
    # Neoplastic
    ## Differentiated-like"
    "AC-like" = "#F5DEB3",
    "NPC-like" = "#20B2AA",
    "MES-like" = "#DDA0DD",
    ## stem-like
    "OPC-like" = "#F5DEB3",

    # Non-neoplastic
    ## Glial-Neuronal
    "Astrocyte" = "#FFD700",
    "OPC" = "#1E90FF",
    "Oligodendrocyte" = "#A000D6",
    "RG" = "#FF7F50",
    "Neuron" = "#87CEEB",
    ## Myeloid
    "Mono" = "#6B8E23",
    "TAM-MG" = "#20B2AA",
    "TAM-BDM" = "#DA70D6",
    "DC" = "#FF4500",
    "Mast" = "#800080",
    ## Lymphoid
    "B cell" = "#FFA500",
    "Plasma B" = "#DC143C",
    "CD4/CD8" = "#FF8C00",
    "NK" = "#A9A9A9",
    ## Vascular
    "Mural cell" = "#FF4500",
    "Endothelial" = "#800000"
  ),

  # GB MAP CORE level 4
  annotation_level_4 = c(
    # Neoplastic
    ## Differentiated-like
    ### AC-like
    "AC-like" = "#8B0000",
    "AC-like Prolif" = "#8B4513",
    ### NPC like
    "NPC-like OPC" = "#F5DEB3",
    "NPC-like Prolif" = "#FFA500",
    "NPC-like neural" = "#FF8C00",
    ### MES-like
    "MES-like hypoxia independent" = "#DDA0DD",
    "MES-like hypoxia/MHC" = "#DA70D6",
    ## Stem-like
    ### OPC-like
    "OPC-like" = "#FF97C4",
    "OPC-like Prolif" = "#FF478C",

    # Non-Neoplastic
    ## Glial-Neuronal
    ### Astrocyte
    "Astrocyte" = "#FFD700",
    ### OPC
    "OPC" = "#FF00C3",
    ### Oligodendrocyte
    "Oligodendrocyte" = "#7700FF",
    ### RG
    "RG" = "#FF7F50",
    ### Neuron
    "Neuron" = "#ADFF2F",

    ## Myeloid
    ### Mono
    "Mono anti-infl" = "#006400",
    "Mono hypoxia" = "#6B8E23",
    "Mono naive" = "#32CD32",  # Changed to avoid duplicate
    ### TAM-MG
    "TAM-MG aging sig" = "#556B2F",  # Changed to avoid duplicate
    "TAM-MG pro-infl I" = "#A9A9A9",
    "TAM-MG pro-infl II" = "#696969",
    "TAM-MG prolif" = "#20B2AA",
    ### TAM-BDM
    "TAM-BDM MHC" = "#800080",
    "TAM-BDM anti-infl" = "#228B22",  # Changed to avoid duplicate
    "TAM-BDM hypoxia/MES" = "#808000",  # Changed to avoid duplicate
    "TAM-BDM INF" = "#BA55D3",  # Changed to avoid duplicate
    ### DC
    "DC1" = "#00242D",
    "DC2" = "#004F63",
    "DC3" = "#007A99",
    "cDC1" = "#1E90FF",
    "cDC2" = "#F5DEB3",
    "pDC" = "#4682B4",  # Changed to avoid duplicate
    ### Mast
    "Mast" = "#8B008B",  # Changed to avoid duplicate

    ## Lymphoid
    ### B cell
    "B cell" = "#DC143C",
    ### Plasma B
    "Plasma B" = "#FF4500",
    ### CD4/CD8
    "CD4 INF" = "#038DBE",
    "CD4 rest" = "#00CED1",  # Changed to avoid duplicate
    "CD8 EM" = "#4682B4",  # Changed to avoid duplicate
    "CD8 NK sig" = "#5F9EA0",  # Changed to avoid duplicate
    "CD8 cytotoxic" = "#05AFF2",
    "Prolif T" = "#1E90FF",
    "Reg T" = "#4B0082",
    "Stress sig" = "#DDA0DD",
    ### NK
    "NK" = "#00008B",

    ## Vascular
    ### Mural cell
    "Pericyte" = "#B22222",  # Changed to avoid duplicate
    "Perivascular fibroblast" = "#D6BB3B",
    "SMC" = "#FF4500",
    "SMC COL" = "#FF69B4",
    "SMC prolif" = "#FF1493",
    "Scavenging endothelial" = "#FF00FF",
    "Scavenging pericyte" = "#8A2BE2",
    "VLMC" = "#87CEEB",
    ### Endothelial
    "Endo arterial" = "#FF6347",  # Changed to avoid duplicate
    "Endo capilar" = "#8B4513",  # Changed to avoid duplicate
    "Tip-like" = "#3CB371"
  ),


  ## ---- GBMap annotation hierarchy ----

  "gbmap.hierarchy.annotation_level_1" = list(
    "Neoplastic",
    "Non-Neoplastic"
  ),
  "gbmap.hierarchy.annotation_level_2" = list(
    "Neoplastic" = list(
      "Differentiated-like",
      "Stem-like"
    ),
    "Non-Neoplastic" = list(
      "Glial-Neuronal",
      "Myeloid",
      "Lymphoid",
      "Vascular"
    )
  ),

  "gbmap.hierarchy.annotation_level_3" = list(
    "Neoplastic" = list(
      "Differentiated-like" = list(
        "AC-like",
        "NPC-like",
        "MES-like"
      ),
      "Stem-like" = list(
        "OPC-like"
      )
    ),
    "Non-Neoplastic" = list(
      "Glial-Neuronal" = list(
        "Astrocyte",
        "OPC",
        "Oligodendrocyte",
        "RG",
        "Neuron"
      ),
      "Myeloid" = list(
        "Mono",
        "TAM-MG",
        "TAM-BDM",
        "DC",
        "Mast"
      ),
      "Lymphoid" = list(
        "B cell",
        "Plasma B",
        "CD4/CD8",
        "NK"
      ),
      "Vascular" = list(
        "Mural cell",
        "Endothelial"
      )
    )
  ),

  "gbmap.hierarchy.annotation_level_4" = list(
    "Neoplastic" = list(
      "Differentiated-like" = list(
        "AC-like" = list("AC-like", "AC-like Prolif"),
        "NPC-like" = list("NPC-like OPC", "NPC-like Prolif", "NPC-like neural"),
        "MES-like" = list("MES-like hypoxia independent", "MES-like hypoxia/MHC")
      ),
      "Stem-like" = list(
        "OPC-like" = list("OPC-like", "OPC-like Prolif")
      )
    ),
    "Non-Neoplastic" = list(
      "Glial-Neuronal" = list(
        "Astrocyte" = list("Astrocyte"),
        "OPC" = list("OPC"),
        "Oligodendrocyte" = list("Oligodendrocyte"),
        "RG" = list("RG"),
        "Neuron" = list("Neuron")
      ),
      "Myeloid" = list(
        "Mono" = list("Mono anti-infl", "Mono hypoxia", "Mono naive"),
        "TAM-MG" = list("TAM-MG aging sig", "TAM-MG pro-infl I", "TAM-MG pro-infl II", "TAM-MG prolif"),
        "TAM-BDM" = list("TAM-BDM MHC", "TAM-BDM anti-infl", "TAM-BDM hypoxia/MES", "TAM-BDM INF"),
        "DC" = list("DC1", "DC2", "DC3", "cDC1", "cDC2", "pDC"),
        "Mast" = list("Mast")
      ),
      "Lymphoid" = list(
        "B cell" = list("B cell"),
        "Plasma B" = list("Plasma B"),
        "CD4/CD8" = list("CD4 INF", "CD4 rest", "CD8 EM", "CD8 NK sig", "CD8 cytotoxic", "Prolif T", "Reg T", "Stress sig"),
        "NK" = list("NK")
      ),
      "Vascular" = list(
        "Mural cell" = list("Pericyte", "Perivascular fibroblast", "SMC", "SMC COL", "SMC prolif", "Scavenging endothelial", "Scavenging pericyte", "VLMC"),
        "Endothelial" = list("Endo arterial", "Endo capilar", "Tip-like")
      )
    )
  ),


  ## ---- Allen ref large ----
  allen.ref.large = c(
    # Non-Neuronal
    "Astro" = "#15F1DA",
    "Endo" = "#A77C70",
    "Micro-PVM" = "#A6666F",
    "Oligo" = "#765FD2",
    "SMC-Peri" = "#71238C",
    "VLMC" = "#9A85EC",

    # Glutamatergic
    "CR" = "#00FF66",


    "CA1-ProS" = "#873C46",
    "CA2-IG-FC" = "#30E6BA",
    "CA3" = "#3842EC",
    "Car3" = "#297F98",
    "CT SUB" = "#42EC04",
    "DG" = "#8059FF",

    "L2 IT ENTl" = "#2FA4EB", # x
    "L2 IT ENTm" = "#1AD475",
    "L2-3 IT ENTl" = "#4E9E9E",

    "L2-3 IT CTX" = "#C0F27F",
    "L2-3 IT PPP" = "#21EC1D",
    "L2-3 IT RHP" = "#48CB80",
    "L3 IT ENT" = "#513577",

    "L4 RSP-ACA" = "#5BC6A1",
    "L4-5 IT CTX" = "#09CCC6",
    "L5 IT CTX" = "#4D483C", # x
    "L5 PPP" = "#0E8C30",
    "L5 PT CTX" = "#0EED60",
    "L5-6 IT TPE-ENT" = "#873C46",
    "L5-6 NP CTX" = "#BD76DC",
    "L6 CT CTX" = "#3E98A5",
    "L6 IT CTX" = "#8D8C20",
    "L6 IT ENTl" = "#C4EC04",
    "L6b CTX" = "#C0F27F",
    "L6b-CT ENT" = "#613D90",
    "NP PPP" = "#6E7C6F",
    "NP SUB" = "#503C40",
    "SUB-ProS" = "#613D90",

    # GABAergic
    "Lamp5" = "#520CB9",
    "Pvalb" = "#3F8FA1",
    "Sncg" = "#8C620E",
    "Sst" = "#ECD704",
    "Sst Chodl" = "#68621B",
    "Vip" = "#4D483C"
  ),

  allen.ref.large.level2 = c(
    # Non-Neuronal
    "Astro" = "#FF9900",
    "Endo" = "#990000",
    "Micro-PVM" = "#A6666F",
    "Oligo" = "#CC66FF",
    "SMC-Peri" = "#e31a1c",
    "VLMC" = "#9A85EC",
    "Glutaminergic" = "#00FF66",
    "GABAergic" = "#C4EC04"
  ),

  allen.ref.small = c(
    "Astro" = "#15F1DA",
    "CR" = "#00FF66",
    "Endo" = "#A77C70",
    "L2-3 IT" = "#C0F27F",
    "L4" = "#5BC6A1",
    "L5 IT" = "#4D483C",
    "L5 PT" = "#0EED60",
    "L6 CT" = "#3E98A5",
    "L6 IT" = "#8D8C20",
    "L6b" = "#C0F27F",
    "Lamp5" = "#520CB9",
    "Macrophage" = "#A6666F",
    "Meis2" = "#68621B",
    "NP" = "#503C40",
    "Oligo" = "#765FD2",
    "Peri" = "#71238C",
    "Pvalb" = "#3F8FA1",
    "Serpinf1" = "#ECD704",
    "SMC" = "#71238C",
    "Sncg" = "#8C620E",
    "Sst" = "#ECD704",
    "Vip" = "#4D483C",
    "VLMC" = "#9A85EC"
  ),

  ## ---- UCSC gbm data ----
  gbm.ucsc.Level2 = c(
    "Tumor cell" = "#8B0000",          # Inspired by "AC-like"

    "Astrocyte" = "#FFD700",          # Inspired by "Astrocyte" from previous palette
    "Oligodendrocyte" = "#7700FF",    # Inspired by "Oligodendrocyte"
    "OPC" = "#FF00C3",                # Inspired by "OPC"
    "Neuron" = "#ADFF2F",             # Inspired by "Neuron"

    "Microglia" = "#006400",          # Inspired by "Microglia-H-6a"
    "Macrogphage" = "#4682B4",        # Inspired by "Macrophage-TA-4d"
    "DC" = "#FF4500",                 # Inspired by "Dendritic-17a"
    "Monocyte" = "#228B22",           # Inspired by "Microglia-H-6b"

    "Granulocyte" = "#00BFFF",        # Inspired by "Macrophage-TA-4c"
    "Lymphocyte" = "#1E90FF",         # Inspired by "Macrophage-TA-4b"
    "RBC" = "#DC143C",                # Inspired by "B cell"

    "Endothelial" = "#8A2BE2",        # Inspired by "Scavenging pericyte"
    "gbmEndo" = "#FF00FF",            # Inspired by "Scavenging endothelial"

    "gbmEndoPeri" = "#3CB371",        # Inspired by "Microglia-TA-5b"
    "gbmPeri" = "#2E8B57",            # Inspired by "Microglia-TA-5c"
    "Pericyte" = "#B22222",           # Inspired by "Pericyte"

    "Fibroblast" = "#D6BB3B",         # Inspired by "Perivascular fibroblast"
    "gbmFib" = "#FF69B4",             # Inspired by "SMC COL"
    "SMC" = "#FF4500"                # Inspired by "SMC"
  ),

  gbm.ucsc.celltype = c(
    # Tumor cells
    "TC_AC" = "#8B0000",          # Inspired by "AC_like"
    "TC_mesh" = "#8B4513",        # Inspired by "AC_like Prolif"
    "TC_mesnh" = "#DDA0DD",       # Inspired by "MES_like hypoxia independent"
    "TC_MTC" = "#DA70D6",         # Inspired by "MES_like hypoxia/MHC"
    "TC_NPC" = "#FFA500",         # Inspired by "NPC_like Prolif"
    "TC_oligo" = "#FF8C00",       # Inspired by "NPC_like neural"
    "TC_OPC" = "#FF97C4",         # Inspired by "OPC_like"
    "TC_prolif" = "#FF478C",      # Inspired by "OPC_like Prolif"

    # glial
    "AstroAQP1" = "#FFD700",      # Inspired by "Astrocyte"
    "AstroPLCG1" = "#FFD700",     # Same as "AstroAQP1" to keep consistency
    "AstroSERPINI2" = "#FF97C4",  # Inspired by "OPC_like"

    "Oligo_2_1" = "#7700FF",      # Inspired by "Oligodendrocyte"
    "Oligo_2_2" = "#7700FF",      # Same as "Oligo_2_1" to keep consistency
    "Oligo_2_3_1" = "#7700FF",    # Same as "Oligo_2_1" to keep consistency
    "Oligo_2_3_2" = "#7700FF",    # Same as "Oligo_2_1" to keep consistency
    "OPC" = "#FF00C3",            # Inspired by "OPC"

    "Neurons_Ex" = "#ADFF2F",     # Inspired by "Neuron"
    "Neurons_In" = "#32CD32",     # Inspired by "Mono naive"

    # lmphocytes
    "Bcell" = "#DC143C",          # Inspired by "BCell_26a"
    "Bcell_plasma" = "#B22222",   # Inspired by "BCell_26b"
    "RBC" = "#DC143C",            # Inspired by "B cell"
    "CD4" = "#038DBE",            # Inspired by "CD4 INF"
    "CD8" = "#05AFF2",            # Inspired by "CD8 cytotoxic"
    "CD8_TRM" = "#1E90FF",        # Inspired by "Prolif T"
    "Tcell_prolif" = "#1E90FF",   # Inspired by "Prolif T"
    "Treg" = "#4B0082",           # Inspired by "Reg T"
    "Granulocyte" = "#00BFFF",    # Inspired by "Macrophage_TA_4c"
    "NKT" = "#00008B",            # Inspired by "NKT_TAH_14a"

    # myeloid
    "DC" = "#FF4500",             # Inspired by "Dendritic_17a"
    "DC_mig" = "#FF6347",         # Inspired by "Dendritic_17b"
    "Mac_5_1_2" = "#4682B4",      # Inspired by "Macrophage_TA_4d"
    "Mac_5_1_3" = "#5F9EA0",      # Inspired by "Macrophage_TA_4e"
    "Mac_prolif" = "#6495ED",     # Inspired by "Macrophage_TA_4f"
    "Mg_1_1" = "#006400",         # Inspired by "Microglia_H_6a"
    "Mg_1_3" = "#228B22",         # Inspired by "Microglia_H_6b"
    "Mg_prolif" = "#3CB371",      # Inspired by "Microglia_TA_5b"
    "Monocyte" = "#228B22",       # Inspired by "Microglia_H_6b"
    "Monocyte_Nc" = "#006400",    # Inspired by "Microglia_H_6a"

    # endo
    "Art1" = "#FF4500",           # Inspired by "Dendritic_17a"
    "Art2" = "#FF6347",           # Inspired by "Dendritic_17b"
    "Art3" = "#FF7F50",           # Inspired by "Dendritic_17c"
    "Endo_cap" = "#8A2BE2",       # Inspired by "Scavenging pericyte"
    "gbmEndo_9_1" = "#FF00FF",    # Inspired by "Scavenging endothelial"
    "gbmEndo_9_1a" = "#3CB371",   # Inspired by "Microglia_TA_5b"
    "gbmEndo_9_3" = "#2E8B57",    # Inspired by "Microglia_TA_5c"
    "gbmEndo_9_4" = "#66CDAA",    # Inspired by "Microglia_TA_10c"
    "Venule" = "#87CEEB",          # Inspired by "Macrophage_TA_4a"
    "Venous" = "#8A2BE2",         # Inspired by "Scavenging pericyte"

    "Pericyte" = "#B22222",       # Inspired by "Pericyte"
    "gbmEndoPeri_1" = "#8FBC8F",  # Inspired by "Microglia_TA_10d"
    "gbmEndoPeri_2" = "#006400",  # Inspired by "Microglia_H_6a"
    "gbmPeri" = "#2E8B57",        # Inspired by "Microglia_TA_5c"

    # fibroblast
    "FBMC" = "#D6BB3B",           # Inspired by "Perivascular fibroblast"
    "Fibroblast" = "#D6BB3B",     # Same as "FBMC" to keep consistency
    "gbmFib_4" = "#FF69B4",       # Inspired by "SMC COL"
    "gbmFib_5" = "#FF1493",       # Inspired by "SMC prolif"
    "gbmFib_6" = "#FF00C3",       # Inspired by "OPC"
    "SMC" = "#FF4500"            # Inspired by "SMC"

  )
)

# --- wrappers ----

## Wrapper for running FGSEA analysis of all clusters in a seurat object
## Input is a seurat object with pre-defined clusters and a list of gene sets
# Run FindAllaMarkers
# Then run fgsea::fgseaMultilevel from the fgea package using curated gene list with relevant cell types/signatures

.wrapper_fgsea <- function(
    seurat.object, # pre-filtered data
  cluster.label, # vector
  gene.sets.list,
  results.path,
  analysis.name
){
  require(fgsea)
  require(Seurat)
  require(tidyr)
  require(dplyr)

  stopifnot(analysis.name != "")
  stopifnot(file.exists(results.path))

  ## set ident
  seurat.object <- Seurat::SetIdent(
    seurat.object, value = cluster.label
  )
  print(table(Seurat::Idents(seurat.object)))
  seurat.object <- seurat.object %>%
    Seurat::NormalizeData() %>%
    Seurat::FindVariableFeatures(nfeatures = Inf)

  ## Run DGE
  cat(" ... Finding DGEs for each cluster")
  DGE_table <- Seurat::FindAllMarkers(
    seurat.object,
    logfc.threshold = 0,
    test.use = "wilcox",
    min.pct = 0.1,
    min.diff.pct = 0,
    only.pos = TRUE,
    # max.cells.per.ident = 25,
    return.thresh = 1
    #assay = "RNA"
  )
  DGE_list <- split(DGE_table, DGE_table$cluster)
  # remove cluster if no DEG
  u <- unlist(lapply(DGE_list, nrow))==0
  if(any(u)){
    cat(" ... removing cluster with no DEGs: ",names(DGE_list)[u])
    DGE_list <- DGE_list[which(u == FALSE)]
  }
  # run fgsea simple for each of the clusters in the list
  cat(" ... Running fgsea on DGEs for each cluster")

  res <- lapply(
    DGE_list, function(x) {
      gene_rank <- setNames(x$avg_log2FC, x$gene)
      fgseaRes <- fgsea::fgseaMultilevel(
        pathways = gene.sets.list,
        scoreType = "pos",
        eps = 0,
        stats = gene_rank
      )
      return(fgseaRes)
    })
  names(res) <- names(DGE_list)
  #names(res2) <- names(DGE_list)

  res <- lapply(res, function(x) {
    x[x$padj < 0.01, ]
  })
  res <- lapply(res, function(x) {
    x[x$size > 3, ]
  })
  res <- lapply(res, function(x) {
    x[order(x$padj, decreasing = F), ]
  })

  ## Write results to csv
  file_name <- file.path(
    results.path,
    paste0("fgseaWrapper-res-",cluster.label,"-",analysis.name,".tsv")
  )
  write.append <- FALSE
  for ( i in 1:length(res)){
    # cat(paste0("\n",names(res)[i]),file=file_name, sep = "\t", append = T)
    if( nrow(res[[i]]) >0){
      res_df <- data.frame(res[[i]])
      res_df <- cbind(
        data.frame(
          snn_cluster = rep(names(res)[i]),
          nrow(res_df)), res_df)
      res_df$leadingEdge <- unlist(
        lapply(res_df$leadingEdge,
          function(x) paste(x, collapse = "|")))
      # cnames <- ifelse(i==1, TRUE, FALSE)
      # app <- ifelse(i==1, FALSE, TRUE)
      write.table(
        res_df,
        file=file_name,
        row.names=FALSE,
        col.names = !write.append,
        append = write.append,
        sep="\t",
        quote=F
      )
      write.append <- TRUE
    }
  }

} # end fgsea wrapper


.wrapper_fgsea_plot <- function(fgsea.dir, annotation, seurat.object, analysis.name){
  # raed output
  file.name <- file.path(
    fgsea.dir,  paste0("fgseaWrapper-res-",annotation,"-",analysis.name,".tsv"))
  if (!file.exists(file.name)){
    cli::cli_alert("No FGSEA tsv results file found: {.file {file.name}}")
    next
  }

  fgsea.df <- read.delim(file = file.name, sep='\t')
  str(fgsea.df)

  ## order annotations on pathway NES similarity (hclust)
  matrix.df <- fgsea.df %>%
    dplyr::select(snn_cluster, pathway, NES) %>%
    spread(key = snn_cluster, value = NES, fill = 0)

  # Set the row names to the celltypes column
  rownames(matrix.df) <- matrix.df$pathway
  # Remove the celltypes column
  matrix.df <- select(matrix.df, -pathway)
  str(matrix.df)
  # Calculate the distance matrix
  dist.mat <- stats::dist(t(matrix.df), method = "euclidean")
  dist.mat.2 <- stats::dist(matrix.df, method = "euclidean")

  ## select top 5 fgsea NES results for each cluster
  top.n <- 5
  grep("microglia", fgsea.df$pathway, ignore.case = T, value = T)

  my.df <- fgsea.df %>%
    dplyr::select(snn_cluster, pathway, NES, padj, size) %>%
    group_by(snn_cluster) %>%
    top_n(top.n, NES) %>%
    ungroup() %>%
    arrange(snn_cluster, desc(NES))


  # plot
  # Create a data frame of tiles
  tiles.df <- expand.grid(
    snn_cluster = unique(my.df$snn_cluster),
    pathway = unique(my.df$pathway))

  ## add UMAP for reference
  seurat.object <- Seurat::SetIdent(seurat.object, value = annotation)
  table(Idents(seurat.object), useNA = "always")

  p2 <- SCpubr::do_DimPlot(
    colors.use =
      .colorPalDiscrete(sort(unique(seurat.object[[annotation]][,1]))),
    # legend.position = "none",
    # cells.highlight = cell.barcodes,
    seurat.object,
    label = FALSE,
    plot_cell_borders = FALSE
    #plot_marginal_distributions = TRUE
  ) + ggtitle(paste0(active.ident))


  # Create the plot
  p <- ggplot() +
    geom_tile(data = tiles.df,
      aes(x = as.character(snn_cluster), y = pathway),
      color = "black", fill = NA, width = 1, height = 1) +
    geom_point(data = my.df,
      aes(x = snn_cluster, y = pathway, fill = padj, size = NES),
      shape = 21, colour = "black", stroke = 0.5) +
    scale_fill_gradientn(colors = RColorBrewer::brewer.pal(9, "YlOrRd")) +
    theme_minimal() +
    # scale_x_discrete(position = "top") +
    labs(
      fill = "Adj pvalue",
      title = paste0("fGSEA\n", annotation),)

  # Modify the theme
  # p <- p + theme(
  #   axis.text.x = element_text(size = rel(1),  vjust = -45),
  #   axis.ticks.x = element_blank(),
  #   panel.grid.major = element_blank(),
  #   panel.grid.minor = element_blank(),
  #   panel.background = element_blank(),
  #   plot.margin = margin(5, 1, 1, 1, "cm")
  # )
  p <- p + theme(axis.text.x = element_text(angle = 90))

  # set a file name of pdf and plot to pdf

  plot(p2 + p)
  rm(p, p2, fgsea.df, my.df, tiles.df, dist.mat, dist.mat.2, matrix.df)
  #}
  #}

}




.wrapper_mca <- function(
    results.dir, clustersets, sd
){
  require(scMCA)
  cli::cli_h2("Running MCA")
  # mca
  # mca.dir <- results.dir # file.path(results.dir, "mca_results")
  # dir.create(mca.dir, showWarnings = FALSE, recursive = TRUE)


  # clustersets <- grep("snn", colnames(clustersets.table), value = TRUE)
  # datasets <- c("sd.subset", "sd_raw_sub")
  # dataset <- c("seurat.reference")
  # cli::cli_alert("Running MCA: {.var {dataset}}")
  # sd <- get(dataset)
  normalized_data <- LayerData(object = sd, assay="RNA", layer = "data")
  head(rownames(normalized_data))
  str(normalized_data)
  #cli::cli_alert("Running MCA: {.var {dataset}}")
  mca_result <- scMCA::scMCA(scdata = normalized_data, numbers_plot = 3)

  # save mca_result as rds object
  saveRDS(mca_result, file = file.path(results.dir,"mca_results.rds"))
}


.wrapper_mca_plot <- function(
    results.dir, clustersets,
  sd.ref
){
  mca.dir <- file.path(results.dir, "mca_results")

  require(scMCA)
  mca_result <- readRDS(file = file.path(results.dir,"mca_results.rds"))
  fname <- paste0("09b-MCA_results_DotPlots.pdf")
  file.name <- file.path(results.dir, fname)
  cli::cli_alert("Plotting to file: {.file {file.name}}")

  cli::cli_h2("Running MCA per clusterset")
  # mca


  fname <- paste0("09b-MCA_results_DotPlots.pdf")
  file.name <- file.path(results.dir, fname)
  cli::cli_alert("Plotting to file: {.file {file.name}}")

  pdf(file.name, width = 17, height = 17)

  # Plot results for clustersset annotations
  # clustersets <- c("snn025_sd.subset","ident_w_subclusters")
  for ( clusterset in clustersets ){
    # set active ident
    Idents(object = sd.ref) <- clusterset
    table(Idents(sd.ref))
    mca.df <- data.frame(
      Cell = colnames(sd.ref),
      cluster = Idents(sd.ref))

    # get the mca result
    df <- mca_result$scMCA_probility
    str(df)

    # merge the two dataframes on cell name (get Cluster assigmnent for all cells)
    df <- base::merge(df, mca.df, by = "Cell")
    head(df)

    # get min max and median for the df$score column
    # hist(df$Score)

    # loop all unique clusters in mca results df
    # find all unique Cell types in cluster and calculate fraction

    celltypes.list <- list()
    for (i in unique(df$cluster)){
      # n.cells
      n.cells <- length(unique(df[df$cluster == i,]$Cell))
      # get all cell types in cluster
      celltypes <- unique(df[df$cluster == i,]$`Cell type`)
      celltypes.list[[i]] <- data.frame(cluster=i, celltypes=celltypes, frac=NA)
      # calulate fraction of cells in cluster that are of cell type j
      celltypes.list[[i]]$frac <-  sapply(celltypes, function(x) {
        nrow(df[df$cluster == i & df$`Cell type` == x,])/n.cells
      })
      # sort table on frac column
      celltypes.list[[i]] <- celltypes.list[[i]] [ order(celltypes.list[[i]]$frac, decreasing = T),]
    }
    # catenate all celltypes.list[[i]] into one dataframe and save as csv
    celltypes.df <- do.call(rbind, celltypes.list)
    str(celltypes.df) # 802 obs.

    # Histogrram of the cettype assignment fractions
    # hist(celltypes.df$frac)
    # how many unique cell types in list
    length(unique(celltypes.df$celltypes)) #  311)
    # write to csv
    if(length(unique(celltypes.df$celltypes)) == 0 ){
      cli::cli_alert("No cell types identified")
      next
    }
    write.csv(
      celltypes.df,
      file = file.path(results.dir,paste0(clusterset, "-mca_results-celltypes.csv"))
    )

    my.df <- celltypes.df[celltypes.df$frac > 0.1,]
    str(my.df)
    # how many unique cell types in list
    if( nrow( my.df ) == 0 ){
      cli::cli_alert("No cell types with fraction > 0.1")
      next
    }
    length(unique(my.df$celltypes)) #  83 left after filtering


    # ----- build wide matrix from the FILTERED df and cluster rows + cols -----

    # write filtered table (you were writing the unfiltered one by mistake)
    write.csv(
      my.df,
      file = file.path(results.dir, paste0(clusterset, "-mca_results-celltypes-filtered.csv"))
    )

    # Wide matrix for clustering (use filtered my.df to match what you plot)
    matrix.df <- my.df %>%
      tidyr::spread(key = cluster, value = frac, fill = 0)

    # keep rownames as cell types, then drop the label column
    rownames(matrix.df) <- matrix.df$celltypes
    matrix.df <- dplyr::select(matrix.df, -celltypes)

    # Cluster columns (clusters)
    dist_cols <- stats::dist(t(matrix.df), method = "euclidean")
    hc_cols   <- stats::hclust(dist_cols, method = "average")
    clusters.ordered <- colnames(matrix.df)[hc_cols$order]

    # Cluster rows (cell types)
    dist_rows <- stats::dist(matrix.df, method = "euclidean")
    hc_rows   <- stats::hclust(dist_rows, method = "average")
    celltypes.ordered <- rownames(matrix.df)[hc_rows$order]

    # Apply factor levels to match clustered orders
    my.df$cluster   <- factor(my.df$cluster,   levels = clusters.ordered)
    my.df$celltypes <- factor(my.df$celltypes, levels = rev(celltypes.ordered))  # reverse if you prefer "top-heavy" heatmaps

    # Tiles grid must use the SAME levels to preserve order
    tiles.df <- tidyr::expand_grid(
      cluster   = factor(levels(my.df$cluster),   levels = levels(my.df$cluster)),
      celltypes = factor(levels(my.df$celltypes), levels = levels(my.df$celltypes))
    )

    ## add UMAP for reference
    sd.ref <- Seurat::SetIdent(sd.ref, value = clusterset)
    table(Idents(sd.ref), useNA = "always")

    p2 <- SCpubr::do_DimPlot(
      sd.ref,
      colors.use =
        .palette_discrete_bp(sort(unique(sd.ref[[clusterset]][,1]))),
      label = FALSE,
      plot_cell_borders = FALSE
      #plot_marginal_distributions = TRUE
    ) + ggtitle(paste0(clusterset))


    # Create the plot
    p <- ggplot() +
      geom_tile(data = tiles.df,
        aes(x = cluster, y = celltypes),
        color = "black", fill = NA, width = 1, height = 1) +
      geom_point(data = my.df,
        aes(x = cluster, y = celltypes, fill = frac),
        shape = 21, colour = "black", stroke = 0.5, size = 3) +
      scale_fill_gradientn(colors = RColorBrewer::brewer.pal(9, "YlOrRd")) +
      theme_minimal() +
      scale_x_discrete(position = "top") +
      labs(
        fill = "Fraction of \ncells in clusterset",
        title = paste0("MCA assignment\n", clusterset))
    p <- p + theme(axis.text.x = element_text(angle = 90))


    plot(p2 + p)

  } # end clusterset  MCA loop
  dev.off()
} # end wrapper mca

###############
#' @description #' wrapper for looping subclustering
#' @param seurat.reference seurat object, can be same as sd_raw_sub, but also downsampled data etc.
#' Gene expression DEGs will be calculated on the sd_raw_sub data
#' @param subset.name Name of the input data set
#' @param umap.subset name of reference umap from the sd.subset object
#' the cluster to be subclustered
#' @param active.ident Name of the active.ident column in the metadata
#' @param results.main.dir main output directory

#' @param seurat.reference seurat object, reference data to compare against when comparing ALL input samples vs rest
#' @param umap.reference name of reference umap from the seurat.reference object
#' @param seurat.full seurat object, reference data to compare against when comparing ALL input samples vs rest
#' @param umap.full name of reference umap from the seurat.reference object


#' @param topN number of top genes to use for clustering
#' @param louvain.res resolution for clustering
#' @param do.liana.full logical, do liana plots on FULL data
#' @param downsample.n.heatmap number of cells to downsample for heatmap
#' @param downsample.n.liana number of cells to downsample for liana
#' @param dims vector of dimensions to use for PCA in clustering
#' @param nfeatures number of features to use for clustering
#' @return write data, results and plots to a defined output folder
#' @export
.subClusteringWrap <- function(

  sd.subset,
  subset.name,
  umap.subset,
  results.main.dir,

  seurat.reference,
  umap.reference,

  seurat.full,
  umap.full,

  dims = c(1, 30),
  nfeatures = 2000,
  active.ident,
  topN = 20,
  louvain.res = 0.25,
  k.param = 30,
  snn.prune.factor = 15,

  downsample.n.heatmap = 5000,
  min.cells.per.cluster = 15

){
  gc()
  ##
  # options(future.globals.maxSize = 4 * 1024^3)
  # library(future)
  # num_cores <- parallel::detectCores()
  # plan("multisession", workers = 3)

  #  environment
  require(Seurat)
  require(dplyr)
  require(ggplot2)
  require(Matrix)
  require(ggalluvial)
  require(cli)
  require(conflicted)
  require(glue)
  conflicted::conflicts_prefer(dplyr::filter)
  set.seed(169999)


  cli::cli_h1(" =========================== ")
  cli::cli_h1("Starting subClusteringWrap: {.var {subset.name}}")
  cli::cli_h1(" =========================== ")

  #  check input

  cli::cli_h2("Checking input")
  # results.dir must be a valid dir
  # if(!dir.exists(results.main.dir)){
  #   stop("results.main.dir does not exist")
  # }
  # subset.name and subset.name must be valid
  if(is.na(subset.name) | is.na(subset.name) | is.na(subset.name)){
    stop("subset.name and subset.name and subset.name must be defined")
  }

  # stop if umap.subset is not in Reductions(seurat.reference)
  if(!(umap.subset %in% Seurat::Reductions(seurat.reference))){
    stop("umap.subset must be in Reductions(seurat.reference)")
  }
  # stop if active.ident is not a column in the metadata of seurat.reference
  if(!(active.ident %in% colnames(seurat.reference@meta.data))){
    stop("active.ident must be a column in the metadata of seurat.reference")
  }


  #    Internal functions
  .addPrefix <- function(num) {
    if (num < 10) {
      return(paste0("0", num))
    } else {
      return(as.character(num))
    }
  }

  #    Set gene filters
  # Basic gene filter G1 for the seurat.reference
  G1 <- .CreateSeuratFilter(
    min.cells.expressed = 25,
    features.remove = c(
      gene.list$mitochondria,
      gene.list$haemoglobin,
      gene.list$rcas_flex_probes,
      gene.list$y_chr,
      "Cst3")
  )

  annot_columns <- c("sample_type","sample_id","nCount_RNA","nFeature_RNA","rcas_module_bin")


  #    Prep folders and names
  cli::cli_h2("Prep folders and names")
  analysis.name <- subset.name
  # start by creating folder structure
  out.dir <- file.path(results.main.dir, analysis.name)
  dir.create(out.dir, recursive = TRUE, showWarnings = FALSE)

  #   Idents and Subset data
  cli::cli_h2("Idents and Subset data")
  # subset the seurat object with cell.barcodes
  seurat.reference <- Seurat::SetIdent(seurat.reference, value = active.ident)
  sd.subset <- Seurat::SetIdent(sd.subset, value = active.ident)
  table(Idents(seurat.reference))
  table(Idents(sd.subset))

  #  Print params
  params.file <- file.path(out.dir, "00-params.txt")
  input.params <- list(

    subset.name = subset.name,
    analysis.name = analysis.name,
    sd.subset = glue("ngenes: {dim(sd.subset)[1]}, ncells: {dim(sd.subset)[2]}"),
    umap.subset = umap.subset,

    seurat.reference = glue("ngenes: {dim(seurat.reference)[1]}, ncells: {dim(seurat.reference)[2]}"),
    umap.reference = umap.reference,

    results.main.dir = results.main.dir,
    out.dir = out.dir,

    #reduction.names = reduction.names,
    dims = glue("{dims[1]}, {dims[2]}"),
    nfeatures = nfeatures,
    active.ident = active.ident,
    topN = topN,

    downsample.n.heatmap = downsample.n.heatmap,
    min.cells.per.cluster = min.cells.per.cluster,
    k.param = k.param,
    louvain.res = louvain.res
  )

  # Function to pretty print parameters
  .pretty_print_params <- function(params, file_path) {
    # Open a connection to the file
    file_conn <- file(file_path, "w")

    # Write the parameters to the file
    for (param_name in names(params)) {
      param_value <- params[[param_name]]
      if (is.null(param_value)) {
        param_value <- "NULL"
      } else if (is.vector(param_value)) {
        param_value <- paste(param_value, collapse = ", ")
      }
      cat(sprintf("%s: %s\n", param_name, param_value), file = file_conn)
    }

    # Close the file connection
    close(file_conn)
  }
  # Call the function to pretty print parameters to the file
  pretty_print_params(input.params, params.file)
  # Verify the file contents
  cat(readLines(params.file), sep = "\n")



  #  Add subcluster assignments to metadata
  # set idents for the full datasets `input_subset (subset.name)
  cli::cli_alert("Adding subcluster assignments to metadata")


  seurat.reference <- SetIdent(seurat.reference, value = "other")
  seurat.reference <- SetIdent(seurat.reference, cells = colnames(sd.subset), value = subset.name)
  seurat.reference[["input_subset"]] <- Idents(seurat.reference)
  table(seurat.reference[["input_subset"]])

  sd.subset <- SetIdent(sd.subset, value = "other")
  sd.subset <- SetIdent(sd.subset, cells = colnames(sd.subset), value = subset.name)
  sd.subset[["input_subset"]] <- Idents(sd.subset)
  table(sd.subset[["input_subset"]])

  seurat.full <- SetIdent(seurat.full, value = "other")
  seurat.full <- SetIdent(seurat.full, cells = colnames(sd.subset), value = subset.name)
  seurat.full[["input_subset"]] <- Idents(seurat.full)
  table(seurat.full[["input_subset"]])


  #  FindAllMarkers Input
  cli::cli_h2("FindAllMarkers - Input datasets")
  # set ident of seurat.reference and sd_raw to the input_subset
  datasets <- c("seurat.reference")

  # loop the full data sets to get DEG genes using FindAllMarkers
  genelists_findallmarkers <- list()
  for (dataset in datasets){

    sd <- get(dataset)
    sd <- sd %>%
      Seurat::SetIdent(value = "input_subset") %>%
      NormalizeData() %>%
      FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
      ScaleData()
    table(Idents(sd))

    genelists_findallmarkers[[dataset]] <- FindAllMarkers(
      sd,
      min.pct = 0.5,
      logfc.threshold = 1,
      only.pos = TRUE,
      verbose = T)

    # sort on avg_log2FC
    deg_df <- genelists_findallmarkers[[dataset]] %>%
      dplyr::filter(avg_log2FC > 1) %>%
      dplyr::arrange( desc(cluster), desc(avg_log2FC) )

    # save  deg_df as csv
    fname <- file.path(out.dir, "01-FindAllMarkers_Input_vs_Ref.csv")
    write.csv(deg_df, file = fname, row.names = FALSE)

    # create a tmp downsampled object (=fewer cells for plotting)
    sd.downsampled <- subset(x = sd, downsample = downsample.n.heatmap)
    table(Idents(sd.downsampled))

    # Plot heatmap on downsampled data since on full data takes TOO Much time ...
    plot_n_genes <- 25
    topGenes <- deg_df %>%
      # group_by(cluster) %>%
      dplyr::filter(cluster != "other") %>%
      dplyr::filter(avg_log2FC > 2) %>%
      slice_head(n = plot_n_genes) %>%
      ungroup()
    annot_by <- c(annot_columns, "input_subset")
    if(nrow(topGenes) > 5){
      p <- suppressWarnings(dittoSeq::dittoHeatmap(
        drop_levels = T,
        sd.downsampled[ topGenes$gene,],
        # filename = plot.name,
        cluster_cols = TRUE,
        # scaled.to.max=TRUE,
        complex = TRUE,
        annot.by = annot_by
      )) + ggtitle(glue::glue("FindAllMarker DEGs - Input data - dataset: {dataset}"))

      plot.name <- file.path(out.dir, glue("01-Heatmap_FindAllMarkers_Input_vs_Ref.pdf"))
      cli::cli_alert("Saving plot to {.file {plot.name}}")
      pdf(plot.name, width=12, height=10)
      plot(p)
      dev.off()
    } else {
      cli::cli_alert_warning("Not enough DEG genes (input vs other) to plot heatmap. skipping")
    }
    # assign back to dataset
    assign(dataset, sd)
    rm(sd, sd.downsampled)
  }


  #  01f: Barplots Input Healthy & RCASneg
  cli::cli_h2("01 - BarPlots - input_subset")
  seurat.reference <- Seurat::SetIdent(seurat.reference, value = "input_subset")
  seurat.reference <- .seuratFactorizeMdata(seurat.reference)

  file.name <- file.path(
    out.dir, paste0("01-BarPlots_rawProportions_inputSubset.pdf")
  )


  plist <- list()
  group.by.vec <- c(
    "sample_id", "sample_type", "disease_state","rcas_module_bin", "rcas_counts_bin"
  )
  for ( group.by in group.by.vec ){
    my.pal <-  .color_pal[[group.by]][
      levels(seurat.reference@meta.data[,group.by])]
    plist[[group.by]] <- SCpubr::do_BarPlot(
      seurat.reference,
      split.by = "input_subset",
      group.by = group.by,
      facet.by = "input_subset",
      colors.use = my.pal,
      add.n = F,
      position = 'fill'
    )
  }
  require(patchwork)
  # Combine the plots
  p_combined <- plist[[1]]
  for (i in 2:length(plist)) {
    p_combined <- p_combined + plist[[i]]
  }
  pdf(file.name, width = 18, height = 18)
  # Print the combined plot
  print(p_combined)
  dev.off()



  #   02a UMAP Input
  # UMAP DimPlots Input (active.ident and subset data using cell.barcodes)

  # mostly to check distributions etc
  cli::cli_h2("UMAP DimPlots - Input data - Idents and Subset")
  file.name <- file.path(out.dir, paste0("02-UMAPs-input.pdf"))
  cli::cli_alert("Plotting UMAP DimPlots for input_subset to {.file {file.name}}")

  # do_DimPlots
  p0 <- SCpubr::do_DimPlot(
    reduction = umap.reference,
    group.by = "input_subset",
    seurat.reference,
    label = TRUE,
    plot_cell_borders = FALSE
  ) + ggtitle(paste0("input_subset"))

  p1 <- SCpubr::do_DimPlot(
    seurat.full,
    reduction = umap.full,
    group.by = "input_subset",
    label = TRUE,
    plot_cell_borders = FALSE
  ) + ggtitle(paste0("input_subset"))

  p2 <- SCpubr::do_DimPlot(
    reduction = umap.subset,
    group.by = active.ident,
    colors.use = .color_pal[[active.ident]],
    seurat.reference,
    label = TRUE,
    plot_cell_borders = FALSE
  ) + ggtitle(paste0(active.ident))

  pdf(file.name, width = 24, height = 12)
  plot(p0 + p1 +  p2 )
  dev.off()




  #  subclustering: sd.subset
  # Process the subseted datasets
  #   Louvain SNN clusterings
  cli::cli_h2("Do Louvain SNN clusterings. Res: {.var {louvain.res}}")

  reduction.name <- umap.subset

  #  Create reults table
  clustersets.table <- data.frame(
    barcode = colnames(sd.subset),
    analysis = analysis.name
  )
  clusterset.names <- list()


  dataset <- c("sd.subset")
  cli::cli_h2(glue::glue("Processing {dataset} object"))
  sd <- get(dataset)
  sd <- .seuratFilterFoo(sd, filter.object = G1)

  # Core preprocessing
  sd <- Seurat::NormalizeData(sd)
  sd <- Seurat::FindVariableFeatures(sd, nfeatures = nfeatures)
  sd <- Seurat::ScaleData(sd, verbose = TRUE)
  sd <- Seurat::RunPCA(sd, verbose = TRUE, npcs = dims[2])
  sd <- Seurat::RunUMAP(sd, reduction.name = reduction.name, dims = dims[1]:dims[2], verbose = TRUE)

  # Graph + neighbors
  cli::cli_alert("Adding graphs to sd.subset")
  snn.graph.name <- paste0("SNN_vst_p", snn.prune.factor)
  sd@graphs <- list()
  sd <- Seurat::FindNeighbors(sd,
    dims = 1L:ncol(sd@reductions$pca@cell.embeddings),
    k.param = k.param,
    prune.SNN = 1 / snn.prune.factor,
    graph.name = snn.graph.name)

  # Clustering
  cli::cli_alert(glue::glue("Doing Louvain SNN clustering on {dataset}"))
  subclustering.name <- glue::glue("snn{gsub('[.]','', louvain.res)}_{dataset}")
  clusterset.names[[dataset]] <- subclustering.name

  sd <- Seurat::FindClusters(sd,
    graph.name = names(sd@graphs)[1],
    resolution = louvain.res,
    algorithm = 1)

  # Safer: read cluster IDs directly from Idents, then zero-pad
  cl_raw <- as.character(Seurat::Idents(sd))
  cluster.assignment.vec <- sprintf("%02d", as.integer(cl_raw))

  print(table(cluster.assignment.vec))
  clustersets.table[[subclustering.name]] <- paste0(cluster.assignment.vec, "_", dataset)

  # Write cluster sizes
  dims_line <- paste(as.integer(table(cluster.assignment.vec)), collapse = ",")
  writeLines(dims_line, con = file.path(out.dir, glue::glue("data-dims-{subclustering.name}.csv")))

  assign(dataset, sd)
  rm(sd)


  names(sd.subset@graphs)
  #  Do ElbowPlot
  cli::cli_alert("Doing ElbowPlot")
  p1 <- ElbowPlot(sd.subset, ndims = length(sd.subset@reductions$pca)) + ggtitle("sd.subset")
  # p2 <- ElbowPlot(sd_raw_sub, ndims = 30) + ggtitle("sd_raw_sub")
  p <- p1
  # plot to pdf
  pdf(file.path(out.dir, "02-ElbowPlot.pdf"), width=12, height=12)
  plot(p)
  dev.off()



  #  calculate number of clusters
  #  calculate number of identified clusters

  clustersets <- grep("snn", colnames(clustersets.table), value = TRUE)
  oneCluster <- FALSE
  for(clusterset in clustersets){
    if(length(unique(clustersets.table[[clusterset]])) == 1){
      oneCluster <- TRUE
    } else {
      cli::cli_alert("Number of clusters identified: {length(unique(clustersets.table[[clusterset]]))}")
    }
  }


  #  exit function if only one cluster was identified
  # exit function if only one cluster was identified.
  if (oneCluster){
    cli::cli_alert("Only one cluster was identified. Exiting function.")
    return(clustersets.table)
  }


  #  save the clustersets table as csv
  cli::cli_alert("Saving clustersets table assignments as csv")
  write.csv(clustersets.table, file.path(out.dir, "02-clustersets_assignments_table.csv"), row.names=F)




  #  03c Barplot of cluster assignments
  ## Make a plain barplot of cluster assignments

  cli::cli_alert("Making a barplot of cluster assignments")
  str(clustersets.table)
  # melt clustersets.table, id.vars = "barcode", grep the 'subclusters'. Set value name to 'subcluster'
  subclusters_n <- reshape2::melt(
    clustersets.table, value.name = "snn_cluster", variable.name = "clusterset",
    id.vars = "barcode",
    measure.vars = grep("snn", colnames(clustersets.table), value = TRUE))
  str(subclusters_n)
  p0 <- .barplot_stacked(
    plot_df = subclusters_n, group.var = "clusterset", color.var = "snn_cluster"
  )
  file.name <- file.path(out.dir, "03-Barplot_ClusterAssignments.pdf")
  pdf(file.name, width = 4, height = 4)
  print(p0)
  dev.off()

  rm(file.name, p0, subclusters_n)


  #  Add SubClustering assignments to metadata (for ditto heatmaps)

  cli::cli_alert("Adding subcluster assignments to metadata")
  # clustersets <- c("subcluster025_ds", "subcluster05_ds", "subcluster025_raw", "subcluster05_raw")
  clustersets <- grep("snn", colnames(clustersets.table), value = TRUE)
  datasets <- c("sd.subset", "seurat.reference", "seurat.full")

  for (clusterset in clustersets){
    for (dataset in datasets){
      sd <- get(dataset)
      sd <- AddMetaData(sd, metadata = clustersets.table %>%  dplyr::select(barcode,!!clusterset)  %>% tibble::column_to_rownames("barcode"))
      sd[[clusterset]][is.na(sd[[clusterset]])] <- "other"
      assign(dataset, sd)
    }
  }
  rm(clusterset, clustersets, dataset, datasets, mdata, sd)

  # For the two reference sets, seurat.reference and sd full,
  # Create and alternative annoation - active.ident - but where the subcusters are introduced
  cli::cli_alert("Creating alternative annotation - active.ident - with subclusters introduced")
  datasets <- c("seurat.reference", "seurat.full")
  clusterset <- "snn025_sd.subset"
  for (dataset in datasets){
    sd <- get(dataset)
    sd[["ident_w_subclusters"]] <- as.character(sd@meta.data[[active.ident]])
    v <- sd[[clusterset]]
    idx <- !is.na(v) & v != "other"
    sd[["ident_w_subclusters"]][idx] <- v[idx]
    table(sd[["ident_w_subclusters"]])

    assign(dataset, sd)
  }


  #   DEG Genes sd.subset
  # Loop all variants of defining clustersets
  # Loop both datasets
  # Find DEG genes between clustersets
  cli::cli_h1("Finding DEG genes between clustersets")
  cli::cli_h2("Loop all cluster assignments and datasets")
  # set Ident to the cluster column (graph_cluster)
  clustersets <- grep("snn", colnames(clustersets.table), value = TRUE)
  datasets <- c("sd.subset")

  allMarkers_list <- list()

  # loop clusters and datasets
  for( clusterset in clustersets ){
    cli::cli_alert(paste0("Looping clusterset assignment ", clusterset))

    for( dataset in datasets ){
      cli::cli_alert(paste0("Looping dataset ", dataset))
      # get data into temp dataset
      sd <- get(dataset)
      sd <- Seurat::SetIdent(sd, value = clusterset)
      table(Idents(sd))
      cli::cli_alert(paste0("Finding DEG genes for ", clusterset, " in ", dataset))
      annot_by <- c(annot_columns, grep("snn", colnames(sd@meta.data), value = TRUE))

      # Seurat FindAllMarkers to get DEG genes
      table(Idents(sd))
      allMarkers_list[[clusterset]][[dataset]] <- FindAllMarkers(sd, only.pos = TRUE)

      df <- allMarkers_list[[clusterset]][[dataset]] %>%
        arrange(cluster, p_val_adj)

      # Save DEG genes to file
      write.csv(df, file.path(out.dir, paste0("05-DEG_", clusterset, "_", dataset, ".csv")))

      #  Select topN genes from each Subcluster (From this dataset and clusterset definition)
      topGenes <- df %>%
        group_by(cluster) %>%
        arrange(p_val_adj, .by_group = TRUE) %>%
        slice_min(n=topN, p_val_adj ) %>%
        # ungroup() %>%
        pull(gene)

      # Heatmap of DEG genes - topN per clusterset
      for ( plot_n_genes in c(15)){
        plot.name <- file.path(out.dir, paste0("05-Heatmap_DEG_", clusterset, "_", dataset, "_top", plot_n_genes,".pdf"))
        # ditto Heatmap of DEG genes - top 10 per clusterset
        p <- suppressWarnings(dittoSeq::dittoHeatmap(
          genes = topGenes,
          drop_levels = T,
          sd,
          # filename = plot.name, width=12, height=24,
          cluster_cols = TRUE,
          # scaled.to.max=TRUE,
          complex = TRUE,
          annot.by = annot_by
        ))
        # plot to pdf
        pdf(plot.name, width=16, height=4*length(unique(clustersets.table[,clusterset])))
        plot(p)
        dev.off()
      }
    } # end dataset
  } # end clusterset

  rm(clusterset, dataset, datasets, topGenes, df, p, plot.name, sd)


  #  05b UMAP SubClusters
  #  UMAP variations for this dataset/clusterset assignment combo
  cli::cli_h1("Plotting UMAP SubClusters. 05b")
  clustersets <- grep("snn", colnames(clustersets.table), value = TRUE)
  datasets <- c("sd.subset")
  annot_by <- c(grep("snn", colnames(sd.subset@meta.data), value = TRUE))

  # loop clustersets and datasets
  fname <- file.path(out.dir, paste0("05-UMAP_SubClusters", ".pdf"))
  pdf(fname, width=24, height=24)

  for( clusterset in clustersets ){ # clusterset <- clustersets[1]
    cli::cli_alert(paste0("Looping clusterset assignment ", clusterset))

    for( dataset in datasets ){ # dataset <- datasets[1]
      cli::cli_alert(paste0("Looping dataset ", dataset))
      # get data into temp dataset
      sd <- get(dataset)
      sd <- Seurat::SetIdent(sd, value = clusterset)
      #table(Idents(sd))
      # Visualize the sub-cluster assignments on UMAP reductiions
      # Plot UMAPS generatad on subset data as well as reductions on full data
      cli::cli_alert(paste0("Plotting UMAPs for cluster assignment ", clusterset, " in ", dataset))

      p1 <- SCpubr::do_DimPlot(
        #reduction = "umap_qF03_dbl10_ds",
        sd,
        # colors.use = .colorPalDiscrete(unique(sort(as.character(Idents(sd))))),
        colors.use = .color_pal[["subcluster_sd_sub"]],
        label = TRUE,
        plot_cell_borders = FALSE
        #plot_marginal_distributions = TRUE
      ) + ggtitle(paste0(clusterset, " in ", dataset))
      # set background to black

      p2 <- SCpubr::do_DimPlot(
        seurat.reference,
        group.by = clusterset,
        colors.use = .color_pal[["subcluster_sd_sub"]],
        label = FALSE,
        plot_cell_borders = FALSE
        #plot_marginal_distributions = TRUE
      ) + ggtitle(paste0(clusterset, " in ", dataset))

      p3 <- SCpubr::do_DimPlot(
        group.by = "rcas_counts_bin",
        colors.use = .color_pal[["rcas_counts_bin"]],
        sd,
        label = FALSE,
        plot_cell_borders = FALSE
        #plot_marginal_distributions = TRUE
      ) + ggtitle("rcas_bin", subtitle = paste0(clusterset, " in ", dataset))

      #  p3 + theme(panel.background = element_rect(fill = "gray95"))
      p4 <- SCpubr::do_DimPlot(
        group.by = "sample_type",
        colors.use = .color_pal[["sample_type"]],
        # legend.position = "none",
        # cells.highlight = cell.barcodes,
        sd,
        label = FALSE,
        plot_cell_borders = FALSE
        #plot_marginal_distributions = TRUE
      ) + ggtitle("sample_type", subtitle = paste0(clusterset, " in ", dataset))

      p5 <- SCpubr::do_DimPlot(
        group.by = "rcas_module_bin",
        colors.use = .color_pal[["rcas_module_bin"]],
        # legend.position = "none",
        # cells.highlight = cell.barcodes,
        sd,
        label = FALSE,
        plot_cell_borders = FALSE
        #plot_marginal_distributions = TRUE
      ) + ggtitle("rcas_both_3bin", subtitle = paste0(clusterset, " in ", dataset))

      p6 <- SCpubr::do_DimPlot(
        group.by = "sample_id",
        colors.use = .color_pal[["sample_id"]],
        # legend.position = "none",
        # cells.highlight = cell.barcodes,
        sd,
        label = FALSE,
        plot_cell_borders = FALSE
        #plot_marginal_distributions = TRUE
      ) + ggtitle("sample_id", subtitle = paste0(clusterset, " in ", dataset))



      plot(p2 + p1 +  p6 + p4 + p3 + p5 )
    } # end loop over datasets
  } # end clustersets
  dev.off()
  rm(clusterset, dataset, datasets, p1, p2, sd)



  #  06 Beeswarms of Signatures
  #  GENE SIGNATURE SCORES OVER ALL DATASETS n clustersetS

  clustersets <- grep("snn", colnames(clustersets.table), value = TRUE)


  # loop clustersets and datasets
  do.beeswarms = FALSE
  if ( do.beeswarms ){

    for ( dataset in datasets){ # dataset <- datasets[1]
      for ( clusterset in clustersets ){ # clusterset <- clustersets[1]
        cli::cli_h3(paste0("Looping clusterset assignment: ", clusterset))
        cli::cli_h3(paste0("Looping dataset: ", dataset))

        file.name <- file.path(out.dir, paste0("06_ClusterModuleScores_Beeswarm_", clusterset, "_", dataset, ".pdf"))
        cli::cli_alert("Plotting to {.file {file.name}}")
        pdf(file.name, width=27, height=18)

        # a tmp object to caclulate cell signature scores over all data

        # calculate genscores over all data - blot beewarm

        # read DEG genes from csv for this clusterset
        fname.deg <- file.path(out.dir, paste0("DEG_", clusterset, "_", dataset, ".csv"))
        cli::cli_alert(paste0("Reading DEG genes from {.file {fname.deg}}"))
        df <- read_csv(fname.deg)
        subclusters <- sort(unique( df$cluster ))
        subclusters.deg.list <- list()
        # loop all subclusters
        for ( subcluster in subclusters ){ # subcluster <- subclusters[1]
          cli::cli_alert(paste0("Calculating gene scores for clustersets: ", clusterset, " in ", dataset, " subclusters"))

          # table(sd_raw.tmp[["input_subset"]])
          # get the topN genes for the subcluster
          topGenes <- df %>%
            dplyr::filter(cluster == subcluster) %>%
            arrange(p_val_adj) %>%
            slice_min(n=topN, p_val_adj ) %>%
            # ungroup() %>%
            pull(gene) %>%
            # remove genes not present in rownames(sd_raw)
            intersect(rownames(sd_raw))

          subclusters.deg.list[[dataset]][[clusterset]][[subcluster]] <- topGenes

          # if number (length) of topGenes is less than 5, skip and mext subcluster
          if ( length(topGenes) < 5 ){
            cli::cli_alert(paste0("Skipping subcluster ", subcluster, " in ", clusterset, " in ", dataset, " as less than 5 genes"))
            next
          }
          # print(sort(topGenes))
          # calculate the gene scores (full data)
          cli::cli_alert(paste0("Calculating gene scores for ", subcluster, " in ", clusterset, " in sd_raw.tmp"))

          sd_raw.tmp <- Seurat::AddModuleScore(
            sd_raw,
            list(topGenes),
            name = subcluster)

          # caclulate gene scores (subset data sd.subset)
          cli::cli_alert(paste0("Calculating gene scores for ", subcluster, " in ", clusterset, " in sd.subset"))

          sd.subset.tmp <- Seurat::AddModuleScore(
            sd.subset,
            list(topGenes),
            name = subcluster)

          # caclulate gene scores (subset data sd.subset)
          cli::cli_alert(paste0("Calculating gene scores for ", subcluster, " in ", clusterset, " in sd_raw_sub"))

          sd_raw_sub.tmp <- Seurat::AddModuleScore(
            sd_raw_sub,
            list(topGenes),
            name = subcluster)

          # plot the beeswarm plot
          cli::cli_alert("Plotting beeswarm plot for clusterset {.var {clusterset}} and subcluster: {.var {subcluster}}")
          module.name <- paste0(subcluster,"1")

          p1 <- SCpubr::do_BeeSwarmPlot(
            group.by = active.ident,
            pt.size = 0.1,
            sd_raw.tmp,
            feature_to_rank = module.name,
            order=TRUE,
            continuous_feature = TRUE,
            raster = TRUE,
            verbose = FALSE,
            plot_cell_borders = FALSE
          )
          p2 <- SCpubr::do_BoxPlot(
            group.by = active.ident,
            sd_raw.tmp,
            feature = module.name,
            use_silhouette = TRUE
          )
          p3 <- SCpubr::do_BoxPlot(
            group.by = clusterset,
            sd_raw.tmp,
            feature = module.name,
            use_silhouette = TRUE
          )
          p4 <- SCpubr::do_FeaturePlot(reduction = umap.subset,
            sd_raw.tmp,
            features = module.name,
            pt.size = 2,
            raster = TRUE,
            verbose = FALSE,
            plot_cell_borders = TRUE
          )
          p5 <- SCpubr::do_FeaturePlot(
            #reduction = "umap_qF03_dbl10_ds",
            sd.subset.tmp,
            features = module.name,
            #colors.use = .color_pal[[clusterset]],
            pt.size = 8,
            raster = TRUE,
            verbose = FALSE,
            plot_cell_borders = TRUE
            #plot_marginal_distributions = TRUE
          ) + ggtitle(paste0(subcluster, " in sd.subset.tmp"))
          # p6 <- SCpubr::do_FeaturePlot(
          #   sd_raw_sub.tmp,
          #   features = module.name,
          #   pt.size = 8,
          #   raster = TRUE,
          #   verbose = FALSE,
          #   plot_cell_borders = TRUE
          # )+ ggtitle(paste0(subcluster, " in sd_raw_sub.tmp"))


          plot(p1 + p2 + p3 + p4 + p5 )
          rm(sd_raw.tmp, p1, p2, p3, p4, p5, topGenes, module.name, sd.subset.tmp, sd_raw_sub.tmp)
        } # end loop module scores over subclusters
        dev.off()
        rm(df, subclusters)
      } # end datasets
    } # end clusters


    #  01d beeswarm plot for input vs rest
    # Also add a beeswarm for input vs rest
    # datasets   <- c("seurat.reference")
    cli::cli_h1("Plotting beeswarm plot for input vs rest")

    datasets   <- c("sd_raw")
    subcluster <- c("input_subset")
    clustersets <- grep("snn", colnames(clustersets.table), value = TRUE)
    # loop datasets
    for ( dataset in datasets){ # dataset <- datasets[1]
      for ( clusterset in clustersets ){ # clusterset <- clustersets[1]

        cli::cli_h3(paste0("Looping dataset: ", dataset, " ans clusterset ",clusterset))

        # Get DEG genws from genelists_findallmarkers
        topGenes <- genelists_findallmarkers[[dataset]] %>%
          dplyr::arrange( desc(avg_log2FC) ) %>%
          dplyr::filter(cluster != "other") %>%
          dplyr::filter(avg_log2FC > 2) %>%
          slice_head(n = topN) %>%
          pull(gene)
        #ungroup()
        if ( length(topGenes) < 5 ){
          cli::cli_alert(paste0("Not enough genes - Skipping subcluster ", subcluster, " in ", clusterset, " in ", dataset, " as less than 5 genes"))
          next
        }

        # calculate genscores over all data - blot beewarm
        cli::cli_alert(paste0("Calculating gene scores for input_samples ", " in ",clusterset,"_", dataset, " subclusters"))

        # print(sort(topGenes))
        # calculate the gene scores (full data)
        cli::cli_alert(paste0("Calculating gene scores for input_samples",clusterset,"_"," in sd_raw.tmp"))

        sd_raw.tmp <- Seurat::AddModuleScore(
          seurat.full,
          list(topGenes),
          name = subcluster)

        # caclulate gene scores (subset data sd.subset)
        cli::cli_alert(paste0("Calculating gene scores for input_samples",clusterset,"_", " in sd.subset"))

        sd.subset.tmp <- Seurat::AddModuleScore(
          sd.subset,
          list(topGenes),
          name = subcluster)

        # caclulate gene scores (subset data sd.subset)
        cli::cli_alert(paste0("Calculating gene scores for input_samples  ",clusterset,"_", " in sd_raw_sub"))

        sd_raw_sub.tmp <- Seurat::AddModuleScore(
          sd_raw_sub,
          list(topGenes),
          name = subcluster)

        file.name <- file.path(out.dir, paste0("01d_ClusterModuleScore_Beeswarm_input_samples_",clusterset,"_", dataset, ".pdf"))
        cli::cli_alert("Plotting to {.file {file.name}}")
        pdf(file.name, width=18, height=18)

        # plot the beeswarm plot
        cli::cli_alert("Plotting beeswarm plot for input_cluster")
        module.name <- paste0(subcluster,"1")
        p1 <- SCpubr::do_BeeSwarmPlot(
          group.by = active.ident,
          pt.size = 0.1,
          sd_raw.tmp,
          feature_to_rank = module.name,
          order=TRUE,
          continuous_feature = TRUE,
          raster = TRUE,
          verbose = FALSE,
          plot_cell_borders = FALSE
        )
        p2 <- SCpubr::do_BoxPlot(
          group.by = active.ident,
          sd_raw.tmp,
          feature = module.name,
          use_silhouette = TRUE
        )
        p3 <- SCpubr::do_BoxPlot(
          group.by = clusterset,
          sd_raw.tmp,
          feature = module.name,
          use_silhouette = TRUE
        )
        p4 <- SCpubr::do_FeaturePlot(reduction = umap.subset,
          sd_raw.tmp,
          features = module.name,
          pt.size = 2,
          raster = TRUE,
          verbose = FALSE,
          plot_cell_borders = TRUE
        )

        plot(p1 + p2 + p3 + p4 )
        rm(sd_raw.tmp, p1, p2, p3, p4, topGenes, module.name, sd.subset.tmp, sd_raw_sub.tmp)

        dev.off()
      } # end loop clusterset
    } # end datasets

  } # end if do.beeswarms


  #





  #  UMAPS

  ## UMAPS
  # nFeatures
  # nCount_RNA
  # cell cycle
  cli::cli_h2("UMAPs")

  # mostly to check distributions etc
  file.name <- file.path(out.dir, paste0("08a-UMAPs_SubcClusters.pdf"))
  pdf(file.name, width = 18, height = 9)

  # loop clusters and datasets
  clustersets <- grep("snn", colnames(clustersets.table), value = TRUE)
  # datasets <- c("sd.subset", "sd_raw_sub")

  for( clusterset in clustersets ){
    cli::cli_alert("Plotting UMAPs for clusterset {.var {clusterset}}")
    # for( dataset in datasets ){

    p1 <- SCpubr::do_DimPlot(
      group.by = clusterset,
      colors.use =
        .palette_discrete_bp(sort(unique(seurat.reference[[clusterset]][,1]))),
      seurat.reference,
      label = TRUE,
      plot_cell_borders = FALSE
      #plot_marginal_distributions = TRUE
    ) + ggtitle(paste0(clusterset))
    p2 <- SCpubr::do_DimPlot(
      group.by = clusterset,
      colors.use =
        .palette_discrete_bp(sort(unique(seurat.reference[[clusterset]][,1]))),
      seurat.reference,
      label = FALSE,
      plot_cell_borders = FALSE
    ) + ggtitle(paste0(clusterset))
    plot( p1 + p2 )
    # dev.off()
  } # end loop over clusters

  dev.off()



  #  BARPLOTS (subclusters)

  ## UMAPS
  do.barplots.8b <- F
  if ( do.barplots.8b ){

    # Now plot normalized values
    file.name <- file.path(
      out.dir, paste0("08b-BarPlots_RawProportions_subClusters.pdf")
    )
    datasets <- c("seurat.reference")
    pdf(file.name, width = 14, height = 14)
    for( clusterset in clustersets ){
      for ( dataset in datasets ){
        cli::cli_alert("Plotting normalized barplots for dataset {.var {dataset}} and clusterset {.var {clusterset}}")
        sd <- get( dataset )
        mdata <- sd@meta.data
        table(mdata$sample_type, mdata[[clusterset]])
        df <- mdata %>%
          dplyr::filter(input_subset != "other")
        p0 <- .barplot_stacked(
          my.pal = .color_pal[["sample_type"]],
          plot_df = df,
          group.var = clusterset,
          color.var = "sample_type",
          scaled.y = F,
          scaled.y.wihtin.color = F
        )

        p1 <-  .barplot_stacked(
          my.pal = .color_pal[["sample_type"]],
          plot_df = mdata,
          group.var = clusterset,
          color.var = "sample_type",
          scaled.y = T,
          scaled.y.wihtin.color = F
        )

        p2 <-  .barplot_stacked(
          my.pal = .color_pal[["rcas_counts_bin"]],
          plot_df = mdata,
          group.var = clusterset,
          color.var = "rcas_counts_bin",
          scaled.y = T,
          scaled.y.wihtin.color = F
          #scaled.y = T
        )
        p3 <-  .barplot_stacked(
          my.pal = .color_pal[["sample_id"]],
          plot_df = mdata,
          group.var = clusterset,
          color.var = "sample_id",
          scaled.y = T,
          scaled.y.wihtin.color = F
          #scaled.y = T
        )

        # plot color.var = "rcas_both"
        pX <-  .barplot_stacked(
          plot_df = mdata,
          my.pal = .color_pal[["rcas_module_bin"]],
          group.var = clusterset,
          color.var = "rcas_module_bin",
          scaled.y = T,
          scaled.y.wihtin.color = F
        )

        #
        #         # Add subseteted plots -  rcas (tumor only)
        #
        #         df <- mdata %>%
        #           dplyr::filter(sample_type %in% c("Primary","Recurrent"))
        #
        #         p5 <-  .barplot_stacked(
        #           plot_df = df,
        #           my.pal = .color_pal[["sample_type"]],
        #           group.var = clusterset,
        #           color.var = "sample_type",
        #           scaled.y = T,
        #           scaled.y.wihtin.color = F
        #         ) + ggtitle(glue::glue("Primary and Recurrent Samples only
        #               (Non-healthy) n={nrow(df)}"))
        #
        #
        plot( p0 + p1 + p2 + p3 + pX )
      }
    }
    dev.off()

    # plot normalized values

  }


  #  FGSEA PER clusterset
  # See if any of the recursive subclusters identify as a known cell type
  # subclusters are compared to other cells (As a group) using FindAllMarkers
  cli::cli_h2("Starting FGSEA")


  celltype_list <-readRDS(file = "./data/metadata/genesets/fgsea-celltypes.Rds")
  #  grep("microglia", names(celltype_list), ignore.case = T, value = T)

  fgsea.dir <- file.path(out.dir, "fgsea_results")
  dir.create(fgsea.dir, showWarnings = FALSE, recursive = TRUE)

  #fname <- paste0("09a-FGSEA_dotplots_GeneSets.1.2_all.pdf")
  # file.name <- file.path(out.dir, fname)
  #pdf(file.name, width = 10, height = 7)

  cli::cli_alert("Plotting to file: {.file {file.name}}")
  # loop clusters and datasets
  clustersets <- grep("snn", colnames(clustersets.table), value = TRUE)
  # datasets <- c("sd.subset", "sd_raw_sub")
  dataset <- c("seurat.reference")
  clusterset <- "snn025_sd.subset"

  #for ( clusterset in clustersets ){
  cli::cli_alert("Starting clusterset: {.info {clusterset}}")
  #  for (dataset in datasets){
  cli::cli_alert("Starting dataset: {.info {dataset}}")
  sd <- get(dataset)
  table(sd[[clusterset]])
  table(sd[["ident_w_subclusters"]])

  .wrapper_fgsea(
    seurat.object =  sd, # pre-filtered data
    # cluster.label =  clusterset, # vector
    cluster.label =  "ident_w_subclusters", # vector
    gene.sets.list = celltype_list,
    results.path = fgsea.dir,
    analysis.name = glue("brain-signatures")
  )

  fgsea.plot.file <- file.path(out.dir,  "09a-fgsea-DotPlots.pdf")
  pdf(fgsea.plot.file, width = 20, height = 20)
  .wrapper_fgsea_plot(
    fgsea.dir = fgsea.dir,  annotation = "ident_w_subclusters", seurat.object = sd)
  dev.off()



  #  MCA PER clusterset
  .wrapper_mca(
    out.dir,
    sd = sd.subset
    #sd = seurat.reference
  )
  .wrapper_mca_plot(
    out.dir,
    clustersets =  c("snn025_sd.subset"),
    #clustersets =  c("snn025_sd.subset"-,"ident_w_subclusters"),
    sd = sd.subset
  )


  cli::cli_h1(" =========================== ")
  cli::cli_h1("All Done! subClusteringWrap: {.var {subset.name}}")
  cli::cli_h1(" =========================== ")
  # sink()
} # End function

