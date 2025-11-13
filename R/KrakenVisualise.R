#' Kraken Report Visualise Distributions
#'
#' @inheritParams kraken_report_add_zscore
#' @param taxids_of_interest taxonomy ids to show as distinct facets (numeric)
#' @param y_aesthetic variable representing y axis. Typically "ScientificName".
#' @param metric Name of metric to use: c("Zscore", "ZscoreRobust", "MedianIQR", "MedianQn", "MedianTau2"). Note depending on how reports have been parsed, you may only have access to ZscoreRobust - which is the recommended choice (string)
#' @param show_y_axis_labels show y axis labels? (boolean)
#' @param show_sample_vlines show vertical lines indicating positions of samples of interest (boolean)
#' @param supporting_reads_threshold minimum number of reads that must be assigned to a taxid for us to be confident its unlikely to be noise (int)
#' @param sample_of_interest which samples are we most interested in - use show_sample_vlines flag to mark these samples (bool)
#' @param use_loggable_version should we use the loggable version of the metric.
#' @param seed random seed.
#' @param pointsize how large should the points be
#' @param ... other aesthetics passed to aes_string. Mostly used to show more properties in hoverbox once plots are made interactive
#'
#' @return ggplot
#' @export
#'
kraken_report_visualise_distributions_old <- function(kraken_report_df, taxids_of_interest, y_aesthetic = "ScientificName", metric = "ZscoreRobust", show_y_axis_labels = FALSE,show_sample_vlines = FALSE, supporting_reads_threshold = 50, sample_of_interest = NA, use_loggable_version = TRUE, seed = 1, pointsize = c(1, 2), ...){
  is_db = "tbl_sql" %in% class(kraken_report_df)

  # Defensive Assertions
  standardised_metrics = c("Zscore", "ZscoreRobust", "MedianIQR", "MedianQn", "MedianTau2")
  assertthat::assert_that(assertthat::is.string(metric))
  assertthat::assert_that(metric %in% standardised_metrics, msg = paste0("metric (", metric,") must be one of: [", paste0(standardised_metrics, collapse = ", "), "]"))
  #assertthat::assert_that(is.data.frame(kraken_report_df))
  assertthat::assert_that(y_aesthetic %in% colnames(kraken_report_df)) #Change this once we add the option to y_aesthetic by metadata

  if(use_loggable_version) metric <- paste0(metric, "Loggable")
  assertthat::assert_that(metric %in% colnames(kraken_report_df), msg = paste0("kraken_report_df must contain the column: [", metric, "]. Please see ?kraken_info_standardised_metrics for information on how to add standard metrics to your dataframe."))
  columns = c("SampleID", "TaxonomyID", "ScientificName", "ReadsCoveredByClade", "RPM", metric, "ZscoreRobust")
  assertthat::assert_that(all(columns %in% colnames(kraken_report_df)), msg = paste0("kraken_report_df must contain the columns [", paste0(columns, collapse=", "), "]"))
  assertthat::assert_that(length(pointsize) == 2, msg = "pointsize must a vector of length 2. First number is size of most samples and second number for size of 'samples of interest' ")
  #assertthat::assert_that(color_aesthetic %in% colnames(kraken_report_df), msg = paste0("kraken_report_df must contain the column [", paste0(color_aesthetic, collapse=", "), "], to use it as a color aesthetic"))

  # Warning messages
  if(y_aesthetic != "ScientificName"){ message("You've chosen to set a Y aesthetic - you will probably want to also set show_y_axis_labels = TRUE ")}

  # Subset the data
  subset_df = kraken_report_df %>%
    dplyr::filter(TaxonomyID %in% taxids_of_interest) #%>%
    #dplyr::select(SampleID, ScientificName, ZscoreRobust,ZscoreRobustLoggable, TaxonomyID, ReadsCoveredByClade, RPM)

  if(is_db)
    subset_df = dplyr::collect(subset_df)

  subset_df <- subset_df %>%
    dplyr::mutate(
      Infinite = is.infinite(subset_df[[metric]]),
      SampleOfInterest = SampleID %in% sample_of_interest,
      Shape = data.table::fcase(SampleOfInterest, "Sample of Interest", is.infinite(get(metric)), "Infinite", default = "General Sample")
      )

  # Relevel factors so order of taxids_of_interest means something
  scientific_names = subset_df$ScientificName[match(taxids_of_interest, subset_df$TaxonomyID)]
  #browser()
  subset_df <- subset_df %>%
    dplyr::mutate(ScientificName = forcats::fct_relevel(ScientificName, scientific_names))


  set.seed(seed)

  # d = subset_df[SampleID == sample_of_interest,]
  pos_jitter = ggplot2::position_jitter(width = 0, height = 0.5, seed = seed)

  #browser()
  plot = subset_df %>%
    ggplot2::ggplot(ggplot2::aes_string(x = metric, y= y_aesthetic, SampleID = "SampleID", ReadsCoveredByClade = "ReadsCoveredByClade", RPM = "RPM", Metric = metric, ZscoreRobust = "ZscoreRobust", customdata = "SampleID", ...), alpha = 0.5) +
    ggplot2::geom_jitter(
      data=~dplyr::filter(.x, !SampleOfInterest),
      position = pos_jitter,
      ggplot2::aes(
        shape = Shape,
        color =  ReadsCoveredByClade > supporting_reads_threshold,
        size = SampleOfInterest)
      ) +
    ggplot2::scale_x_continuous(oob = scales::oob_squish_infinite, trans="log10") +
    ggplot2::facet_wrap(~ScientificName, ncol = 1, scales = "free_y") +
    #ggplot2::theme(strip.text.x = ggplot2::element_blank()) +
    ggplot2::scale_shape_manual(values = c("General Sample" = 21, "Infinite" = 24, "Sample of Interest" = 17)) +
    #ggplot2::scale_alpha_manual(values = c(FALSE=0.7, TRUE=1)) +
    ggplot2::scale_size_manual(values = pointsize) +
    utilitybeltgg::theme_fivethirtyeight_two() +
    utilitybeltgg::theme_legend_right() +
    ggplot2::scale_color_manual(values = c("FALSE" = "red", "TRUE" = "blue")) +
    ggplot2::guides(color = ggplot2::guide_legend(paste0("ReadsCoveredByClade > ", supporting_reads_threshold)))


  if(sum(subset_df[["SampleOfInterest"]]) > 0  ){
    message("Drawing Samples Of Interest ")
    plot <- plot + ggplot2::geom_jitter(
        data=~dplyr::filter(.x, SampleOfInterest),
        position = pos_jitter,
        ggplot2::aes(
          shape = Shape,
          color =  ReadsCoveredByClade > supporting_reads_threshold,
          size = SampleOfInterest)
        ) +
      ggplot2::geom_vline(
        data=~dplyr::filter(.x, SampleOfInterest),
        color = "green",
        size = 2,
        ggplot2::aes_string(
          linetype = "SampleOfInterest",
          xintercept = metric)
        )
  }

  if(!show_y_axis_labels){
    plot = plot + ggplot2::theme(axis.text.y = ggplot2::element_blank(), axis.ticks.y = ggplot2::element_blank(), axis.title.y  = ggplot2::element_blank())
  }

  if(show_sample_vlines){
    plot = plot + ggplot2::geom_vline(data = function(df) {df[SampleID %in% sample_of_interest,]} ,ggplot2::aes_string(xintercept = metric, linetype = "SampleID"))
  }

  return(plot)
}


kraken_report_visualise_signal_focus <- function(kraken_report_df, parent_taxids){

  lapply(
    parent_taxids, kraken_calculate_proportion_of_signal_explained_by_n_strongest_taxids(kraken_report_df = as.data.frame(kraken_report_df), n = 1, parent_taxid)
    )

  return(results)
}

#' Confidence Palette
#'
#' @return colours for confidence levels (character)
#' @export
#'
#' @examples
#' kraken_confidence_palette()
kraken_confidence_palette  <- function(){
  c(
  'No Hit' = "#B3B3B3",
  'Very Low' = '#E5C494',
  'Low' =  '#FFD92F' ,
  'Medium' = "#8DA0CB",
  'High' = "#E78AC3"
  )
}

#' Kraken Report Visualise Distributions
#'
#' @inheritParams kraken_report_add_zscore
#' @param taxids_of_interest taxonomy ids to show as distinct facets (numeric)
#' @param show_y_axis_labels show y axis labels? (boolean)
#' @param show_sample_vlines show vertical lines indicating positions of samples of interest (boolean)
#' @param sample_of_interest which samples are we most interested in - use show_sample_vlines flag to mark these samples (bool)
#' @param ... other aesthetics passed to aes_string. Mostly used to show more properties in hoverbox once plots are made interactive
#'
#' @return ggplot
#' @export
#'
kraken_report_visualise_distributions <- function(kraken_report_df, taxids_of_interest, show_y_axis_labels = TRUE,show_sample_vlines = FALSE, sample_of_interest = NA, pointsize = 2, pointstroke = 1.5){
  is_db = "tbl_sql" %in% class(kraken_report_df)

  columns = c("SampleID", "TaxonomyID", "ScientificName", "ReadsCoveredByClade", "RPM", "ZscoreRobustLoggable", "Confidence")
  assertthat::assert_that(all(columns %in% colnames(kraken_report_df)), msg = paste0("kraken_report_df must contain the columns [", paste0(columns, collapse=", "), "]"))


  # Subset the data
  subset_df = kraken_report_df %>%
    dplyr::filter(TaxonomyID %in% taxids_of_interest)

  if(is_db)
    subset_df = dplyr::collect(subset_df)


  if(nrow(subset_df) == 0 )
    stop('None of the taxids of interest are in the kraken data')

  subset_df <- subset_df %>%
    dplyr::mutate(
      SampleOfInterest = SampleID %in% sample_of_interest
    )

  # Relevel factors so order of taxids_of_interest means something
  scientific_names = subset_df$ScientificName[match(taxids_of_interest, subset_df$TaxonomyID)]

  subset_df <- subset_df %>%
    dplyr::mutate(ScientificName = forcats::fct_relevel(ScientificName, scientific_names))


  plot = subset_df %>%
    ggplot2::ggplot(ggplot2::aes(x = ZscoreRobustLoggable, y = as.double(ReadsCoveredByClade), colour = Confidence, tooltip = SampleID)) +
    ggiraph::geom_point_interactive(shape = 21, size = pointsize, stroke = pointstroke) +
    ggplot2::facet_wrap(~ScientificName, ncol = 1, scales = "free_y") +
    ggplot2::scale_x_continuous(trans = "log10", oob = scales::oob_squish_infinite) +
    ggplot2::scale_y_continuous(trans = "log10", oob = scales::oob_squish_infinite, expand = ggplot2::expansion(c(0.1, 0.1))) +
    ggplot2::theme_bw() +
    ggplot2::ylab("Supporting Reads") +
    ggplot2::xlab("Microbial Signal\n(robust z-score, loggable)") +
    ggplot2::scale_color_manual(values = kraken_confidence_palette()) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(size = 20),
      axis.title.x = ggplot2::element_text(size = 22, face = "bold"),
      axis.title.y = ggplot2::element_text(size = 22, face = "bold")
    )

  if(!show_y_axis_labels){
    plot = plot + ggplot2::theme(axis.text.y = ggplot2::element_blank(), axis.ticks.y = ggplot2::element_blank(), axis.title.y  = ggplot2::element_blank())
  }

  if(show_sample_vlines){
    plot = plot + ggplot2::geom_vline(data = function(df) {df[SampleID %in% sample_of_interest,]} ,ggplot2::aes_string(xintercept = ZscoreRobustLoggable, linetype = "SampleID"))
  }

  return(plot)
}


kraken_report_visualise_signal_focus <- function(kraken_report_df, parent_taxids){

  lapply(
    parent_taxids, kraken_calculate_proportion_of_signal_explained_by_n_strongest_taxids(kraken_report_df = as.data.frame(kraken_report_df), n = 1, parent_taxid)
  )

  return(results)
}


#' kraken_report_visualise_single_sample
#'
#' @param kraken_report_df the kraken dataframe produced by \link{kraken_reports_parse} or link{kraken_report_parse}
#' @param samples_of_interest Sample IDs whose kraken report you want to visualise (character)
#' @param taxonomy_id_to_highlight taxonomy id to highlight (int)
#' @param supporting_reads_threshold taxons with greater than this number of reads will be indicated using the color aesthetic (int)
#' @param type_of_taxids_to_show one of c("all", "bacterial", "viral", "fungal", "sar")
#' @inheritParams kraken_calculate_proportion_of_signal_explained_by_n_strongest_taxids
#'
#' @return ggplot object
#' @export
#'
kraken_report_visualise_single_sample <- function(kraken_report_df, samples_of_interest, rank = "S",taxonomy_id_to_highlight = NULL, supporting_reads_threshold = 50, type_of_taxids_to_show = c("all", "bacterial", "viral", "fungal", "sar")){
  assertthat::assert_that("ZscoreRobustLoggable" %in% colnames(kraken_report_df), msg = "Report dataframe must include the column [ZscoreRobustLoggable]")
  type_of_taxids_to_show = rlang::arg_match(type_of_taxids_to_show)

  is_db = "tbl_sql" %in% class(kraken_report_df)

  taxonomy_name_to_highlight <- kraken_report_df[["ScientificName"]][match(taxonomy_id_to_highlight, kraken_report_df[["TaxonomyID"]])]

  data = kraken_report_df %>%
    dplyr::filter(SampleID %in% samples_of_interest & ReadsCoveredByClade > 0 & Rank == rank)

  if(is_db)
    data <- dplyr::collect(data)

  if(type_of_taxids_to_show == "bacterial"){
    message("Including only ", type_of_taxids_to_show, " taxids in plot")
    data = data %>%
      dplyr::filter(TaxonomyID %in% ncbitaxids::taxids_load_list("all_bacterial_taxids.05_05_2022"))
  }
  else if(type_of_taxids_to_show == "viral"){
    message("Including only ", type_of_taxids_to_show, " taxids in plot")
    data = data %>%
      dplyr::filter(TaxonomyID %in% ncbitaxids::taxids_load_list("all_viral_taxids.05_05_2022"))
  }
  else if(type_of_taxids_to_show == "fungal"){
    message("Including only ", type_of_taxids_to_show, " taxids in plot")
    data = data %>%
      dplyr::filter(TaxonomyID %in% ncbitaxids::taxids_load_list("all_fungal_taxids.05_05_2022"))
  }
  else if(type_of_taxids_to_show == "sar"){
    message("Including only ", type_of_taxids_to_show, " taxids in plot")
    data = data %>%
      dplyr::filter(TaxonomyID %in% ncbitaxids::taxids_load_list("all_sar_taxids.05_05_2022"))
  }
  else if(type_of_taxids_to_show == "sar"){
    message("Including only ", type_of_taxids_to_show, " taxids in plot")
    data = data %>%
      dplyr::filter(TaxonomyID %in% ncbitaxids::taxids_load_list("all_sar_taxids.05_05_2022"))
  }
  else if(type_of_taxids_to_show == "all"){
   message("Including all taxids in plot")
  }

  #browser()
  #browser()
  gg = data %>%
    ggplot2::ggplot(
      ggplot2::aes(
        color = ReadsCoveredByClade > supporting_reads_threshold,
        ReadsCoveredByClade = ReadsCoveredByClade,
        x=ZscoreRobustLoggable,
        ZscoreRobust = ZscoreRobust,
        y=RPM,
        ScientificName = ScientificName,
        TaxonomyID = TaxonomyID,
        #ReadsCoveredByClade = ReadsCoveredByClade,
        shape = is.infinite(ZscoreRobust)
        #color = TaxonomyID %in% taxonomy_id_to_highlight,
        )) +
    ggplot2::scale_shape_manual(values = c("FALSE" = 21,"TRUE" = 24)) +
    ggplot2::scale_x_continuous(trans = "log10", oob = scales::oob_squish_infinite) +
    ggplot2::scale_y_continuous(trans = "log10", oob = scales::oob_squish_infinite) +
    ggplot2::ggtitle(paste0(samples_of_interest, collapse = ",")) +
    ggplot2::guides(shape=ggplot2::guide_legend("Infinite")) +
    ggplot2::geom_point() +
    #ggthemes::theme_fivethirtyeight() +
    #ggthemes::theme_calc() +
    #utilitybeltgg::theme_common_adjustments() +
    utilitybeltgg::theme_fivethirtyeight_two() +
    ggplot2::scale_color_brewer(palette = "Set2") +
    ggplot2::guides(color = ggplot2::guide_legend(paste0("ReadsCoveredByClade  >", supporting_reads_threshold)))

  if (is.null(taxonomy_id_to_highlight)){
    #gg = gg + ggplot2::guides(color="none")
  }
  else
  {
    gg = gg + ggrepel::geom_text_repel(data = function(df) {df[(TaxonomyID %in% taxonomy_id_to_highlight),]}, ggplot2::aes(label=ScientificName), color="black", min.segment.length = 0, nudge_y=3)
    #gg = gg + ggplot2::guides(color = ggplot2::guide_legend(taxonomy_name_to_highlight))
  }

  return(gg)

}


#' Interactive Single Sample Visualisation
#'
#' @param tooltip aesthetics to include in tooltip c("ScientificName", "TaxonomyID", "ReadsCoveredByClade", "ZscoreRobust", "RPM")
#' @param ... arguments passed to \strong{kraken_report_visualise_single_sample}
#' @return plotly graph
#' @export
#'
kraken_report_visualise_single_sample_interactive <- function(tooltip = c("ScientificName", "TaxonomyID", "ReadsCoveredByClade", "ZscoreRobust", "RPM"), ...){
  plotly::ggplotly(kraken_report_visualise_single_sample(...), tooltip = tooltip)
}

kraken_report_visualise_two_group_comparision <- function(kraken_report_df, taxonomy1 = 2, taxonomy2 = 10239, sample_of_interest = NA){

  subset = kraken_report_df[(TaxonomyID %in% c(taxonomy1, taxonomy2)),]
  taxname1 = subset[["ScientificName"]][match(taxonomy1, subset[["TaxonomyID"]])]
  taxname2 = subset[["ScientificName"]][match(taxonomy2, subset[["TaxonomyID"]])]

  gg = subset %>%
    tidyr::pivot_wider(id_cols = SampleID, names_from = ScientificName, values_from = ZscoreRobustLoggable) %>%
    ggplot2::ggplot(ggplot2::aes_string(x=taxname1, y=taxname2, SampleID = "SampleID")) +
    ggplot2::geom_point() +
    ggplot2::scale_x_continuous(trans = "log10", oob = scales::oob_squish_infinite)+
    ggplot2::scale_y_continuous(trans = "log10", oob = scales::oob_squish_infinite) +
    utilitybeltgg::theme_fivethirtyeight_two()

  if(!is.na(sample_of_interest)){
    gg = gg + ggplot2::geom_point(data = function(df) {dplyr::filter(df, SampleID %in% sample_of_interest)} ,ggplot2::aes(color = "Sample of Interest"))
  }

  return(gg)
}


kraken_report_umap <- function(kraken_report_df){

}

kraken_report_tsne <- function(kraken_report_df){
  #assertthat::assert_that("ZscoreRobust" %in% colnames(kraken_report_df), msg = "Report dataframe must include the column [ZscoreRobustLoggable]")
  assertthat::assert_that("ZscoreRobust" %in% colnames(kraken_report_df), msg = "Report dataframe must include the column [ZscoreRobustLoggable]")
  assertthat::assert_that(nchar(system.file(package = "Rtsne")) > 0, msg = "Tsne plots require that the package `Rtsne` is installed. Run install.packages('Rtsne')")

  message("Using ZscoreRobust column in TSNE dimensionality reduction")

  #Filter for rank of interest
  tsne_input_dt = kraken_report_df[Rank=="S", .(SampleID, ZscoreRobust, ScientificName)]

  #Filter zerovariance groups
  tsne_input_dt[, `:=`(novariance = data.table::uniqueN(RPM) == 1), by = ScientificName]
  tsne_input_dt = tsne_input_dt[novariance == FALSE]

  # Remove samples with infinite values ((alternative - maybe turn those inf values into the max Zscore observed in DF))
  tsne_input_noinfinites_dt = tsne_input_dt %>%
    #dtplyr::lazy_dt(immutable = FALSE) %>%
    dplyr::group_by(ScientificName) %>%
    dplyr::filter(!any(is.infinite(ZscoreRobust))) %>%
    data.table::as.data.table()

  # convert to wide form
  tsne_input_wide_dt = data.table::dcast(tsne_input_dt, SampleID ~ ScientificName, value.var = "ZscoreRobust")


  # Convert to matrix
  tsne_input_wide_mx = as.matrix(tsne_input_wide_dt[,!"SampleID"])

  # Add SampleID as rownmaes
  rownames(tsne_input_wide_mx) <- tsne_input_wide_dt$SampleID
  #return(tsne_input_dt)
}


kraken_report_visualise_single_sample_sunbursts <- function(kraken_report_df, sample_id, ranks = c("G", "S")){
  # rlang::check_installed("taxizedbextra", reason = "to use `taxizedbextra::taxid2lineage()`")
  # rlang::check_installed("sunburst", reason = "to use `lineage2sunburst()`")


  kraken_report_df |>
    dplyr::filter(SampleID == sample_id, Rank %in% ranks) |>
    dplyr::collect() |>
    dplyr::mutate(
      parent =
    )
}

#' Annotate kraken dataframes with parent taxids
#'
#' @param kraken_report_df the kraken dataframe produced by \link{kraken_reports_parse} or link{kraken_report_parse}.
#' Can be filtered to include only ranks of interest (e.g. just genus and species) but must not be sorted in any way. (relies on default sample order).
#'
#' @returns a kraken data.frame with ParentTaxonomyID and ParentScientificName columns added
#' @export
#'
#' @examples
#' # Read in all kraken reports from a directory into a data.frame
#' kreport_dir <- system.file(package="krakenR", "simulated_data/simulated_kraken_reports_inc_zero_counts/")
#' df_kreports <- kraken_reports_parse(kreport_dir)
#'
#' # Annotate parent taxids (useful for sunburst visualisation)
#' df_kreports <- kraken_annotate_parents(df_kreports)
#'
kraken_annotate_parents <- function(kraken_report_df, ancestor = "Ancestor"){
  Level <- kraken_report_df$Level
  SampleID <- kraken_report_df$SampleID
  TaxonomyID <- kraken_report_df$TaxonomyID
  # ScientificName <- kraken_report_df$ScientificName
  prev_level = Level[1]
  prev_sample = SampleID[1]

  nrows = nrow(kraken_report_df)
  prev_taxid_at_level = numeric()
  parent_taxids = numeric(nrows)

  for (i in seq_len(nrows)){
    # message("Completed: ", i, "/", nrows)

    # When sample Changes reset prev_taxid_levels
    if(SampleID[i] != prev_sample){ prev_taxid_at_level <- numeric()} # If sample changed

    if(i == 1) {
      # If first row just set parent taxid to Inf
      parent_taxids[i] <- Inf
    }
    else{
      # Grab Current Level
      current_level = Level[i]

      # Lookup the last taxid which was higher ranked (lower level)
      parent_taxid = fetch_last_index_with_lower_level(prev_taxid_at_level, current_level=current_level)

      # If NULL replace with Inf
      parent_taxid <- if(is.na(parent_taxid)) Inf else parent_taxid

      # Update ls_parent_taxids
      parent_taxids[i] <- parent_taxid
    }

    # Update
    prev_taxid_at_level[Level[i]] <- TaxonomyID[i]
    prev_sample = SampleID[i]
  }

  kraken_report_df$ParentTaxonomyID <- parent_taxids

  # Add Scientific Name Lookup
  ## A) Build a unique lookup of TaxonomyID -> ScientificName
  idx_unique <- !duplicated(kraken_report_df$TaxonomyID)
  lookup_taxid <- kraken_report_df$TaxonomyID[idx_unique]
  lookup_name  <- kraken_report_df$ScientificName[idx_unique]

  # B) Match parent_taxid to this lookup
  m <- match(kraken_report_df$ParentTaxonomyID, lookup_taxid)
  kraken_report_df$ParentScientificName <- lookup_name[m]
  kraken_report_df$ParentScientificName[is.na(kraken_report_df$ParentScientificName)] <- ancestor
  return(kraken_report_df)
}

fetch_last_index_with_lower_level <- function(prev_taxid_at_level, current_level){
  levels_observed = seq_along(prev_taxid_at_level)

  # If no lower level (higher ranked) observations were seen - return NULL
  if(all(levels_observed >= current_level)) return(NaN)

  # Get the highest level thats still lower than current
  highest_valid_level <- max(levels_observed[levels_observed <  current_level])

  # Return the taxonomy ID of that level
  taxid=prev_taxid_at_level[highest_valid_level]
  return(taxid)
}


#' Visualise Kraken classifications as a sunburst plot
#'
#' Creates an interactive sunburst plot of Kraken classification results for one
#' or more samples. Each ring represents a taxonomic rank and segment sizes are
#' proportional to \code{ReadsCoveredByClade}. Parent–child relationships are
#' inferred using \code{\link{kraken_annotate_parents}}.
#'
#' @param kraken_report_df A Kraken report data.frame produced by
#'   \code{\link{kraken_reports_parse}} or \code{\link{kraken_report_parse}}.
#'   The data must contain at least the columns \code{SampleID}, \code{Rank},
#'   \code{ScientificName}, \code{TaxonomyID}, \code{Level} and
#'   \code{ReadsCoveredByClade}. The data can be pre-filtered to ranks of
#'   interest but must not be resorted from its original Kraken order.
#' @param sample Character vector giving one or more sample IDs (matching the
#'   \code{SampleID} column) to include in the plot.
#' @param ranks Character vector of rank codes to include (as stored in the
#'   \code{Rank} column), for example \code{c("F", "G", "S")} for family, genus
#'   and species.
#' @param ancestor Character scalar used as the label for the inferred root
#'   ancestor when no parent taxon exists (passed to
#'   \code{\link{kraken_annotate_parents}}). Defaults to \code{"Ancestor"}.
#' @param insidetextorientation Text orientation for labels inside the
#'   sunburst sectors. Must be one of \code{"horizontal"}, \code{"radial"} or
#'   \code{"tangential"}. Passed to \code{\link[plotly]{plot_ly}}.
#' @param return_data Instead of returning the htmlwidget plot, return the data
#' @returns A \code{plotly} object representing the sunburst plot.
#'
#' @export
#'
#' @examples
#' # Read in all kraken reports from a directory into a data.frame
#' kreport_dir <- system.file(
#'   package = "krakenR",
#'   "simulated_data/simulated_kraken_reports_inc_zero_counts/"
#' )
#' df_kreports <- kraken_reports_parse(kreport_dir)
#'
#' # Visualise the first sample as a sunburst at family, genus and species level
#' first_sample <- unique(df_kreports$SampleID)[1]
#' kraken_visualise_sunburst(
#'   df_kreports,
#'   sample = first_sample,
#'   ranks = c("F", "G", "S")
#' )
kraken_visualise_sunburst <- function(kraken_report_df, sample, ranks = c("F", "G", "S"), ancestor = "Ancestor", insidetextorientation = c("horizontal", "radial", "tangential"), return_data = FALSE){
  insidetextorientation <- rlang::arg_match(insidetextorientation)

  # Subset and Annotate
  kraken_report_annotated <- kraken_report_df |>
    dplyr::filter(SampleID %in% sample, Rank %in% ranks, ReadsCoveredByClade > 0) |>
    kraken_annotate_parents(ancestor = ancestor)

  if(return_data)
    return(kraken_report_annotated)

  # Step2: Create sunburst
  plotly::plot_ly(
    type = "sunburst",
    labels = kraken_report_annotated$ScientificName,
    parents = kraken_report_annotated$ParentScientificName,
    values = kraken_report_annotated$ReadsCoveredByClade,
    insidetextorientation = insidetextorientation
  )

}


