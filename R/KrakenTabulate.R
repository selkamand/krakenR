#' Title
#'
#' @inheritParams kraken_report_add_zscore
#' @param rank one of: (S, G, F, O, C, P, K). S=speies, G=genus, etc. Determines which rank in calculation of binning 'diffuseness'. You will almost always want to use species (S), maybe genus.
#' @param zscore_robust_min minimum robust zscore to be considered a hit at the lowest confidence level
#' @param min_reads minimum reads classified under clade required to be considered a hit at the lowest confidence level
#'
#' @return dataframe describing all hits passing \code{zscore_robust_min} and \code{min_reads} filters, with levels of confidence added (based on robust zscore and the number of reads assigned under clade
#' @export
#'
kraken_report_tabulate_hits <- function(kraken_report_df, rank = "S", zscore_robust_min = 5, min_reads = 100){
  is_db = "tbl_sql" %in% class(kraken_report_df)

  # Filter for everythin from very-low confidence to high.
  kreport_hits_df <- kraken_report_df |>
    dplyr::filter(Rank == rank, ZscoreRobust >= zscore_robust_min, ReadsCoveredByClade >= min_reads)


  #browser()
  if(is_db)
    kreport_hits_df = dplyr::collect(kreport_hits_df)

  cli::cli_alert_info(
    'Classifying hits as follows:\n
    ZscoreRobust >= 100 & ReadsCoveredByClade >= 100000 ~ "Very High",
    ZscoreRobust >= 100 & ReadsCoveredByClade >= 1000 ~ "High",
    ZscoreRobust >= 10 & ReadsCoveredByClade >= 1000 ~ "Medium",
    ZscoreRobust >= 5 & ReadsCoveredByClade >= 100 ~ "Low",
    ZscoreRobust >= {zscore_robust_min} & ReadsCoveredByClade > {min_reads} ~ "Very Low"'

    )

  kreport_hits_df <- kreport_hits_df |>
    dplyr::mutate(Confidence = dplyr::case_when(
      ZscoreRobust >= 100 & ReadsCoveredByClade >= 100000 ~ "Very High",
      ZscoreRobust >= 100 & ReadsCoveredByClade >= 1000 ~ "High",
      ZscoreRobust >= 10 & ReadsCoveredByClade >= 1000 ~ "Medium",
      ZscoreRobust >= 5 & ReadsCoveredByClade >= 100 ~ "Low",
      TRUE ~ "Very Low"
      ))

  kreport_hits_df <- kreport_hits_df |>
    dplyr::add_count(TaxonomyID, name = "Recurrence") |>
    dplyr::select(
      ScientificName, TaxonomyID, Confidence, ZscoreRobust, ReadsCoveredByClade,Recurrence,SampleID, Rank
      )


  #kraken_db |> kraken_report_tabulate_hits() |> dplyr::count(ScientificName, Confidence, sort = TRUE, name = "Recurrence")
  return(kreport_hits_df)
}


#
# kreports_to_genus_species <- function(kraken_report_df){
#   df_genus_counts <- kraken_report_df |>
#     dplyr::filter(Rank == "G") |>
#     dplyr::summarise(
#       .by = c(SampleID, TaxonomyID, ScientificName), Reads = ReadsCoveredByClade
#     )
#
#   df_species_counts <- kraken_report_df |>
#     dplyr::filter(Rank == "G") |>
#     dplyr::summarise(
#       .by = c(SampleID, TaxonomyID, ScientificName), Reads = ReadsCoveredByClade
#     )
# }
