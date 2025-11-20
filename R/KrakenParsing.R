
# Reports -----------------------------------------------------------------


#' Parse Kraken Reports
#'
#' @param path_to_kreport path to kraken report (string)
#' @param sample_id  sample identifier. By default will guess from the filename (takes everything before the first '.' as sample name)
#' @param verbose print more informative messages (boolean)
#'
#' @return a dataframe describing all samples kraken reports
#' @export
#'
kraken_report_parse <- function(path_to_kreport, sample_id = NULL, verbose = TRUE) {

  kreport_basename = basename(path_to_kreport)

  assertthat::assert_that(file.exists(path_to_kreport), msg = paste0("Could not find file: ", path_to_kreport))

  kreport_headings <- c("PercentReadsCoveredByCladeLowResolution","ReadsCoveredByClade", "ReadsDirectlyAssigned", "Rank", "TaxonomyID", "ScientificName")

  # Read Files into data.table
  kraken_reports_df <- data.table::fread(path_to_kreport, col.names = kreport_headings, sep="\t", strip.white = FALSE)

  #Use number of indents in scientific name to extrapolate depth
  kraken_reports_df[, `:=` (Level = stringr::str_count(ScientificName, pattern = "  "))]

  # Strip tabs from scientific name column
  kraken_reports_df[, `:=` (ScientificName = stringr::str_replace(string = ScientificName, pattern = "^ +", replacement = ""))]
  #kraken_df$ScientificName <- sub(x = kraken_df$ScientificName, pattern = "^ +", replacement = "")

  # Simplify Rank code - In case we don't care about differentiating between S, S1, or S2 ranks - theyre all species level -- most of the time original rankings should be fine
  kraken_reports_df[, `:=` (RankSimple = stringr::str_replace_all(string =  Rank, pattern = "[0-9]", replacement = ""))]

  if(is.null(sample_id)){
    if (verbose) message("No Sample ID supplied (sample_id = NULL) ... Guessing SampleID from filename (assuming everything before first '.' is the sample ID)")
    #Add SampleID column derived from report filename
    kraken_reports_df[, `:=`(SampleID = stringr::str_replace(string = kreport_basename, pattern = '\\..*', replacement = ""))]
  }
  else{
    assertthat::assert_that(assertthat::is.string(sample_id), msg = paste0("User-specified sample identifier must be a string, not a [", class(sample_id) ,"]"))
    kraken_reports_df[, `:=`(SampleID = sample_id)]
  }

  # Calculate RPM (reads covered by clade per million total reads
  kraken_reports_df[, `:=`(TotalReadsInSample = sum(ReadsDirectlyAssigned)), by = .(SampleID)]
  kraken_reports_df[, `:=`(RPM = ReadsCoveredByClade * 1e+06/TotalReadsInSample)]

  return(kraken_reports_df)
}

#' Parse Kraken Reports
#'
#' To parse a single kraken report - see
#'
#' @param kraken2directory path to a directory filled with ONLY kraken2 reports
#'
#' @return a dataframe describing all samples kraken reports
#' @export
#'
#' @examples
#' directory = system.file("simulated_data/simulated_kraken_reports", package = "krakenR")
#' kraken_reports_parse(directory)
#'
kraken_reports_parse <- function(kraken2directory){

  assertthat::assert_that(dir.exists(kraken2directory), msg = paste0("Could not find directory: ", kraken2directory))
  kreport_headings <- c("PercentReadsCoveredByCladeLowResolution","ReadsCoveredByClade", "ReadsDirectlyAssigned", "Rank", "TaxonomyID", "ScientificName")
  krakenreportpaths <- dir(path = kraken2directory, full.names = TRUE, recursive = FALSE, include.dirs = FALSE)

  # Read Files into one huge df
  kraken_reports_df <- data.table::rbindlist(
    lapply(krakenreportpaths, data.table::fread, col.names = kreport_headings, sep="\t", strip.white = FALSE) %>%
      magrittr::set_names(basename(krakenreportpaths)),
    idcol = "Filename")


  #Use number of indents in scientific name to extrapolate depth
  kraken_reports_df[, `:=` (Level = stringr::str_count(ScientificName, pattern = "  "))]


  # Strip tabs from scientific name column
  kraken_reports_df[, `:=` (ScientificName = stringr::str_replace(string = ScientificName, pattern = "^ +", replacement = ""))]
  #kraken_df$ScientificName <- sub(x = kraken_df$ScientificName, pattern = "^ +", replacement = "")

  # Simplify Rank code - In case we don't care about differentiating between S, S1, or S2 ranks - theyre all species level -- most of the time original rankings should be fine
  kraken_reports_df[, `:=` (RankSimple = stringr::str_replace_all(string =  Rank, pattern = "[0-9]", replacement = ""))]

  #Add SampleID column derived from report filename
  kraken_reports_df[, `:=`(SampleID = stringr::str_replace(string = Filename, pattern = '\\..*', replacement = ""))]

  # Assert that theres no files describing the same sample.
  n_samples = dplyr::n_distinct(kraken_reports_df[,SampleID])
  n_files = dplyr::n_distinct(kraken_reports_df[,Filename])
  assertthat::assert_that(n_files == n_samples, msg = paste0("The number of files [", n_files ,"] is not the same as the number of distinct sample IDs [", n_samples,"].  Sample IDs are extrapolated from filenames, so please ensure files are named appropriately (i.e. filenames start with a unique sampleID, where the end of the sampleID is indicated by a period. e.g. `sample1.kreport`)"))

  # Calculate RPM (reads covered by clade per million total reads
  kraken_reports_df[, `:=`(TotalReadsInSample = sum(ReadsDirectlyAssigned)), by = .(SampleID)]
  kraken_reports_df[, `:=`(RPM = ReadsCoveredByClade * 1e+06/TotalReadsInSample)]


  return(kraken_reports_df[])
}

#' Identify Taxid Descendancy Status
#'
#' Identify taxids in your \strong{kraken_report_df} that are either equal to the user specified taxonomy ID or one of its descendants/children
#' This is useful for telling what species are bacterial / viral / belong to a particular genus etc.
#' \strong{WARNING: } this function relies on inherent properties of kraken reports - make sure you run this ONLY on unsorted dataframes produced by \link{kraken_reports_parse} or link{kraken_report_parse}
#'
#'
#' @param kraken_report_df an UNSORTED kraken report dataframe produced by \link{kraken_reports_parse} or link{kraken_report_parse}
#' @param taxonomy a taxonomy object returned by  [kraken_taxonomy_parse()]
#' @param taxid ncbi taxonomic id (integer). If you supply a numeric vector then we check whether taxids belong to ANY of the ones supplied.
#' @param colname name of the column to add to the dataframe
#' @param inclusive include the supplied taxid.
#'
#' @return This function returns the same dataframe from \strong{kraken_report_df} produced by \link{kraken_reports_parse} or link{kraken_report_parse}  but with  a new columns describing inclusive descendency status to a particular taxid (logical vector)
#'
#' @examples
#' # Read in all kraken reports from a directory into a data.frame
#' kreport_dir <- system.file(package="krakenR", "simulated_data/simulated_kraken_reports_inc_zero_counts/")
#' df_kreports <- kraken_reports_parse(kreport_dir)
#'
#' # Read in Kraken Taxonomy
#' path_ktaxonomy <- system.file("example_data/ktaxonomy.tsv", package = "krakenR")
#' taxonomy <- kraken_taxonomy_parse(path_ktaxonomy)
#'
#' # Add descendency status
#' kraken_report_add_descendancy_status(df_kreports, taxonomy, taxid=2, colname = "is_bacterial", inclusive = TRUE)
#'
#' @export
kraken_report_add_descendancy_status <- function(kraken_report_df, taxonomy, taxid, colname="is_descendent", inclusive = TRUE){

  ls_child_taxids <- lapply(unique(taxid), function(t){ kraken_fetch_child_taxids(taxonomy, taxid = taxid, inclusive=inclusive)})
  child_taxids <- unique(unlist(ls_child_taxids))

  kraken_report_df[[colname]] <- kraken_report_df[["TaxonomyID"]] %in% child_taxids

  return(kraken_report_df)
}


# kraken_fetch_child_taxids <- function(kraken_report_df, taxid, inclusive = TRUE){
#
#   levels  <- kraken_report_df$Level
#   taxids  <- kraken_report_df$TaxonomyID
#
#   # Locate first occurrence
#   start <- match(taxid, taxids)
#   if (is.na(start)) return(numeric(0))
#
#   target_level <- levels[start]
#   last_low_level_was_target <- TRUE
#
#   n <- length(taxids)
#   is_child <- logical(n)   # mask to mark children
#
#   for (i in start:n) {
#
#     if (taxids[i] == taxid) {
#       # reset the chain whenever we see the target taxid again
#       last_low_level_was_target <- TRUE
#       next
#     }
#
#     lvl <- levels[i]
#
#     if (lvl > target_level && last_low_level_was_target) {
#       # mark as child
#       is_child[i] <- TRUE
#     } else if (lvl <= target_level) {
#       # chain broken
#       last_low_level_was_target <- FALSE
#     }
#   }
#
#   out <- taxids[is_child]
#   if (inclusive) out <- c(taxid, out)
#   return(out)
# }


#' Quantifying Signal Spread Across Taxid
#'
#' This function tells you the proportion of reads belonging to a taxid (or its descendents) whose classification within the larger group is driven by classification into the most commonly classified species/genus/family (user can choose).
#' Say we notice a high number of bacterial reads and we want to know if this phenomina is driven largely by a single species or a more diffuse binning of reads into lots of different bacterial species.
#' The low spread binning (most in single species) would increase our confidence that the sample actually has the in it. If binning is more diffuse its more likely to have some other explanation.
#'
#' @param parent_taxid the taxid representing a larger group within which we want to know whether binning of reads is diffuse or focused (integer)
#' @param n  use the n most binned species/genus'/etc to calculate 'proportion of reads in parent_taxid lineage explained by a focused group' (integer)
#' @param rank one of: (S, G, F, O, C, P, K). S=speies, G=genus, etc. Determines which rank in calculation of binning 'diffuseness'. You will almost always want to use species (S), maybe genus.
#' @inheritParams kraken_report_add_zscore
#'
#' @return Dataframe describing the percentage of reads in parent clade that are explained by the n most frequently classified taxids
#' @export
kraken_calculate_proportion_of_signal_explained_by_n_strongest_taxids <- function(kraken_report_df,parent_taxid, n=1, rank = "S"){
  assertthat::assert_that(assertthat::is.number(parent_taxid))
  assertthat::assert_that(is.data.frame(kraken_report_df))
  assertthat::assert_that(assertthat::has_name(kraken_report_df, c("TaxonomyID", "Level", "SampleID", "ReadsCoveredByClade", "ScientificName")))
  assertthat::assert_that(parent_taxid %in% kraken_report_df[["TaxonomyID"]])
  assertthat::assert_that(dplyr::n_distinct(table(kraken_report_df[["TaxonomyID"]])) == 1, msg = "Some Taxonomy Ids appear more than others. Its likely you haven't used --report-zero-counts. I would advise you add this flag as it makes it very easy to catch problems arising from different databases. All the functions in this package will still work fine - just run this code again with the option careful=FALSE")

  taxid_name = kraken_report_df[["ScientificName"]][match(parent_taxid, kraken_report_df[["TaxonomyID"]])]

  kraken_report_with_descendancy_status = kraken_report_add_descendancy_status(kraken_report_df = kraken_report_df, taxid = parent_taxid, columname = "InParentClade", verbose=F)
  reads_covered_by_clade <- kraken_report_with_descendancy_status %>%
    dplyr::filter(TaxonomyID == parent_taxid) %>%
    dplyr::distinct(SampleID, ReadsCoveredByClade,.keep_all = TRUE) %>%
    dplyr::select(SampleID, ReadsCoveredByParentClade = ReadsCoveredByClade)

  kraken_report_with_descendancy_status %>%
    dplyr::filter(InParentClade==TRUE,Rank == rank) %>%
    dplyr::group_by(SampleID) %>%
    dplyr::mutate(ContributionRank = bit64::rank(ties.method="first", -ReadsCoveredByClade)) %>%
    dplyr::summarise(
      TopRankClade = paste0(ScientificName[ContributionRank<=n], collapse = ","),
      TopRankCladeTaxid = paste0(TaxonomyID[ContributionRank<=n], collapse = ","),
      ReadsCoveredByTopRankClade = sum(ReadsCoveredByClade[ContributionRank<=n]),
      #ReadsCoveredByParentClade = reads_covered_by_clade$ReadsCoveredByClade[match(SampleID, reads_covered_by_clade$SampleID)],
      #PercentageOfReadsInParentCladeCoveredByTopRankClade = ReadsCoveredByTopRankClade*100/ReadsCoveredByParentClade
      ) %>%
    dplyr::left_join(reads_covered_by_clade, by="SampleID") %>%
    dplyr::mutate(
      ParentCladeQueried = parent_taxid,
      Nused = n,
      PercentageOfReadsInParentCladeCoveredByTopRankClade = ReadsCoveredByTopRankClade*100/ReadsCoveredByParentClade,
      RPM = ReadsCoveredByTopRankClade*1000000/ReadsCoveredByParentClade
      )
}


# Inspect -----------------------------------------------------------------


#' Parse a kraken inspect file
#'
#' @param path Path to file produced by k2 inspect
#'
#' @returns a data.table representation of the kraken inspect file.
#' @export
#'
#' @examples
#' path <- system.file("example_data/inspect_pluspf_20210517.txt", package = "krakenR")
#' kraken_inspect_parse(path)
kraken_inspect_parse <- function(path){

  dt <- data.table::fread(
    file = path,
    header = FALSE,
    sep="\t", strip.white = FALSE,
    col.names = c("PercentMinimisersCoveredByClade", "MinimisersCoveredByClade", "MinimisersDirectlyAssigned", "Rank", "TaxonomyID", "ScientificName")
  )



  #Use number of indents in scientific name to extrapolate depth
  dt[, `:=` (Level = stringr::str_count(ScientificName, pattern = "  "))]


  # Strip tabs from scientific name column
  dt[, `:=` (ScientificName = stringr::str_replace(string = ScientificName, pattern = "^ +", replacement = ""))]

  dt[]
}

#' Taxids of Interest
#'
#'
#' @param purpose one of ("info","general", "Detecting Infection from Human Samples", "Detecting Viruses Potentially Related to Paediatric Cancer")
#'
#' @return Named numeric vector with some commonly used taxids (returns invisible vector)
#' @export
#'
#' @examples
#' kraken_info_taxids_of_interest()
kraken_info_taxids_of_interest <- function(purpose = "info"){
  purpose_possible_values = c("info","general", "Detecting Infection from Human Samples", "Detecting Viruses Potentially Related to Paediatric Cancer")
  assertthat::assert_that(purpose %in% purpose_possible_values, msg = paste0("[purpose] must be one of [", paste0(purpose_possible_values, collapse = ", "), "]"))

  if(purpose == "info"){
    message("Please run one of the following based on your goal: \n\n",
            'taxids <- kraken_info_taxids_of_interest("general")\n\n',
            'taxids <- kraken_info_taxids_of_interest("Detecting Infection from Human Samples")\n\n',
            'taxids <- kraken_info_taxids_of_interest("Detecting Viruses Potentially Related to Paediatric Cancer")\n\n'
            )
  }
  else if(purpose == "general"){
    message("Listing General Purpose Taxids")
    general_vector <- c(
      "unclassified" = 0,
      "root (a.k.a classified)" = 1,
      "Bacteria" = 2,
      "Viruses" = 10239,
      "Fungi" = 4751,
      "Archaea" = 2157,
      "SAR" = 2698737,
      "Homo sapiens" = 9606,
      "Homo" = 9605,
      "Hominidae" = 9604,
      "Chordata" = 7711,
      "Other Sequences" = 28384,
      "Escherichia virus phiX174" = 10847
      #"Eukaryota" = 2759,
      )

    message(paste0("[",general_vector, "]", "\t\t",names(general_vector), collapse = "\n"))
    message("\nNotes:\n> SAR is a subdomain of eukaryotic microbes that includes many eukaryotic microbes including protists")
    message("> Bacteria, Viruses, and Archaea are Domains, while Fungi is a kingdom within the eukarya domain")
    message("> Chordata is a phylum containing humans and many other animals - no microbes")
    message("> If Using PlusPF DB: Other sequences refers to UniVecCore sequences (common artificial contaminants, plasmids, etc)\n  and includes PhiX")
    return(invisible(general_vector))
  }
  else if(purpose == "Detecting Viruses Potentially Related to Paediatric Cancer"){
    pedcanviral = c(
      "Human gammaherpesvirus 4 (EBV)" = 10376,
      "Human betaherpesvirus 6B (HHV-6B)" = 32604,
      "Hepatitis B virus (Hep B)" = 10407,
      "Torque teno virus (TTV)" = 68887,
      "Human papillomavirus (HPV)" = 10566,
      "Cytomegalovirus (CMV)" = 10358,
      "Macaca mulatta polyomavirus 1 (SV40)" = 1891767,
      "Human polyomavirus 1 (BKV)" = 1891762
      )
    message(paste0("[",pedcanviral, "]", "\t\t",names(pedcanviral), collapse = "\n"))
    return(invisible(pedcanviral))
  }
  else if (purpose == "Detecting Infection from Human Samples"){
    humaninfectiondetection = c(
      "unclassified" = 0,
      "Homo" = 9605,
      "Viruses" = 10239,
      "Bacteria"  = 2,
      "Fungi" = 4751,
      "SAR" = 2698737,
      message("\nNotes:\n> SAR is a subdomain of eukaryotic microbes that includes many eukaryotic microbes including protists (but also includes non-microbial organisms like kelp)")
      )
    message(paste0("[",humaninfectiondetection, "]", "\t\t",names(humaninfectiondetection), collapse = "\n"))
    return(invisible(humaninfectiondetection))
  }
  else
    stop("purpose must be one of: ",paste0(purpose_possible_values, collapse=", "))
}





# Taxonomy ----------------------------------------------------------------


#' Parse a KrakenTools ktaxonomy file
#'
#' Parse the tabular taxonomy file produced by
#' \code{KrakenTools::make_ktaxonomy.py}, or the equivalent
#' \code{ktaxonomy.tsv} files distributed with Kraken databases (e.g. Ben
#' Langmead's Kraken databases built after 2023).
#'
#' The ktaxonomy format uses a multi-character field separator
#' \code{"\\t|\\t"}. This function normalises that separator to a single
#' tab and then reads the file via \code{data.table::fread()}.
#'
#' @param path_ktaxonomy Path to a ktaxonomy file (typically named
#'   \code{"ktaxonomy.tsv"} or similar).
#'
#' @returns
#' A \code{data.table} with one row per taxon, containing the columns:
#' \itemize{
#'   \item \code{TaxonomyID}: NCBI taxonomy ID.
#'   \item \code{ParentTaxonomyID}: parent NCBI taxonomy ID.
#'   \item \code{Rank}: taxonomic rank code/string (e.g. \code{"S"},
#'     \code{"G"}, \code{"F"}, etc.).
#'   \item \code{Level}: an integer or numeric level/depth indicator as
#'     stored in the ktaxonomy file.
#'   \item \code{ScientificName}: scientific name of the taxon.
#' }
#'
#' @export
#'
#' @examples
#' path <- system.file("example_data/ktaxonomy.tsv", package = "krakenR")
#' taxonomy_table <- kraken_taxonomy_parse_to_dt(path)
#' head(taxonomy_table)
kraken_taxonomy_parse_to_dt <- function(path_ktaxonomy){
  assertions::assert_file_exists(path_ktaxonomy)
  vec_ktaxonomy <- readLines(path_ktaxonomy)
  vec_ktaxonomy <- gsub(x=vec_ktaxonomy, pattern = "\t\\|\t", replacement = "\t")
  dt <- data.table::fread(text = vec_ktaxonomy, sep = "\t", col.names = c("TaxonomyID", "ParentTaxonomyID", "Rank", "Level", "ScientificName"))

  # dt$ParentTaxonomyID <- ifelse(dt$ParentTaxonomyID == dt$TaxonomyID, Inf, dt$ParentTaxonomyID)
  return(dt[])
}


#' Parse a KrakenTools ktaxonomy file
#'
#' Parse the tabular taxonomy file produced by
#' \code{KrakenTools::make_ktaxonomy.py}, or the equivalent
#' \code{ktaxonomy.tsv} files distributed with Kraken databases (e.g. Ben
#' Langmead's Kraken databases built after 2023).
#'
#' The ktaxonomy format uses a multi-character field separator
#' \code{"\\t|\\t"}. This function normalises that separator to a single
#' tab and then reads the file via \code{data.table::fread()}.
#'
#' @param path_ktaxonomy Path to a ktaxonomy file (typically named
#'   \code{"ktaxonomy.tsv"} or similar).
#'
#' @returns
#' A directed \code{igraph} network (parent->child) with node attributes:
#' \itemize{
#'   \item \code{TaxonomyID}: NCBI taxonomy ID.
#'   \item \code{ParentTaxonomyID}: parent NCBI taxonomy ID.
#'   \item \code{Rank}: taxonomic rank code/string (e.g. \code{"S"},
#'     \code{"G"}, \code{"F"}, etc.).
#'   \item \code{Level}: an integer or numeric level/depth indicator as
#'     stored in the ktaxonomy file.
#'   \item \code{ScientificName}: scientific name of the taxon.
#' }
#'
#' @export
#'
#' @examples
#' path <- system.file("example_data/ktaxonomy.tsv", package = "krakenR")
#' taxonomy <- kraken_taxonomy_parse(path)
#' print(taxonomy)
kraken_taxonomy_parse <- function(path_ktaxonomy){
  kraken_taxonomy_to_igraph(kraken_taxonomy_parse_to_dt(path_ktaxonomy))
}

kraken_taxonomy_to_igraph <- function(taxonomy_table){
  # taxonomy_table$TaxonomyID <- as.numeric(taxonomy_table$TaxonomyID)
  # taxonomy_table$ParentTaxonomyID <- as.numeric(taxonomy_table$ParentTaxonomyID)

  # Drop any cases where parent taxonomy = taxonomy ID (just removes the root parent)
  dt_edges <- taxonomy_table[taxonomy_table$ParentTaxonomyID != taxonomy_table$TaxonomyID, c("ParentTaxonomyID", "TaxonomyID")]

  igraph::graph_from_data_frame(
    d = dt_edges,
    vertices = taxonomy_table[,c("TaxonomyID", "Rank", "Level", "ScientificName")],
    directed = TRUE
  )
}

igraph_to_taxonomy_table <- function(taxonomy){
  ls_graph <- igraph::as_data_frame(taxonomy, what = "both")
  dt_edges <- data.table::as.data.table(ls_graph$edges)
  dt_edges <- dt_edges[,c(2, 1)]

  colnames(ls_graph$vertices) <- c("ParentTaxonomyID", "ParentRank", "ParentLevel", "ParentScientificName")
  df_vertex_annotation_direct <- ls_graph$vertices
  colnames(df_vertex_annotation_direct) <- c("TaxonomyID", "Rank", "Level", "ScientificName")

  colnames(dt_edges) <- c("TaxonomyID", "ParentTaxonomyID")
  dt_edges <- dt_edges |>
    dplyr::left_join(df_vertex_annotation_direct, by = "TaxonomyID") |>
    dplyr::left_join(ls_graph$vertices, by = "ParentTaxonomyID")

  return(dt_edges)


}

#' Fetch Child Taxids
#'
#' Identify taxids in your \strong{kraken_report_df} that are either equal to the user specified taxonomy ID or one of its descendants/children
#' This is useful for telling what species are bacterial / viral / belong to a particular genus etc.
#' \strong{WARNING: } this function relies on inherent properties of kraken reports - make sure you run this ONLY on unsorted dataframes produced by \link{kraken_reports_parse} or link{kraken_report_parse}
#'
#'
#' @param taxonomy a taxonomy object produce by [kraken_taxonomy_parse()]
#' @param taxid ncbi taxonomic id (integer)
#' @param inclusive include the supplied taxid.
#' @param verbose print more informative messages
#'
#' @return This function returns the same dataframe from \strong{kraken_report_df} produced by \link{kraken_reports_parse} or link{kraken_report_parse}  but with  a new columns describing inclusive descendancey status to a particular taxid
#' @export
kraken_fetch_child_taxids <- function(taxonomy, taxid, inclusive = FALSE){
  # Find all descendants of taxid
  indexes <- as.vector(igraph::subcomponent(taxonomy, v=taxid, mode = "out"))
  descendants <- igraph::vertex_attr(graph = taxonomy, name = "names", index = indexes)

  if(!inclusive) descendants <- descendants[descendants!=taxid]
  return(descendants)
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
#' taxonomy_path <- system.file(package="krakenR", "example_data/ktaxonomy.tsv")
#'
#' # Read in taxonomy
#' df_kreports <- kraken_reports_parse(kreport_dir)
#' taxonomy <- kraken_taxonomy_parse(taxonomy_path)
#'
#' # Filter for a specific sample
#' df_kreports <- df_kreports[df_kreports$SampleID == "e_coli_1",]
#'
#' # Annotate parent taxids (useful for sunburst visualisation)
#' df_kreports_annotated <- kraken_annotate_parents(df_kreports, taxonomy)
#'
kraken_annotate_parents <- function(kraken_report_df, taxonomy, ancestor = "Ancestor"){
  df_taxonomies <- igraph_to_taxonomy_table(taxonomy)
  df_taxonomies$TaxonomyID <-  as.numeric(df_taxonomies$TaxonomyID)
  df_taxonomies$Rank <- NULL
  df_taxonomies$Level <- NULL
  df_taxonomies$ScientificName <- NULL

  df_annotated <- dplyr::left_join(kraken_report_df, df_taxonomies, by = "TaxonomyID")
  df_annotated$ParentScientificName[is.na(df_annotated$ParentScientificName)] <- ancestor
  return(df_annotated)
}


#' Map taxonomy IDs to parents at specified ranks
#'
#' Given a taxonomy table parsed by \code{\link{kraken_taxonomy_parse}},
#' compute, for each taxon, the closest ancestor whose rank is in a
#' user-supplied set of ranks.
#'
#' @export
#'
#' @examples
#' # Read in all kraken reports from a directory into a data.frame
#' kreport_dir <- system.file(package="krakenR", "simulated_data/simulated_kraken_reports_inc_zero_counts/")
#' taxonomy_path <- system.file(package="krakenR", "example_data/ktaxonomy.tsv")
#'
#' # Read in taxonomy
#' df_kreports <- kraken_reports_parse(kreport_dir)
#' taxonomy <- kraken_taxonomy_parse(taxonomy_path)
#'
#' # Filter for a specific sample
#' df_kreports <- df_kreports[df_kreports$SampleID == "e_coli_1",]
#'
#' # Annotate parent taxids (useful for sunburst visualisation)
#' df_kreports_annotated <- kraken_annotate_parents_(df_kreports, taxonomy)
#'
kraken_annotate_parents_at_ranks <- function(kraken_report_df,  taxonomy,
                                            ranks = c("S", "G", "F", "K"), ancestor = "Ancestor") {

  kraken_report_df$ParentTaxonomyID <- find_parent_with_rank(
    taxonomy,
    taxids = kraken_report_df$TaxonomyID,
    ranks = ranks,
    noparent = NA_character_
  )

  df_taxonomy <- igraph_to_taxonomy_table(taxonomy)
  df_taxonomy <- df_taxonomy[,c("TaxonomyID", "Rank", "Level", "ScientificName")]
  colnames(df_taxonomy) <- paste0("Parent",colnames(df_taxonomy))
  kraken_report_df <- dplyr::left_join(kraken_report_df, df_taxonomy, by = "ParentTaxonomyID")
  kraken_report_df$ParentScientificName[is.na(kraken_report_df$ParentScientificName)] <- ancestor

  # taxonomy_reduced <- collapse_taxonomy_to_ranks(taxonomy, ranks)
  # df_annotated <- kraken_annotate_parents(kraken_report_df, taxonomy_reduced)
  # df_annotated$ParentTaxonomyID[is.na(df_annotated$ParentTaxonomyID)] <- ancestor

  return(kraken_report_df)
}

#' Collapse a Taxonomy Graph to Selected Ranks
#'
#' Given a full taxonomy stored as a directed \code{igraph} where edges point
#' \emph{parent → child}, this function constructs a reduced taxonomy graph
#' containing only vertices whose rank is in \code{valid_ranks}. Intermediate
#' vertices (those with ranks not in \code{valid_ranks}) are removed and their
#' structural role is preserved by connecting each retained vertex directly to
#' its nearest retained ancestor.
#'
#' This effectively creates a "rank-skeletonized" taxonomy: a smaller graph
#' that preserves ancestor–descendant relationships among selected ranks.
#'
#' @param g A directed \code{igraph} object representing a taxonomy. Vertices
#'   must have a rank attribute (default: \code{"rank"}), and vertex names
#'   (\code{V(g)$name}) are treated as taxonomic IDs.
#' @param valid_ranks A character vector of ranks to retain (e.g. \code{c("S","F")}).
#'   Only these vertices remain in the reduced graph.
#' @param name_attr Name of the vertex attribute storing taxids (default \code{"name"}).
#' @param rank_attr Name of the vertex attribute storing taxonomic rank codes
#'   (default \code{"rank"}).
#'
#' @details
#' The algorithm assumes the taxonomy is a tree or forest (each vertex has at
#' most one parent). It works in \eqn{O(V + E)} time and is suitable for large
#' taxonomies such as NCBI.
#'
#' Steps:
#' \enumerate{
#'   \item Identify all vertices whose rank is in \code{valid_ranks}.
#'   \item For each vertex, precompute its nearest ancestor whose rank is valid.
#'   \item For each retained vertex, determine its retained parent.
#'   \item Build a new directed graph containing only retained vertices and
#'     edges between them.
#'   \item Preserve all vertex attributes for retained vertices.
#'   \item Store \code{parent_taxid} and \code{child_taxid} as edge attributes
#'     in the reduced graph.
#' }
#'
#' @return A reduced \code{igraph} object containing:
#' \itemize{
#'   \item only vertices with ranks in \code{valid_ranks},
#'   \item edges representing direct relationships between nearest retained ancestors,
#'   \item vertex attributes copied from the original graph,
#'   \item edge attributes \code{parent_taxid} and \code{child_taxid}.
#' }
#'
#' @examples
#' library(igraph)
#'
#' # Tiny toy taxonomy:
#' edges <- data.frame(
#'   from = c("1","2","2","3","3","4"),
#'   to   = c("2","3","4","5","6","7")
#' )
#'
#' ranks <- data.frame(
#'   name = as.character(1:7),
#'   rank = c("R","F","G","G","S","S","S")
#' )
#'
#' g <- graph_from_data_frame(edges, directed = TRUE, vertices = ranks)
#'
#' # Collapse to Family (F) and Species (S)
#' g_red <- collapse_taxonomy_to_ranks(g, valid_ranks = c("S","F"))
#'
#' V(g_red)$name
#' E(g_red)$parent_taxid
#' E(g_red)$child_taxid
#'
#' @export
collapse_taxonomy_to_ranks <- function(g, valid_ranks,
                                       taxid_attr = "name",
                                       rank_attr  = "Rank") {
  stopifnot(igraph::is_igraph(g))

  # Extract attributes
  ranks  <- igraph::vertex_attr(g, rank_attr)
  taxids <- igraph::vertex_attr(g, taxid_attr)

  if (is.null(ranks) || is.null(taxids)) {
    stop("Graph must have vertex attributes '", taxid_attr, "' and '", rank_attr, "'.")
  }

  # Which vertices we keep in the reduced graph
  valid_rank <- ranks %in% valid_ranks
  if (!any(valid_rank)) {
    stop("No vertices have ranks in 'valid_ranks'.")
  }

  N <- igraph::vcount(g)
  keep_idx <- which(valid_rank)

  ## 1. Build parent index: parent_idx[child] = parent (vertex IDs, 1..N)
  el <- igraph::as_edgelist(g, names = FALSE)   # matrix: parent, child
  parent_idx <- rep(NA_integer_, N)
  if (nrow(el) > 0) {
    parent_idx[el[, 2]] <- el[, 1]
  }

  ## 2. Topological order (root -> leaves)
  topo <- as.integer(igraph::topo_sort(g, mode = "out"))

  ## 3. For each vertex, record the nearest valid-rank vertex
  ##    on the path from the root to that vertex (including itself).
  nearest_valid <- rep(NA_integer_, N)

  for (v in topo) {
    p <- parent_idx[v]
    if (is.na(p)) {
      # root
      nearest_valid[v] <- if (valid_rank[v]) v else NA_integer_
    } else {
      if (valid_rank[v]) {
        nearest_valid[v] <- v
      } else {
        nearest_valid[v] <- nearest_valid[p]
      }
    }
  }

  ## 4. For each kept vertex v, find its parent in the reduced graph:
  ##    the nearest valid ancestor of its original parent.
  red_parent <- rep(NA_integer_, N)
  for (v in keep_idx) {
    p <- parent_idx[v]
    red_parent[v] <- if (!is.na(p)) nearest_valid[p] else NA_integer_
  }

  # Edges in terms of original vertex IDs (old indices)
  from_old <- red_parent[keep_idx]
  to_old   <- keep_idx
  sel <- !is.na(from_old)
  from_old <- from_old[sel]
  to_old   <- to_old[sel]

  # Remove duplicates just in case
  if (length(from_old)) {
    edges_old <- unique(cbind(from_old, to_old))
  } else {
    edges_old <- matrix(integer(0), ncol = 2)
  }

  ## 5. Map old vertex IDs to new vertex IDs
  old_to_new <- rep(NA_integer_, N)
  old_to_new[keep_idx] <- seq_along(keep_idx)

  if (nrow(edges_old) > 0) {
    el_new <- cbind(old_to_new[edges_old[, 1]],
                    old_to_new[edges_old[, 2]])
    g_red <- igraph::graph_from_edgelist(el_new, directed = TRUE)
  } else {
    # No edges, just isolated kept vertices
    g_red <- igraph::make_empty_graph(n = length(keep_idx), directed = TRUE)
  }

  ## 6. Copy vertex attributes for the kept vertices
  for (attr_name in igraph::vertex_attr_names(g)) {
    igraph::vertex_attr(g_red, attr_name) <- igraph::vertex_attr(g, attr_name, index = keep_idx)
  }

  ## 7. Build edge data.frame: parent_taxid, child_taxid
  if (igraph::ecount(g_red) > 0) {
    el_new_numeric <- igraph::as_edgelist(g_red, names = FALSE)
    parent_taxid   <- igraph::vertex_attr(g_red, taxid_attr)[el_new_numeric[, 1]]
    child_taxid    <- igraph::vertex_attr(g_red, taxid_attr)[el_new_numeric[, 2]]
    edges_df <- data.frame(
      parent_taxid = parent_taxid,
      child_taxid  = child_taxid,
      stringsAsFactors = FALSE
    )
  } else {
    edges_df <- data.frame(
      parent_taxid = integer(0),
      child_taxid  = integer(0)
    )
  }

  # Return reduced igraph
  return(g_red)
}


find_parent_with_rank <- function(taxonomy, taxids, ranks = c("S", "G", "F", "O"), noparent = NA){
  taxids_uniq <- unique(taxids)

  taxids_in_index_order <- igraph::vertex_attr(graph = taxonomy, name = "name")
  ranks_in_index_order <- igraph::vertex_attr(graph = taxonomy, name = "Rank")
  sciname_in_index_order <- igraph::vertex_attr(graph = taxonomy, name = "ScientificName")

  vertex_ids <- match(taxids_uniq, taxids_in_index_order)
  names(vertex_ids) <- as.character(taxids_uniq)

  # Drop NAs
  vertex_ids <- na.omit(vertex_ids)

  # If no taxids are in taxonomy - retun NA for all
  if(length(vertex_ids) == 0){
    values <- rep(NA_character_, times = length(taxids))
    names(values) <- taxids
    return(values)
  }

  # Lookup first parent with valid rank
  parent_taxids <- vapply(vertex_ids, function(id){
    parent_indexes <- as.numeric(igraph::subcomponent(graph = taxonomy, v = id, mode = "in"))
    parents <- tail(x = parent_indexes, n=-1)
    ranks_of_parents <- ranks_in_index_order[parents]
    first_parent_with_valid_rank <- parents[which(ranks_of_parents %in% ranks)[1]]
    taxids_in_index_order[first_parent_with_valid_rank]
    }, FUN.VALUE = character(1)
  )
  names(parent_taxids) <- names(vertex_ids)

  # Expand back out
  parent_taxids_expanded <- parent_taxids[match(taxids, names(vertex_ids))]
  names(parent_taxids_expanded) <- taxids
  return(parent_taxids_expanded)
}


# Writing ----------------------------------------------------------------


#' Write Kraken TSV
#'
#' @inheritParams kraken_report_add_descendancy_status
#' @param filepath filename/path to write tab-separated dataframe to
#'
#' @return Run for its side effects
#' @export
#'
kraken_write_tsv <- function(kraken_report_df, filepath = paste0("kraken_report_database_",Sys.Date(), ".tsv")){
  data.table::fwrite(file = filepath, sep="\t")
}
