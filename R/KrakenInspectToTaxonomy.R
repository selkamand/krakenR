

# Kraken Taxanomic Operations ---------------------------------------------
kraken_taxonomy_get_ancestors <- function(taxonomy_database, taxids, ranks){
  taxonomy_db_conn <- DBI::dbConnect(RSQLite::SQLite(), taxonomy_database)

  #nodes_db <- dplyr::tbl(taxonomy_db_conn, "nodes")
  heirarchy_db <- dplyr::tbl(taxonomy_db_conn, "heirarchy")

  heirarchy_df <- heirarchy_db |>
    dplyr::filter(TaxonomyID %in% taxids) |>
    dplyr::collect()

  df <- heirarchy_df |>
    dplyr::group_by(TaxonomyID) |>
    dplyr::summarise(
      Lineage = paste0(ScientificNameAncestor, collapse = ">"),
      ScientificName = unique(ScientificName)
      ) |>
    dplyr::mutate(Lineage = paste0(Lineage, ">", ScientificName)) |>
    dplyr::select(-ScientificName)

  DBI::dbDisconnect(taxonomy_db_conn)
  return(df)
}

# Add a 'get all child taxids' function so we can more easily graph distributions of reads in specific clades


ancestor_table_to_lineage <- function(df){
  #TODO: Inplement functionality
}

# Generate Taxonomy Sqlite DB ---------------------------------------------

#' Parse Kraken Inspect File
#'
#' @param file path to the file produced by running 'kraken2 inspect' on database
#'
#' @return tibble representation of kraken inspect.
#' @export
#'
#' @examples
#' inspect_file = system.file(package='krakenR', "example_data/inspect_pluspf_20210517.txt")
#' kraken_inspect_parse(inspect_file)
kraken_inspect_parse_to_sqlite_db <- function(file, dbname_taxonomy = paste0(getwd(),"/kraken_taxonomic_database.sqlite"), dbname_heirarchy = paste0(getwd(),"/kraken_taxonomic_heirarchy.sqlite")){
  assertthat::assert_that(assertthat::is.string(file))
  assertthat::assert_that(file.exists(file))

  column_names = c("PercentMinimisersMappingToClade", "MinimisersMappedToClade", "MinimisersDirectlyMapped", "Rank", "TaxonomyID", "ScientificNameRaw")

  inspect_df <- utils::read.csv(
    file = file,
    sep = "\t",
    comment.char = "#",
    col.names = column_names
)
  inspect_df <- dplyr::as_tibble(inspect_df)

  required_columns = 6
  if(ncol(inspect_df) != required_columns)
    cli::cli_abort("Expected {required_columns} columns. Found {ncol(inspect_df)}")

  inspect_df[["Level"]] <- stringr::str_count(string = inspect_df[["ScientificNameRaw"]], pattern = "  ")
  inspect_df[["ScientificName"]] <- stringr::str_replace(string = inspect_df[["ScientificNameRaw"]], pattern = "^ +", replacement = "")


  inspect_df[["ParentRow"]] <- find_parent_rownumber(Lvls = inspect_df[["Level"]], Rows = seq_along(inspect_df[["ScientificName"]]))
  inspect_df[["ParentScientificName"]] <- inspect_df[["ScientificName"]][inspect_df[["ParentRow"]]]
  inspect_df[["ParentTaxonomyID"]] <-  inspect_df[["TaxonomyID"]][inspect_df[["ParentRow"]]]

  # Create heirarchy table
  heirarchy_df <- network_style_to_heirarchy_table(
    taxonomy_ids = inspect_df$TaxonomyID,
    parent_taxonomy_ids = inspect_df$ParentTaxonomyID,
    rank = inspect_df$Rank,
    level = inspect_df$Level,
    sciname = inspect_df$ScientificName
    )
  # cli::cli_alert_info("Adding taxonomyID lineage strings")
  # inspect_df[["LineageTaxonomyID"]] <- network_style_to_linage_strings(inspect_df[["TaxonomyID"]], inspect_df[["ParentTaxonomyID"]])

  # WriteDatabase
  kraken_write_taxonomy_database(
    taxonomy_df = inspect_df[c('ScientificName', 'TaxonomyID', 'Rank')],
    heirarchy_df = heirarchy_df,
    dbname = dbname_taxonomy
    )

  # add ancestors table

  return(inspect_df)
}

kraken_write_taxonomy_database <- function(taxonomy_df, heirarchy_df,dbname = paste0(getwd(),"/kraken_taxonomic_database")){
  assertthat::assert_that(all(c('ScientificName', 'TaxonomyID', 'Rank') %in% names(taxonomy_df)))
  assertthat::assert_that(all(c('TaxonomyID', 'Ancestor', 'Rank', 'Level', 'ScientificName', 'ScientificNameAncestor') %in% names(heirarchy_df)))

  assertthat::assert_that(anyDuplicated(taxonomy_df$TaxonomyID)==0)

  duplicated_name_indexes <- anyDuplicated(taxonomy_df$ScientificName)
  if(duplicated_name_indexes != 0){
    duplicated_names <- taxonomy_df[['ScientificName']][duplicated(taxonomy_df[['ScientificName']])]
    cli::cli_alert_info("Duplicate Scientific Names found. Adding suffix to make the following names unique: ")
    cli::cli_alert("{duplicated_names}")
    taxonomy_df[["ScientificName"]] <- make.unique(taxonomy_df[["ScientificName"]])
  }


  cli::cli_h1("Creating Taxonomy Database")
  cli::cli_alert_info("Writing Taxonomic Table sqlite database table to {.path {dbname}}")
  kraken_taxonomy_db <- DBI::dbConnect(RSQLite::SQLite(), dbname)

  node_table_name = "nodes"
  cli::cli_alert_info("Adding table {.strong {node_table_name}} connecting Taxonomy,ScientificName, and Rank")
  DBI::dbWriteTable(conn = kraken_taxonomy_db, name = node_table_name, value = taxonomy_df, overwrite = TRUE)

  cli::cli_alert_info("Indexing table {.strong {node_table_name}} based on TaxonomyID")
  DBI::dbSendQuery(conn = kraken_taxonomy_db, statement = paste0(
    "CREATE INDEX nodes_taxonomy_index
    ON ",node_table_name," (TaxonomyID);
    "))

  cli::cli_alert_info("Indexing table {.strong {node_table_name}} based on ScientifcName")
  DBI::dbSendQuery(conn = kraken_taxonomy_db, statement = paste0(
    "CREATE INDEX nodes_sciname_index
    ON ",node_table_name," (ScientificName);
    "))

  heirarchy_table_name = "heirarchy"
  cli::cli_alert_info("Adding table {.strong {heirarchy_table_name}} descibing all Ancestors of each Taxid (including their Ranks)")
  DBI::dbWriteTable(conn = kraken_taxonomy_db, name = heirarchy_table_name, value = heirarchy_df, overwrite = TRUE)

  cli::cli_alert_info("Indexing table {.strong {heirarchy_table_name}} based on TaxonomyID")
  DBI::dbSendQuery(conn = kraken_taxonomy_db, statement = paste0(
    "CREATE INDEX heirarchy_taxonomy_index
    ON ",heirarchy_table_name," (TaxonomyID);
    "))

  DBI::dbDisconnect(kraken_taxonomy_db)

}

find_parent_rownumber <- function(Lvls, Rows){
  rlang::check_installed("progress")

  cli::cli_alert_info("Identifying parents of each taxid in kraken database ... ")
  possible_parents = list()
  parent_rows = numeric(length(Rows))

  pb <- progress::progress_bar$new(total = length(Rows))

  for (i in seq_along(Rows)){
    #browser()
    possible_parent_levels = names(possible_parents)[as.numeric(names(possible_parents)) < Lvls[i]]

    if(length(possible_parent_levels) == 0)
      parent_row = NA_integer_
    else
      parent_row = max(unlist(possible_parents[possible_parent_levels]))



    possible_parents[[as.character(Lvls[i])]] <- Rows[i]
    parent_rows[[i]] <- parent_row

    pb$tick()
  }

  return(parent_rows)
}

network_style_to_linage_strings <- function(node_label, parent_label){
  rlang::check_installed("progress")

  assertthat::assert_that(length(node_label) == length(parent_label))
  assertthat::assert_that(anyDuplicated(node_label)==0)

  pb <- progress::progress_bar$new(total = length(node_label))

  lineage_string <- character(length(node_label))

  for (i in seq_along(node_label)){
    lineage = node_label[i]
    parent_taxid = parent_label[i]


    while (!is.na(parent_taxid) & !parent_taxid < 0){
      lineage <- c(lineage,parent_taxid)
      parent_taxid_index = match(parent_taxid, node_label)
      parent_taxid = parent_label[parent_taxid_index]
    }

    lineage_string[[i]] <- paste0(rev(lineage), collapse=">")
    pb$tick()
  }

  return(lineage_string)
}

#' Title
#'
#' @param taxonomy_ids taxid
#' @param parent_taxonomy_ids taxid parent
#' @param ranks taxid rank
#' @param level taxid level
#'
#' @return heirarchy table. For each taxid list all their ancestors
network_style_to_heirarchy_table <- function(taxonomy_ids, parent_taxonomy_ids, ranks, level, sciname){

  assertthat::assert_that(length(taxonomy_ids) == length(parent_taxonomy_ids))
  assertthat::assert_that(length(taxonomy_ids) == length(ranks))
  assertthat::assert_that(length(taxonomy_ids) == length(level))
  assertthat::assert_that(length(taxonomy_ids) == length(sciname))

  sciname <- make.unique(sciname)

  cli::cli_alert_info("Generating Heirarchy Table")
  pb <- progress::progress_bar$new(total = length(taxonomy_ids))

  heirarchy_df <- purrr::map_dfr(seq_along(taxonomy_ids), .f = function(i){

    lineage = numeric()
    parent_taxid = parent_taxonomy_ids[i]

    while (!is.na(parent_taxid) & !parent_taxid < 0){
      lineage <- c(lineage,parent_taxid)
      parent_taxid_index = match(parent_taxid, taxonomy_ids)
      parent_taxid = parent_taxonomy_ids[parent_taxid_index]

    }
    pb$tick()
    n_ancestors =  length(lineage)
    data.frame(TaxonomyID=rep(taxonomy_ids[i], times=n_ancestors), Ancestor=lineage, Rank = rep(ranks[i], times=n_ancestors))
    #lineage_string[[i]] <- paste0(rev(lineage), collapse=">")
    })

  heirarchy_df[["Level"]] <-level[match(heirarchy_df[["Ancestor"]], taxonomy_ids)]
  heirarchy_df[["ScientificName"]] <-sciname[match(heirarchy_df[["TaxonomyID"]], taxonomy_ids)]
  heirarchy_df[["ScientificNameAncestor"]] <-sciname[match(heirarchy_df[["Ancestor"]], taxonomy_ids)]

  return(heirarchy_df)
}

taxid_lineage_string_to_scientific_name <- function(taxid_lineage_string, taxonomy_ids, scinames, ranks,  ranks_to_include = c("Unclassified",  "Root",  "Domain",  "Kingdom",  "Phylum",  "Class",  "Order",  "Family",  "Genus",  "Species")){
  assertthat::assert_that(length(taxonomy_ids) == length(scinames))
  assertthat::assert_that(length(scinames) == length(ranks))


  valid_ranks = c(
    'U'="Unclassified",
    'R'="Root",
    'D'="Domain",
    'K'="Kingdom",
    'P'="Phylum",
    'C'="Class",
    'O'="Order",
    'F'="Family",
    'G'="Genus",
    'S'="Species"
  )


  # Assertions
  assertthat::assert_that(all(ranks_to_include %in% valid_ranks),
                          msg = utilitybeltfmt::fmterror(
                            "Ranks must be one of ", paste0("[", valid_ranks, "]", collapse = " "),
                            "\n\nInvalid Ranks: \n", paste0(ranks_to_include[!ranks_to_include %in% valid_ranks])
                          ))

  ranks_to_include_shorthand = names(valid_ranks)[match(ranks_to_include, valid_ranks)]


  lineage_string_split <- strsplit(x = taxid_lineage_string, split = ">")

  pb = progress::progress_bar$new(total = length(taxid_lineage_string))

  purrr::map_chr(lineage_string_split, .f = function(taxids){
    taxids <- as.numeric(taxids)
    taxid_rownumbers = match(taxids, taxonomy_ids)
    scinames_ = scinames[taxid_rownumbers]
    ranks_ = gsub(ranks[taxid_rownumbers], pattern = "[0-9]", replacement = "")
    #ranks_to_incluede = ranks_[ranks_ %in% ranks_to_include_shorthand]
    scinames_to_include <- scinames_[ranks_ %in% ranks_to_include_shorthand | taxids == taxids[length(taxids)]]
    pb$tick()

    paste0(scinames_to_include, collapse=">")
  })

}
