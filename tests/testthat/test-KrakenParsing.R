test_that("kraken_reports_parse works when reports include zero counts", {

  kraken_report_dir = system.file(package = "krakenR", "simulated_data/simulated_kraken_reports_inc_zero_counts/")

  # Runs without error
  expect_error(
    kraken_reports_parse(kraken2directory = kraken_report_dir),
    NA
  )

  # Throws error if path doesnt exist
  expect_error(
    kraken_reports_parse(kraken2directory = "asdaksdajwlkdjalwkjdawjdlakwjldawk"),
    "Could not find directory: asdaksdajwlkdjalwkjdawjdlakwjldawk"
  )


  kraken_report_dt=kraken_reports_parse(kraken2directory = kraken_report_dir)

  expect_setequal(
    object = colnames(kraken_report_dt),
    expected = c(
      "Filename",
      "PercentReadsCoveredByCladeLowResolution",
      "ReadsCoveredByClade",
      "ReadsDirectlyAssigned",
      "Rank",
      "TaxonomyID",
      "ScientificName",
      "Level",
      "RankSimple",
      "SampleID",
      "TotalReadsInSample",
      "RPM")
    )

  # Test SampleID is solid
  expect_equal(unique(kraken_report_dt$SampleID), c("e_coli_1", "e_coli_10", "e_coli_2", "e_coli_3", "e_coli_4",
                                                    "e_coli_5", "e_coli_6", "e_coli_7", "e_coli_8", "e_coli_9", "l_monocytogenes_1",
                                                    "l_monocytogenes_10", "l_monocytogenes_2", "l_monocytogenes_3",
                                                    "l_monocytogenes_4", "l_monocytogenes_5", "l_monocytogenes_6",
                                                    "l_monocytogenes_7", "l_monocytogenes_8", "l_monocytogenes_9"
  ))

  # Make sure the kraken report dataframe doesn't change between updates
  expect_snapshot_value(x = kraken_report_dt, style = "serialize")

  # ensure entries per taxid == 1 per sample (should be 20 for our simulated dataset)
  entries_per_taxid = kraken_report_dt |>
    dplyr::count(TaxonomyID) |>
    dplyr::pull(n) |>
    unique()

  expect_equal(entries_per_taxid, dplyr::n_distinct(kraken_report_dt$SampleID))

  #Expected number of rows = sample nuber * number of taxids
  expect_equal(nrow(kraken_report_dt), 20*21112)


  # ensure df is ungrouped
  expect_equal(kraken_report_dt, dplyr::ungroup(kraken_report_dt))

  # RPM calculation is correct
  expect_equal(
    kraken_report_dt[["RPM"]],
    kraken_report_dt[["ReadsCoveredByClade"]] * 1e6 / kraken_report_dt[["TotalReadsInSample"]]
  )

  # Ensure each sample has one unique value for 'TotalReadsInSample'
  expect_setequal(
    kraken_report_dt |>
      dplyr::ungroup() |>
      dplyr::group_by(SampleID) |>
      dplyr::summarise(n_different_total_reads_in_sample = dplyr::n_distinct(TotalReadsInSample)) |>
      dplyr::pull(n_different_total_reads_in_sample),
    1
  )


  # Test ZscoreRobust Calculation

})
