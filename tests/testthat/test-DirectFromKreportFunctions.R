test_that("kraken_report_get_child_taxids works", {

  testkreport_path <- system.file(package = "krakenR", "testfiles/test_ecoli_paeruginosa.kreport")

  # Expect function to run without error
  expect_error(
    kraken_report_get_child_taxids(testkreport_path, 286),
    NA
  )

  # Expect function to run without warning
expect_warning(
  kraken_report_get_child_taxids(testkreport_path, 286),
  NA
  )
#Basic inclusive query works
  expect_equal(
   object = kraken_report_get_child_taxids(testkreport_path, 286, inclusive = TRUE),
   expected = c(286, 136841, 287)
   )
  # Basic non-inclusive query works
  expect_equal(
   object = kraken_report_get_child_taxids(testkreport_path, 286, inclusive = FALSE),
   expected = c(136841, 287)
  )


  # If taxid supplied has no children in kreport, returns only itself if inclusive=TRUE
  expect_equal(
   object = kraken_report_get_child_taxids(testkreport_path, 287, inclusive = TRUE),
   expected = c(287)
  )

  # If taxid supplied has no children in kreport, throws error if inclusive=FALSE
  expect_error(
   object = kraken_report_get_child_taxids(testkreport_path, 287, inclusive = FALSE)
  )


  # [last_taxid_in_report] If taxid supplied has no children in kreport, returns only itself if inclusive=TRUE
  expect_equal(
   object = kraken_report_get_child_taxids(testkreport_path, 562, inclusive = TRUE),
   expected = c(562)
  )

  # [last_taxid_in_report] If taxid supplied has no children in kreport, throws error if inclusive=FALSE
  expect_error(
   object = kraken_report_get_child_taxids(testkreport_path, 562, inclusive = FALSE)
  )

  # [first_taxid_in_report] works as expected (inclusive)
  expect_equal(
    object = kraken_report_get_child_taxids(testkreport_path, 1, inclusive = TRUE),
    expected = c(1, 131567, 2, 1224, 1236, 72274, 135621, 286, 136841, 287, 91347, 543, 561, 562)
  )

  # [first_taxid_in_report] works as expected (not inclusive)
  expect_equal(
    object = kraken_report_get_child_taxids(testkreport_path, 1, inclusive = FALSE),
    expected = c(131567, 2, 1224, 1236, 72274, 135621, 286, 136841, 287, 91347, 543, 561, 562)
  )

  # [first_taxid_in_report] works as expected (not inclusive)
  expect_equal(
    object = kraken_report_get_child_taxids(testkreport_path, 1236, inclusive = FALSE),
    expected = c(72274, 135621, 286, 136841, 287, 91347, 543, 561, 562)
  )

  # [first_taxid_in_report] works as expected (inclusive)
  expect_equal(
    object = kraken_report_get_child_taxids(testkreport_path, 1236, inclusive = TRUE),
    expected = c(1236, 72274, 135621, 286, 136841, 287, 91347, 543, 561, 562)
  )

  # exluding children without reads directly classified works
  expect_equal(
    object = kraken_report_get_child_taxids(testkreport_path, 543, inclusive = TRUE, exclude_children_without_reads_directly_classified = TRUE),
    expected = c(543,562)
  )
  # exluding children without reads directly classified works noninclusive
  expect_equal(
    object = kraken_report_get_child_taxids(testkreport_path, 543, inclusive = FALSE, exclude_children_without_reads_directly_classified = TRUE),
    expected = c(562)
  )
})



