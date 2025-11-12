# kraken_extract_reads ----------------------------------------------------
test_that("kraken_extract_reads works", {

  # Toggle whether to run tests pertaining to kraken_extract_reads.
  # I often disable these tests when developing new ones since running this particular set takes a while
  skip_tests=FALSE
  if(skip_tests) { message("\nSkipping kraken_extract_reads tests"); return(TRUE)}



  # Setup options for maximum reproducibility
  withr::local_options(.new = list(crayon.enabled = FALSE, cli.unicode = FALSE, width = 80))

  # Make a temporary_directory in which to store test ouptuts that will be deleted after test
  self_deleting_tempdir0 = withr::local_tempdir(pattern = "kraken_extract_reads")
  self_deleting_tempdir = withr::local_tempdir(pattern = "kraken_extract_reads")
  self_deleting_tempdir2 = withr::local_tempdir(pattern = "kraken_extract_reads")
  self_deleting_tempdir3 = withr::local_tempdir(pattern = "kraken_extract_reads")
  self_deleting_tempdir4 = withr::local_tempdir(pattern = "kraken_extract_reads")
  self_deleting_tempdir5 = withr::local_tempdir(pattern = "kraken_extract_reads")
  self_deleting_tempdir6 = withr::local_tempdir(pattern = "kraken_extract_reads")

  kreport_path = system.file(package="krakenR", "testfiles/e_coli_l_monocytogones_s_enterica_50seqs_each.minikraken8GB.kreport")
  kraken_stdout_path = system.file(package="krakenR", "testfiles/e_coli_l_monocytogones_s_enterica_50seqs_each.minikraken8GB.stdout")

  kraken_input_seqs_fasta = system.file(package="krakenR", "testfiles/e_coli_l_monocytogones_s_enterica_50seqs_each.fasta")
  kraken_input_seqs_invalid = "aajksdlajdlksajkldjawkljdakljwaadlkwjadwawddsalkfdsjkflenjcne"
  kraken_input_seqs_fastq_gzip = system.file(package="krakenR", "testfiles/e_coli_l_monocytogones_s_enterica_50seqs_each.fastq.gz")

  sequence_outfile_prefix0 = paste0(self_deleting_tempdir0,"/", "kraken_extract_reads")
  sequence_outfile_prefix = paste0(self_deleting_tempdir,"/", "kraken_extract_reads")
  sequence_outfile_prefix2 = paste0(self_deleting_tempdir2,"/", "kraken_extract_reads")
  sequence_outfile_prefix3 = paste0(self_deleting_tempdir3,"/", "kraken_extract_reads")
  sequence_outfile_prefix4 = paste0(self_deleting_tempdir4,"/", "kraken_extract_reads")
  sequence_outfile_prefix5 = paste0(self_deleting_tempdir5,"/", "kraken_extract_reads")
  sequence_outfile_prefix6 = paste0(self_deleting_tempdir6,"/", "kraken_extract_reads")



  # Only run full test suite if all dependencies are present ---------------------------------------
  # We only want to test read extraction on systems with the relevent dependencies:
  # (seqkit, parallel, and awk) otherwise all results will throw the 'can't find dependencies error'
  dependencies = c("awk", "parallel", "seqkit")
  if(!all(dependencies_present(dependencies))){

    #browser()
    expect_error(
      kraken_extract_reads(
        kreport_path = kreport_path,
        kraken_input_seqs = kraken_input_seqs_fasta,
        outfile_prefix = sequence_outfile_prefix,
        compress_output_seqfile = TRUE,
        taxid = 561,
        include_children = TRUE,
        kraken_stdout_path = kraken_stdout_path
      ),
      "Missing Dependencies"
      ) |> suppressMessages()

    message("\nSince dependencies for kraken_extract_reads [",paste0(dependencies[!dependencies_present(dependencies)], collapse = ", "),"] are missing we'll make sure that the functions throws the missing dependency error but skip the rest of the tests, since they'll all throw the same error")
    return(TRUE)
  }

  ## Unpaired Fasta Input ----------------------------------------------------
  # Runs without error
  #browser()
  expect_error(
    kraken_extract_reads(
      kreport_path = kreport_path,
      kraken_input_seqs = kraken_input_seqs_fasta,
      outfile_prefix = sequence_outfile_prefix,
      compress_output_seqfile = TRUE,
      taxid = 561,
      include_children = TRUE,
      kraken_stdout_path = kraken_stdout_path
    ),
    NA
  ) |> suppressMessages()

  # Produces Expected Files
  expected_outseq_path = paste0(sequence_outfile_prefix, ".reads_classified_as_taxid561_or_children.unpaired.1.fasta.gz")
  expected_readnames_path = paste0(sequence_outfile_prefix, ".reads_classified_as_taxid561_or_children.txt")
  expect_true(file.exists(expected_outseq_path))
  expect_true(file.exists(expected_readnames_path))

  # Check outseq file contents
  outseq_file_contents = readLines(expected_outseq_path)
  expected_readnames_contents = readLines(expected_readnames_path)

  # --> Is Fasta
  expect_true(grepl(outseq_file_contents[1], pattern = "^>"))

  # --> Expected Number of Reads in fasta
  expect_true(sum(grepl(outseq_file_contents, pattern = "^>")) == 23)

  # --> Expected Number of Read names classified as taxid
  expect_length(expected_readnames_contents, n = 23)

  # --> Expected contents of read names classified as taxid
  expect_equal(expected_readnames_contents, c("coli_NC_002127.1-347", "coli_NC_002695.2-213", "coli_NC_002695.2-375",
                                              "coli_NC_002695.2-860", "coli_NC_002128.1-69", "coli_NC_002127.1-931",
                                              "coli_NC_002128.1-54", "coli_NC_002127.1-625", "coli_NC_002128.1-782",
                                              "coli_NC_002127.1-954", "coli_NC_002128.1-930", "coli_NC_002695.2-111",
                                              "coli_NC_002128.1-210", "coli_NC_002695.2-673", "coli_NC_002128.1-585",
                                              "coli_NC_002695.2-1", "coli_NC_002695.2-896", "coli_NC_002695.2-174",
                                              "coli_NC_002695.2-215", "coli_NC_002127.1-166", "coli_NC_002128.1-64",
                                              "coli_NC_002695.2-796", "coli_NC_002127.1-619"))

  # Produces Expected result
  # add a browser() call - check file contents against what we expect - then write new tests to make sure everything works

  # Correctly identifies child taxids with reads directly assigned
  expect_message(
    kraken_extract_reads(
      kreport_path = kreport_path,
      kraken_input_seqs = kraken_input_seqs_fasta,
      outfile_prefix = sequence_outfile_prefix,
      taxid = 561,
      include_children = TRUE,
      kraken_stdout_path = kraken_stdout_path,
      compress_output_seqfile = TRUE
    ),
    "Taxids: 561, 562, 83334"
  ) |> suppressMessages()

  # Detects filetype
  expect_message(
    kraken_extract_reads(
      kreport_path = kreport_path,
      kraken_input_seqs = kraken_input_seqs_fasta,
      outfile_prefix = sequence_outfile_prefix3,
      compress_output_seqfile = TRUE,
      taxid = 561,
      include_children = TRUE,
      kraken_stdout_path = kraken_stdout_path
    ),
    "input sequence type: fasta"
  ) |> suppressMessages()

  # Detects compression
  expect_message(
    kraken_extract_reads(
      kreport_path = kreport_path,
      kraken_input_seqs = kraken_input_seqs_fasta,
      outfile = sequence_outfile_prefix,
      compress_output_seqfile = FALSE,
      taxid = 561,
      include_children = TRUE,
      kraken_stdout_path = kraken_stdout_path
    ),
    "input_seqs_compressed: FALSE"
  ) |> suppressMessages()


  #Throw Error if no reads were classified as taxid
  expect_error(
    kraken_extract_reads(
      kreport_path = kreport_path,
      kraken_input_seqs = kraken_input_seqs_fasta,
      outfile = sequence_outfile_prefix2,
      taxid = 561,
      include_children = FALSE,
      kraken_stdout_path = kraken_stdout_path
    ),
    "No reads were classified as taxid/s: 561"
  ) |> suppressMessages()

  # Unpaired Compressed Fastq Input --------------------------------------------------
  # Runs without error
  expect_error(
    kraken_extract_reads(
      kreport_path = kreport_path,
      kraken_input_seqs = kraken_input_seqs_fastq_gzip,
      outfile = sequence_outfile_prefix3,
      taxid = 561,
      include_children = TRUE,
      kraken_stdout_path = kraken_stdout_path,
      compress_output_seqfile = FALSE
    ),
    NA
  ) |> suppressMessages()

  # Produces Expected Files
  expected_outseq_path_fastq = paste0(sequence_outfile_prefix3, ".reads_classified_as_taxid561_or_children.unpaired.1.fastq")
  expected_readnames_path_fastq = paste0(sequence_outfile_prefix3, ".reads_classified_as_taxid561_or_children.txt")
  expect_true(file.exists(expected_outseq_path_fastq))
  expect_true(file.exists(expected_readnames_path_fastq))

  # Check outseq file contents
  outseq_file_contents_fastq = readLines(expected_outseq_path_fastq)
  expected_readnames_contents_fastq = readLines(expected_readnames_path_fastq)

  # --> Is Fastq
  expect_true(grepl(outseq_file_contents_fastq[1], pattern = "^@"))

  # --> Expected Number of Reads in fasta
  expect_true(length(outseq_file_contents_fastq)/4 == 23)

  # --> Expected Number of Read names classified as taxid
  expect_length(expected_readnames_contents_fastq, n = 23)

  # --> Expected contents of read names classified as taxid
  expect_equal(expected_readnames_contents_fastq, c("coli_NC_002127.1-347", "coli_NC_002695.2-213", "coli_NC_002695.2-375",
                                                    "coli_NC_002695.2-860", "coli_NC_002128.1-69", "coli_NC_002127.1-931",
                                                    "coli_NC_002128.1-54", "coli_NC_002127.1-625", "coli_NC_002128.1-782",
                                                    "coli_NC_002127.1-954", "coli_NC_002128.1-930", "coli_NC_002695.2-111",
                                                    "coli_NC_002128.1-210", "coli_NC_002695.2-673", "coli_NC_002128.1-585",
                                                    "coli_NC_002695.2-1", "coli_NC_002695.2-896", "coli_NC_002695.2-174",
                                                    "coli_NC_002695.2-215", "coli_NC_002127.1-166", "coli_NC_002128.1-64",
                                                    "coli_NC_002695.2-796", "coli_NC_002127.1-619"))

  # Detects filetype
  expect_message(
    kraken_extract_reads(
      kreport_path = kreport_path,
      kraken_input_seqs = kraken_input_seqs_fastq_gzip,
      outfile = sequence_outfile_prefix3,
      taxid = 561,
      include_children = TRUE,
      kraken_stdout_path = kraken_stdout_path
    ),
    "input sequence type: fastq"
  ) |> suppressMessages()

  # Detects Compression
  expect_message(
    kraken_extract_reads(
      kreport_path = kreport_path,
      kraken_input_seqs = kraken_input_seqs_fastq_gzip,
      outfile = sequence_outfile_prefix3,
      taxid = 561,
      include_children = TRUE,
      kraken_stdout_path = kraken_stdout_path,
      compress_output_seqfile = TRUE
    ),
    "input_seqs_compressed: TRUE"
  ) |> suppressMessages()

  # Detects Inputs as Unpaired
  expect_message(
    kraken_extract_reads(
      kreport_path = kreport_path,
      kraken_input_seqs = kraken_input_seqs_fastq_gzip,
      outfile = sequence_outfile_prefix3,
      taxid = 561,
      include_children = TRUE,
      kraken_stdout_path = kraken_stdout_path,
      compress_output_seqfile = TRUE
    ),
    "paired inputs: FALSE"
  ) |> suppressMessages()

  # Threads
  expect_message(
    kraken_extract_reads(
      kreport_path = kreport_path,
      kraken_input_seqs = kraken_input_seqs_fastq_gzip,
      outfile = sequence_outfile_prefix3,
      taxid = 561,
      include_children = TRUE,
      kraken_stdout_path = kraken_stdout_path,
      compress_output_seqfile = TRUE, threads = 2
    ),
    "threads: 2"
  ) |> suppressMessages()


  # Paired Inputs -----------------------------------------------------------
  # Warns user in inputing 2 seqs to kraken_input_seqs to use kraken_input_seqs2 if their data is paired
  expect_warning(
    kraken_extract_reads(
      kreport_path = kreport_path,
      kraken_input_seqs = c(kraken_input_seqs_fastq_gzip, kraken_input_seqs_fastq_gzip),
      outfile = sequence_outfile_prefix4,
      taxid = 561,
      include_children = TRUE,
      kraken_stdout_path = kraken_stdout_path,
      compress_output_seqfile = TRUE,
      threads = 1
    ),
    "Supplied 2 sequences but have not indicated the data is paired"
  ) |> suppressMessages()


  # Expect error if different numbers of forward and reverse end seqs are provided
  expect_error(
    kraken_extract_reads(
      kreport_path = kreport_path,
      kraken_input_seqs = c(kraken_input_seqs_fastq_gzip, kraken_input_seqs_fastq_gzip),
      kraken_input_seqs2 = kraken_input_seqs_fastq_gzip,
      outfile = sequence_outfile_prefix5,
      taxid = 561,
      include_children = TRUE,
      kraken_stdout_path = kraken_stdout_path,
      compress_output_seqfile = TRUE,
      threads = 1
    ),
    "length of kraken_input_seqs [2] should match kraken_input_seqs2 [1]", fixed=TRUE
  ) |> suppressMessages()

  # Expect error if more than one forward and reverse seqfiles are provided
  expect_error(
    kraken_extract_reads(
      kreport_path = kreport_path,
      kraken_input_seqs = c(kraken_input_seqs_fastq_gzip, kraken_input_seqs_fastq_gzip),
      kraken_input_seqs2 = c(kraken_input_seqs_fastq_gzip, kraken_input_seqs_fastq_gzip),
      outfile = sequence_outfile_prefix5,
      taxid = 561,
      include_children = TRUE,
      kraken_stdout_path = kraken_stdout_path,
      compress_output_seqfile = TRUE,
      threads = 1
    ),
    "Too many kraken_input_sequences", fixed=TRUE
  ) |> suppressMessages()

  expect_error(
    kraken_extract_reads(
      kreport_path = kreport_path,
      kraken_input_seqs = kraken_input_seqs_fastq_gzip,
      kraken_input_seqs2 = kraken_input_seqs_fastq_gzip,
      outfile = sequence_outfile_prefix6,
      taxid = 561,
      include_children = TRUE,
      kraken_stdout_path = kraken_stdout_path,
      compress_output_seqfile = FALSE,
      threads = 1
    ),
    NA
  ) |> suppressMessages()

  # Produces Expected Files
  expected_outseq_path_fastq_forward = paste0(sequence_outfile_prefix6, ".reads_classified_as_taxid561_or_children.forward.1.fastq")
  expected_outseq_path_fastq_reverse = paste0(sequence_outfile_prefix6, ".reads_classified_as_taxid561_or_children.reverse.2.fastq")
  expected_readnames_path_fastq_paired = paste0(sequence_outfile_prefix6, ".reads_classified_as_taxid561_or_children.txt")
  expect_true(file.exists(expected_outseq_path_fastq_forward))
  expect_true(file.exists(expected_outseq_path_fastq_reverse))
  expect_true(file.exists(expected_readnames_path_fastq_paired))

  # Check outseq file contents
  outseq_file_contents_fastq_forward = readLines(expected_outseq_path_fastq_forward)
  outseq_file_contents_fastq_reverse = readLines(expected_outseq_path_fastq_reverse)
  expected_readnames_contents_fastq_paired = readLines(expected_readnames_path_fastq_paired)

  # --> Is Fastq
  expect_true(grepl(outseq_file_contents_fastq_forward[1], pattern = "^@"))
  expect_true(grepl(outseq_file_contents_fastq_reverse[1], pattern = "^@"))

  # --> Expected Number of Reads in fasta
  expect_true(length(outseq_file_contents_fastq_forward)/4 == 23)
  expect_true(length(outseq_file_contents_fastq_reverse)/4 == 23)

  # --> Expected Number of Read names classified as taxid
  expect_length(expected_readnames_contents_fastq_paired, n = 23)

  # --> Expected contents of read names classified as taxid
  expect_equal(expected_readnames_contents_fastq_paired, c("coli_NC_002127.1-347", "coli_NC_002695.2-213", "coli_NC_002695.2-375",
                                                    "coli_NC_002695.2-860", "coli_NC_002128.1-69", "coli_NC_002127.1-931",
                                                    "coli_NC_002128.1-54", "coli_NC_002127.1-625", "coli_NC_002128.1-782",
                                                    "coli_NC_002127.1-954", "coli_NC_002128.1-930", "coli_NC_002695.2-111",
                                                    "coli_NC_002128.1-210", "coli_NC_002695.2-673", "coli_NC_002128.1-585",
                                                    "coli_NC_002695.2-1", "coli_NC_002695.2-896", "coli_NC_002695.2-174",
                                                    "coli_NC_002695.2-215", "coli_NC_002127.1-166", "coli_NC_002128.1-64",
                                                    "coli_NC_002695.2-796", "coli_NC_002127.1-619"))

})
