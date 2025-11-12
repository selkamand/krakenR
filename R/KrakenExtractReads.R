#' Extract Reads
#'
#' Extract reads from a fastq/fasta that kraken has classified as belonging to some species + optionally children.
#'
#' @param kreport_path path to kraken report (string)
#' @param taxid Taxonomy ID of reads to extract from input seqs (numeric)
#' @param include_children include children of supplied taxid (boolean)
#' @param kraken_stdout_path path to stdout (string)
#' @param kraken_input_seqs path to fasta/fastq (optionally compressed) files that were the inputs for the kraken. If input was paired-end sequences, include only paths to forward reads here. Supply any reverse reads to kraken_input_seqs2 (character)
#' @param kraken_input_seqs2  path to reverse reads. Leave null if input was unpaired (character)
#' @param outfile_prefix outfile prefix (string)
#' @param threads number of threads to use (string)
#' @param compress_output_seqfile should output sequence files be gzipped? (flag)
#'
#' @return invisibly returns paths to fastq/fasta sequences classifed as \strong{taxid} and (optionally) child taxids if \strong{include_children = TRUE}. (character)
#' @export
#'
#' @examples
#' \dontrun{
#' kreport_path = system.file(package="krakenR", "testfiles/DRR002014_1.100sequences.minikraken8GB.kreport")
#' kraken_stdout_path = system.file(package="krakenR", "testfiles/DRR002014_1.100sequences.minikraken8GB.stdout.tsv")
#' kraken_input_seqs = c(
#'   system.file(package="krakenR", "testfiles/DRR002014_1.100sequences.fastq.gz"),
#'   system.file(package="krakenR", "testfiles/DRR002014_2.100sequences.fastq.gz")
#' )
#' sequence_outfile_prefix = "mykrakenrun"
#'
#' kraken_extract_reads(
#'   kreport_path = kreport_path,
#'   kraken_input_seqs = kraken_input_seqs,
#'   outfile_prefix = sequence_outfile_prefix,
#'   taxid = 562,
#'   include_children = TRUE,
#'   kraken_stdout_path = kraken_stdout_path,
#'   threads = 1
#' )
#' }
#'
kraken_extract_reads <- function(kreport_path, kraken_stdout_path, kraken_input_seqs, kraken_input_seqs2 = NULL, outfile_prefix, taxid, include_children, threads = 1, compress_output_seqfile = TRUE){

  # Validating Inputs --------------------------------------------------------------
  paired = !is.null(kraken_input_seqs2)


  ## Check Inputs: Assertions --------------------------------------------------------------
  cli::cli_h1("Checking Inputs")
  #assertthat::assert_that(assertthat::is.number(taxid))
  if(!assertthat::is.number(taxid)) stop(glue::glue("Expected taxid [{taxid}] to be a number, not a {class(taxid)}"))

  if(!assertthat::is.string(kreport_path)) stop(glue::glue("Expected kreport_path [{kreport_path}] to be a string, not a {class(kreport_path)}"))

  #assertthat::assert_that(file.exists(kreport_path))
  if(!file.exists(kreport_path)) stop(glue::glue("kreport_path: File doesn't exist [{kreport_path}]"))

  #assertthat::assert_that(assertthat::is.string(kraken_stdout_path))
  if(!assertthat::is.string(kraken_stdout_path)) stop(glue::glue("Expected kraken_stdout_path to be a string, not a {class(kraken_stdout_path)}"))

  #assertthat::assert_that(file.exists(kraken_stdout_path))
  if(!file.exists(kraken_stdout_path)) stop(glue::glue("kraken_stdout_path: File doesn't exist [{kraken_stdout_path}]"))

  #assertthat::assert_that(is.character(kraken_input_seqs))
  if(!all(is.character(kraken_input_seqs))) stop(glue::glue("Expected kraken_input_seqs to be a character vector, not a {class(kraken_input_seqs[!is.character(kraken_input_seqs)])}"))

  #assertthat::assert_that(all(file.exists(kraken_input_seqs)))
  if(!all(file.exists(kraken_input_seqs))) stop(glue::glue("kraken_input_seqs: File doesn't exist [{kraken_input_seqs[!file.exists(kraken_input_seqs)]}]"))

  if(!is.null(kraken_input_seqs2)){
    if(!all(is.character(kraken_input_seqs2))) stop(glue::glue("Expected kraken_input_seqs2 to be a character vector, not a {class(kraken_input_seqs2[!is.character(kraken_input_seqs2)])}"))
    if(!all(file.exists(kraken_input_seqs2))) stop(glue::glue("kraken_input_seqs2: File doesn't exist [{kraken_input_seqs2[!file.exists(kraken_input_seqs2)]}]"))
  }

  #assertthat::assert_that(assertthat::is.flag(include_children))
  if(!assertthat::is.flag(include_children)) stop(glue::glue("Expected include_children to be a character vector, not a {class(include_children)}"))

  #assertthat::assert_that(assertthat::is.string(outfile_prefix))
  if(!assertthat::is.string(outfile_prefix)) stop(glue::glue("Expected outfile_prefix to be a string, not a {class(outfile_prefix)}"))

  # Pairing Assertions
  if (paired){
    if(length(kraken_input_seqs) != length(kraken_input_seqs2)){
      #stop(glue::glue("length of kraken_input_seqs [{length(kraken_input_seqs)}] should match kraken_input_seqs2 [{length(kraken_input_seqs2)}] for paired data. If your data is unpaired, supply path to all fastqs/fasta files to kraken_input_seqs as a vector"))
      stop(glue::glue("length of kraken_input_seqs [{length(kraken_input_seqs)}] should match kraken_input_seqs2 [{length(kraken_input_seqs2)}] for paired data. If your data is unpaired, supply path to all fastqs/fasta files to kraken_input_seqs as a vector"))
    }

    if(length(kraken_input_seqs) != 1)
      stop(glue::glue("Too many kraken_input_sequences. Expected 1, got [{length(kraken_input_seqs)}]. If input is paired please supply only 1 input seq file to kraken_input_seqs, and 1 to kraken_input_seqs2"))

    names(kraken_input_seqs) <- rep("forward", times=length(kraken_input_seqs))
    names(kraken_input_seqs2) <- rep("reverse", times=length(kraken_input_seqs2))
  }
  else{
    if(length(kraken_input_seqs) == 2){
      cli::cli_warn(
        "Supplied 2 sequences but have not indicated the data is paired.
       Treating data as unpaired. If data is paired please rerun with {.code kraken_input_seqs = 'path/to/forward_reads.fa' AND kraken_input_seqs2 = 'path/to/reverse_reads.fa'}"
      )
    }
    names(kraken_input_seqs) <- rep("unpaired", times=length(kraken_input_seqs))
  }

  input_seqs <- c(kraken_input_seqs, kraken_input_seqs2)


  # Ensure extension is appropriate
  input_seqs_compressed = grepl(input_seqs, pattern = ".gz$")
  input_seq_class = filepath_get_sequence_type(input_seqs)

  extension = paste0(input_seq_class, ifelse(input_seqs_compressed, yes = ".gz", no=""))
  assertthat::assert_that(all(input_seq_class != "other"), msg = "Input Sequences had")



  ## Print Config Info -------------------------------------------------------
  cli::cli_h2("File Inputs")
  cli::cli_bullets("{.strong kraken report}")
  cli::cli_bullets(c(v = "{kreport_path}", ""))


  cli::cli_bullets("{.strong kraken stdout}")
  cli::cli_bullets(c(v = "{kraken_stdout_path}", ""))

  cli::cli_bullets("{.strong input sequences}")
  cli::cli_bullets(c(v = "{paste0(input_seqs, ' (',names(input_seqs), ')')}", ""))
  cli::cli_bullets(c(">" = "input_seqs_compressed: {input_seqs_compressed}",
                     ">" = "input sequence type: {input_seq_class}"
  ))


  cli::cli_h2("Other")
  #cli::cli_rule("{.strong Other}")
  cli::cli_bullets(c(
    v ="{.strong taxid}: {taxid}",
    v= "{.strong include children?}: {include_children}",
    v= "{.strong paired inputs}: {paired}",
    v= "{.strong threads}: {threads}"
  ))



  ## Expected Output Files ---------------------------------------------------
  cli::cli_h2("Expected Outputs")

  # Reads classified as taxid
  outfile_readnames = paste0(outfile_prefix, ".reads_classified_as_taxid", taxid, ifelse(include_children, yes = "_or_children", no = ""), ".txt")

  # Fastqs/Fasta output files
  outfile_sequences = paste0(outfile_prefix, ".reads_classified_as_taxid", taxid, ifelse(include_children, yes = "_or_children", no = ""), ".", names(input_seqs), ".", seq_along(input_seqs),".", input_seq_class, ifelse(compress_output_seqfile, yes = ".gz", no = ""))

  if(anyDuplicated(outfile_sequences)){
    duplicated_outfile_sequences = outfile_sequences[duplicated(outfile_sequences)]
    stop(glue::glue(
      "Cannot proceed since outfiles containing classified reads includes duplicates. \n{duplicated_outfile_sequences} are duplicated."))
  }

  # CLI messages()
  cli::cli_bullets("{.strong readnames classified as taxid}")
  cli::cli_bullets(c(">"="{outfile_readnames}", ""))

  cli::cli_bullets("{.strong fasta/fastqs classified as taxid}")
  cli::cli_bullets(c(">"="{outfile_sequences}", ""))


  ## Check Dependencies ------------------------------------------------------
  cli::cli_h1("Checking Dependencies")
  assert_dependencies_present(c("awk", "parallel", "seqkit"))

  # Extract Reads by Taxid --------------------------------------------------
  cli::cli_h1("Extract Reads by Taxid")
  cli::cli_alert_info("Extracting taxid {taxid} {ifelse(include_children, 'and children ', '')}from {length(input_seqs)} input file{?s}:")

  # Get taxids
  if(include_children) {
    cli::cli_h2("Getting child taxids")
    cli::cli_alert_info("Extracting child taxids from kreport: {basename(kreport_path)}")

    taxids = kraken_report_get_child_taxids(path_to_kraken_report = kreport_path, taxid = taxid, inclusive = TRUE, exclude_children_without_reads_directly_classified = TRUE)

    # Will only return child taxids which have at least 1 read directly assigned to them according to kraken report
    cli::cli_alert_success("Found {length(taxids)} taxid{?s} that represent {taxid} + any children with at least 1 read directly assigned to them according to kraken")
    cli::cli_alert(c("Taxids: ", paste0(taxids, collapse = ", ")))
  }
  else
  {
    taxids = taxid
  }


  # Identify Readnames ------------------------------------------------------
  ## Build Awk Command -------------------------------------------------------
  cli::cli_h1("Build Awk command")
  cli::cli_alert_info("We use a parallel awk to create make a file of readnames belonging to our taxids of interest{ifelse(include_children, yes=' and its children', no='')}")

  awk_begin = paste0(
    "BEGIN{",
    paste0("arr[",taxids, "]=0", collapse = "; "),
    "}"
  )
  awk_main = paste0(
    "{ if(\\$3 in arr) print \\$2 }"
  )
  awk_command = glue::glue("awk '{awk_begin} {awk_main}'")

  #cli::cli_alert_info("Awk Command: {awk_command}")

  parallel_command = glue::glue("parallel --pipepart -j {threads} -a {kraken_stdout_path}")

  parallel_awk_command = glue::glue('{parallel_command} "{awk_command}" > {outfile_readnames}') # we don't need to do a sort and uniq to deduplicate readnames (duplicates from same readname in forward and reverse reads) since seqkit handles this well
  cli::cli_bullets("")
  cli::cli_alert_info("Parallel Awk Command:\n {.code {parallel_awk_command}}")


  ## Run Parallel Awk ----------------------------------------------------------------
  cli::cli_h2("Running Parallel Awk (this may take a while) ...")
  parallel_awk_exit_code <- system(parallel_awk_command)

  if(parallel_awk_exit_code != 0 || !file.exists(outfile_readnames)){
    stop(glue::glue("Failed to produce the file listing read names of the taxid of interest {ifelse(include_children, ' and its children', '')}"))
  }
  else{
    cli::cli_alert_success("Names of reads classified as taxids of interest children saved to: {outfile_readnames}")
  }

  # Make sure file containin readnames isn't empty
  if(length(readLines(outfile_readnames, n = 50)) == 0) {
    stop(glue::glue("Terminating early: No reads were classified as taxid/s: {taxids}"))
  }

  # Extract Sequences -------------------------------------------------------
  cli::cli_h1("Extract Reads")
  seqkit_commands <- character(0)
  cli::cli_h2("Building Seqkit Command:")
  for (i in seq_along(input_seqs)){
    input_seq = input_seqs[[i]]
    outfile_sequence = outfile_sequences[[i]]
    seqkit_command = glue::glue("seqkit grep --threads {threads} --pattern-file {outfile_readnames} {input_seq} -o {outfile_sequence}")
    cli::cli_alert("`[{i}] {seqkit_command}`")
    cli::cli_bullets("")
    seqkit_commands <- c(seqkit_commands, seqkit_command)
  }

  cli::cli_h2("Extracting Reads: ")
  seqkit_exit_codes = numeric(0)
  for (i in seq_along(seqkit_commands)){
    cli::cli_inform("running seqkit command [{i}]")
    seqkit_exit_code = system(seqkit_commands[i])
    if (seqkit_exit_code == 0)  cli::cli_alert_success("success")
    else stop(glue::glue("Seqkit command threw exit code [{seqkit_exit_code}]"))

    seqkit_exit_codes = c(seqkit_exit_codes, seqkit_exit_code)
  }

  cli::cli_bullets("")
  if(all(seqkit_exit_code == 0) & all(file.exists(outfile_sequences))) {
    cli::cli_alert_success("Successully extracted reads to {outfile_sequences}")
  }


  return(invisible(outfile_sequences))
}


dependencies_present <- function(program_names){
  program_paths <- Sys.which(program_names)
  return(!(program_paths==""))
}

assert_dependencies_present <- function(program_names){
  #program_names <- c("awk", "parallel", "seqkit")

  names(program_names) <- ifelse(dependencies_present(program_names), no = "x", yes = "v")
  programs_not_found <- program_names[!dependencies_present(program_names)]
  programs_not_found_string <- paste0(names(programs_not_found), collapse = ", ")

  cli::cli_bullets(c("Checking Dependencies: ", program_names))
  if(!all(dependencies_present(program_names))){
    stop(glue::glue("Missing Dependencies: ", paste0(programs_not_found_string, collapse = ",")))
  }
}
