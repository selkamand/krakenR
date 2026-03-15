

#' Extract total primary alignments from flagstats file
#'
#' Since most aligners keep unmapped reads and label them as primary -
#' this value is usually a good proxy
#' for the total number of reads fed into the aligner.
#'
#' @param flagstats path to a 'flagstats' file created by \code{samtools flagstat <bam>}
#'
#' @returns total number of primary alignments (whether or not they pass or fail QC)
#' @export
#'
flagstats_to_total_primary_reads <- function(flagstats){
  flagstat_line2 <- readLines(flagstats, n=2)[2]
  csv_string <- sub(x = flagstat_line2, pattern = "([0-9]+) \\+ ([0-9]+) .*", replacement = "\\1,\\2")
  total_primary <- utils::read.csv(text = csv_string, header = FALSE, col.names = c("qc_passed", "qc_failed"))
  total_primary$total <- total_primary$qc_passed + total_primary$qc_failed
  return(total_primary$total)
}
