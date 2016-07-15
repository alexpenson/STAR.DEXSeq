#' Test taken from DEXSeq documentation
#'
#' @return NULL
#' @export test_from_the_docs
#'
## #' @examples
#' @import DEXSeq
#' @importFrom S4Vectors DataFrame
test_from_the_docs <- function(){

  S4Vectors::DataFrame()
  if(packageVersion("DEXSeq") < "1.16.10"){
    warning("may need to update DEXSeq")
  }

  #######################################
  ### From elementary data structures ###
  #######################################
  countData <- matrix( rpois(10000, 100), nrow=1000 )
  sampleData <- data.frame(
    condition=rep( c("untreated", "treated"), each=5 ) )
  design <- formula( ~ sample + exon + condition:exon )
  groupID <- rep(
    paste0("gene", 1:10),
    each= 100 )
  featureID <- rep(
    paste0("exon", 1:100),
    times= 10 )
  dds <- DEXSeq::DEXSeqDataSet( countData, sampleData, design,
                                featureID, groupID )
  dds <- DEXSeq::DEXSeq(dds)
  dds
}


getSampleData <- function(sampleTable,
                          read_depth_threshold = 10 ### remove junctions covered by less than this amount summed over all samples
){
  ### read splice junction data
  SJ <- sampleTable[, fread(fileName), list(sampleName)]
  setnames(SJ, c("sampleName", "chr", "start", "end", "strand", "motif",
                 "annotated", "Nunique", "Nmulti", "maxOverhang"))
  ### create table and filter
  SJc <- dcast.data.table(SJ[, sum_Nunique := sum(Nunique),
                             list(chr, start, end, strand, motif, annotated) ### junction identifier
                             ][sum_Nunique >= read_depth_threshold],
                          chr + start + end + strand + motif + annotated + sum_Nunique ~ sampleName,
                          value.var = "Nunique", fill = 0)
  SJc
}


### extract splice junctions from gtf
extract_junctions <- function(
  gtf_filename = 'gencode.v19.annotation.gtf',
  gene_name_regexp = '(?<=gene_name \")[A-Z0-9-.]*',
  junctions_output_file = gsub(".gtf$", ".junctions.rds", gtf_filename)
){

  if(file.exists(junctions_output_file)){
    junctions <- readRDS(junctions_output_file)
  } else {
    gtf <- fread(gtf_filename)
    setnames(gtf, c("chr", "source", "feature", "start", "end",
                    "score", "strand", "frame", "attribute"))
    gtf[, gene_name := stringr::str_extract(attribute,
                                            gene_name_regexp)]
    setkeyv(gtf, c("chr", "start", "end"))
    junctions <- gtf[!duplicated(gtf) & feature == "exon"]
    setnames(junctions, c("chr", "source", "feature",
                          "end", "start", ### reversed!
                          "score", "strand", "frame", "attribute",
                          "gene_name"))
    junctions[, start := start + 1]
    junctions[, end := end - 1]

    setkeyv(junctions, c("chr", "start", "end"))
    saveRDS(junctions, junctions_output_file)
  }
}

annotate_known_splice_sites <- function(SJ, junctions){
  setkeyv(junctions, c("chr", "end"))
  setkeyv(SJ, c("chr", "end"))
  SJend <- merge(SJ,
                 junctions[, gene_name, keyby = list(chr, end)],
                 all.x = T,
                 allow.cartesian=TRUE)

  setkeyv(junctions, c("chr", "start"))
  setkeyv(SJend, c("chr", "start"))
  SJnovel <- merge(SJend,
                   junctions[, list(chr, start, gene_name)],
                   all.x = T,
                   suffixes = c("_end", "_start"))

  ### silently remove splice junctions with two novel endpoints
  SJnovel <- SJnovel[!( is.na(gene_name_end) & is.na(gene_name_start) )]

  ### silently remove splice junctions joining two genes
  SJnovel <- SJnovel[!( !is.na(gene_name_end) & !is.na(gene_name_start) & gene_name_end != gene_name_start )]

  ### annotate as known-known or novel splice site at the 3' or 5' end.
  SJnovel[, ann := ifelse(!is.na(gene_name_end) & !is.na(gene_name_start),
                          "kk",
                          ifelse((strand == 1 & is.na(gene_name_end) |
                                    strand == 2 & is.na(gene_name_start)),
                                 "n3",
                                 "n5")
  )
  ]

  ### consolidate gene name
  SJnovel[, gene_name := ifelse(is.na(gene_name_end), gene_name_start, gene_name_end)]
  SJnovel[, gene_name_end := NULL]
  SJnovel[, gene_name_start := NULL]

  ### remove duplicate lines
  setkey(SJnovel, NULL)
  SJnovel <- SJnovel[!duplicated(SJnovel)]

  SJnovel
}

getDEXSeqDataSet <- function(sampleTable, junctions){
  SJ <- getSampleData(sampleTable)
  SJ_ann <- annotate_known_splice_sites(SJ, junctions = junctions)
  dds <- DEXSeqDataSet(countData = as.matrix(SJ_ann[, sampleTable$sampleName, with = F]),
                       sampleData = data.frame(condition = factor(sampleTable$condition)),
                       featureID = SJ_ann[, paste0(chr, "_", start, "-", end)],
                       groupID = SJ_ann$gene_name)
}
