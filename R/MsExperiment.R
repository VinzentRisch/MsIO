#' @include PlainTextParam.R
#' @title Methods to save and load contents of a MsExperiment object
#'
#' @author Philippine Louail
#'
#' @importFrom ProtGenerics spectra
#'
#' @importFrom S4Vectors DataFrame SimpleList
#'
#' @noRd
NULL

#' @rdname PlainTextParam
setMethod("saveMsObject",
          signature(object = "MsExperiment",
                    param = "PlainTextParam"),
          function(object, param){
              dir.create(path = param@path,
                         recursive = TRUE,
                         showWarnings = FALSE)
              ## sample data
              sdata <- object@sampleData
              if (!length(sdata)) # initialize with empty data frame
                  sdata <- DataFrame(sample_name = character())
              write.table(as.data.frame(sdata), sep = "\t",
                          file = file.path(param@path,
                                           "ms_experiment_sample_data.txt"))

              ## sample data links
              sdl <- object@sampleDataLinks
              if (length(sdl) > 0) {
                  lapply(names(sdl), function(x){
                      fl <- file.path(
                          param@path,
                          paste0("ms_experiment_sample_data_links_", x, ".txt"))
                      write.table(sdl[[x]], file = fl, row.names = FALSE,
                                  col.names = FALSE, sep = "\t")
                  })
                  write.table(
                      sdl@elementMetadata, sep = "\t", quote = TRUE,
                      file = file.path(param@path,
                                       "ms_experiment_link_mcols.txt"))
              }
              ## call export of individual other objects (not MsExperiment data)
              if (length(spectra(object)))
                  saveMsObject(spectra(object), param)
              ## at some point also chromatograms, etc.
          }
)

#' @rdname PlainTextParam
#'
#' @importMethodsFrom S4Vectors mcols<-
setMethod("readMsObject",
          signature(object = "MsExperiment",
                    param = "PlainTextParam"),
          function(object, param, ...) {
              ## read sample data
              fl <- file.path(param@path, "ms_experiment_sample_data.txt")
              if (!file.exists(fl))
                  stop("No 'ms_experiment_sample_data.txt' file found in ",
                       "the provided path.")
              sd <- read.table(fl, sep = "\t", header = TRUE)
              object@sampleData <- DataFrame(sd, row.names = NULL)

              ## read spectra
              if (file.exists(file.path(param@path, "spectra_slots.txt")))
                  object@spectra <- readMsObject(Spectra::Spectra(), param, ...)
              ## sample data links
              fl <- list.files(
                  param@path,
                  pattern = "ms_experiment_sample_data_links_.*\\.txt",
                  full.names = TRUE)
              if (length(fl) > 0) {
                  n <- gsub("ms_experiment_sample_data_links_|\\.txt", "",
                            basename(fl))
                  sdl <- lapply(fl, function(x) {
                      unname(as.matrix(read.table(x, sep = "\t")))
                  })
                  names(sdl) <- n
                  object@sampleDataLinks <- SimpleList(sdl)
                  em <- read.table(file.path(param@path,
                                             "ms_experiment_link_mcols.txt"),
                                   sep = "\t", header = TRUE)
                  mcols(object@sampleDataLinks) <- DataFrame(
                      em, row.names = NULL)
              }

              validObject(object)
              object
          })

################################################################################
##
## MetaboLights readMsObject
##
################################################################################
#' @rdname MetaboLightsParam
#' @importFrom utils menu
setMethod("readMsObject",
          signature(object = "MsExperiment",
                    param = "MetaboLightsParam"),
          function(object, param, keepOntology = TRUE, keepProtocol = TRUE,
                   simplify = TRUE, ...) {
              if (!.is_ms_backend_metabo_lights_installed())
                  stop("Required package 'MsBackendMetaboLights' is missing. ",
                       "Please install it and try again.", call. = FALSE)
              if (!.is_spectra_installed())
                  stop("Required package 'Spectra' is missing. ",
                       "Please install and try again.", call. = FALSE)
              pth <- MsBackendMetaboLights::mtbls_ftp_path(param@mtblsId)
              all_fls <- MsBackendMetaboLights::mtbls_list_files(param@mtblsId)

              ## Extract and read assay files
              assays <- all_fls[grepl("^a_", all_fls)]
              if (length(param@assayName) > 0) {
                  selected_assay <- param@assayName
                  if (!selected_assay %in% assays)
                      stop("Specified assay \"", selected_assay, "\" does ",
                           "not exist.", call. = FALSE)
              }
              else {
                  if (length(assays) == 1) {
                      selected_assay <- assays
                      message("Only one assay file found: ", selected_assay)
                  } else {
                      message("Multiple assay files found:\n")
                      selection <- menu(assays,
                                        title = paste("Please choose the assay",
                                                      "file you want to use:"))
                      selected_assay <- assays[selection]
                  }
              }

              assay_data <- read.table(paste0(pth, selected_assay),
                                       header = TRUE, sep = "\t",
                                       check.names = FALSE)

              ## Extract and read sample info files
              s_files <- all_fls[grepl("^s_", all_fls)]
              sample_info <- read.table(paste0(pth, s_files),
                                        header = TRUE, sep = "\t",
                                        check.names = FALSE)

              ## merging
              ord <- match(assay_data$`Sample Name`, sample_info$`Sample Name`)
              merged_data <- cbind(assay_data, sample_info[ord, ])
              if (keepProtocol || keepOntology || simplify)
                  merged_data <- MsIO:::.clean_merged(x = merged_data,
                                               keepProtocol = keepProtocol,
                                               keepOntology = keepOntology,
                                               simplify = simplify)
              ## Assemble object
              b <- MsBackendMetaboLights::MsBackendMetaboLights()
              object@spectra <- Spectra::Spectra(
                                             mtblsId = param@mtblsId,
                                             source = b,
                                             assayName = selected_assay,
                                             filePattern = param@filePattern)

              ## sample to spectra link
              fl <- object@spectra@backend@spectraData$derived_spectral_data_file[1L]
              ## identify the column in the assay/sample description containing
              ## the file name information.
              nms <- c("Raw Spectral Data File", "Derived Spectral Data File")
              nms <- nms[nms %in% colnames(merged_data)]
              nme <- nms[vapply(nms, function(z)
                  any(merged_data[, z] %in% fl), NA)]
              merged_data <- merged_data[grepl(param@filePattern,
                                               merged_data[, nme]), ,
                                         drop = FALSE]
              nnme <- gsub(" ", "_", nme, fixed = TRUE)
              colnames(merged_data)[colnames(merged_data) == nme] <- nnme
              object@sampleData <- DataFrame(merged_data, check.names = FALSE,
                                             row.names = NULL)
              w <- paste0("sampleData.", nnme,
                          "= spectra.derived_spectral_data_file")
              object <- MsExperiment::linkSampleData(object, with = w)
              validObject(object)
              object
          })


#####HELPERS

#' Function that takes the extra parameters and clean the metadata if asked by
#' the user.
#'
#' @noRd
.clean_merged <- function(x, keepProtocol, keepOntology, simplify) {
    # remove ontology
    if (!keepOntology)
        x <- x[, -which(grepl("Term", names(x))), drop = FALSE]

    # remove protocol
    if (!keepProtocol)
        x <- x[, -which(grepl("Protocol|Parameter", names(x))),  drop = FALSE]

    # remove duplicated columns contents and NAs
    if (simplify) {
        x <- x[, !duplicated(as.list(x)), drop = FALSE]
        x <- x[, colSums(is.na(x)) != nrow(x), drop = FALSE]
    }
    x
}
