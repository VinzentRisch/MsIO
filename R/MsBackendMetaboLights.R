#' @rdname PlainTextParam
setMethod("readMsObject", signature(object = "MsBackendMetaboLights",
                                    param = "PlainTextParam"),
          function(object, param, offline = FALSE) {
              fl <- file.path(param@path, "ms_backend_data.txt")
              if (!file.exists(fl))
                  stop("No 'ms_backend_data.txt' file found in ",
                       "the provided path.")
              l2 <- readLines(fl, n = 2)
              if (l2[1] != "# MsBackendMetaboLights")
                  stop("Invalid class in 'ms_backend_data.txt' file.")
              if (length(l2) > 1L) {
                  data <- read.table(file = fl, sep = "\t", header = TRUE)
                  rownames(data) <- NULL
                  slot(object, "spectraData", check = FALSE) <- DataFrame(data)
                  MsBackendMetaboLights::mtbls_sync(object, offline = offline)
              }
              validObject(object)
              object
          })
