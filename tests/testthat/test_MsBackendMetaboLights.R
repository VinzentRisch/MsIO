## Tests for the MsBackendMetaboLights backend.

library(Spectra)
library(MsBackendMetaboLights)

be_mtbls <- backendInitialize(MsBackendMetaboLights(), mtblsId = "MTBLS39",
                              filePattern = "63A.cdf", offline = FALSE)

test_that("readMsObject,PlainTextParam,MsBackendMetaboLights works", {
    pth <- file.path(tempdir(), "test_backend_ml")
    saveMsObject(be_mtbls, PlainTextParam(pth))

    ## read
    res <- readMsObject(MsBackendMetaboLights(), PlainTextParam(pth),
                        offline = TRUE)
    expect_true(validObject(res))
    expect_equal(rtime(be_mtbls), rtime(res))
    expect_equal(mz(be_mtbls), mz(res))

    ## Clean cache first to check errors etc.
    bfc <- BiocFileCache::BiocFileCache()
    BiocFileCache::cleanbfc(bfc, days = -10, ask = FALSE)

    ## read again offline throws error
    expect_error(readMsObject(MsBackendMetaboLights(), PlainTextParam(pth),
                              offline = TRUE), "No locally cached")

    ## read re-downloading the data
    res <- readMsObject(MsBackendMetaboLights(), PlainTextParam(pth),
                        offline = FALSE)
    expect_true(validObject(res))
    expect_equal(rtime(be_mtbls), rtime(res))
    expect_equal(mz(be_mtbls), mz(res))

    unlink(pth, recursive = TRUE)

    pth <- file.path(tempdir(), "test_backend_ml")
    dir.create(pth)
    expect_error(readMsObject(MsBackendMetaboLights(), PlainTextParam(pth)),
                 "found in the provided path.")
    writeLines("# Some line\nnext line\nthird line\n",
               con = file.path(pth, "ms_backend_data.txt"))
    expect_error(readMsObject(MsBackendMetaboLights(), PlainTextParam(pth)),
                 "Invalid class in")
})
