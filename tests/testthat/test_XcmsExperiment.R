library(xcms)

xmse <- loadXcmsData()
xmseg_filt <- filterMzRange(xmse, c(200, 500))
xmseg_filt <- filterRt(xmseg_filt, c(3000, 4000))


test_that("saveMsObject,readMsObject,PlainTextParam,XcmsExperiment works", {
    pth = file.path(tempdir(), "test_xcmsexp")
    param <- PlainTextParam(path = pth)
    saveMsObject(xmseg_filt, param = param)
    expect_true(dir.exists(pth))
    expect_true(file.exists(
        file.path(param@path, "ms_experiment_sample_data.txt")))
    expect_true(file.exists(
        file.path(param@path, "ms_backend_data.txt")))
    expect_true(file.exists(
        file.path(param@path, "spectra_slots.txt")))
    expect_true(file.exists(
        file.path(param@path, "spectra_processing_queue.json")))
    expect_true(file.exists(
        file.path(param@path, "xcms_experiment_process_history.json")))
    expect_true(file.exists(
        file.path(param@path, "xcms_experiment_chrom_peaks.txt")))
    expect_true(file.exists(
        file.path(param@path, "xcms_experiment_chrom_peak_data.txt")))
    expect_true(file.exists(
        file.path(param@path, "xcms_experiment_feature_definitions.txt")))
    expect_true(file.exists(
        file.path(param@path, "xcms_experiment_feature_peak_index.txt")))

    ## load data again
    ## This error is not thrown with rcmdcheck::rcmdcheck()
    ## expect_error(readMsObject(new("XcmsExperiment"), param), "load the library")
    ## library(Spectra)
    load_xmse <- readMsObject(new("XcmsExperiment"), param)
    expect_true(inherits(load_xmse, "XcmsExperiment"))
    expect_equal(xmseg_filt@featureDefinitions,
                 load_xmse@featureDefinitions)
    expect_equal(featureValues(xmseg_filt), featureValues(load_xmse))
    expect_equal(adjustedRtime(xmseg_filt), adjustedRtime(load_xmse))
    expect_no_error(filterRt(load_xmse, c(3000, 3500)))
    expect_equal(xmseg_filt@chromPeaks, load_xmse@chromPeaks)
    expect_equal(xmseg_filt@chromPeakData, load_xmse@chromPeakData)
    expect_equal(xmseg_filt@sampleData, load_xmse@sampleData)
    expect_equal(length(xmseg_filt@processHistory), length(load_xmse@processHistory))
    expect_equal(xmseg_filt@processHistory[[1L]], load_xmse@processHistory[[1L]])
    expect_equal(xmseg_filt@processHistory[[2L]], load_xmse@processHistory[[2L]])
    expect_equal(xmseg_filt@processHistory[[3L]], load_xmse@processHistory[[3L]])
    expect_equal(xmseg_filt@processHistory[[4L]], load_xmse@processHistory[[4L]])
    expect_equal(xmseg_filt@processHistory[[5L]], load_xmse@processHistory[[5L]])
    ## The 6th param object contains functions for which the comparison fails
    ## because of the name/namespace mentioned. See e.g.
    ## xmseg_filt@processHistory[[6]]@param load_xmse@processHistory[[6]]@param

    ## Check the spectraPath parameter.
    bp <- dataStorageBasePath(xmseg_filt@spectra)
    ## manually change dataStorage path of backend
    sd <- read.table(file.path(param@path, "ms_backend_data.txt"),
                     header = TRUE)
    sd$dataStorage <- sub("faahKO", "other", sd$dataStorage)
    writeLines(
        "# MsBackendMzR", con = file.path(param@path, "ms_backend_data.txt"))
    write.table(sd,
                file = file.path(param@path, "ms_backend_data.txt"),
                sep = "\t", quote = TRUE, append = TRUE)
    expect_error(readMsObject(new("XcmsExperiment"), param), "invalid class")
    expect_no_error(readMsObject(XcmsExperiment(),
                                 param, spectraPath = bp))

    param <- PlainTextParam(tempdir())
    expect_error(readMsObject(XcmsExperiment(), param),
                 "No 'ms_experiment_sample_data")

    ## Export an empty object.
    a <- XcmsExperiment()
    pth = file.path(tempdir(), "test_xcmsexp_empty")
    param <- PlainTextParam(path = pth)
    saveMsObject(a, param)
    res <- readMsObject(XcmsExperiment(), param)
    expect_equal(nrow(chromPeaks(a)), nrow(chromPeaks(res)))
    expect_equal(colnames(chromPeaks(a)), colnames(chromPeaks(res)))
    expect_equal(nrow(chromPeakData(a)), nrow(chromPeakData(res)))
    expect_equal(colnames(chromPeakData(a)), colnames(chromPeakData(res)))
    expect_equal(a@featureDefinitions, res@featureDefinitions)
    expect_equal(a@processHistory, res@processHistory)
})

test_that(".import_chrom_peaks works", {
    pth <- tempdir()
    expect_error(.import_chrom_peaks(xcmse, pth), "chrom_peaks.txt")
    write.table(chromPeaks(xmse),
                file = file.path(pth, "xcms_experiment_chrom_peaks.txt"),
                sep = "\t")
    expect_error(.import_chrom_peaks(xmse, pth),
                 "chrom_peak_data.txt")
    file.remove(file.path(pth, "xcms_experiment_chrom_peaks.txt"))
})

test_that(".import_features works", {
    pth <- tempdir()
    write.table(
        featureDefinitions(xmse)[, 1:8],
        file = file.path(pth, "xcms_experiment_feature_definitions.txt"),
        sep = "\t")
    expect_error(.import_features(xmse, pth), "feature_peak_index.txt")
})

test_that(".import_process_history works", {
    pth <- tempdir()
    expect_error(.import_process_history(xmse, pth), "process_history.json")
})

test_that("saveMsObject,mzTabParam works", {
    faahko <- loadXcmsData("faahko_sub2")
    faahko <- groupChromPeaks(
        faahko, PeakDensityParam(sampleGroups = rep(1, length(faahko))))

    d <- file.path(tempdir(), "mzt_test")
    dir.create(d, recursive = TRUE)

    ## errors
    expect_error(
        saveMsObject(faahko, mzTabParam(studyId = "test_study", path = d,
                                        sampleDataColumn = "sample_name")),
        "has to correspond to column names of the sampleData()")
    expect_error(
        saveMsObject(faahko, mzTabParam(studyId = "test_study", path = d,
                                        sampleDataColumn = "sample_index",
                                        optionalFeatureColumns = "other")),
        "'optionalFeatureColumns' have to correspond")

    p <- mzTabParam(studyId = "test_study", path = d,
                    sampleDataColumn = "sample_index",
                    optionalFeatureColumns = "peakidx")
    saveMsObject(faahko, p)
    expect_true(file.exists(file.path(d, "test_study.mztab")))
    res <- readLines(file.path(d, "test_study.mztab"))
    expect_true(length(res) > 0L)
    expect_true(length(grep("^MTD", res)) > 0)
    expect_true(length(grep("^SML", res)) > 0)
    expect_true(length(grep("^SMF", res)) > 0)
    ## Check for empty lines
    expect_true(length(grep(c("^MTD|SML|SMF"), res, invert = TRUE)) == 2)

    expect_error(
        saveMsObject(faahko, p), "File \"test_study.mztab\" already exists")

    unlink(d, recursive = TRUE)
})
