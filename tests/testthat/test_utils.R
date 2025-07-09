test_that(".is_spectra_installed works", {
    expect_true(.is_spectra_installed())
})

test_that(".is_ms_experiment_installed works", {
    expect_true(.is_ms_experiment_installed())
})

test_that(".is_ms_backend_metabo_lights_installed works", {
    expect_true(.is_ms_backend_metabo_lights_installed())
})

test_that(".is_xcms_installed works", {
    expect_true(.is_xcms_installed())
})
