test_that("an exact 0.50 tie is non-directional Instability", {
  skip("Phase 2A contract: FSF computation is not implemented")

  identity <- data.frame(
    feature_id = "tie",
    p_up = 0.50,
    p_down = 0.50,
    p_const = 0
  )
  observed <- fsf_classify(identity)
  expect_equal(observed$dominant_state, "tied")
  expect_equal(observed$ssi, 0.50)
  expect_equal(as.character(observed$stability_region), "Instability")
  expect_equal(as.character(observed$signal_class), "Instability")
})

test_that("an exact 0.50 unique plurality remains Instability", {
  skip("Phase 2A contract: FSF computation is not implemented")

  identity <- data.frame(
    feature_id = "plurality",
    p_up = 0.50,
    p_down = 0.30,
    p_const = 0.20
  )
  observed <- fsf_classify(identity)
  expect_equal(observed$dominant_state, "up")
  expect_equal(observed$ssi, 0.50)
  expect_equal(as.character(observed$stability_region), "Instability")
  expect_equal(as.character(observed$signal_class), "Instability")
})

test_that("strict majority produces a unique Transitional class", {
  skip("Phase 2A contract: FSF computation is not implemented")

  identity <- data.frame(
    feature_id = "majority",
    p_up = 0.51,
    p_down = 0.49,
    p_const = 0
  )
  observed <- fsf_classify(identity)
  expect_equal(observed$dominant_state, "up")
  expect_equal(as.character(observed$stability_region), "Transitional")
  expect_equal(as.character(observed$signal_class), "Transitional Up")
})


test_that("raw near-half probabilities do not create a fuzzy tie", {
  skip("Phase 2A contract: FSF computation is not implemented")

  identity <- data.frame(
    feature_id = "raw_majority",
    p_up = 0.5000000001,
    p_down = 0.4999999999,
    p_const = 0
  )
  observed <- fsf_classify(identity)
  expect_equal(observed$dominant_state, "up")
  expect_equal(observed$ssi, 0.5000000001)
  expect_equal(as.character(observed$stability_region), "Transitional")
  expect_equal(as.character(observed$signal_class), "Transitional Up")
})
