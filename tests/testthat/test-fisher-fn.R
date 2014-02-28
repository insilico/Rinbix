testValues <- seq(from=-1, to=1, by=0.25)

require(Rinbix)
expect_that(fisherRtoZ(testValues[0]), equals(0.5 * log((1+testValues[0])/(1-testValues[0]))))
