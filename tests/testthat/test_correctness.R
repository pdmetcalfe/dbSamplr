context("testing correctness")

test_that("matches gsDesign #1",
          {
    prob  <- 0.6
    sizes <- c(10, 20, 30)
    crits <- c(4, 12, 16)

    standard <-
        gsDesign::gsBinomialExact(3,
                                  theta=1 - prob,
                                  n.I=sizes,
                                  a=c(-1, -1, -1),
                                  b=sizes - crits)

    internal <- gsProbs(prob, sizes, crits)

    expect_equal(as.numeric(standard$upper$prob),
                 internal)
})

test_that("matches gsDesign #2",
          {
    prob  <- 0.7
    sizes <- c(10, 20, 30)
    crits <- c(8, 18, 25)

    standard <-
        gsDesign::gsBinomialExact(3,
                                  theta=1 - prob,
                                  n.I=sizes,
                                  a=c(-1, -1, -1),
                                  b=sizes - crits)

    internal <- gsProbs(prob, sizes, crits)

    expect_equal(as.numeric(standard$upper$prob),
                 internal)
})
