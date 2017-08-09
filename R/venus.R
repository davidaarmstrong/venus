venus <- function(formula, nuisance = NULL, data = NULL,
                  method = c("earth", "lm"), nuisanceArgs = NULL, mainArgs = NULL) {
  if (is.null(nuisance)) do.call(earth, c(formula = formula, mainArgs))
  else {
    stopifnot(inherits(nuisance, "formula"),
              inherits(formula, "formula"),
              length(nuisance) == 3,
              identical(nuisance[[2]], formula[[2]]),
              is.character(method),
              length(method) == 2)

    # method <- match.arg(method, several.ok = TRUE)

    # Do the fit of the nuisance parameters
    nuisanceFit <- do.call(method[1], c(list(formula = nuisance, data = data), nuisanceArgs))
    yResids <- residuals(nuisanceFit)
    nuisanceModelmatrix <- model.matrix(nuisanceFit)
    nuisanceIntercept <- attr(terms(nuisance, data = data), "intercept") > 0

    # Use the same linear model to correct the predictor variables
    mainModelmatrix <- model.matrix(formula, data = data)
    mainIntercept <- attr(terms(formula, data = data), "intercept") > 0
    if (nuisanceIntercept && mainIntercept)
      mainModelmatrix <- mainModelmatrix[, -1, drop = FALSE]

    mainModelmatrix <- residuals(lm(mainModelmatrix ~ nuisanceModelmatrix))

    # Now fit the residuals of y to the residuals of the predictors
    if (method[2] == "lm")
      mainFormula <- yResids ~ mainModelmatrix - 1
    else
      mainFormula <- yResids ~ mainModelmatrix
    mainFit <- do.call(method[2], c(mainFormula, mainArgs))
    structure(list(mainFit = mainFit, nuisanceFit = nuisanceFit), class = "venus")
  }
}

print.venus <- function(x, ...) {
  cat("The nuisance fit:\n")
  print(x$nuisanceFit, ...)
  cat("\nThe main fit:\n")
  print(x$mainFit, ...)
}

plot.venus <- function(x, ...) {
  plot(x$nuisanceFit, ...)
  par(mfrow=c(2,2), oma = c(0,0,2,0))
  plot(x$mainFit, ...)
}

