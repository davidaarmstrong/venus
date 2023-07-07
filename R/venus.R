 venus <- function(formula, 
                   nuisance = NULL, 
                   data = NULL,
                   method = c("earth", "lm"), 
                   nuisanceArgs = NULL, 
                   mainArgs = NULL, 
                   svd.thresh=.99, 
                   ps_method = c("po", "ds"), 
                   ...) {
  if (is.null(nuisance)) do.call(earth, c(formula = formula, mainArgs))
  else {
    stopifnot(inherits(nuisance, "formula"),
              inherits(formula, "formula"),
              length(nuisance) == 3,
              identical(nuisance[[2]], formula[[2]]),
              is.character(method),
              length(method) == 2)

    method <- match.arg(method, several.ok = TRUE)

    # Do the fit of the nuisance parameters
    nuisanceCall <- quote(earth(formula = nuisance, data = data))
    if (method[1] == "lm")
      nuisanceCall[[1]] <- as.name("lm")
    nm <- names(nuisanceArgs)
    for (i in seq_along(nuisanceArgs)) {
      if (!is.null(nm) && !is.na(nm[i]) && nchar(nm[i]))
        nuisanceCall[[nm[i]]] <- nuisanceArgs[[i]]
      else
        nuisanceCall <- c(nuisanceCall, pairlist(nuisanceArgs[[i]]))
    }
    nuisanceFit <- eval(nuisanceCall)
    yResids <- residuals(nuisanceFit)
    nuisanceModelmatrix <- model.matrix(nuisanceFit)
    nuisanceIntercept <- attr(terms(nuisance, data = data), "intercept") > 0
    ctrlMat <- nuisanceModelmatrix
    avform <- all.vars(formula)[-1]
    Xfits <- vector(mode="list", length=length(avform))
    for(j in 1:length(avform)){
      f <- update(nuisance, as.formula(paste0(avform[j], "~ .")))
      nxCall <- quote(earth(formula = f, data = data))
      nm <- names(nuisanceArgs)
      for (i in seq_along(nuisanceArgs)) {
        if (!is.null(nm) && !is.na(nm[i]) && nchar(nm[i]))
          nxCall[[nm[i]]] <- nuisanceArgs[[i]]
        else
          nxCall <- c(nxCall, pairlist(nuisanceArgs[[i]]))
      }
      Xfits[[j]] <- eval(nxCall)
     ndiff <- setdiff(colnames(model.matrix(Xfits[[j]])), colnames(ctrlMat))
     ctrlMat <- cbind(ctrlMat, model.matrix(Xfits[[j]])[, ndiff])
    }
    # Use the same linear model to correct the predictor variables
    mainModelmatrix <- model.matrix(formula, data = data)
    mainIntercept <- attr(terms(formula, data = data), "intercept") > 0
    if (nuisanceIntercept && mainIntercept)
      mainModelmatrix <- mainModelmatrix[, -1, drop = FALSE]
    novars <- apply(ctrlMat, 1, sd)
    if(any(novars == 0)){
      ctrlMat <- ctrlMat[,-which(novars == 0)]
    }
    scm <- svd(ctrlMat)
    scm.pct <- cumsum(scm$d/sum(scm$d))
    scm.cut <- min(which(scm.pct >= svd.thresh))
    sMat <- scm$u[,1:scm.cut, drop=FALSE] %*% diag(scm$d[1:scm.cut])[,,drop=FALSE]
    mainCall <- list()
    mnrg <- names(mainArgs)
    for (i in seq_along(mainArgs)) {
      if (!is.null(mnrg) && !is.na(mnrg[i]) && nchar(mnrg[i])) 
        mainCall[[nm[i]]] <- mainArgs[[i]]
      else mainCall <- c(mainCall, pairlist(mainArgs[[i]]))
    }
    if(psm == "po"){
      predictorFit <- lm(mainModelmatrix ~ sMat)
      mainModelmatrix <- residuals(predictorFit)
      mainFormula <- yResids ~ mainModelmatrix
      mainCall$formula <- mainFormula
      mainFit <- do.call(method[2], mainCall)
      res <- structure(list(mainFit = mainFit, nuisanceFit = nuisanceFit, 
                            predictorFit = predictorFit), class = "venus")
    }
    if(psm == "ds"){
      maindf <- get_all_vars(formula, data=data)
      sMat <- as.data.frame(sMat)
      maindf <- cbind(maindf, sMat)
      mainFormula <- reformulate(names(maindf)[-1], 
                                 names(maindf)[1])
      mainCall$formula <- mainFormula
      mainCall$data <- maindf
      mainFit <- do.call(method[2], mainCall)
      res <- structure(list(mainFit = mainFit, nuisanceFit = nuisanceFit), class = "venus")
    }
  }
   res 
 }
 
 
print.venus <- function(x, types = c("nuisance", "main"), ...) {
  types <- match.arg(types, choices = c("nuisance", "predictor", "main"), several.ok = TRUE)
  for (i in seq_along(types)) {
    type <- types[i]
    if (i > 1)
      cat("\n")
    if (type == "nuisance") {
      cat("The nuisance fit:\n")
      print(x$nuisanceFit, ...)
    }
    if (type == "predictor" & "predictorFit" %in% names(x)) {
      cat("The predictor fit:\n")
      print(x$predictorFit, ...)
    }
    if (type == "main") {
      cat("The main fit:\n")
      print(x$mainFit, ...)
    }
  }
}

plot.venus <- function(x, types = c("nuisance", "main"), ...) {
  types <- match.arg(types, choices = c("nuisance", "predictor", "main"), several.ok = TRUE)
  for (type in types) {
    if (type == "nuisance") 
      plot(x$nuisanceFit, ...)
    else if (type == "main") 
      plot(x$mainFit, ...)
    else if (type == "predictor" & "predictorFit" %in% names(x))
      plot(x$predictorFit)
  }
}

summary.venus <- function(object, ...) 
  structure(list(nuisanceSummary = summary(object$nuisanceFit, ...),
      mainSummary = summary(object$mainFit, ...),
      predictorSummary = summary(object$predictorFit, ...)),
      class = "summary.venus")

print.summary.venus <- function(x, types = c("main"), ...) {
  types <- match.arg(types, choices = c("nuisance", "predictor", "main"), several.ok = TRUE)
  for (i in seq_along(types)) {
    type <- types[i]
    if (i > 1)
      cat("\n")
    if (type == "nuisance") {
      cat("Nuisance summary:\n")
      print(x$nuisanceSummary, ...)
    }
    if (type == "predictor" & "predictorSummary" %in% names(x)) {
      cat("Predictor summary:\n")
      print(x$predictorSummary, ...)
    }
    if (type == "main") {
      cat("Main summary:\n")
      print(x$mainSummary, ...)
    }
  }
}

residuals.venus <- function(object, type = "main", ...) {
  type <- match.arg(type, choices = c("nuisance", "predictor", "main"))
  if (type == "main")
    residuals(object$mainFit, ...)
  else if (type == "nuisance")
    residuals(object$nuisanceFit, ...)
  else if (type == "predictor" & "predictorFit" %in% names(object))
    residuals(object$predictorFit, ...)
}
 
coef.venus <- function(object, type = "main", ...) {
  type <- match.arg(type, choices = c("nuisance", "predictor", "main"))
  if (type == "main")
    coef(object$mainFit, ...)
  else if (type == "nuisance")
    coef(object$nuisanceFit, ...)
  else if (type == "predictor" & "predictorFit" %in% names(object))
    coef(object$predictorFit, ...)
}

effects.venus <- function(object, type = "main", ...) {
  type <- match.arg(type, choices = c("nuisance", "predictor", "main"))
  if (type == "main")
    effects(object$mainFit, ...)
  else if (type == "nuisance")
    effects(object$nuisanceFit, ...)
  else if (type == "predictor" & predictorFit %in% names(object))
    effects(object$predictorFit, ...)
}

fitted.venus <- function(object, type = "main", ...) {
  type <- match.arg(type, choices = c("nuisance", "predictor", "main"))
  if (type == "main")
    fitted(object$mainFit, ...)
  else if (type == "nuisance")
    fitted(object$nuisanceFit, ...)
  else if (type == "predictor" & predictorFit %in% names(object))
    fitted(object$predictorFit, ...)
}

vcov.venus <- function(object, type = "main", ...) {
  type <- match.arg(type, choices = c("nuisance", "predictor", "main"))
  if (type == "main")
    vcov(object$mainFit, ...)
  else if (type == "nuisance")
    vcov(object$nuisanceFit, ...)
  else if (type == "predictor" & "predictorFit" %in% names(object))
    vcov(object$predictorFit, ...)
}