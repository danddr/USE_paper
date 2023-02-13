Dsquared<-function (model = NULL, obs = NULL, pred = NULL, family = NULL, 
                    adjust = FALSE, npar = NULL) {
  model.provided <- ifelse(is.null(model), FALSE, TRUE)
  if (model.provided) {
    if (!("glm" %in% class(model))) 
      stop("'model' must be of class 'glm'.")
    if (!is.null(pred)) 
      message("Argument 'pred' ignored in favour of 'model'.")
    if (!is.null(obs)) 
      message("Argument 'obs' ignored in favour of 'model'.")
    obs <- model$y
    pred <- model$fitted.values
  }
  else {
    if (is.null(obs) | is.null(pred)) 
      stop("You must provide either 'obs' and 'pred', or a 'model' object of class 'glm'.")
    if (length(obs) != length(pred)) 
      stop("'obs' and 'pred' must have the same number of values (and in the same order).")
    if (is.null(family)) 
      stop("With 'obs' and 'pred' arguments (rather than a model object), you must also specify one of two model family options: 'binomial' or 'poisson' (in quotes).")
    else if (!is.character(family)) 
      stop("Argument 'family' must be provided as character (i.e. in quotes: 'binomial' or 'poisson').")
    else if (length(family) != 1 | !(family %in% c("binomial", 
                                                   "poisson"))) 
      stop("'family' must be either 'binomial' or 'poisson' (in quotes).")
    pred[pred == 0] <- 2e-16
    pred[pred == 1] <- 1 - 2e-16
    dat <- data.frame(obs, pred)
    n.in <- nrow(dat)
    dat <- na.omit(dat)
    n.out <- nrow(dat)
    if (n.out < n.in) 
      warning(n.in - n.out, " observations removed due to missing data; ", 
              n.out, " observations actually evaluated.")
    obs <- dat$obs
    pred <- dat$pred
    if (family == "binomial") {
      if (any(!(obs %in% c(0, 1)) | pred < 0 | pred > 1)) 
        stop("'binomial' family implies that 'obs' data should be binary (with values 0 or 1) and 'pred' data should be bounded between 0 and 1.")
      link <- log(pred/(1 - pred))
    }
    else if (family == "poisson") {
      if (any(obs%%1 != 0)) 
        stop("'poisson' family implies that 'obs' data should consist of whole numbers.")
      link <- log(pred)
    }
    model <- glm(obs ~ link, family = family)
  }
  D2 <- (model$null.deviance - model$deviance)/model$null.deviance
  if (adjust) {
    if (model.provided) {
      n <- length(model$y)
      p <- attributes(logLik(model))$df
    }
    else {
      if (is.null(npar)) 
        stop("Adjusted D-squared from 'obs' and 'pred' values (rather than a model object) requires specifying the number of parameters in the underlying model ('npar').")
      n <- length(na.omit(obs))
      p <- npar
    }
    D2 <- 1 - ((n - 1)/(n - p)) * (1 - D2)
  }
  return(D2)
}