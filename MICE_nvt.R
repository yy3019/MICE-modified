

combine_foo1 <- function(coef, vars, m, n){
  
  # coefficient of interest (gamma_x_hat)
  gamma_hat_MI <- mean(coef)
  n = n
  # The following code gives non-sensical results
  #Calculating W,B,U
  gamma_matrix <- matrix(coef, ncol= n, byrow=TRUE)
  mean_n_gamma_hat <- apply(gamma_matrix, 1, mean)
  
  W <- sum(sapply(1:m, function(x){
    sum((gamma_matrix[x,] - mean_n_gamma_hat[x])^2)
  }))/(m*(n-1))
  B <- sum((mean_n_gamma_hat - gamma_hat_MI)^2)/(m-1)
  
  U <- mean(vars)
  
  # Calculating variance of our coefficient of interest
  T_MI <- U - W + (1+1/m)*B - W/n
  if(T_MI < 0) T_MI <- (1+1/m)*B
  
  # Confidence intervals are done with t-distribution with these
  # degrees of freedom
  df <- 1/(((((1+1/m)*B)^2)/((m-1)*T_MI^2)) + 
             ((((1+1/n)*W)^2)/(m*(n-1)*T_MI^2))) #Note that here was the mistake
  if(T_MI < 0) df <- m-1
  
  CI_low <- gamma_hat_MI + qt(0.025, df=df)*sqrt(T_MI)
  CI_upp <- gamma_hat_MI + qt(0.975, df=df)*sqrt(T_MI)
  
  return(c(gamma_hat_MI, CI_low, CI_upp, T_MI, W, B, U, df))
}



dd1 = function(data1, true_x1, true_x2, true_z){
  
  result_x1 = data1[,1]
  result_x2 = data1[,4]
  result_z = data1[,7]
  
  CI_lowx1 = data1[,2]
  CI_upx1 = data1[,3]
  CI_lowx2 = data1[,5]
  CI_upx2 = data1[,6]
  CI_lowz = data1[,8]
  CI_upz = data1[,9]
  
  bias_x1 = round(abs(sum(result_x1 - true_x1)/S_times)*1000)
  rmse_x1 = round(abs(sqrt(sum((result_x1 - true_x1)^2)/S_times))*1000)
  n_x1 = round(mean(((true_x1 <= CI_lowx1) | (true_x1 >= CI_upx1))) * 1000)
  aw_x1 = round(mean(CI_upx1-CI_lowx1)*1000)
  
  
  bias_x2 = round(abs(sum(result_x2 - true_x2)/S_times)*1000)
  rmse_x2 = round(abs(sqrt(sum((result_x2 - true_x2)^2)/S_times))*1000)
  n_x2 = round(mean(((true_x2 <= CI_lowx2) | (true_x2 >= CI_upx2))) * 1000)
  aw_x2 = round(mean(CI_upx2-CI_lowx2)*1000)
  
  bias_z = round(abs(sum(result_z - true_z)/S_times)*1000)
  rmse_z = round(abs(sqrt(sum((result_z - true_z)^2)/S_times))*1000)
  n_z = round(mean(((true_z <= CI_lowz) | (true_z >= CI_upz))) * 1000)
  aw_z = round(mean(CI_upz-CI_lowz)*1000)
  
  return(cbind(
    c(bias_x1, rmse_x1, n_x1, aw_x1),
    c(bias_x2, rmse_x2, n_x2, aw_x2),
    c(bias_z, rmse_z, n_z, aw_z)))
}


combine_rubin <- function(coef, vars){
  
  # coefficient of interest (gamma_x_hat)
  gamma_thetad <- mean(coef)
  d = length(coef)                 
  # The following code gives non-sensical results
  #Calculating W,B,U
  
  W <- mean(vars)
  B <- sum((coef - gamma_thetad)^2)/(d-1)
  
  T_d <- W + (d+1)/d*B
  
  gamma_d = (1+1/d)*B/T_d
  
  df = (d-1)*(1+1/(d+1)*W/B)^2
  
  CI_low <- gamma_thetad + qt(0.025, df=df)*sqrt(T_d)
  CI_upp <- gamma_thetad + qt(0.975, df=df)*sqrt(T_d)
  
  return(c(gamma_thetad, CI_low, CI_upp))
}

combine_foo <- function(coef, vars, m, n){
  
  # coefficient of interest (gamma_x_hat)
  gamma_hat_MI <- mean(coef)
  n = n
  # The following code gives non-sensical results
  #Calculating W,B,U
  gamma_matrix <- matrix(coef, ncol= n, byrow=TRUE)
  mean_n_gamma_hat <- apply(gamma_matrix, 1, mean)
  
  W <- sum(sapply(1:m, function(x){
    sum((gamma_matrix[x,] - mean_n_gamma_hat[x])^2)
  }))/(m*(n-1))
  B <- sum((mean_n_gamma_hat - gamma_hat_MI)^2)/(m-1)
  
  U <- mean(vars)
  
  # Calculating variance of our coefficient of interest
  T_MI <- U - W + (1+1/m)*B - W/n
  if(T_MI < 0) T_MI <- (1+1/m)*B
  
  # Confidence intervals are done with t-distribution with these
  # degrees of freedom
  df <- 1/(((((1+1/m)*B)^2)/((m-1)*T_MI^2)) + 
             ((((1+1/n)*W)^2)/(m*(n-1)*T_MI^2))) #Note that here was the mistake
  if(T_MI < 0) df <- m-1
  
  CI_low <- gamma_hat_MI + qt(0.025, df=df)*sqrt(T_MI)
  CI_upp <- gamma_hat_MI + qt(0.975, df=df)*sqrt(T_MI)
  
  return(c(gamma_hat_MI, CI_low, CI_upp))
}


check.deprecated <- function(...) {
  # print warnings for deprecated argument names
  nms <- names(list(...))
  replace.args <- list(
    imputationMethod = "method",
    defaultImputationMethod = "defaultMethod",
    form = "formulas"
  )
  
  wrn <- names(replace.args) %in% nms
  if (any(wrn)) {
    for (i in which(wrn)) {
      msg <- paste0(
        "The '", names(replace.args)[i],
        "' argument is no longer supported. Please use '",
        replace.args[i], "' instead."
      )
      warning(msg)
    }
  }
  
  invisible(NULL)
}

check.data <- function(data, method) {
  check.dataform(data)
}


check.dataform <- function(data) {
  if (!(is.matrix(data) || is.data.frame(data))) {
    stop("Data should be a matrix or data frame", call. = FALSE)
  }
  if (ncol(data) < 2) {
    stop("Data should contain at least two columns", call. = FALSE)
  }
  data <- as.data.frame(data)
  mat <- sapply(data, is.matrix)
  df <- sapply(data, is.data.frame)
  if (any(mat)) {
    stop(
      "Cannot handle columns with class matrix: ",
      colnames(data)[mat]
    )
  }
  if (any(df)) {
    stop(
      "Cannot handle columns with class data.frame: ",
      colnames(data)[df]
    )
  }
  
  dup <- duplicated(colnames(data))
  if (any(dup)) {
    stop(
      "Duplicate names found: ",
      paste(colnames(data)[dup], collapse = ", ")
    )
  }
  data
}


check.m <- function(m) {
  m <- m[1L]
  if (!is.numeric(m)) {
    stop("Argument m not numeric", call. = FALSE)
  }
  m <- floor(m)
  if (m < 1L) {
    stop("Number of imputations (m) lower than 1.", call. = FALSE)
  }
  m
}


check.cluster <- function(data, predictorMatrix) {
  # stop if the cluster variable is a factor
  isclassvar <- apply(predictorMatrix == -2, 2, any)
  for (j in colnames(predictorMatrix)) {
    if (isclassvar[j] && lapply(data, is.factor)[[j]]) {
      stop("Convert cluster variable ", j, " to integer by as.integer()")
    }
  }
  TRUE
}

check.ignore <- function(ignore, data) {
  if (is.null(ignore)) {
    return(rep(FALSE, nrow(data)))
  }
  if (!is.logical(ignore)) {
    stop("Argument ignore not a logical.")
  }
  if (length(ignore) != nrow(data)) {
    stop(
      "length(ignore) (", length(ignore),
      ") does not match nrow(data) (", nrow(data), ")."
    )
  }
  if (sum(!ignore) < 10L) {
    warning(
      "Fewer than 10 rows for fitting the imputation model. Are you sure?",
      call. = FALSE
    )
  }
  ignore
}

check.newdata <- function(newdata, data) {
  if (is.null(newdata)) {
    stop("No newdata found.")
  }
  if (!is.data.frame(newdata)) {
    stop("newdata not a data.frame.")
  }
  newdata
}

check.data <- function(data, method) {
  check.dataform(data)
}


check.dataform <- function(data) {
  if (!(is.matrix(data) || is.data.frame(data))) {
    stop("Data should be a matrix or data frame", call. = FALSE)
  }
  if (ncol(data) < 2) {
    stop("Data should contain at least two columns", call. = FALSE)
  }
  data <- as.data.frame(data)
  mat <- sapply(data, is.matrix)
  df <- sapply(data, is.data.frame)
  if (any(mat)) {
    stop(
      "Cannot handle columns with class matrix: ",
      colnames(data)[mat]
    )
  }
  if (any(df)) {
    stop(
      "Cannot handle columns with class data.frame: ",
      colnames(data)[df]
    )
  }
  
  dup <- duplicated(colnames(data))
  if (any(dup)) {
    stop(
      "Duplicate names found: ",
      paste(colnames(data)[dup], collapse = ", ")
    )
  }
  data
}


check.m <- function(m) {
  m <- m[1L]
  if (!is.numeric(m)) {
    stop("Argument m not numeric", call. = FALSE)
  }
  m <- floor(m)
  if (m < 1L) {
    stop("Number of imputations (m) lower than 1.", call. = FALSE)
  }
  m
}


check.cluster <- function(data, predictorMatrix) {
  # stop if the cluster variable is a factor
  isclassvar <- apply(predictorMatrix == -2, 2, any)
  for (j in colnames(predictorMatrix)) {
    if (isclassvar[j] && lapply(data, is.factor)[[j]]) {
      stop("Convert cluster variable ", j, " to integer by as.integer()")
    }
  }
  TRUE
}

check.ignore <- function(ignore, data) {
  if (is.null(ignore)) {
    return(rep(FALSE, nrow(data)))
  }
  if (!is.logical(ignore)) {
    stop("Argument ignore not a logical.")
  }
  if (length(ignore) != nrow(data)) {
    stop(
      "length(ignore) (", length(ignore),
      ") does not match nrow(data) (", nrow(data), ")."
    )
  }
  if (sum(!ignore) < 10L) {
    warning(
      "Fewer than 10 rows for fitting the imputation model. Are you sure?",
      call. = FALSE
    )
  }
  ignore
}

check.newdata <- function(newdata, data) {
  if (is.null(newdata)) {
    stop("No newdata found.")
  }
  if (!is.data.frame(newdata)) {
    stop("newdata not a data.frame.")
  }
  newdata
}

make.blocks <- function(data,
                        partition = c("scatter", "collect", "void"),
                        calltype = "type") {
  if (is.vector(data) && !is.list(data)) {
    v <- as.list(as.character(data))
    names(v) <- as.character(data)
    ct <- rep(calltype, length(v))
    names(ct) <- names(v)
    attr(v, "calltype") <- ct
    return(v)
  }
  if (is.list(data) && !is.data.frame(data)) {
    v <- name.blocks(data)
    if (length(calltype) == 1L) {
      ct <- rep(calltype, length(v))
      names(ct) <- names(v)
      attr(v, "calltype") <- ct
    }
    else {
      ct <- calltype
      names(ct) <- names(v)
      attr(v, "calltype") <- ct
    }
    return(v)
  }
  data <- as.data.frame(data)
  partition <- match.arg(partition)
  switch(partition,
         scatter = {
           v <- as.list(names(data))
           names(v) <- names(data)
         },
         collect = {
           v <- list(names(data))
           names(v) <- "collect"
         },
         void = {
           v <- list()
         },
         {
           v <- as.list(names(data))
           names(v) <- names(data)
         }
  )
  if (length(calltype) == 1L) {
    ct <- rep(calltype, length(v))
    names(ct) <- names(v)
    attr(v, "calltype") <- ct
  }
  else {
    ct <- calltype
    names(ct) <- names(v)
    attr(v, "calltype") <- ct
  }
  v
}


name.blocks <- function(blocks, prefix = "B") {
  if (!is.list(blocks)) {
    return(make.blocks(blocks))
  }
  if (is.null(names(blocks))) names(blocks) <- rep("", length(blocks))
  inc <- 1
  for (i in seq_along(blocks)) {
    if (names(blocks)[i] != "") next
    if (length(blocks[[i]]) == 1) {
      names(blocks)[i] <- blocks[[i]][1]
    } else {
      names(blocks)[i] <- paste0(prefix, inc)
      inc <- inc + 1
    }
  }
  blocks
}

check.blocks <- function(blocks, data, calltype = "type") {
  data <- check.dataform(data)
  blocks <- name.blocks(blocks)
  
  # check that all variable names exists in data
  bv <- unique(unlist(blocks))
  notFound <- !bv %in% colnames(data)
  if (any(notFound)) {
    stop(paste(
      "The following names were not found in `data`:",
      paste(bv[notFound], collapse = ", ")
    ))
  }
  
  if (length(calltype) == 1L) {
    ct <- rep(calltype, length(blocks))
    names(ct) <- names(blocks)
    attr(blocks, "calltype") <- ct
  }
  else {
    ct <- calltype
    names(ct) <- names(blocks)
    attr(blocks, "calltype") <- ct
  }
  
  blocks
}

construct.blocks <- function(formulas = NULL, predictorMatrix = NULL) {
  blocks.f <- blocks.p <- NULL
  if (!is.null(formulas)) {
    if (!all(sapply(formulas, is.formula))) {
      return(NULL)
    }
    blocks.f <- name.blocks(lapply(name.formulas(formulas), lhs))
    ct <- rep("formula", length(blocks.f))
    names(ct) <- names(blocks.f)
    attr(blocks.f, "calltype") <- ct
    if (is.null(predictorMatrix)) {
      return(blocks.f)
    }
  }
  
  if (!is.null(predictorMatrix)) {
    if (is.null(row.names(predictorMatrix))) {
      stop("No row names in predictorMatrix", call. = FALSE)
    }
    blocks.p <- name.blocks(row.names(predictorMatrix))
    ct <- rep("type", length(blocks.p))
    names(ct) <- names(blocks.p)
    attr(blocks.p, "calltype") <- ct
    if (is.null(formulas)) {
      return(blocks.p)
    }
  }
  
  # combine into unique blocks
  blocknames <- unique(c(names(blocks.f), names(blocks.p)))
  keep <- setdiff(blocknames, names(blocks.f))
  blocks <- c(blocks.f, blocks.p[keep])
  ct <- c(
    rep("formula", length(formulas)),
    rep("type", length(keep))
  )
  names(ct) <- names(blocks)
  attr(blocks, "calltype") <- ct
  blocks
}


make.predictorMatrix <- function(data, blocks = make.blocks(data)) {
  data <- check.dataform(data)
  predictorMatrix <- matrix(1, nrow = length(blocks), ncol = ncol(data))
  dimnames(predictorMatrix) <- list(names(blocks), colnames(data))
  for (i in row.names(predictorMatrix)) {
    predictorMatrix[i, colnames(predictorMatrix) %in% i] <- 0
  }
  predictorMatrix
}

check.predictorMatrix <- function(predictorMatrix,
                                  data,
                                  blocks = NULL) {
  data <- check.dataform(data)
  
  if (!is.matrix(predictorMatrix)) {
    stop("predictorMatrix not a matrix", call. = FALSE)
  }
  if (any(dim(predictorMatrix) == 0L)) {
    stop("predictorMatrix has no rows or columns", call. = FALSE)
  }
  
  # if we have no blocks, restrict to square predictorMatrix
  if (is.null(blocks)) {
    if (nrow(predictorMatrix) != ncol(predictorMatrix)) {
      stop(paste(
        "If no blocks are specified, predictorMatrix must",
        "have same number of rows and columns"
      ),
      call. = FALSE
      )
    }
    if (is.null(dimnames(predictorMatrix))) {
      if (ncol(predictorMatrix) == ncol(data)) {
        dimnames(predictorMatrix) <- list(colnames(data), colnames(data))
      } else {
        stop("Missing row/column names in predictorMatrix", call. = FALSE)
      }
    }
    for (i in row.names(predictorMatrix)) {
      predictorMatrix[i, grep(paste0("^", i, "$"), colnames(predictorMatrix))] <- 0
    }
    return(predictorMatrix)
  }
  
  # check conforming arguments
  if (nrow(predictorMatrix) > length(blocks)) {
    stop(paste0(
      "predictorMatrix has more rows (", nrow(predictorMatrix),
      ") than blocks (", length(blocks), ")"
    ),
    call. = FALSE
    )
  }
  
  # borrow rownames from blocks if needed
  if (is.null(rownames(predictorMatrix)) &&
      nrow(predictorMatrix) == length(blocks)) {
    rownames(predictorMatrix) <- names(blocks)
  }
  if (is.null(rownames(predictorMatrix))) {
    stop("Unable to set row names of predictorMatrix", call. = FALSE)
  }
  
  # borrow blocknames from predictorMatrix if needed
  if (is.null(names(blocks)) &&
      nrow(predictorMatrix) == length(blocks)) {
    names(blocks) <- rownames(predictorMatrix)
  }
  if (is.null(names(blocks))) {
    stop("Unable to set names of blocks", call. = FALSE)
  }
  
  # check existence of row names in blocks
  found <- rownames(predictorMatrix) %in% names(blocks)
  if (!all(found)) {
    stop("Names not found in blocks: ",
         paste(rownames(predictorMatrix)[!found], collapse = ", "),
         call. = FALSE
    )
  }
  
  # borrow colnames from data if needed
  if (is.null(colnames(predictorMatrix)) &&
      ncol(predictorMatrix) == ncol(data)) {
    colnames(predictorMatrix) <- names(data)
  }
  if (is.null(colnames(predictorMatrix))) {
    stop("Unable to set column names of predictorMatrix", call. = FALSE)
  }
  
  # check existence of variable names on data
  found <- colnames(predictorMatrix) %in% names(data)
  if (!all(found)) {
    stop("Names not found in data: ",
         paste(colnames(predictorMatrix)[!found], collapse = ", "),
         call. = FALSE
    )
  }
  
  list(
    predictorMatrix = predictorMatrix,
    blocks = blocks
  )
}

make.formulas <- function(data, blocks = make.blocks(data),
                          predictorMatrix = NULL) {
  data <- check.dataform(data)
  formulas <- as.list(rep("~ 0", length(blocks)))
  names(formulas) <- names(blocks)
  
  for (h in names(blocks)) {
    y <- blocks[[h]]
    if (is.null(predictorMatrix)) {
      predictors <- colnames(data)
    } else {
      type <- predictorMatrix[h, ]
      predictors <- names(type)[type != 0]
    }
    x <- setdiff(predictors, y)
    formulas[[h]] <- paste(
      paste(y, collapse = "+"), "~",
      paste(c("0", x), collapse = "+")
    )
  }
  
  formulas <- lapply(formulas, as.formula)
  formulas
}


name.formulas <- function(formulas, prefix = "F") {
  if (!is.list(formulas)) {
    stop("Argument `formulas` not a list", call. = FALSE)
  }
  if (!all(sapply(formulas, is.formula) | sapply(formulas, is.list))) {
    stop("Not all elements in `formulas` are a formula or a list")
  }
  if (is.null(names(formulas))) names(formulas) <- rep("", length(formulas))
  inc <- 1
  for (i in seq_along(formulas)) {
    if (names(formulas)[i] != "") next
    # if (hasdot(formulas[[i]]) && is.null(data))
    #  stop("Formula with dot requires `data` argument", call. = FALSE)
    y <- lhs(formulas[[i]])
    if (length(y) == 1) {
      names(formulas)[i] <- y
    } else {
      names(formulas)[i] <- paste0(prefix, inc)
      inc <- inc + 1
    }
  }
  formulas
}


check.formulas <- function(formulas, data) {
  formulas <- name.formulas(formulas)
  formulas <- handle.oldstyle.formulas(formulas, data)
  formulas <- lapply(formulas, expand.dots, data)
  # escape if formula is list of two formula's
  if (any(sapply(formulas, is.list))) {
    return(formulas)
  }
  formulas <- lapply(formulas, as.formula)
  formulas
}

extend.formulas <- function(formulas, data, blocks, predictorMatrix = NULL,
                            auxiliary = TRUE,
                            include.intercept = FALSE,
                            ...) {
  # Extend formulas with predictorMatrix
  if (is.null(predictorMatrix)) {
    return(formulas)
  }
  for (h in names(blocks)) {
    type <- predictorMatrix[h, ]
    predictors <- names(type)[type != 0]
    ff <- extend.formula(
      formula = formulas[[h]],
      predictors = predictors,
      auxiliary = auxiliary,
      include.intercept = include.intercept
    )
    formulas[[h]] <- ff
  }
  formulas
}

extend.formula <- function(formula = ~0,
                           predictors = NULL,
                           auxiliary = TRUE,
                           include.intercept = FALSE, ...) {
  if (!is.formula(formula)) formula <- ~0
  
  # handle dot in RHS
  if (hasdot(formula)) {
    if (length(predictors) > 1) {
      fr <- as.formula(c("~", paste(predictors, collapse = "+")))
    } else {
      fr <- ~0
    }
  } else {
    fr <- reformulate(c(".", predictors))
  }
  
  if (auxiliary) formula <- update(formula, fr, ...)
  if (include.intercept) formula <- update(formula, ~ . + 1, ...)
  formula
}

make.visitSequence <- function(data = NULL, blocks = NULL) {
  if (!is.null(blocks)) {
    blocks <- name.blocks(blocks)
    return(names(blocks))
  }
  
  data <- check.dataform(data)
  blocks <- make.blocks(data)
  names(blocks)
}

check.visitSequence <- function(visitSequence = NULL,
                                data, where = NULL, blocks) {
  if (is.null(names(blocks)) || any(is.na(names(blocks)))) {
    stop("Missing names in `blocks`.")
  }
  
  if (is.null(visitSequence)) {
    return(make.visitSequence(data, blocks))
  }
  
  if (is.null(where)) where <- is.na(data)
  nimp <- nimp(where, blocks)
  if (length(nimp) == 0) visitSequence <- nimp
  
  if (length(visitSequence) == 1 && is.character(visitSequence)) {
    code <- match.arg(visitSequence,
                      choices = c("roman", "arabic", "monotone", "revmonotone")
    )
    visitSequence <- switch(
      code,
      roman = names(blocks)[nimp > 0],
      arabic = rev(names(blocks)[nimp > 0]),
      monotone = names(blocks)[order(nimp)],
      revmonotone = rev(names(blocks)[order(nimp)])
    )
  }
  
  # legacy handling
  if (is.numeric(visitSequence)) {
    visitSequence <- colnames(data)[visitSequence]
  }
  
  # check against names(blocks)
  visitSequence <- visitSequence[is.element(visitSequence, names(blocks))]
  
  # remove any blocks without missing data
  visitSequence <- names((nimp > 0L)[visitSequence])
  visitSequence
}

handle.oldstyle.formulas <- function(formulas, data) {
  # converts old-style character vector to formula list
  oldstyle <- length(formulas) == ncol(data) && is.vector(formulas) &&
    is.character(formulas)
  if (!oldstyle) {
    return(formulas)
  }
  formulas[formulas != ""] <- "~ 0"
  fl <- as.list(formulas)
  names(fl) <- names(formulas)
  fl
}


is.empty.model.data <- function(x, data) {
  tt <- terms(x, data = data)
  (length(attr(tt, "factors")) == 0L) & (attr(tt, "intercept") == 0L)
}

lhs <- function(x) all.vars(update(x, . ~ 1))

is.formula <- function(x) {
  inherits(x, "formula")
}

hasdot <- function(f) {
  if (is.recursive(f)) {
    return(any(sapply(as.list(f), hasdot)))
  } else {
    f == as.symbol(".")
  }
}

expand.dots <- function(formula, data) {
  if (!is.formula(formula)) {
    return(formula)
  }
  if (!hasdot(formula)) {
    return(formula)
  }
  
  y <- lhs(formula)
  x <- setdiff(colnames(data), y)
  fs <- paste(paste(y, collapse = "+"), "~", paste(x, collapse = "+"))
  as.formula(fs)
}

make.where <- function(data,
                       keyword = c("missing", "all", "none", "observed")) {
  keyword <- match.arg(keyword)
  
  data <- check.dataform(data)
  where <- switch(keyword,
                  missing = is.na(data),
                  all = matrix(TRUE, nrow = nrow(data), ncol = ncol(data)),
                  none = matrix(FALSE, nrow = nrow(data), ncol = ncol(data)),
                  observed = !is.na(data)
  )
  
  dimnames(where) <- dimnames(data)
  where
}


check.where <- function(where, data, blocks) {
  if (is.null(where)) {
    where <- make.where(data, keyword = "missing")
  }
  
  if (!(is.matrix(where) || is.data.frame(where))) {
    if (is.character(where)) {
      return(make.where(data, keyword = where))
    } else {
      stop("Argument `where` not a matrix or data frame", call. = FALSE)
    }
  }
  if (!all(dim(data) == dim(where))) {
    stop("Arguments `data` and `where` not of same size", call. = FALSE)
  }
  
  where <- as.logical(as.matrix(where))
  if (anyNA(where)) {
    stop("Argument `where` contains missing values", call. = FALSE)
  }
  
  where <- matrix(where, nrow = nrow(data), ncol = ncol(data))
  dimnames(where) <- dimnames(data)
  where[, !colnames(where) %in% unlist(blocks)] <- FALSE
  where
}

make.method <- function(data,
                        where = make.where(data),
                        blocks = make.blocks(data),
                        defaultMethod = c("pmm", "logreg", "polyreg", "polr")) {
  method <- rep("", length(blocks))
  names(method) <- names(blocks)
  for (j in names(blocks)) {
    yvar <- blocks[[j]]
    y <- data[, yvar]
    def <- sapply(y, assign.method)
    k <- ifelse(all(diff(def) == 0), k <- def[1], 1)
    method[j] <- defaultMethod[k]
  }
  nimp <- nimp(where, blocks)
  method[nimp == 0] <- ""
  method
}


check.method <- function(method, data, where, blocks, defaultMethod) {
  if (is.null(method)) {
    return(make.method(
      data = data,
      where = where,
      blocks = blocks,
      defaultMethod = defaultMethod
    ))
  }
  nimp <- nimp(where, blocks)
  
  # expand user's imputation method to all visited columns
  # single string supplied by user (implicit assumption of two columns)
  if (length(method) == 1) {
    if (is.passive(method)) {
      stop("Cannot have a passive imputation method for every column.")
    }
    method <- rep(method, length(blocks))
    method[nimp == 0] <- ""
  }
  
  # check the length of the argument
  if (length(method) != length(blocks)) {
    stop("Length of method differs from number of blocks", call. = FALSE)
  }
  
  # add names to method
  names(method) <- names(blocks)
  
  # check whether the requested imputation methods are on the search path
  active.check <- !is.passive(method) & nimp > 0 & method != ""
  passive.check <- is.passive(method) & nimp > 0 & method != ""
  check <- all(active.check) & any(passive.check)
  if (check) {
    fullNames <- rep.int("mice.impute.passive", length(method[passive.check]))
  } else {
    fullNames <- paste("mice.impute", method[active.check], sep = ".")
    if (length(method[active.check]) == 0) fullNames <- character(0)
  }
  
  # type checks on built-in imputation methods
  for (j in names(blocks)) {
    vname <- blocks[[j]]
    y <- data[, vname, drop = FALSE]
    mj <- method[j]
    mlist <- list(
      m1 = c("logreg", "logreg.boot", "polyreg", "lda", "polr"),
      m2 = c(
        "norm", "norm.nob", "norm.predict", "norm.boot",
        "mean", "2l.norm", "2l.pan",
        "2lonly.norm", "2lonly.pan",
        "quadratic", "ri"
      ),
      m3 = c(
        "norm", "norm.nob", "norm.predict", "norm.boot",
        "mean", "2l.norm", "2l.pan",
        "2lonly.norm", "2lonly.pan",
        "quadratic", "logreg", "logreg.boot"
      )
    )
    cond1 <- sapply(y, is.numeric)
    cond2 <- sapply(y, is.factor) & sapply(y, nlevels) == 2
    cond3 <- sapply(y, is.factor) & sapply(y, nlevels) > 2
    if (any(cond1) && mj %in% mlist$m1) {
      warning("Type mismatch for variable(s): ",
              paste(vname[cond1], collapse = ", "),
              "\nImputation method ", mj, " is for categorical data.",
              call. = FALSE
      )
    }
    if (any(cond2) && mj %in% mlist$m2) {
      warning("Type mismatch for variable(s): ",
              paste(vname[cond2], collapse = ", "),
              "\nImputation method ", mj, " is not for factors.",
              call. = FALSE
      )
    }
    if (any(cond3) && mj %in% mlist$m3) {
      warning("Type mismatch for variable(s): ",
              paste(vname[cond3], collapse = ", "),
              "\nImputation method ", mj, " is not for factors with >2 levels.",
              call. = FALSE
      )
    }
  }
  method[nimp == 0] <- ""
  unlist(method)
}


# assign methods based on type,
# use method 1 if there is no single method within the block
assign.method <- function(y) {
  if (is.numeric(y)) {
    return(1)
  }
  if (nlevels(y) == 2) {
    return(2)
  }
  if (is.ordered(y) && nlevels(y) > 2) {
    return(4)
  }
  if (nlevels(y) > 2) {
    return(3)
  }
  if (is.logical(y)) {
    return(2)
  }
  1
}

make.post <- function(data) {
  post <- vector("character", length = ncol(data))
  names(post) <- colnames(data)
  post
}

check.post <- function(post, data) {
  if (is.null(post)) {
    return(make.post(data))
  }
  
  # check
  if (length(post) != ncol(data)) {
    stop("length(post) does not match ncol(data)", call. = FALSE)
  }
  
  # change
  if (is.null(names(post))) names(post) <- colnames(data)
  
  post
}

make.blots <- function(data, blocks = make.blocks(data)) {
  data <- check.dataform(data)
  blots <- vector("list", length(blocks))
  for (i in seq_along(blots)) blots[[i]] <- alist()
  names(blots) <- names(blocks)
  blots
}

check.blots <- function(blots, data, blocks = NULL) {
  data <- check.dataform(data)
  
  if (is.null(blots)) {
    return(make.blots(data, blocks))
  }
  
  blots <- as.list(blots)
  for (i in seq_along(blots)) blots[[i]] <- as.list(blots[[i]])
  
  if (length(blots) == length(blocks) && is.null(names(blots))) {
    names(blots) <- names(blocks)
  }
  blots
}

edit.setup <- function(data, setup,
                       allow.na = FALSE,
                       remove.constant = TRUE,
                       remove.collinear = TRUE,
                       remove_collinear = TRUE,
                       ...) {
  # legacy handling
  if (!remove_collinear) remove.collinear <- FALSE
  
  # edits the imputation model setup
  # When it detec constant or collinear variables, write in loggedEvents
  # and continues imputation with reduced model
  
  pred <- setup$predictorMatrix
  meth <- setup$method
  vis <- setup$visitSequence
  post <- setup$post
  
  # FIXME: this function is not yet adapted to blocks
  if (ncol(pred) != nrow(pred) || length(meth) != nrow(pred)
      || ncol(data) != nrow(pred)) {
    return(setup)
  }
  
  varnames <- colnames(data)
  
  # remove constant variables but leave passive variables untouched
  for (j in seq_len(ncol(data))) {
    if (!is.passive(meth[j])) {
      d.j <- data[, j]
      v <- if (is.character(d.j)) NA else var(as.numeric(d.j), na.rm = TRUE)
      constant <- if (allow.na) {
        if (is.na(v)) FALSE else v < 1000 * .Machine$double.eps
      } else {
        is.na(v) || v < 1000 * .Machine$double.eps
      }
      didlog <- FALSE
      if (constant && any(pred[, j] != 0) && remove.constant) {
        out <- varnames[j]
        pred[, j] <- 0
        updateLog(out = out, meth = "constant")
        didlog <- TRUE
      }
      if (constant && meth[j] != "" && remove.constant) {
        out <- varnames[j]
        pred[j, ] <- 0
        if (!didlog) {
          updateLog(out = out, meth = "constant")
        }
        meth[j] <- ""
        vis <- vis[vis != j]
        post[j] <- ""
      }
    }
  }
  
  ## remove collinear variables
  ispredictor <- apply(pred != 0, 2, any)
  if (any(ispredictor)) {
    droplist <- find.collinear(data[, ispredictor, drop = FALSE], ...)
  } else {
    droplist <- NULL
  }
  if (length(droplist) > 0) {
    for (k in seq_along(droplist)) {
      j <- which(varnames %in% droplist[k])
      didlog <- FALSE
      if (any(pred[, j] != 0) && remove.collinear) {
        # remove as predictor
        out <- varnames[j]
        pred[, j] <- 0
        updateLog(out = out, meth = "collinear")
        didlog <- TRUE
      }
      if (meth[j] != "" && remove.collinear) {
        out <- varnames[j]
        pred[j, ] <- 0
        if (!didlog) {
          updateLog(out = out, meth = "collinear")
        }
        meth[j] <- ""
        vis <- vis[vis != j]
        post[j] <- ""
      }
    }
  }
  
  if (all(pred == 0L)) {
    stop("`mice` detected constant and/or collinear variables. No predictors were left after their removal.")
  }
  
  setup$predictorMatrix <- pred
  setup$visitSequence <- vis
  setup$post <- post
  setup$method <- meth
  setup
}

is.mira <- function(x) {
  inherits(x, "mira")
}



is.mipo <- function(x) {
  inherits(x, "mipo")
}



is.mitml.result <- function(x) {
  inherits(x, "mitml.result")
}


is.passive <- function(string) {
  "~" == substring(string, 1, 1)
}

is.mads <- function(x) {
  inherits(x, "mads")
}

keep.in.model <- function(y, ry, x, wy) {
  (complete.cases(y, x) & ry) | (complete.cases(x) & wy)
}


impute.with.na <- function(x, wy) !complete.cases(x) & wy


check.df <- function(x, y, ry) {
  # if needed, writes the df warning message to the log
  df <- sum(ry) - ncol(x) - 1
  mess <- paste("df set to 1. # observed cases:", sum(ry), " # predictors:", ncol(x) + 1)
  if (df < 1 && sum(ry) > 0) {
    updateLog(out = mess, frame = 4)
  }
}


remove.lindep <- function(x, y, ry, eps = 1e-04, maxcor = 0.99,
                          allow.na = TRUE, frame = 4, ...) {
  # returns a logical vector of length ncol(x)
  
  if (ncol(x) == 0) {
    return(NULL)
  }
  
  # setting eps = 0 bypasses remove.lindep()
  if (eps == 0) {
    return(rep.int(TRUE, ncol(x)))
  }
  if (eps < 0) {
    stop("\n Argument 'eps' must be positive.")
  }
  
  # Keep all predictors if we allow imputation of fully missing y
  if (allow.na && sum(ry) == 0) {
    return(rep.int(TRUE, ncol(x)))
  }
  
  xobs <- x[ry, , drop = FALSE]
  yobs <- as.numeric(y[ry])
  if (var(yobs) < eps) {
    return(rep(FALSE, ncol(xobs)))
  }
  
  keep <- unlist(apply(xobs, 2, var) > eps)
  keep[is.na(keep)] <- FALSE
  highcor <- suppressWarnings(unlist(apply(xobs, 2, cor, yobs) < maxcor))
  keep <- keep & highcor
  if (all(!keep)) {
    updateLog(
      out = "All predictors are constant or have too high correlation.",
      frame = frame
    )
  }
  
  # no need to calculate correlations, so return
  k <- sum(keep)
  if (k <= 1L) {
    return(keep)
  } # at most one TRUE
  
  # correlation between x's
  cx <- cor(xobs[, keep, drop = FALSE], use = "all.obs")
  eig <- eigen(cx, symmetric = TRUE)
  ncx <- cx
  while (eig$values[k] / eig$values[1] < eps) {
    j <- seq_len(k)[order(abs(eig$vectors[, k]), decreasing = TRUE)[1]]
    keep[keep][j] <- FALSE
    ncx <- cx[keep[keep], keep[keep], drop = FALSE]
    k <- k - 1
    eig <- eigen(ncx)
  }
  if (!all(keep)) {
    out <- paste(dimnames(x)[[2]][!keep], collapse = ", ")
    updateLog(out = out, frame = frame)
  }
  return(keep)
}


## make list of collinear variables to remove
find.collinear <- function(x, threshold = 0.999, ...) {
  nvar <- ncol(x)
  x <- data.matrix(x)
  r <- !is.na(x)
  nr <- apply(r, 2, sum, na.rm = TRUE)
  ord <- order(nr, decreasing = TRUE)
  xo <- x[, ord, drop = FALSE] ## SvB 24mar2011
  varnames <- dimnames(xo)[[2]]
  z <- suppressWarnings(cor(xo, use = "pairwise.complete.obs"))
  hit <- outer(seq_len(nvar), seq_len(nvar), "<") & (abs(z) >= threshold)
  out <- apply(hit, 2, any, na.rm = TRUE)
  return(varnames[out])
}


updateLog <- function(out = NULL, meth = NULL, frame = 1) {
  
  # find structures defined a mice() level
  pos_state <- ma_exists("state", frame)$pos
  pos_loggedEvents <- ma_exists("loggedEvents", frame)$pos
  
  s <- get("state", pos_state)
  r <- get("loggedEvents", pos_loggedEvents)
  
  rec <- data.frame(
    it = s$it,
    im = s$im,
    dep = s$dep,
    meth = if (is.null(meth)) s$meth else meth,
    out = if (is.null(out)) "" else out
  )
  
  if (s$log) {
    rec <- rbind(r, rec)
  }
  s$log <- TRUE
  assign("state", s, pos = pos_state, inherits = TRUE)
  assign("loggedEvents", rec, pos = pos_loggedEvents, inherits = TRUE)
  return()
}


sym <- function(x) {
  (x + t(x)) / 2
}


# This helper function was copied from
# https://github.com/alexanderrobitzsch/miceadds/blob/master/R/ma_exists.R
ma_exists <- function(x, pos, n_index = 1:8) {
  n_index <- n_index + 1
  is_there <- exists(x, where = pos)
  obj <- NULL
  if (is_there) {
    obj <- get(x, pos)
  }
  if (!is_there) {
    for (nn in n_index) {
      pos <- parent.frame(n = nn)
      is_there <- exists(x, where = pos)
      if (is_there) {
        obj <- get(x, pos)
        break
      }
    }
  }
  #--- output
  res <- list(is_there = is_there, obj = obj, pos = pos)
  return(res)
}

initialize.chain <- function(blocks, maxit, m) {
  vars <- unique(unlist(blocks))
  chain <- array(NA, dim = c(length(vars), maxit, m))
  dimnames(chain) <- list(
    vars,
    seq_len(maxit),
    paste("Chain", seq_len(m))
  )
  chain
}

handles.arg <- function(f, a = "data") {
  # determine whether function f handles argument a
  if (!is.function(f)) {
    return(FALSE)
  }
  a %in% names(formals(f))
}


handles.format <- function(fn) {
  # determine whether function fn handles the `format` argument
  f <- get(fn)
  handles.arg(f, "format")
}


norm.draw <- function(y, ry, x, rank.adjust = TRUE, ...) {
  return(.norm.draw(y, ry, x, rank.adjust = TRUE, ...))
}


.norm.draw <- function(y, ry, x, rank.adjust = TRUE, ...) {
  p <- estimice(x[ry, , drop = FALSE], y[ry], ...)
  sigma.star <- sqrt(sum((p$r)^2) / rchisq(1, p$df))
  #print(chol(p$v))
  beta.star <- p$c + (t(chol(p$v)) %*% rnorm(ncol(x))) * sigma.star
  #print(beta.star)
  parm <- list(p$c, beta.star, sigma.star, p$ls.meth)
  names(parm) <- c("coef", "beta", "sigma", "estimation")
  if (any(is.na(parm$coef)) & rank.adjust) {
    parm$coef[is.na(parm$coef)] <- 0
    parm$beta[is.na(parm$beta)] <- 0
  }
  parm
}



estimice <- function(x, y, ls.meth = "qr", ridge = 1e-05, ...) {
  df <- max(length(y) - ncol(x), 1)
  if (ls.meth == "qr") {
    qr <- lm.fit(x = x, y = y)
    c <- t(qr$coef)
    f <- qr$fitted.values
    r <- t(qr$residuals)
    v <- try(solve(as.matrix(crossprod(qr.R(qr$qr)))), silent = TRUE)
    if (inherits(v, "try-error")) {
      xtx <- as.matrix(crossprod(qr.R(qr$qr)))
      # calculate ridge penalty
      pen <- diag(xtx) * ridge
      # add ridge penalty to allow inverse of v
      v <- solve(xtx + diag(pen))
      mess <- paste0(
        "mice detected that your data are (nearly) multi-collinear.\n",
        "It applied a ridge penalty to continue calculations, but the results can be unstable.\n",
        "Does your dataset contain duplicates, linear transformation, or factors with unique respondent names?"
      )
      updateLog(out = mess, frame = 6)
      if (get.printFlag()) {
        cat("*")
      } # indicator of added ridge penalty in the printed iteration history
    }
    return(list(c = t(c), r = t(r), v = v, df = df, ls.meth = ls.meth))
  }
  if (ls.meth == "ridge") {
    xtx <- crossprod(x)
    pen <- ridge * diag(xtx)
    if (length(pen) == 1) {
      pen <- matrix(pen)
    }
    v <- solve(xtx + diag(pen))
    c <- t(y) %*% x %*% v
    r <- y - x %*% t(c)
    return(list(c = t(c), r = r, v = v, df = df, ls.meth = ls.meth))
  }
  if (ls.meth == "svd") {
    s <- svd(x)
    c <- s$v %*% ((t(s$u) %*% y) / s$d)
    f <- x %*% c
    r <- f - y
    v <- try(solve(s$v %*% diag(s$d)^2 %*% t(s$v)), silent = TRUE)
    if (inherits(v, "try-error")) {
      xtx <- s$v %*% diag(s$d)^2 %*% t(s$v)
      # calculate ridge penalty
      pen <- diag(xtx) * ridge
      # add ridge penalty to allow inverse of v
      v <- solve(xtx + diag(pen))
      mess <- paste0(
        "mice detected that your data are (nearly) multi-collinear.\n",
        "It applied a ridge penalty to continue calculations, but the results can be unstable.\n",
        "Does your dataset contain duplicates, linear transformation, or factors with unique respondent names?"
      )
      updateLog(out = mess, frame = 6)
      if (get.printFlag()) {
        cat("*")
      } # indicator of added ridge penalty in the printed iteration history
    }
    return(list(c = c, r = r, v = v, df = df, ls.meth = ls.meth))
  }
}

get.printFlag <- function(start = 4) {
  while (inherits(
    try(get("printFlag", parent.frame(start)), silent = TRUE),
    "try-error"
  )) {
    start <- start + 1
  }
  get("printFlag", parent.frame(start))
}

obtain.design <- function(data, formula = ~.) {
  mf <- model.frame(formula, data = data, na.action = na.pass)
  model.matrix(formula, data = mf)
}

update.design <- function(design, data, varname = ".") {
  # Updates columns of the design matrix related to variable
  # varname in data
  
  varname <- as.character(varname[1])
  idx <- attr(design, "assign") %in% grep(varname, names(data))
  
  # variable j not found
  if (varname == "" || !any(idx)) {
    return(design)
  }
  
  # create model frame of variable j only
  fj <- as.formula(paste("~", varname))
  mfj <- model.frame(fj, data = data, na.action = na.pass)
  design[, idx] <- model.matrix(fj, data = mfj)[, -1, drop = FALSE]
  design
}

initialize.imp <- function(data, m, ignore, where, blocks, visitSequence,
                           method, nmis, data.init) {
  imp <- vector("list", ncol(data))
  names(imp) <- names(data)
  r <- !is.na(data)
  for (h in visitSequence) {
    for (j in blocks[[h]]) {
      y <- data[, j]
      ry <- r[, j] & !ignore
      wy <- where[, j]
      imp[[j]] <- as.data.frame(matrix(NA, nrow = sum(wy), ncol = m))
      dimnames(imp[[j]]) <- list(row.names(data)[wy], 1:m)
      if (method[h] != "") {
        for (i in seq_len(m)) {
          #print(i)
          if (nmis[j] < nrow(data)) {
            if (is.null(data.init)) {
              imp[[j]][, i] <- mice.impute.sample(y, ry, wy = wy)
              #print(imp[[j]][, i])
            } else {
              imp[[j]][, i] <- data.init[wy, j]
            }
          } else {
            imp[[j]][, i] <- rnorm(nrow(data))
          }
        }
      }
    }
  }
  imp
}

mv = function(x, wy, parm){
  x[wy, ] %*% parm$beta + rnorm(sum(wy)) * parm$sigma
}

mice.impute.normv <- function(y, ry, x, wy = NULL,...) {
  if (is.null(wy)) wy <- !ry
  x <- cbind(1, as.matrix(x))
  parm <- .norm.draw(y, ry, x, ...)
  r1 = x[wy, ] %*% parm$beta + rnorm(sum(wy)) * parm$sigma
  
  r2 = c()
  
  for(z in 1:(n-1)){
    #set.seed(n)
    r2 = cbind(r2, mv(x, wy, parm))
  }
  r2 = cbind(r1, r2)
  return(list(r1 = r1, r2 = r2))
}

sampler <- function(data, m, ignore, where, imp, blocks, method, 
                    visitSequence, predictorMatrix, formulas, blots, 
                    post, fromto, printFlag, ...) {
  from <- fromto[1]
  to <- fromto[2]
  maxit <- to - from + 1
  r <- !is.na(data)
  
  # set up array for convergence checking
  chainMean <- chainVar <- initialize.chain(blocks, maxit, m)
  
  ## THE MAIN LOOP: GIBBS SAMPLER ##
  if (maxit < 1) iteration <- 0
  if (maxit >= 1) {
    if (printFlag) {
      cat("\n iter imp variable")
    }
    for (k in from:to) {
      # begin k loop : main iteration loop
      rr1 = c()
      rr2 = c()
      rr11 = c()
      rr22 = c()
      bbz = 1
      iteration <- k
      for (i in seq_len(m)) {
        # begin i loop: repeated imputation loop
        
        if (printFlag) {
          cat("\n ", iteration, " ", i)
        }
        
        # prepare the i'th imputation
        # do not overwrite any observed data
        for (h in visitSequence) {
          for (j in blocks[[h]]) {
            y <- data[, j]
            ry <- r[, j]
            wy <- where[, j]
            data[(!ry) & wy, j] <- imp[[j]][(!ry)[wy], i]
          }
        }
        
        # impute block-by-block
        for (h in visitSequence) {
          ct <- attr(blocks, "calltype")
          calltype <- ifelse(length(ct) == 1, ct[1], ct[h])
          
          b <- blocks[[h]]
          if (calltype == "formula") ff <- formulas[[h]] else ff <- NULL
          if (calltype == "type") type <- predictorMatrix[h, ] else type <- NULL
          
          user <- blots[[h]]
          
          # univariate/multivariate logic
          theMethod <- method[h]
          empt <- theMethod == ""
          univ <- !empt && !is.passive(theMethod) &&
            !handles.format(paste0("mice.impute.", theMethod))
          mult <- !empt && !is.passive(theMethod) &&
            handles.format(paste0("mice.impute.", theMethod))
          pass <- !empt && is.passive(theMethod) && length(blocks[[h]]) == 1
          if (printFlag & !empt) cat(" ", b)
          
          ## store current state
          oldstate <- get("state", pos = parent.frame())
          newstate <- list(
            it = k, im = i,
            dep = h,
            meth = theMethod,
            log = oldstate$log
          )
          assign("state", newstate, pos = parent.frame(), inherits = TRUE)
          
          # (repeated) univariate imputation - type method
          if (univ) {
            
            for (j in b) {
              #print(3)
              aa = sampler.univ(
                data = data, r = r, where = where,
                type = type, formula = ff,
                method = theMethod,
                yname = j, k = k,
                calltype = calltype,
                user = user, ignore = ignore,
                ...
              )
              
              imp[[j]][,i] <-
                aa$imputes
              
              if(nrow(aa$nn) == n_main){
                
                
                
                if(bbz %% 2 == 1){ 
                  #colnames(bb) = c(paste0(i,"_",1), paste0(i,"_",2), paste0(i,"_",3), paste0(i,"_",4))
                  rr1 = cbind(rr1, aa$nn)
                  #rr11 = cbind(rr11, aa$nn2)
                  }
                
                
                if(bbz %% 2 == 0){ 
                  #colnames(bb) = c(paste0(i,"_",1), paste0(i,"_",2), paste0(i,"_",3), paste0(i,"_",4))
                  rr2 = cbind(rr2, aa$nn)
                  #rr22 = cbind(rr22, aa$nn2)
                  }
                
                bbz = bbz + 1}
              
              
              
              data[(!r[, j]) & where[, j], j] <-
                imp[[j]][(!r[, j])[where[, j]], i]
              
              # optional post-processing
              cmd <- post[j]
              if (cmd != "") {
                eval(parse(text = cmd))
                data[where[, j], j] <- imp[[j]][, i]
              }
            }
          }
          
          # multivariate imputation - type and formula
          if (mult) {
            mis <- !r
            mis[, setdiff(colnames(data), b)] <- FALSE
            data[mis] <- NA
            
            fm <- paste("mice.impute", theMethod, sep = ".")
            if (calltype == "formula") {
              imputes <- do.call(fm, args = list(
                data = data,
                formula = ff, ...
              ))
            } else if (calltype == "type") {
              imputes <- do.call(fm, args = list(
                data = data,
                type = type, ...
              ))
            } else {
              stop("Cannot call function of type ", calltype,
                   call. = FALSE
              )
            }
            if (is.null(imputes)) {
              stop("No imputations from ", theMethod,
                   h,
                   call. = FALSE
              )
            }
            for (j in names(imputes)) {
              print(1)
              imp[[j]][, i] <- imputes[[j]]
              data[!r[, j], j] <- imp[[j]][, i]
            }
          }
          
          # passive imputation
          # applies to all rows, so no ignore needed
          if (pass) {
            for (j in b) {
              print(2)
              wy <- where[, j]
              ry <- r[, j]
              imp[[j]][, i] <- model.frame(as.formula(theMethod), data[wy, ],
                                           na.action = na.pass
              )
              data[(!ry) & wy, j] <- imp[[j]][(!ry)[wy], i]
            }
          }
        } # end h loop (blocks)
      } # end i loop (imputation number)
      
      # store means and sd of m imputes
      k2 <- k - from + 1L
      if (length(visitSequence) > 0L) {
        for (h in visitSequence) {
          for (j in blocks[[h]]) {
            if (!is.factor(data[, j])) {
              chainVar[j, k2, ] <- apply(imp[[j]], 2L, var, na.rm = TRUE)
              chainMean[j, k2, ] <- colMeans(as.matrix(imp[[j]]), na.rm = TRUE)
            }
            if (is.factor(data[, j])) {
              for (mm in seq_len(m)) {
                nc <- as.integer(factor(imp[[j]][, mm], levels = levels(data[, j])))
                chainVar[j, k2, mm] <- var(nc, na.rm = TRUE)
                chainMean[j, k2, mm] <- mean(nc, na.rm = TRUE)
              }
            }
          }
        }
      }
    } # end main iteration
    
    if (printFlag) {
      r <- get("loggedEvents", parent.frame(1))
      ridge.used <- any(grepl("A ridge penalty", r$out))
      if (ridge.used) {
        cat("\n * Please inspect the loggedEvents \n")
      } else {
        cat("\n")
      }
    }
  }
  list(iteration = maxit, imp = imp, chainMean = chainMean, chainVar = chainVar, rr1 = rr1, rr2 = rr2)
}


sampler.univ <- function(data, r, where, type, formula, method, yname, k,
                         calltype = "type", user, ignore, ...) {
  j <- yname[1L]
  
  if (calltype == "type") {
    vars <- colnames(data)[type != 0]
    pred <- setdiff(vars, j)
    if (length(pred) > 0L) {
      formula <- reformulate(pred, response = j)
      formula <- update(formula, ". ~ . ")
    } else {
      formula <- as.formula(paste0(j, " ~ 1"))
    }
  }
  
  if (calltype == "formula") {
    # move terms other than j from lhs to rhs
    ymove <- setdiff(lhs(formula), j)
    formula <- update(formula, paste(j, " ~ . "))
    if (length(ymove) > 0L) {
      formula <- update(formula, paste("~ . + ", paste(ymove, collapse = "+")))
    }
  }
  
  # get the model matrix
  x <- obtain.design(data, formula)
  
  # expand type vector to model matrix, remove intercept
  if (calltype == "type") {
    type <- type[labels(terms(formula))][attr(x, "assign")]
    x <- x[, -1L, drop = FALSE]
    names(type) <- colnames(x)
  }
  if (calltype == "formula") {
    x <- x[, -1L, drop = FALSE]
    type <- rep(1L, length = ncol(x))
    names(type) <- colnames(x)
  }
  
  # define y, ry and wy
  y <- data[, j]
  ry <- complete.cases(x, y) & r[, j] & !ignore
  wy <- complete.cases(x) & where[, j]
  
  # nothing to impute
  if (all(!wy)) {
    return(numeric(0))
  }
  
  cc <- wy[where[, j]]
  if (k == 1L) check.df(x, y, ry)
  
  # remove linear dependencies
  keep <- remove.lindep(x, y, ry, ...)
  x <- x[, keep, drop = FALSE]
  type <- type[keep]
  if (ncol(x) != length(type)) {
    stop("Internal error: length(type) != number of predictors")
  }
  
  # here we go
  f <- paste("mice.impute", method, sep = ".")
  imputes <- data[wy, j]
  imputes[!cc] <- NA

  y_1 = y
  ry_1 = ry
  x_1 = x
  wy_1 = wy
  user_1 = user
  l1 = list(...)
  args <- c(list(y = y, ry = ry, x = x, wy = wy, type = type), user, list(...))
  #args2 <- c(list(y = y_1, ry = ry_1, x = x_1, wy = wy_1, type = type), user_1, l1)
  rr = do.call(f, args = args)
  imputes[cc] <- rr$r1
  #print("a")
  #print(imputes[cc])
  nn = rr$r2
  
  #print(nn) 
  return(list(imputes = imputes, nn = nn))
}




mice_nv <- function(data,
                    m = 5,
                    method = NULL,
                    predictorMatrix,
                    ignore = NULL,
                    where = NULL,
                    blocks,
                    visitSequence = NULL,
                    formulas,
                    blots = NULL,
                    post = NULL,
                    defaultMethod = c("pmm", "logreg", "polyreg", "polr"),
                    maxit = 5,
                    printFlag = TRUE,
                    seed = NA,
                    data.init = NULL,
                    ...) {
  call <- match.call()
  check.deprecated(...)
  if (!is.na(seed)) set.seed(seed)
  
  # check form of data and m
  data <- check.dataform(data)
  m <- check.m(m)
  
  # determine input combination: predictorMatrix, blocks, formulas
  mp <- missing(predictorMatrix)
  mb <- missing(blocks)
  mf <- missing(formulas)
  
  # case A
  if (mp & mb & mf) {
    # blocks lead
    blocks <- make.blocks(colnames(data))
    predictorMatrix <- make.predictorMatrix(data, blocks)
    formulas <- make.formulas(data, blocks)
  }
  # case B
  if (!mp & mb & mf) {
    # predictorMatrix leads
    predictorMatrix <- check.predictorMatrix(predictorMatrix, data)
    blocks <- make.blocks(colnames(predictorMatrix), partition = "scatter")
    formulas <- make.formulas(data, blocks, predictorMatrix = predictorMatrix)
  }
  
  # case C
  if (mp & !mb & mf) {
    # blocks leads
    blocks <- check.blocks(blocks, data)
    predictorMatrix <- make.predictorMatrix(data, blocks)
    formulas <- make.formulas(data, blocks)
  }
  
  # case D
  if (mp & mb & !mf) {
    # formulas leads
    formulas <- check.formulas(formulas, data)
    blocks <- construct.blocks(formulas)
    predictorMatrix <- make.predictorMatrix(data, blocks)
  }
  
  # case E
  if (!mp & !mb & mf) {
    # predictor leads
    blocks <- check.blocks(blocks, data)
    z <- check.predictorMatrix(predictorMatrix, data, blocks)
    predictorMatrix <- z$predictorMatrix
    blocks <- z$blocks
    formulas <- make.formulas(data, blocks, predictorMatrix = predictorMatrix)
  }
  
  # case F
  if (!mp & mb & !mf) {
    # formulas lead
    formulas <- check.formulas(formulas, data)
    predictorMatrix <- check.predictorMatrix(predictorMatrix, data)
    blocks <- construct.blocks(formulas, predictorMatrix)
  }
  
  # case G
  if (mp & !mb & !mf) {
    # blocks lead
    blocks <- check.blocks(blocks, data, calltype = "formula")
    formulas <- check.formulas(formulas, blocks)
    predictorMatrix <- make.predictorMatrix(data, blocks)
  }
  
  # case H
  if (!mp & !mb & !mf) {
    # blocks lead
    blocks <- check.blocks(blocks, data)
    formulas <- check.formulas(formulas, data)
    predictorMatrix <- check.predictorMatrix(predictorMatrix, data, blocks)
  }
  
  chk <- check.cluster(data, predictorMatrix)
  where <- check.where(where, data, blocks)
  visitSequence <- check.visitSequence(visitSequence,
                                       data = data,
                                       where = where, blocks = blocks
  )
  method <- check.method(
    method = method, data = data, where = where,
    blocks = blocks, defaultMethod = defaultMethod
  )
  post <- check.post(post, data)
  blots <- check.blots(blots, data, blocks)
  ignore <- check.ignore(ignore, data)
  
  # data frame for storing the event log
  state <- list(it = 0, im = 0, dep = "", meth = "", log = FALSE)
  loggedEvents <- data.frame(it = 0, im = 0, dep = "", meth = "", out = "")
  
  # edit imputation setup
  setup <- list(
    method = method,
    predictorMatrix = predictorMatrix,
    visitSequence = visitSequence,
    post = post
  )
  setup <- edit.setup(data, setup, ...)
  method <- setup$method
  predictorMatrix <- setup$predictorMatrix
  visitSequence <- setup$visitSequence
  post <- setup$post
  
  # initialize imputations
  nmis <- apply(is.na(data), 2, sum)
  imp <- initialize.imp(
    data, m, ignore, where, blocks, visitSequence,
    method, nmis, data.init
  )
  
  # and iterate...
  from <- 1
  to <- from + maxit - 1
  q <- sampler(
    data, m, ignore, where, imp, blocks, method,
    visitSequence,predictorMatrix, formulas, blots,
    post, c(from, to), printFlag, ...
  )
  
  if (!state$log) loggedEvents <- NULL
  if (state$log) row.names(loggedEvents) <- seq_len(nrow(loggedEvents))
  
  ## save, and return
  midsobj <- list(
    data = data,
    imp = q$imp,
    rr1 = q$rr1,
    rr2 = q$rr2,
    m = m,
    where = where,
    blocks = blocks,
    call = call,
    nmis = nmis,
    method = method,
    predictorMatrix = predictorMatrix,
    visitSequence = visitSequence,
    formulas = formulas,
    post = post,
    blots = blots,
    ignore = ignore,
    seed = seed,
    iteration = q$iteration,
    lastSeedValue = .Random.seed,
    chainMean = q$chainMean,
    chainVar = q$chainVar,
    loggedEvents = loggedEvents,
    version = packageVersion("mice"),
    date = Sys.Date()
  )
  oldClass(midsobj) <- "mids"
  
  if (!is.null(midsobj$loggedEvents)) {
    warning("Number of logged events: ", nrow(midsobj$loggedEvents),
            call. = FALSE
    )
  }
  midsobj
}
