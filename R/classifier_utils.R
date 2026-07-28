require_classifier_backend <- function() {
  if (requireNamespace("ranger", quietly = TRUE)) return("ranger")
  if (requireNamespace("randomForest", quietly = TRUE)) return("randomForest")
  stop("Install either 'ranger' or 'randomForest' for classification.", call. = FALSE)
}

binary_auc <- function(truth, score) {
  if (requireNamespace("pROC", quietly = TRUE)) {
    return(as.numeric(pROC::auc(truth, score)))
  }
  NA_real_
}

make_stratified_folds <- function(y, k = 5, seed = 1) {
  set.seed(seed)
  y <- as.factor(y)
  idx_by_class <- split(seq_along(y), y)
  folds <- vector("list", k)
  for (i in seq_len(k)) folds[[i]] <- integer()
  for (cls in names(idx_by_class)) {
    idx <- sample(idx_by_class[[cls]])
    chunks <- split(idx, rep(seq_len(k), length.out = length(idx)))
    for (i in seq_len(k)) folds[[i]] <- c(folds[[i]], chunks[[i]])
  }
  folds
}

run_binary_cv_classifier <- function(df, meta, outcome = "Target_Condition", k = 5, seed = 1) {
  backend <- require_classifier_backend()
  dat <- df %>% dplyr::inner_join(meta, by = "sample_id")
  dat <- dat %>% dplyr::filter(!is.na(.data[[outcome]]))
  y <- as.factor(dat[[outcome]])
  if (nlevels(y) != 2) stop("Classifier currently expects a binary outcome", call. = FALSE)

  x <- dat %>% dplyr::select(-dplyr::all_of(c("sample_id", outcome)))
  folds <- make_stratified_folds(y, k = k, seed = seed)
  res <- list()

  for (i in seq_along(folds)) {
    test_idx <- folds[[i]]
    train_idx <- setdiff(seq_len(nrow(dat)), test_idx)
    x_train <- x[train_idx, , drop = FALSE]
    x_test <- x[test_idx, , drop = FALSE]
    y_train <- y[train_idx]
    y_test <- y[test_idx]

    if (backend == "ranger") {
      fit <- ranger::ranger(
        y = y_train,
        x = as.data.frame(x_train),
        probability = TRUE,
        num.trees = 500,
        seed = seed + i
      )
      pred <- predict(fit, data = as.data.frame(x_test))$predictions[, 2]
    } else {
      fit <- randomForest::randomForest(x = as.data.frame(x_train), y = y_train, ntree = 500)
      pred <- predict(fit, newdata = as.data.frame(x_test), type = "prob")[, 2]
    }

    pred_class <- levels(y)[1 + as.integer(pred >= 0.5)]
    acc <- mean(pred_class == as.character(y_test))
    auc <- binary_auc(y_test, pred)
    res[[i]] <- tibble::tibble(fold = i, accuracy = acc, auc = auc)
  }

  dplyr::bind_rows(res)
}
