library(caret)

decide.selection <- function(selection, boots, train.x, train.y, n.cores, positive.class=NULL) {
  if (selection == "lasso") {
    mclapply(seq_len(boots), one.lasso.boot, x=train.x, y=train.y, mc.cores=n.cores)
  } else if (selection == "mrmr") {
    mclapply(seq_len(boots), one.mrmr.boot, x=train.x, y=train.y, positive.class=positive.class, mc.cores=n.cores)
  } else if (selection == "t-test") {
    mclapply(seq_len(boots), one.t.boot, x=train.x, y=train.y, mc.cores=n.cores)
  } else if (selection == "rf-rfe") {
    mclapply(seq_len(boots), one.rf.rfe.boot, x=train.x, y=train.y, mc.cores=n.cores)
  }
  
  
}


one.lasso.boot <- function(b, x, y) {
  # set.seed(boot.seeds[b])
  
  x <- as.matrix(x)
  
  col.sel.data <- rep(F, ncol(x))
  names(col.sel.data) <- colnames(x)
  
  boot.i <- sample(1:nrow(x), nrow(x),
                   replace=T)
  sel.model <- cv.glmnet(x[boot.i,,drop=F],
                         y[boot.i],
                         alpha=1,
                         family=binomial)
  # Extract coefficients, only obtain those that are non-zero
  all.coefs <- coef(sel.model, s="lambda.1se")
  sel.coefs <- rownames(all.coefs)[which(all.coefs != 0 & rownames(all.coefs) != "(Intercept)")]
  col.sel.data[names(col.sel.data) %in% sel.coefs] <- T
  col.sel.data
}

one.mrmr.boot <- function(b, x, y, positive.class) {
  x <- as.matrix(x)
  col.sel.data <- rep(F, ncol(x))
  names(col.sel.data) <- colnames(x)
  
  boot.i <- sample(1:nrow(x), nrow(x),
                   replace=T)
  
  mrmr.dat <- data.frame(outcome=ifelse(y[boot.i] == positive.class, 1, 0),
                         x[boot.i,,drop=F])
  
  mrmr.dat <- mRMR.data(data=mrmr.dat)
  mrmr.fit <- mRMR.classic(
    "mRMRe.Filter",
    data=mrmr.dat,
    target_indices=1,
    feature_count=25
  )
  
  selected.idx <- solutions(mrmr.fit)[[1]][,1]
  selected.names <- featureNames(mrmr.dat)[selected.idx]
  
  col.sel.data[names(col.sel.data) %in% selected.names] <- T
  col.sel.data
  
  
}


one.t.boot <- function(b, x, y) {
  
  col.sel.data <- rep(F, ncol(x))
  names(col.sel.data) <- colnames(x)
  
  boot.i <- sample(1:nrow(x), nrow(x),
                   replace=T)
  
  p.val.list <- c()
  for (j in 1:ncol(x)) {
    t.matrix <- data.frame("y"=y[boot.i],
                           "x"=x[boot.i,j])
    
    lm.model <- lm(x~y, data=t.matrix)
    p.val <- summary(lm.model)$coefficients[2, "Pr(>|t|)"]
    p.val.list <- c(p.val.list, p.val)
  }
  p.val.list <- p.adjust(p.val.list, method="BH")
  
  selected.names <- colnames(x)[which(p.val.list <= 0.05)]
  
  col.sel.data[names(col.sel.data) %in% selected.names] <- T
  col.sel.data
}

one.rf.rfe.boot <- function(b, x, y, positive.class) {
  x <- as.matrix(x)
  col.sel.data <- rep(F, ncol(x))
  names(col.sel.data) <- colnames(x)
  
  boot.i <- sample(1:nrow(x), nrow(x),
                   replace=T)
  
  control.rfe <- rfeControl(
    functions=rfFuncs,
    method="repeatedcv",
    number=10,
    repeats=5
  )
  
  set.seed(42)
  
  rfe.profile <- rfe(
    x=x,
    y=y,
    rfeControl=control.rfe
  )
  
  selected.names <- rfe.profile$optVariables
  
  col.sel.data[names(col.sel.data) %in% selected.names] <- T
  col.sel.data
}




