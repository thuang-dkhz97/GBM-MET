
factorise <- function(y, comparison) {
  if (comparison == "gbm.v.met") {
    factor(y, levels=c("MET", "GBM"))
  } else if (comparison == "gbm.v.hc") {
    factor(y, levels=c("HC", "GBM"))
  } else if (comparison == "met.v.hc") {
    factor(y, levels=c("HC", "MET"))
  }
}

get.positive.class <- function(comparison) {
  if (comparison == "gbm.v.met") {
    "GBM"
  } else if (comparison == "gbm.v.hc") {
    "GBM"
  } else if (comparison == "met.v.hc") {
    "MET"
  }
}

get.negative.class <- function(comparison) {
  if (comparison == "gbm.v.met") {
    "MET"
  } else if (comparison == "gbm.v.hc") {
    "HC"
  } else if (comparison == "met.v.hc") {
    "HC"
  }
}