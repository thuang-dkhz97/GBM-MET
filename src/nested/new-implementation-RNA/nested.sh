# ======================== Protein
# ======== FSM=Lasso, Model=Ridge =====================

Rscript './src/nested/new-implementation-RNA/nested-pipeline.R' -s 1000 -S "lasso" -B 100 -P 0 -k 25 -f 1.2 -c "gbm.v.met"
Rscript './src/nested/new-implementation-RNA/nested-pipeline.R' -s 1000 -S "lasso" -B 100 -P 0 -k 25 -f 1.2 -c "gbm.v.hc"
Rscript './src/nested/new-implementation-RNA/nested-pipeline.R' -s 1000 -S "lasso" -B 100 -P 0 -k 25 -f 1.2 -c "met.v.hc"


# ======== FSM=mrmr, Model=Ridge =====================


Rscript './src/nested/new-implementation-RNA/nested-pipeline.R' -s 1000 -S "mrmr" -B 100 -P 0 -k 25 -f 1.2 -c "gbm.v.met"
Rscript './src/nested/new-implementation-RNA/nested-pipeline.R' -s 1000 -S "mrmr" -B 100 -P 0 -k 25 -f 1.2 -c "gbm.v.hc"
Rscript './src/nested/new-implementation-RNA/nested-pipeline.R' -s 1000 -S "mrmr" -B 100 -P 0 -k 25 -f 1.2 -c "met.v.hc"

# ======== FSM=t-test, Model=Ridge =====================


Rscript './src/nested/new-implementation-RNA/nested-pipeline.R' -s 1000 -S "t-test" -B 100 -P 0 -k 25 -f 1.2 -c "gbm.v.met"
Rscript './src/nested/new-implementation-RNA/nested-pipeline.R' -s 1000 -S "t-test" -B 100 -P 0 -k 25 -f 1.2 -c "gbm.v.hc"
Rscript './src/nested/new-implementation-RNA/nested-pipeline.R' -s 1000 -S "t-test" -B 100 -P 0 -k 25 -f 1.2 -c "met.v.hc"

# ======== FSM=rf_rfe, Model=Ridge =====================


Rscript './src/nested/new-implementation-RNA/nested-pipeline.R' -s 1000 -S "rf-rfe" -B 100 -P 0 -k 25 -f 1.2 -c "gbm.v.met"
Rscript './src/nested/new-implementation-RNA/nested-pipeline.R' -s 1000 -S "rf-rfe" -B 100 -P 0 -k 25 -f 1.2 -c "gbm.v.hc"
Rscript './src/nested/new-implementation-RNA/nested-pipeline.R' -s 1000 -S "rf-rfe" -B 100 -P 0 -k 25 -f 1.2 -c "met.v.hc"

# ======== FSM=Lasso, Model=Rf =====================

Rscript './src/nested/new-implementation-RNA/nested-pipeline.R' -s 1000 -S "lasso" -m "rf" -B 100 -P 0 -k 25 -f 1.2 -c "gbm.v.met"
Rscript './src/nested/new-implementation-RNA/nested-pipeline.R' -s 1000 -S "lasso" -m "rf" -B 100 -P 0 -k 25 -f 1.2 -c "gbm.v.hc"
Rscript './src/nested/new-implementation-RNA/nested-pipeline.R' -s 1000 -S "lasso" -m "rf" -B 100 -P 0 -k 25 -f 1.2 -c "met.v.hc"


# ======== FSM=mrmr, Model=Rf =====================


Rscript './src/nested/new-implementation-RNA/nested-pipeline.R' -s 1000 -S "mrmr" -m "rf" -B 100 -P 0 -k 25 -f 1.2 -c "gbm.v.met"
Rscript './src/nested/new-implementation-RNA/nested-pipeline.R' -s 1000 -S "mrmr" -m "rf" -B 100 -P 0 -k 25 -f 1.2 -c "gbm.v.hc"
Rscript './src/nested/new-implementation-RNA/nested-pipeline.R' -s 1000 -S "mrmr" -m "rf" -B 100 -P 0 -k 25 -f 1.2 -c "met.v.hc"

# ======== FSM=t-test, Model=Rf =====================


Rscript './src/nested/new-implementation-RNA/nested-pipeline.R' -s 1000 -m "rf" -S "t-test" -B 100 -P 0 -k 25 -f 1.2 -c "gbm.v.met"
Rscript './src/nested/new-implementation-RNA/nested-pipeline.R' -s 1000 -m "rf" -S "t-test" -B 100 -P 0 -k 25 -f 1.2 -c "gbm.v.hc"
Rscript './src/nested/new-implementation-RNA/nested-pipeline.R' -s 1000 -m "rf" -S "t-test" -B 100 -P 0 -k 25 -f 1.2 -c "met.v.hc"

# ======== FSM=rf_rfe, Model=Rf =====================

Rscript './src/nested/new-implementation-RNA/nested-pipeline.R' -s 1000 -m "rf" -S "rf-rfe" -B 100 -P 0 -k 25 -f 1.2 -c "gbm.v.met"
Rscript './src/nested/new-implementation-RNA/nested-pipeline.R' -s 1000 -m "rf" -S "rf-rfe" -B 100 -P 0 -k 25 -f 1.2 -c "gbm.v.hc"
Rscript './src/nested/new-implementation-RNA/nested-pipeline.R' -s 1000 -m "rf" -S "rf-rfe" -B 100 -P 0 -k 25 -f 1.2 -c "met.v.hc"

