library(survivalROC)

## Define a function
fun.survivalROC <- function(lp, t) {
  res <- with(pbc,
              survivalROC(Stime        = time,
                          status       = event,
                          marker       = get(lp),
                          predict.time = t,
                          method       = "KM"))       # KM method without smoothing
  
  ## Plot ROCs
  with(res, plot(TP ~ FP, type = "l", main = sprintf("t = %.0f, AUC = %.2f", t, AUC)))
  abline(a = 0, b = 1, lty = 2)
  
  res
}

## 2 x 5 layout
layout(matrix(1:10, byrow = T, ncol = 5))

## Model with age and sex
res.survivalROC.age.sex <- lapply(1:10 * 365.25, function(t) {
  fun.survivalROC(lp = "lp.age.sex", t)
})