# first two biased ES, var, last two unbiased

f.Q.NA = function(dadj, varadj) {
  w <- 1/varadj
  tmp1 <- w * dadj
  mu <- rowSums(tmp1, na.rm = TRUE)/rowSums(w, na.rm = TRUE)
  Q <- rowSums(w * (dadj - mu)^2, na.rm = TRUE)
}
