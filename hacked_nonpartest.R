#added the progress bar...

hackNonParTest<-function (formula, data, permtest = TRUE, permreps = 10000, plots = TRUE, 
          tests = c(1, 1, 1, 1), releffects = TRUE) 
{
  if (!is(formula, "formula")) {
    return("Error: Please give a formula")
  }
  formula = Formula(formula)
  frame = model.frame(formula, data = data)
  groupvar.location = length(frame[1, ])
  groupvar = names(frame)[groupvar.location]
  vars = names(frame)[1:(groupvar.location - 1)]
  if (sum(is.na(frame)) > 0) {
    return("Error: Missing Data")
  }
  if (plots) {
    if (length(vars) > 1) {
      oask <- devAskNewPage(TRUE)
      on.exit(devAskNewPage(oask))
    }
    for (i in 1:length(vars)) {
      boxplot(frame[, vars[i]] ~ frame[, groupvar], data = frame, 
              main = vars[i], ylab = vars[i], xlab = names(frame)[1])
    }
  }
  if (!is.factor(frame[, groupvar])) {
    frame[, groupvar] <- factor(frame[, groupvar])
  }
  o <- order(frame[, groupvar])
  frame <- frame[o, ]
  p <- length(vars)
  a <- length(levels(frame[, groupvar]))
  N <- length(frame[, 1])
  ssize <- array(NA, a)
  lims <- matrix(NA, 2, a)
  for (i in 1:a) {
    ssize[i] <- length(frame[frame[, groupvar] == levels(frame[, 
                                                               groupvar])[i], 1])
    lims[1, i] <- min(which(frame[, groupvar] == levels(frame[, 
                                                              groupvar])[i]))
    lims[2, i] <- max(which(frame[, groupvar] == levels(frame[, 
                                                              groupvar])[i]))
  }
  if (sum(ssize < 2) > 0) {
    return("Error: Each group must have sample size of at least 2")
  }
  Rmat <- matrix(NA, N, p)
  for (j in 1:p) {
    Rmat[, j] <- rank(frame[, vars[j]], ties.method = "average")
  }
  Rbars <- matrix(NA, a, p)
  for (i in 1:a) {
    for (j in 1:p) {
      Rbars[i, j] <- mean(Rmat[(lims[1, i]:lims[2, i]), 
                               j])
    }
  }
  Rtilda <- (1/a) * colSums(Rbars)
  Rbarovr <- (1/N) * colSums(Rmat)
  H1 <- matrix(0, p, p)
  H2 <- matrix(0, p, p)
  G1 <- matrix(0, p, p)
  G2 <- matrix(0, p, p)
  G3 <- matrix(0, p, p)
  for (i in 1:a) {
    H1 <- H1 + ssize[i] * (Rbars[i, ] - Rbarovr) %*% t(Rbars[i, 
                                                             ] - Rbarovr)
    H2 <- H2 + (Rbars[i, ] - Rtilda) %*% t(Rbars[i, ] - Rtilda)
    for (j in 1:ssize[i]) {
      G1 <- G1 + (((Rmat[(lims[1, i] + j - 1), (1:p)]) - 
                     Rbars[i, ]) %*% t((Rmat[(lims[1, i] + j - 1), 
                                             (1:p)]) - Rbars[i, ]))
      G2 <- G2 + (1 - (ssize[i]/N)) * (1/(ssize[i] - 1)) * 
        (((Rmat[(lims[1, i] + j - 1), (1:p)]) - Rbars[i, 
                                                      ]) %*% t((Rmat[(lims[1, i] + j - 1), (1:p)]) - 
                                                                 Rbars[i, ]))
      G3 <- G3 + (1/(ssize[i] * (ssize[i] - 1))) * (((Rmat[(lims[1, 
                                                                 i] + j - 1), (1:p)]) - Rbars[i, ]) %*% t((Rmat[(lims[1, 
                                                                                                                      i] + j - 1), (1:p)]) - Rbars[i, ]))
    }
  }
  H1 <- (1/(a - 1)) * H1
  H2 <- (1/(a - 1)) * H2
  G1 <- (1/(N - a)) * G1
  G2 <- (1/(a - 1)) * G2
  G3 <- (1/a) * G3
  if (det(H1) == 0 | det(H2) == 0 | det(G1) == 0 | det(G2) == 
        0 | det(G3) == 0) {
    warning("Rank Matrix is Singular, only ANOVA test returned")
    tests = c(1, 0, 0, 0)
  }
  if (releffects) {
    rel <- as.data.frame(matrix(NA, a, p))
    names(rel) <- vars
    rownames(rel) <- levels(frame[, groupvar])
    for (i in 1:a) {
      for (j in 1:p) {
        rel[i, j] <- (1/N) * (mean(Rmat[which(frame[, 
                                                    groupvar] == levels(frame[, groupvar])[i]), 
                                        j]) - 0.5)
      }
    }
    if (a == 2) {
      origrel <- rel
      for (j in 1:p) {
        rel[1, j] <- origrel[1, j] - origrel[2, j] + 
          0.5
        rel[2, j] <- origrel[2, j] - origrel[1, j] + 
          0.5
      }
    }
  }
  base <- basenonpartest(frame, groupvar, vars, tests)
  pvalanova <- base$pvalanova
  pvalLH <- base$pvalLH
  FLH <- base$FLH
  pvalBNP <- base$pvalBNP
  FBNP <- base$FBNP
  Fanova <- base$Fanova
  FWL <- base$FWL
  pvalWL <- base$pvalWL
  if (permtest) {
    perms <- matrix(NA, permreps, 4)
    tempframe <- frame
    pb<-txtProgressBar(0,permreps,style=3)
    for (i in 1:permreps) {
      setTxtProgressBar(pb,i)
      tempframe[groupvar] <- sample(as.vector(t((tempframe[groupvar]))))
      permout <- basenonpartest(tempframe, groupvar, vars, 
                                tests)
      perms[i, 1] <- permout$Fanova
      perms[i, 2] <- permout$FLH
      perms[i, 3] <- permout$FBNP
      perms[i, 4] <- permout$FWL
    }
    pvalanovperm <- sum(as.numeric(perms[, 1] > Fanova))/permreps
    pvalLHperm <- sum(as.numeric(perms[, 2] > FLH))/permreps
    pvalBNPperm <- sum(as.numeric(perms[, 3] > FBNP))/permreps
    pvalWLperm <- sum(as.numeric(perms[, 4] > FWL))/permreps
    results = matrix(c(Fanova, FLH, FBNP, FWL, pvalanova, 
                       pvalLH, pvalBNP, pvalWL, pvalanovperm, pvalLHperm, 
                       pvalBNPperm, pvalWLperm), ncol = 3)
    results = data.frame(results, row.names = c("ANOVA type test p-value", 
                                                "McKeon approx. for the Lawley Hotelling Test", "Muller approx. for the Bartlett-Nanda-Pillai Test", 
                                                "Wilks Lambda"))
    colnames(results) = c("Test Statistic", "P-value", "Permutation Test p-value")
    results[, 3] <- round(as.numeric(results[, 3]), digits = 3)
  }
  else {
    results = matrix(c(Fanova, FLH, FBNP, FWL, pvalanova, 
                       pvalLH, pvalBNP, pvalWL), ncol = 2)
    results = data.frame(results, row.names = c("ANOVA type test p-value", 
                                                "McKeon approx. for the Lawley Hotelling Test", "Muller approx. for the Bartlett-Nanda-Pillai Test", 
                                                "Wilks Lambda"))
    colnames(results) = c("Test Statistic", "P-value")
  }
  results[, 1] <- round(as.numeric(results[, 1]), digits = 3)
  results[, 2] <- round(as.numeric(results[, 2]), digits = 3)
  if (releffects) {
    if (a == 2) {
      return(list(results = results, twogroupreleffects = rel))
    }
    else {
      return(list(results = results, releffects = rel))
    }
  }
  return(results)
}

environment(hackNonParTest)<-as.environment(package:npmv)