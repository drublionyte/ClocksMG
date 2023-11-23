install.packages("ClocksMG", repos = "https://github.com/drublionyte/ClocksMG")

anti.trafo <- function(x, adult.age=20) {
  ifelse(
    x < 0,
    (1 + adult.age) * exp(x) - 1,
    (1 + adult.age) * x + adult.age
  )
}

horwath <- function(betas) {
  coefs <- setNames(Horvath1_CpGs$CoefficientTraining, Horvath1_CpGs$CpGmarker)
  CpGs <- intersect(names(coefs), rownames(betas))

  betas <- betas[CpGs, ]
  coefs <- coefs[CpGs]

  tt <- betas * coefs
  anti.trafo(colSums(tt,na.rm=TRUE)+0.696)
}

zhang_en_2019 <- function(betas) {
  betas.norm <- apply(betas,1,scale)
  rownames(betas.norm) <- colnames(betas)

  coefs <- setNames(encoef$CoefficientTraining, encoef$CpGmarker)
  CpGs <- intersect(names(coefs), colnames(betas))

  tt <- betas.norm[CpGs, ] * coefs[CpGs]
  colSums(tt,na.rm=TRUE)+65.79295
}

zhang_blup_2019 <- function(betas) {
  betas.norm <- apply(betas,1,scale)
  rownames(betas.norm) <- colnames(betas)

  coefs <- setNames(blupcoef$CoefficientTraining, blupcoef$CpGmarker)
  CpGs <- intersect(names(coefs), colnames(betas))

  tt <- betas.norm[CpGs, ] * coefs[CpGs]
  colSums(tt,na.rm=TRUE)+91.15396
}

Wu <- function(betas) {
  coefs <- setNames(Horvath1_CpGs$CoefficientTraining, Horvath1_CpGs$CpGmarker)
  CpGs <- intersect(names(coefs), rownames(betas))

  betas <- betas[CpGs, ]
  coefs <- coefs[CpGs]

  tt <- betas * coefs
  anti.trafo(colSums(tt,na.rm=TRUE)+2.376853787, adult.age = 48)/12
}
