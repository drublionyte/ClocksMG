anti.trafo <- function(x, adult.age=20) {
  ifelse(
    x < 0,
    (1 + adult.age) * exp(x) - 1,
    (1 + adult.age) * x + adult.age
  )
}

epi_horvath <- function(betas) {
  data_path <- system.file("data", "Horvath_CpGs.coef", package = "ClocksMG")
  Horvath_CpGs <- read.table(data_path, stringsAsFactors = FALSE, header = TRUE)
  intercept <- Horvath_CpGs[1, 2]
  coefs <- setNames(Horvath_CpGs$CoefficientTraining, Horvath_CpGs$CpGmarker)
  CpGs <- intersect(names(coefs), rownames(betas))

  betas <- betas[CpGs, ]
  coefs <- coefs[CpGs]

  tt <- betas * coefs
  anti.trafo(colSums(tt,na.rm=TRUE)+intercept)
}

epi_horvath_snb <- function(betas) {
  data_path <- system.file("data", "Horvath_snb_CpGs.coef", package = "ClocksMG")
  Horvath_CpGs <- read.table(data_path, stringsAsFactors = FALSE, header = TRUE)
  intercept <- Horvath_CpGs[1, 2]
  coefs <- setNames(Horvath_CpGs$CoefficientTraining, Horvath_CpGs$CpGmarker)
  CpGs <- intersect(names(coefs), rownames(betas))

  betas <- betas[CpGs, ]
  coefs <- coefs[CpGs]

  tt <- betas * coefs
  anti.trafo(colSums(tt,na.rm=TRUE)+intercept)
}

epi_zhang_en <- function(betas) {
  data_path <- system.file("data", "Zhang_en_CpGs.coef", package = "ClocksMG")
  Zhang_en_CpGs <- read.table(data_path, header = TRUE)
  intercept <- Zhang_en_CpGs[1, 2]
  betas.norm <- apply(t(betas),1,scale)
  rownames(betas.norm) <- rownames(betas)

  coefs <- setNames(Zhang_en_CpGs$CoefficientTraining, Zhang_en_CpGs$CpGmarker)
  CpGs <- intersect(names(coefs), rownames(betas))

  tt <- betas.norm[CpGs, ] * coefs[CpGs]
  colSums(tt,na.rm=TRUE) + intercept
}

epi_zhang_blup <- function(betas) {
  data_path <- system.file("data", "Zhang_blup_CpGs.coef", package = "ClocksMG")
  Zhang_blup_CpGs <- read.table(data_path, header = TRUE)
  intercept <- Zhang_blup_CpGs[1, 2]
  betas.norm <- apply(t(betas),1,scale)
  rownames(betas.norm) <- colnames(betas)

  coefs <- setNames(Zhang_blup_CpGs$CoefficientTraining, Zhang_blup_CpGs$CpGmarker)
  CpGs <- intersect(names(coefs), rownames(betas))

  tt <- betas.norm[CpGs, ] * coefs[CpGs]
  colSums(tt,na.rm=TRUE) + intercept
}

epi_wu <- function(betas) {
  data_path <- system.file("data", "Wu_CpGs.coef", package = "ClocksMG")
  Wu_CpGs <- read.table(data_path, header = TRUE)
  intercept <- Wu_CpGs[1, 2]
  coefs <- setNames(Wu_CpGs$CoefficientTraining, Wu_CpGs$CpGmarker)
  CpGs <- intersect(names(coefs), rownames(betas))

  betas <- betas[CpGs, ]
  coefs <- coefs[CpGs]

  tt <- betas * coefs
  anti.trafo(colSums(tt,na.rm=TRUE) + intercept, adult.age = 48)/12
}

epi_bocklandt <- function(betas) {
  data_path <- system.file("data", "Bocklandt_CpGs.coef", package = "ClocksMG")
  Bocklandt_CpGs <- read.table(data_path, header = TRUE)
  coefs <- setNames(Bocklandt_CpGs$CoefficientTraining, Bocklandt_CpGs$CpGmarker)
  CpGs <- intersect(names(coefs), rownames(betas))

  betas <- betas[CpGs, ]
  coefs <- coefs[CpGs]

  tt <- betas * coefs
  colSums(tt, na.rm = TRUE)
}

epi_garagnani <- function(betas) {
  data_path <- system.file("data", "Garagnani_CpGs.coef", package = "ClocksMG")
  Garagnani_CpGs <- read.table(data_path, stringsAsFactors = FALSE, header = TRUE)
  coefs <- setNames(Garagnani_CpGs$CoefficientTraining, Garagnani_CpGs$CpGmarker)
  CpGs <- intersect(names(coefs), rownames(betas))

  betas <- betas[CpGs, ]
  coefs <- coefs[CpGs]

  tt <- betas * coefs
  colSums(tt, na.rm = TRUE)
}

epi_hannum <- function(betas) {
  data_path <- system.file("data", "Hannum_CpGs.coef", package = "ClocksMG")
  Hannum_CpGs <- read.table(data_path, header = TRUE)
  coefs <- setNames(Hannum_CpGs$CoefficientTraining, Hannum_CpGs$CpGmarker)
  CpGs <- intersect(names(coefs), rownames(betas))

  betas <- betas[CpGs, ]
  coefs <- coefs[CpGs]

  tt <- betas * coefs
  colSums(tt, na.rm = TRUE)
}

epi_lin <- function(betas) {
  data_path <- system.file("data", "Lin_CpGs.coef", package = "ClocksMG")
  Lin_CpGs <- read.table(data_path, header = TRUE)
  intercept <- Lin_CpGs[1, 2]
  coefs <- setNames(Lin_CpGs$CoefficientTraining, Lin_CpGs$CpGmarker)
  CpGs <- intersect(names(coefs), rownames(betas))

  betas <- betas[CpGs, ]
  coefs <- coefs[CpGs]

  tt <- betas * coefs
  colSums(tt, na.rm = TRUE) + intercept
}

epi_vidalbralo <- function(betas) {
  data_path <- system.file("data", "VidalBralo_CpGs.coef", package = "ClocksMG")
  VidalBralo_CpGs <- read.table(data_path, header = TRUE)
  intercept <- VidalBralo_CpGs[1, 2]
  coefs <- setNames(VidalBralo_CpGs$CoefficientTraining, VidalBralo_CpGs$CpGmarker)
  CpGs <- intersect(names(coefs), rownames(betas))

  betas <- betas[CpGs, ]
  coefs <- coefs[CpGs]

  tt <- betas * coefs
  colSums(tt, na.rm = TRUE) + intercept
}

epi_weidner <- function(betas) {
  data_path <- system.file("data", "Weidner_CpGs.coef", package = "ClocksMG")
  Weidner_CpGs <- read.table(data_path, header = TRUE)
  intercept <- Weidner_CpGs[1, 2]
  coefs <- setNames(Weidner_CpGs$CoefficientTraining, Weidner_CpGs$CpGmarker)
  CpGs <- intersect(names(coefs), rownames(betas))

  betas <- betas[CpGs, ]
  coefs <- coefs[CpGs]

  tt <- betas * coefs
  colSums(tt, na.rm = TRUE) + intercept
}
