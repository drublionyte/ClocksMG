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

#zhang <- function(betas) {
 # coefs <- setNames(Zhang_CpGs$CoefficientTraining, Zhang_CpGs$CpGmarker)
  #CpGs <- intersect(names(coefs), rownames(betas))

  #betas <- betas[CpGs, ]
  #coefs <- coefs[CpGs]

  #tt <- betas * coefs
  #colSums(tt,na.rm=TRUE)+65.79295
#}

#zhang(betas)

#brain <- function(betas) {
 # coefs <- setNames(Brain_CpGs$Coefficient, Brain_CpGs$IlluminaProbeID)
  #CpGs <- intersect(names(coefs), rownames(betas))

  #betas <- betas[CpGs, ]
  #coefs <- coefs[CpGs]

  #tt <- betas * coefs
  #anti.trafo(colSums(tt,na.rm=TRUE)+0.57768257)
#}
