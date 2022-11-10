#Phenome-wide association study of genetically predicted B vitamins and homocysteine biomarkers with multiple health and disease outcomes: analysis of the UK Biobank by Lijuan Wang et al.
#PheWAS analysis
library(PheWAS)
id.vocab.code.count <- data$id.vocab.code.count
genotypes <- data$genotypes
phenotypes <- createPhenotypes(id.vocab.code.count, aggregate.fun=sum)
dat <- inner_join(inner_join(phenotypes,genotypes))
results <- phewas(dat, phenotypes=names(phenotypes)[-1], genotypes, covariates, cores=1)

#MR
library("TwoSampleMR")
exp_dat<-read_exposure_data(
  filename= "exposure.txt",
  sep = "\t",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  eaf_col = "eaf",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele"
)

outcome_dat <- read_outcome_data(
  snps = exp_dat$SNP,
  filename = "outcome.txt",
  sep = "\t",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  eaf_col = "eaf",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
)

dat <- harmonise_data(exp_dat, outcome_dat)
mr(dat)
mr_heterogeneity(dat)
mr_pleiotropy_test(dat)

#dose-response analysis
library(rms)
describe(data)
dd <- datadist(data)
options(datadist='dd')
fit <-lrm(outcome~rcs(exposure,5)+covariates, data=data)
anova(fit)
fit <- update(fit)
anova(fit)
OR <- Predict(fit, exposure, fun=exp, ref.zero=TRUE)

#mediation analysis
library("mediation")
b<-lm(outcome~mediator+covariates, data=data)
c<-glm(outcome~exposure+mediator+covariates,data=data, family = binomial("probit"))
contcont<-mediate(b,c,sims=100,treat="exposure",mediator="mediator")
summary(contcont)

#bioinformatics analysis was conducted based on the pathway analysis panel implemented in the STRING database.
