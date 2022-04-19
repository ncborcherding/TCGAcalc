library(survival)
library(dplyr)

TCGAcalc.Data <- readRDS("TCGAcalc.Data.rds")

cancer.types.available <- c("BLCA", "BRCA", "CESC", "COAD", "ESCA", "GBM",  "HNSC", "KIRC", 
                            "KIRP", "LAML", "LGG", "LIHC", "LUAD", "LUSC", "OV",   "PAAD", 
                            "PCPG", "PRAD", "READ", "SARC", "SKCM", "STAD", "TCGT", "THCA", 
                            "THYM", "UCEC")

checkGene <- function(cancer.types = NULL, 
                      gene = NULL,
                      survival.type = "overall",
                      cut.off = "median") {
dat <- TCGAcalc.Data
survival.columns <- survival.parser(survival.type)
gene.values <- lapply(seq_len(length(dat[[1]])), function(x) {
    y <- data.frame(t(dat[[1]][[x]][gene,]))
    y$cancer.type <- names(dat[[1]])[x]
    y
  })
gene.values <- dplyr::bind_rows(gene.values) 

survival.values <- lapply(seq_len(length(dat[[2]])), function(x) {
  y <- data.frame(dat[[2]][[x]][,survival.columns])
  y[,1] <- suppressWarnings(as.integer(y[,1]))
  y[,2] <- suppressWarnings(as.numeric(y[,2]))
  y
})
survival.values <- dplyr::bind_rows(survival.values) 
joined <- merge(gene.values, survival.values, by = 0)
joined <- split(joined, joined$cancer.type)
cancers.calculated <- names(joined)

survival.stats <- lapply(seq_len(length(joined)), function(x) {
  calculate.survival(joined[[x]], cut.off, gene, x)
})
results <- dplyr::bind_rows(survival.stats) 

plot <- ggplot(results, aes(x = HR, y = -log10(p))) + 
  geom_point(shape= 21, aes(fill = HR)) +
  geom_hline(yintercept = -log10(0.05), lty = 2) + 
  geom_vline(xintercept = 1) + 
  geom_text_repel(data = subset(results, p < 0.05), aes(label = cancer.type)) + 
  theme_classic() + 
  scale_fill_distiller(palette = "Spectral") + 
  guides(fill = "none")
return(plot)
  
}

survival.parser <- function(survival.type) {
  if(tolower(survival.type) %in%  c("overall", "os")) {
    col.surv <- c("OS", "OS.time")
  } else if(tolower(survival.type) %in% c("progresion-free", "progression free", "pfi")) {
    col.surv <- c("PFI", "PFI.time")
  } else if(tolower(survival.type) %in% c("disease-free", "disease free", "dfi")) {
    col.surv <- c("DFI", "DFI.time")
  } else if(tolower(survival.type) %in% c("disease-specific", "disease specific", "dss")) {
    col.surv <- c("DSS", "DSS.time")
  }
  return(col.surv)
}

calculate.survival <- function(data.prime, cut.off, gene, x) {
  slice <- translate.cut.off(cut.off)
  if(slice %in% c(2,3,4,5)) {
    data.prime <- data.prime %>%
                  mutate(group = ntile(data.prime[,gene], n = slice))
    if(slice %in% c(2,3,4)) {
      new.levels <- relevel.cut.off(slice)
      data.prime[,"group"] <- unname(new.levels[data.prime[,"group"]])
    }
    OS.calc <- survfit(Surv(data.prime[,survival.columns[2]], 
                            data.prime[,survival.columns[1]]) ~ data.prime[,"group"])
    cox.calc<- coxph(Surv(data.prime[,survival.columns[2]], 
                          data.prime[,survival.columns[1]]) ~ data.prime[,"group"])
  } else if(slice == "optimum") {
    opt.cut <- survminer::surv_cutpoint(data.prime, survival.columns[2], survival.columns[1], 
                             variables = gene, minprop = 0.15)
    opt.cat <- survminer::surv_categorize(opt.cut)
    OS.calc <- survfit(Surv(opt.cat[,1], 
                            opt.cat[,2]) ~ opt.cat[,3])
    cox.calc<- coxph(Surv(opt.cat[,1], 
                          opt.cat[,2]) ~ opt.cat[,3])
  }
  res.summary <- summary(cox.calc)
  cox.HR = 1/(res.summary$coefficients[,2])      #HR
  p.value <- res.summary$coefficients[,5]
  summary <- data.frame(HR = cox.HR, 
                        p = p.value, 
                        comparator = str_split(row.names(res.summary$coefficients), "]", simplify = TRUE)[,2], 
                        cancer.type = cancers.calculated[x])
  summary 
}

translate.cut.off <- function(cut.off) {
  list.cut.off <- list(tertile = 3, quartile =4, quintile = 5, median = 2, optimum = "optimum")
  slice <- unname(unlist(list.cut.off[cut.off]))
  return(slice)
}

relevel.cut.off <- function(slice) {
  likert <- list(c("High", "Low"), 
              c("High", "Med", "Low"), 
              c("High", "Med-High", "Med-Low", "Low"))
  relevels <- unlist(likert[slice-1])
  names(relevels) <- rev(seq_len(slice))
  return(relevels)
}

