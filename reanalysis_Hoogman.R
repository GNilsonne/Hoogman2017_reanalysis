# By Gustav Nilsonne 2017-02-25
# Free to use with attribution

# Load packages
require(metafor)

# Read data
s_table_2 <- read.csv2("C:/Users/Gustav Nilsonne/Box Sync/Gustavs_arbete/Pek/2017_LPsych_letter_Hoogman/s_table_2.csv", skip = 1)
s_table_4 <- read.csv2("C:/Users/Gustav Nilsonne/Box Sync/Gustavs_arbete/Pek/2017_LPsych_letter_Hoogman/s_table_4.csv", skip = 1)

names(s_table_2) <- c(
  "Cohort",
  "ADHD symptom score in cases only",
  "N Lifetime co-morbid disorder present",
  "N Lifetime co-morbid depression present",
  "N Lifetime co-morbid anxiety",
  "N Lifetime co-morbid SUD",
  "IQ")
names(s_table_4) <- c(
  "Site",
  "Dx",
  "Accumbens",
  "Accumbens_SD",
  "Amygdala",
  "Amygdala_SD",
  "Caudate",
  "Caudate_SD",
  "Hippocampus",
  "Hippocampus_SD",
  "Pallidum",
  "Pallidum_SD",
  "Hippocampus",
  "Hippocampus_SD",
  "Thalamus",
  "Thalamus_SD",
  "ICV",
  "ICV_SD")

# Change IQ data to numeric and merge into data frame with regional brain volumes
IQ <- data.frame(strsplit(as.character(s_table_2$IQ), "[^0-9]+"))
IQ2 <- data.frame(t(IQ))
s_table_4$IQ <- as.integer(substr(s_table_2$IQ, 1, 3))
s_table_4$IQ_SD <- as.integer(as.character(IQ2$X2))
s_table_4$n <- as.integer(as.character(IQ2$X3))

# Test plots to verify that data are accessible
plot(Accumbens ~ IQ, data = s_table_4[s_table_4$Dx == "Controls", ], col = "blue", cex = sqrt(s_table_4[s_table_4$Dx == "Controls", ]$n/pi), frame.plot = F, xlim = c(90, 120), ylim = c(400, 900), main = "Accumbens", ylab = "Volume (mm3)")
points(Accumbens ~ IQ, data = s_table_4[s_table_4$Dx == "Cases", ], col = "red", cex = sqrt(s_table_4[s_table_4$Dx == "Cases", ]$n/pi))
legend("bottomleft", col = c("blue", "red"), legend = c("Controls", "Cases"), pch = 1)

plot(ICV ~ IQ, data = s_table_4[s_table_4$Dx == "Controls", ], col = "blue", cex = sqrt(s_table_4[s_table_4$Dx == "Controls", ]$n/pi), frame.plot = F, xlim = c(90, 120), ylim = c(1200000, 1800000), main = "ICV", ylab = "Volume (mm3)")
points(ICV ~ IQ, data = s_table_4[s_table_4$Dx == "Cases", ], col = "red", cex = sqrt(s_table_4[s_table_4$Dx == "Cases", ]$n/pi))
legend("bottomleft", col = c("blue", "red"), legend = c("Controls", "Cases"), pch = 1)

# Perform meta-regression by brain region
# THis involves making a new data frame for each region in the format suitable for the metafor package
accumbens <- data.frame(study = c(1:23),
                  source = s_table_4$Site[s_table_4$Dx == "Cases"],
                  n1i = s_table_4$n[s_table_4$Dx == "Cases"],
                  m1i = s_table_4$Accumbens[s_table_4$Dx == "Cases"],
                  sd1i = s_table_4$Accumbens_SD[s_table_4$Dx == "Cases"],
                  n2i = s_table_4$n[s_table_4$Dx == "Controls"],
                  m2i = s_table_4$Accumbens[s_table_4$Dx == "Controls"],
                  sd2i = s_table_4$Accumbens_SD[s_table_4$Dx == "Controls"],
                  delta_IQ = s_table_4$IQ[s_table_4$Dx == "Controls"] - s_table_4$IQ[s_table_4$Dx == "Cases"])
accumbens <- escalc(measure="SMD", m1i=m1i, sd1i=sd1i, n1i=n1i, m2i=m2i, sd2i=sd2i, n2i=n2i, data=accumbens)
res_accumbens <- rma(accumbens, slab = source)
forest(res_accumbens, main = "Accumbens, without meta-regression")
res_accumbens
res_accumbens_IQ <- rma(accumbens, mods = ~ delta_IQ, slab = source)
forest(res_accumbens_IQ, main = "Accumbens, with meta-regression")
res_accumbens_IQ

IQ_range <- -18:2
preds <- predict(res_accumbens_IQ, newmods = IQ_range)
plot(accumbens$delta_IQ, accumbens$yi, cex=accumbens$vi*60, frame.plot = F, xlab = "delta IQ", ylab = "SMD", main = "Accumbens, meta-regression")
lines(IQ_range, preds$pred)
lines(IQ_range, preds$ci.lb, lty="dashed")
lines(IQ_range, preds$ci.ub, lty="dashed")


amygdala <- data.frame(study = c(1:23),
                        source = s_table_4$Site[s_table_4$Dx == "Cases"],
                        n1i = s_table_4$n[s_table_4$Dx == "Cases"],
                        m1i = s_table_4$Amygdala[s_table_4$Dx == "Cases"],
                        sd1i = s_table_4$Amygdala_SD[s_table_4$Dx == "Cases"],
                        n2i = s_table_4$n[s_table_4$Dx == "Controls"],
                        m2i = s_table_4$Amygdala[s_table_4$Dx == "Controls"],
                        sd2i = s_table_4$Amygdala_SD[s_table_4$Dx == "Controls"],
                        delta_IQ = s_table_4$IQ[s_table_4$Dx == "Controls"] - s_table_4$IQ[s_table_4$Dx == "Cases"])
amygdala <- escalc(measure="SMD", m1i=m1i, sd1i=sd1i, n1i=n1i, m2i=m2i, sd2i=sd2i, n2i=n2i, data=amygdala)
res_amygdala <- rma(amygdala, slab = source)
forest(res_amygdala, main = "amygdala, without meta-regression")
res_amygdala
res_amygdala_IQ <- rma(amygdala, mods = ~ delta_IQ, slab = source)
forest(res_amygdala_IQ, main = "amygdala, with meta-regression")
res_amygdala_IQ


caudate <- data.frame(study = c(1:23),
                       source = s_table_4$Site[s_table_4$Dx == "Cases"],
                       n1i = s_table_4$n[s_table_4$Dx == "Cases"],
                       m1i = s_table_4$Caudate[s_table_4$Dx == "Cases"],
                       sd1i = s_table_4$Caudate_SD[s_table_4$Dx == "Cases"],
                       n2i = s_table_4$n[s_table_4$Dx == "Controls"],
                       m2i = s_table_4$Caudate[s_table_4$Dx == "Controls"],
                       sd2i = s_table_4$Caudate_SD[s_table_4$Dx == "Controls"],
                       delta_IQ = s_table_4$IQ[s_table_4$Dx == "Controls"] - s_table_4$IQ[s_table_4$Dx == "Cases"])
caudate <- escalc(measure="SMD", m1i=m1i, sd1i=sd1i, n1i=n1i, m2i=m2i, sd2i=sd2i, n2i=n2i, data=caudate)
res_caudate <- rma(caudate, slab = source)
forest(res_caudate, main = "caudate, without meta-regression")
res_caudate
res_caudate_IQ <- rma(caudate, mods = ~ delta_IQ, slab = source)
forest(res_caudate_IQ, main = "caudate, with meta-regression")
res_caudate_IQ


hippocampus <- data.frame(study = c(1:23),
                      source = s_table_4$Site[s_table_4$Dx == "Cases"],
                      n1i = s_table_4$n[s_table_4$Dx == "Cases"],
                      m1i = s_table_4$Hippocampus[s_table_4$Dx == "Cases"],
                      sd1i = s_table_4$Hippocampus_SD[s_table_4$Dx == "Cases"],
                      n2i = s_table_4$n[s_table_4$Dx == "Controls"],
                      m2i = s_table_4$Hippocampus[s_table_4$Dx == "Controls"],
                      sd2i = s_table_4$Hippocampus_SD[s_table_4$Dx == "Controls"],
                      delta_IQ = s_table_4$IQ[s_table_4$Dx == "Controls"] - s_table_4$IQ[s_table_4$Dx == "Cases"])
hippocampus <- escalc(measure="SMD", m1i=m1i, sd1i=sd1i, n1i=n1i, m2i=m2i, sd2i=sd2i, n2i=n2i, data=hippocampus)
res_hippocampus <- rma(hippocampus, slab = source)
forest(res_hippocampus, main = "hippocampus, without meta-regression")
res_hippocampus
res_hippocampus_IQ <- rma(hippocampus, mods = ~ delta_IQ, slab = source)
forest(res_hippocampus_IQ, main = "hippocampus, with meta-regression")
res_hippocampus_IQ


pallidum <- data.frame(study = c(1:23),
                          source = s_table_4$Site[s_table_4$Dx == "Cases"],
                          n1i = s_table_4$n[s_table_4$Dx == "Cases"],
                          m1i = s_table_4$Pallidum[s_table_4$Dx == "Cases"],
                          sd1i = s_table_4$Pallidum_SD[s_table_4$Dx == "Cases"],
                          n2i = s_table_4$n[s_table_4$Dx == "Controls"],
                          m2i = s_table_4$Pallidum[s_table_4$Dx == "Controls"],
                          sd2i = s_table_4$Pallidum_SD[s_table_4$Dx == "Controls"],
                          delta_IQ = s_table_4$IQ[s_table_4$Dx == "Controls"] - s_table_4$IQ[s_table_4$Dx == "Cases"])
pallidum <- escalc(measure="SMD", m1i=m1i, sd1i=sd1i, n1i=n1i, m2i=m2i, sd2i=sd2i, n2i=n2i, data=pallidum)
res_pallidum <- rma(pallidum, slab = source)
forest(res_pallidum, main = "pallidum, without meta-regression")
res_pallidum
res_pallidum_IQ <- rma(pallidum, mods = ~ delta_IQ, slab = source)
forest(res_pallidum_IQ, main = "pallidum, with meta-regression")
res_pallidum_IQ


putamen <- data.frame(study = c(1:23),
                          source = s_table_4$Site[s_table_4$Dx == "Cases"],
                          n1i = s_table_4$n[s_table_4$Dx == "Cases"],
                          m1i = s_table_4$Putamen[s_table_4$Dx == "Cases"],
                          sd1i = s_table_4$Putamen_SD[s_table_4$Dx == "Cases"],
                          n2i = s_table_4$n[s_table_4$Dx == "Controls"],
                          m2i = s_table_4$Putamen[s_table_4$Dx == "Controls"],
                          sd2i = s_table_4$Putamen_SD[s_table_4$Dx == "Controls"],
                          delta_IQ = s_table_4$IQ[s_table_4$Dx == "Controls"] - s_table_4$IQ[s_table_4$Dx == "Cases"])
putamen <- escalc(measure="SMD", m1i=m1i, sd1i=sd1i, n1i=n1i, m2i=m2i, sd2i=sd2i, n2i=n2i, data=putamen)
res_putamen <- rma(putamen, slab = source)
forest(res_putamen, main = "putamen, without meta-regression")
res_putamen
res_putamen_IQ <- rma(putamen, mods = ~ delta_IQ, slab = source)
forest(res_putamen_IQ, main = "putamen, with meta-regression")
res_putamen_IQ


thalamus <- data.frame(study = c(1:23),
                          source = s_table_4$Site[s_table_4$Dx == "Cases"],
                          n1i = s_table_4$n[s_table_4$Dx == "Cases"],
                          m1i = s_table_4$Thalamus[s_table_4$Dx == "Cases"],
                          sd1i = s_table_4$Thalamus_SD[s_table_4$Dx == "Cases"],
                          n2i = s_table_4$n[s_table_4$Dx == "Controls"],
                          m2i = s_table_4$Thalamus[s_table_4$Dx == "Controls"],
                          sd2i = s_table_4$Thalamus_SD[s_table_4$Dx == "Controls"],
                          delta_IQ = s_table_4$IQ[s_table_4$Dx == "Controls"] - s_table_4$IQ[s_table_4$Dx == "Cases"])
thalamus <- escalc(measure="SMD", m1i=m1i, sd1i=sd1i, n1i=n1i, m2i=m2i, sd2i=sd2i, n2i=n2i, data=thalamus)
res_thalamus <- rma(thalamus, slab = source)
forest(res_thalamus, main = "thalamus, without meta-regression")
res_thalamus
res_thalamus_IQ <- rma(thalamus, mods = ~ delta_IQ, slab = source)
forest(res_thalamus_IQ, main = "thalamus, with meta-regression")
res_thalamus_IQ

# Write results
out <- data.frame(Area = c("Accumbens", "Amygdala", "Caudate", "Pallidum", "Putamen", "Thalamus"),
                  Unadjusted_SMD_confint = c(paste(round(res_accumbens$b, 2), "[", round(res_accumbens$ci.lb, 2), round(res_accumbens$ci.ub, 2), "]"),
                                 paste(round(res_amygdala$b, 2), "[", round(res_amygdala$ci.lb, 2), round(res_amygdala$ci.ub, 2), "]"),
                                 paste(round(res_caudate$b, 2), "[", round(res_caudate$ci.lb, 2), round(res_caudate$ci.ub, 2), "]"),
                                 paste(round(res_pallidum$b, 2), "[", round(res_pallidum$ci.lb, 2), round(res_pallidum$ci.ub, 2), "]"),
                                 paste(round(res_putamen$b, 2), "[", round(res_putamen$ci.lb, 2), round(res_putamen$ci.ub, 2), "]"),
                                 paste(round(res_thalamus$b, 2), "[", round(res_thalamus$ci.lb, 2), round(res_thalamus$ci.ub, 2), "]")),
                  Adjusted_SMD_confint = c(paste(round(res_accumbens_IQ$b[1], 2), "[", round(res_accumbens_IQ$ci.lb[1], 2), round(res_accumbens_IQ$ci.ub[1], 2), "]"),
                                           paste(round(res_amygdala_IQ$b[1], 2), "[", round(res_amygdala_IQ$ci.lb[1], 2), round(res_amygdala_IQ$ci.ub[1], 2), "]"),
                                           paste(round(res_caudate_IQ$b[1], 2), "[", round(res_caudate_IQ$ci.lb[1], 2), round(res_caudate_IQ$ci.ub[1], 2), "]"),
                                           paste(round(res_pallidum_IQ$b[1], 2), "[", round(res_pallidum_IQ$ci.lb[1], 2), round(res_pallidum_IQ$ci.ub[1], 2), "]"),
                                           paste(round(res_putamen_IQ$b[1], 2), "[", round(res_putamen_IQ$ci.lb[1], 2), round(res_putamen_IQ$ci.ub[1], 2), "]"),
                                           paste(round(res_thalamus_IQ$b[1], 2), "[", round(res_thalamus_IQ$ci.lb[1], 2), round(res_thalamus_IQ$ci.ub[1], 2), "]")),
                  Delta_IQ_confint = c(paste(round(res_accumbens_IQ$b[2], 2), "[", round(res_accumbens_IQ$ci.lb[2], 2), round(res_accumbens_IQ$ci.ub[2], 2), "]"),
                                       paste(round(res_amygdala_IQ$b[2], 2), "[", round(res_amygdala_IQ$ci.lb[2], 2), round(res_amygdala_IQ$ci.ub[2], 2), "]"),
                                       paste(round(res_caudate_IQ$b[2], 2), "[", round(res_caudate_IQ$ci.lb[2], 2), round(res_caudate_IQ$ci.ub[2], 2), "]"),
                                       paste(round(res_pallidum_IQ$b[2], 2), "[", round(res_pallidum_IQ$ci.lb[2], 2), round(res_pallidum_IQ$ci.ub[2], 2), "]"),
                                       paste(round(res_putamen_IQ$b[2], 2), "[", round(res_putamen_IQ$ci.lb[2], 2), round(res_putamen_IQ$ci.ub[2], 2), "]"),
                                       paste(round(res_thalamus_IQ$b[2], 2), "[", round(res_thalamus_IQ$ci.lb[2], 2), round(res_thalamus_IQ$ci.ub[2], 2), "]")))

write.table(out, file = "IQ_correction.txt", sep = "\t", row.names = F)
                  