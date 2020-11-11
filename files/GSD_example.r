# Gallium numbers
alphaGallium <- 0.05
betaGallium <- 0.2
m1 <- 12
hr <- 0.75
m2 <- m1 / hr

# Required events without interim
nevent <- getSampleSizeSurvival(lambda1 = getLambdaByMedian(m2), lambda2 = getLambdaByMedian(m1), 
                                sided = 2, alpha = alphaGallium, beta = betaGallium)
nevent <- ceiling(nevent$maxNumberOfEvents)

# Events with interim (as per Gallium protocol)
nevents_i1 <- c(111, 248, 370)


# ==============================================
# slides1.rnw
# ==============================================

## Replicate Gallium design with interims as closely as possible. 
## Two tricks are applied to do this:
## 1) Designs with futility are always one-sided in rpact. To mimic the "two-sided" Gallium design,
##    a futility boundary at the efficacy interim is added which would stop the trial if it went 
##    "significantly" in the wrong direction. This should get stopping probabilities correct.
## 2) Beta is slightly increased (i.e. power slightly decreased from 80%) to replicate the exactly 
##    370 events from Gallium. [For exactly 80% power ]

# OBF design without interims (to get futility bound at second interim as per 1) above 
designOBF <- getDesignGroupSequential(informationRates = nevents_i1 / max(nevents_i1),
                                      typeOfDesign = "asOF", sided = 1, 
                                      alpha = alphaGallium / 2)

# Gallium GS design 
# PS: slight beta increase as per 2) above calculated using getPowerSurvival
designGallium <- getDesignGroupSequential(informationRates = nevents_i1 / max(nevents_i1),
                                          typeOfDesign = "asOF", sided = 1, 
                                          alpha = alphaGallium / 2, 
                                          beta = betaGallium + 0.0050968,
                                          futilityBounds = c(0, -designOBF$criticalValues[2]), 
                                          bindingFutility = FALSE)

samplesizeGallium <- getSampleSizeSurvival(design = designGallium,
                                           lambda2 = log(2) / m1, hazardRatio = hr,
                                           dropoutRate1 = 0.025, dropoutRate2 = 0.025, dropoutTime = 12,
                                           accrualTime = 0:6, accrualIntensity = seq(6, 42, by = 6),
                                           maxNumberOfSubjects = 1200)

# Local significance levels (two-sided)
alphas_i1 <- designGallium$stageLevels * 2

# compute MDD 
hrMDD_i1 <- as.vector(samplesizeGallium$criticalValuesEffectScale)

# stopping probabilities at futility and efficacy interim under H0 and H1
designChar <- getDesignCharacteristics(designGallium)
stopProbsH0 <- getPowerAndAverageSampleNumber(designGallium, theta = 0, nMax = designChar$shift)
stopProbsH1 <- getPowerAndAverageSampleNumber(designGallium, theta = 1, nMax = designChar$shift)

stopFutIA_H0 <- stopProbsH0$earlyStop["stage = 1", 1]
stopEffIA_H0 <- stopProbsH0$earlyStop["stage = 2", 1]

stopFutIA_H1 <- stopProbsH1$earlyStop["stage = 1", 1]
stopEffIA_H1 <- stopProbsH1$earlyStop["stage = 2", 1]

# Expected number of events under H0 and H1 
expH0 <- samplesizeGallium$expectedEventsH0

# Alternative calculation (gives same result):
# expH0 <- nevents_i1[1] * stopFutIA_H0+nevents_i1[2] * stopEffIA_H0 + nevents_i1[3] * (1 - stopFutIA_H0-stopEffIA_H0)

expH1 <- samplesizeGallium$expectedEventsH1

# ==============================================
# for simulations in slides2.rnw
# ==============================================
loghr <- log(hr)
sd <- sqrt(4)
fut_hr <- 1
eff_hr <- hrMDD_i1[2]
final_hr <- hrMDD_i1[3]

