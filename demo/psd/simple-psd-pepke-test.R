source("~/projects/mlm/KaSimTutorial/analyse.R")
dat <- run.kappa("simple-psd-pepke.ka", t=1000)

vol <- 0.5 # um3
## Concentration of 1 molecule in uM
yscale <- 1E18/6.02E23/vol*1000

par(mfrow=c(2,2))
plot(dat[,c("time",  "CaCaM1N",   "CaCaM2N",   "CaCaM1C",   "CaCaM2C",   "CaCaM4")], xlab="Time (ms)", ylab="conc (uM)", yscale=yscale)
plot(dat[,c("time", "KCaCaM1N",  "KCaCaM2N",  "KCaCaM1C",  "KCaCaM2C",  "KCaCaM4")], xlab="Time (ms)", ylab="conc (uM)", yscale=yscale)
plot(dat[,c("time", "KCaCaM1Np", "KCaCaM2Np", "KCaCaM1Cp", "KCaCaM2Cp", "KCaCaM4p")], xlab="Time (ms)", ylab="conc (uM)", yscale=yscale)
plot(dat[,c("time", "CaMKIIp")], xlab="Time (ms)", ylab="conc (uM)", yscale=yscale)
dev.print(pdf, "simple-psd-pepke.pdf")
