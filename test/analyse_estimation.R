library(abcrf)

# entry data
pos = read.table("/home/croux/Programmes/msms/ABCsweep/test/selection/rep1/outputABC_positions.txt", h=T)
param = read.table("/home/croux/Programmes/msms/ABCsweep/test/selection/rep1/outputABC_prior.txt", h=T)
stats = read.table("/home/croux/Programmes/msms/ABCsweep/test/selection/rep1/outputABC_sumStats.txt", h=T)

pos = rbind(pos, read.table("/home/croux/Programmes/msms/ABCsweep/test/selection/rep2/outputABC_positions.txt", h=T))
param = rbind(param, read.table("/home/croux/Programmes/msms/ABCsweep/test/selection/rep2/outputABC_prior.txt", h=T))
stats = rbind(stats, read.table("/home/croux/Programmes/msms/ABCsweep/test/selection/rep2/outputABC_sumStats.txt", h=T))
#
Sp = regAbcrf(param$Sp, stats, paral=T, ntree=2000)
SF = regAbcrf(param$SF, stats, paral=T, ntree=2000)
SAA = regAbcrf(param$SAA, stats, paral=T, ntree=2000)

#
obs_stats = read.table("/home/croux/Programmes/msms/ABCsweep/test/selection/rep1/obs/outputABC_sumStats.txt", h=T)
obs_param = read.table("/home/croux/Programmes/msms/ABCsweep/test/selection/rep1/obs/outputABC_prior.txt", h=T)

obs = cbind(obs_stats, obs_param)
obs = na.omit(obs)

Sp_prediction = predict(Sp, obs[, 1:133])
SF_prediction = predict(SF, obs[, 1:133])
SAA_prediction = predict(SAA, obs[, 1:133])

# plot 1
dev.new(width=12, height=5)
par(mfrow=c(1,3))
plot(obs$Sp, Sp_prediction$med, xlim=c(0,1), ylim=c(0,1), pch=16, col=pouet, xlab = "real values", ylab = "ABC estimates", main = "position of the\ntarget of selection\nalong a 100kb fragment");abline(a=0, b=1, col="red", lwd=3, lty=2)
tmp = loess(Sp_prediction$med~obs$Sp)
lines(0:100/100, predict(tmp, 0:100/100), col="red", lwd=3)

plot(obs$SF, SF_prediction$med, xlim=c(0,1), ylim=c(0,1), pch=16, col=pouet, xlab = "real values", ylab = "ABC estimates", main = "Fixation time\n(in 400,000 generations)");abline(a=0, b=1, col="red", lwd=3, lty=2)
tmp = loess(SF_prediction$med~obs$SF)
lines(0:100/100, predict(tmp, 0:100/100), col="red", lwd=3)

#plot(obs$SAA/200000, SAA_prediction$med/200000, xlim=c(0,10000/200000), ylim=c(0,10000/200000), pch=16, col=pouet, xlab = "real values", ylab = "ABC estimates", main = "s_AA");abline(a=0, b=1, col="red", lwd=2, lty=2)
#tmp = loess(SAA_prediction$med~obs$SAA)
#lines(0:1000000/200000/100, predict(tmp, 0:1000000/200000/100), col="red", lwd=2)

a=obs$SAA/200000
b=SAA_prediction$med/200000
plot(a, b, xlim=c(0,0.05), ylim=c(0,0.05), pch=16, col=pouet, xlab = "real values", ylab = "ABC estimates", main = "s_AA");abline(a=0, b=1, col="red", lwd=3, lty=2)
tmp = loess(b~a)
lines(0:5000/10000, predict(tmp, 0:5000/10000), col="red", lwd=3)


# plot 2 
dev.new(width=12, height=5)
par(mfrow=c(1,3))
# A
bin1 = which(obs$Sp < 0.2)
bin2 = which(obs$Sp >= 0.2 & obs$Sp < 0.4)
bin3 = which(obs$Sp >= 0.4 & obs$Sp < 0.6)
bin4 = which(obs$Sp >= 0.6 & obs$Sp < 0.8)
bin5 = which(obs$Sp >= 0.8 & obs$Sp <= 1.0)
tmp = Sp_prediction$med/obs$Sp
pouet = c("<0.2", "0.2-0.4", "0.4-0.6", "0.6-0.8", "0.8-1")
boxplot(tmp[bin1],tmp[bin2],tmp[bin3],tmp[bin4],tmp[bin5], outline=F, xlab = "position Sp", ylab= "ratio [estimated / real]", names=pouet, cex.lab=1.2); abline(h=1, col="red", lwd=2)

# B
bin1 = which(obs$SF < 0.2)
bin2 = which(obs$SF >= 0.2 & obs$SF < 0.4)
bin3 = which(obs$SF >= 0.4 & obs$SF < 0.6)
bin4 = which(obs$SF >= 0.6 & obs$SF < 0.8)
bin5 = which(obs$SF >= 0.8 & obs$SF <= 1.0)
tmp = SF_prediction$med/obs$SF
pouet = c("<0.2", "0.2-0.4", "0.4-0.6", "0.6-0.8", "0.8-1")
boxplot(tmp[bin1],tmp[bin2],tmp[bin3],tmp[bin4],tmp[bin5], outline=F, xlab = "fixation time SF", ylab= "ratio [estimated / real]", names=pouet, cex.lab=1.2); abline(h=1, col="red", lwd=2)

# C
a=obs$SAA/200000
b=SAA_prediction$med/200000
tmp=b/a
bin1 = which(a < 0.2*0.05)
bin2 = which(a >= 0.2*0.05 & a < 0.4*0.05)
bin3 = which(a >= 0.4*0.05 & a < 0.6*0.05)
bin4 = which(a >= 0.6*0.05 & a < 0.8*0.05)
bin5 = which(a >= 0.8*0.05 & a <= 1.0*0.05)
pouet = c("<0.01", "0.01-0.02", "0.02-0.03", "0.03-0.04", "0.04-05")
boxplot(tmp[bin1],tmp[bin2],tmp[bin3],tmp[bin4],tmp[bin5], outline=F, xlab = "coefficient SAA", ylab= "ratio [estimated / real]", names=pouet, cex.lab=1.2); abline(h=1, col="red", lwd=2)

