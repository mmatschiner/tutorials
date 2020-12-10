
# read in the results with 50 SNP windows and a step of 25 SNPs
bigStep <- read.table("mbuna_deep_Diplotaxodon_localFstats__50_25.txt",as.is=T,header=T)

# plot D in windows along scaffold 18:
plot(bigStep$windowStart/1000000, bigStep$D,type="l",xlab="scaffold 18 coordinate (Mb)",ylab="D (ABBA-BABA)")
# plot f_dM in windows along scaffold 18:
plot(bigStep$windowStart/1000000, bigStep$f_dM,type="l",xlab="scaffold 18 coordinate (Mb)",ylab="f_dM")
# plot f_d in windows along scaffold 18:
plot(bigStep$windowStart/1000000, bigStep$f_d,type="l",xlab="scaffold 18 coordinate (Mb)",ylab="f_d")
# use f_dM and zoom-in on the region of interest 
plot(bigStep$windowStart/1000000, bigStep$f_dM,type="l",xlim=c(4,4.5),xlab="scaffold 18 coordinate (Mb)",ylab="f_dM",ylim=c(-0.2,0.8))

# Same window: 50SNPs, and the smallest 1 SNP step
smallestStep <- read.table("mbuna_deep_Diplotaxodon_localFstats__50_1.txt",as.is=T,header=T)
# plot f_dM again:
plot(smallestStep$windowStart/1000000, smallestStep$f_dM,type="l",xlim=c(4,4.5),xlab="scaffold 18 coordinate (Mb)",ylab="f_dM",ylim=c(-0.2,0.8))

# Smaller 10SNP window, and the smallest 1 SNP step
smallerWindow <- read.table("mbuna_deep_Diplotaxodon_localFstats__10_1.txt",as.is=T,header=T)
# plot f_dM again:
plot(smallerWindow$windowStart/1000000, smallerWindow$f_dM,type="l",xlim=c(4.0,4.5),xlab="scaffold 18 coordinate (Mb)",ylab="f_dM",ylim=c(-0.2,0.9),pch=16)

# Finally, a tiny window with 2 SNPs and a step of 1 SNP
smallWindow <- read.table("mbuna_deep_Diplotaxodon_localFstats__2_1.txt",as.is=T,header=T)
plot(smallWindow$windowStart, smallWindow$f_dM,type="l",xlim=c(4000000,4500000),xlab="scaffold 18 coordinate (Mb)",ylab="f_dM",ylim=c(-0.2,0.8))
