#library(ggplot2)
library(plotrix)
setwd("/home/roustique/Документы/Курсач/Прога")

Data1 <- read.table ("result.dat",sep=",", header=T)

#pdf(file = "Trajectory0.pdf")

plot(Data1$x, Data1$y, xlab = "x, cm", ylab = "y, cm", type="l", col="white", asp=1)

grid (NULL,NULL, lty = 6, col = "gray")
lines(0, 0, type="p")
lines(Data1$x, Data1$y, type="l")

draw.circle(0,0,3e+11,density=5,
            angle=45)

#dev.off()

#plot(Data1$x, Data1$y, xlab = "x, cm", ylab = "y, cm", type="l", asp=1, xlim=c(-5e+9, 5e+9), ylim=c(-5e+9, 5e+9))
#lines(0, 0, type="p")
#draw.circle(0,0,2.95e+9,density=5,
#            angle=45)

plot(Data1$E/50762.62-1, type="l", xlab="index", ylab=expression('E/E'[0]-1))
sd(Data1$E/50762.62-1)
plot(Data1$E, type="l")

#Data1["r"] <- sqrt(Data1$x^2+Data1$y^2)
#Data1["v"] <- sqrt(Data1$vx^2+Data1$vy^2)
