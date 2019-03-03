#library(ggplot2)
library(plotrix)
setwd("/home/roustique/Документы/Курсач/Прога")

Data1 <- read.table ("result.dat",sep=",", header=T)
Data2 <- read.table ("resultN.dat", sep=",", header=T)

#pdf(file = "Trajectory0.pdf")

plot(Data1$x, Data1$y, xlab = "x, 10*km", ylab = "y, 10*km", type="l", col="white", asp=1)

grid (NULL,NULL, lty = 6, col = "gray")
lines(0, 0, type="p")
lines(Data1$x, Data1$y, type="l")

#draw.circle(0,0,2.95e+5,density=5,angle=45)

#dev.off()

#plot(Data1$x, Data1$y, xlab = "x, cm", ylab = "y, cm", type="l", asp=1, xlim=c(-5e+9, 5e+9), ylim=c(-5e+9, 5e+9))
#lines(0, 0, type="p")
#draw.circle(0,0,2.95e+9,density=5,
#            angle=45)

plot(Data1$E/Data1$E[1]-1, type="l", xlab="index", ylab=expression('E/E'[0]-1))
sd(Data1$E)
plot(Data1$E, type="l")

Data1$m <- Data1$x*Data1$vy-Data1$y*Data1$vx

plot(Data1$m/Data1$m[1]-1, type="l", xlab="index", ylab="angular moment")

#Data1["r"] <- sqrt(Data1$x^2+Data1$y^2)
#Data1["v"] <- sqrt(Data1$vx^2+Data1$vy^2)

# Обрабатываем данные для нулевого эксцентриситета

Data0N1 <- read.table ("result/result(e0,1p)N.dat",sep=",", header=T)
Data0N1["v"] <- sqrt(Data0N1$vx^2+Data0N1$vy^2)
Data0P1 <- read.table ("result/result(e0,1p)P.dat",sep=",", header=T)
Data0P1["v"] <- sqrt(Data0P1$vx^2+Data0P1$vy^2)

plot(Data0N1$v, type="l", xlab="index", ylab="v, cm/s", ylim=c(42e+8, 44e+8), col="white")
grid (NULL,NULL, lty = 6, col = "gray")
lines(Data0N1$v, type="l", col = "black")
lines(Data0P1$v, type="l", col = "red")
sd(Data0N1$v)
sd(Data0P1$v)

Data0N10 <- read.table ("result/result(e0,10p)N.dat",sep=",", header=T)
Data0N10["v"] <- sqrt(Data0N10$vx^2+Data0N10$vy^2)
Data0P10 <- read.table ("result/result(e0,10p)P.dat",sep=",", header=T)
Data0P10["v"] <- sqrt(Data0P10$vx^2+Data0P10$vy^2)

plot(Data0N10$v, type="l", xlab="index", ylab="v, cm/s", ylim=c(42e+8, 44e+8), col="white")
grid (NULL,NULL, lty = 6, col = "gray")
lines(Data0N10$v, type="l", col = "black")
lines(Data0P10$v, type="l", col = "red")
sd(Data0N10$v)
sd(Data0P10$v)

Data0N100 <- read.table ("result/result(e0,100p)N.dat",sep=",", header=T)
Data0N100["v"] <- sqrt(Data0N100$vx^2+Data0N100$vy^2)
Data0P100 <- read.table ("result/result(e0,100p)P.dat",sep=",", header=T)
Data0P100["v"] <- sqrt(Data0P100$vx^2+Data0P100$vy^2)

plot(Data0N100$v, type="l", xlab="index", ylab="v, cm/s", ylim=c(42e+8, 44e+8), col="white")
grid (NULL,NULL, lty = 6, col = "gray")
lines(Data0N100$v, type="l", col = "black")
lines(Data0P100$v, type="l", col = "red")
sd(Data0N100$v)
sd(Data0P100$v)

# Обрабатываем данные для эксцентриситета 0.1

Data01N1 <- read.table ("result/result(e01,1p)N.dat",sep=",", header=T)
Data01N1["v"] <- sqrt(Data01N1$vx^2+Data01N1$vy^2)
Data01P1 <- read.table ("result/result(e01,1p)P.dat",sep=",", header=T)
Data01P1["v"] <- sqrt(Data01P1$vx^2+Data01P1$vy^2)

plot(Data01N1$v, type="l", xlab="index", ylab="v, cm/s", ylim=c(35e+8, 50e+8), col="white")
grid (NULL,NULL, lty = 6, col = "gray")
lines(Data01N1$v, type="l", col = "black")
lines(Data01P1$v, type="l", col = "red")
sd(Data01N1$v)
sd(Data01P1$v)

Data01N10 <- read.table ("result/result(e01,10p)N.dat",sep=",", header=T)
Data01N10["v"] <- sqrt(Data01N10$vx^2+Data01N10$vy^2)
Data01P10 <- read.table ("result/result(e01,10p)P.dat",sep=",", header=T)
Data01P10["v"] <- sqrt(Data01P10$vx^2+Data01P10$vy^2)

plot(Data01N10$v, type="l", xlab="index", ylab="v, cm/s", ylim=c(35e+8, 50e+8), col="white")
grid (NULL,NULL, lty = 6, col = "gray")
lines(Data01N10$v, type="l", col = "black")
lines(Data01P10$v, type="l", col = "red")
sd(Data01N10$v)
sd(Data01P10$v)

Data01N100 <- read.table ("result/result(e01,100p)N.dat",sep=",", header=T)
Data01N100["v"] <- sqrt(Data01N100$vx^2+Data01N100$vy^2)
Data01P100 <- read.table ("result/result(e01,100p)P.dat",sep=",", header=T)
Data01P100["v"] <- sqrt(Data01P100$vx^2+Data01P100$vy^2)

plot(Data01N100$v, type="l", xlab="index", ylab="v, cm/s", ylim=c(35e+8, 50e+8), col="white")
grid (NULL,NULL, lty = 6, col = "gray")
lines(Data01N100$v, type="l", col = "black")
lines(Data01P100$v, type="l", col = "red")
sd(Data01N100$v)
sd(Data01P100$v)

# Обрабатываем данные для эксцентриситета 0.6

Data06N1 <- read.table ("result/result(e06,1p)N.dat",sep=",", header=T)
Data06N1["v"] <- sqrt(Data06N1$vx^2+Data06N1$vy^2)
Data06P1 <- read.table ("result/result(e06,1p)P.dat",sep=",", header=T)
Data06P1["v"] <- sqrt(Data06P1$vx^2+Data06P1$vy^2)

plot(Data06N1$v, type="l", xlab="index", ylab="v, cm/s", ylim=c(20e+8, 90e+8), col="white")
grid (NULL,NULL, lty = 6, col = "gray")
lines(Data06N1$v, type="l", col = "black")
lines(Data06P1$v, type="l", col = "red")
sd(Data06N1$v)
sd(Data06P1$v)

Data06N10 <- read.table ("result/result(e06,10p)N.dat",sep=",", header=T)
Data06N10["v"] <- sqrt(Data06N10$vx^2+Data06N10$vy^2)
Data06P10 <- read.table ("result/result(e06,10p)P.dat",sep=",", header=T)
Data06P10["v"] <- sqrt(Data06P10$vx^2+Data06P10$vy^2)

plot(Data06N10$v, type="l", xlab="index", ylab="v, cm/s", ylim=c(20e+8, 90e+8), col="white")
grid (NULL,NULL, lty = 6, col = "gray")
lines(Data06N10$v, type="l", col = "black")
lines(Data06P10$v, type="l", col = "red")
sd(Data06N10$v)
sd(Data06P10$v)

Data06N100 <- read.table ("result/result(e06,100p)N.dat",sep=",", header=T)
Data06N100["v"] <- sqrt(Data06N100$vx^2+Data06N100$vy^2)
Data06P100 <- read.table ("result/result(e06,100p)P.dat",sep=",", header=T)
Data06P100["v"] <- sqrt(Data06P100$vx^2+Data06P100$vy^2)

plot(Data06N100$v, type="l", xlab="index", ylab="v, cm/s", ylim=c(20e+8, 90e+8), col="white")
grid (NULL,NULL, lty = 6, col = "gray")
lines(Data06N100$v, type="l", col = "black")
lines(Data06P100$v, type="l", col = "red")
sd(Data06N100$v)
sd(Data06P100$v)

# Рассмотрим случаи с одним оборотом:

rap0N=sqrt(Data0N1$x[50000]^2+Data0N1$y[50000]^2)
rpe0N=sqrt(Data0N1$x[100000]^2+Data0N1$y[100000]^2)
exc0N=(rap0N-rpe0N)/(rap0N+rpe0N)
exc0N

rap0P=sqrt(Data0P1$x[50000]^2+Data0P1$y[50000]^2)
rpe0P=sqrt(Data0P1$x[100000]^2+Data0P1$y[100000]^2)
exc0P=(rap0P-rpe0P)/(rap0P+rpe0P)
exc0P

rap01N=sqrt(Data01N1$x[50000]^2+Data01N1$y[50000]^2)
rpe01N=sqrt(Data01N1$x[100000]^2+Data01N1$y[100000]^2)
exc01N=(rap01N-rpe01N)/(rap01N+rpe01N)
exc01N

rap01P=sqrt(Data01P1$x[50000]^2+Data01P1$y[50000]^2)
rpe01P=sqrt(Data01P1$x[100000]^2+Data01P1$y[100000]^2)
exc01P=(rap01P-rpe01P)/(rap01P+rpe01P)
exc01P

# А теперь с 10 оборотами:

Deltaexc10 <- data.frame(matrix(ncol = 6, nrow = 10))
colnames(Deltaexc10) <- c("e0N","e0P","e01N","e01P","e06N","e06P")

for (i in 1:10)
{ rap0N[i]=sqrt(Data0N10$x[10000*i-5000]^2+Data0N10$y[10000*i-5000]^2)
  rpe0N[i]=sqrt(Data0N10$x[10000*i]^2+Data0N10$y[10000*i]^2)
  exc0N[i] = (rap0N[i]-rpe0N[i])/(rap0N[i]+rpe0N[i])
  Deltaexc10$e0N[i] = exc0N[i]}
for (i in 1:10)
{ rap0P[i]=sqrt(Data0P10$x[10000*i-5000]^2+Data0P10$y[10000*i-5000]^2)
  rpe0P[i]=sqrt(Data0P10$x[10000*i]^2+Data0P10$y[10000*i]^2)
  exc0P[i] = (rap0P[i]-rpe0P[i])/(rap0P[i]+rpe0P[i])
  Deltaexc10$e0P[i] = exc0P[i]}
plot(Deltaexc10$e0P, type="l", col="red", ylab="Excentricity")
lines(Deltaexc10$e0N, type="l", col="black")

for (i in 1:10)
{ rap01N[i]=sqrt(Data01N10$x[10000*i-5000]^2+Data01N10$y[10000*i-5000]^2)
  rpe01N[i]=sqrt(Data01N10$x[10000*i]^2+Data01N10$y[10000*i]^2)
  exc01N[i] = (rap01N[i]-rpe01N[i])/(rap01N[i]+rpe01N[i])
  Deltaexc10$e01N[i] = exc01N[i]}
for (i in 1:10)
{ rap01P[i]=sqrt(Data01P10$x[10000*i-5000]^2+Data01P10$y[10000*i-5000]^2)
  rpe01P[i]=sqrt(Data01P10$x[10000*i]^2+Data01P10$y[10000*i]^2)
  exc01P[i] = (rap01P[i]-rpe01P[i])/(rap01P[i]+rpe01P[i])
  Deltaexc10$e01P[i] = exc01P[i]}
plot(Deltaexc10$e01P, type="l", col="red", ylab="Excentricity", ylim=c(0.1, 0.101))
lines(Deltaexc10$e01N, type="l", col="black")

