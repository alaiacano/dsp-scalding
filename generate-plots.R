setwd("~/Dropbox/ml-presentation/")
library(ggplot2)
library(reshape)
rm(list=ls())

# weekly trend data generation

T = 30              # number of periods to use
Ts = 1              # sampling period
f <- 1 / (Ts * 7)   # signal frequency (weekly)

d <- data.frame(t = seq(0, 5*T, by=Ts))

d$y <- sin( 2 * pi * f * d$t) + # variation
       4 * sin(d$t / 15)        # trend

# take 7-day moving average
d$trend7 <- filter(d$y, rep(1/7, 7))
d$trend28 <- filter(d$y, rep(1/28, 28))

p <- ggplot(d, aes(x=t)) +
  geom_line(aes(y=y)) +
  scale_x_continuous("Time") +
  scale_y_continuous("Signal") +
  ggtitle("Daily Users")
ggsave("figures/01-users-weekly.png", plot=p, height=7, width=12, dpi=72)

p2 <- p + geom_line(aes(y=trend7), color='red')
ggsave("figures/02-users-weekly-trend.png", plot=p2, height=7, width=12, dpi=72)

p3 <- p2 + geom_line(aes(y=trend28), color='blue')
ggsave("figures/02a-users-weekly-monthly-trend.png", plot=p3, height=7, width=12, dpi=72)
## picture of ideal low-pass filter
rm(list=ls())
d <- data.frame(
  x = seq(0,1,length.out=10000),
  y = c(rep(1, 2000), rep(0, 8000)),
  gap = c(rep(1,1500), seq(1,0, length.out=1000), rep(0, 7500))
)
p <- ggplot(d, aes(x=x)) + 
  geom_line(aes(y=y)) +
  scale_x_continuous("Frequency") +
  scale_y_continuous("Magnitude")
ggsave("figures/03-lowpass-diagram.png", plot=p, height=7, width=12, dpi=72)

p2 <- p + geom_line(aes(y=gap), size=1.5)
ggsave("figures/04-lowpass-diagram-2.png", plot=p2, height=7, width=12, dpi=72)

###################################### fourier transform demo
rm(list=ls())
Ts = .001                # sampling period
fn = 1 / (2 * Ts)        # nyquist frequency
t <- seq(0,5,by=Ts)      # generate a time vector

filt <- c(-0.020104118828857324, -0.05842798004352509, -0.06117840364782199, -0.01093939338533898, 0.0512509644353497, 0.033220867678947864, -0.056552769718339224, -0.08565500737264507, 0.06337959966054496, 0.31085440365663597, 0.4344309124179416, 0.31085440365663597, 0.06337959966054496, -0.08565500737264507, -0.056552769718339224, 0.033220867678947864, 0.0512509644353497, -0.01093939338533898, -0.06117840364782199, -0.05842798004352509, -0.020104118828857324)

d <- data.frame(
  t = t,
  x = sin(2*pi*.5*t) + .5 * sin(2*pi*250*t) + .1 * sin(2*pi*400*t),
  x2 = 10 * sin(2*pi*.5*t) + .5 * sin(2*pi*3*t)# + .1 * sin(2*pi*10*t)
)

d$filt <- filter(d$x, filt)
d <- subset(d, !is.na(d$filt))

p <- qplot(t, x, data=subset(d, t < 2), geom='line') +
  scale_x_continuous("Time (seconds)", limits=c(0, 2)) +
  scale_y_continuous("Amplitude")
ggsave("figures/08-fft-signal-time-domain.png", plot=p, height=3.5, width=6, dpi=72)

X <- abs(fft(d$x))
Xf <- abs(fft(d$filt))
D <- data.frame(
  f = seq(0, fn, length.out=floor(length(X)/2)),
  X = X[1:floor(length(X)/2)],
  Xf = Xf[1:floor(length(Xf)/2)]
)

p2 <- qplot(f, X/max(X), data=D, geom='line') +
  scale_x_continuous("Frequency (Hz)") +
  scale_y_continuous("Magnitude")
ggsave("figures/09-fft-signal-frequency-domain.png", plot=p2, height=3.5, width=6, dpi=72)

p3 <- qplot(f, Xf/max(Xf), data=D, geom='line') +
  scale_x_continuous("Frequency (Hz)") +
  scale_y_continuous("Magnitude")
ggsave("figures/10-fft-signal-frequency-domain-filtered.png", plot=p3, height=3.5, width=6, dpi=72)

p4 <- qplot(t, filt, data=subset(d, t < 2), geom='line') +
  scale_x_continuous("Time (seconds)") +
  scale_y_continuous("Amplitude")
ggsave("figures/11-fft-signal-time-domain-filtered.png", plot=p4, height=3.5, width=6, dpi=72)

########### white noise example.
rm(list=ls())
# filt <- c(-0.020104118828857324, -0.05842798004352509, -0.06117840364782199, -0.01093939338533898, 0.0512509644353497, 0.033220867678947864, -0.056552769718339224, -0.08565500737264507, 0.06337959966054496, 0.31085440365663597, 0.4344309124179416, 0.31085440365663597, 0.06337959966054496, -0.08565500737264507, -0.056552769718339224, 0.033220867678947864, 0.0512509644353497, -0.01093939338533898, -0.06117840364782199, -0.05842798004352509, -0.020104118828857324)
Ts = .001                # sampling period
fn = 1 / (2 * Ts)        # nyquist frequency
t <- seq(0,5,by=Ts)      # generate a time vector

# make a data.frame with all of our signals
# the runif function will generate an approximately normal distribution of *frequencies*.
# Just take a histogram of the fft to see what I mean.
d <- data.frame(
  t  = t,
  v  = runif(length(t))
)
d$filtered <- filter(d$v, rep(1/7, 7))    # 7-point average
d$filtered2 <- filter(d$v, filt) # 28-point average

# There are some NA values left over from the filter.
# It's OK to just drop them before going into the fft
d <- subset(d, !is.na(filtered) & !is.na(filtered2))

# take the FFT of the four signals
V = abs(fft(d$v))
Vf = abs(fft(d$filtered))
Vf2 = abs(fft(d$filtered2))

# put all of them into a data.frame for plotting
D <- data.frame(
  f  = seq(0, fn, length.out=floor(length(V)/2)),
  V  = V[1:floor(length(V)/2)],
  Vf = Vf[1:floor(length(Vf)/2)],
  Vf2 = Vf2[1:floor(length(Vf2)/2)]
)

# one of these has a really huge spike at DC (zero Hz). Drop it for plotting.
D <- D[2:nrow(D),]
p <- ggplot(D, aes(f,V)) + 
  geom_line() +
  geom_smooth() +
  scale_x_continuous("Frequency (Hz)") +
  scale_y_continuous("Magnitude")
ggsave("figures/05-white-noise-fft.png", plot=p, height=7, width=12, dpi=72)

p2 <- ggplot(D, aes(x=f)) + 
  geom_line(aes(y=V), alpha=.2) +
  geom_line(aes(y=Vf)) +
  geom_smooth(aes(y=Vf)) +
  scale_x_continuous("Frequency (Hz)") +
  scale_y_continuous("Magnitude")
ggsave("figures/06-white-noise-fft-filtered7.png", plot=p2, height=7, width=12, dpi=72)

p3 <- ggplot(D, aes(x=f)) + 
  geom_line(aes(y=V), alpha=.2) +
  geom_line(aes(y=Vf), alpha=.3) +
  geom_line(aes(y=Vf2)) +
  geom_smooth(aes(y=Vf2)) +
  scale_x_continuous("Frequency (Hz)") +
  scale_y_continuous("Magnitude")
ggsave("figures/07-white-noise-fft-filtered7,30.png", plot=p3, height=7, width=12, dpi=72)


############ Impulse response (hard-coded, sorry!)
rm(list=ls())
t <- -1:20
fir <- c(0, rep(1/7,7), rep(0, 14))

stemm <- function(x, y, impulse=FALSE) {
  d <- data.frame(x=x, y=y)
  p <- ggplot(d, aes(x=x))
  # add the impulse line before the other data.
  if (impulse == TRUE) {
    p <- p + geom_linerange(aes(x=0, ymin=0, ymax=1), color='red', alpha=.5)
  }
  p <- p +
    geom_point(aes(y=y)) +
    geom_linerange(aes(ymin=sapply(y, min, 0), ymax=sapply(y, max, 0)))

  p
}

p <- stemm(t, fir, impulse=TRUE) +
  scale_x_continuous("Time") +
  scale_y_continuous("Impulse Response") +
  ggtitle("FIR filter - 7 day MA")
ggsave("figures/11-fir-impulse-response.png", plot=p, height=3.5, width=6, dpi=72)

iir <- c(0.0, 0.5, 0.25, 0.125, 0.0625, 0.03125, 0.015625, 0.0078125, 0.00390625, 0.001953125, 9.765625E-4, 4.8828125E-4, 2.44140625E-4, 1.220703125E-4, 6.103515625E-5, 3.0517578125E-5, 1.52587890625E-5, 7.62939453125E-6, 3.814697265625E-6, 1.9073486328125E-6, 9.5367431640625E-7, 4.76837158203125E-7)
p2 <- stemm(t, iir, impulse=TRUE) +
  scale_x_continuous("Time") +
  scale_y_continuous("Impulse Response") +
  ggtitle("IIR filter - Exponential MA (.5, .5)")
ggsave("figures/12-iir-impulse-response.png", plot=p2, height=3.5, width=6, dpi=72)

iir.unstable <- c(0.0, 0.5, 0.55, 0.6050000000000001, 0.6655000000000002, 0.7320500000000003, 0.8052550000000004, 0.8857805000000005, 0.9743585500000007, 1.0717944050000008, 1.178973845500001, 1.2968712300500012, 1.4265583530550014, 1.5692141883605017, 1.7261356071965521, 1.8987491679162074, 2.088624084707828, 2.297486493178611, 2.527235142496472, 2.7799586567461194, 3.0579545224207316, 3.363749974662805)
p3 <- stemm(t, iir.unstable, impulse=TRUE) +
  scale_x_continuous("Time") +
  scale_y_continuous("Impulse Response") +
  ggtitle("IIR filter - Exponential MA (1.1, .5)")
ggsave("figures/13-iir-impulse-response-unstable.png", plot=p3, height=3.5, width=6, dpi=72)



############### sampling
rm(list=ls())
generate.samples <- function(Tsamp, Npoints, offset=0) {
  ((0:(Npoints-1)) %% round(1/Tsamp) == offset)
}
Ts = .001                # sampling period
fn = 1 / (2 * Ts)        # nyquist frequency
t <- seq(0,5,by=Ts)      # generate a time vector

d <- data.frame(
  t = t,
  x1 = sin(2*pi*t)
)

# original graph (no sampling)
p <- qplot(t, x1, data=d, geom='line') +
  scale_x_continuous("Periods") +
  scale_y_continuous("")
ggsave("figures/14-sampling-raw.png", plot=p, height=300/72, width=500/72, dpi=72)

make.graph <- function(d, div, offset=0) {
  d.samp <- subset(d, generate.samples(Ts/div, nrow(d), offset))
  p <- ggplot() +
  geom_line(data=d, aes(x=t, y=x1),alpha=.5) +
  geom_linerange(data=d.samp, aes(x=t, ymin=0, ymax=x1), alpha=.3) +
  geom_point(data=d.samp, aes(x=t, y=x1), color='red') +
  geom_line(data=d.samp, aes(x=t, y=x1), color='red') +
  scale_x_continuous("Periods") +
  scale_y_continuous("")
  p
}

p <- make.graph(d, .1)
ggsave("figures/15-sampling-.1.png", plot=p, height=300/72, width=500/72, dpi=72)

p <- make.graph(d, .5, 100)
ggsave("figures/16-sampling-.5.png", plot=p, height=300/72, width=500/72, dpi=72)

p <- make.graph(d, .88)
ggsave("figures/17-sampling-.8.png", plot=p, height=300/72, width=500/72, dpi=72)

p <- make.graph(d, .47)
ggsave("figures/18-sampling-..47.png", plot=p, height=300/72, width=500/72, dpi=72)

p <- make.graph(d, 1, 100)
ggsave("figures/19-sampling-1.png", plot=p, height=300/72, width=500/72, dpi=72)



################### Toeplitz style
rm(list=ls())
t <- seq(0, 6*pi, length.out=1000)

x <- data.frame(
  x1 = sin(2*pi*t),
  x2 = 1/3 * sin(3 * 2*pi*t),
  x3 = 1/5 * sin(5 * 2*pi*t)
)
x$x4 <- with(x, x1+x2+x3)  # square wave!

X <- as.matrix(x)

# WARNING! This apparently doesn't build the right type of matrix for what we
# want. See toeplitz_filter.jl for a better implementation.
H <- toeplitz(c(rep(1/7, 7), rep(0, nrow(x)-7)))

Xf <- t(t(X) %*% H)

names(Xf) <- paste(names(x), "f", sep="")
Xf$t <- t
x <- cbind(x, Xf)

x.m <- melt.data.frame(x, id.vars="t")
x.m$filtered <- factor(ifelse(!grepl("f", x.m$variable), "Original", "Filtered"), levels=c("Original", "Filtered"))
x.m$p <- substr(x.m$variable, 2, 2)
x.m <- subset(x.m, t <= pi)
p <- qplot(t, value, data=x.m, color=p, geom='line') +
  facet_grid(p ~ filtered) +
  theme(legend.position="none")
ggsave("figures/20-matrix-filter.png", plot=p, height=7, width=12, dpi=72)