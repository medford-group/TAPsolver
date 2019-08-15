pracma::erf(sin(pi/(pi + 3/2)))

pracma::erf(sin(pi / 2) * 3/2)

#pracma::erf(sin(pi) * 3/2)
# half of the react or scaled by 3/2 based on the three dimensions?
tempVal = seq(.5, 1, length.out = 100)
test = pracma::erf(sin(pi * tempVal) * 3/2)
plot(test)
test[80]
pracma::erf(sin(pi/2) * pi)

gasProfile = pracma::erf(sqrt(tempVal))
plot(gasProfile)
min(gasProfile)
plot(pracma::erf((tempVal)))

tempVal = seq(0,1,length.out = 100) #/ (2 * sqrt(.01 * 1))
plot(pracma::erfc(tempVal))
gasProfile = 1-   exp(-tempVal^2 * 3 / 2  / (8) / 6)
plot(gasProfile)
min(gasProfile)


gasProfile = 1- exp(-tempVal^2 * 3 / 2  / (8) * 1.5)
plot(gasProfile)
min(gasProfile)

(pracma::erf( (3/4) ) )



1 - exp(-tempVal)
