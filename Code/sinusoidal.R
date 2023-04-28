# Trying a new idea

x = runif(100)
y = sin(x * 4 * pi) + rnorm(100, sd = 0.2)
plot(x, y)
title("Noisy sin wave")
