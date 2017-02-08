#install.packages('DiceKriging')
#install.packages('DiceOptim')
library(DiceKriging)
library(DiceOptim)

# 1D Example with known parameters -> simple kriging
inputs = c(-1,-0.5,0,0.5,1)
output = c(-9,-5,-1,9,11)
theta = 0.4
sigma = 5
trend = c(0,11,2)
model = km(formula = ~x+I(x^2), design = data.frame(x=inputs), 
           response = output, covtype = "matern5_2",
           coef.trend = trend, coef.cov = theta, coef.var = sigma^2)
model

t = seq(from = -2, to = 2, length = 200)
p = predict(model, newdata = data.frame(x=t), type = "SK")

plot(t, p$mean, type = "l", xlim = c(-2, 2), ylim = c(-30, 30),
     xlab = "x", ylab = "y")
lines(t, p$lower95, col = "black", lty = 2)
lines(t, p$upper95, col = "black", lty = 2)
points(inputs, output, col = "red", pch = 19)
abline(h = 0)

#Simulations of Gaussian processes underlying Kringing
 #Unconditional
y = simulate(model, nsim = 5, newdata = data.frame(x=t))

ytrend = trend[1]+trend[2]*t+trend[3]*t^2
par(mfrow = c(1, 1))
plot(t, ytrend, type = "l", col = "black", ylab = "y", lty = "dashed",
     ylim = c(min(ytrend) - 2 * sigma, max(ytrend) + 2 * sigma))
for(i in 1:5) lines(t, y[i, ], col = i)

 #Conditional (known points)
y = simulate(model, nsim = 5, newdata = data.frame(x=t), cond = TRUE)

plot(t, ytrend, type = "l", col = "black", ylab = "y", lty = "dashed",
     ylim = c(min(ytrend) - 2 * sigma, max(ytrend) + 2 * sigma))
for(i in 1:5) lines(t, y[i, ], col = i)
points(inputs, output, col = "red", pch = 19)


## 2D data
X <- expand.grid(x1 = seq(0, 1, length = 4), x2 = seq(0, 1, length = 4))
y <- apply(X, 1, branin)
m = km(~., design = X, response = y, covtype = "gauss")

 #concentrates likelihood
n.grid = 30
x.grid = seq(0.01,2,length = n.grid)
X.grid = expand.grid(x.grid,x.grid)
logLik.grid = apply(X.grid,1,logLikFun,m)

contour(x.grid, x.grid, matrix(logLik.grid, n.grid, n.grid), 40,
            xlab = expression(theta[1]), ylab = expression(theta[2]))
opt <- m@covariance@range.val
points(opt[1], opt[2], pch = 19, col = "red")

 #draw kriging mean and visualize pediction accuracy
n.grid = 50
x.grid = seq(0,1,length = n.grid)
X.grid = expand.grid(x1=x.grid, x2 = x.grid)
y.grid = apply(X.grid,1,branin)
pred.m = predict(m,X.grid,"UK")
par(mfrow = c(1,3))

contour(x.grid, x.grid, matrix(y.grid, n.grid, n.grid), 50,
            main = "Branin")
points(X[ , 1], X[ , 2], pch = 19, cex = 1.5, col = "red")
contour(x.grid, x.grid, matrix(pred.m$mean, n.grid, n.grid), 50,
               main = "Kriging mean")
points(X[ , 1], X[ , 2], pch = 19, cex = 1.5, col = "red")
contour(x.grid, x.grid, matrix(pred.m$sd^2, n.grid, n.grid), 15,
               main = "Kriging variance")
points(X[ , 1], X[ , 2], pch = 19, cex = 1.5, col = "red")

plot(m)



## My data
library(R.matlab)

residuals = readMat('/Users/Alexandra/Documents/SFU/Data/SnowDepth/Kriging/residuals.mat')
res = residuals$res
utm = residuals$utm

X <- data.frame(utm)
y <- res
m = km(~.,design = X, response = y, covtype = "matern5_2")
m = km(~., design = X, response = y, covtype = "gauss", nugget = 1e-8 * var(y))
plot(m)
m

 #concentrates likelihood
n.grid = 30+1
x.grid = seq(0,30,length = n.grid)
X.grid = expand.grid(x.grid,x.grid)
logLik.grid = apply(X.grid,1,logLikFun,m)

contour(x.grid, x.grid, matrix(logLik.grid, n.grid, n.grid), 40,
        xlab = expression(theta[1]), ylab = expression(theta[2]))
opt <- m@covariance@range.val
points(opt[1], opt[2], pch = 19, col = "red")

#draw kriging mean and visualize pediction accuracy
n.grid = 51
x.grid = seq(0,3000,length = n.grid)
X.grid = expand.grid(X1=x.grid, X2 = x.grid)
pred.m = predict(m,X.grid,"SK")
par(mfrow = c(1,3))


#contour(x.grid, x.grid, matrix(y.grid, n.grid, n.grid), 50,
 #       main = "Branin")
#points(X[ , 1], X[ , 2], pch = 19, cex = 1.5, col = "red")
contour(x.grid, x.grid, matrix(pred.m$mean, n.grid, n.grid), 20,
        main = "Kriging mean")
points(X[ , 1], X[ , 2], pch = 19, cex = 1.5, col = "red")
#contour(x.grid, x.grid, matrix(pred.m$sd^2, n.grid, n.grid), 15,
  #      main = "Kriging variance")
#points(X[ , 1], X[ , 2], pch = 19, cex = 1.5, col = "red")


### 1D attempt
# 1D Example with known parameters -> simple kriging
inputs = data.frame(utm[1:5 ,1])
colnames(inputs) = c("V1")
output = res[1:5]
model = km(formula = ~., design = inputs, 
           response = output, covtype = "gauss", nugget = 1e-8 * var(output))
model

t = seq(from = 0, to = max(inputs), length = 5)
new.data = data.frame(x=utm[5:9 ,1])
colnames(new.data) <- c("V1")
p = predict(model, newdata = new.data, type = "SK")

plot(t, p$mean, type = "l", xlab = "x", ylab = "y", ylim = c(0, 0.3))
  lines(t, p$lower95, col = "black", lty = 2)
  lines(t, p$upper95, col = "black", lty = 2)
  points(utm[1:5 ,1], output, col = "red", pch = 19)

#Simulations of Gaussian processes underlying Kringing
#Unconditional
y = simulate(model, nsim = 5, newdata = data.frame(x=t))

ytrend = trend[1]+trend[2]*t+trend[3]*t^2
par(mfrow = c(1, 1))
plot(t, ytrend, type = "l", col = "black", ylab = "y", lty = "dashed",
     ylim = c(min(ytrend) - 2 * sigma, max(ytrend) + 2 * sigma))
for(i in 1:5) lines(t, y[i, ], col = i)

#Conditional (known points)
y = simulate(model, nsim = 5, newdata = data.frame(x=t), cond = TRUE)

plot(t, ytrend, type = "l", col = "black", ylab = "y", lty = "dashed",
     ylim = c(min(ytrend) - 2 * sigma, max(ytrend) + 2 * sigma))
for(i in 1:5) lines(t, y[i, ], col = i)
points(inputs, output, col = "red", pch = 19)
