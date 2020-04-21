segmenTree package
================
Iyar Lin
21 April, 2020

The segmenTree package implements a decision tree based algorithm for
exploration of heterogeneous treatment effects in Randomised Controlled
Trial data. It uses the rpart package as the backend with a predifined
method on top.

You can install and load the package from github using:

``` r
pacman::p_load_gh("IyarLin/segmenTree")
```

Let’s generate a dataset:

``` r
set.seed(1) # vary this and the sample size to get a sense of the model performance
p_x <- function(Tr, X1, X2, X3){
  lp <- -X1 + 0.2*X2 + as.numeric(X3)/6
  0.25 + 2*Tr*X1^3 + 0.5*exp(lp)/(1+exp(lp))
}

n <- 10000
Tr <- rbinom(n, 1, 0.3)
X1 <- runif(n, -0.5, 0.5)
X2 <- rnorm(n)
X3 <- factor(sample(LETTERS[1:3], size = n, replace = T))
p <- p_x(Tr, X1, X2, X3)
y <- sapply(p, function(x) rbinom(1, 1, x))
y_mat <- cbind(y, Tr)
dat <- data.frame(y = I(y_mat), X1, X2, X3)
```

Next, let’s import the lift method and use rpart to fit a causal tree:

``` r
lift_method <- import_lift_method()
causal_tree <- rpart(y ~ ., data = dat,
              method = lift_method,
              control = rpart.control(cp = 0.001),
              parms = list(), y = T)
```

We used `cp = 0.001` above based on some experimentation. A function for
tuning the `cp` paramter using cross validation is in the works.

Let’s explore the resulting segments:

``` r
segments <- extract_segments(causal_tree, alpha = 0.15)
pander::pandoc.table(segments)
```

|       |         segment         |  n   |   lift    | lift\_lower | lift\_upper |
| :---: | :---------------------: | :--: | :-------: | :---------: | :---------: |
| **3** |   X1\>=0.3228,X1\<Inf   | 1781 |  0.1152   |   0.07867   |   0.1518    |
| **2** | X1\>=-0.4742,X1\<0.3228 | 7969 | \-0.02376 |  \-0.04125  | \-0.006277  |
| **1** |  X1\>=-Inf,X1\<-0.4742  | 250  | \-0.3055  |  \-0.3985   |  \-0.2126   |

Let’s compare the estimated lift (in black) vs the true lift (in red):

``` r
tau <- predict(causal_tree, dat)
p_treat <- p_x(rep(1, n), X1, X2, X3)
p_cont <- p_x(rep(0, n), X1, X2, X3)
cate <- p_treat - p_cont

y_lim <- c(min(tau, cate), max(tau, cate))
plot(c(min(dat$X1), max(dat$X1)), y_lim, type = "n", main = "segmenTree",
     xlab = "X1", ylab = "true (red) vs predicted (black) lift")
points(dat$X1, cate, col = "red")
points(dat$X1, tau)
```

![](README_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

Just to get a sense of why we need this kind of model, let’s try to
estimate lift by modelling both the treatment and covariates:

``` r
dat2 <- data.frame(y, X1, X2, X3, Tr)
fit2 <- rpart(y ~ ., data = dat2, control = rpart.control(cp = 0))
fit2 <- rpart(y ~ ., data = dat,
                     method = lift_method,
                     control = rpart.control(cp = fit2$cptable[max(which(fit2$cptable[, 2] < 3)), 1]),
                     parms = list())

dat2_treat <- dat2; dat2_cont <- dat2
dat2_treat$Tr <- 1L; dat2_cont$Tr <- 0L
tau2 <- predict(fit2, dat2_treat) - predict(fit2, dat2_cont)

y_lim <- c(min(tau2, cate), max(tau2, cate))
plot(c(min(dat$X1), max(dat$X1)), y_lim, type = "n", main = "regular model",
     xlab = "X1", ylab = "")
points(dat$X1, cate, col = "red")
points(dat$X1, tau2)
```

![](README_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

Let’s also try fitting a seperate model for treatment and control:

``` r
dat3_treat <- dat2[dat2$Tr == 1, -5]
dat3_cont <- dat2[dat2$Tr == 0, -5]

fit3_treat <- rpart(y ~ ., data = dat3_treat)
fit3_cont <- rpart(y ~ ., data = dat3_cont)
dat2_treat <- dat2; dat2_cont <- dat2
tau3 <- predict(fit3_treat, dat2) - predict(fit3_cont, dat2)

y_lim <- c(min(tau3, cate), max(tau3, cate))
plot(c(min(dat$X1), max(dat$X1)), y_lim, type = "n", main = "seperate models",
     xlab = "X1", ylab = "")
points(dat$X1, cate, col = "red")
points(dat$X1, tau3)
```

![](README_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->
