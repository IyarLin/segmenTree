segmenTree package
================
Iyar Lin
16 July, 2020

The `segmenTree` package implements a decision tree based algorithm for
exploration of heterogeneous treatment effects in Randomized Controlled
Trial data. It supports continuous and binary response variable, and
binary treatments. It uses the rpart package as the back-end with a
predefined method on top.

You can install and load the package from github using:

``` r
pacman::p_load_gh("IyarLin/segmenTree")
```

Let’s generate a dataset with a binary response:

``` r
set.seed(1) # vary seed, n and effect_size below to get a sense of the model performance sensetivity
effect_size <- 0.25
p_x <- function(Tr, X1, X2, X3){
  lp <- 2*X1 + 0.2*X2 + as.numeric(X3)/6
  effect_size + (effect_size/0.125)*Tr*X1^3 + 
    (1 - 2*effect_size)*exp(lp)/(1+exp(lp))
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

Next, let’s import the lift method and use rpart to fit a segment tree:

``` r
lift_method <- import_lift_method()
segment_tree <- rpart(y ~ ., data = dat,
                      method = lift_method, 
                      control = rpart.control(cp = 0, minbucket = 1000),
                      x = T)
```

One way of exploring the resulting tree is by printing the rpart object:

``` r
segment_tree
```

    ## n= 10000 
    ## 
    ## node), split, n, deviance, yval
    ##       * denotes terminal node
    ## 
    ##   1) root 10000 1.000000e+16 -0.001051900  
    ##     2) X1< -0.4005301 1003 1.012054e+12 -0.150025300 *
    ##     3) X1>=-0.4005301 8997 6.552256e+15  0.013425090  
    ##       6) X1< 0.4016537 7983 4.061295e+15 -0.002719803  
    ##        12) X1< 0.2939269 6939 2.318396e+15 -0.019487100  
    ##          24) X1< -0.2554171 1467 4.631487e+12 -0.078092060 *
    ##          25) X1>=-0.2554171 5472 8.965703e+14 -0.004552848  
    ##            50) X1< 0.1890516 4430 3.851367e+14 -0.014668450  
    ##             100) X3=B,C 2992 8.013945e+13  0.008535898 *
    ##             101) X3=A 1438 4.275979e+12 -0.060894750 *
    ##            51) X1>=0.1890516 1042 1.178883e+12  0.043870860 *
    ##        13) X1>=0.2939269 1044 1.187960e+12  0.096902740 *
    ##       7) X1>=0.4016537 1014 1.057187e+12  0.152712600 *

The `segmenTree` package also contains a utility function to extract the
resulting segments:

``` r
segments <- extract_segments(segment_tree, alpha = 0.15)
print(segments)
```

    ##                            segment    n         lift   lift_lower  lift_upper
    ## 7  X1>=0.4017,X1<0.499855161411688 1014  0.152712621  0.109679301  0.19574594
    ## 6             X1>=0.2939,X1<0.4017 1044  0.096902744  0.052632159  0.14117333
    ## 5             X1>=0.1891,X1<0.2939 1042  0.043870858 -0.004383936  0.09212565
    ## 3  X1>=-0.2554,X1<0.1891, X3={B,C} 2992  0.008535898 -0.019895457  0.03696725
    ## 4    X1>=-0.2554,X1<0.1891, X3={A} 1438 -0.060894751 -0.101565692 -0.02022381
    ## 2           X1>=-0.4005,X1<-0.2554 1467 -0.078092056 -0.118100169 -0.03808394
    ## 1 X1>=-0.49972922494635,X1<-0.4005 1003 -0.150025305 -0.196506492 -0.10354412

Next, let’s compare how the segment tree predictions (black) compare
with the true treatment effect (red):

``` r
tau <- predict(segment_tree, dat)
p_treat <- p_x(rep(1, n), X1, X2, X3)
p_cont <- p_x(rep(0, n), X1, X2, X3)
cate <- p_treat - p_cont

y_lim <- c(min(tau, cate), max(tau, cate))
plot(c(min(dat$X1), max(dat$X1)), y_lim, type = "n", main = "segmenTree",
     xlab = "X1", ylab = "true (red) vs predicted (black) lift")
points(dat$X1, cate, col = "red")
points(dat$X1, tau)
```

![](README_files/figure-gfm/predict%20treatment%20effect%20and%20compare%20with%20actual%20treatment%20effect-1.png)<!-- -->

We can see that the algorithm predicts well most of the curve, but did
some over-fitting in the middle section.

The `segmenTree::prune_tree` function uses cross validation to find the
optimal `cp` hyper parameter value.

Below we can see the effective lift and error as a function of the `cp`
parameter:

``` r
optimal_cp_cv <- prune_tree(segment_tree)
```

![](README_files/figure-gfm/prune%20tree%20using%20prune_tree-1.png)<!-- -->

``` r
optimal_cp_cv <- optimal_cp_cv$optimal_cp
pruned_segment_tree <- prune(segment_tree, cp = optimal_cp_cv)
```

The optimal `cp` maximizes the ratio between effective lift and error.

Let’s check the new pruned tree predictions:

``` r
tau <- predict(pruned_segment_tree, dat)

y_lim <- c(min(tau, cate), max(tau, cate))
plot(c(min(dat$X1), max(dat$X1)), y_lim, type = "n", main = "segmenTree",
     xlab = "X1", ylab = "true (red) vs predicted (black) lift")
points(dat$X1, cate, col = "red")
points(dat$X1, tau)
```

![](README_files/figure-gfm/predict%20treatment%20effect%20pruned%20tree-1.png)<!-- -->

Looks like the tree got pruned a little too much.

One other way we can reduce over-fitting is by adjusting how the leaf
size is weighted.

``` r
lift_method <- import_lift_method(f_n = function(x) x)

weighted_segment_tree <- rpart(y ~ ., data = dat,
                      method = lift_method, 
                      control = rpart.control(cp = 0, minbucket = 1000),
                      x = T)
```

Let’s see the resulting predictions:

``` r
tau <- predict(weighted_segment_tree, dat)

y_lim <- c(min(tau, cate), max(tau, cate))
plot(c(min(dat$X1), max(dat$X1)), y_lim, type = "n", main = "segmenTree",
     xlab = "X1", ylab = "true (red) vs predicted (black) lift")
points(dat$X1, cate, col = "red")
points(dat$X1, tau)
```

![](README_files/figure-gfm/predict%20treatment%20effect%20weighted%20tree-1.png)<!-- -->

To demonstrate why we even need a specialized algorithm we’ll compare
the segmenTree model with 2 other approaches that use plain ML models:

1.  Model jointly the treatment and covariates. Predict:
    tau(x)=f(x,T=1)-f(x,T=0)  
2.  Train a model on the treatment units f\_{T=1} and a separate model
    f\_{T=0} on the control units and predict: tau(x) = f\_{T=1}(x) -
    f\_{T=0}(x). This is also called “Two model approach”.

We’ll use a decision tree in both of the above.

``` r
par(mfrow = c(1, 3))
# segmenTree pruned model
tau <- predict(pruned_segment_tree, dat)
p_treat <- p_x(rep(1, n), X1, X2, X3)
p_cont <- p_x(rep(0, n), X1, X2, X3)
cate <- p_treat - p_cont

y_lim <- c(min(tau, cate), max(tau, cate))
plot(c(min(dat$X1), max(dat$X1)), y_lim, type = "n", main = "segmenTree",
     xlab = "X1", ylab = "true (red) vs predicted (black) lift")
points(dat$X1, cate, col = "red")
points(dat$X1, tau)

# Model jointly the treatment and covariates
dat2 <- data.frame(y, X1, X2, X3, Tr)
fit2 <- rpart(y ~ ., data = dat2)

dat2_treat <- dat2; dat2_cont <- dat2
dat2_treat$Tr <- 1L; dat2_cont$Tr <- 0L
tau2 <- predict(fit2, dat2_treat) - predict(fit2, dat2_cont)

y_lim <- c(min(tau2, cate), max(tau2, cate))
plot(c(min(dat$X1), max(dat$X1)), y_lim, type = "n", main = "regular model",
     xlab = "X1", ylab = "")
points(dat$X1, cate, col = "red")
points(dat$X1, tau2)

# Train a model on the treatment units and a seperate model on the control units

dat3_treat <- dat2[dat2$Tr == 1, -5]
dat3_cont <- dat2[dat2$Tr == 0, -5]

fit3_treat <- rpart(y ~ ., data = dat3_treat)
fit3_cont <- rpart(y ~ ., data = dat3_cont)
dat2_treat <- dat2; dat2_cont <- dat2
tau3 <- predict(fit3_treat, dat2) - predict(fit3_cont, dat2)

y_lim <- c(min(tau3, cate), max(tau3, cate))
plot(c(min(dat$X1), max(dat$X1)), y_lim, type = "n", main = "Two models",
     xlab = "X1", ylab = "")
points(dat$X1, cate, col = "red")
points(dat$X1, tau3)
```

![](README_files/figure-gfm/compare%20segmenTree%20with%202%20other%20approches-1.png)<!-- -->
