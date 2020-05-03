library(segmenTree)
# In the script below you can vary the seed, the effect_size and sample size 
# to get a sense of model performance sensitivity
set.seed(1) 


# Generate a dataset
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

# Fit a segment tree
lift_method <- import_lift_method()
segment_tree <- rpart(y ~ ., data = dat,
                      method = lift_method, 
                      control = rpart.control(cp = 0, minbucket = 1000),
                      x = T)

# explore the resulting segments
segments <- extract_segments(segment_tree, alpha = 0.15)
print(segments)

# Predict treatment effect and compare with actual treatment effect
tau <- predict(segment_tree, dat)
p_treat <- p_x(rep(1, n), X1, X2, X3)
p_cont <- p_x(rep(0, n), X1, X2, X3)
cate <- p_treat - p_cont

y_lim <- c(min(tau, cate), max(tau, cate))
plot(c(min(dat$X1), max(dat$X1)), y_lim, type = "n", main = "segmenTree",
     xlab = "X1", ylab = "true (red) vs predicted (black) lift")
points(dat$X1, cate, col = "red")
points(dat$X1, tau)

# find optimal cp using cross validation
optimal_cp_cv <- tune_cp(segment_tree, cp_num = 6, train_frac = 0.5, M = 100)
optimal_cp_cv <- optimal_cp_cv$optimal_cp
pruned_segment_tree <- prune(segment_tree, cp = optimal_cp_cv)

tau <- predict(pruned_segment_tree, dat)
y_lim <- c(min(tau, cate), max(tau, cate))

par(mfrow = c(1, 3))
plot(c(min(dat$X1), max(dat$X1)), y_lim, type = "n", main = "segmenTree - pruned",
     xlab = "X1", ylab = "true (red) vs predicted (black) lift")
points(dat$X1, cate, col = "red")
points(dat$X1, tau)

# Compare to a regular classfication model
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

# Compare to training seperately on control and treatment
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

# Preferred lift
p_x <- function(Tr, X1, X2, X3, X4){
  lp <- X1 - 0.2*X2 + as.numeric(X3)/6 - X4/3
  0.25 + 
    Tr*0.25*(1-X2)*X1^2 -
    Tr*0.25*(1-X1)*X2^2 + 
    0.5*exp(lp)/(1+exp(lp))
}

n <- 10000
Tr <- rbinom(n, 1, 0.3)
X1 <- runif(n, 0, 1)
X2 <- runif(n, 0, 1)
X3 <- factor(sample(LETTERS[1:3], size = n, replace = T))
X4 <- rnorm(n)
p <- p_x(Tr, X1, X2, X3, X4)
y <- sapply(p, function(x) rbinom(1, 1, x))
y_mat <- cbind(y, Tr)
dat <- data.frame(y = I(y_mat), X1, X2, X3)

lift_method <- import_lift_method()
segment_tree <- rpart(y ~ ., data = dat,
                      method = lift_method, 
                      parms = list(preferred_lift = "lower"), 
                      control = rpart.control(cp = 0, minbucket = 1000),
                      x = T)

# explore the resulting segments
segments <- extract_segments(segment_tree, alpha = 0.15)
print(segments)
