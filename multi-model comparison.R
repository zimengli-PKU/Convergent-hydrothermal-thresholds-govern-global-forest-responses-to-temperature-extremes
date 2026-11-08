data<-read.csv('response.csv')
data$order<-1:dim(data)[1]
load ('TP60.RData')
data$MAT<-T60$MAT
data$MAP<-P60$MAP

#install.packages("modelr",repos="http://mirrors.tuna.tsinghua.edu.cn/CRAN/")

library(nlme)  #nonlinear regression
library(ggplot2)
library(AICcmodavg)
library(modelr)
dataset<-na.omit(data[,c('MAP','MAT','hotr')])
colnames(dataset)[3]<-'response'

zerovalue<-data.frame(x=NA,y=NA,p=NA,p_inc=NA,k=NA,std=NA,num=NA,uniq=NA)
for (i in 1:50){
  findzero<-which(dataset$MAP>100*(i-1)&dataset$MAP<=100*i)
  datazero<-dataset[findzero,c('MAT','response')]
  datazero<-na.omit(datazero)
  zerovalue[i,8]<-length(unique(datazero$MAT))
  if (length(findzero)>20)
  {
    lmzero<-lm(response~MAT,data=datazero)
    pzero<-summary(lmzero)$coefficients[,4]
    zerovalue[i,3:6]<-c(pzero[2],pzero[1],summary(lmzero)$coefficients[2,1:2])
    if (pzero[2]<0.05)
    {
      zerovalue[i,1]<-i
      zerovalue[i,2]<-lmzero$coefficients[1]/lmzero$coefficients[2]*(-1)
      
    }
  }
}

datah <- data.frame(x = zerovalue$x, y = zerovalue$y)
datah<-na.omit(datah)
# 1. linear
linear_model <- lm(y ~ x, data = datah)
summary(linear_model)
pred <- predict(linear_model)
rmse <- sqrt(mean((datah$y - pred)^2))
mae  <- mean(abs(datah$y - pred))
# 2. quadratic
quadratic_model <- lm(y ~ x + I(x^2), data = datah)
summary(quadratic_model)
pred <- predict(quadratic_model)
rmse <- sqrt(mean((datah$y - pred)^2))
mae  <- mean(abs(datah$y - pred))
# 3. linear-plateau
linear_plateau <- function(x, a, b, c) {
  pmin(a + b * x, a + b * c)  # 当 x >= c 时，y = a + b * c（高原）
}

linear_plateau_model <- nls(y ~ linear_plateau(x, a, b, c),
                            data = datah,
                            start = list(a = -20, b = 4,
                                         c = 8))#-5,1,10;-10,0.01,10
library(rcompanion)

plotPredy(data  = datah,
          x     = x,
          y     = y,
          model = linear_plateau_model)

a<-zerovalue$y
a<-a[!is.na(a)]
cor(fitted(linear_plateau_model),a)^2
coef(linear_plateau_model)["a"]+coef(linear_plateau_model)["b"]*coef(linear_plateau_model)["c"]

pred <- predict(linear_plateau_model)
rmse <- sqrt(mean((datah$y - pred)^2))
mae  <- mean(abs(datah$y - pred))

SSE <- sum((datah$y - pred)^2)
SST <- sum((datah$y - mean(datah$y))^2)

R2_linear_plateau <- 1 - SSE / SST
adj_R2_linear_plateau <- 1 -
  (SSE / (nrow(datah) - 3)) / (SST / (nrow(datah) - 1))
# 4. qua-plateau

fit.lm<-lm(y ~ poly(x,2, raw=TRUE),data=datah)
a_parameter<-fit.lm$coefficients[1]
b_parameter<-fit.lm$coefficients[2]
c_parameter<-fit.lm$coefficients[3]
x_mean<-mean(datah$x)


qp <- function(x, a, b, jp) {
  c <- -0.5 * b / jp
  if_else(condition = x < jp,
          true  = a + (b * x) + (c * x * x),
          false = a + (b * jp) + (c * jp * jp))
}

library(dplyr)
quadratic_plateau_model <- nls(y ~ qp(x, a, b, jp),
                               data = datah,
                               start = list(a = a_parameter, b = b_parameter, jp = x_mean))#,nls.control(maxiter=1000

summary(quadratic_plateau_model)
coef(quadratic_plateau_model)["b"]/2/coef(quadratic_plateau_model)["jp"]*(-1) #c value
coef(quadratic_plateau_model)["a"]+coef(quadratic_plateau_model)["b"]*coef(quadratic_plateau_model)["jp"]/2


nullfunct = function(x, m){m}
m.ini= mean(datah$y)
null = nls(y ~ nullfunct(x, m),
           data = datah,
           start = list(m = m.ini),
           trace = FALSE)

a<-datah$y
a<-a[!is.na(a)]
cor(fitted(quadratic_plateau_model),a)^2  # 对于 nls 模型，近似计算

pred <- predict(quadratic_plateau_model)
rmse <- sqrt(mean((datah$y - pred)^2))
mae  <- mean(abs(datah$y - pred))

datah$pred_linear <- predict(linear_model)
datah$pred_quadratic <- predict(quadratic_model)
datah$pred_linear_plateau <- predict(linear_plateau_model)
datah$pred_quadratic_plateau <- predict(quadratic_plateau_model)

p_models<-ggplot(datah, aes(x = x, y = y)) +
  geom_point(color = "black", alpha = 0.5) +
  geom_line(aes(y = pred_linear, color = "Linear"), size = 1) +
  geom_line(aes(y = pred_quadratic, color = "Quadratic"), size = 1) +
  geom_line(aes(y = pred_linear_plateau, color = "Linear-Plateau"), size = 1) +
  geom_line(aes(y = pred_quadratic_plateau, color = "Quadratic-Plateau"), size = 1) +
  labs(title = "",
       x = " Binned MAP(mm)", y = expression('T'[threshold])) +
  scale_color_manual(values = c("Linear" = "blue", "Quadratic" = "red",
                                "Linear-Plateau" = "green", "Quadratic-Plateau" = "purple")) +
  theme_minimal()


model_list <- list(
  linear = linear_model,
  quadratic = quadratic_model,
  linear_plateau = linear_plateau_model,
  quadratic_plateau = quadratic_plateau_model
)
print(sapply(model_list, nobs))


####2. define residual diagnosis function####

residual_diagnosis <- function(model, data, model_name = "model") {
  res <- as.numeric(residuals(model))
  fit <- as.numeric(fitted(model))
  std_res <- as.numeric(scale(res))
  old_par <- par(mfrow = c(2, 2))
  on.exit(par(old_par), add = TRUE)
  plot(fit, res,
       xlab = "Fitted values", ylab = "Residuals",
       main = paste(model_name, "Residuals vs Fitted"),
       pch = 19, col = "grey40")
  abline(h = 0, lty = 2, col = "red")
  if (length(unique(fit)) > 1) {
    lines(lowess(fit, res), col = "blue", lwd = 2)
  }
  
  # (2) Q-Q plot:test residual normality
  qqnorm(std_res, main = paste(model_name, "Normal Q-Q"),
         pch = 19, col = "grey40")
  qqline(std_res, col = "red", lwd = 2)
  
  # (3) Scale-location
  plot(fit, sqrt(abs(std_res)),
       xlab = "Fitted values",
       ylab = expression(sqrt("|Standardized residuals|")),
       main = paste(model_name, "Scale-Location"),
       pch = 19, col = "grey40")
  if (length(unique(fit)) > 1) {
    lines(lowess(fit, sqrt(abs(std_res))), col = "blue", lwd = 2)
  }
  
  # (4) residual~x
  plot(data$x, res,
       xlab = "x", ylab = "Residuals",
       main = paste(model_name, "Residuals vs x"),
       pch = 19, col = "grey40")
  abline(h = 0, lty = 2, col = "red")
  if (length(unique(data$x)) > 1) {
    lines(lowess(data$x, res), col = "blue", lwd = 2)
  }
  
  # Shapiro-Wilk
  shapiro_result <- shapiro.test(res)

}

data_diag <- datah[complete.cases(datah[, c("x", "y")]), , drop = FALSE]

diag_linear <- residual_diagnosis(
  linear_model, data_diag, "Linear model"
)

diag_quadratic <- residual_diagnosis(
  quadratic_model, data_diag, "Quadratic model"
)

diag_linear_plateau <- residual_diagnosis(
  linear_plateau_model, data_diag, "Linear-plateau model"
)

diag_quadratic_plateau <- residual_diagnosis(
  quadratic_plateau_model, data_diag, "Quadratic-plateau model"
)


diagnostic_results <- rbind(
  diag_linear,
  diag_quadratic,
  diag_linear_plateau,
  diag_quadratic_plateau
)
print(diagnostic_results)


####for cold extremes####
dataset<-na.omit(data[,c('MAP','MAT','coldr')])
colnames(dataset)[3]<-'response'
zerovalue<-data.frame(x=NA,y=NA,sd=NA,cnt=NA)
for (i in 1:50){
  findzero<-which(dataset$MAP>100*(i-1)&dataset$MAP<=100*i)
  datazero<-dataset[findzero,c('MAT','response')]
  datazero<-na.omit(datazero)
  if (length(findzero)>20)
  {
    tmp<-array(NA)
    cnt=0
    for (iseed in 1:100){
      set.seed(iseed)
      x_seed<-sample(1:length(findzero),size=floor(length(findzero)*0.9))
      dataseed<-datazero[x_seed,]
      lmzero<-lm(response~MAT,data=dataseed)
      pzero<-summary(lmzero)$coefficients[,4]
      if (pzero[2]<0.05)
      {
        tmp[iseed]<-lmzero$coefficients[1]/lmzero$coefficients[2]*(-1)
        cnt<-cnt+1
      }
    }
    zerovalue[i,3]<-sd(tmp,na.rm=TRUE)
    zerovalue[i,1]<-i
    zerovalue[i,2]<-mean(tmp,na.rm=TRUE)
    zerovalue[i,4]<-cnt
  }
}

dataset<-na.omit(data[,c('MAP','MAT','coldr')])
colnames(dataset)[3]<-'response'
zerovalue<-data.frame(x=NA,y=NA,p=NA,p_inc=NA,k=NA,std=NA,SD_y=NA,uniq=NA)
for (i in 1:50){
  findzero<-which(dataset$MAP>100*(i-1)&dataset$MAP<=100*i)
  datazero<-dataset[findzero,c('MAT','response')]
  datazero<-na.omit(datazero)
  zerovalue[i,8]<-length(unique(datazero$MAT))
  if (length(findzero)>20)
  {
    lmzero<-lm(response~MAT,data=datazero)
    pzero<-summary(lmzero)$coefficients[,4]
    zerovalue[i,3:6]<-c(pzero[2],pzero[1],summary(lmzero)$coefficients[2,1:2])
    if (pzero[2]<0.05)
    {
      zerovalue[i,1]<-i
      zerovalue[i,2]<-lmzero$coefficients[1]/lmzero$coefficients[2]*(-1)
      lm_xy <- lm(MAT~response,data=datazero)
      a<-predict(lm_xy,newdata=data.frame(response=0),se.fit = TRUE)
      zerovalue[i, 7]<-a$se.fit
    }
  }
}
datac <- data.frame(x = zerovalue$x, y = zerovalue$y)
datac<-na.omit(datac)

linear_model <- lm(y ~ x, data = datac)
summary(linear_model)
pred <- predict(linear_model)
rmse <- sqrt(mean((datac$y - pred)^2))
mae  <- mean(abs(datac$y - pred))


quadratic_model <- lm(y ~ x + I(x^2), data = datac)
summary(quadratic_model)

pred <- predict(quadratic_model)
rmse <- sqrt(mean((datac$y - pred)^2))
mae  <- mean(abs(datac$y - pred))

linear_plateau <- function(x, a, b, c) {
  pmin(a + b * x, a + b * c)  # 当 x >= c 时，y = a + b * c（高原）
}

linear_plateau_model <- nls(y ~ linear_plateau(x, a, b, c),
                            data = datac,
                            start = list(a = -13, b = 3,
                                         c = 8))#-5,1,10
library(rcompanion)

plotPredy(data  = datac,
          x     = x,
          y     = y,
          model = linear_plateau_model)

a<-zerovalue$y
a<-a[!is.na(a)]
cor(fitted(linear_plateau_model),a)^2
coef(linear_plateau_model)["a"]+coef(linear_plateau_model)["b"]*coef(linear_plateau_model)["c"]

pred <- predict(linear_plateau_model)
rmse <- sqrt(mean((datac$y - pred)^2))
mae  <- mean(abs(datac$y - pred))

fit.lm<-lm(y ~ poly(x,2, raw=TRUE),data=datac)
a_parameter<-fit.lm$coefficients[1]
b_parameter<-fit.lm$coefficients[2]
c_parameter<-fit.lm$coefficients[3]
x_mean<-mean(datac$x)


qp <- function(x, a, b, jp) {
  c <- -0.5 * b / jp
  if_else(condition = x < jp,
          true  = a + (b * x) + (c * x * x),
          false = a + (b * jp) + (c * jp * jp))
}
# 使用 nls 拟合
library(dplyr)
quadratic_plateau_model <- nls(y ~ qp(x, a, b, jp),
                               data = datac,
                               start = list(a = a_parameter, b = b_parameter, jp = x_mean))#,nls.control(maxiter=1000

summary(quadratic_plateau_model)

coef(quadratic_plateau_model)["b"]/2/coef(quadratic_plateau_model)["jp"]*(-1) #c value
coef(quadratic_plateau_model)["a"]+coef(quadratic_plateau_model)["b"]*coef(quadratic_plateau_model)["jp"]/2

nullfunct = function(x, m){m}
m.ini= mean(datac$y)
null = nls(y ~ nullfunct(x, m),
           data = datac,
           start = list(m   = m.ini),
           trace = FALSE)
pred <- predict(quadratic_plateau_model)
rmse <- sqrt(mean((datac$y - pred)^2))
mae  <- mean(abs(datac$y - pred))

###  Find p-value and pseudo R-squared
nagelkerke(quadratic_plateau_model,null)

a<-zerovalue$y
a<-a[!is.na(a)]
cor(fitted(quadratic_plateau_model),a)^2 

datac$pred_linear <- predict(linear_model)
datac$pred_quadratic <- predict(quadratic_model)
datac$pred_linear_plateau <- predict(linear_plateau_model)
datac$pred_quadratic_plateau <- predict(quadratic_plateau_model)

p_models<-ggplot(datac, aes(x = x, y = y)) +
  geom_point(color = "black", alpha = 0.5) +
  geom_line(aes(y = pred_linear, color = "Linear"), size = 1) +
  geom_line(aes(y = pred_quadratic, color = "Quadratic"), size = 1) +
  geom_line(aes(y = pred_linear_plateau, color = "Linear-Plateau"), size = 1) +
  geom_line(aes(y = pred_quadratic_plateau, color = "Quadratic-Plateau"), size = 1) +
  labs(title = "",
       x = " Binned MAP(mm)", y = expression('T'[threshold])) +
  scale_color_manual(values = c("Linear" = "blue", "Quadratic" = "red",
                                "Linear-Plateau" = "green", "Quadratic-Plateau" = "purple")) +
  theme_minimal()



##residual diagnostics
data_diag <- datah[complete.cases(datac[, c("x", "y")]), , drop = FALSE]

diag_linear <- residual_diagnosis(
  linear_model, data_diag, "Linear model"
)

diag_quadratic <- residual_diagnosis(
  quadratic_model, data_diag, "Quadratic model"
)

diag_linear_plateau <- residual_diagnosis(
  linear_plateau_model, data_diag, "Linear-plateau model"
)

diag_quadratic_plateau <- residual_diagnosis(
  quadratic_plateau_model, data_diag, "Quadratic-plateau model"
)


diagnostic_results <- rbind(
  diag_linear,
  diag_quadratic,
  diag_linear_plateau,
  diag_quadratic_plateau
)
print(diagnostic_results)
