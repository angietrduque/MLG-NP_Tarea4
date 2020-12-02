#----------------------------------------------------------------------------------------#
# MLGNP Tarea 3
# Cesar A Saavedra
# Angie Rodriguez uque
#----------------------------------------------------------------------------------------#
suppressMessages(library(dplyr))
suppressMessages(library(tidyverse))
suppressMessages(library(factoextra))
suppressMessages(library(corrplot))
suppressMessages(library(polycor))
suppressMessages(library(psych))
suppressMessages(library(gplots))
suppressMessages(library(gridExtra))
suppressMessages(library(nlme))
suppressMessages(library(MASS))
suppressMessages(library(reshape))
suppressMessages(library(homals))
suppressMessages(library(broom))
suppressMessages(library(readr))
suppressMessages(library(lubridate))
suppressMessages(library(purrr))
suppressMessages(library(VGAM))
#----------------------------------------------------------------------------------------#
setwd("/Users/cesar.saavedra/Documents/GitHub/MLG-NP_Tarea4")
#----------------------------------------------------------------------------------------#
Datos <- read.table("Datos.txt",header=T,sep = ",")
Datos
#----------------------------------------------------------------------------------------#
n <- 60
set.seed(12345)
muestra <- Datos %>% sample_n(size=n,replace=FALSE)
muestra <- muestra %>% arrange(pH)
muestra
#----------------------------------------------------------------------------------------#
x <- muestra  %>% dplyr::select(fixed.acidity, pH)
x
ggplot() + geom_point(data = x, aes(x = pH, y = fixed.acidity))
#----------------------------------------------------------------------------------------#
n = nrow(x)
y  = pull(x, fixed.acidity)
sigma.rice <- 1/(2*(n-1))*sum((y - lag(y, k=1))^2, na.rm = T)
sigma.rice
#----------------------------------------------------------------------------------------#
base_cons <- function(x,j){
  sqrt(2)*cos((j-1)*pi*x)
}

lambda.select <- function(x, lambda, salida=1){
  df <-   dplyr::select(x, fixed.acidity)
  hora_normada <- x$pH; i <- list(); y <- list()
  
  for (i in 2:lambda) {
    y[[i]] <- base_cons(hora_normada, i)
  }
  
  y <- data.frame(matrix(unlist(y), ncol = lambda-1))
  df <-  bind_cols(df, y)
  
  f.i <- df %>% 
    dplyr::select(contains("x")) %>% 
    colnames()
  
  lambda <- lambda %>%
    tibble(lambda = . )
  
  f.i_sum <- paste(f.i, collapse = "+")
  mdl_formula <- as.formula(paste("fixed.acidity", f.i_sum, sep = "~"))
  
  mdl <- lm(mdl_formula, data = df)
  
  fitted <- predict(mdl) %>%
    tibble(fitted = . )
  
  S = lm.influence(mdl)$hat
  tr = sum(S)
  
  UBRE <-  (1/n) * sum(resid(mdl)^2) + (2/n) * sigma.rice*tr - sigma.rice %>% 
    tibble(ubre = .)
  
  CV <- (1/n) * sum(((residuals(mdl) / (1-S))^2)) %>%
    tibble(cv = .)
  
  GCV <- (1/n) * ( sum(resid(mdl)^2) / (1 - (1/n) * tr)^2 ) %>% 
    tibble(GCV = .)
  
  R <- bind_cols(UBRE, CV, GCV, lambda)
  
  if(salida == 1) {
    return(R)
  } else {
    return(fitted)
  } 
}

all.R <- function(x, lambda){
  R <- list()
  for(i in 2:lambda){
    R[[i]] <- lambda.select(x,i)
  }
  R <- as.data.frame(t(matrix(unlist(R), ncol = lambda-1)))
  names(R) <- c("UBRE","CV","GCV","LAMBDA")
  return(R)
}
#----------------------------------------------------------------------------------------#
lambda <- 30
all.R <- all.R(x, lambda)
all.R
UBRE <- cbind(all.R$UBRE,all.R$LAMBDA)
UBRE <- data.frame(UBRE)
names(UBRE)[1] = "UBRE";names(UBRE)[2] = "LAMBDA"
UBRE
#----------------------------------------------------------------------------------------#
plot1 <- ggplot()+
  geom_point(data = all.R, aes(x = LAMBDA, y = UBRE)) +
  theme_bw() +
  labs(x =  expression(lambda), y = expression(hat(R)(lambda)))

plot2 <- ggplot()+
  geom_point(data = all.R, aes(x = LAMBDA, y = CV)) +
  theme_bw() +
  labs(x =  expression(lambda), y = expression(hat(CV)(lambda)))

plot3 <- ggplot()+
  geom_point(data = all.R, aes(x = LAMBDA, y = GCV)) +
  theme_bw() +
  labs(x =  expression(lambda), y = expression(hat(GCV)(lambda)))

grid.arrange(plot1, plot2, plot3, ncol=2)
#----------------------------------------------------------------------------------------#
lambda <- 3
salida <- 2 # 1: Riesgo, 2: fitted values

fitted <- lambda.select(x, lambda, salida)
x <- bind_cols(x, fitted)
x
#----------------------------------------------------------------------------------------#
ggplot()+ geom_point(data = x, aes(x = pH, y = fixed.acidity)) +
  geom_line(data = x, aes(x =pH, y = fitted), col="red") +
  labs(subtitle = expression(lambda==3))
#----------------------------------------------------------------------------------------#
library(readxl)
library(dplyr)
library(lubridate)
library(ggplot2)
library(purrr)

# Cargar los datos
Datos2 <- read_excel("OzonoCompartir2019.xlsx")
Datos2  


x <- Datos2  %>% 
  dplyr::select(`Fecha & Hora`, Datos2$O3) %>% 
  mutate(dia = day(`Fecha & Hora`)) %>%
  filter(dia %in% c(9,10,11,12,13)) %>%  
  group_by(dia) %>%
  
  mutate(hora = row_number()) %>% 
  ungroup() %>%
  
  mutate(hora_normada = (2*hora - 1)/(2*24)) %>%
  mutate(dia = as.factor(dia))
#----------------------------------------------------------------------------------------#
ggplot() +
  geom_line(data = x, aes(x = hora_normada, y = O3, color = dia))+
  labs(y = expression(O[3]))

ggplot()+
  geom_point(data = x, aes(x = hora_normada, y = O3)) +
  labs(y = expression(O[3]))+ facet_wrap(~dia)
#----------------------------------------------------------------------------------------#
base_cons <- function(x,j){
  sqrt(2)*cos((j-1)*pi*x)
}

df <-  x %>% 
  dplyr::select(O3) 

hora_normada <- x$hora_normada; lambda <- 7; i <- NULL; aux <- NULL; y <- NULL

for (i in 2:lambda) {
  aux <- base_cons(hora_normada, i)
  y[[i]] <- aux
}

y <- data.frame(matrix(unlist(y), ncol = lambda-1)) # Variables predictoras
df <-  bind_cols(df, y)
#----------------------------------------------------------------------------------------#
lm_grupo <- function(x){
  model <- lm(O3 ~  X1+X2+X3+X4+X5+X6, data = x)
  fitted_NO2 <- predict(model) %>% 
    as_tibble
  bind_cols(x, fitted_NO2)
}

df <- df  %>% 
  mutate(Dia = as.factor(x$dia)) %>%
  group_split(Dia) %>% 
  map(.x = ., .f = lm_grupo) %>% 
  bind_rows()
#----------------------------------------------------------------------------------------#
fitted <- as.data.frame(df$value)
names (fitted) = "fitted"
df = bind_cols(x, fitted)

#----------------------------------------------------------------------------------------#
ggplot()+
  geom_point(data = df, aes(x = hora_normada, y = O3)) +
  geom_line(data = df, aes(x =hora_normada, y = fitted)) +
  labs(y = expression(O[3])) +
  labs(title = "Estimacion con base de cosenos") +
  labs(subtitle = expression(lambda==)) +
  facet_wrap(~dia)


#----------------------------------------------------------------------------------------#