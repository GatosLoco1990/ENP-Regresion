library(dplyr)
library(fpc)
library(GSE)

# 
# library(ggplot2)
# library(tidyverse)
# library(corrplot)
# library(mlbench)
# library(caret)
# library(lattice)

#CSV generado desde python
datos <- read.csv("C:/Users/usuario/Desktop/abalone/abalone.csv", header = FALSE)
head(datos)
names(datos) <- c("Sex", "Length", "Diameter", "Height", 
                  "Whole_Weight", "Shucked_Weight", "Viscera_Weight", 
                  "Shell_Weight", "Rings")
#Quitar primer dato que está mal por la conversión
datos <- datos[-1, ]

datos$Sex <- factor(datos$Sex, levels = c("M", "F", "I"))

# Aplicar One-Hot Encoding a la columna 'Sex' usando model.matrix
onehot_sex <- model.matrix(~ Sex - 1, data = datos)
onehot_sex
# Combinar las nuevas columnas de One-Hot Encoding con los datos originales
datos_final <- cbind(onehot_sex, datos[, -1])  # Eliminamos la columna 'Sex' original
datos_final
summary(datos_final)

#Datos genericos

# Establecer la semilla para reproducibilidad
#set.seed(123)

# Definir el tamaño de la muestra
n <- 500

# Generar variables aleatorias normales
x1 <- rnorm(n, mean = 1, sd = 1)
x2 <- rnorm(n, mean = 2, sd = 2)
x3 <- rnorm(n, mean = 3, sd = 3)
e <- rnorm(n, mean = 0, sd = sqrt(2))

# Calcular y
y <- 1.20 + 9.76 * x1 + 6.89 * x2 + 3.53 * x3 + e

# Crear una matriz para almacenar los datos originales
matr <- cbind(x1, x2, x3)

# Contaminar 9% de los individuos de la muestra
n_contaminados <- floor(n * 0.09)

contaminados <- sample(1:n, n_contaminados, replace = FALSE)

# Contaminar los datos seleccionados
for (i in contaminados) {
  matr[i, ] <- rnorm(3, mean = 22, sd = sqrt(1.5))
}

plot(matr[, 1],
     main = "Gráfico de dispersión de x1",
     xlab = "Índice",
     ylab = "Valor de x1",
     pch = 19,
     col = "blue")


# Función para aplicar DBSCAN a cada columna
apply_dbscan_to_column <- function(data_column, eps = 0.5, minPts = 2) {
  # Aplicar DBSCAN
  clusters <- dbscan(data_column, eps = eps, MinPts = minPts)
  
  # Devolver los labels
  return(clusters$cluster)
}

# Aplicar DBSCAN a cada columna
outliers_x1 <- apply_dbscan_to_column(matr[, 1])  # Parámetros ajustados según datos
outliers_x2 <- apply_dbscan_to_column(matr[, 2])
outliers_x3 <- apply_dbscan_to_column(matr[, 3])


convertir_valores <- function(labels) {
  # Encuentra el valor más repetido
  contar <- table(labels)
  # Convierte los valores que no son el más repetido a 1 y el más repetido a 0
  return(ifelse(labels != names(contar)[which.max(contar)], 1, 0))
}

outliers_x1_n <- convertir_valores(outliers_x1)
outliers_x2_n <- convertir_valores(outliers_x2)
outliers_x3_n <- convertir_valores(outliers_x3)

# Añadir resultados a los datos originales
matr_outliers <- cbind(matr, outliers_x1 = outliers_x1_n, 
                       outliers_x2 = outliers_x2_n, 
                       outliers_x3 = outliers_x3_n)

# Visualizar los resultados
print("Resultados de DBSCAN en cada columna:")
head(matr_outliers)


# Gráfico de dispersión coloreado según outliers
plot(matr_outliers[, 1],
     
     col = ifelse(outliers_x1_n == 0, "blue", "red"),
     main = "x1, coloreado según outliers",
     xlab = "índice", ylab = "x1")

legend("right", legend = c("No outlier", "Outlier"),
       fill = c("blue", "red"))

# Crear una lista de las columnas de outliers
outliers_list <- list(outliers_x1_n, outliers_x2_n, outliers_x3_n)

# Usar un ciclo for para recorrer las columnas y asignar NaN donde corresponda
for (i in 1:3) {
  # Asignar NaN donde el valor del outlier es 1 en cada variable x1, x2, x3
  matr_outliers[outliers_list[[i]] == 1, i] <- NaN
}

# Visualizar los resultados actualizados
print("Matriz con NaN donde outliers es 1 en x1, x2, x3:")
head(matr_outliers)

# Eliminar las columnas de outliers

matr_outliers_f <- matr_outliers[, -c(4, 5, 6)]
matr_outliers_f

matr_outliers_f_y <- cbind(matr_outliers_f, y)
matr_outliers_f_y

resultados_GSE <- GSE(matr_outliers_f_y, tol = 1e-5, maxiter = 3000, method = "rocke", init = "qc")
sigma <- getScatter(resultados_GSE)

sigma_xx <- sigma[1:3, 1:3]
sigma_xy <- sigma[4, 1:3]

beta_sin_cero <- solve(sigma_xx) %*% sigma_xy
beta_sin_cero

matr_outliers_f_y

# Calcular el promedio de las primeras 3 columnas
mux <- colMeans(matr_outliers_f_y[, 1:3], na.rm = TRUE)
muy <- mean(matr_outliers_f_y[, 4], na.rm = TRUE)

beta_con_cero <- muy - t(mux) %*% beta_sin_cero
beta_con_cero

beta_final <- cbind(beta_con_cero, t(beta_sin_cero))
beta_final

calcular_betas <- function(matr_s, y_s){
  aba_outliers_x1 <- apply_dbscan_to_column(matr_s[, 1])  # Parámetros ajustados según datos
  aba_outliers_x2 <- apply_dbscan_to_column(matr_s[, 2])
  aba_outliers_x3 <- apply_dbscan_to_column(matr_s[, 3])
  
  aba_outliers_x1_n <- convertir_valores(aba_outliers_x1)
  aba_outliers_x2_n <- convertir_valores(aba_outliers_x2)
  aba_outliers_x3_n <- convertir_valores(aba_outliers_x3)
  
  
  aba_matr_outliers <- cbind(matr, aba_outliers_x1 = aba_outliers_x1_n, 
                             aba_outliers_x2 = aba_outliers_x2_n, 
                             aba_outliers_x3 = aba_outliers_x3_n)
  
  #print(aba_matr_outliers)
  
  # Crear una lista de nombres de las columnas de outliers
  outliers_list <- list(aba_outliers_x1_n, aba_outliers_x2_n, aba_outliers_x3_n)
  
  
  # Usar un ciclo for para recorrer las columnas y asignar NaN donde corresponda
  for (k in 1:3) {
    # Asignar NaN donde el valor del outlier es 1 en cada variable x1, x2, x3
    print(length(table(outliers_list[k])))
    if (length(table(outliers_list[k])) > 1){
      print(aba_matr_outliers[outliers_list[[k]]] == 1)
      aba_matr_outliers[outliers_list[[k]] == 1, k] <- NaN
    }
  }
  
  aba_matr_outliers_f <- aba_matr_outliers[, -c(4, 5, 6)]
  aba_matr_outliers_f_y <- cbind(aba_matr_outliers_f, y_s)
  
  resultados_GSE <- GSE(aba_matr_outliers_f_y, tol = 1e-5, maxiter = 500, method = "rocke", init = "emve")
  sigma <- getScatter(resultados_GSE)
  
  sigma_xx <- sigma[1:3, 1:3]
  sigma_xy <- sigma[4, 1:3]
  
  beta_sin_cero <- solve(sigma_xx) %*% sigma_xy
  
  # Calcular el promedio de las primeras 3 columnas
  mux <- colMeans(aba_matr_outliers_f_y[, 1:3], na.rm = TRUE)
  muy <- mean(aba_matr_outliers_f_y[, 4], na.rm = TRUE)
  
  beta_con_cero <- muy - t(mux) %*% beta_sin_cero
  
  beta_final <- cbind(beta_con_cero, t(beta_sin_cero))
  return(beta_final)
}


# Prueba del algoritmo

max_iter <- 100
tamaños <- cbind(100, 500, 1000)


resultados_totales <- list()
MSE = list()
MAE = list()
AIC = list()
BIC = list()

for (j in 1:length(tamaños)){
  resultados <- matrix(nrow = max_iter, ncol = 4, data = NA)
  mse <- 0
  mae <- 0
  aic_total  <- 0
  bic_total  <- 0
  for (i in 1:max_iter){
    n <- tamaños[j]
    x1 <- rnorm(n, mean = 1, sd = 1)
    x2 <- rnorm(n, mean = 2, sd = 2)
    x3 <- rnorm(n, mean = 3, sd = 3)
    e <- rnorm(n, mean = 0, sd = sqrt(2))
    
    # Calcular y
    y <- 1.20 + 9.76 * x1 + 6.89 * x2 + 3.53 * x3 + e
    
    # Crear una matriz para almacenar los datos originales
    matr <- cbind(x1, x2, x3)
    
    # Contaminar 9% de los individuos de la muestra
    n_contaminados <- floor(n * 0.09)
    
    contaminados <- sample(1:n, n_contaminados, replace = FALSE)
    
    # Contaminar los datos seleccionados
    for (k in contaminados) {
      matr[k, ] <- rnorm(3, mean = 22, sd = sqrt(1.5))
    }
    
    # Aplicar DBSCAN a cada columna
    outliers_x1 <- apply_dbscan_to_column(matr[, 1])  # Parámetros ajustados según datos
    outliers_x2 <- apply_dbscan_to_column(matr[, 2])
    outliers_x3 <- apply_dbscan_to_column(matr[, 3])
    
    outliers_x1_n <- convertir_valores(outliers_x1)
    outliers_x2_n <- convertir_valores(outliers_x2)
    outliers_x3_n <- convertir_valores(outliers_x3)
    
    matr_outliers <- cbind(matr, outliers_x1 = outliers_x1_n, 
                           outliers_x2 = outliers_x2_n, 
                           outliers_x3 = outliers_x3_n)
    
    # Crear una lista de nombres de las columnas de outliers
    outliers_list <- list(outliers_x1_n, outliers_x2_n, outliers_x3_n)
    
    # Usar un ciclo for para recorrer las columnas y asignar NaN donde corresponda
    for (k in 1:3) {
      # Asignar NaN donde el valor del outlier es 1 en cada variable x1, x2, x3
      matr_outliers[outliers_list[[k]] == 1, k] <- NaN
    }
    
    matr_outliers_f <- matr_outliers[, -c(4, 5, 6)]
    matr_outliers_f_y <- cbind(matr_outliers_f, y)
    
    resultados_GSE <- GSE(matr_outliers_f_y, tol = 1e-5, maxiter = 3000, method = "rocke", init = "emve_c")
    sigma <- getScatter(resultados_GSE)
    
    sigma_xx <- sigma[1:3, 1:3]
    sigma_xy <- sigma[4, 1:3]
    
    beta_sin_cero <- solve(sigma_xx) %*% sigma_xy
    
    # Calcular el promedio de las primeras 3 columnas
    mux <- colMeans(matr_outliers_f_y[, 1:3], na.rm = TRUE)
    muy <- mean(matr_outliers_f_y[, 4], na.rm = TRUE)
    
    beta_con_cero <- muy - t(mux) %*% beta_sin_cero
    
    beta_final <- cbind(beta_con_cero, t(beta_sin_cero))
    y_est <- beta_final[1] + beta_final[2] * x1 + beta_final[3] * x2 + beta_final[4] * x3 + e 
    error <- y - y_est
    
    # Cálculo de MSE y MAE
    error_mse <- mean(error**2)
    error_mae <- mean(abs(error))
    mse <- mse + error_mse
    mae <- mae + error_mae
    
    # RSS para AIC y BIC
    RSS <- sum(error^2)
    k <- 4  # Número de parámetros (intercepto + 3 coeficientes)
    
    # Cálculo de AIC y BIC
    aic <- 2 * k + n * log(RSS / n)
    bic <- k * log(n) + n * log(RSS / n)
    
    aic_total <- aic_total + aic
    bic_total <- bic_total + bic
    
    resultados[i, 1:4] <- beta_final
  }
  #Matriz de betas de todas las iteraciones
  resultados_totales<- append(resultados_totales, list(resultados))
  
  mse <- mse / max_iter
  mae <- mae / max_iter
  aic_prom <- aic_total / max_iter
  bic_prom <- bic_total / max_iter
  
  MSE <- append(MSE, mse)
  MAE <- append(MAE, mae)
  AIC <- append(AIC, aic_prom)
  BIC <- append(BIC, bic_prom)
  
  cat("terminado el tamaño ", tamaños[j], "\n")
}

b0 <- colMeans(resultados, na.rm = T)

cat(paste("Mae: ", paste(MAE, collapse = " "), "\n"))
cat(paste("Mse: ", paste(MSE, collapse = " "), "\n"))
cat(paste("Aic: ", paste(AIC, collapse = " "), "\n"))
cat(paste("Bic: ", paste(BIC, collapse = " "), "\n"))

metricas <- rbind(MAE, MSE,AIC,BIC)
`colnames<-`(metricas, c("100", "500", "1000"))


resultados_totales[3]

#Usando datos de Abalone

aba_x1 <- datos_final$Shell_Weight
aba_x2 <- datos_final$Height
aba_x3 <- datos_final$Whole_Weight

y_aba <- datos_final$Rings

# Crear una matriz para almacenar los datos originales
aba_matr <- cbind(aba_x1, aba_x2, aba_x3)
f_aba <- cbind(aba_matr, y_aba)

outliers_x1 <- apply_dbscan_to_column(aba_matr[, 1])  # Parámetros ajustados según datos
outliers_x2 <- apply_dbscan_to_column(aba_matr[, 2])
outliers_x3 <- apply_dbscan_to_column(aba_matr[, 3])

outliers_x1_n <- convertir_valores(outliers_x1)
outliers_x2_n <- convertir_valores(outliers_x2)
outliers_x3_n <- convertir_valores(outliers_x3)

matr_outliers <- cbind(aba_matr, outliers_x1 = outliers_x1_n, 
                       outliers_x2 = outliers_x2_n, 
                       outliers_x3 = outliers_x3_n)

# Crear una lista de nombres de las columnas de outliers
outliers_list <- list(outliers_x1_n, outliers_x2_n, outliers_x3_n)

# Usar un ciclo for para recorrer las columnas y asignar NaN donde corresponda
for (k in 1:3) {
  # Asignar NaN donde el valor del outlier es 1 en cada variable x1, x2, x3
  matr_outliers[outliers_list[[k]] == 1, k] <- NaN
}

matr_outliers_f <- matr_outliers[, -c(4, 5, 6)]
matr_outliers_f_y <- cbind(matr_outliers_f, y_boot)


sigma_1 <- getScatter(GSE(f_aba, tol = 1e-5, maxiter = 500, method = "rocke", init = "emve_c"))

sigma_xx <- sigma_1[1:3, 1:3]
sigma_xy <- sigma_1[4, 1:3]

beta_sin_cero <- solve(sigma_xx) %*% sigma_xy
# Calcular el promedio de las primeras 3 columnas
mux <- colMeans(matr_outliers_f_y[, 1:3], na.rm = TRUE)
muy <- mean(matr_outliers_f_y[, 4], na.rm = TRUE)

beta_con_cero <- muy - t(mux) %*% beta_sin_cero

beta_final <- cbind(beta_con_cero, t(beta_sin_cero))
beta_final
#Remuestreo para intervalos de confianza

n <- nrow(matr)
B <- 10000 # Número de bootstrap samples
#Se puede cambiar B por los N datos

betas_bootstrap <- matrix(nrow = B, ncol = length(c(beta_con_cero, t(beta_sin_cero))))
for (i in 1:B) {
  rm(outliers_x1, outliers_x2, outliers_x3,
     sample_indices, X_boot, y_boot,
     outliers_x1_n, outliers_x2_n, outliers_x3_n,
     matr_outliers, outliers_list, matr_outliers_f,
     matr_outliers_f_y, sigma_1, sigma_xx, sigma_xy,
     beta_sin_cero, mux, muy, beta_con_cero, beta_final)
  
  sample_indices <- sample(1:n, n-1, replace = TRUE)
  X_boot <- aba_matr[sample_indices, ]
  y_boot <- y_aba[sample_indices]
  
  # Aplicar DBSCAN a cada columna
  outliers_x1 <- apply_dbscan_to_column(X_boot[, 1])  # Parámetros ajustados según datos
  outliers_x2 <- apply_dbscan_to_column(X_boot[, 2])
  outliers_x3 <- apply_dbscan_to_column(X_boot[, 3])
  
  outliers_x1_n <- convertir_valores(outliers_x1)
  outliers_x2_n <- convertir_valores(outliers_x2)
  outliers_x3_n <- convertir_valores(outliers_x3)
  
  matr_outliers <- cbind(X_boot, outliers_x1 = outliers_x1_n, 
                         outliers_x2 = outliers_x2_n, 
                         outliers_x3 = outliers_x3_n)
  
  # Crear una lista de nombres de las columnas de outliers
  outliers_list <- list(outliers_x1_n, outliers_x2_n, outliers_x3_n)
  
  # Usar un ciclo for para recorrer las columnas y asignar NaN donde corresponda
  for (k in 1:3) {
    # Asignar NaN donde el valor del outlier es 1 en cada variable x1, x2, x3
    matr_outliers[outliers_list[[k]] == 1, k] <- NaN
  }
  
  matr_outliers_f <- matr_outliers[, -c(4, 5, 6)]
  matr_outliers_f_y <- cbind(matr_outliers_f, y_boot)
  
  sigma_1 <- getScatter(GSE(matr_outliers_f_y, tol = 1e-5, maxiter = 500, method = "rocke", init = "emve_c"))
  
  sigma_xx <- sigma_1[1:3, 1:3]
  sigma_xy <- sigma_1[4, 1:3]
  
  beta_sin_cero <- solve(sigma_xx) %*% sigma_xy
  # Calcular el promedio de las primeras 3 columnas
  mux <- colMeans(matr_outliers_f_y[, 1:3], na.rm = TRUE)
  muy <- mean(matr_outliers_f_y[, 4], na.rm = TRUE)
  
  beta_con_cero <- muy - t(mux) %*% beta_sin_cero
  
  beta_final <- cbind(beta_con_cero, t(beta_sin_cero))
  
  
  betas_bootstrap[i, ] <- beta_final
  
  if (i %% 10 == 0){cat("Estamos en: ", i, "\n")}
}

quantile(betas_bootstrap[,1], c(0.025, 0.975), na.rm = TRUE)

IC95 <- list()

for (i in 1:4) {
  # Calcula los cuantiles para la i-ésima columna
  intervalo <- quantile(betas_bootstrap[, i], c(0.025, 0.975), na.rm = TRUE)
  
  # Guarda el resultado en la lista IC95
  IC95[[i]] <- intervalo
  
  # Imprime el intervalo de confianza para la i-ésima columna
  cat("Intervalo de confianza del 95% para la variable", i, ":\n")
  print(intervalo)
}

df_IC95 <- do.call(rbind.data.frame, IC95)
colnames(df_IC95) <- c("Limite inferior", "Limite superior")
rownames(df_IC95) <- paste("Beta", 1:4)

# Imprime el data.frame resultante
print(df_IC95)

hist(betas_bootstrap[,4])
