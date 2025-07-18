#library
#install.packages("dplyr")
#install.packages("readxl")
#install.packages("forecast")
#install.packages("tseries")
#install.packages("MASS")
library(dplyr)
library(readxl)
library(forecast)
library(tseries)
library(MASS)

#regresi polinomial
poly_reg = function(x, y, degree, weight) {
  if (length(x) != length(y) || length(x) != length(weight)) {
    stop("x, y, dan bobot harus memiliki jumlah data yang sama")
  }
  #matriks x
  X = sapply(0:degree, function(d) x^d)
  #weighted least squares
  W = diag(weight)  # diagonal matriks untuk bobot
  beta = ginv(t(X) %*% W %*% X) %*% (t(X) %*% W %*% y)
  return(beta)
}

#regresi LOESS
loess_reg = function(x, y, a=0.5,b=NULL,degree=1) {
  n = length(x)
  fitted_values = numeric(n)
  if (is.null(b)){
    span=ceiling(a*n)
  } else {
    span=b
  }
  for (i in 1:n) {
    #jarak euclidean satu dimensi
    distances = abs(x - x[i])
    x = x - x[i]
    #dataframe yang diurutkan berdasarkan distances
    df_span = data.frame(x,y,distances)
    df_span_sorted = arrange(df_span,distances)
    #variabel yang sudah dipotong sesuai span
    x_span = df_span_sorted$x[1:span]
    y_span = df_span_sorted$y[1:span]
    distances_span = df_span_sorted$distances[1:span]
    #normalisasi dengan  simple feature scalling
    scaled_dist = distances_span/max(distances_span)
    #pembobotan tricube
    weight = ifelse(scaled_dist < 1, (1-(scaled_dist)^3)^3, 0)
    weight[weight < 0] = 0
    #pemodelan dengan regresi polinomial
    local_model = poly_reg(x_span,y_span,degree=degree,
                           weight=weight)
    #nilai prediksi
    fitted_values[i]=sum(local_model*((x[i]-x[i])^(0:degree)))
  }
  return(list(fitted_values=fitted_values,span_in=span))
}

#evaluasi model
model_eval = function(y,y_pred,degree) {
  #R-Squared
  ss_total = sum((y-mean(y))^2)
  ss_residual = sum((y-y_pred)^2)
  r_squared = (1-(ss_residual/ss_total))*100
  r_squared = round(r_squared,3)
  #MAPE
  mape = mean(abs((y-y_pred)/y))*100
  mape = round(mape,3)
  
  return(list(r_squared=r_squared,mape=mape))
}

#inverse stasionerity
inverse_stasionerity=function(y_pred,y_trans,diff_level,lambda){
  m = length(y_pred)
  y_pred_ori = numeric(m)
  
  #inverse differencing
  for (j in 1:m) {
    if (diff_level == 2) {
      y_pred_ori[j] = y_pred[j]+2*y_trans[j+1]-y_trans[j]  
    } else if (diff_level == 1) {
      y_pred_ori[j] = y_pred[j]+y_trans[j]
    } else if (diff_level == 0) {
      y_pred_ori[j] = y_pred[j]
    }
  }
  
  #inverse transformasi box-cox
  if (lambda != 1) {
    y_pred_orig <- ((y_pred_ori*lambda)+1)^(1/lambda)
  } else {
    y_pred_orig <- y_pred_ori
  }
  
  return(y_pred_orig)
}

#visualisasi data awal
#data
datapl = read_excel("D:/Pelajaran/Kuliah/Semester 7 (Skripsi)/Coding/Dataset/DataTA.xlsx", 
                    sheet = "Data Inflasi")
#plot data
par(mar = c(8, 4, 4, 2) + 0.1) 
plot(1:length(datapl$Periode), datapl$Data_Inflasi,
     col="blue",pch=20,lwd = 2,xaxt = "n", 
     xlab = "", ylab = "Data Inflasi (%)")
selected_indices = seq(1, length(datapl$Periode), by = 5)
axis(1, at = selected_indices, labels = datapl$Periode[selected_indices], las = 2, cex.axis = 0.8)
grid(col = "gray", lty = "dotted")
mtext("Periode", side = 1, line = 7, cex = 1.2)

###ANALISIS DATA
#data
data = datapl$Data_Inflasi
prop_in=0.80 #proporsi in-sample
N=length(data)
data_in_sample = data[1:round(prop_in*N)]  
data_out_sample = data[round(prop_in*N-2):N] 
k = length(data_in_sample)
x_ori = data_in_sample[1:k-1]
y_ori = data_in_sample[2:k]

#statistik deskriptif
desc_stats <- data.frame(
  Mean = mean(data),Varians = var(data),Std.Dev = sd(data),
  Minimum = min(data),Maksimum = max(data))
print(desc_stats)

##uji stasioneritas
#uJi stasioneritas dalam varians
lambda=BoxCox.lambda(data_in_sample)
round(lambda,3)
if (lambda!=1) {
  data_in_sample = ((data_in_sample^lambda)-1)/lambda
  x_trans = data_in_sample[1:k-1]
  y_trans = data_in_sample[2:k]
  round(BoxCox.lambda(data_in_sample))
} 
data_in_sample
#Uji Stasioneritas dalam Mean
p_value=adf.test(data_in_sample)$p.value
diff_level=0
while (p_value > 0.05) {
  cat("Differencing Level =",diff_level,";","P-Value =",p_value,"\n")
  data_in_sample = diff(data_in_sample, differences = 1)  
  p_value = suppressWarnings(adf.test(data_in_sample)$p.value)  
  diff_level=diff_level+1
  print(data_in_sample)
}
cat("Differencing Level =",diff_level,";","P-Value =",p_value,"\n")
#plot data stasioner
l=length(data_in_sample)
t=seq(1,l,1)
plot(t,data_in_sample,col="blue",pch=20,lwd=2,
     xlab="t",ylab="Data Stasioner")
grid(col = "gray", lty = "dotted")

#pendefinisian variabel prediktor dan respon
x_stas = data_in_sample[1:l-1]
y_stas = data_in_sample[2:l]
x_stas
y_stas

#plot data x dan y yang Sudah Stasioner
plot(x_stas,y_stas,col="blue",pch=20,lwd=2,xlab="X",ylab="Y")
grid(col = "gray", lty = "dotted")

##optimasi span dan degree
results = data.frame()
span = seq(0.10,0.95,0.05)
degree = 2

for (i in span) {
  model_in=loess_reg(x=x_stas,y=y_stas,a=i,degree=degree)
  y_pred=model_in$fitted_values
  
  y_pred_orig=inverse_stasionerity(y_pred,y_trans,
                                    diff_level,lambda)
  
  #evaluasi model data in-sample
  n=length(x_ori)
  eval = model_eval(y_ori[(diff_level+1):n],y_pred_orig,degree)
  
  results = rbind(results, data.frame(
    Span = i*100,
    R_Squared = eval$r_squared
  ))
}
sorted_results = results[order(-results$R_Squared), ]
span_opt=sorted_results[1,1]/100
cat("======================================= \n",
    "span optimal =",span_opt,"dengan degree =",degree,"\n",
    "======================================= \n")
print(sorted_results)

#pemodelan regresi LOESS dengan span dan degree optimal beserta visualisasinya
model_in = loess_reg(x=x_stas,y=y_stas,a=span_opt,degree=degree)
y_pred=model_in$fitted_values
y_pred

y_pred_orig=inverse_stasionerity(y_pred, y_trans, diff_level, lambda)
y_pred_orig

eval_insample = model_eval(y_ori[(diff_level+1):n],y_pred_orig,degree)
eval_insample$mape

plot((2+diff_level):(length(y_pred)+1+diff_level), y_pred,col = "red",
     pch = 20, lty = 1,xlab = "t", ylab = "Data Stasioner")
grid(col = "gray", lty = "dotted")

plot((2+diff_level):(length(y_pred_orig)+1+diff_level),
     y_ori[(1+diff_level):length(y_ori)],col = "blue", pch = 20,
     lty = 1,xlab = "t", ylab = "Data Inflasi (%)",
     ylim=c((min(y_pred_orig)-1),(max(y_pred_orig)+1)))
lines((2+diff_level):(length(y_pred_orig)+1+diff_level), y_pred_orig, type = "o",
      col = "red", pch = 20, lty = 1)
grid(col = "gray", lty = "dotted")
legend("topleft", legend = c("Aktual", "Prediksi"),
       col = c("blue", "red"), pch = 20, lty = 1)

##evaluasi model data out-sample dengan MAPE
o = length(data_out_sample)
x.out_ori = data_out_sample[1:o-1]
y.out_ori = data_out_sample[2:o]

data_out_sample = ((data_out_sample^lambda)-1)/lambda
x.out_trans = data_out_sample[1:o-1]
y.out_trans = data_out_sample[2:o]

data_out_sample = diff(data_out_sample,differences = diff_level)
p=length(data_out_sample)
x.out_stas = data_out_sample[1:p-1]
y.out_stas = data_out_sample[2:p]

g=length(y_pred)

x.out_stas_uji=append(x_stas,x.out_stas[1])
y.out_stas_uji=append(y_stas,y_pred[g])

z=length(x.out_stas)
y.out_pred=numeric(z)
for (c in 2:(z+1)){
  model_out=loess_reg(x=x.out_stas_uji,y=y.out_stas_uji,
                      b=model_in$span_in,degree=degree)
  y.out_pred_temp=model_out$fitted_values
  
  x.out_stas_uji[g+1]=x.out_stas[c]
  y.out_stas_uji[g+1]=tail(y.out_pred_temp,1)
  y.out_pred[(c-1)]=tail(y.out_pred_temp,1)
}

y.out_pred_orig = inverse_stasionerity(y.out_pred,y.out_trans,
                                        diff_level,lambda)

r=length(x.out_ori)
eval_outsample=model_eval(y.out_ori[(diff_level+1):r],
                          y.out_pred_orig,degree)
eval_outsample$mape

##visualisasi data aktual vs prediksi
#data frame aktual vs prediksi
compare=data.frame(Periode=tail(datapl$Periode,length(y.out_pred_orig)),
                   Aktual=y.out_ori[(diff_level+1):r],Prediksi=y.out_pred_orig)
compare
#plot aktual vs prediksi
par(mar = c(8, 4, 4, 2) + 0.1) 
plot(1:length(compare$Aktual), compare$Aktual, col = "blue",
     pch = 20, lty = 1,xlab = "", ylab = "Data Inflasi (%)",xaxt = "n",
     ylim=c((min(y_pred_orig)-1),(max(y_pred_orig)+1)))
lines(1:length(compare$Prediksi), compare$Prediksi, type = "o",
      col = "red", pch = 20, lty = 1)
axis(1,at=1:length(compare$Prediksi),labels=compare$Periode,
     las = 2, cex.axis = 0.8)
grid(col = "gray", lty = "dotted")
mtext("Periode", side = 1, line = 7, cex = 1.2)
legend("topright", legend = c("Aktual", "Prediksi"),
       col = c("blue", "red"), pch = 20, lty = 1)

##prediksi beberapa bulan ke depan
x_future=append(x_stas,tail(y.out_stas,1)) 
y_future=append(y_stas,tail(y.out_pred,1)) 

month=4
y_future_pred=numeric(month)
f=length(x_future)
for (s in 1:month) {
  model = loess_reg(x=x_future,y=y_future,
                    b=model_in$span_in,degree=degree)
  y_future_pred_temp = model$fitted_values
  x_future[f]=tail(y_future,1)
  y_future[f]=tail(y_future_pred_temp,1)
  y_future_pred[s]=tail(y_future_pred_temp,1)
}

y_future_pred_ori=numeric(month)
y_future_pred_ori[1]=y_future_pred[1]+2*y.out_trans[r]-y.out_trans[(r-1)]
y_future_pred_ori[2]=y_future_pred[2]+2*y_future_pred_ori[1]-y.out_trans[r]
  
for (u in 3:month){
  y_future_pred_ori[u]=y_future_pred[u]+2*y_future_pred_ori[u-1]-y_future_pred_ori[u-2] 
}
y_future_pred_orig=((y_future_pred_ori*lambda)+1)^(1/lambda)

data.frame(Prediksi=y_future_pred_orig)
