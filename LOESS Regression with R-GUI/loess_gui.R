#==========================
# Import Library  
#==========================
library(shiny)
library(shinythemes)
library(readxl)
library(DT)
library(dplyr)
library(forecast)
library(tseries)
library(MASS)

#==========================
# Inputasi Data  
#==========================
cekdata <- function(dt, sep){
  if (endsWith(dt$name, ".csv")) {
    return(read.csv(dt$datapath, sep = sep)) 
  }
  else if (endsWith(dt$name, ".txt")) { 
    return(read.table(dt$datapath, sep = sep, header = TRUE))
  }
  else if (endsWith(dt$name, ".xls") || endsWith(dt$name, ".xlsx")) { 
    return(read_excel(dt$datapath))
  }
}

#==========================
# Regresi Polinomial  
#==========================
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

#==========================
# Regresi LOESS
#==========================
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

#================
# Evaluasi Model
#================
model_eval = function(y,y_pred,degree) {
  n = length(y)
  k = degree + 1
  #R-Squared
  ss_total = sum((y - mean(y))^2)
  ss_residual = sum((y - y_pred)^2)
  r_squared = (1 - (ss_residual / ss_total)) * 100
  r_sq_adj = 1 - (((1 - r_squared) * (length(y) - 1)) / (length(y) - degree - 1))
  r_squared = round(r_squared, 3)
  r_sq_adj = round(r_sq_adj, 3)
  #MAPE
  mape = mean(abs((y - y_pred) / y)) * 100
  mape = round(mape, 3)
  #sMAPE
  smape = mean(2 * abs(y - y_pred) / (abs(y) + abs(y_pred))) * 100
  smape = round(smape, 3)
  #MSE
  mse = mean((y - y_pred)^2)
  mse = round(mse,3)
  #AIC
  aic = n * log(mse) + 2 * k
  aic = round(aic, 3)
  
  return(list(r_squared = r_squared, 
              r_squared_adj = r_sq_adj, 
              mape = mape, 
              smape = smape, 
              aic = aic,
              mse = mse))
}

#======================
# Invers Stasioneritas
#======================
inverse_stasionerity=function(y_pred,y_trans,diff_level,lambda){
  m = length(y_pred)
  y_pred_ori = numeric(m)
  
  #Inverse Differencing
  for (j in 1:m) {
    if (diff_level == 2) {
      y_pred_ori[j] = y_pred[j]+2*y_trans[j+1]-y_trans[j]  
    } else if (diff_level == 1) {
      y_pred_ori[j] = y_pred[j]+y_trans[j]
    } else if (diff_level == 0) {
      y_pred_ori[j] = y_pred[j]
    }
  }
  
  #Inverse Transformasi Box-Cox
  if (lambda != 1) {
    y_pred_orig <- ((y_pred_ori*lambda)+1)^(1/lambda)
  } else {
    y_pred_orig <- y_pred_ori
  }
  
  return(y_pred_orig)
}

#======================
# Komponen UI
#======================
ui <- fluidPage(
  theme = shinytheme("lumen"),
  titlePanel("PEMODELAN REGRESI LOCALLY ESTIMATED SCATTERPLOT SMOOTHING (LOESS) UNIVARIABEL DENGAN BEBERAPA METODE OPTIMASI (R-squared, R-squared Adjusted, AIC, MSE)",),
  tags$p("Oleh Faoza Akbar, dkk.", style = "font-size: 20px; color: black"),
  navbarPage("GUI-R",
             tabPanel("Pendahuluan", icon = icon("house"),
                      tabsetPanel(type = "pills", id = "navbar",
                                  tabPanel("Halaman Depan", icon = icon("user"), hr(),
                                           h2("Pemodelan Regresi LOCALLY ESTIMATED SCATTERPLOT SMOOTHING (LOESS)", style = "text-align: center;font-weight: bold"),
                                           h2("Dengan Beberapa Metode Optimasi (R-squared, R-squared Adjusted, AIC, MSE)", style = "text-align: center;font-weight: bold"),
                                           br(),
                                           img(src = "https://upload.wikimedia.org/wikipedia/id/thumb/2/2d/Undip.png/180px-Undip.png",
                                               width = "270px", height = "316.5px", style = "display:block;margin:0 auto;"),
                                           br(),
                                           h5("Inventor :", style = "text-align:center"),
                                           h5("Faoza Akbar (24050121140113)", style = "text-align:center"),
                                           h5("Dra. Suparti, M.Si. (196509131990032001)", style = "text-align:center"),
                                           h5("Puspita Kartikasari, S.Si., M.Si. (199105212019032021)", style = "text-align:center"),
                                           br(),
                                           h4("DEPARTEMEN STATISTIKA", style = "text-align:center"),
                                           h4("FAKULTAS SAINS DAN MATEMATIKA", style = "text-align:center"),
                                           h4("UNIVERSITAS DIPONEGORO", style = "text-align:center"),
                                           h4("2025", style = "text-align:center")
                                  ),
                                  tabPanel("Petunjuk", icon = icon("info-circle"), hr(),
                                           h3("Panduan Penggunaan GUI-R"),
                                           h4("Aplikasi ini khusus diperuntukkan menganalisis Regresi Locally Estimated Scatterplot Smoothing (LOESS) dengan 1 variabel prediktor dan 1 variabel respon."),
                                           h4("Analisis hanya dapat dilakukan untuk data time series."), br(),
                                           h4("1. Mengunggah Data"),
                                           h5("- Unggah file data pada panel 'Data' dengan format CSV, XLSX, XLS, atau TXT."),
                                           h5("- Pilih separator yang sesuai untuk file data. Pilihan yang tersedia adalah koma (,), titik koma (;), tab (\t), atau spasi ( )."),
                                           h5("- Pada bagian 'Data Time Series' pilih data time series (Zt) yang akan digunakan dalam analisis data."),
                                           h5("- Pada bagian 'Periode Data' pilih data yang berisi periode (t) data time series yang berupa hari, bulan, tahun, atau berupa angka dari 1 hingga banyaknya data."),
                                           h5("- Pada bagian 'Pilih Variabel X' pilih berapa periode sebelumnya yang akan digunakan sebagai variabel X. Pilihan yang tersedia 1 periode sebelumnya (Zt-1), 2 periode sebelumnya (Zt-2), atau 3 periode sebelumnya (Zt-3)."),
                                           h5("- Tentukan rasio untuk membagi data menjadi data in-sample dan out-sample."),
                                           h5("- Klik tombol 'Jalankan' untuk memuat data."),
                                           h5("- Data times series dengan periodenya akan ditampilkan pada tabel yang tersedia di subpanel 'Tabel Data'."),
                                           h5("- Plot data time series dapat dilihat pada subpanel 'Plot Data'."),
                                           h5("- Statistik deskriptif dari data time series dapat dilihat pada subpanel 'Statistik Deskriptif'."),
                                           h4("2. Pemodelan Data In-Sample"),
                                           h5("- Pada panel 'Pemodelan Data In-Sample' berisi hasil analisis data time series untuk data in-sample."),
                                           h5("- Masukkan minimum span, maksimal span, dan penambahan span untuk parameter span yang dibutuhkan dalam metode optimasi."),
                                           h5("- Pilih parameter degree yang akan digunakan dalam analisis. Pilihan yang tersedia adalah 1, 2, 3, 4, atau 5."),
                                           h5("- Pilih evaluasi yang digunakan untuk pengurutan hasil regresi terbaik dalam metode optimasi. Pilihan yang tersedia adalah R-squared, R-squared Adjusted, AIC, atau MSE."),
                                           h5("- Klik tombol 'Jalankan' untuk melakukan analisis data untuk data in-sample."),
                                           h5("- Pada subpanel 'Stasioneritas', ditampilkan hasil uji stasioneritas terhadap varians dan mean beserta plotnya."),
                                           h5("- Pada subpanel 'Variabel X dan Y', ditampilkan hasil tabel data yang berisi data X (Periode Sebelumnya) dan Y (Periode Sekarang) beserta plotnya."),
                                           h5("- Pada subpanel 'Optimasi', ditampilkan hasil span dan degree optimal beserta tabel data yang berisi hasil optimasi berdasarkan span dan degree yang diujikan."),
                                           h5("- Pada subpanel 'Estimasi', ditampilkan hasil estimasi model dengan span dan degree optimal. Hasil yang tertera adalah nilai MAPE dan SMAPE serta tabel perbandingan aktual dan prediksi."),
                                           h5("- Pada subpanel 'Plot Prediksi', ditampilkan hasil plot prediksi dari hasil estimasi menggunakan span dan degree optimal."),
                                           h4("3. Evaluasi Data Out-Sample"),
                                           h5("- Pada panel 'Evaluasi Data Out-Sample', berisi hasil dari evaluasi data out-sample berdasarkan span dan degree optimal."),
                                           h5("- Pada subpanel 'Evaluasi', ditampilkan nilai MAPE dan SMAPE serta tabel perbandingan aktual dan prediksi."),
                                           h5("- Pada subpanel 'Plot Prediksi', ditampilkan plot perbandingan aktual dan prediksi."),
                                           h4("4. Prediksi Beberapa Periode ke Depan"),
                                           h5("- Pada panel 'Prediksi Beberapa Periode ke Depan', berisi hasil prediksi beberapa periode ke depan setelah data time series yang digunakan."),
                                           h5("- Masukkan seberapa banyak periode ke depan yang akan ditampilkan pada hasil prediksi."),
                                           h5("- Klik tombol 'Jalankan' untuk melihat hasil prediksi."),
                                           h5("- Subpanel 'Hasil prediksi' menampilkan hasil dari prediksi beberapa periode ke depan.")
                                           
                                  )
                      )
             ),
             tabPanel("Data", icon = icon("upload"),
                      sidebarLayout(
                        sidebarPanel(
                          fileInput("file", "Input File Here:", multiple = FALSE, accept = c(".csv", ".xlsx", ".xls", ".txt")),
                          radioButtons(inputId = 'sep', label = 'Separator', choices = c(Comma = ',', Semicolon = ';', Tab = '\t', Space = ' '), selected = ','),
                          uiOutput("pilih"),
                          sliderInput(inputId = "ratio", label = "Rasio In-Sample dan Out-Sample", min = 0.5, max = 0.9, value = 0.8, step = 0.05),
                          actionButton(inputId = "run_data", label = "Jalankan")
                        ),
                        mainPanel(
                          tabsetPanel(type = "tabs",
                                      tabPanel("Tabel Data",br(),DTOutput("datahead")),
                                      tabPanel("Plot Data",plotOutput("plot_data")),
                                      tabPanel("Statistik Deskriptif",br(),DTOutput("statdesk"))
                          )
                        )
                      )
             ),
             tabPanel("Pemodelan Data In-Sample", icon = icon("chart-area"),
                      sidebarLayout(
                        sidebarPanel(
                          numericInput(inputId = "sf", label = "Minimial Span", value = 0.1, min = 0, max = 1, step = 0.01),
                          numericInput(inputId = "st", label = "Maksimal Span", value = 0.9, min = 0, max = 1, step = 0.01),
                          numericInput(inputId = "sb", label = "Penambahan Span", value = 0.01, min = 0.01, max = 0.1, step = 0.01),
                          sliderInput(inputId = "p", label = "Derajat", min = 1, max = 5, value = 1, step = 1),
                          selectInput(inputId = "sort_by",label = "Diurutkan Berdasarkan",choices=c("R-squared"="Rsq","R-squared Adjusted"="Radj","Akaike Information Criterion"="AIC","Mean Squared Error"="MSE")),
                          actionButton(inputId = "optimize", "Jalankan")
                        ),
                        mainPanel(
                          tabsetPanel(type = "tabs",
                                      tabPanel("Stasioneritas",
                                               br(),h4("Hasil Stasioneritas terhadap Varians dengan Uji Box-Cox"),
                                               verbatimTextOutput("lambda_text"),br(),
                                               h4("Hasil Stasioneritas terhadap Mean dengan Uji ADF"),
                                               verbatimTextOutput("pvalue_text"),br(),
                                               DTOutput('table_stas'),br(),
                                               plotOutput("plot_stas")
                                      ),
                                      tabPanel("Variabel X dan Y", 
                                               br(),DTOutput('varxy'), br(),
                                               plotOutput("plot_xy")
                                      ),
                                      tabPanel("Optimasi",
                                               br(),verbatimTextOutput("hasil_span"),br(),
                                               DTOutput('tabel_optimasi')
                                      ),
                                      tabPanel("Estimasi",
                                               br(),verbatimTextOutput('evalin'), br(),
                                               DTOutput("tabel_estimasi")
                                      ),
                                      tabPanel("Plot Prediksi",
                                               br(),
                                               plotOutput('plotpred'), br(),
                                               plotOutput('plotpredinv')
                                      )
                          )
                        )
                      )
             ),
             tabPanel("Evaluasi Data Out-Sample", icon = icon("chart-area"),
                        mainPanel(
                          tabsetPanel(type = "tabs",
                                      tabPanel("Evaluasi", 
                                               br(),verbatimTextOutput('evalout'),br(),
                                               DTOutput('table_evalout')
                                      ),
                                      tabPanel("Plot Prediksi",
                                               br(),plotOutput('plotevalout')      
                                      )
                          )
                        )
             ),
             tabPanel("Prediksi Beberapa Periode ke Depan", icon = icon("chart-line"),
                      sidebarLayout(
                        sidebarPanel(
                          numericInput(inputId = "predakt", label = "Prediksi Berapa Periode", value = 10, min = 1, step = 1),
                          actionButton(inputId = "predict", "Jalankan")
                        ),
                        mainPanel(
                          tabsetPanel(type = "tabs",
                                      tabPanel("Hasil Prediksi",
                                               br(),DTOutput('table_pred')
                                      )
                          )
                        )
                      )
             )
  )
)

#==========================
# Server 
#==========================
server <- function(input, output) {
  datapl <- reactive({
    req(input$file)
    cekdata(input$file, input$sep)
  })
  opts <- reactive({
    req(datapl())
    unique(colnames(datapl()))
  })
  
  output$pilih <- renderUI({
    req(datapl())
    tagList(
      selectInput(inputId = "datats", label = "Data Time Series", choices = opts(), selected = opts()[2]),
      selectInput(inputId = "period", label = "Periode Data", choices = opts(), selected = opts()[1]),
      selectInput(inputId = "hub", label = "Pilih Variabel X", choices = c("1 Periode Sebelumnya" = "cone", "2 Periode Sebelumnya" = "ctwo",                                                                    "3 Periode Sebelumnya" = "ctree"))
    )
  })
  
  data = reactiveVal(NULL)
  data_in_sample <- reactiveVal(NULL)
  data_out_sample <- reactiveVal(NULL)
  prop_in <- reactiveVal(NULL)
  minhub <- reactiveVal(NULL)
  x_ori <- reactiveVal(NULL)
  y_ori <- reactiveVal(NULL)
  x_stas = reactiveVal(NULL)
  y_stas = reactiveVal(NULL)
  N = reactiveVal(NULL)
  periode = reactiveVal(NULL)
  x.out_trans = reactiveVal(NULL)
  y.out_trans = reactiveVal(NULL)
  x.out_stas = reactiveVal(NULL)
  y.out_stas = reactiveVal(NULL)
  y.out_pred_f = reactiveVal(NULL)
  span_insample = reactiveVal(NULL)
  degree_insample = reactiveVal(NULL)
  r_f = reactiveVal(NULL)
  lambda_before_f = reactiveVal(NULL)
  diff_level_f = reactiveVal(NULL)
  
  observeEvent(input$run_data, {
    req(datapl(), input$datats, input$period, input$hub)
   
    if (input$hub == "cone") {
    minhub(1)  
    } else if (input$hub == "ctwo") {
    minhub(2)
    } else if (input$hub == "ctree") {
    minhub(3)
    } 
    
    data(datapl()[[input$datats]])
    periode(datapl()[[input$period]])
    data(as.numeric(data()))
    periode(as.character(periode()))
    prop_in(input$ratio) 
    N(length(data()))
    data_in_sample(data()[1:round(prop_in()*N())])  
    k = length(data_in_sample())
    x_ori(data_in_sample()[1:(k-minhub())])
    y_ori(data_in_sample()[(1+minhub()):k])
    
    output$datahead <- renderDT({
      req(datapl())
      datatable(datapl())
    })
    
    output$plot_data <- renderPlot({
      req(datapl(), input$datats, input$period)
      
      n <- length(periode())
      if (n > 10) { selected_indices <- seq(1, n, length.out = 10)  # Ambil sekitar 10 titik yang merata
      } else { selected_indices <- seq_along(periode()) }
      avg_length <- mean(nchar(periode()[selected_indices]))
      if (avg_length > 10) { line_pos <- 9  
      } else { line_pos <- 7 }
      par(mar = c(line_pos + 2, 4, 4, 2) + 0.1)
      
      plot(seq_along(data()), data(), type = "o", col = "blue", pch = 20, lwd = 1,
           xaxt = "n", xlab = "", ylab = "Data Time Series", main = "Plot Data Time Series")
      axis(1, at = selected_indices, labels = periode()[selected_indices], las = 2, cex.axis = 0.8)
      mtext("Periode", side = 1, line = line_pos, cex = 1, adj = 0.5)
      grid(col = "gray", lty = "dotted")
    })
    
    output$statdesk <- renderDT({
      data.frame(
       Mean = round(mean(data(), na.rm = TRUE),3),
       Varians = round(var(data(), na.rm = TRUE),3),
       Std.Dev = round(sd(data(), na.rm = TRUE),3),
       Minimum = round(min(data(), na.rm = TRUE),3),
       Maksimum = round(max(data(), na.rm = TRUE),3)
        )
    })
  }) #end observe
  
  observeEvent(input$optimize, {
    req(data_in_sample(),minhub(),x_ori(),y_ori(),input$sf,input$st,input$sb,input$p,input$sort_by,N(),prop_in())
    
    #uJi stasioneritas dalam varians (box cox test)
    k = length(data_in_sample())
    lambda_before <- BoxCox.lambda(data_in_sample())  # Simpan lambda sebelum transformasi
    lambda_before_f (lambda_before)
    
    if (lambda_before != 1) {
      data_in_sample (((data_in_sample()^lambda_before) - 1)/lambda_before)
      lambda_after <- BoxCox.lambda(data_in_sample())  # Simpan lambda setelah transformasi
    }
    x_trans = data_in_sample()[1:(k-minhub())]
    y_trans = data_in_sample()[(1+minhub()):k]
    
    output$lambda_text <- renderPrint({
      cat("Lambda sebelum transformasi:", round(lambda_before, 3), "\n")
      if (lambda_before != 1) {
      cat("Lambda setelah transformasi:", round(lambda_after))
      }
    })
    
    #Uji Stasioneritas dalam Mean (ADF test)
    reactive_data <- reactiveValues(data = data_in_sample())
    
    output$pvalue_text <- renderPrint({
      initial_data <- reactive_data$data
      p_value_print <- adf.test(initial_data)$p.value  
      diff_level <- 0
      cat("P-value sebelum differencing:", round(p_value_print, 3), "\n")
      while (p_value_print > 0.05) {
        initial_data <- diff(initial_data, differences = 1)  # Perbarui data
        p_value_print <- suppressWarnings(adf.test(initial_data)$p.value)  # Update p-value
        diff_level <- diff_level + 1  
        cat("Differencing Level", diff_level, "dengan P-value:", round(p_value_print, 3), "\n")
      }
    })
    
    p_value = adf.test(data_in_sample())$p.value
    diff_level = 0
    while (p_value > 0.05) {
      data_in_sample(diff(data_in_sample(),differences = 1))
      p_value <- suppressWarnings(adf.test(data_in_sample())$p.value)  
      diff_level <- diff_level + 1
    }
    diff_level_f(diff_level)

    #plot data stasioner
    l=length(data_in_sample())
    t=seq_along(data_in_sample())
    output$plot_stas <- renderPlot({
      plot(t,data_in_sample(),col="blue",pch=20,lwd=2,
           xlab="t",ylab="Data Stasioner",main = "Plot Data Stasioner")
      grid(col = "gray", lty = "dotted")
    })
    
    #pendefinisian variabel prediktor dan respon
    x_stas(data_in_sample()[1:(l-minhub())])
    y_stas(data_in_sample()[(1+minhub()):l])
    
    output$varxy <- renderDT({
      datatable(data.frame(
        Variabel_X = round(x_stas(),3),
        Variabel_Y = round(y_stas(),3)
      ),colnames=c("X (Periode Sebelumnya)","Y (Periode Sekarang)"))
    })
   
    #plot data x dan y yang Sudah Stasioner
    output$plot_xy <- renderPlot({
      plot(x_stas(),y_stas(),col="blue",pch=20,lwd=2,xlab="X",ylab="Y",main="Plot Data X dan Y")
      grid(col = "gray", lty = "dotted")
    })
    
    results <- data.frame(Span=numeric(),R_Squared=numeric(),R2_Adj=numeric(),
                          AIC=numeric(),MSE=numeric(),stringsAsFactors=FALSE)
    span_seq <- seq(input$sf, input$st, input$sb)
    degree <- input$p
    degree_insample(degree)
    
    for (i in span_seq) {
      model_in <- loess_reg(x = x_stas(), y = y_stas(), a = i, degree = degree)
      y_pred <- model_in$fitted_values
      
      y_pred_orig <- inverse_stasionerity(y_pred, y_trans, diff_level, lambda_before)
      
      n <- length(x_ori())
      eval <- model_eval(y_ori()[(diff_level+1):n], y_pred_orig, degree)
      
    results <- rbind(results,data.frame(
        Span = i * 100,
        R_Squared = as.numeric(eval$r_squared),
        R2_Adj = as.numeric(eval$r_squared_adj),
        AIC = as.numeric(eval$aic),
        MSE = as.numeric(eval$mse)
     ))
    }
    
    #urutkan berdasarkan metrik yang dipilih
    sorting_col <- switch(input$sort_by,
                          "Rsq" = "R_Squared",
                          "Radj" = "R2_Adj",
                          "AIC" = "AIC",
                          "MSE" = "MSE")
    
    if (sorting_col %in% c("R_Squared", "R2_Adj")) {
      sorted_results <- results[order(-results[[sorting_col]]), ]
    } else {
      sorted_results <- results[order(results[[sorting_col]]), ]
    }
    
    span_opt <- sorted_results[1, 1] / 100
    
    output$hasil_span <- renderPrint({
      cat("Span optimal =", span_opt, "dengan degree =", degree)
    })
    
    output$tabel_optimasi <- renderDT({ datatable(sorted_results,colnames = c("Span (%)", "R-squared (%)", "R-squared Adjusted (%)", "AIC", "MSE")) }) 
    
    #estimasi
    model_in = loess_reg(x=x_stas(),y=y_stas(),a=span_opt,degree=degree)
    y_pred=model_in$fitted_values
    y_pred_orig=inverse_stasionerity(y_pred, y_trans, diff_level, lambda_before)
    
    span_insample(model_in$span_in)
    
    eval_insample = model_eval(y_ori()[(diff_level+1):n],y_pred_orig,degree)
    
    output$evalin = renderPrint({
      cat("Nilai MAPE data in-sample sebesar",as.numeric(eval_insample$mape))
      cat("\nNilai SMAPE data in-sample sebesar",as.numeric(eval_insample$smape))
    })
    output$tabel_estimasi <- renderDT({ datatable(data.frame(round(y_pred,3),round(y_pred_orig,3),y_ori()[(diff_level+1):n]),colnames = c("Prediksi","Prediksi Setelah Inverse Stasioneritas", "Aktual")) }) 
    
    output$plotpred = renderPlot({
      plot((2+diff_level):(length(y_pred)+1+diff_level), y_pred,col = "red",
           pch = 20, lty = 1,xlab = "t", ylab = "Data Stasioner",main="Plot Prediksi untuk Data In-Sample")
      grid(col = "gray", lty = "dotted")
    })
    
    output$plotpredinv = renderPlot({
    plot((2+diff_level):(length(y_pred_orig)+1+diff_level),
         y_ori()[(1+diff_level):length(y_ori())],col = "blue", pch = 20,
         lty = 1,xlab = "t", ylab = "Data Setelah Inverse Stasioneritas",main="Plot Prediksi untuk Data In-Sample Setelah Inverse Stasioneritas",
         ylim=c((min(y_pred_orig)-1),(max(y_pred_orig)+1)))
    lines((2+diff_level):(length(y_pred_orig)+1+diff_level), y_pred_orig, type = "o",
          col = "red", pch = 20, lty = 1)
    grid(col = "gray", lty = "dotted")
    legend("topleft", legend = c("Aktual", "Prediksi"),
           col = c("blue", "red"), pch = 20, lty = 1)
    })
    
    #evaluasi data out sample
    data_out_sample(data()[round(prop_in()*N()-diff_level):N()]) 
    
    o = length(data_out_sample())
    x.out_ori = data_out_sample()[1:(o-minhub())]
    y.out_ori = data_out_sample()[(1+minhub()):o]
    
    if (lambda_before!=1){
      data_out_sample(((data_out_sample()^lambda_before)-1)/lambda_before)
    }
    x.out_trans(data_out_sample()[1:(o-minhub())])
    y.out_trans(data_out_sample()[(1+minhub()):o])
    
    if (diff_level!=0){
      data_out_sample(diff(data_out_sample(),differences = diff_level))
    }
    p=length(data_out_sample())
    x.out_stas(data_out_sample()[1:(p-minhub())])
    y.out_stas(data_out_sample()[(1+minhub()):p])
    
    g=length(y_pred)
  
    x.out_stas_uji=append(x_stas(),x.out_stas()[1])
    y.out_stas_uji=append(y_stas(),y_pred[g])
    
    z=length(x.out_stas())
    y.out_pred=numeric(z)
    for (c in 2:(z+1)){
      model_out=loess_reg(x=x.out_stas_uji,y=y.out_stas_uji,
                          b=model_in$span_in,degree=degree)
      y.out_pred_temp=model_out$fitted_values
      
      x.out_stas_uji[g+1]=x.out_stas()[c]
      y.out_stas_uji[g+1]=tail(y.out_pred_temp,1)
      y.out_pred[(c-1)]=tail(y.out_pred_temp,1)
    }
    y.out_pred_f(y.out_pred)
    
    y.out_pred_orig = inverse_stasionerity(y.out_pred,y.out_trans(),
                                           diff_level,lambda_before)
    
    r=length(x.out_ori)
    r_f(r)
    eval_outsample=model_eval(y.out_ori[(diff_level+1):r],
                              y.out_pred_orig,degree)
    
    output$evalout = renderPrint({
      cat("Nilai MAPE data out-sample sebesar",as.numeric(eval_outsample$mape))
      cat("\nNilai SMAPE data out-sample sebesar",as.numeric(eval_outsample$smape))
    })
    
    #data aktual vs prediksi
    compare=data.frame(Periode=tail(periode(),length(y.out_pred_orig)),
                       Aktual=y.out_ori[(diff_level+1):r],Prediksi=round(y.out_pred_orig,3))
    
    output$table_evalout=renderDT({datatable(compare)})
    
    #plot aktual vs prediksi
    output$plotevalout=renderPlot({
      par(mar = c(8, 4, 4, 2) + 0.1) 
      plot(1:length(compare$Aktual), compare$Aktual, col = "blue",
           pch = 20, lty = 1,xlab = "", ylab = "Data  Out-Sample",main="Plot Aktual vs Prediksi Data Out-Sample",xaxt = "n",
           ylim=c((min(y_pred_orig)-1),(max(y_pred_orig)+1)))
      lines(1:length(compare$Prediksi), compare$Prediksi, type = "o",
            col = "red", pch = 20, lty = 1)
      axis(1,at=1:length(compare$Prediksi),labels=compare$Periode,
           las = 2, cex.axis = 0.8)
      grid(col = "gray", lty = "dotted")
      mtext("Periode", side = 1, line = 7, cex = 1.2)
      legend("topright", legend = c("Aktual", "Prediksi"),
             col = c("blue", "red"), pch = 20, lty = 1)
    })
  }) #end observe
  
  
  observeEvent(input$predict, {
    req(input$predakt)
    
    #PREDIKSI BEBERAPA PERIODE KE DEPAN
    x_future=append(x_stas(),tail(y.out_stas(),1)) 
    y_future=append(y_stas(),tail(y.out_pred_f(),1)) 
    
    preddepan=input$predakt
    y_future_pred=numeric(preddepan)
    f=length(x_future)
    for (s in 1:preddepan) {
      model = loess_reg(x=x_future,y=y_future,
                        b=span_insample(),degree_insample())
      y_future_pred_temp = model$fitted_values
      x_future[f]=tail(y_future,1)
      y_future[f]=tail(y_future_pred_temp,1)
      y_future_pred[s]=tail(y_future_pred_temp,1)
    }
    
    y_future_pred_ori=numeric(preddepan)
    
    if (diff_level_f() == 2){
      y_future_pred_ori[1]=y_future_pred[1]+2*y.out_trans()[r_f()]-y.out_trans()[(r_f()-1)]
      y_future_pred_ori[2]=y_future_pred[2]+2*y_future_pred_ori[1]-y.out_trans()[r_f()]
      for (u in 3:preddepan){
        y_future_pred_ori[u]=y_future_pred[u]+2*y_future_pred_ori[u-1]-y_future_pred_ori[u-2] 
      } 
    } 
    else if (diff_level_f() == 1){
      y_future_pred_ori[1]=y_future_pred[1]+y.out_trans()[(r_f()-1)]
      for (u in 2:preddepan){
      y_future_pred_ori[u]=y_future_pred[u]+y_future_pred_ori[u-1]
      }
    } 
    else if (diff_level_f() == 0){
      y_future_pred_ori[u]=y_future_pred[u]
    }
    
    if (lambda_before_f() != 1){
    y_future_pred_orig=((y_future_pred_ori*lambda_before_f())+1)^(1/lambda_before_f())
    } else {
    y_future_pred_orig=y_future_pred_ori
    }
    
    output$table_pred= renderDT({
      datatable(data.frame(Prediksi=round(y_future_pred_orig,3)))
    })
    
  }) #end observe
  
}

#==========================
# Jalankan Aplikasi 
#==========================
shinyApp(ui = ui, server = server)

