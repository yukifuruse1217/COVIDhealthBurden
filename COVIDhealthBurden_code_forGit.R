# このコードは、新規陽性者数の変化と、それに応じて「酸素投与を要する人」と「重症者」が生じる関係をモデル化したものです。
# ワクチンの接種率や有効性がどのように影響するのかや、必要な確保病床数の推定部分については、本コードに含まれていません。
# それらは予測ツールのEXCEL内で計算されています。EXCELにある非表示の計算式は「567567567」で保護を解除して表示させることができます。


require(deSolve) # for the "ode" function

#Equation#######################

sir_1 <- function(status1, status2, hospRate, sevRate, delta1, delta2,
                  I0, Ha0, Hb0, Hc0, Da0, Db0, Dc0,R0,
                  times) {
  
  # the differential equations:
  sir_equations <- function(time, variables, parameters) {
    with(as.list(c(variables, parameters)), {
      
      if (time<30) {
        dI <- I * (status1^(1/7)-1)
      } else {
        if (status2 == 5) {
          dI <- I * (status1^(1/7)-1)
        } else if (status2 == 6){
          dI <- 0
        } else if (status2 == 7){
          dI <- I * (0.85^(1/5)-1)
        }
      }
      
      dHa <- hospRate * I - delta1 * Ha
      dHb <- delta1 * Ha - delta1 * Hb
      dHc <- delta1 * Hb - delta1 * Hc
      
      dDa <- (sevRate / hospRate) * delta1 * Hc - delta2 * Da
      dDb <- delta2 * Da - delta2 * Db
      dDc <- delta2 * Db - delta2 * Dc
      
      dR <- (1 - (sevRate / hospRate)) * delta1 * Hc + delta2 * Dc
      
      return(list(c(dI, dHa, dHb, dHc, dDa, dDb, dDc, dR)))
    })
  }
  
  # the parameters values:
  parameters_values <- c(status1=status1, status2=status2, hospRate=hospRate, sevRate=sevRate, delta1=delta1, delta2=delta2)
  
  # the initial values of variables:
  initial_values <- c(I=I0, Ha=Ha0, Hb=Hb0, Hc=Hc0, Da=Da0, Db=Db0, Dc=Dc0, R=R0)
  
  # solving
  out <- ode(initial_values, times, sir_equations, parameters_values)

  #out <- ode(method="rk4",initial_values, times, sir_equations, parameters_values)
  
  # returning the output:
  simulationresult <- as.data.frame(out)
  return(simulationresult)
}



#Simulations

list_status1 <- c(1) #increasing rate per week
list_status2 <- c(7) #"same", "flat", "decrease2"

list_delta_check <- c(1)

list_I0 <- c(100)

list_hospRate <- c(0.15) #ex data for 50s
list_sevRate <- c(0.04) #ex data for 50s

list_delta1 <- c(10) #ex data for 50s
list_delta2 <- c(5) #ex data for 50s

i <- 0
for (status1_temp in list_status1) {
  for (status2_temp in list_status2) {
    for (I0_temp in list_I0) {
      for (hospRate_temp in list_hospRate) {
        for (sevRate_temp in list_sevRate) {
          for (delta_check in list_delta_check) {
          
        
      i = i+1
      
      delta1_temp <- list_delta1[delta_check]/3
      delta2_temp <- list_delta2[delta_check]/3
      
      result <- sir_1(status1=status1_temp, status2=status2_temp, hospRate=hospRate_temp, sevRate=sevRate_temp, delta1=1/delta1_temp, delta2=1/delta2_temp,
                      I0=1*I0_temp, Ha0=0, Hb0=0, Hc0=0, Da0=0, Db0=0, Dc0=0, R0=0,
                      times = seq(0, 60, by=1))

  var_status1 <- rep(status1_temp, 61)
  var_status2 <- rep(status2_temp, 61)
  var_I0 <- rep(I0_temp, 61)
  var_hospRate <- rep(hospRate_temp, 61)
  var_sevRate <- rep(sevRate_temp, 61)
  var_deltaCheck <- rep(delta_check, 61)

  result <- cbind(result,var_status1)
  result <- cbind(result,var_status2)
  result <- cbind(result,var_I0)
  result <- cbind(result,var_hospRate)
  result <- cbind(result,var_sevRate)
  result <- cbind(result,var_deltaCheck)

  if (i == 1) {
    comb_result <- result
  } else {
    comb_result <- rbind(comb_result,result)
  }
    
}}}}}}

hosp <- comb_result$Ha + comb_result$Hb + comb_result$Hc + comb_result$Da + comb_result$Db + comb_result$Dc
severe <- comb_result$Da + comb_result$Db + comb_result$Dc + (comb_result$Hc * sevRate_temp / hospRate_temp)

comb_result <- cbind(comb_result, hosp)
comb_result <- cbind(comb_result, severe)

comb_result[is.na(comb_result)] <- 0
