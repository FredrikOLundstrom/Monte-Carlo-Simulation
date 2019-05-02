#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------
#A Monte Carlo simulation and drawing the efficient frontier of portfolios containing stocks
#McDonalds, Bank of America, IBM, Chevron, Coca-Cola, Novartis and AT&T
#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------
# IMPORT DATA
data<-read.table("closes.dat.txt", col.names = c("McDonalds", "BankofAmerica", "IBM", "Chevron", "CocaCola", "Novartis", "ATnT"))
#--------------------------------------------------------------------------------------------------
# Monte Carlo Simulation
#--------------------------------------------------------------------------------------------------
# Here a Monte Carlo Simulation is used to find the efficient frontier 
# and to be able to select the optimal portfolio with corresponding alphas
# Expected return E[Ri]=1/T*SumRi(t), A measure of risk Var[Ri]=1/T*Sum(1:T) (Ri(t)-E[Ri])
Returns <- function(y){
  y_returns<-c(0, diff(y)/y[-length(y)])
  return(y_returns)
}
return_McDonalds <- Returns(data$McDonalds) 
return_BankofAmerica <- Returns(data$BankofAmerica)
return_IBM <- Returns(data$IBM)
return_Chevron <- Returns(data$Chevron)
return_CocaCola <- Returns(data$CocaCola)
return_Novartis <- Returns(data$Novartis)
return_ATnT <- Returns(data$ATnT)
return_STOCKS <- matrix(c(return_McDonalds, return_BankofAmerica, return_IBM, return_Chevron, return_CocaCola, return_Novartis, return_ATnT), nrow = 252, ncol = 7)
avg_return_STOCKS <- matrix(c(mean(return_McDonalds), mean(return_BankofAmerica), mean(return_IBM), mean(return_Chevron), mean(return_CocaCola), mean(return_Novartis), mean(return_ATnT)), nrow = 1, ncol = 7)
# RETURNS
cov_matrix <- cov(return_STOCKS, y = NULL, use = "all.obs")
# COVARIANCE MATRIX OF RETURNS
Monte_Carlo <- function(i_max=100, short_sell, avg_return_STOCKS){
  result_list_x = c()
  result_list_y = c()
  result_list_xx = c()
  for(i in 1:i_max){
    w = c()
    set.seed(i)
    if(short_sell){ w = runif(n = 7,min = -1, max = 1) }
    else { w = runif(n = 7) }
    w = w/sum(w)
    weighted_stock_returns <- avg_return_STOCKS * w
    portfolio_expected_returns <- sum(weighted_stock_returns)
    portfolio_variance <- t(w) %*% cov_matrix %*% w
    portfolio_standard_deviation <- portfolio_variance^(.5)
    result_list_x = c(result_list_x, portfolio_variance)
    result_list_xx = c(result_list_xx, portfolio_standard_deviation)
    result_list_y = c(result_list_y, portfolio_expected_returns)
  }
  return(list(result_list_x, result_list_y, result_list_xx))
}
#--------------------------------------------------------------------------------------------------
# Analytical Expression
#--------------------------------------------------------------------------------------------------
# Here the analytical expression is used to find the optimal portfolio and effient frontier
# and to be able to select the optimal portfolio with corresponding alphas
C <- cov_matrix
C_inv <- solve(C)
vector_of_ones <- matrix(rep(1, 7), nrow = 7, ncol = 1)
my <- t(avg_return_STOCKS)#set of returns
# a, b, c & d are scalar
a = as.numeric( t(vector_of_ones) %*% C_inv %*% vector_of_ones )
b = as.numeric( t(vector_of_ones) %*% C_inv %*% my )
c = as.numeric( t(my) %*% C_inv %*% my )
d = a * c - b^2
Markowitz <- function(my_p, my, a, b, c, d){
  result_efficient_frontier = c()
  for (x in my_p){
    lambda1 <- (c - b * x)/d
    lambda2 <- (a * x - b)/d
    result_efficient_frontier <-c(result_efficient_frontier, ((lambda1 + lambda2*x)^0.5) )
  }
  return(result_efficient_frontier)
}
#--------------------------------------------------------------------------------------------------
# Running the functions
#--------------------------------------------------------------------------------------------------
result = c()
result <- Monte_Carlo(10000, F, avg_return_STOCKS)
res_y <- unlist(result[2])
res_xx <- unlist(result[3])
result_2 = c()
result_2 <- Monte_Carlo(30000, T, avg_return_STOCKS)
res_y_2 <- unlist(result_2[2])
res_xx_2 <- unlist(result_2[3])
my_p <- seq(-1, 1, 0.00001)
variance_efficient_frontier <- Markowitz(my_p, my, a, b, c, d)
#--------------------------------------------------------------------------------------------------
# Minumum risk portfolio
minimum_risk_portfolio <- min(variance_efficient_frontier)
minimum_risk_portfolio
x_min <- sqrt(1/a)
y_min <- b/a

lambda1_efficient <- (c - b * y_min)/d
lambda2_efficient <- (a * y_min - b)/d

w_efficient <- C_inv %*% (lambda1_efficient*vector_of_ones + lambda2_efficient*my)

#--------------------------------------------------------------------------------------------------
# Plotting
#--------------------------------------------------------------------------------------------------
plot(res_xx, res_y, pch = 20, col = "BLUE", ylim = c(-6.4e-4, 6e-4), xlim = c(0.006, 0.015), xlab = "Risk (Standard Deviation)", ylab = "Expected Return", main = "Monte Carlo Simulation, with Markowitz's frontier")
points(res_xx_2, res_y_2, pch = 20, col = "BLACK")
points(res_xx, res_y, pch = 20, col = "BLUE")
lines(variance_efficient_frontier, my_p, col="RED", type = "l", lwd = 2)
points(x_min, y_min, pch = 21, cex = 2, bg = "RED",col = "BLACK")
legend_titles = c("Porfolios", "Portfolios allowing short-sell", "Markowitz's efficent frontier", "The minimum risk portfolio")
col_for_lines = c("BLUE","BLACK", "RED", "BLACK")
legend("topleft", legend = legend_titles, col = col_for_lines, 
       lty = c(0, 0, 1, 0), 
       lwd = c(0, 0, 2, 0), 
       pch = c(20, 20, NA, 21),
       pt.bg = c("BLUE","BLACK", NA, "RED"),
       pt.cex = 2)




