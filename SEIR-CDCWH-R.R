#PCL implementation of SEIR for CDC@SZ
#Authors: WH Li, XY Wei
#Jan 27, 2020

remove (list = objects() )
library (deSolve) 

#the function for SEIR
seir_model = function (current_timepoint, state_values, parameters)
{
  # create state variables (local variables)
  S = state_values [1]        # susceptibles
  E = state_values [2]        # exposed
  I = state_values [3]        # infectious
  R = state_values [4]        # recovered
  
  with ( 
    as.list (parameters),     # variable names within parameters can be used 
    {
      # compute derivatives
      dS = (-beta * S * I)
      dE = (beta * S * I) - (delta * E)
      dI = (delta * E) - (gamma * I)
      dR = (gamma * I)
      
      # combine results
      results = c (dS, dE, dI, dR)
      list (results)
    }
  )
}


#set parameters
contact_rate = 5                     # number of contacts per day
transmission_probability = 0.04       # transmission probability
infectious_period = 14                # infectious period
latent_period = 10                    # latent period

beta_value = contact_rate * transmission_probability
#beta_value = 0.8
gamma_value = 1 / infectious_period
delta_value = 1 / latent_period
Ro = beta_value / gamma_value

parameter_list = c (beta = beta_value, gamma = gamma_value, delta = delta_value)

#Set initial S,E,I,R values
W = 40000        # susceptible hosts
X = 1000           # infectious hosts
Y = 500           # recovered hosts
Z = 2000           # exposed hosts

N = W + X + Y + Z

initial_values = c (S = W/N, E = X/N, I = Y/N, R = Z/N)

timepoints = seq (0, 200, by=1)

output = lsoda (initial_values, timepoints, seir_model, parameter_list)

# plot (I * N ~ time, data = output, type='b', col = 'blue', ylab = 'Infectious', main = 'Wuhan infectious')

# susceptible hosts over time plot (S * N ~ time, data = output, type='l', col = 'blue', ylab = 'S, E, I, R', main = 'SEIR epidemic') 
plot (S * N ~ time, data = output, type='l', ylim = c(0,50000), col = 'blue', ylab = 'S, E, I, R', main = '武汉管制变化预测') 
# remain on same frame
par (new = TRUE)    
# exposed hosts over time
plot (E * N ~ time, data = output, type='l', ylim = c(0,50000), col = 'pink', ylab = '', axes = FALSE)
# remain on same frame
par (new = TRUE) 
# infectious hosts over time
plot (I * N ~ time, data = output, type='l', ylim = c(0,50000), col = 'red', ylab = '', axes = FALSE) 
# remain on same frame
par (new = TRUE)  
# recovered hosts over time
plot (R * N ~ time, data = output, type='l', ylim = c(0,50000), col = 'green', ylab = '', axes = FALSE)
# remain on same frame
par (new = TRUE) 
legend(x=150,y=20000,legend=c("Susceptible","Exposed","Infectious","Recovered"),col=c("blue","pink","red","green"), lty=1, cex=0.8)
