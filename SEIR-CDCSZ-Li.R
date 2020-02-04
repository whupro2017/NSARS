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
contact_rate = 10                     # number of contacts per day
transmission_probability = 0.2       # transmission probability
infectious_period = 14                # infectious period
latent_period = 10                    # latent period

beta_value = contact_rate * transmission_probability
#beta_value = 0.8
gamma_value = 1 / infectious_period
delta_value = 1 / latent_period
Ro = beta_value / gamma_value

parameter_list = c (beta = beta_value, gamma = gamma_value, delta = delta_value)

#Set initial S,E,I,R values
W = 13026619        # susceptible hosts
X = 1           # infectious hosts
Y = 0           # recovered hosts
Z = 9           # exposed hosts

N = W + X + Y + Z

initial_values = c (S = W/N, E = X/N, I = Y/N, R = Z/N)

timepoints = seq (0, 200, by=1)

output = lsoda (initial_values, timepoints, seir_model, parameter_list)

# susceptible hosts over time
plot (S ~ time, data = output, type='l', ylim = c(0,1), col = 'blue', ylab = 'S, E, I, R', main = 'SEIR epidemic') 
# remain on same frame
par (new = TRUE)    
# exposed hosts over time
plot (E ~ time, data = output, type='l', ylim = c(0,1), col = 'pink', ylab = '', axes = FALSE)
# remain on same frame
par (new = TRUE) 
# infectious hosts over time
plot (I ~ time, data = output, type='l', ylim = c(0,1), col = 'red', ylab = '', axes = FALSE) 
# remain on same frame
par (new = TRUE)  
# recovered hosts over time
plot (R ~ time, data = output, type='l', ylim = c(0,1), col = 'green', ylab = '', axes = FALSE)
# remain on same frame
par (new = TRUE) 
legend(180,0.7,legend=c("S","E","I","R"),col=c("blue","pink","red","green"), lty=1, cex=0.8)

