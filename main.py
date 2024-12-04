import numpy as np
import matplotlib.pyplot as plt
import time
import random
from statsmodels.tsa.stattools import adfuller
from tqdm import tqdm

# Set the parameters
Omega = 10
mu = 1
epsilon = 0.1
gamma = 0.2

# Set the initial condition
N_0 = 1

# Set the maximum simulation time
t_max = 100

# Plot the simulated trajectory# Initialize time and population arrays
start_time=time.time()
time_points = [0]
population = [N_0]



def test_stationarity(data):
    # Perform Augmented Dickey-Fuller test
    result = adfuller(data)
    
    # Extract test statistics and p-value
    test_statistic = result[0]
    p_value = result[1]
    
    # Set the significance level
    significance_level = 0.00001
    print(p_value)
    # Compare test statistics with critical values and p-value with significance level
    if test_statistic < result[4]['5%'] and p_value < significance_level:
        
        print("The process is likely stationary")
        return True
    else:
        print("The process is likely non-stationary")
        return False


key = False
# Implement the Gillespie algorithm
while time_points[-1] < t_max  and key ==False :
    N = population[-1]  # Current population

    # Calculate the birth and death rates
    birth_rate = mu * N + epsilon
    death_rate = gamma * N * N / Omega

    # Calculate the total rate
    total_rate = birth_rate + death_rate

    # Generate a random exponential time step
    delta_t = np.random.exponential(1 / total_rate)
    time_points.append(time_points[-1] + delta_t)

    # Determine the type of event (birth or death)
    if np.random.rand() < birth_rate / total_rate:
        population.append(N + 1)
    else:
        population.append(N - 1)

    if len(time_points)%1500 ==0:
        key = test_stationarity(population[-1500:])


# Convert the population list to a NumPy array
population = np.array(population)

# Calculate the mean and variance in the steady-state
steady_state_population = population[-1500:]
mean_population = np.mean(steady_state_population)
var_population = np.var(steady_state_population)

plt.plot(time_points, population)
plt.xlabel('Time')
plt.ylabel('N')
plt.title('Simulated Trajectory of N(t)')

plt.savefig('Traj.png')
plt.show()
print('Mean in steady-state:', mean_population)
print('Variance in steady-state:', var_population)

end_time = time.time()
print('the code run in '+str(end_time-start_time))


# Set the parameters
mu_A = 1
mu_B = 1

# Set the maximum simulation time
t_max = 100

# Set the initial conditions
N_0 = 50  # Initial number of A-bacteria
M_0 = 50  # Initial number of B-bacteria

# Initialize time, A-bacteria population, and B-bacteria population arrays
time_points = [0]
populations_A =[]
populations_B = []


# Implement the Gillespie algorithm for five trajectories
for i in range(5):
    # Reset the population arrays for each trajectory
    population_A = [N_0]
    population_B = [M_0]
    time_points = [0]

    # Implement the Gillespie algorithm
    while time_points[-1] < t_max:
        N = population_A[-1]  # Current A-bacteria population
        M = population_B[-1]  # Current B-bacteria population

        # Calculate the birth and death rates for A-bacteria
        birth_rate_A = mu_A * N + epsilon
        death_rate_A = gamma * (M + N / 2) * N / Omega

        # Calculate the birth and death rates for B-bacteria
        birth_rate_B = mu_B * M + epsilon
        death_rate_B = gamma * (N + M / 2) * M / Omega

        # Calculate the total rates for A-bacteria and B-bacteria
        total_rate_A = birth_rate_A + death_rate_A
        total_rate_B = birth_rate_B + death_rate_B

        # Generate a random exponential time step
        delta_t = np.random.exponential(1 / (total_rate_A + total_rate_B))
        time_points.append(time_points[-1] + delta_t)
        # Determine the type of event (birth or death) for A-bacteria and B-bacteria
        p_A = birth_rate_A / total_rate_A
        p_B = birth_rate_B / total_rate_B
        if np.random.rand() < p_A:
            population_A.append(N + 1)
        else:
            population_A.append(N - 1)

        if np.random.rand() < p_B:
            population_B.append(M + 1)
        else:
            population_B.append(M - 1)
    populations_A.append(population_A)
    populations_B.append(population_B)
    

for i in range(5):
        
    # Plot the simulated trajectories of A-bacteria and B-bacteria
    plt.plot(range(len(populations_A[i])), populations_A[i], label='A-bacteria')
    plt.plot(range(len(populations_B[i])), populations_B[i], label='B-bacteria')

    # Set the plot labels and title
    plt.xlabel('Time')
    plt.ylabel('Population')
    plt.title('Simulated Trajectories of A-bacteria and B-bacteria')

    # Show the legend
    plt.legend()
    stri = str(i)+".png"

    plt.savefig(stri)

    # Show the plot
    plt.show()
plt.close()



#---------------Exercice 2------------------
c=1
k=1
D=0.1
eta=2
delta_t = 0.01
time_steps = 1000
iters = 400 

def Euler_step(X_last,c=1,k=1,D=0.1,eta=2):
    delta_x = -1* delta_t * (k * X_last - c) / eta + random.uniform(0, 1)
    return X_last + delta_x


X = np.zeros([time_steps])

for i in range(time_steps):
    X[i]=Euler_step(X[-1])
plt.plot(range(time_steps),X)
plt.title('Random trajectory')
plt.ylabel("Displacement")
plt.xlabel("timestep")
plt.show()
plt.show()
plt.close()

iters = 10000
X_ = np.zeros([iters,time_steps])

for j in tqdm(range(iters)):
    for i in range(time_steps):
        X_[j,i]=Euler_step(X_[j,-1])
X_avg = np.mean(X_, axis = 0)

plt.plot(range(time_steps),X_avg)
plt.title('Mean trajectory over 10000 trajectories')
plt.ylabel("Displacement")
plt.xlabel("timestep")
plt.show()
plt.close()

X_var = X_.var(0)
plt.plot(range(time_steps),X_var)
plt.title('Mean Variance over 10000 trajectories')
plt.ylabel("Variance")
plt.xlabel("timestep")
plt.show()
plt.close()


