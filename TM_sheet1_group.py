#Group
import numpy as np
import matplotlib.pyplot as plt

#Q3
def construct_hamiltonian(N: int, t: float, t_prime: float, mu_A=None, mu_B=None, potentials=False) -> np.ndarray:
    H  = np.zeros((N, N))
    for i in range(N):
        for j in range(N):
            if i + 1 == j and i % 1 == 0 and i > 0:
                H[i, j] = t_prime
            if i + 1 == j and i % 2 == 0:
                H[i, j] = t
            if i - 1 == j and i % 1 == 0 and i > 0:
                H[i, j] = t
            if i - 1 == j and i % 2 == 0:
                H[i, j] = t_prime

            if potentials:
                if i == j and i % 2 == 0:
                    H[i, j] = mu_A
                if i == j and i % 1 == 0 and i > 0:
                    H[i, j] = mu_B
    return H

def plot_eigenvals(H: np.ndarray, t: float, t_prime: float) -> tuple:
    eigenvals, eigenvecs = np.linalg.eig(H)
    plt.figure(figsize=(6, 4))
    plt.title(f'Eigenvalues for $t=${t}, $t^*=${t_prime}')
    plt.plot(sorted(eigenvals), 'o',color='red')
    plt.axhline(0, color='black', lw=1, ls='--')
    plt.xlabel('States')
    plt.ylabel('Energy')
    plt.show()
    return eigenvals, eigenvecs

N = 33

for t, t_prime in [(1, 0), (1, 1), (1, 2), (1, 3),(100,1), (1, 100)]:
    H = construct_hamiltonian(N, t, t_prime)
    eigenvals, eigenvecs = plot_eigenvals(H, t, t_prime)
    print(f'the approx 0 eigenvalue value is:{eigenvals[np.isclose(eigenvals, 0, 1e-10)]}')
    


#Q4
def plot_eigenvec(H:np.ndarray,t:float,t_prime:float)-> None:
    eigenvals, eigenvecs =  np.linalg.eig(H)
    zero_value_eigenvec = eigenvecs[:,np.isclose(eigenvals,0,1e-10)]
    
    #Finding the localized state
    abs_zve = abs(zero_value_eigenvec)
    maximum = max(abs_zve)
    index_maximum = np.argmax(abs_zve)
    
    plt.figure(figsize=(6,4))
    plt.title(f'Eigenvector for $t=${t}, $t^*=${t_prime}')
    x_s = range(1,len(zero_value_eigenvec)+1)
    plt.plot(x_s, np.squeeze(zero_value_eigenvec))#,label=f't={t},t_prime={t_prime}')
    if maximum>0.5:
        plt.scatter(index_maximum+1,zero_value_eigenvec[index_maximum], color = 'red',marker = 'x')
    plt.xlabel('Atom/site')
    plt.ylabel('Amplitude')
    plt.grid()
    #plt.legend()
    plt.show()

    
for t, t_prime in [(1, 0), (1, 1), (1, 2), (1, 3),(100,1), (1, 100)]:
    H = construct_hamiltonian(N,t,t_prime)
    plot_eigenvec(H,t,t_prime)
