import numpy as np

class Schrodinger:
    def __init__ (self, N, L): # N = number of points, L = size along x
        self.N = N
        self.L = L
        self.dx = L/(N-1)
        self.dt = 2*(self.dx)**2
        self.alpha = 2*(self.dx)**2/(self.dt) # hence alpha = 1
        # Arrays for tridiagonal matrix
        self.b = np.zeros(N, dtype=complex)
        self.e = np.zeros(N, dtype=complex)
        # Recursion coefficients (the other are a=1=c and d=-1=f)
        self.beta = np.zeros(N, dtype=complex)
        self.gamma = np.zeros(N, dtype=complex)
        self.x = np.zeros(N, dtype=complex)

    def psi0 (self, x0, k0, sigma0): # Wave function at t=0
        x = np.zeros(self.N)
        psi0 = np.zeros(self.N, dtype=complex)
        # Boundary condition
        psi0[0] = 0.0
        psi0[self.N-1] = 0.0
        x[self.N-1] = self.L # x[0] is already 0 thanks to np.zeros(N)
        for i in range(1, self.N-1): # [1, N-1[
            x[i] = i*self.dx
            psi0[i] = np.exp(1j*k0*x[i])*np.exp(-(x[i]-x0)**2/(2*sigma0**2))
        return x, psi0

    def initialize (self, V): # V is a list containing the value of the potential
        # Boundary condition (for e[0] and e[N-1] this is already set to 0.0 in the constructor)
        self.b[0] = 1.0
        self.b[self.N-1] = 1.0
        for k in range(1, self.N-1): # [1, N-1[
            self.b[k] = 1j*self.alpha - (self.dx**2)*V[k] - 2.0
            self.e[k] = 1j*self.alpha + (self.dx**2)*V[k] + 2.0
        # Set up recursion coefficients
        self.beta[0] = self.b[0]
        self.gamma[0] = 1.0/self.beta[0] # Because c[k] = 1.0 for all k
        for k in range(1, self.N): # [1, N[
            self.beta[k] = self.b[k] - 1.0*self.gamma[k-1] # Because a[k] = 1.0 for all k
            self.gamma[k] = 1.0/self.beta[k] # Because c[k] = 1.0 for all k

    def update (self, psi0, t_f):
        psi = psi0.copy()
        t = 0.0
        while t < t_f:
            self.x[0] = (self.e[0]*psi[0] - psi[1])/(self.beta[0]) # Because psi[k-1] = 0 for k=0
            for k in range(1, self.N-1): # [1, N-1[
                R_k = -psi[k-1] + self.e[k]*psi[k] - psi[k+1]
                self.x[k] = (R_k-1.0*self.x[k-1])/(self.beta[k]) # Because a[k] = 1.0 for all k 
            # Last element
            R_N_minus_1 = -psi[self.N-2] + self.e[self.N-1]*psi[self.N-1] # Because psi[k+1] = 0 for k=N-1
            self.x[self.N-1] = (R_N_minus_1 - 1.0*self.x[self.N-2])/(self.beta[self.N-1]) # Because a[k] = 1.0 for all k
            psi[self.N-1] = self.x[self.N-1]
            for k in reversed(range(self.N-1)): # ]N-1, 0] = N-2, N-3, ..., 2, 1, 0
                psi[k] = self.x[k] - self.gamma[k]*psi[k+1]
            # Increment time
            t += self.dt

        return psi
    
    def step_potential (self, V0): # V0 is the step height
        V = [0.0 if i <= int(self.N/2) else V0 for i in range(self.N)]
        return np.array(V)
    
    def barrier_potential (self, V0, width):
        n = int(width/(self.dx/2)) # Number of points in the half width of the barrier potential
        middle = int(self.N/2)
        V = np.zeros(self.N)
        V[middle-n:middle+n] = V0
        return V
    
    def well_potential (self, V0, width):
        n = int(width/self.dx/2) # Number of points in the half width of the barrier potential
        middle = int(self.N/2)
        V = np.zeros(self.N)
        V[:middle-n] = V0
        V[middle+n:] = V0
        return V
    
    def gaussian_potential (self, x0, sigma):
        V = np.zeros(self.N)
        A = 1.0/(np.sqrt(2*np.pi)*sigma) # Normalization constant
        for i in range(self.N):
            x = i*self.dx
            V[i] = A*np.exp(-(x-x0)**2/(2*sigma**2))/8.0 # Not optimized : the division by 8.0 works only for sigma = 0.05

        # Allow to reduce noise before and after the potential
        n = int(6*sigma/self.dx)
        middle = int(self.N/2)
        V[:middle-n] = 0.0
        V[middle+n:] = 0.0
        
        return V


