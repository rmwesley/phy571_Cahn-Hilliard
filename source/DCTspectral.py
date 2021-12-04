import numpy as np
from scipy.fft import dctn, idctn, fftfreq

class Spectral:
    """Spectral stepping method for solving the Cahn-Hilliard equation"""
    """The class provides some built-in stepping methods to evolve the problem in time. Default is step_method5."""

    def __init__(self, Lx, Ly, Nx, Ny, eps, initialization='random', step=None):
        """Initialize the spatial and spectral grids based on passed arguments. Initializes the desired stepping function as well."""
        self.Lx, self.Ly = Lx, Ly
        self.Nx, self.Ny = Nx, Ny
        self.x, self.dx = np.linspace(0, Lx, Nx, retstep=True)
        self.y, self.dy = np.linspace(0, Ly, Ny, retstep=True)
        self.kx = 2 * np.pi * np.linspace(0, (Nx-1)/Lx, Nx).reshape(-1, 1)
        self.ky = 2 * np.pi * np.linspace(0, (Ny-1)/Ly, Ny).reshape(1, -1)

        self.k2 = self.kx**2 * np.ones([1, Ny]) + np.ones([Nx, 1]) * self.ky**2
        if initialization == 'trivial_x':
            solution = lambda x: np.tanh(x/(np.sqrt(2)*eps))
            self.u = ( solution(np.linspace(-.5,.5, Nx)).reshape(-1, 1) *
                np.ones([1, Ny]))
        elif initialization == 'trivial_y':
            solution = lambda y: np.tanh(y/(np.sqrt(2)*eps))
            self.u = ( solution(np.linspace(-.5,.5, Ny)).reshape(1, -1) *
                np.ones([Nx, 1]) )
        elif initialization == 'separated_x':
            self.u = np.concatenate( (np.ones([(Nx+1)//2, Ny]), 
                -np.ones([Nx//2, Ny])) )
        elif initialization == 'separated_y':
            self.u = np.concatenate( (np.ones([Nx, (Ny+1)//2]),
                -np.ones([Nx, Ny//2])), axis=1)
        else:
            self.u = np.random.uniform(-1, 1, [Nx, Ny])

        """The homogeneous free energy per molecule."""
        """This function isn't the total density of free energy! It can be though of as its first-order term, though."""
        self.func = lambda t: t**3 - t
        self.step = step 
        if step==None: self.step = self.step_method5
        self.eps = eps

    def advance(self, u_hat, dt, Nt):
        """Advance a spectral grid temporally. Does not update any os the object's attributes."""
        for i in range(Nt):
            u_hat = self.step(u_hat, dt)
        return u_hat

    def evolve(self, reps, dt, Nt):
        """Evolves the system temporally, updating self.u 'reps' times."""
        """To calculate the updated value, we evolve the wavelength spectral decomposition (Fourier transform) of self.u with the advance procedure. We then go back to the space domain through an inverse Fourier transform."""
        for i in range(reps):
            u_hat = dctn(self.u, norm='ortho')
            u_hat = self.advance(u_hat, dt, Nt)
            self.u = idctn(u_hat, norm='ortho')
    
    def step_method2(self, u_hat, dt):
        """Spectral version of semi-implicit Euler's scheme"""
        return ((u_hat -
            dt * self.k2 * dctn(self.func(self.u), norm='ortho')) /
            (1 + dt * self.eps**2 * self.k2**2))
    def step_method5(self, u_hat, dt):
        """Spectral version of linearly stabilized splitting scheme"""
        return ((u_hat -
            dt * self.k2 * (dctn(self.u**3, norm='ortho')) +
            3*dt * self.k2 * u_hat ) /
            (1 + dt * (self.eps**2 * self.k2**2 + 2*self.k2)))
