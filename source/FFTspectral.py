import numpy as np

class FFTspectral:
    """Class representing the Cahn-Hilliard problem system.
    Evolved in time through a FFT scheme.
    The class provides some built-in stepping methods for time evolution.
    Default is step_method5.
    """

    def __init__(self, Lx, Ly, Nx, Ny, eps, initialization='random'):
        """Initialize the spatial and spectral grids with the given arguments."""
        self.Lx, self.Ly = Lx, Ly
        self.Nx, self.Ny = Nx, Ny
        self.x, self.dx = np.linspace(0, Lx, Nx, retstep=True, endpoint=False)
        self.y, self.dy = np.linspace(0, Ly, Ny, retstep=True, endpoint=False)
        self.kx = 2 * np.pi * np.fft.fftfreq(Nx, d=self.dx).reshape(-1, 1)
        self.ky = 2 * np.pi * np.fft.fftfreq(Ny, d=self.dy).reshape(1, -1)
        self.t = 0

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
        self.u_hat = np.fft.fft2(self.u)

        """The homogeneous free energy per molecule.
        This function isn't the total density of free energy!
        It can be though of as its first-order term, though.
        """
        self.func = lambda t: t**3 - t
        self.step = self.step_method5
        self.eps = eps

    def advance(self, dt, Nt):
        """Evolves a spectral grid temporally without updating it."""
        for i in range(Nt):
            self.u_hat = self.step(dt)
        return self.u_hat

    def evolve(self, reps, dt, Nt):
        """Evolves the system temporally, updating self.u 'reps' times.
        To calculate the updated values, we work in the wavelength domain.
        That is, we evolve the spectral decomposition (FT) of self.u, self.u_hat, temporally .
        We then go back to the space domain through an inverse FT to update self.u.
        """
        for i in range(reps):
            self.u_hat = self.advance(dt, Nt)
            self.u = np.fft.ifft2(self.u_hat).real
            self.t += Nt*dt
    
    def step_method2(self, dt):
        """Spectral version of semi-implicit Euler's scheme"""
        return ((self.u_hat -
            dt * self.k2 * np.fft.fft2(self.func(self.u))) /
            (1 + dt * self.eps**2 * self.k2**2))
    def step_method5(self, dt):
        """Spectral version of linearly stabilized splitting scheme"""
        return ((self.u_hat -
            dt * self.k2 * (np.fft.fft2(self.u**3)) +
            3*dt * self.k2 * self.u_hat ) /
            (1 + dt * (self.eps**2 * self.k2**2 + 2*self.k2)))
