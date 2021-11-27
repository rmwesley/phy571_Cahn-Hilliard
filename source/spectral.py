import numpy as np

class Spectral:
    """Spectral stepping method for solving the Cahn-Hilliard equation"""
    """The class provides some built-in stepping methods to evolve the system in time. Default is step_method5."""
    def step_method5(self, u_hat, dt):
        """Spectral version of linearly stabilized splitting scheme"""
        return ((u_hat +
            dt * self.k2 * (np.fft.fft2(self.u**3)) -
            3*dt * self.k2 * u_hat ) /
            (1 + dt * (self.eps**2 * self.k2**2 - 2*self.k2)))

    def __init__(self, Lx, Ly, Nx, Ny, eps, initialization='random', step=step_method5):
        """Initialize the spatial and spectral grids based on passed arguments. Initializes the desired stepping function as well."""
        self.Lx, self.Ly = Lx, Ly
        self.Nx, self.Ny = Nx, Ny
        self.x, self.dx = np.linspace(0, Lx, Nx, retstep=True)
        self.y, self.dy = np.linspace(0, Ly, Ny, retstep=True)
        self.kx = 2 * np.pi * np.fft.fftfreq(Nx, d=self.dx).reshape(-1, 1)
        self.ky = 2 * np.pi * np.fft.fftfreq(Ny, d=self.dy).reshape(1, -1)

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

        """The homogeneous free energy per molecule function."""
        self.func = lambda t: t**3 - t
        """This function isn't the total density of free energy! It can be though of as only its first-order term. The second-order term is dependent on the squared 2-norm of the gradient of the concentration, that is, this one accounts for the increment in the free energy due to the variation of the concentration or heterogeneity in the system (Higher variation, higher concentration gradient). For a high epsilon value, to minimize the Hemholtz free energy we must make the second order term smaller overall (tendence towards zero). This gives rise to a weaker phase separation and thus results in a more homogeneous, miscible, monophasic fluid"""
    #TODO
    #Call step methods "recursion", "recursive_step", "step", "evolution" ...?
        self.step = step 
        self.eps = eps
    
    def advance(self, dt, Nt):
        u_hat = np.fft.fft2(self.u)
        for i in range(Nt):
            u_hat = self.step(u_hat, dt)
        self.u = np.fft.ifft2(u_hat).real
    
    def step_method2(self, u_hat, dt):
        """Spectral version of semi-implicit Euler's scheme"""
        return ((u_hat -
            dt * self.k2 * np.fft.fft2(self.func(self.u))) /
            (1 + dt * self.eps**2 * self.k2**2))
        #u_hat = ((u_hat - dt * self.k2 * np.fft.fft2(func(self.u))) /
        #    (1 + 2*dt*self.k2  + dt * self.eps**2 * self.k2**2))
