import numpy as np

class Spectral:
    """Spectral stepping method for solving PDEs"""
    def __init__(self, Lx, Ly, Nx, Ny, eps, initialization='random'):

        """Initialize the spatial and spectral grids"""
        self.Lx, self.Ly = Lx, Ly
        self.Nx, self.Ny = Nx, Ny
        self.x, self.dx = np.linspace(0, Lx, Nx, retstep=True)
        self.y, self.dy = np.linspace(0, Ly, Ny, retstep=True)
        self.kx = 2 * np.pi * np.fft.fftfreq(Nx, d=self.dx).reshape(-1, 1)
        self.ky = 2 * np.pi * np.fft.fftfreq(Ny, d=self.dy).reshape(1, -1)
        
        self.k2 = self.kx**2 * np.ones([1, Ny]) + np.ones([Nx, 1]) * self.ky**2
        if initialization == 'trivial_x':
            solution = lambda x: np.tanh(x/(np.sqrt(2)*eps))
            self.u = solution(np.linspace(-.5,.5, Nx)).reshape(-1, 1) * np.ones([1, Ny])
        elif initialization == 'trivial_y':
            solution = lambda y: np.tanh(y/(np.sqrt(2)*eps))
            self.u = solution(np.linspace(-.5,.5, Ny)).reshape(1, -1) * np.ones([Nx, 1])
        elif initialization == 'separated_x':
            self.u = np.concatenate((np.ones([(Nx+1)//2, Ny]), -np.ones([Nx//2, Ny])))
        elif initialization == 'separated_y':
            self.u = np.concatenate((np.ones([Nx, (Ny+1)//2]), -np.ones([Nx, Ny//2])), axis=1)
        else:
            self.u = np.random.uniform(-1, 1, [Nx, Ny])
        self.f = lambda t: t**3 - t
        self.eps = eps

    def advance(self, dt, Nt):
        u_hat = np.fft.fft2(self.u)
        #print(u_hat.shape, self.k2.shape, self.f(self.u).shape, np.fft.fft2(self.f(self.u)).shape)
        for i in range(Nt):
            u_hat = ((u_hat - dt * self.k2 * np.fft.fft2(self.f(self.u)))/
                (1 + dt * self.eps**2 * self.k2**2))
            #u_hat = ((u_hat - dt * self.k2 * np.fft.fft2(self.f(self.u)))/
            #    (1 + 2*dt*self.k2  + dt * self.eps**2 * self.k2**2))
        self.u = np.fft.ifft2(u_hat).real

    def step(self, dt):
        """Make a step dt forward in time"""
        self.advance(dt, 1)
