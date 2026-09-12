"""A two-dimensional Kelvin wave on affine Couette flow, with exact viscosity."""
from dataclasses import dataclass
import numpy as np


@dataclass(frozen=True)
class KelvinWave:
    shear: float = 1.0
    viscosity: float = 0.005
    kx: float = 1.0
    ky0: float = 4.0
    initial_speed: float = 1.0

    def __post_init__(self):
        for key in ("shear", "viscosity", "kx", "ky0", "initial_speed"):
            if np.ndim(getattr(self,key)) or not np.isfinite(getattr(self,key)):
                raise ValueError("parameters must be finite scalars")
        if self.viscosity < 0 or self.kx == 0 or self.initial_speed <= 0:
            raise ValueError("viscosity >=0, nonzero kx, and positive initial speed required")

    def state(self, time):
        t=np.asarray(time,float)
        if np.any(~np.isfinite(t)) or np.any(t<0):
            raise ValueError("time must be finite and nonnegative")
        ky=self.ky0-self.shear*self.kx*t
        k2=self.kx**2+ky**2
        # Equivalent to the expanded cubic, but explicitly positive: integrate
        # the square of the affine ky using midpoint plus variance on [0,t].
        integral=t*(self.kx**2+(self.ky0-self.shear*self.kx*t/2)**2
                    +(self.shear*self.kx*t)**2/12)
        zeta0=self.initial_speed*np.hypot(self.kx,self.ky0)
        zeta=zeta0*np.exp(-self.viscosity*integral)
        amplitude=np.stack((-ky*zeta/k2,self.kx*zeta/k2),axis=-1)
        speed=zeta/np.sqrt(k2)
        qxy=amplitude[...,0]*amplitude[...,1]/2
        energy=speed**2/4
        production=-self.shear*qxy
        dissipation=self.viscosity*zeta**2/2
        pressure_amplitude=2*self.shear*self.kx*amplitude[...,1]/k2
        return dict(ky=ky,k2=k2,damping_integral=integral,vorticity_amplitude=zeta,
            velocity_amplitude=amplitude,speed=speed,energy=energy,qxy=qxy,
            production=production,dissipation=dissipation,energy_rate=production-dissipation,
            logarithmic_speed_rate=self.shear*self.kx*ky/k2-self.viscosity*k2,
            wavelength=2*np.pi/np.sqrt(k2),pressure_amplitude=pressure_amplitude,
            heat_reference_speed=self.initial_speed*np.exp(-self.viscosity*(self.kx**2+self.ky0**2)*t))

    def field(self,x,y,time):
        state=self.state(time)
        phase=self.kx*np.asarray(x)+state["ky"]*np.asarray(y)
        return state["velocity_amplitude"]*np.sin(phase)[...,None]

    def pressure(self,x,y,time):
        state=self.state(time)
        return state["pressure_amplitude"]*np.cos(self.kx*np.asarray(x)+state["ky"]*np.asarray(y))

    def amplitude_rhs(self,time,amplitude):
        """Projected momentum ODE; independent of the vorticity solution formula."""
        ky=self.ky0-self.shear*self.kx*time
        wavevector=np.array([self.kx,ky])
        k2=float(wavevector@wavevector)
        return (-self.shear*amplitude[1]*np.array([1.,0.])
                +2*self.shear*self.kx*amplitude[1]*wavevector/k2
                -self.viscosity*k2*amplitude)


def integrate_amplitude(model,end,step):
    """Fixed-step classical RK4, without projecting away constraint errors."""
    steps=int(round(end/step))
    if steps<=0 or not np.isclose(steps*step,end):
        raise ValueError("end must be an integer multiple of positive step")
    times=np.linspace(0,end,steps+1)
    values=np.empty((steps+1,2))
    values[0]=model.initial_speed*np.array([-model.ky0,model.kx])/np.hypot(model.kx,model.ky0)
    for i,t in enumerate(times[:-1]):
        b=values[i]
        k1=model.amplitude_rhs(t,b)
        k2=model.amplitude_rhs(t+step/2,b+step*k1/2)
        k3=model.amplitude_rhs(t+step/2,b+step*k2/2)
        k4=model.amplitude_rhs(t+step,b+step*k3)
        values[i+1]=b+step*(k1+2*k2+2*k3+k4)/6
    return times,values
