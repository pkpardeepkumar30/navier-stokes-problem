"""Exact periodic covariance toys; not the paper's viscous pulse construction."""
import numpy as np


def positive_integer(value, name):
    if not isinstance(value, (int, np.integer)) or isinstance(value, bool) or value <= 0:
        raise ValueError(f"{name} must be a positive integer for a periodic carrier")
    return value


def periodic_grid(size):
    positive_integer(size, "grid size")
    points = 2*np.pi*np.arange(size)/size
    return points, *np.meshgrid(points, points, indexing="ij")


def derivative(values, axis, length=2*np.pi):
    """Real Fourier derivative; even-grid Nyquist mode has derivative zero.

    This convention cannot recover an under-resolved carrier. Validate bandwidth
    separately; small covariance error alone does not demonstrate resolution.
    """
    values = np.asarray(values, float)
    size = values.shape[axis]
    modes = 2*np.pi*np.fft.fftfreq(size, d=length/size)
    if size % 2 == 0:
        modes[size//2] = 0
    shape = [1]*values.ndim
    shape[axis] = size
    return np.fft.ifft(1j*modes.reshape(shape)*np.fft.fft(values, axis=axis), axis=axis).real


def two_mode_velocity(x, y, lambda_plus=3.0, lambda_minus=1.0, n=8, m=9):
    positive_integer(n, "n")
    positive_integer(m, "m")
    for value in (lambda_plus, lambda_minus):
        if np.ndim(value) or not np.isfinite(value) or value < 0:
            raise ValueError("covariance eigenvalues must be finite and nonnegative")
    first = np.sqrt(lambda_plus)*np.sin(n*(x-y))
    second = np.sqrt(lambda_minus)*np.sin(m*(x+y))
    return np.stack((first+second, first-second), axis=-1)


def two_mode_target(lambda_plus=3.0, lambda_minus=1.0):
    return 0.5*np.array([[lambda_plus+lambda_minus, lambda_plus-lambda_minus],
                         [lambda_plus-lambda_minus, lambda_plus+lambda_minus]])


def covariance(velocity, axes):
    products = velocity[..., :, None]*velocity[..., None, :]
    return products.mean(axis=axes)


def divergence(velocity):
    return derivative(velocity[..., 0], 0)+derivative(velocity[..., 1], 1)


def envelope(y, epsilon=0.4):
    if np.ndim(epsilon) or not np.isfinite(epsilon) or not 0 <= epsilon < 1:
        raise ValueError("epsilon must be in [0,1) for this positive envelope")
    return 1+epsilon*np.cos(y), -epsilon*np.sin(y), -epsilon*np.cos(y)


def enveloped_velocity(x, y, epsilon=0.4, k=8, slope=1):
    """w=(psi_y,-psi_x), psi=a(y)*cos(k*(x+slope*y))/k.

    The a'/k term is essential for exact incompressibility.
    """
    positive_integer(k, "k")
    if not isinstance(slope, (int, np.integer)):
        raise ValueError("slope must be an integer for periodicity in y")
    a, ap, _ = envelope(y, epsilon)
    phase = k*(x+slope*y)
    return np.stack((ap/k*np.cos(phase)-slope*a*np.sin(phase), a*np.sin(phase)), axis=-1)


def envelope_covariance(y, epsilon=0.4, k=8, slope=1):
    a, ap, _ = envelope(y, epsilon)
    direction = np.array([[slope*slope, -slope], [-slope, 1]])
    leading = a[..., None, None]**2/2*direction
    exact = leading.copy()
    exact[..., 0, 0] += ap*ap/(2*k*k)
    return exact, leading


def envelope_mean_force(y, epsilon=0.4, slope=1):
    """Minus divergence of Q: horizontal shear force plus vertical gradient."""
    a, ap, _ = envelope(y, epsilon)
    return np.stack((slope*a*ap, -a*ap), axis=-1)


def averaged_transport(velocity):
    """Mean over x of (w dot grad)w, computed before averaging."""
    value = (velocity[..., 0, None]*derivative(velocity, 0)
             + velocity[..., 1, None]*derivative(velocity, 1))
    return value.mean(axis=0)


def force_from_mean_covariance(mean_covariance):
    """For Q(y), -div Q_i = -partial_y Q_iy."""
    return -derivative(mean_covariance[..., :, 1], 0)


def project_mean_force(force):
    """Periodic Leray projection for a vector F(y); preserve its zero mode.

    All nonconstant vertical components are gradients and removed. The horizontal
    component is already divergence free. This is an exact restricted projection,
    not a general two-dimensional Poisson solver.
    """
    projected = np.array(force, copy=True)
    projected[:, 1] = force[:, 1].mean()
    return projected
