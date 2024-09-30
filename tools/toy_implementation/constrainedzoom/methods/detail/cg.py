import numpy as np

def conjugate_gradient(Q, b, rtol=1e-6):
    """Solve Qx=b using conjugate gradient. Return x.

    Q should be a function such that for any vector a it returns Qa.

    If mask is not None, it is a function that masks (zeros) pixels which are in the null space of Q.
    b is assumed to lie outside the null space of Q."""


    Ndim = len(b)
    x = np.zeros(Ndim)

    residual = - b  # + Q(x), but that's zero
    direction = -residual

    scale = b.norm()

    for i in range(Ndim + 1):
        # find distance to travel in direction d
        alpha = -np.dot(residual, direction) / np.dot(direction, Q(direction))
        # update x
        x = x + alpha * direction

        # displacement vector from correct solution
        residual = Q(x) - b

        if (residual.norm() < rtol * scale):
            break

        # now figure out the direction of the next update
        # it must be Q-orthogonal to all previous updates
        beta = np.dot(residual, Q(direction)) / np.dot(direction, Q(direction))
        direction = -residual + beta * direction

    print("Converged after %d iterations" % i)

    return x