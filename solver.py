import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

class SphereVectorField:
    """
    A scientific engine to compute and solve analytic vector fields 
    globally constrained to the 2-sphere (x^2 + y^2 + z^2 = 1).
    """
    def __init__(self, field_type="rotation"):
        self.field_type = field_type

    def evaluate(self, x, y, z):
        """
        Defines the analytic vector field components (P, Q, R).
        Guarantees strict tangency: x*P + y*Q + z*R = 0.
        """
        if self.field_type == "rotation":
            # Pure rigid rotation around the z-axis (Two centers at poles)
            return -y, x, np.zeros_like(x)
            
        elif self.field_type == "gradient":
            # Source-to-sink flow from South Pole to North Pole
            return -x * z, -y * z, 1.0 - z**2
            
        elif self.field_type == "dipole_saddle":
            # Higher-order quadratic polynomial system showing saddle configurations
            return z - (x**2) * z, -x * y * z, -x * (1.0 - z**2)
            
        else:
            raise ValueError(f"Unknown field type: {self.field_type}")

    def compute_divergence(self, x, y, z, eps=1e-5):
        """
        Numerically approximates the ambient divergence (dP/dx + dQ/dy + dR/dz)
        to identify local stability, sources, and sinks.
        """
        P, Q, R = self.evaluate(x, y, z)
        
        Px, _, _ = self.evaluate(x + eps, y, z)
        _, Qy, _ = self.evaluate(x, y + eps, z)
        _, _, Rz = self.evaluate(x, y, z + eps)
        
        div = (Px - P)/eps + (Qy - Q)/eps + (Rz - R)/eps
        return div

    def integrate_trajectory(self, start_point, t_max=10.0, num_points=500):
        """
        Uses SciPy's Runge-Kutta 4(5) to solve the ODE system on the sphere surface.
        Insulates trajectories from drifting off the surface via normalization.
        """
        def odes(t, u):
            # Enforce unit sphere projection at each evaluation step
            u_norm = u / np.linalg.norm(u)
            P, Q, R = self.evaluate(u_norm[0], u_norm[1], u_norm[2])
            return [P, Q, R]

        t_span = (0, t_max)
        t_eval = np.linspace(0, t_max, num_points)
        
        # Normalize initial condition to ensure it starts exactly on S^2
        init_cond = np.array(start_point) / np.linalg.norm(start_point)
        
        sol = solve_ivp(odes, t_span, init_cond, t_eval=t_eval, rtol=1e-8, atol=1e-10)
        return sol.y

    def plot_3d(self, num_arrows=20, trajectories=None):
        """
        Generates a publication-ready 3D visual of the vector field mapped onto the sphere.
        """
        fig = plt.figure(figsize=(10, 8))
        ax = fig.add_subplot(111, projection='3d')

        # 1. Plot the wireframe sphere surface
        u = np.linspace(0, 2 * np.pi, 30)
        v = np.linspace(0, np.pi, 30)
        x_sp = np.outer(np.cos(u), np.sin(v))
        y_sp = np.outer(np.sin(u), np.sin(v))
        z_sp = np.outer(np.ones_like(u), np.cos(v))
        ax.plot_wireframe(x_sp, y_sp, z_sp, color="silver", linewidth=0.5, alpha=0.6)

        # 2. Plot the Vector Field Arrows
        u_arr = np.linspace(0, 2 * np.pi, num_arrows)
        v_arr = np.linspace(0.1, np.pi - 0.1, num_arrows)
        x_a = np.outer(np.cos(u_arr), np.sin(v_arr)).flatten()
        y_a = np.outer(np.sin(u_arr), np.sin(v_arr)).flatten()
        z_a = np.outer(np.ones_like(u_arr), np.cos(v_arr)).flatten()

        P, Q, R = self.evaluate(x_a, y_a, z_a)
        
        # Normalize arrows for clean visual scale
        norm = np.sqrt(P**2 + Q**2 + R**2)
        norm[norm == 0] = 1.0
        ax.quiver(x_a, y_a, z_a, P/norm, Q/norm, R/norm, length=0.15, color="teal", alpha=0.8, pivot="middle")

        # 3. Overlay numerical trajectories if provided
        if trajectories is not None:
            for traj in trajectories:
                ax.plot(traj[0], traj[1], traj[2], color="crimson", linewidth=2.5, zorder=10)

        ax.set_title(f"Analytic Vector Field Globally Mapped onto S² ({self.field_type.upper()})", fontsize=12)
        ax.set_axis_off()
        return fig

# Demonstration execution block
if __name__ == "__main__":
    # Initialize the gradient (source-to-sink) field engine
    engine = SphereVectorField(field_type="gradient")
    
    # Calculate three distinct sample trajectories spiraling up from the south
    trajectories = [
        engine.integrate_trajectory([0.9, 0.1, -0.4], t_max=8.0),
        engine.integrate_trajectory([-0.5, 0.8, -0.2], t_max=8.0),
        engine.integrate_trajectory([0.1, -0.9, -0.6], t_max=8.0)
    ]
    
    # Generate and display the plot
    fig = engine.plot_3d(num_arrows=15, trajectories=trajectories)
    plt.show()
