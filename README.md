# sphere-vector-fields
SPHERE VECTOR FIELDS
# Sphere Vector Fields Engine 🌌

An optimized, publication-ready Python visualization and numerical simulation toolkit designed to solve and analyze **real-analytic vector fields globally constrained to the 2-sphere (\(S^2\))**. 

This engine enables researchers, students, and mathematical hobbyists to instantly model continuous dynamical systems on compact surfaces without numerical drift, mapping trajectories exactly along the manifold curvature.

---

## 🚀 Key Features

* **Strict Tangency Enforcement:** All built-in equations are strictly bound by the geometric constraint \(xP(x,y,z) + yQ(x,y,z) + zR(x,y,z) = 0\).
* **Symplectic Stability:** The Runge-Kutta 4(5) solver features dynamic algebraic normalization at each integration step to prevent numerical leakage off the sphere surface.
* **Vector Field Topology Analysis:** Includes a numerical divergence estimator (\(dP/dx + dQ/dy + dR/dz\)) to quickly verify topological singularities (saddles, sources, sinks, and centers).
* **3D Matplotlib Renderings:** Generates wireframe models wrapped in vector grids with multi-trajectory tracking overlay capabilities.

---

## 📐 Mathematical Presets & Topology

The toolkit comes pre-programmed with three classic analytical presets illustrating the **Poincaré-Hopf Theorem** (\(\sum \text{Index}(x_i) = \chi(S^2) = 2\)):

### 1. Pure Rotational Field (`"rotation"`)
Models a rigid body spinning around the polar axis. Lines of flow exactly match lines of latitude.
* **Vector System:** 
  \[\mathbf{V}(x, y, z) = \begin{pmatrix} -y \\ x \\ 0 \end{pmatrix}\]
* **Topology:** 2 Centers located at the poles \((0,0,1)\) and \((0,0,-1)\), each carrying a topological index of \(+1\).

### 2. Gradient Flow Field (`"gradient"`)
Models an irreversible flow running down a scalar height landscape from an unstable point to a stable point.
* **Vector System:** 
  \[\mathbf{V}(x, y, z) = \begin{pmatrix} -xz \\ -yz \\ 1 - z^2 \end{pmatrix}\]
* **Topology:** 1 Source at the South Pole (Index \(+1\)) and 1 Sink at the North Pole (Index \(+1\)).

### 3. Quadratic Dipole & Saddle Config (`"dipole_saddle"`)
A more complex, higher-order polynomial vector field showcasing topological deformation.
* **Vector System:** 
  \[\mathbf{V}(x, y, z) = \begin{pmatrix} z - x^2z \\ -xyz \\ -x(1-z^2) \end{pmatrix}\]

---

## 🛠️ Installation & Requirements

Ensure you have a modern Python installation (Python 3.8+) along with standard scientific computing packages:

```bash
pip install numpy scipy matplotlib
```

---

## 💻 Quickstart Usage Example

To generate an interactive 3D model tracking custom particle trajectories, simply run the engine as shown below:

```python
from solver import SphereVectorField
import matplotlib.pyplot as plt

# 1. Initialize the solver with your chosen flow topology
engine = SphereVectorField(field_type="gradient")

# 2. Integrate particle trajectories by supplying a 3D starting vector
# The system automatically normalizes the input to sit exactly on S^2
trajectories = [
    engine.integrate_trajectory([0.9, 0.1, -0.4], t_max=8.0),
    engine.integrate_trajectory([-0.5, 0.8, -0.2], t_max=8.0),
    engine.integrate_trajectory([0.1, -0.9, -0.6], t_max=8.0)
]

# 3. Compile the global 3D visualization window
fig = engine.plot_3d(num_arrows=15, trajectories=trajectories)
plt.show()
```

---

## 📝 Ongoing Academic Exploration

This open-source engine serves as a foundation for exploring the **Qualitative Theory of Planar Dynamical Systems** and variants of **Hilbert's 16th Problem**. Current frontiers include calculating the structural boundaries of algebraic limit cycles and determining the local cyclicity bounds of analytic families under compact parameter deformations.

Contributions, feature extensions, and numerical optimization pull requests are welcome!
