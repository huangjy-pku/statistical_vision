# Statistical Vision

This is a repo for the course *Computer Vision: Statistical Models for Marr's Paradigm* in PKU, which is also known as [Stat 232A](http://www.stat.ucla.edu/~sczhu/Courses/UCLA/Stat_232A/Stat_232A.html) in UCLA. This repo contains the implementation of projects in the course, and belows are summary notes of key ideas in this course.

#### Marr's paradigm
- Disentangle three components
  - Computation on top-level, problem formulation
  - Representation and algorithm
  - Implementation about hardware
- Representation framework
  - 2D image -> primal sketch -> 2.1D sketch -> 2.5D sketch -> 3D model representation
- Contribute to two fields
  - Computer Vision
  - Computational Neuroscience
#### Unsolved problem with current vision models
  - Continuous computing process
  - Top-down path, dark matter awareness
  - Flexibly jumping over the heterogeneous solution space
  - Definition of concepts (explicit grounding rather than simple classifiers)
  - Unification/interpretation of representations
#### Statistics of natural image
- Power law
  - Apply 2D Fourier Transform on image
  $$\hat{I}(\xi,\eta) = \sum\limits_{(x,y)\in D_1 \times D_2} I(x,y) e^{-i2\pi(\frac{x\xi}{W}+\frac{y\eta}{H})}$$
  - Inverse relation between amplitude and frequency
  $$\log{A(f)} = C - \log{f}, f=\sqrt{\xi^2+\eta^2}, A(f)=|\hat{I}(\xi,\eta)|^2$$
  - Invariance of power integrated on frequency band
  $$\iint\limits_{f_0<f<2f_0} A^2(f) d\xi d\eta = C, \forall f_0$$
- Kurtosis and sparsity
  - Histogram response
  $$h(a) = \int_\Omega f(I) \delta(<F,I>-a) dI$$
  - The point cloud formed by natural images contains many low-dimensional spikes, reflective of the high kurtosis of natural images, which may be used as guides, or viable directions, to traverse and better learn or represent the image space.
- Gradient histogram scale-invariance, as well as similar kurtosis
- Toy example for scale-invariance: line segment world, crops from various scale appear similarly
- The transition among low/middle/high entropy regimes
#### Texture
- Markov Random Fields (MRF) $\Longleftrightarrow$ Clique-based Gibbs Models
$$p(I_s  \vert I_{-s}) = p(I_s  \vert I_{N_s} \Longleftrightarrow p(I) = \frac{1}{Z}\exp\left(-\sum\limits_{C\in \mathcal{C}} \lambda_C\left(   I(C) \right) \right)$$
- Ising model: one of MRF models with pair potentials, $H$ denoting Hamiltonian Energy
$$H(I) = - \sum\limits_t \alpha I_t - \sum\limits_{<s,t>} \beta I_s I_t$$
- Potts model: generalization of Ising model on lattice with each particle having $n$ choices
$$H(I) = - \sum\limits_t \alpha I_t - \sum\limits_{<s,t>} \beta \delta(I_s,I_t)$$
- Gaussian Markov Random Fields (GMRF)
  - Grid structure: $\mathit{Neighbor}\left(I(x,y)\right) = \{I(x\pm1,y\pm1)\}$
  $$
  I \sim \mathcal{N}(0, \beta^{-1}Q^{-1}), Q_{(x_1,y_1),(x_2,y_2)} = \begin{cases}
  1, &(x_1,y_1)=(x_2,y_2) \\
  -\frac{1}{4}, &(x_2,y_2)\in \mathit{Neighbor\left(I(x_1,y_1)\right)} \\
  0, &\text{else} \end{cases}
  $$
  - Discrete form
  $$H(I) = -\beta \sum\limits_{x,y} \left[I(x+1,y)-I(x,y)\right]^2 + \left[I(x,y+1)-I(x,y)\right]^2$$
  - Continuous form
  $$H(I) = -\beta \iint \left[\nabla_x I\right]^2 + \left[\nabla_y I\right]^2 dxdy$$
  - Property 1: scale-invariance
    - Fourier Transform on first-order gradient
    $$F(I) = \hat{I}(\xi,\eta) = \iint I(x,y)e^{-2{\pi}i(\xi x+\eta y)} dxdy \\
    F(\nabla_x I) = 2{\pi}i{\xi}\hat{I}, \ \ \ (\nabla_y I) = 2{\pi}i\eta\hat{I}$$
    - Plancherel Theorem (power preservation)
    $$\beta \iint \left[\nabla_x I\right]^2 + \left[\nabla_y I\right]^2 dxdy = 4\pi^2\beta \iint (\xi^2+\eta^2)|\hat{I}(\xi,\eta)|^2 d{\xi}d\eta$$
    - Factorization: form in Gaussian
    $$p\left(\hat{I}(\xi,\eta)\right) \propto \exp\left(-4\pi^2\beta(\xi^2+\eta^2) |\hat{I}(\xi,\eta)|^2\right), \hat{I} \sim \mathcal{N}\left(0, \frac{1}{8\pi^2\beta(\xi^2+\eta^2)}\right)$$
  - Property 2: connection to heat equation
    - Specify dynamics
    $$\frac{dI(x,y,t)}{dt} = -\frac{\delta H\left(I(x,y,t)\right)}{\delta I}$$
    - Form in heat diffusion equation
    $$\frac{dI(x,y,t)}{dt} = \Delta I(x,y)$$
    - Variational method
    $$H\left(I(x,y,t)\right) = \iint L(I,I_x,I_y,x,y) dxdy \\
    \frac{\delta H}{\delta I} = \frac{\partial L}{\partial f} - \frac{d}{dx}(\frac{\partial L}{\partial I_x}) - \frac{d}{dy}(\frac{\partial L}{\partial I_y})$$
  - Limitation of GMRF: can capture the local regularity in natural images, but the joint pixel density is unimodal Gaussian; can be adapted as a filter for realistic image synthesis.
