
**************************************
 The Chemical Concentration Evolution
**************************************


.. math:: \frac{d c_i^\gamma}{dt}=\sum_{j\neq i} w_{ij}^{\gamma},

.. math::
   F=&\frac{1}{2}\sum_{i,\gamma,j,\delta}c_i^\gamma c_j^\delta \omega^{\gamma \delta}_{ij}+\sum_{i,\gamma} c_i^\gamma E_i^\gamma \left(\sum_{j,\delta} c_j^{\delta}\psi_{ji}^{\delta\gamma}\right)-3k_BT\sum_{i,\gamma} c_i^\gamma \log\left(\frac{\alpha_i^\gamma}{\Lambda^\gamma}\sqrt{\pi e}\right)\\
   &+k_BT\sum_{i,\gamma} c_i^\gamma\log c_i^\gamma+k_BT\sum_{i}c_i^v\log (c_i^v),

.. math:: \Lambda^\gamma=\sqrt{\frac{h^2}{2\pi m^\gamma k_BT}}, \quad c_i^v=1-\sum_{\gamma}c_i^\gamma,

.. math:: \frac{\partial F}{\partial c_i^\gamma}=E_i^{\gamma}+\sum_{j,\delta}c_j^{\delta}(\omega_{ij}^{\gamma\delta}+\psi_{ij}^{\gamma\delta}{E'}_j^{\delta})-3k_BT\log\left(\frac{\alpha_i^\gamma}{\Lambda^\gamma}\sqrt{\pi e}\right)+k_BT\log\left( \frac{c_i^\gamma}{c_i^v}\right)


.. math:: \mu_i^{\gamma}=f_i^{\gamma}-3k_BT\log\left(\frac{\alpha_i^\gamma}{\Lambda^\gamma}\sqrt{\pi e}\right)+k_BT\log\left( \frac{c_i^\gamma}{c_i^v}\right)


.. math::
   Z_i^\gamma &= \int d^3 \mathbf{x} \exp \left(-\beta f_i^{\gamma} - \left(\frac{\mathbf{x}-\mathbf{x}_i}{\alpha_i^\gamma}\right)^2 \right)\\
   &= \left( \alpha_i^\gamma \sqrt{\pi} \right)^3 e^{ -\beta f_i^{\gamma}}

.. math:: Q_i^\gamma &= \exp \left(-\beta f_i^{\gamma} - \left(\frac{\mathbf{x}-\mathbf{x}_i}{\alpha_i^\gamma}\right)^2 \right)

.. math:: Q_j^\gamma &= \exp \left(-\beta f_j^{\gamma} - \left(\frac{\mathbf{x}-\mathbf{x}_j}{\alpha_j^\gamma}\right)^2 \right)


.. math:: \nabla u(0,0) = 0, \quad \nabla u(0,x_{ij}) = 0

.. math:: \nabla \nabla u(0,0) = \frac{2}{{\alpha_i^{\gamma}}^2} \left( \mathbf{e}_r \otimes\mathbf{e}_r + \mathbf{e}_\theta \otimes\mathbf{e}_\theta + \mathbf{e}_z \otimes\mathbf{e}_z\right)

.. math:: \nabla \nabla u(0,x_{ij}) = \frac{2}{{\alpha_j^{\gamma}}^2} \left( \mathbf{e}_r \otimes\mathbf{e}_r + \mathbf{e}_\theta \otimes\mathbf{e}_\theta + \mathbf{e}_z \otimes\mathbf{e}_z\right)


.. math::
   u(r,z) = U\left(\frac{z}{x_{ij}}\right) +\left(\frac{r}{\alpha_j^\gamma}\right)^2 \frac{z}{x_{ij}}
   +\left(\frac{r}{\alpha_i^\gamma}\right)^2 \left(1 -\frac{z}{x_{ij}}\right)

.. math:: U(0) = 0, \quad U(1) = \beta \left(f_j^{\gamma}-f_i^{\gamma}\right), \quad U'(0) = U'(1) = 0,

.. math:: U''(0) = 2 \left(\frac{x_{ij}}{\alpha_i^\gamma}\right)^2, \quad U''(1) &= 2 \left(\frac{x_{ij}}{\alpha_j^\gamma}\right)^2,


.. math::
   U(\zeta)=&\zeta^2 (1-\zeta)^2\biggl[\left(1-\zeta\right)\left(\frac{x_{ij}}{\alpha_i^\gamma}\right)^2+\zeta\left(\frac{x_{ij}}{\alpha_j^\gamma}\right)^2\biggr]\\
   &+\zeta^3\left(6 \zeta^2-15 \zeta+10\right)\beta \left(f_j^{\gamma}-f_i^{\gamma}\right)

.. math::
   \frac{U'(\zeta)}{\zeta(1-\zeta)}=&\frac{5}{2}\biggl[6\beta\left(f_j^\gamma-f_i^\gamma\right)+\left(\frac{x_{ij}}{\alpha_j}\right)^2
   -\left(\frac{x_{ij}}{\alpha_i}\right)^2\biggr]\zeta\left(1-\zeta\right)\\
   &+\left(\frac{x_{ij}}{\alpha_i}\right)^2\left(1-\zeta\right)-\left(\frac{x_{ij}}{\alpha_j}\right)^2\zeta

.. math:: U'(\zeta_0)=0

.. math::
   \frac{\partial^2 u}{\partial r^2}\bigg|_{z=z_0,r=0}=\frac{2}{{\alpha_j^\gamma}^2}\zeta_0+
   \frac{2}{{\alpha_i^\gamma}^2}\left(1-\zeta_0\right)


.. math::
   w_{ij}^{\gamma}-w_{ji}^{\gamma}=\frac{-x_{ij}^4}{\pi\tilde{\alpha}_{ij}^2(\zeta_0)}\sqrt{\frac{k_B T}{2 m_\gamma}}\left(\frac{c_i^{\gamma}c_j^v}{{\alpha_i^{\gamma}}^3}
   e^{-u(\zeta_0)}-\frac{c_j^{\gamma}c_i^v}{{\alpha_j^{\gamma}}^3}
   e^{-u(\zeta_0)+\beta(f_j^\gamma-f_i^\gamma)}\right)



















