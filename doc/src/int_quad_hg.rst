******************************
Numerical Integration Details
******************************

The objective of this section is to outline a numerical approach to calculate the integrations of the form

.. math:: g(\mathbf{x},\alpha)=\frac{1}{\left(\alpha\sqrt{\pi}\right)^{3}}\int d^3\mathbf{r}e^{-\left(\mathbf{y}-\mathbf{x}\right)^2/\alpha^2} f(|\mathbf{y}|),

rotating the coordinate system such that :math:`\mathbf{x}` is positioned on :math:`z` axis and using the spherical coordinate system,

.. math:: g(|\mathbf{x}|,\alpha)=\frac{1}{\left(\alpha\sqrt{\pi}\right)^{3}}\int_0^{\infty}\!\!\!\int_0^{\pi}\!\!\!\int_0^{2\pi}\!\!\!dr d\phi d\theta r^2\sin\phi e^{-\left(r^2+x^2-2rx\cos\phi\right)/\alpha^2} f(r),

One can reduce the three dimensional integration to a one dimensional one

.. math:: g(x,\alpha)=\frac{1}{x\alpha\sqrt{\pi}}\int_0^{\infty}\!\!\!dr rf(r) \biggl[e^{-(x-r)^2/\alpha^2}-e^{-(x+r)^2/\alpha^2} \biggr],

In our case the potential functions normal have a cutoff radius namely, :math:`r_c`. Therefore, we can limit our integration domain from zero to :math:`r_c`

.. math:: g(x,\alpha)=\frac{1}{x\alpha\sqrt{\pi}}\int_0^{r_c}\!\!\!dr rf(r) \biggl[e^{-(x-r)^2/\alpha^2}-e^{-(x+r)^2/\alpha^2} \biggr],

Employing two change of variables the integration simplifies to

.. math:: g(x,\alpha)=\frac{1}{x\sqrt{\pi}}\int_{\left(x-r_c\right)/\alpha}^{\left(x+r_c\right)/\alpha}\!\!\!dt e^{-t^2} \left(x-\alpha t\right)f(|x-\alpha t |),

where :math:`|.|` denotes absolute value. Before going further it is instructive to show the derivative of :math:`g(x,\alpha)` with respect to :math:`x` and :math:`\alpha`, in terms of this simplified integration

.. math::
   \frac{\partial}{\partial x}g(x,\alpha)=&\frac{1}{x\sqrt{\pi}}\int_{\left(x-r_c\right)/\alpha}^{\left(x+r_c\right)/\alpha}\!\!\!dt e^{-t^2} \biggl[f'(|x-\alpha t |)+|x-\alpha t| f(|x-\alpha t |)\biggr]\\
   &-\frac{g(x,\alpha)}{x}

.. math:: \frac{\partial}{\partial \alpha}g(x,\alpha)=-\frac{1}{x\sqrt{\pi}}\int_{\left(x-r_c\right)/\alpha}^{\left(x+r_c\right)/\alpha}\!\!\!dt t e^{-t^2} \biggl[f'(|x-\alpha t |)+|x-\alpha t| f(|x-\alpha t |)\biggr],

Now we are ready employ Hermite-Gauss quadrature to calculate such integrations

.. math:: g(x,\alpha)\approx\frac{1}{x\sqrt{\pi}}\sum_{i=1}^m \left(x-\alpha x_i\right)f\left(|x-\alpha x_i|\right)w_i,

where :math:`m` is the number of quadrature points and :math:`x_i` and :math:`w_i` are the quadrature abscissas and weight at :math:`i` th point, respectively. Of course if :math:`|x-\alpha x_i|\ge r_c` that term would be ignored. Considering the fact that Hermite-Gauss quadrature abscissas range from a negative value to absolute value of the same value i.e.

.. math:: -x_{\mathrm{max}}\le x_i\le x_{\mathrm{max}}, \quad i=1, \cdots, m

it can be realized that the result of this numerical scheme is :math:`0` whenever

.. math:: r_c\le x-\alpha x_{\mathrm{max}}

thus giving us a closed form for the effictive cuttoff radius to be used in the neighbor lists. In other words two atoms (sites) would be interacting if and only if

.. math:: x_{ij}-\alpha_{ij} x_{\mathrm{max}} \le r_c



