******************************
Numerical Integration Details
******************************

The objectuve of this section is to outline a numerical approach to calculate the integrations of the form 

.. math::
   g(x,\alpha)=\left(\frac{\alpha}{\pi}\right)^{3/2}\int d^3\mathbf{r}e^{-\alpha(\mathbf{r}-\mathbf{x})^2} f(|\mathbf{r}|)

Rotating the coordinate system such that :math:`\mathbf{x}` is positioned on :math:`z` axis and using the spherical coordinate system,

.. math::
   g(x,\alpha)=\left(\frac{\alpha}{\pi}\right)^{3/2}\int_0^{\infty}\!\!\!\int_0^{\pi}\!\!\!\int_0^{2\pi}\!\!\!dr d\theta d\phi r^2\sin\theta e^{-\alpha (r^2+x^2-2rx\cos\theta)} f(r),

One can reduce the three dimensional integration to a one dimensional one

.. math::
   g(x,\alpha)&=2\pi\left(\frac{\alpha}{\pi}\right)^{3/2}\int_0^{\infty}\!\!\!\int_0^{\pi}\!\!\!dr d\theta r^2\sin\theta e^{-\alpha (r^2+x^2-2rx\cos\theta)} f(r),\\
   &=\frac{1}{x}\left(\frac{\alpha}{\pi}\right)^{1/2}\int_0^{\infty}\!\!\!dr rf(r) \biggl[e^{-\alpha (x-r)^2}-e^{-\alpha (x+r)^2} \biggr]

In our case the potential functions have cutoff :math:`r_c`. Therefore,

.. math::
   g(x,\alpha)=\frac{1}{x}\left(\frac{\alpha}{\pi}\right)^{1/2}\int_0^{r_c}\!\!\!dr rf(r) \biggl[e^{-\alpha (x-r)^2}-e^{-\alpha (x+r)^2} \biggr].

One can reduce the equation to

.. math::
   g(x,\alpha)=\frac{1}{x\sqrt{\pi}}\int_{(x-r_c)\sqrt{\alpha}}^{(x+r_c)\sqrt{\alpha}}\!\!\!\!\!\!\!\!\!dt e^{-t^2} \left(x-\frac{t}{\sqrt{\alpha}}\right)f\left(|x-\frac{t}{\sqrt{\alpha}}|\right),

Now we can employ Hermite-Gauss quadrature to calculate such integrations

.. math::
   g(x,\alpha)\approx\frac{1}{x\sqrt{\pi}}\sum_{i=1}^M \left(x-\frac{t_i}{\sqrt{\alpha}}\right)f\left(|x-\frac{t_i}{\sqrt{\alpha}}|\right)W_i,

where :math:`M`, and :math:`t_i` is the number of quadrature points and quadrature abscissas, respectively; and

.. math::
   W_i=\left\{\begin{array}{lr}
   \text{Hermite-Gauss Weight}\quad& (x-r_c)\sqrt{\alpha}<t_i<(x+r_c)\sqrt{\alpha}\\
   0 & \text{otherwise}
   \end{array}\right.

in addition to calculating :math:`g(x,\alpha)`, one needs to calculate the values of its derivative with respect to :math:`x` and :math:`\alpha`. Using the second fundamental theorem of calculus and the gaussian quadrature theory


.. math::
   \frac{\partial g(x,\alpha)}{\partial x}\approx &\frac{1}{x\sqrt{\pi}}\sum_{i=1}^m W_i\bigg[\frac{-1}{x}\left(x-\frac{t_i}{\sqrt{\alpha}}\right)f\left(|x-\frac{t_i}{\sqrt{\alpha}}|\right)\\
   &+f\left(|x-\frac{t_i}{\sqrt{\alpha}}|\right)+|x-\frac{t_i}{\sqrt{\alpha}}|f'\left(|x-\frac{t_i}{\sqrt{\alpha}}|\right)\bigg]


.. math::
   \frac{\partial g(x,\alpha)}{\partial \alpha}\approx&\frac{1}{2\alpha^{3/2}x\sqrt{\pi}}\sum_{i=1}^m W_it_i\bigg[f\left(|x-\frac{t_i}{\sqrt{\alpha}}|\right)\\
   &+|x-\frac{t_i}{\sqrt{\alpha}}|f'\left(|x-\frac{t_i}{\sqrt{\alpha}}|\right)\bigg]


.. math::
   \frac{\partial g(x,\alpha)}{\partial x}\approx &\frac{1}{x\sqrt{\pi}}\sum_{i=1}^m W_i\bigg[\frac{t_i}{x\sqrt{\alpha}}f\left(|x-\frac{t_i}{\sqrt{\alpha}}|\right)\\
   &+|x-\frac{t_i}{\sqrt{\alpha}}|f'\left(|x-\frac{t_i}{\sqrt{\alpha}}|\right)\bigg]


.. math::
   g(x,\alpha)&=\frac{\sqrt{\alpha}}{x} \frac{r_c}{\sqrt{\pi}} \int_{-1}^{1}dt tr_c f(|t r_c|)e^{-\alpha(x-tr_c)^2}

.. math::
   -\frac{\partial}{\partial x}g(x,\alpha)
   &=\frac{\sqrt{\alpha}}{x} \frac{r_c}{\sqrt{\pi}}\int_{-1}^{1}dt tr_cf(|t r_c|)e^{-\alpha(x-tr_c)^2}\bigg[\frac{1}{x}+2\alpha(x-tr_c)\bigg]

.. math::
   -\frac{\partial}{\partial \alpha}g(x,\alpha)
   &=\frac{\sqrt{\alpha}}{x} \frac{r_c}{\sqrt{\pi}}\int_{-1}^{1}dt tr_cf(|t r_c|)e^{-\alpha(x-tr_c)^2}\bigg[(x-tr_c)^2-\frac{1}{2\alpha}\bigg]



