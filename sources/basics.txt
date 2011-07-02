
Relation between Source Excitation and Far-Field Pattern
========================================================

:Release: |version|
:Date: |today|

Many students either dislike or afraid of Electromagnetic theory! Obviously, the
reason is the mathematical complexity involved in the EM theory. At the same time,
nobody has that much trouble with subjects such as signal processing, communication
theory, etc. But, if observed closely, the fundamental theory behind all these 
subjects is nothing but the `Fourier analysis <http://en.wikipedia.org/wiki/Fourier_analysis>`_.
In fact, Fourier analysis plays crucial role in understanding many areas of science and
engineering.

So, it is the authors' opinion that if the books on antenna theory starts with a
slightly different approach using **Fourier analysis as the starting point**, then
the learning process becomes easier and more fun. In this brief introductory 
article, authors use this approach to explain the antenna theory  (and of course, 
array antenna theory too).

Symmetric Maxwell's Equations
-----------------------------

Symmetric Maxwell's equations are given as

.. math::

        \nabla \cdot \mathbf{D} &= \rho_e
        \\
        \nabla \cdot \mathbf{B} &= \rho_m
        \\
        \nabla \times \mathbf{E} &= -\mathbf{J}_m-\frac{\partial \mathbf{B}}{\partial t}
        \\
        \nabla \times \mathbf{H} &= +\mathbf{J}_e+\frac{\partial \mathbf{D}}{\partial t}
        
where :math:`\mathbf{D}=\epsilon\mathbf{E}` and :math:`\mathbf{B}=\mu\mathbf{H}`.

Vector Potentials and Helmholtz Wave Equations
----------------------------------------------

Assuming the entire system is `linear`, any source can be separated into 
electric and magnetic currents (sources). The vector potentials :math:`\mathbf{A}`
and :math:`\mathbf{F}` corresponding to these given electric and magnetic currents
:math:`\mathbf{J}_e` and :math:`\mathbf{J}_m` are given as

.. math::

        \nabla^{2}\mathbf{A}+k^{2}\mathbf{A} &= -\mu\mathbf{J}_{e}
        \\
        \nabla^{2}\mathbf{F}+k^{2}\mathbf{F} &= -\epsilon\mathbf{J}_{m}.
        
The above equations are derived from the Maxwell's equations using simple
`vector identities <http://en.wikipedia.org/wiki/Vector_identities#Summary_of_all_identities>`_. 
Using the `Sommerfeld radiation condition <http://en.wikipedia.org/wiki/Sommerfeld_radiation_condition>`_,
solutions to the above `inhomogeneous Hemlholtz equations <http://en.wikipedia.org/wiki/Helmholtz
_equation#Inhomogeneous_Helmholtz_equation>`_ are given as

.. math::

        \mathbf{A} &= \frac{\mu}{4\pi}\int_\mathrm{volume}\mathbf{J}_e\frac{e^{-jk_0|\mathbf{r-r^\prime}|}}{|\mathbf{r-r^\prime}|}dv^\prime
        \\
        \mathbf{F} &= \frac{\epsilon}{4\pi}\int_\mathrm{volume}\mathbf{J}_m\frac{e^{-jk_0|\mathbf{r-r^\prime}|}}{|\mathbf{r-r^\prime}|}dv^\prime.
        
And of course, from the definitions, electric and magnetic fields can be written in terms of
the vector potentials as

.. math::

        \mathbf{H}_A &= \frac{1}{\mu}\left(\nabla\times\mathbf{A}\right)
        \\
        \mathbf{E}_A &= \frac{1}{j\omega\epsilon\mu}\nabla\left(\nabla.\mathbf{A}\right)-j\omega\mathbf{A}
        \\
        \mathbf{E}_F &= -\frac{1}{\epsilon}\left(\nabla\times\mathbf{F}\right)
        \\
        \mathbf{H}_F &= \frac{1}{j\omega\epsilon\mu}\nabla\left(\nabla.\mathbf{F}\right)-j\omega\mathbf{F}
        
where :math:`\mathbf{E}_{\mathrm{tot}}=\mathbf{E}_A+\mathbf{E}_F`  and
:math:`\mathbf{H}_{\mathrm{tot}}=\mathbf{H}_A+\mathbf{H}_F`. For further details, 
please refer to (page.135, [Balanis]).
        
Far-field Approximations
------------------------

In rectangular co-ordinate system, when :math:`r \gg r'`, the term :math:`|\mathbf{r-r^\prime}|` can be approximated as

.. math::

        |\mathbf{r-r^{\prime}}| &\approx \left(r-\hat{\mathbf{r}}\cdot\mathbf{r^{\prime}}\right)
        \\
        &\approx r-\left(\sin\theta\cos\phi x'+\sin\theta\sin\phi y'+\cos\theta z'\right)
        \\
        &\approx r-\left(ux'+vy'+wz'\right).
        
So, far-field vector potentials are given as

.. math::

        \mathbf{A} &\approx \frac{\mu}{4\pi}\frac{e^{-jk_0r}}{r}\int_\mathrm{volume}\mathbf{J}_e{e^{+jk_0(ux'+vy'+wz')}}dv^\prime
        \\
        \mathbf{F} &\approx \frac{\epsilon}{4\pi}\frac{e^{-jk_0r}}{r}\int_\mathrm{volume}\mathbf{J}_m{e^{+jk_0(ux'+vy'+wz')}}dv^\prime.
        
From the above equations, it is evident that :math:`\mathbf{A} \& \mathbf{J}_e`
and :math:`\mathbf{F} \& \mathbf{J}_m` form **Fourier transform pairs**.
In comparison to the signal processing terminology, :math:`(x,y,z)` and 
:math:`(u,v,w)` are analogous to time :math:`t` and frequency :math:`f`, respectively.

Also, far-field electric and magnetic field components can be approximated as (for the 
:math:`\theta` and :math:`\phi` components only since :math:`E_r\simeq0` and :math:`H_r\simeq0`)

.. math::

        \mathbf{E}_A &\simeq -j\omega\mathbf{A}
        \\
        \mathbf{H}_F &\simeq -j\omega\mathbf{F}.
        
The far-field components :math:`\mathbf{H}_A` and :math:`\mathbf{E}_F` are related
to the above components as :math:`\mathbf{H}_A\simeq\frac{1}{\eta}\left(\hat{\mathbf{a}}_r\times\mathbf{E}_A\right)`
and :math:`\mathbf{E}_F\simeq-\eta\left(\hat{\mathbf{a}}_r\times\mathbf{H}_F\right)`,
where :math:`\eta=\sqrt{\mu / \epsilon}` is the free space wave impedance.

Far-field Green’s Functions of Infinitesimal Dipoles
----------------------------------------------------

Now, a simple example will be considered. This example deals with evaluation of
the far-field components corresponding to a infinitesimal electric dipole placed
at the origin and oriented along the :math:`z`-axis. The corresponding vector potential
is given as

.. math::

        \mathbf{A} = \left(\frac{\mu}{4\pi}\frac{e^{-jk_0r}}{r}\right)\hat{\mathbf{z}}
        
`Converting <http://en.wikipedia.org/wiki/Vector_fields_in_cylindrical_and_spherical_coordinates#Spherical_coordinate_system>`_
the above equation into spherical co-ordinate system gives

.. math::

        \mathbf{A} = \left(-\frac{\mu}{4\pi}\frac{e^{-jk_0r}}{r}\sin\theta\right)\hat{\mathbf{\theta}}.
        
In deriving the above equation, radial component of the vector potential is neglected.
Finally, far-field electric field is given as

.. math::

        \mathbf{E} = -j\omega\mathbf{A} = \left(j\omega\frac{\mu}{4\pi}\frac{e^{-jk_0r}}{r}\sin\theta\right)\hat{\mathbf{\theta}}.
        
Similar far-field Green’s functions corresponding to infinitesimal dipoles oriented along various
directions are given below.

.. image:: _static/Green.png
   :align: center
