---

---

# Langevin equation

The *Langevin equation* describes the irregular Brownian motion of particles in _an equilibrium system_. In such a system, the particles are under the influence of a slowly varying _drag_ force ($-\zeta \bm{v}$) and a rapidly varying _random_ force ($\bm{F}_\text{R}$) due to the thermal fluctuations

$$
\frac{d \bm{v}(t)}{dt} = - \zeta \bm{v}(t) + \bm{F}_\text{R}(t).
$$

The random force means that we have the following condition

$$
\langle \bm{v}(0).\bm{F}_\text{R}(t) \rangle = 0 \quad , \quad \forall t,
$$

where the average is taken over all particles in ensemble.

> The interesting point is that although the macroscopic properties of the ensemble could reach steady state, the individual particles in the ensemble will not be in steady state.

The _non-Markovian_ shape fo the Langeving equation is

$$
\frac{d\bm{v}(t)}{dt} = - \int_0^t dt \zeta(t-t') \bm{v}(t') + \bm{F}_\text{R}(t).
$$

## Generalized Langevin equation

Many transport processes are discribed by an equation with a similar form, therefore, we call the following equation the _generalized Langeving equation_ for phase variable $A$

$$
\frac{dA(t)}{dt} = - \int_0^t dt K(t-t') A(t') + F(t),
$$

where $K(t)$ is the _time dependent_ transport coefficient.

> We usually seek to evaluate the transport coefficient, $K$.

> Also, note that the transport coefficient of the whole ensemble appears in the Langevin equation which describes the dynamics of individual particles in the ensemble.

Further more we assume that the random force, $F(t)$, is random at all times

$$
\langle A(t_0) F(t) \rangle = \langle F(t) A^*(t_0) \rangle = 0 \quad , \quad \forall t~\text{and}~t_0
$$

We define the _equilibrium autocorrelation function_

$$
C(t) = \langle A(t) A^*(0)\rangle,
$$

and  can show that 

$$
\frac{d C(t)}{dt} = - \int_0^t dt ~ K(t-t')~ C(t).
$$

<details>
<summary>Show details</summary>

We multiply the complex conjugate of the phase variable from the right to the generalized Langevin equation and take the ensemble average

$$
\langle \frac{dA(t)}{dt} A^*(0) \rangle = - \int_0^t dt ~ K(t-t')~ \langle A(t') A^*(0) \rangle + \langle F(t) A^*(0) \rangle.
$$

Using the randomness property of the force $F$, we get

$$
\frac{d\langle A(t) A^*(0)\rangle}{dt} = - \int_0^t dt ~ K(t-t')~ \langle A(t') A^*(0) \rangle.
$$

</details>

Taking the Laplace transform of this equation we get

$$
s \tilde{C}(s) - C(0) = - \tilde{K}(s) \tilde{C}(s).
$$

<details>
<summary>Show details</summary>

> Details is missing!

</details>

We can also define _flux autocorrelation function_

$$
\phi(t) = \langle \dot{A}(t) \dot{A}^*(0) \rangle.
$$

We can show the following relation between equilibrium autocorrelation function and flux autocorrelation function

$$
\frac{d^2 C(t)}{dt^2} = - \phi(t).
$$

<details>
  <summary>Show details</summary>

> The details of the appearance of the minus sign in the end result is still unclear.

$$
\frac{d^2 C(t)}{dt^2} = \frac{d}{dt} \langle \frac{dA(t)}{dt} A^*(0) \rangle = \frac{d}{dt} \langle [i\mathcal{L}A(t)] A^*(0) \rangle = \frac{d}{dt} \langle A(t) [-i\mathcal{L}A^*(0)]\rangle = \langle [i\mathcal{L}A(t)] [-i\mathcal{L}A^*(0)] \rangle = - \langle \dot{A}(t) \dot{A}^*(0) \rangle = -\phi(t).
$$

</details>

Therefore, taking the Laplace transform on both sides and assuming $\frac{dC}{dt}\big|_{t=0} = 0 $, we get

$$
-\tilde{\phi}(s) = s^2 \tilde{C}(s) - s C(0).
$$

Markdown converts text to HTML.

*[HTML]: HyperText Markup Language
<!--stackedit_data:
eyJoaXN0b3J5IjpbLTYzNDAxNTkwMSwxODE2NjM5NDM3LC0xMT
Y0NTk2MjM2XX0=
-->