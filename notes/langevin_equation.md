# Langevin equation

The Langevin equation describes the motion of an _ensemble_ of particles under influence of a slowly varying _drag_ force ($-\zeta \bm{v}$) and a rapidly varying _random_ force ($\bm{F}_\text{R}$) due to the thermal fluctuations

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
