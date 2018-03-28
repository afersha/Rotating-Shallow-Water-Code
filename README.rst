======================
Rotating Shallow Water 
======================

Solves rotating shallow water equations 

.. math::

   u_t =  v (f + zeta) - B_x + nu del^(a) u
   v_t = -u (f + zeta) - B_y + nu del^(a) v
   h_t = -(u h)_x - (v h)_y -(u_x + v_y) 

where

.. math::

   zeta = v_x - u_y 
   B = (u^2+v^2)/2 + Cg^2 h    
   H = 1 + h

Inputs

::

   Sin:       N x N x 3 array containing initial u, v, h, respectively
   f:         Nondim Coriolis [ie inverse Rossby number f_0*L/U]
   Cg:        Nondim GW speed  [sqrt(g*H_0)/U]
   numsteps:  Total number of timesteps 
   savestep:  Frequency, in timesteps, to save output    
   Xpin:      Structure Xpin.x, Xpin.y containing intial particle 
              positions (optional)

Outputs

::

   Sout:      Arranged as Sin, but with 4th dimension for time
   time:      Times at which output is saved
   ke:        Time series of KE
   pe:        Time series of PE
   hmov:      Movie of h field (if hmov included in output list)
   Xp:        Structure with coordinates Xp(j).x, Xp(j).y of particles,
              where j is timestep
    
Numerical details

Model is spectal, in square domain of size 2*pi x 2*pi.  Input fields must have N = 2^n, where n is an integer.  Nonlinear terms are done in physical space using dealiased product via Orszag method.  Uses AB3 timestepping with trapezoidal hyperviscosity of order a.  Timestep and hyperviscosity are set adaptively, via::

   dt = dttune dx/max(|u|) 
   nu = nutune dx^a zeta_rms.
