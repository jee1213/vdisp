# vdisp
vdisp_grid_apt4.py takes the solution of the anisotropic spherical Jeans equation as input.
The code computes the luminosity weighted (by Hernquist profile), aperture-averaged 
projected velocity dispersion through fft-convolving the 2D Gaussian seeing.
grid_create.py takes this aperture-averaged velocity dispersions calculated with a
grid of parameters, (gamma, rani, beta_in, beta_out), and creates grid of velocity
dispersion values. The grid then can be interpolated at various parameter combinations,
enabling the computation of aperture-averaged, luminosity weighted projected velocity dispersion
fast for the importance sampling that comes later. 
