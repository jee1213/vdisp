# This python script generates mock parameters of (Dd, Dds/Ds, kext, gamma, reff, rani, theta_E)to feed in to the jeans equation for B1608+656 analysis #
# among these parameters
# read in parameter file and grid velocity dispersions #
# calculate the normalized vdist #
# create grid form and make a interpolation function #
from numpy import *
import sys,math,os
from math import log10,pi
from vdmodel import sigma_model as sm
def main():
	theta_E = 1.25
	# check if all the output files are created #
	str = []
	nl = 0
	for line in open('sig_files/list','r'):
		str.append(line)
		nl += 1
	print nl
	for j in range(nl):
		string_ = str[j]
		str[j] = string_[41:-5]
	print len(str)
	num = map(int,str)
	num = sort(num)
	print num
	# finding the missing file index #
	nmiss = 0
	missings = []
	for i in range(nl):
		if i != num[i]: 
			print i
			missings.append(i)
			num = insert(num,i,i)
			nmiss += 1
	f = open('missing_idx','w')
	for i in range(nmiss):
		f.write('%s\n'%(missings[i]))
	f.close()	
	print nmiss
	print nmiss + nl
#	return
	ngam = 121
	nbout = 13
	nbin = 13
	nrani = 121
	re = 1.85
        gamma_arr = linspace(1.01,2.99,ngam)
        rani_arr = logspace(log10(0.5*re),log10(5.*re),nrani)
        bin_arr = linspace(-0.6,0.6,nbin)
        bout_arr = linspace(-0.6,0.6,nbout)
	np = 32
	nsamp = 13*13
	data_grid = zeros((ngam,nrani,nbin,nbout),'float64')
	asec_to_rad = 4.84e-6
	# i is the index of vdisp file #
	for i in range(77323):
		nline = 0
		if i%1000 ==0:
			print '...reading %d th file...'%i 	
		vdisp_file = 'sig_files/log_grid_spl_idx_fft_vdisp_out_fine_full_%d.txt'%(i)
		for line in open(vdisp_file,'r'):nline += 1
#		if i%100 ==0:
#			print '...number of lines : %d...'%nline 	
		vdisp_data = zeros((nline,3),'float64')
		count = 0
		for line in open(vdisp_file,'r'):
			if(line.find('#')==-1):
				vdisp_data[count] = line.split()
				count += 1
		gal_idx = vdisp_data[:,0]
		vdisp_norm = vdisp_data[:,1]
		vdisp_raw = vdisp_data[:,2]
		idxs = i*np
		idxf = (i+1)*np-1
		if idxf > ngam*nrani*nbin*nbout:
			idxf = ngam*nrani*nbin*nbout-1
		# nfiles is the index of param file #
		nfiles = idxs / nsamp
		nfilef = idxf / nsamp
		iis = nfiles / ngam
		jjs = nfiles % ngam 
		iif = nfilef / ngam
		jjf = nfilef % ngam
		idx_offset = idxs-nfiles*nsamp
	        file = 'param/grid_gam%03drani%03d.dat'%(iis,jjs)
		data2 = zeros((nsamp,8),'float64')
		data = zeros((nsamp,8),'float64')
		count = 0
		for line in open(file,'r'):
			if(line.find('#')==-1):
				data[count] = line.split()
				count += 1
	        file2 = 'param/grid_gam%03drani%03d.dat'%(iif,jjf)
		data2 = zeros((nsamp,8),'float64')
		count = 0
		for line in open(file2,'r'):
			if(line.find('#')==-1):
				data2[count] = line.split()
				count += 1
		if jjs != jjf:
			data = concatenate((data,data2),axis = 0)
		Dd = data[:,0]*1000000.
		gamma = data[:,1]
		kext = data[:,2]
		theta_E = data[:,3]
		mE = data[:,4]
		mass = mE*(1.-kext)
		if iis == 0 and jjs == 0:
			print mE
		from scipy.special import gamma as gfunc
		rho = -mass*theta_E**(gamma-3)*gfunc(gamma/2.)*pi**(-1.5)/gfunc(0.5*(gamma-3))
		norm = 4.*pi*rho/(3.-gamma)
		
#	        data = zeros((nsamp,8),'float64')
#	        count = 0
#	        for line in open(file,'r'):
 #      	        	if(line.find('#')==-1):
  #                      	data[count] = line.split()
   #                     	count += 1
		gal_idx = gal_idx.astype(int)
	        gamma_idx = gal_idx / (nbout*nbin*nrani)
	        rani_idx = (gal_idx / (nbout*nbin))%nrani
	        bin_idx = (gal_idx / nbout)%nbin
	        bout_idx = gal_idx % nbout
		# data grid is guaranteed to be positive #
		# check what's happening with interpolation #
		for m in range(nline):
			#data_grid[gamma_idx[m],rani_idx[m],bin_idx[m],bout_idx[m]] = (vdisp_raw[m])**2./norm[idx_offset+m]/(asec_to_rad*Dd[idx_offset+m])	
			data_grid[gamma_idx[m],rani_idx[m],bin_idx[m],bout_idx[m]] = (vdisp_norm[m])**2./norm[idx_offset+m]*(Dd[idx_offset+m])	
	print 'data grid created'
	import ndinterp
	from scipy import interpolate
	axes = {}
	axes[0] = interpolate.splrep(gamma_arr,arange(ngam),k=3,s=0)	
	axes[1] = interpolate.splrep(rani_arr,arange(nrani),k=3,s=0)	
	axes[2] = interpolate.splrep(bin_arr,arange(nbin),k=3,s=0)	
	axes[3] = interpolate.splrep(bout_arr,arange(nbout),k=3,s=0)	
	model = ndinterp.ndInterp(axes,data_grid)
	
	filesigv = 'filesigv_adriparam'
	import cPickle
	f = open(filesigv,'wb')
	cPickle.dump(model,f,2)
	f.close()
	print '1 grid interpolated and saved'

	model2 = interpolate.RegularGridInterpolator((gamma_arr,rani_arr,bin_arr,bout_arr),data_grid,method='nearest')

	filesigv = 'filesigv_adriparam_reg_grid_norm_1115'
	import cPickle
	f = open(filesigv,'wb')
	cPickle.dump(model2,f,2)
	f.close()
	print '2 grid interpolated and saved'
	return
	numpy.random.seed(10)

	c = 8.39e-10
	Grav = 6.67408e-11*(3.086e22)**-3.*1.989e30/(1.157e-5)**2.

	zl = 0.6304
	zs = 1.394

	reff = 0.58
######## number of beta_in and beta_out : each has 13 sample points in [-0.6,0.6] ########
	nbin = 13
	nbout = 13
######## number of samples in each parameter file ########
	nsamp = nbin*nbout
######## for grid creation, pick a single value ########
	Ddmax = 2100
	Ddmin = 700
	Dd_samp = numpy.random.random(nsamp)*(Ddmax-Ddmin)+Ddmin

	DdsDsmax = 0.70
	DdsDsmin = 0.30
	DdsDs_samp = numpy.random.random(nsamp)*(DdsDsmax-DdsDsmin) + DdsDsmin

	SigCrit_samp = c**2./(4.*pi*Grav*Dd_samp*DdsDs_samp)
#	SigCrit = c**2./(4.*pi*Grav*Dd*DdsDs)

	
	thE_samp = [0.826 for i in range(nsamp)]
#	thE = 0.826

	mE_samp = SigCrit_samp * pi * thE_samp * thE_samp * Dd_samp * Dd_samp * asec_to_rad * asec_to_rad 
#	mE = SigCrit * pi * thE * thE * Dd * Dd * asec_to_rad * asec_to_rad 
	# number of desired gamma / rani for the grid #
	ngam = 121 
	nrani = 121
	count = 0 
	kextfile = 'kextcounts_for_Inh_16-0623.dat'
	for line in open(kextfile,'r'): count+=1
	kext = numpy.zeros(count)
	kext_samp = numpy.zeros(nsamp)
	m = 0
	for line in open(kextfile,'r'):
		kext[m] = line
	       	m+= 1


	gam_samp = numpy.linspace(1.01,2.99,ngam)
	rani_samp = numpy.logspace(log10(0.5*reff),log10(5*reff),nrani)
	bi_samp = numpy.linspace(-0.6,0.6,nbin)
	bo_samp = numpy.linspace(-0.6,0.6,nbout)
	# give a random seed to numpy random so that one can reproduce the same results #
	for i2 in range(nsamp):
		idx = numpy.random.randint(low=0,high=m)
		kext_samp[i2] = kext[idx]
		mass_samp = (1-kext_samp)*mE_samp
	for i1 in range(ngam):
		gam = gam_samp[i1]
		for j1 in range(nrani):
			rani = rani_samp[j1]
			paramfile = 'param/grid_gam%03drani%03d.dat'%(i1,j1)

			f = open(paramfile,'w')
			f.write("# Dd, gamma, kext, theta_E, m_E, rani, beta_in, beta_out\n")

#	model_Inh = numpy.load('filesigv_Inh_100CF')
			for i3 in range(nbin):
				bin = bi_samp[i3]
				for j3 in range(nbout):
					bout = bo_samp[j3]
					idxp = i3*nbout+j3
#	sig_ref = numpy.sqrt(2.)*sm.getSigmaFromGrid2(mass_samp,thE_samp, gam_samp, rani_samp, model_Inh, Dd_samp)
					f.write("%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\n"%(Dd_samp[idxp],gam,kext_samp[idxp],thE_samp[idxp],mE_samp[idxp],rani,bin,bout))
			f.close()	
main()	
