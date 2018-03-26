import os
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from numpy import *
from numpy.fft import fft2,ifft2
import scipy
import time
from scipy import integrate
from scipy import fftpack
from scipy import interpolate
from scipy import signal
from vdmodel import sigma_model as sm
from vdmodel import profiles
import copy
import math
import multiprocessing
from multiprocessing import pool
def asec_to_rad(asec):
	rad = asec*4.8481*10.**(-6.)
	return rad
def itp_func(x,r,sig):
	mod = interpolate.splrep(r,sig,s=0,k=1)
	val = interpolate.splev(x,mod)
	return val
def beta(rev,rani,bin,bout):
	beta = (bin*rev*rev+bout*rani*rani)/(rev*rev+rani*rani)
	return beta
def quadint(x, R, Ih):
	rx = R*cos(x)
	ry = R*sin(x)
	integrand = Ih(rx,ry)
	return integrand
def intfunc(y,R,rani,bin,bout,rsig):
	x = exp(y)
	rev = R*x
	intf = (1.-beta(rev,rani,bin,bout)*R*R/(rev*rev))*x*interpolate.splev(rev,rsig)*x/sqrt(x*x-1.)
	return 2.*R*intf
def intfunc_lin(rev,R,rani,r,sig):
	intf = (1.-beta(rev,rani)*R*R/(rev*rev))*rev*itp_func(rev,r,sig)/sqrt(rev*rev-R*R)
	return 2.*intf
def l_intfunc(x,R,rani,r,sig):
	l = exp(x)
	rev = sqrt(l*l+R*R)
	intf = l*(1.-beta(rev,rani)*R*R/(rev*rev))*itp_func(rev,r,sig)
	return 2.*intf
def printIsigma2OM(R,r,M,light_profile,lp_args,riso,infile):
    from numpy import arctan
    if type(R)==type(1.):
        R = numpy.asarray([R])
    light = light_profile(r,lp_args)
    a = lp_args

    model = interpolate.splrep(r,M*light,k=3,s=0)
    result = R*0.
    ua2 = (riso/R)**2
    t1 = (ua2+0.5)/(ua2+1.)**1.5
    for i in range(R.size):
	outfile = infile+'_%03d.txt'%(i)
        reval = logspace(log10(R[i]),log10(r[-1]),301)
        reval[0] = R[i] # Avoid sqrt(-epsilon)
        Mlight = interpolate.splev(reval,model)
        u = reval/R[i]
        u2 = u*u
        ua = ua2[i]
        K = t1[i]*(u2+ua)*arctan(((u2-1.)/(ua+1.))**0.5)/u - 0.5*(1.-1./u2)**0.5/(ua+1.)
	integrand = K*Mlight/reval
        mod = interpolate.splrep(reval,K*Mlight/reval,k=3,s=0)
        result[i] = 2.*interpolate.splint(R[i],reval[-1],mod)
	
	f = open(outfile,'w')
	for (ip,jp) in zip(reval,integrand):
	    f.write("%f\t%.10g\n"%(ip,jp))
	f.close()	
    return result

def generate_stuff(i):
    global np
#    np = 9
    for foo in range(i*np,(i+1)*np):
#    for foo in range(i*np,(i)*np+25):
        yield foo
def project_conv_avr(i):
#	fout = open('%dred.txt'%(i),'w')	
#	fout.write('inside parallel session')
#	fout.flush()
	global arc_sig, Dd, gamma, kext, theta_E, m_E, rani, n_idx, spc, fx_red,fy_red, fgauss, r_ap
	global P, ar, rr_ref_1D, x, y, fx, fy, frr,r_fft, fxx,fyy,fr, P_red
	global red_itv, eps, r_mesh,costheta_mesh, f_pos, beta_i, beta_o
	global gamma_arr, rani_arr,bin_arr,bout_arr
#	x = arange(1,n_idx,spc)/(float(200.))
#	y = arange(1,n_idx,spc)/(float(200.))
#	xx,yy = meshgrid(x,y)
#	if beta_i < 0:
#		bin = '-.%02d'%(-100*beta_i)
#	else:
#		bin = '0.%02d'%(100*beta_i)
#	if beta_o < 0:
#		bout = '-.%02d'%(-100*beta_o)
#	else:
#		bout = '0.%02d'%(100*beta_o)
#	infile = 'sol_jeans_adri/rhosig_%03dgamma_%6.4frani_%4.2fbi_%4.2fbo_%4.2f.txt'%(i+1,gamma[i],rani[i],beta_i[i],beta_o[i])
#	count = 0
#	for line in open(infile): count +=1
#	r_sig = zeros((count,2),'float64')
#	fout.write('%s'%(infile))
#	fout.flush()
#	Rani = rani[i]
	### read all the lines ###
#	m = count-1
#	for line in open(infile,'r'):
#		r_sig[m] = line.split()
#		m -= 1
	# bring the unit to arcsec, so that the projection and convolution calculation can be done in units of angle up to 100''#
#	r = r_sig[:count-10,0]/asec_to_rad(Dd[i])
#	sig = r_sig[:count-10,1]
#	sig_interp = interpolate.splrep(r,sig,s=0,k=1)
	# single array that contains ngam * nrani * betain * betao information #
	gamma_idx = i / (nbout*nbin*nrani)
	rani_idx = (i / (nbout*nbin))%nrani
	bin_idx = (i / nbout)%nbin
	bout_idx = i % nbout
	print gamma_idx
	print rani_idx
	print bin_idx
	print bout_idx
	
	gamma = gamma_arr[gamma_idx]
	rani = rani_arr[rani_idx]
	bin = bin_arr[bin_idx]
	bout = bout_arr[bout_idx]
	arc_sig_interp = arc_sig[i]
	Ih_sig_s_int = []
	l_Ih_sig_s_int = []
	delta_r = f_pos[1]-f_pos[0]
	print delta_r
	print 'projection starts\n'
	logre = log10(re)
	R_arr = logspace(-3,3,90)
	r_eval = logspace(logre-3,logre+3,300)
	x_eval = linspace(logre-3,logre+3,300)
#	for j in range(r_eval.size):
#		outfile = 'integrand/Inh_Isigsq_%03d_gamma_%6.4fk_%6.4ftE_%6.4fmE_%4.2frani_%4.2f_%03d.txt'%(i+1,gamma[i],kext[i]*10,theta_E[i],m_E[i]/10**11,rani[i],j)
#		integ = intfunc(r_eval[j:],r_eval[j],Rani,r,sig)
#		f = open(outfile,'w')
#		for (ii,jj) in zip(r_eval[j:],integ):
#			f.write("%f\t%.10g\n"%(ii,jj))
#		f.close()
#	tb_proj = time.time()
	r_ref_1D = []
	for k in range(90):
		RR = R_arr[k]
		r_ref_1D.append(RR)
		a = integrate.quad(intfunc,log(float128(1.+eps)),float128(10),args=(RR,rani,bin,bout,arc_sig_interp),epsabs=1.0e-10)
		#a = integrate.romberg(intfunc,log(float64(1.+eps)),float64(10),args=(RR,Rani,r,sig),tol=1.0e-15,divma=40)
#		b = integrate.quad(l_intfunc,-(float128(20)),float128(20),args=(RR,Rani,r,sig),epsabs=1.0e-10)
		Ih_sig_s_int.append(a[0])
#	ta_proj = time.time()
#		l_Ih_sig_s_int.append(b[0])
#	out = 'time/time_proj%d.txt'%(i)
#	f = open(out,'w')
#	f.write('%.12g'%(ta_proj-tb_proj))
#	f.flush()
#	f.close()
#	out = 'num_conv/l_num_Ih_sig_s_%d.txt'%(i)
#	f = open(out,'w')
#	for (kk,mm) in zip(r_ref_1D,l_Ih_sig_s_int):
#		f.write('%.12g\t%.12g\n'%(kk,mm))
#	f.flush()
#	f.close()
	l = len(f_pos)
	Ih_sig_s_1D = interpolate.splrep(r_ref_1D,Ih_sig_s_int,s=0,k=3)
#	I_sig_s_conv = lambda x,y: interpolate.splrep((x*x+y*y)**0.5, I_sig_s_1D)
	mm,nn = shape(frr)
	Ih_sig_s_1D_ = 0.*frr
	for ii in range(mm):
		for jj in range(nn):
			Ih_sig_s_1D_[ii,jj] = interpolate.splev(frr[ii,jj],Ih_sig_s_1D)
	print "2D array created"
	mini = amin(Ih_sig_s_1D_)
	maxi = amax(Ih_sig_s_1D_)
	#lbl = linspace(log10(mini),log10(maxi),20)
#	CS = plt.contour(fx,fy,log10(Ih_sig_s_1D_),levels = lbl)
#	plt.clabel(CS,inline=1,fontsize=10)
#	plt.colorbar(CS,ticks=lbl)
#	plt.xlim([fx[0],-fx[0]])
#	plt.ylim([fx[0],-fx[0]])
#	plt.savefig('bconv_Ih_sig.png',format='png')
#	plt.close()
	# expand 1D grid to 2D grid and create Ih_sig_s_2D function #
#	Ih_sig_s_2D = interpolate.interp2d(fx,fy,Ih_sig_s_1D_,kind='cubic')
#	Ih_sig_s_2D_ = Ih_sig_s_2D(fx,fy)
	#Ih_sig_s_2D_red = interpolate.interp2d(fx_red,fy_red,Ih_sig_s_1D_,kind='cubic')
	# mapping to a linear grid #
#	Ih_sig_s_2D = interpolate.interp2d(fx,fy,Ih_sig_s_1D_)
#	Ih_sig_s_2D_ = Ih_sig_s_2D(fx,fy)
	print "2D interpolation done"
#	fout.write('2D interpolation done\n')
#	fout.flush()
#	Ih_sig_s_conv_red = fft_convolve2d_fg(Ih_sig_s_2D_red_,fgauss)
	#Ih_sig_s_conv = fft_convolve2d(Ih_sig_s_1D_,P)
	Ih_sig_s_conv = signal.fftconvolve(Ih_sig_s_1D_,P,mode="same")*delta_r*delta_r
	mini = amin(Ih_sig_s_conv)
	maxi = amax(Ih_sig_s_conv)
	#lbl = linspace(log10(mini),log10(maxi),20)
#	CS = plt.contour(fx,fy,log10(Ih_sig_s_conv),levels=lbl)
#	plt.clabel(CS,inline=1,fontsize=10)
#	plt.colorbar(CS,ticks=lbl)
#	plt.xlim([fx[0],-fx[0]])
#	plt.ylim([fx[0],-fx[0]])
#	plt.savefig('conv_Ih_sig.png',format='png')
#	plt.close()
	print "Fourier convolution done"
	Ih_sig_s_conv_2D = interpolate.interp2d(fx,fy,Ih_sig_s_conv)
	
	Ih_sig_s_avr = f_pos*0.
	m,n = shape(Ih_sig_s_conv)
#	for kk in range(l):
#		R = f_pos[kk]
#		a = integrate.quad(quadint,0.,2.*pi,args=(R,Ih_sig_s_conv_2D))
#		Ih_sig_s_avr[kk] = a[0]/(2.*pi)
#	ta_angavr = time.time()	
#	print "angle average done"
	r_ref_1D = asarray(r_ref_1D)
	Ih_sig_s_int = asarray(Ih_sig_s_int)
#	plt.plot(f_pos,f_pos*f_pos*Ih_sig_s_avr,label='ang avr')
#	plt.plot(r_ref_1D,r_ref_1D*r_ref_1D*Ih_sig_s_int,label='before conv')
#	plt.xlim([0.01,1])
#	plt.xscale('log')
#	plt.yscale('log')
#	plt.legend()
#	plt.savefig('rsq_ang_avr_test.png',format='png')
#	plt.close()
#	plt.plot(f_pos,Ih_sig_s_avr,label='ang avr')
#	plt.plot(r_ref_1D,Ih_sig_s_int,label='before conv')
#	plt.xlim([0.01,1])
#	plt.xscale('log')
#	plt.yscale('log')
#	plt.legend()
#	plt.savefig('ang_avr_test.png',format='png')
#	plt.close()
#	Ih_sig_s_2D_ft = fftw(Ih_sig_s_2D)
#	P_ft = fftw(P)
#	Ih_sig_s_conv_ft = Ih_sig_s_2D_ft*P_ft
#	Ih_sig_s_cov = ifft2(Ih_sig_s_conv_ft)
	# do the convolution on 2D grid function, then integrate # 
#	Ih_sig_s_conv = signal.convolve2d(Ih_sig_s_2D_,P,mode='same')
#	fout.write('convolution done\n')
#	fout.flush()
#	Ih_sig_s_conv_ = interpolate.interp2d(fx,fy,Ih_sig_s_conv,kind='cubic')
#	sig_ap = integrate.dblquad(Ih_sig_s_conv_,-n_idx,n_idx,ylower,yupper)
	#Ih_sig_s_conv_red_ = interpolate.interp2d(fx_red,fy_red,Ih_sig_s_conv_red,kind='cubic')

#	tb_2dinterp = time.time()	
	Ih_sig_s_conv_ = interpolate.interp2d(fx,fy,Ih_sig_s_conv)
#	ta_2dinterp = time.time()	
#	out = 'time/time_2dinterp%d.txt'%(i)
#	f = open(out,'w')
#	f.write('%.12g'%(ta_2dinterp-tb_2dinterp))
#	f.flush()
#	f.close()
	print "interpolate again"
#	tb_apavr = time.time()	
	sig_ap = integrate.dblquad(Ih_sig_s_conv_,-r_ap,r_ap,ylower,yupper,epsabs=1e-10)
#	ta_apavr = time.time()	
#	out = 'time/time_apavr%d.txt'%(i)
#	f = open(out,'w')
#	f.write('%.12g'%(ta_apavr-tb_apavr))
#	f.flush()
#	f.close()
	print "calculate aperture average"
	#Ih_sig_s_1D = interpolate.interp1d(r_ref_1D,Ih_sig_s_int)
	# expand 1D grid to 2D grid and create Ih_sig_s_2D function #
	#Ih_sig_s_2D = interpolate.interp2d(x,y,Ih_sig_s_1D(rr),kind='cubic')
	# do the convolution on 2D grid function, then integrate # 
	#Ih_sig_conv = scipy.signal.convolve2d(Ih_sig_s_2D,P,mode='same')
	#sig_ap = integrate.simps(integrate.simps(Ih_sig_cov,x),y)	
#	fout.write('%f'%(sig_ap[0]))
#	fout.close
	print "%d core done"%(i)
	print sig_ap[0]
	return (i,sig_ap[0])
def makeGaussian(fx, steps, sig , center=None):
    """ Make a square gaussian kernel.
    size is the length of a side of the square
    fwhm is full-width-half-maximum, which
    can be thought of as an effective radius.
    """
    fxx, fyy = meshgrid(fx,fx)
#    x = arange(-size, size, steps, float)
#    y = arange(-size, size, steps, float)
#    xx, yy = meshgrid(x,y)
    if center is None:
#        x0 = y0 = size // 2
	x0 = y0 = 0.
    else:
        x0 = center[0]
        y0 = center[1]
    rr = ((fxx-x0)*(fxx-x0) + (fyy-y0)*(fyy-y0))
    G = exp(- (rr) / (2.*sig*sig))/(2.*math.pi*sig*sig)
    return exp(- (rr) / (2.*sig*sig))/(2.*math.pi*sig*sig)
#    return exp(- (rr) / (2.*sig*sig))

def makeFourierGaussian(size, steps, sig , center=None):
    """ Make a square gaussian kernel.
    size is the length of a side of the square
    fwhm is full-width-half-maximum, which
    can be thought of as an effective radius.
    """
    a = 1./(2.*sig*sig)
    x = arange(-size, size, 10, float)
    y = arange(-size, size, 10, float)
    xx, yy = meshgrid(x,y)
    if center is None:
#        x0 = y0 = size // 2
	x0 = y0 = 0.
    else:
        x0 = center[0]
        y0 = center[1]
    rr = ((xx-x0)*(xx-x0) + (yy-y0)*(yy-y0))
    return exp(-pi**2. * (rr) / (a))

def fft_convolve2d(x,y):
    """ 2D convolution, using FFT"""
    fr = fft2(x)
    fr2 = fft2(flipud(fliplr(y)))
    m,n = fr.shape
    cc = real(ifft2(fr*fr2))
    cc = roll(cc, -m/2+1,axis=0)
    cc = roll(cc, -n/2+1,axis=1)
    return cc

def fft_convolve2d_fg(x,y):
    """ 2D convolution, using FFT"""
    global n_fft
    global sig
    fr = fft2(x)
#    fr2 = fft2(flipud(fliplr(y)))
    fr2 = fft.ifftshift(y)
    m,n = fr.shape
    cc = real(ifft2(fr*fr2))
    cc = roll(cc, -m/2+1,axis=0)
    cc = roll(cc, -n/2+1,axis=1)
    return cc


def hernquist(x):
	# x is R/a #
        if x <= 1:
                Ih = 1./(1.-x**2.)**2.*((2.+x**2.)/math.sqrt(-x**2.+1.)*arccosh(1./x)-3.)
        else:
                Ih = 1./(1.-x**2.)**2.*((2.+x**2.)/math.sqrt(x**2.-1.)*arccos(1./x)-3.)
        return Ih

def gaussian(x,mu,sig):
	return 1./((2.*pi)*sig*sig)*exp(((x - mu) / (2. * sig ))** 2.)

def idl_tabulate(x, f, p=5) :
    def newton_cotes(x, f) :
        if x.shape[0] < 2 :
            return 0
        rn = (x.shape[0] - 1) * (x - x[0]) / (x[-1] - x[0])
        weights = scipy.integrate.newton_cotes(rn)[0]
        return (x[-1] - x[0]) / (x.shape[0] - 1) * dot(weights, f)
    ret = 0
    for idx in xrange(0, x.shape[0], p - 1) :
        ret += newton_cotes(x[idx:idx + p], f[idx:idx + p])
    return ret
def ylower(x):
#	global n_ap
	#return -0.35*n_ap
	return -0.35*4
def yupper(x):
#	global n_ap
	return 0.35*4

def main():
	ngam = 121
	nrani = 2
	nbin = 3
	nbout = 3
	re = 1.85
#	n_ap = input('input n_ap : ')
	n_ap = 4
	global ngam,nrani,nbin,nbout,re
	gamma_arr = linspace(1.01,2.99,ngam) 
	rani_arr = logspace(log10(0.5*re),log10(5*re),nrani) 
	bin_arr = linspace(-0.6,0.6,nbin)
	bout_arr = linspace(-0.6,0.6,nbout) 
	input = int(sys.argv[1])
	global np 
	global n_ap 
	np = multiprocessing.cpu_count()
	#npoint = 9
	print np
	## indexs is the starting number of galaxy from the whole sample ##
	## indexf is the finishing number of galaxy from the whole sample ##
	indexs = input * np
	indexf = (input+1) * np
	nsamp = nbin*nbout
	nfiles = indexs / nsamp
	iis = nfiles / nrani
	jjs = nfiles % nrani
	idx_offset = indexs - nfiles*nsamp
	## iis*nsamp + jjs is the starting index of the file #
	## idx_offset is the starting index of the processor in the file #
	nfilef = indexf / nsamp
	dnfile = nfilef-nfiles	
	iif = nfilef / nrani
	jjf = nfilef % nrani
#	print ii
#	jj = int(sys.argv[1])
	if iis+1 == iif:
		file = 'param_apt/grid_gam%03drani%03d.dat'%(iis,jjs)
		data0 = zeros((nsamp,8),'float64')
		count = 0
		for line in open(file,'r'):
			if(line.find('#')==-1):
				data0[count] = line.split()
				count += 1
				Dd0 = data0[:,0]*1000000.
		file1 = 'param_apt/grid_gam%03drani%03d.dat'%(iis,jjs+1)
		data1 = zeros((nsamp,8),'float64')
		count = 0
		for line in open(file1,'r'):
			if(line.find('#')==-1):
				data1[count] = line.split()
				count += 1
				Dd1 = data1[:,0]*1000000.
		file2 = 'param_apt/grid_gam%03drani%03d.dat'%(iif,jjs)
		data2 = zeros((nsamp,8),'float64')
		count = 0
		for line in open(file2,'r'):
			if(line.find('#')==-1):
				data2[count] = line.split()
				count += 1
				Dd2 = data2[:,0]*1000000.
		file3 = 'param_apt/grid_gam%03drani%03d.dat'%(iif,jjf)
		data3 = zeros((nsamp,8),'float64')
		count = 0
		for line in open(file3,'r'):
			if(line.find('#')==-1):
				data3[count] = line.split()
				count += 1
				Dd3 = data3[:,0]*1000000.
	else:
		file = 'param_apt/grid_gam%03drani%03d.dat'%(iis,1)
		data0 = zeros((nsamp,8),'float64')
		count = 0
		for line in open(file,'r'):
			if(line.find('#')==-1):
				data0[count] = line.split()
				count += 1
				Dd0 = data0[:,0]*1000000.
		file1 = 'param_apt/grid_gam%03drani%03d.dat'%(iis+1,0)
		data1 = zeros((nsamp,8),'float64')
		count = 0
		for line in open(file1,'r'):
			if(line.find('#')==-1):
				data1[count] = line.split()
				count += 1
				Dd1 = data1[:,0]*1000000.
		file2 = 'param_apt/grid_gam%03drani%03d.dat'%(iis+1,1)
		data2 = zeros((nsamp,8),'float64')
		count = 0
		for line in open(file2,'r'):
			if(line.find('#')==-1):
				data2[count] = line.split()
				count += 1
				Dd2 = data2[:,0]*1000000.
		file3 = 'param_apt/grid_gam%03drani%03d.dat'%(iif,0)
		data3 = zeros((nsamp,8),'float64')
		count = 0
		for line in open(file3,'r'):
			if(line.find('#')==-1):
				data3[count] = line.split()
				count += 1
				Dd3 = data3[:,0]*1000000.
		if dnfile == 4:
			if jjs == 0:			
				file4 = 'param_apt/grid_gam%03drani%03d.dat'%(iis,jjs)
				data4 = zeros((nsamp,8),'float64')
				count = 0
				for line in open(file4,'r'):
					if(line.find('#')==-1):
						data4[count] = line.split()
						count += 1
						Dd4 = data4[:,0]*1000000.
			if jjs == 1:
				file4 = 'param_apt/grid_gam%03drani%03d.dat'%(iif,jjf)
				data4 = zeros((nsamp,8),'float64')
				count = 0
				for line in open(file4,'r'):
					if(line.find('#')==-1):
						data4[count] = line.split()
						count += 1
						Dd4 = data4[:,0]*1000000.
	print file
	# read data for the individual objects #
	global Dd, gamma, kext, theta_E, m_E, rani, n_idx, spc, beta_i, beta_o
	global P, fx_red,fy_red, P_red, n_fft, sig, fgauss
	global x,y,fx,fy,frr,r_fft, fxx,fyy,fr, r_ap, r_mesh,costheta_mesh
	global red_itv, eps, f_pos
	global arc_sig
	global gamma_arr, rani_arr,bin_arr,bout_arr
	arc_sig = load('jeans_grid')
	## Mpc to pc
	data = row_stack((data0,data1,data2,data3))
	if dnfile== 4:
		if jjs ==0:
			data = row_stack((data4,data))
		if jjs ==1:
			data = row_stack((data,data4))
	Dd = data[:,0]*1000000.
#	Dd = 1400*1000000.
	gamma = data[:,1]
	kext = data[:,2]
	theta_E = data[:,3]
	m_E = data[:,4]
	rani = data[:,5]
	beta_i = data[:,6]
	beta_o = data[:,7]
	# create Gaussian array for convolution #
#	P = linspace(-1,1,2000)
	# seeing in distance unit #
	seeing = 0.7
	# change FWHM to std #
	sig = seeing / 2.355
#	sig = seeing / 3.33
	#sig = seeing 
	# create symmetric 2D Gaussian kernel, with distance unit (Mpc) #
	# in a square grid of 2x2 arcsec, with 0.005 arcsec spacing #
	# create a mesh where we do the 2D convolution of Gaussian seeing, in the first quadrature #
	# n_idx = 0.405/0.05
	# spc = spacing between indices
	reff = 1.85
	a = reff/1.8153
	light = profiles.hernquist
	logre = log10(reff)
	r_eval = logspace(logre-3,logre+3,300)
	limit = r_eval[-1]
	Rvals = r_eval[r_eval<=limit]
	Rvals = r_eval[:-1]	
#	for i1 in range(count):
#		outfile = 'integrand/Matt_Isigsq_%03d_gamma_%6.4fk_%6.4ftE_%6.4fmE_%4.2frani_%4.2f'%(i1+1,gamma[i1],kext[i1]*10,theta_E[i1],m_E[i1]/10**11,rani[i1])
#		print outfile
#		M = r_eval**(3-gamma[i1])
#		sbSigma = printIsigma2OM(Rvals,r_eval,M,light,a,rani[i1],outfile)
		#integrand = interpolate.splev(r_eval,mod0)
	eps = float64(1.0e-7)
	red_itv = 1
	delta_r = 0.01*red_itv
	r_fft = 5.
	r_ap = 0.405*n_ap
	n_fft = int(ceil(r_fft/delta_r))
	n_idx = int(ceil(r_ap/delta_r))
	n_diag = int(ceil(n_idx*1.6))
	delta = 1./500.
	n_sec = 1./delta_r
	spc = 1
	p_arr = linspace(0.5,n_fft+0.5,n_fft+1)*delta_r
	#p_arr = logspace(log10(eps/delta_r),log10(n_fft+eps/delta_r),n_fft+1)*delta_r
	w_arr = append(-p_arr[::-1],p_arr)
        fx = w_arr
        fy = w_arr
	fr = sqrt(fx*fx+fy*fy)
        fxx,fyy = meshgrid(fx,fy)
        frr = sqrt(fxx*fxx+fyy*fyy)
	f_pos = p_arr
#        fgauss = makeFourierGaussian(n_fft,1,sig)
	P = makeGaussian(fx,delta_r,sig)
	print P[0,0]
	print P[n_fft+1,n_fft+1]
	# create Ih = hernquist/a^2, a = 0.551 * reff, reff = 1.85 ''#
	global ar, rr_ref_1D
	global delta_r
	#p_arr = logspace(log10(eps/delta_r),log10(eps/delta_r+2.0*n_fft),2.0*n_fft+1)*delta_r
	p_arr = linspace(0.5,0.5+2.0*n_fft,2.0*n_fft+1)*delta_r
#	Ih = [hernquist(k/(0.551*reff))/(0.551*reff*asec_to_rad(Dd))**2. for k in p_arr] 
	Ih = [hernquist(k/(0.551*reff))/(0.551*reff)**2. for k in p_arr] 
	Ih_rev = Ih[::-1]
	Ih = append(Ih_rev,Ih)
        rr_ref_1D = p_arr
        rr_ref_1D_rev = -rr_ref_1D[::-1]
        rr_ref_1D = append(rr_ref_1D_rev,rr_ref_1D)
	# relate Ih_sig_s_int in rrint grid and create Ih_sig_s_1D function #
#	Ih_1D = interpolate.interp1d(rr_ref_1D,Ih)
#	Ih_1D_ = Ih_1D(frr)
#	Ih_1D_red = interpolate.interp1d(rr_ref_1D_red,Ih_red)
#	Ih_1D_red_ = Ih_1D_red(frr_red)
	Ih_1D = interpolate.splrep(rr_ref_1D,Ih,s=0,k=3)
	Ih_1D_ = 0.*frr
	m,n = shape(frr)
	for i in range(m):
		for j in range(n):
			Ih_1D_[i,j] = interpolate.splev(frr[i,j],Ih_1D)
#	Ih_1D_red_ = interpolate.splev(frr_red,Ih_1D_red)
	# expand 1D grid to 2D grid and create Ih_sig_s_2D function #
#	Ih_2D = interpolate.interp2d(fx,fy,Ih_1D_,kind='cubic')
#	Ih_2D_ = Ih_2D(fx,fy)
	#Ih_2D_red = interpolate.interp2d(fx_red,fy_red,Ih_1D_red_,kind='cubic')
	Ih_2D = interpolate.interp2d(fx,fy,Ih_1D_)
	Ih_2D_ = Ih_2D(fx,fy)
	# do the convolution on 2D grid function, then integrate # 
#	Ih_conv = fft_convolve2d(Ih_2D_,P)
#	Ih_conv_red = fft_convolve2d_fg(Ih_2D_red_,fgauss)
	#Ih_conv = fft_convolve2d(Ih_2D_,P)
	Ih_conv = signal.fftconvolve(Ih_2D_,P,mode="same")*delta_r*delta_r
#	r_mesh = frr*0.
#	costheta_mesh = frr*0.
#	m,n = shape(frr)
#	for iii in range(m):
#		for jj in range(n):
 #			r_mesh[iii,jj] = sqrt(fx[iii]*fx[iii]+fy[jj]*fy[jj])
#			costheta_mesh[iii,jj] = fx[iii]/r_mesh[iii,jj]
#	CS = plt.contour(fx_red,fy_red,fgauss)
#	plt.clabel(CS,inline=1,fontsize=10)
#	plt.savefig('conv_Hernquist.png',format='png')
#	plt.close()
#	Ih_conv = signal.convolve2d(Ih_2D_,P,mode='same')
	#Ih_conv_red_ = interpolate.interp2d(fx_red,fy_red,Ih_conv_red,kind='cubic')
	Ih_conv_ = interpolate.interp2d(fx,fy,Ih_conv)
#	denom_ap = integrate.dblquad(Ih_conv_,-n_idx,n_idx,ylower,yupper)
	denom_ap = integrate.dblquad(Ih_conv_,-r_ap,r_ap,ylower,yupper,epsabs=1e-10)
	print "denominator"
	print denom_ap
	# Ih normalization Dd^-2 : denominator varies depending on gal idx #
	denom_ap_i = denom_ap[0]/(asec_to_rad(Dd))**2.

	# first integration should contain the full array as it requires integration from R to infty # 
	# do a 2d interpolation to create a grid to convolve with 2d Gaussian #
	# units in distance or angular size? #
	pool = multiprocessing.Pool(processes = np)
	print 'starting the parallel session'
	# ii is the one-dimensional index in parameter space #
	results = pool.map(func=project_conv_avr,iterable=generate_stuff(input))
	pool.close()
	print 'finishing the parallel session'
	idx = []
	sig_s = []
	for items in results:
		idx.append(items[0])
		sig_s.append(items[1])
	#sig_ss = [math.sqrt(sig_s[k]/denom_ap[0]*asec_to_rad(Dd[idx_offset+k])**3.) for k in range(np)]
	sig_ss = [math.sqrt(sig_s[k]/denom_ap[0]*asec_to_rad(Dd[idx_offset+k])**3.) for k in range(np)]
	#sig_wodd = [math.sqrt(sig_s[k]/denom_ap[0]) for k in range(np)]
	sig_wodd = [math.sqrt(sig_s[k]/denom_ap[0]) for k in range(np)]
#	sig_norm = [math.sqrt(i/sig_s[0])*sig_ref[ii*np] for i in sig_s]
#	sig_print = []
#	sig_s = [x/denom_ap[0] for x in sig_s]
#    	for i in range(ii*np,(ii+1)*np):
#		sig_print.append(sig_ref[i])
#	fout = open('log_grid_spl_idx_fft_vdisp_out_fine_norm_%d.txt'%(ii),'w')	
#	for i,j,k in zip(idx,sig_norm,sig_print):
#		fout.write("%d\t%.8g\t%.8g\n"%(i,j,k))	
#	fout.close
	fout = open('sig_apt4/log_grid_spl_idx_fft_vdisp_out_fine_full_%d.txt'%(input),'w')	
#	fout = open('sig_files/log_grid_spl_idx_fft_vdisp_out_fine_full_77322.txt','w')	
	for i,j,k in zip(idx,sig_ss,sig_wodd):
		fout.write("%d\t%.8g\t%.8g\n"%(i,j,k))	
	fout.close
	return
main()



