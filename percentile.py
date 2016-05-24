import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import healpy as hp
from sklearn import linear_model

###########################################################################
### read eboss file ###
###########################################################################
hdulist = fits.open('systematics-v0005-lrg.fits')
tb = hdulist[1].data


'''
###########################################################################
##############   		devide galaxies into pixels
###########################################################################
pixel_set = set(tb['IPIX128'])

n_gal = []
n_star = []
i = 0
ra_pixel = []
dec_pixel = []
for x in pixel_set:
	pixel_select = tb['IPIX128']==x
	n_gal.append(tb['TARGET_DENSITY'][pixel_select][0])
	n_star.append(tb['STAR_DENSITY'][pixel_select][0])
	ra_pixel.append(hp.pix2ang(128,x)[0])
	dec_pixel.append(hp.pix2ang(128,x)[1])
	print i*1.0/len(pixel_set)
	i += 1

n_gal = np.array(n_gal)
n_star = np.array(n_star)
ra_pixel = np.array(ra_pixel)*180/np.pi
dec_pixel = np.array(dec_pixel)*180/np.pi

###########################################################################
##############     bin the n_star
###########################################################################
n_gal_star = []
star_range = np.arange(0.0,10000,1000)
for i in range(len(star_range)):
	star_select = (n_star < star_range[i]+1000) & (n_star > star_range[i])
	n_gal_star.append(np.mean(n_gal[star_select]))

###########################################################################
##############     plot
###########################################################################

plt.plot(n_star,n_gal,'.')
plt.plot(star_range+500, n_gal_star,'ro-',linewidth=2)
plt.axis([0.0,10000,0,300])
plt.show()

##### plot the position of n_star bin

star_select = (n_star < 1800) & (n_star > 1500)
plt.plot(dec_pixel,90-ra_pixel,'.')
plt.plot(dec_pixel[star_select],90-ra_pixel[star_select],'r.')
plt.savefig("n_star9000.png")
plt.show()

###### plot percentile 

star10000 = n_star < 10000
hist = plt.hist(n_star[star10000],200)
cusum = np.cumsum(hist[0])
indx = 0
n_gal_star = []
n_gal_star_std = []
x_axis = []
for i in range(20):
	x_axis.append(hist[1][indx])
	p_select = n_star > hist[1][indx]
	while (cusum[indx] < (i+1) * cusum[-1] / 20):	indx += 1
	s_select = n_star < hist[1][indx]
	n_gal_star.append(np.mean(n_gal[p_select & s_select]))
	n_gal_star_std.append(np.std(n_gal[p_select & s_select]) / np.sqrt(len(n_gal[p_select & s_select])))

x_axis.append(10000.0)
x_axis = np.array(x_axis)

xx = (x_axis[:-1] + x_axis[1:])/2


plt.plot(xx[1:],n_gal_star[1:],'o-',linewidth=2)
plt.errorbar(xx[1:],n_gal_star[1:],xerr = 0, yerr = n_gal_star_std[1:])
plt.axis([0,10000,52,65])
plt.savefig("percentile_n_star_n_gal.png")
plt.show()
'''
###########################################################################
##############     use mangle file
###########################################################################

ra = []
dec = []
idd = []
i=0
for line in open('radec_outfile'):
	if i>0:
		columns = line.strip().split()
		ra.append(float(columns[0]))
		dec.append(float(columns[1]))
		if len(columns)>2:
			idd.append(int(columns[2]))
		else:
			idd.append(10000)
	i += 1

ra = np.array(ra)
dec = np.array(dec)
idd = np.array(idd)

idd_select = idd <10000

n_gal_mangle = []
n_star_mangle = []
i = 0
ra_mangle = []
dec_mangle =[]
psf_mangle = []
flux_mangle = []
extinction = []
W1MOON = []
W1cov = []
pixel_set = set(tb['IPIX128'][idd_select])
for x in pixel_set:
	pixel_select = tb['IPIX128']==x
	n_gal_mangle.append(tb['TARGET_DENSITY'][idd_select & pixel_select][0])
	n_star_mangle.append(tb['STAR_DENSITY'][idd_select & pixel_select][0])
	psf_mangle.append(np.mean(tb['PSF_FWHM'][idd_select & pixel_select,3]))
	flux_mangle.append(np.mean(tb['SKYFLUX'][idd_select & pixel_select,4]))
	extinction.append(np.mean(tb['EB_MINUS_V'][idd_select & pixel_select]))
	W1MOON.append(np.mean(tb["WISE_MOONLEV_W1"][idd_select & pixel_select]))
	W1cov.append(np.mean(tb['WISE_W1COVMEDIAN'][idd_select & pixel_select]))
	ra_mangle.append(hp.pix2ang(128,x)[0])
	dec_mangle.append(hp.pix2ang(128,x)[1])
	print i*1.0/len(pixel_set)
	i += 1

psf_mangle = np.array(psf_mangle)
flux_mangle = np.array(flux_mangle)
ra_mangle = np.array(ra_mangle)*180/np.pi
dec_mangle = np.array(dec_mangle)*180/np.pi
n_gal_mangle = np.array(n_gal_mangle)
n_star_mangle = np.array(n_star_mangle)
extinction = np.array(extinction)
W1MOON = np.array(W1MOON)
W1cov = np.array(W1cov)

###########################################################################
##############     plot percentile
###########################################################################
def percentile_average(n_star_mangle,n_gal_mangle):
	star10000 = n_star_mangle < 10000
	hist_mangle = plt.hist(n_star_mangle[star10000],200)
	cusum_mangle = np.cumsum(hist_mangle[0])
	indx = 0
	n_gal_star_mangle = []
	n_gal_star_mangle_std = []
	x_axis = []
	for i in range(20):
		x_axis.append(hist_mangle[1][indx])
		p_select = n_star_mangle > hist_mangle[1][indx]
		while (cusum_mangle[indx] < (i+1) * cusum_mangle[-1] / 20):	indx += 1
		s_select = n_star_mangle < hist_mangle[1][indx]
		n_gal_star_mangle.append(np.mean(n_gal_mangle[p_select & s_select]))
		n_gal_star_mangle_std.append(np.std(n_gal_mangle[p_select & s_select]) / np.sqrt(len(n_gal_mangle[p_select & s_select])))
	x_axis.append(max(n_star_mangle[star10000]))
	x_axis = np.array(x_axis)
	xxx = (x_axis[:-1] + x_axis[1:])/2
	plt.clf()
	return xxx, n_gal_star_mangle, n_gal_star_mangle_std


x_axis, n_gal_star_mangle, n_gal_star_mangle_std = percentile_average(n_star_mangle,n_gal_mangle)

x_psf, n_gal_psf, n_gal_psf_std = percentile_average(psf_mangle, n_gal_mangle)
plt.clf()

plt.plot(x_axis[1:],n_gal_star_mangle[1:],'o-',linewidth=2)
plt.errorbar(x_axis[1:],n_gal_star_mangle[1:], xerr=0, yerr=n_gal_star_mangle_std[1:])
plt.axis([0,10000,52,65])
plt.xlabel("n_star")
plt.ylabel("n_gal")
plt.savefig("percentile_n_star_n_gal_mangle.png")
plt.show()

plt.plot(x_psf[1:],n_gal_psf[1:],'o-',linewidth=2)
plt.errorbar(x_psf[1:],n_gal_psf[1:], xerr=0, yerr=n_gal_psf_std[1:])
plt.axis([1.0,2.2,59,62])
plt.xlabel("psf")
plt.ylabel("n_gal")
plt.savefig("percentile_psf_n_gal10.png")
plt.show()
'''
star_select = (n_star_mangle < 2700) & (n_star_mangle > 2400)
plt.plot(dec_mangle,90-ra_mangle,'.')
plt.plot(dec_mangle[star_select],90-ra_mangle[star_select],'r.')
plt.savefig("n_star2700.png")
plt.clf()
'''
###########################################################################
##############     n_gal distribution
###########################################################################
ra_bin = np.arange(0,361,1)
dec_bin = np.arange(-25,91,1)
ext_grid = []

for i in range(len(dec_bin)-1):
    	psf_aa = []
        for j in range(len(ra_bin)-1):
                psf_select = (dec_mangle <= ra_bin[j+1]) & (dec_mangle > ra_bin[j]) & (90-ra_mangle <= dec_bin[i+1]) & (90-ra_mangle > dec_bin[i])
                psf_aa.append(np.mean(n_gal_mangle[psf_select]))
        ext_grid.append(np.array(psf_aa))
        print i

ext_grid = np.array(ext_grid)
plt.imshow(ext_grid[-1:0:-1],clim=(40.0,90.0), cmap=plt.cm.Reds,extent=[0,360,-25,90])
plt.axis('tight')
plt.colorbar(orientation='horizontal')
plt.xlabel('ra')
plt.ylabel('dec')
plt.title('eBOSS galaxies density distribution')
plt.savefig("n_gal_distribution.png")
plt.show()


###########################################################################
##############     test random file
###########################################################################

hdulist_random = fits.open('eboss_randoms.fits')
tb_random = hdulist_random[1].data
ra_random = tb_random['RA']*np.pi/180
dec_random  = (90-tb_random['DEC'])*np.pi/180

random_IPIX = hp.ang2pix(128,dec_random,ra_random)
n_random = []
i = 0
for x in pixel_set:
	random_select = random_IPIX == x
	n_random.append(len(random_IPIX[random_select]))
	print i*1.0/len(pixel_set)
	i += 1

n_random = np.array(n_random)

###########################################################################
##############     distribution of n_random 
###########################################################################

ra_bin = np.arange(0,361,1)
dec_bin = np.arange(-25,91,1)
ext_grid = []

for i in range(len(dec_bin)-1):
    	psf_aa = []
        for j in range(len(ra_bin)-1):
                psf_select = (dec_mangle <= ra_bin[j+1]) & (dec_mangle > ra_bin[j]) & (90-ra_mangle <= dec_bin[i+1]) & (90-ra_mangle > dec_bin[i])
                psf_aa.append(np.mean(n_random[psf_select]))
        ext_grid.append(np.array(psf_aa))
        print i

ext_grid = np.array(ext_grid)
plt.imshow(ext_grid[-1:0:-1],clim=(0.0,250.0), cmap=plt.cm.Reds,extent=[0,360,-25,90])
plt.axis('tight')
plt.colorbar(orientation='horizontal')
plt.xlabel('ra')
plt.ylabel('dec')
plt.title('eBOSS galaxies density distribution')
#plt.savefig("n_gal_distribution.png")
plt.show()

###########################################################################
##############     n_gal / n_random
###########################################################################

nn_gal = n_gal_mangle[n_random!=0] / n_random[n_random!=0]
nn_gal = nn_gal / np.mean(nn_gal)


xx, n, n_std = percentile_average(n_star_mangle[n_random!=0], nn_gal)

xx_psf, nnpsf, nnpsf_std = percentile_average(psf_mangle[n_random!=0], nn_gal)

plt.plot(xx[1:],n[1:],'o-',linewidth=2)
plt.errorbar(xx[1:],n[1:], xerr=0, yerr=n_std[1:])
plt.axis([0,10000,0.5,1.5])
plt.xlabel("n_star")
plt.ylabel("n_gal / n_random / <n_gal / n_random>")
plt.savefig("percentile_n_star_n_gal_random.png")
plt.show()

plt.plot(xx_psf[1:],nnpsf[1:],'o-',linewidth=2)
plt.errorbar(xx_psf[1:],nnpsf[1:], xerr=0, yerr=nnpsf_std[1:])
plt.axis([1.0,2.2,0.5,1.5])
plt.xlabel("psf")
plt.ylabel("n_gal")
plt.savefig("percentile_psf_n_gal_random.png")
plt.show()


################   flux   ##################
x_flux, n_gal_flux, n_gal_flux_std = percentile_average(flux_mangle[n_random!=0], nn_gal)
plt.clf()

plt.plot(x_flux[1:],n_gal_flux[1:],'o-',linewidth=2)
plt.errorbar(x_flux[1:],n_gal_flux[1:], xerr=0, yerr=n_gal_flux_std[1:])
plt.axis([10.0,70.0,0.5,1.5])
plt.xlabel("z-band skyflux")
plt.ylabel("n_gal")
plt.savefig("percentile_flux_n_gal_random.png")
plt.show()

###########################################################################
##############     linear regression
###########################################################################
NGC_select = (ra_mangle < 300) & (ra_mangle > 60)
SGC_select = (ra_mangle > 300) | (ra_mangle < 60)

clfNGC = linear_model.LinearRegression()
clfSGC = linear_model.LinearRegression()
X_NGC = np.array([n_star_mangle[NGC_select],psf_mangle[NGC_select],flux_mangle[NGC_select],extinction[NGC_select],W1MOON[NGC_select],W1cov[NGC_select]]).transpose()
Y_NGC = n_gal_mangle[NGC_select]
X_SGC = np.array([n_star_mangle[SGC_select],psf_mangle[SGC_select],flux_mangle[SGC_select],extinction[SGC_select],W1MOON[SGC_select],W1cov[SGC_select]]).transpose()
Y_SGC = n_gal_mangle[SGC_select]
clfNGC.fit(X_NGC, Y_NGC)
clfSGC.fit(X_SGC,Y_SGC)

PSD_NGC = clfNGC.predict(X_NGC)
PSD_SGC = clfSGC.predict(X_SGC)
'''
RSD_NGC_star = n_gal_mangle[NGC_select] - PSD_NGC + clfNGC.coef_[0] * X_NGC[:,0]
RSD_SGC_star = n_gal_mangle[SGC_select] - PSD_SGC - clfSGC.coef_[1] * X_SGC[:,0]
RSD_NGC_psf = n_gal_mangle[NGC_select] - clfNGC.intercept_ - clfNGC.coef_[0] * n_star_mangle[NGC_select] - clfNGC.coef_[2] * flux_mangle[NGC_select]
RSD_SGC_psf  =n_gal_mangle[SGC_select] - clfSGC.intercept_ - clfSGC.coef_[0] * n_star_mangle[SGC_select] - clfSGC.coef_[2] * flux_mangle[SGC_select]
RSD_NGC_flux = n_gal_mangle[NGC_select] - clfNGC.intercept_ - clfNGC.coef_[1] * psf_mangle[NGC_select] - clfNGC.coef_[0] * n_star_mangle[NGC_select]
RSD_SGC_flux = n_gal_mangle[SGC_select] - clfSGC.intercept_ - clfSGC.coef_[1] * psf_mangle[SGC_select] - clfSGC.coef_[0] * n_star_mangle[SGC_select]

x_RSD_star_N, RSD_star_N, RSD_star_N_std = percentile_average(n_star_mangle[NGC_select], RSD_NGC_star)
x_RSD_star_S, RSD_star_S, RSD_star_S_std = percentile_average(n_star_mangle[SGC_select], RSD_SGC_star)
x_RSD_psf_N, RSD_psf_N, RSD_psf_N_std = percentile_average(psf_mangle[NGC_select], RSD_NGC_psf)
x_RSD_psf_S, RSD_psf_S, RSD_psf_S_std = percentile_average(psf_mangle[SGC_select], RSD_SGC_psf)
x_RSD_flux_N, RSD_flux_N, RSD_flux_N_std = percentile_average(flux_mangle[NGC_select], RSD_NGC_flux)
x_RSD_flux_S, RSD_flux_S, RSD_flux_S_std = percentile_average(flux_mangle[SGC_select], RSD_SGC_flux)
'''

def ppplot(x_RSD_star_N, RSD_star_N, RSD_star_N_std,n_star_mangle,NGC_select,a,i,clfNGC):
	plt.plot(x_RSD_star_N[1:],RSD_star_N[1:],'o',linewidth=2)
	plt.errorbar(x_RSD_star_N[1:],RSD_star_N[1:], xerr=0, yerr=RSD_star_N_std[1:],linestyle='None',color='b')
	plt.plot(a,clfNGC.coef_[i]*a,'r')
	plt.hist(n_star_mangle[NGC_select],bins=100,histtype='step',color='g',weights=0.005*np.ones(sum(NGC_select)))
	plt.axis([min(a),max(a),-10,15])
	#plt.xlabel("n_star")
	plt.ylabel("Residual Surface Density")

an=[]
bn=[]
cn=[]
sa=[]
sb=[]
sc=[]

for i in range(6):
	a,b,c = percentile_average(X_NGC[:,i],(n_gal_mangle[NGC_select] - PSD_NGC + clfNGC.coef_[i] * X_NGC[:,i]))
	an.append(a)
	bn.append(b)
	cn.append(c)
	a,b,c = percentile_average(X_SGC[:,i],(n_gal_mangle[SGC_select] - PSD_SGC + clfSGC.coef_[i] * X_SGC[:,i]))
	sa.append(a)
	sb.append(b)
	sc.append(c)

#########################       plot 1      ############################
plt.figure(figsize=(20,10))
plt.subplot(321)
plt.title("NGC")
ppplot(an[0],bn[0],cn[0],n_star_mangle,NGC_select,np.arange(0,10000,1000),0,clfNGC)
plt.xlabel("stellar density")

plt.subplot(322)
plt.title("SGC")
ppplot(sa[0],sb[0],sc[0],n_star_mangle,SGC_select,np.arange(0,10000,1000),0,clfSGC)
plt.xlabel("stellar density ")

plt.subplot(323)
ppplot(an[1],bn[1],cn[1],psf_mangle,NGC_select,np.arange(0.8,2.2,0.2),1,clfNGC)
plt.xlabel("i-band PSF_FWHM")

plt.subplot(324)
ppplot(sa[1],sb[1],sc[1],psf_mangle,SGC_select,np.arange(0.8,2.2,0.2),1,clfSGC)
plt.xlabel("i-band PSF_FWHM")

plt.subplot(325)
ppplot(an[2],bn[2],cn[2],flux_mangle,NGC_select,np.arange(10,60,10),2,clfNGC)
plt.xlabel("z-band SKYFLUX")

plt.subplot(326)
ppplot(sa[2],sb[2],sc[2],flux_mangle,SGC_select,np.arange(10,60,10),2,clfSGC)
plt.xlabel("z-band SKYFLUX")

plt.savefig('compare_prakash_6_quantities_1.png')

#########################       plot 2      ############################
plt.figure(figsize=(20,10))
plt.subplot(321)
plt.title("NGC")
ppplot(an[3],bn[3],cn[3],extinction,NGC_select,np.arange(0,0.3,0.03),3,clfNGC)
plt.xlabel("EB_MINUS_V")

plt.subplot(322)
plt.title("SGC")
ppplot(sa[3],sb[3],sc[3],extinction,SGC_select,np.arange(0,0.3,0.03),3,clfSGC)
plt.xlabel("EB_MINUS_V")

plt.subplot(323)
aa = n_gal_mangle[NGC_select] - PSD_NGC + clfNGC.coef_[4] * X_NGC[:,4]
mean_EB = []
mean_std = []
for i in range(5):
	selection = X_NGC[:,4] == i
	mean_EB.append(np.mean(aa[selection]))
	mean_std.append(np.std(aa[selection])/len(aa[selection]))

plt.plot(range(5),mean_EB,'o',linewidth=2)
plt.errorbar(range(5),mean_EB, xerr=0, yerr=mean_std,linestyle='None',color='b')
plt.plot(np.arange(5),clfNGC.coef_[4]*np.arange(5),'r')
plt.hist(W1MOON[NGC_select],bins=100,histtype='step',color='g',weights=0.001*np.ones(sum(NGC_select)))
plt.axis([-0.1,5,-10,15])
plt.xlabel("W1MOON")
plt.ylabel("Residual Surface Density")

plt.subplot(324)
aa = n_gal_mangle[SGC_select] - PSD_SGC + clfSGC.coef_[4] * X_SGC[:,4]
mean_EB = []
mean_std = []
for i in range(5):
	selection = X_SGC[:,4] == i
	mean_EB.append(np.mean(aa[selection]))
	mean_std.append(np.std(aa[selection])/len(aa[selection]))

plt.plot(range(5),mean_EB,'o',linewidth=2)
plt.errorbar(range(5),mean_EB, xerr=0, yerr=mean_std,linestyle='None',color='b')
plt.plot(np.arange(5),clfSGC.coef_[4]*np.arange(5),'r')
plt.hist(W1MOON[SGC_select],bins=100,histtype='step',color='g',weights=0.001*np.ones(sum(SGC_select)))
plt.axis([-0.1,5,-10,15])
plt.xlabel("W1MOON")
plt.ylabel("Residual Surface Density")

plt.subplot(325)
ppplot(an[5],bn[5],cn[5],W1cov,NGC_select,np.arange(15,55,5),5,clfNGC)
plt.xlabel("W1cov")

plt.subplot(326)
ppplot(sa[5],sb[5],sc[5],W1cov,SGC_select,np.arange(15,55,5),5,clfSGC)
plt.xlabel("W1cov")

plt.savefig('compare_prakash_6_quantities_2.png')

plt.show()



