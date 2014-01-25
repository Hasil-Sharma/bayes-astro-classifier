from time import gmtime, strftime, clock
start = clock()
import mixfit as mf
from scipy.stats import norm
from random import seed, sample
from numpy import loadtxt
from scipy.optimize import minimize
import numpy as np
import scipy.stats as st
import csv
from matplotlib import pyplot as plt
import histo as ht
# Model is fit in 4 steps:
# 1. Fit Gaussian mixture model for ClassStat noise distribution by ML
# 2. Use UKIDSS pipeline posterior probability to roughly identify galaxy locus (idxgal = prob_star < 0.1) and fit galaxy locus parameters by LS
# 3. Use UKIDSS pipeline posterior probability to get number counts (idxgal = prob_star < 0.5) and fit number count model by LS
# 4. Fit magnitude distribution by ML

date = strftime("%a %d %b %Y %H:%M:%S", gmtime())
temp_seed = date.split(" ")
temp_seed2 = temp_seed[4].split(":")
temp_seed2 = np.sum([float(temp_seed2[0]),float(temp_seed2[1]),float(temp_seed2[2])])
temp_seed = np.sum([float(temp_seed[1]),float(temp_seed[3])])
temp_seed = temp_seed*temp_seed2
seed(temp_seed)
del(temp_seed,temp_seed2)

f = open("C:\Users\mUlm\Astronomy\Henrion\sample_with-colours_no-miss_std.csv","r")
lasdat = loadtxt(f, delimiter=',',skiprows=1)

# how do I associate a column name with each variable?
n = np.shape(lasdat)[0]
rows = np.arange(0,n,dtype=int)

idx_miss_y = rows[lasdat[:,31]==1]
idx_miss_j = rows[lasdat[:,32]==1]
idx_miss_j = rows[lasdat[:,33]==1]
idx_miss_k = rows[lasdat[:,34]==1]

# Taking a random sample
idx_rv = sample(rows,20000)
lasdatrv = lasdat[idx_rv,:]

#######################################################
# Mixture model for star locus and galaxy locus noise #
#######################################################

starpars = {'y':0.0,'j':0.0,'h':0.0,'k':0.0} # This must be a list

# Not sure where this is used
def dstar(x,starpars):  
	res = (starpars[0]*norm.pdf(x+starpars[1],scale=1) + (1-starpars[0])*norm.pdf(x,loc=starpars[2],scale=starpars[3]))
	return res

## Yband

dtemp = lasdat[((lasdat[:,5] > 13) & (lasdat[:,5] < 17) & (abs(lasdat[:,17]) < 10) & (lasdat[:,31] == 0) & (lasdat[:,17] != 0)),:]

def starloglik(para):
	res = para[0]*norm.pdf(dtemp[:,17]+para[1],scale=1) + (1-para[0])*norm.pdf(dtemp[:,17],loc=para[2],scale=para[3])
	res[res==0] = 10**(-200)
	res = -np.sum(np.log(res))
	return res

init = [0.95,0.1,4,4]
mybounds = [(0,1),(-50,50),(-50,50),(0.01,50)]
optpars = minimize(starloglik,x0=init,bounds=mybounds,method='L-BFGS-B') # gets slightly different results than Henrion's R code
starpars['y'] = optpars 
del(dtemp,optpars)

## Jband

dtemp = lasdat[((lasdat[:,7] > 13) & (lasdat[:,7] < 15) & (abs(lasdat[:,18]) < 12) & (lasdat[:,32] == 0) & (lasdat[:,18] != 0)),:]

def starloglik(para):
	res = para[0]*norm.pdf(dtemp[:,18]+para[1],scale=1) + (1-para[0])*norm.pdf(dtemp[:,18],loc=para[2],scale=para[3])
	res[res==0] = 10**(-200)
	res = -np.sum(np.log(res))
	return res

init = [0.95,0.1,4,4]
mybounds = [(0,1),(-50,50),(-50,50),(0.01,50)]
optpars = minimize(starloglik,x0=init,bounds=mybounds,method='L-BFGS-B') # gets slightly different results than Henrion's R code
starpars['j'] = optpars
del(dtemp,optpars)

## Hband

dtemp = lasdat[((lasdat[:,9] > 12.5) & (lasdat[:,9] < 15) & (abs(lasdat[:,19]) < 12) & (lasdat[:,33] == 0) & (lasdat[:,19] != 0)),:]

def starloglik(para):
	res = para[0]*norm.pdf(dtemp[:,19]+para[1],scale=1) + (1-para[0])*norm.pdf(dtemp[:,19],loc=para[2],scale=para[3])
	res[res==0] = 10**(-200)
	res = -np.sum(np.log(res))
	return res

init = [0.95,0.1,4,4]
mybounds = [(0,1),(-50,50),(-50,50),(0.01,50)]
optpars = minimize(starloglik,x0=init,bounds=mybounds,method='L-BFGS-B') # gets slightly different results than Henrion's R code
starpars['h'] = optpars
del(dtemp,optpars)

vec = lasdat[:,19][(lasdat[:,9] > 12.5) & (lasdat[:,9] < 15) & (abs(lasdat[:,19]) < 35) & (lasdat[:,33] == 0) & (lasdat[:,19] != 0)]
ht.my_hist(vec,2)

## Kband

dtemp = lasdat[((lasdat[:,11] > 12.5) & (lasdat[:,11] < 15) & (abs(lasdat[:,20]) < 12) & (lasdat[:,34] == 0) & (lasdat[:,20] != 0)),:]

def starloglik(para):
	res = para[0]*norm.pdf(dtemp[:,20]+para[1],scale=1) + (1-para[0])*norm.pdf(dtemp[:,20],loc=para[2],scale=para[3])
	res[res==0] = 10**(-200)
	res = -np.sum(np.log(res))
	return res

init = [0.95,0.1,4,4]
mybounds = [(0,1),(-50,50),(-50,50),(0.01,50)]
optpars = minimize(starloglik,x0=init,bounds=mybounds,method='L-BFGS-B') # same results as Henrion's R code
starpars['k'] = optpars
del(dtemp,optpars)

########################
# Estimate Galaxy Locus#
########################

n = np.shape(lasdatrv)[0]
rows = np.arange(0,n,dtype=int)

idx_gal_rv_y = rows[(lasdatrv[:,31]==0) & (lasdatrv[:,27]<0.1)] # This is of length 9906 and should be 9910 and I think the problem with MedPars lies here
idx_gal_rv_j = rows[(lasdatrv[:,32]==0) & (lasdatrv[:,27]<0.1)]
idx_gal_rv_h = rows[(lasdatrv[:,33]==0) & (lasdatrv[:,27]<0.1)]
idx_gal_rv_k = rows[(lasdatrv[:,34]==0) & (lasdatrv[:,27]<0.1)]

SmpSize = 50
Nbins_y = np.trunc(len(idx_gal_rv_y)/SmpSize)
Nbins_j = np.trunc(len(idx_gal_rv_j)/SmpSize)
Nbins_h = np.trunc(len(idx_gal_rv_h)/SmpSize)
Nbins_k = np.trunc(len(idx_gal_rv_k)/SmpSize)

rnk_y = st.mstats.rankdata(lasdatrv[idx_gal_rv_y,5]) # no random tie method for scipy.stats
rnk_j = st.mstats.rankdata(lasdatrv[idx_gal_rv_j,7])
rnk_h = st.mstats.rankdata(lasdatrv[idx_gal_rv_h,9])
rnk_k = st.mstats.rankdata(lasdatrv[idx_gal_rv_k,11])

meanmag = {'y':0,'j':0,'h':0,'k':0}
var_cs = {'y':0,'j':0,'h':0,'k':0}
var_mag = {'y':0,'j':0,'h':0,'k':0}
med = {'y':0,'j':0,'h':0,'k':0}

meanmag['y'] = np.repeat(0.0,Nbins_y+1); meanmag['j'] = np.repeat(0.0,Nbins_j+1); meanmag['h'] = np.repeat(0.0,Nbins_h+1); meanmag['k'] = np.repeat(0.0,Nbins_k+1)
var_cs['y'] = np.repeat(0.0,Nbins_y+1); var_cs['j'] = np.repeat(0.0,Nbins_j+1); var_cs['h'] = np.repeat(0.0,Nbins_h+1); var_cs['k'] = np.repeat(0.0,Nbins_k+1)
var_mag['y'] = np.repeat(0.0,Nbins_y+1); var_mag['j'] = np.repeat(0.0,Nbins_j+1); var_mag['h'] = np.repeat(0.0,Nbins_h+1); var_mag['k'] = np.repeat(0.0,Nbins_k+1)
med['y'] = np.repeat(0.0,Nbins_y+1); med['j'] = np.repeat(0.0,Nbins_j+1); med['h'] = np.repeat(0.0,Nbins_h+1); med['k'] = np.repeat(0.0,Nbins_k+1)

for band in xrange(0,4):
	if band == 0:
		n = np.arange(0,len(idx_gal_rv_y),dtype=int)
		dtemp = lasdatrv[idx_gal_rv_y,:]
		for i in xrange(0,int(Nbins_y)):
			if i==0:
				magdata_idx = n[rnk_y <= SmpSize*(i+1)]
			if 0 < i < (Nbins_y - 1):
				magdata_idx = n[(rnk_y > SmpSize*(i)) & (rnk_y <= SmpSize*(i+1))]
			if i == (Nbins_y -1):
				magdata_idx = n[rnk_y > SmpSize*(i)]
			magvec = dtemp[:,5]
			csvec = dtemp[:,17]
			meanmag['y'][i] = np.mean(magvec[magdata_idx])
			var_cs['y'][i] = np.var(csvec[magdata_idx],ddof=1)
			var_mag['y'][i] = np.var(magvec[magdata_idx],ddof=1)
			med['y'][i] = np.mean(csvec[magdata_idx])
		ntotal = np.arange(0,np.shape(lasdat)[0],dtype=int)
		
		idx_max = ntotal[(lasdat[:,31]==0) & (lasdat[:,5]==np.max(lasdat[:,5]))]
		idx_sim = ntotal[lasdat[:,5]==max(lasdat[:,5])]
		minus1 = np.setdiff1d(ntotal,idx_sim)
		idx_max = np.concatenate((idx_max,ntotal[(lasdat[:,31] == 0) & (lasdat[:,5] == np.max(lasdat[:,5][minus1]))]))  
		conc1 = np.concatenate((idx_sim,idx_max))
		idx_sim = np.unique(conc1) 
		minus2 = np.setdiff1d(ntotal,idx_sim)
		idx_max = np.concatenate((idx_max,ntotal[(lasdat[:,31] == 0) & (lasdat[:,5] == np.max(lasdat[:,5][minus2]))]))
		dtemp = lasdat[idx_max,:]
		meanmag['y'][Nbins_y] = np.mean(dtemp[:,5])
		var_cs['y'][Nbins_y] = np.var(dtemp[:,17],ddof=1)
		var_mag['y'][Nbins_y] = np.var(dtemp[:,5],ddof=1)
		med['y'][Nbins_y] = np.mean(dtemp[:,17])
	if band == 1:
		n = np.arange(0,len(idx_gal_rv_j),dtype=int)
		dtemp = lasdatrv[idx_gal_rv_j,:]
		for i in xrange(0,int(Nbins_j)):
			if i==0:
				magdata_idx = n[rnk_j <= SmpSize*(i+1)]
			if 0 < i < (Nbins_j - 1):
				magdata_idx = n[(rnk_j > SmpSize*(i)) & (rnk_j <= SmpSize*(i+1))]
			if i == (Nbins_j -1):
				magdata_idx = n[rnk_j > SmpSize*(i)]
			magvec = dtemp[:,7]
			csvec = dtemp[:,18]
			meanmag['j'][i] = np.mean(magvec[magdata_idx])
			var_cs['j'][i] = np.var(csvec[magdata_idx],ddof=1)
			var_mag['j'][i] = np.var(magvec[magdata_idx],ddof=1)
			med['j'][i] = np.mean(csvec[magdata_idx])
		ntotal = np.arange(0,np.shape(lasdat)[0],dtype=int)
		
		idx_max = ntotal[(lasdat[:,32]==0) & (lasdat[:,7]==np.max(lasdat[:,7]))]
		idx_sim = ntotal[lasdat[:,7]==max(lasdat[:,7])]
		minus1 = np.setdiff1d(ntotal,idx_sim)
		idx_max = np.concatenate((idx_max,ntotal[(lasdat[:,32] == 0) & (lasdat[:,7] == np.max(lasdat[:,7][minus1]))]))  
		conc1 = np.concatenate((idx_sim,idx_max))
		idx_sim = np.unique(conc1) 
		minus2 = np.setdiff1d(ntotal,idx_sim)
		idx_max = np.concatenate((idx_max,ntotal[(lasdat[:,32] == 0) & (lasdat[:,7] == np.max(lasdat[:,7][minus2]))]))
		dtemp = lasdat[idx_max,:]
		meanmag['j'][Nbins_j] = np.mean(dtemp[:,7])
		var_cs['j'][Nbins_j] = np.var(dtemp[:,18],ddof=1)
		var_mag['j'][Nbins_j] = np.var(dtemp[:,7],ddof=1)
		med['j'][Nbins_j] = np.mean(dtemp[:,18])				
	if band == 2:
		n = np.arange(0,len(idx_gal_rv_h),dtype=int)
		dtemp = lasdatrv[idx_gal_rv_h,:]
		for i in xrange(0,int(Nbins_h)):
			if i==0:
				magdata_idx = n[rnk_h <= SmpSize*(i+1)]
			if 0 < i < (Nbins_h - 1):
				magdata_idx = n[(rnk_h > SmpSize*(i)) & (rnk_h <= SmpSize*(i+1))]
			if i == (Nbins_h -1):
				magdata_idx = n[rnk_h > SmpSize*(i)]
			magvec = dtemp[:,9]
			csvec = dtemp[:,19]
			meanmag['h'][i] = np.mean(magvec[magdata_idx])
			var_cs['h'][i] = np.var(csvec[magdata_idx],ddof=1)
			var_mag['h'][i] = np.var(magvec[magdata_idx],ddof=1)
			med['h'][i] = np.mean(csvec[magdata_idx])
		ntotal = np.arange(0,np.shape(lasdat)[0],dtype=int)
		
		idx_max = ntotal[(lasdat[:,33]==0) & (lasdat[:,9]==np.max(lasdat[:,9]))]
		idx_sim = ntotal[lasdat[:,9]==max(lasdat[:,9])]
		minus1 = np.setdiff1d(ntotal,idx_sim)
		idx_max = np.concatenate((idx_max,ntotal[(lasdat[:,33] == 0) & (lasdat[:,9] == np.max(lasdat[:,9][minus1]))]))  
		conc1 = np.concatenate((idx_sim,idx_max))
		idx_sim = np.unique(conc1) 
		minus2 = np.setdiff1d(ntotal,idx_sim)
		idx_max = np.concatenate((idx_max,ntotal[(lasdat[:,33] == 0) & (lasdat[:,9] == np.max(lasdat[:,9][minus2]))]))
		dtemp = lasdat[idx_max,:]
		meanmag['h'][Nbins_h] = np.mean(dtemp[:,9])
		var_cs['h'][Nbins_h] = np.var(dtemp[:,19],ddof=1)
		var_mag['h'][Nbins_h] = np.var(dtemp[:,9],ddof=1)
		med['h'][Nbins_h] = np.mean(dtemp[:,19])					
	if band == 3:
		n = np.arange(0,len(idx_gal_rv_k),dtype=int)
		dtemp = lasdatrv[idx_gal_rv_k,:]
		for i in xrange(0,int(Nbins_k)):
			if i==0:
				magdata_idx = n[rnk_k <= SmpSize*(i+1)]
			if 0 < i < (Nbins_k - 1):
				magdata_idx = n[(rnk_k > SmpSize*(i)) & (rnk_k <= SmpSize*(i+1))]
			if i == (Nbins_k -1):
				magdata_idx = n[rnk_k > SmpSize*(i)]
			magvec = dtemp[:,11]
			csvec = dtemp[:,20]
			meanmag['k'][i] = np.mean(magvec[magdata_idx])
			var_cs['k'][i] = np.var(csvec[magdata_idx],ddof=1)
			var_mag['k'][i] = np.var(magvec[magdata_idx],ddof=1)
			med['k'][i] = np.mean(csvec[magdata_idx])
		ntotal = np.arange(0,np.shape(lasdat)[0],dtype=int)
		
		idx_max = ntotal[(lasdat[:,34]==0) & (lasdat[:,11]==np.max(lasdat[:,11]))]
		idx_sim = ntotal[lasdat[:,11]==max(lasdat[:,11])]
		minus1 = np.setdiff1d(ntotal,idx_sim)
		idx_max = np.concatenate((idx_max,ntotal[(lasdat[:,34] == 0) & (lasdat[:,11] == np.max(lasdat[:,11][minus1]))]))  
		conc1 = np.concatenate((idx_sim,idx_max))
		idx_sim = np.unique(conc1) 
		minus2 = np.setdiff1d(ntotal,idx_sim)
		idx_max = np.concatenate((idx_max,ntotal[(lasdat[:,34] == 0) & (lasdat[:,11] == np.max(lasdat[:,11][minus2]))]))
		dtemp = lasdat[idx_max,:]
		meanmag['k'][Nbins_k] = np.mean(dtemp[:,11])
		var_cs['k'][Nbins_k] = np.var(dtemp[:,20],ddof=1)
		var_mag['k'][Nbins_k] = np.var(dtemp[:,11],ddof=1)
		med['k'][Nbins_k] = np.mean(dtemp[:,20])

def SSmedCSHyp(para):
	f_idx = np.round(Nbins_y/10)
	check = 0
	para[2] = 142
	if check==0:
		mmm = np.max([(np.max(lasdat[:,5][lasdat[:,31]==0])),(np.max(lasdat[:,7][lasdat[:,32]==0])+para[5]-para[6]),(np.max(lasdat[:,9][lasdat[:,33]==0])+para[5]-para[7]),(np.max(lasdat[:,11][lasdat[:,34]==0])+para[5]-para[8])])
		res1 = (med['y'][f_idx:(Nbins_y+1)] - ( (1-(meanmag['y'][f_idx:(Nbins_y+1)]+(para[5]-para[5]))/mmm)*( (para[0]*(meanmag['y'][f_idx:(Nbins_y+1)]-para[5])**2 + para[1]*(meanmag['y'][f_idx:(Nbins_y+1)]-para[5]) + para[2])**(para[4]) + para[3]) ) )**2
		res2 = (med['j'][f_idx:(Nbins_j+1)] - ( (1-(meanmag['j'][f_idx:(Nbins_j+1)]+(para[5]-para[6]))/mmm)*( (para[0]*(meanmag['j'][f_idx:(Nbins_j+1)]-para[6])**2 + para[1]*(meanmag['j'][f_idx:(Nbins_j+1)]-para[6]) + para[2])**(para[4]) + para[3]) ) )**2
		res3 = (med['h'][f_idx:(Nbins_h+1)] - ( (1-(meanmag['h'][f_idx:(Nbins_h+1)]+(para[5]-para[7]))/mmm)*( (para[0]*(meanmag['h'][f_idx:(Nbins_h+1)]-para[7])**2 + para[1]*(meanmag['h'][f_idx:(Nbins_h+1)]-para[7]) + para[2])**(para[4]) + para[3]) ) )**2
		res4 = (med['k'][f_idx:(Nbins_k+1)] - ( (1-(meanmag['k'][f_idx:(Nbins_k+1)]+(para[5]-para[8]))/mmm)*( (para[0]*(meanmag['k'][f_idx:(Nbins_k+1)]-para[8])**2 + para[1]*(meanmag['k'][f_idx:(Nbins_k+1)]-para[8]) + para[2])**(para[4]) + para[3]) ) )**2
		res = np.sum(res1) + np.sum(res2) + np.sum(res3) + np.sum(res4)
	return res

def gr_SSmedCSHyp(para):
	res = SSmedCSHyp(para)
	res = np.gradient(res)
	return res

# the function in SSmedCSHyp is an empirical function simply chosen to give a good fit; this hunction can change depending on the data
init = [10,8,142,-56,0.8,23,22,21,20]
MedPars = minimize(SSmedCSHyp,x0=init,method='BFGS',options={'maxiter':250}) #would like to use the gradient here
print MedPars

mmagshift = {'y':0.0,'j':0.0,'h':0.0,'k':0.0}
mmagshift['y'] = MedPars.x[5] - MedPars.x[5]
mmagshift['j'] = MedPars.x[5] - MedPars.x[6]
mmagshift['h'] = MedPars.x[5] - MedPars.x[7]
mmagshift['k'] = MedPars.x[5] - MedPars.x[8]

pa = MedPars.x # parameters are not right
mmm = np.max([(np.max(lasdat[:,5][lasdat[:,31]==0])),(np.max(lasdat[:,7][lasdat[:,32]==0])+pa[5]-pa[6]),(np.max(lasdat[:,9][lasdat[:,33]==0])+pa[5]-pa[7]),(np.max(lasdat[:,11][lasdat[:,34]==0])+pa[5]-pa[8])])
step = (23-11)/249.0
m = np.arange(11,(23+step),step)
yval1 = (1-(m+pa[5]-pa[5])/mmm)*((pa[0]*(m-pa[5])**2+pa[1]*(m-pa[5])+pa[2])**(pa[4])+pa[3])
yval2 = (1-(m+pa[5]-pa[6])/mmm)*((pa[0]*(m-pa[6])**2+pa[1]*(m-pa[6])+pa[2])**(pa[4])+pa[3])
yval3 = (1-(m+pa[5]-pa[7])/mmm)*((pa[0]*(m-pa[7])**2+pa[1]*(m-pa[7])+pa[2])**(pa[4])+pa[3])
yval4 = (1-(m+pa[5]-pa[8])/mmm)*((pa[0]*(m-pa[8])**2+pa[1]*(m-pa[8])+pa[2])**(pa[4])+pa[3])

def meanhist():
	plt.scatter(meanmag['y'], med['y'], c='k')
	plt.scatter(meanmag['j'], med['j'], c='r')
	plt.scatter(meanmag['h'], med['h'], c='m')
	plt.scatter(meanmag['k'], med['k'], c='y')
	plt.plot(m,yval1)
	plt.plot(m,yval2)
	plt.plot(m,yval3)
	plt.plot(m,yval4)

	plt.xlabel("Mean Magnitude")
	plt.xlim(13,23) # X axis limits

	plt.ylabel("Median Magnitude")
	plt.ylim(-2,25) # Y axis limits
	plt.title("Distribution of mu'(m)")
	plt.show() # Uncomment to show plot
	#plt.savefig("meanhist.png".format(file))
	
meanhist()

def SSvariCS(para):
	f_idx = round(Nbins_y/70)-1
	check = 0
	if para[1] < (-1):
		res = 10**(12); check = 1
	if check == 0:
		res1 = (var_cs['y'][f_idx:(Nbins_y-1)] -1 - para[0]*(10**5)*((10)**(para[1]*(meanmag['y'][f_idx:(Nbins_y-1)] + mmagshift['y'] - 11))) )**2
		res2 = (var_cs['j'][f_idx:(Nbins_j-1)] -1 - para[0]*(10**5)*((10)**(para[1]*(meanmag['j'][f_idx:(Nbins_j-1)] + mmagshift['j'] - 11))) )**2
		res3 = (var_cs['h'][f_idx:(Nbins_h-1)] -1 - para[0]*(10**5)*((10)**(para[1]*(meanmag['h'][f_idx:(Nbins_h-1)] + mmagshift['h'] - 11))) )**2
		res4 = (var_cs['k'][f_idx:(Nbins_k-1)] -1 - para[0]*(10**5)*((10)**(para[1]*(meanmag['k'][f_idx:(Nbins_k-1)] + mmagshift['k'] - 11))) )**2    
		res = np.sum(res1) + np.sum(res2) + np.sum(res3) + np.sum(res4)
	return res
	
init = [1.5,-0.5]	
VarPars = minimize(SSvariCS,x0=init,method='BFGS')

## NB the model pa[1]*(10^5)*10(pa[2]*(m-11)) is equivalent to the model 10^(pa[2]*(m-pa'[1]))
## where pa'[1]=11-(5+log10(pa[1]))/pa[2]
## we have given the simpler model (with pa[1]) in the thesis...
## NBB whether we specify that sigma=10^pa[2]*(m-pa[1]) or sigma^2=10^pa'[2]*(m-pa[1]) is unimportant: pa'[2]=2*pa[2]

plt.scatter(meanmag['y'] + mmagshift['y'],var_cs['y'] - 1, c='k')
plt.scatter(meanmag['j'] + mmagshift['j'],var_cs['j'] - 1, c='c')
plt.scatter(meanmag['h'] + mmagshift['h'],var_cs['h'] - 1, c='r')
plt.scatter(meanmag['k'] + mmagshift['k'],var_cs['k'] - 1, c='y')
	
def varhist():
	step = (23-14)/249.0
	m = np.arange(14,(23+step),step)
	vcs = VarPars.x[0]*(10**5)*((10)**(VarPars.x[1]*(m-11)))
	plt.plot(m,vcs,c='b')

	plt.xlabel("Mean + Shifted Magnitude")
	plt.xlim(14,23) # X axis limits

	plt.ylabel("CS Variance - 1")
	plt.ylim(-2,200) # Y axis limits
	plt.title("Distribution of sigma'(m)")
	plt.show() # Uncomment to show plot
	#plt.savefig("varhist.png".format(file))
	
varhist()


###############################
#  Fit Magnitude Distribution #
###############################

dtemp = np.concatenate((lasdatrv[lasdatrv[:,31]==0,:][:,5]+mmagshift['y'],lasdatrv[lasdatrv[:,32]==0][:,7]+mmagshift['j'],lasdatrv[lasdatrv[:,33]==0][:,9]+mmagshift['h'],lasdatrv[lasdatrv[:,34]==0][:,11]+mmagshift['k']))
dis = lasdatrv[lasdatrv[:,31]==0,5]
print mmagshift['y'], 'y'
print mmagshift['j'], 'j'
print mmagshift['h'], 'h'
print mmagshift['k'], 'k'

# length of dtemp is not quite right
init = [20,0.5,0.25,12]
#MagParsGl = mf.FitMagInit(dtemp,init) 
#print MagParsGl

#####################
#   Number Counts   #
#####################

# NB as we want to get parameters that hold true simultaneously for all 4 bands, we use stars and galaxies for all bands.
# This means that any object can be counted up to 4 times during the estimation procedure.
# However, this is more a weighting problem than anything else: as we normalize to get the final a(m), the number of objects (and their sum of weights) is nor relevant.

n = np.shape(lasdatrv)[0]
ntotal = np.arange(0,n,1)
idx_gal_y = ntotal[(lasdatrv[:,31]==0) & (lasdatrv[:,27]<0.5)]
idx_gal_j = ntotal[(lasdatrv[:,32]==0) & (lasdatrv[:,27]<0.5)]
idx_gal_h = ntotal[(lasdatrv[:,33]==0) & (lasdatrv[:,27]<0.5)]
idx_gal_k = ntotal[(lasdatrv[:,34]==0) & (lasdatrv[:,27]<0.5)]

step = (23-12)/1999.0
m = np.arange(12,(23+step),step)

dtemp = {'st':0.0,'ga':0.0}
minusy = np.setdiff1d(ntotal,idx_gal_y)
minusj = np.setdiff1d(ntotal,idx_gal_j)
minush = np.setdiff1d(ntotal,idx_gal_h)
minusk = np.setdiff1d(ntotal,idx_gal_k)

sty = lasdatrv[:,5][minusy]; sty = sty[(lasdatrv[:,31][minusy]==0) & (lasdatrv[:,5][minusy]+mmagshift['y'] > 15) & (lasdatrv[:,5][minusy]+mmagshift['y'] < 21)]
stj = lasdatrv[:,7][minusj]; stj = stj[(lasdatrv[:,32][minusj]==0) & (lasdatrv[:,7][minusj]+mmagshift['j'] > 15) & (lasdatrv[:,7][minusj]+mmagshift['j'] < 21)]
sth = lasdatrv[:,9][minush]; sth = sth[(lasdatrv[:,33][minush]==0) & (lasdatrv[:,9][minush]+mmagshift['h'] > 15) & (lasdatrv[:,9][minush]+mmagshift['h'] < 21)]
stk = lasdatrv[:,11][minusk]; stk = stk[(lasdatrv[:,34][minusk]==0) & (lasdatrv[:,11][minusk]+mmagshift['k'] > 15) & (lasdatrv[:,11][minusk]+mmagshift['k'] < 21)]
dtemp['st'] = np.concatenate((sty + mmagshift['y'],stj + mmagshift['j'],sth + mmagshift['h'],stk + mmagshift['k']))

gay = lasdatrv[:,5][idx_gal_y]; gay = gay[(lasdatrv[:,5][idx_gal_y] + mmagshift['y'] > 15) & (lasdatrv[:,5][idx_gal_y] + mmagshift['y'] < 21)]
gaj = lasdatrv[:,7][idx_gal_j]; gaj = gaj[(lasdatrv[:,7][idx_gal_j] + mmagshift['j'] > 15) & (lasdatrv[:,7][idx_gal_j] + mmagshift['j'] < 21)]
gah = lasdatrv[:,9][idx_gal_h]; gah = gah[(lasdatrv[:,9][idx_gal_h] + mmagshift['h'] > 15) & (lasdatrv[:,9][idx_gal_h] + mmagshift['h'] < 21)]
gak = lasdatrv[:,11][idx_gal_k]; gak = gak[(lasdatrv[:,11][idx_gal_k] + mmagshift['k'] > 15) & (lasdatrv[:,11][idx_gal_k] + mmagshift['k'] < 21)]
dtemp['ga'] = np.concatenate((gay + mmagshift['y'],gaj + mmagshift['j'],gah + mmagshift['h'],gak + mmagshift['k']))

print np.shape(sty)[0]+np.shape(gay)[0]; print np.shape(lasdatrv[(lasdatrv[:,31]==0) & (lasdatrv[:,5]+mmagshift['y'] > 15) & (lasdatrv[:,5]+mmagshift['y'] < 21)])[0]
print np.shape(stj)[0]+np.shape(gaj)[0]; print np.shape(lasdatrv[(lasdatrv[:,32]==0) & (lasdatrv[:,7]+mmagshift['j'] > 15) & (lasdatrv[:,7]+mmagshift['j'] < 21)])[0]
print np.shape(sth)[0]+np.shape(gah)[0]; print np.shape(lasdatrv[(lasdatrv[:,33]==0) & (lasdatrv[:,9]+mmagshift['h'] > 15) & (lasdatrv[:,9]+mmagshift['h'] < 21)])[0]
print np.shape(stk)[0]+np.shape(gak)[0]; print np.shape(lasdatrv[(lasdatrv[:,34]==0) & (lasdatrv[:,11]+mmagshift['k'] > 15) & (lasdatrv[:,11]+mmagshift['k'] < 21)])[0]
del(sty,stj,sth,stk,gay,gaj,gah,gak)


step = (20.99-15.01)/299
mids = np.arange(15.01,(20.99),step)
conc = np.concatenate((dtemp['st'],dtemp['ga']))
hall = plt.hist(conc, bins = 299, normed = False)
hst = plt.hist(dtemp['st'], bins = 299,normed = False)
hga = plt.hist(dtemp['ga'], bins = 299, normed = False)

plt.show()

yval1 = np.log(hall[0])
plt.plot(mids, yval1, color = 'k', label = 'Number Counts (Log Scale)')
plt.plot(mids, np.log(hst[0]), color = 'b')
plt.plot(mids, np.log(hga[0]), color = 'r')
plt.xlim(15,22)
plt.ylim(2,7)
plt.xlabel("Bins")
plt.ylabel("Log Counts")
plt.title("Number Counts (Log Scale)")
#plt.show()



def SSst(para):
	check = 0
	if para[0] < 0.0001:
		check = 1; res = 10**20
	if check == 0:
		m = mids[hst[0] > 0]
		y = np.log(hst[0][hst[0] > 0])
		res = np.sum((y - np.log(10**(para[0]*(m - para[1]))*0.5*mf.erfc((m-para[2])/para[3])))**2)
	return res
	
def SSga(para):
	check = 0
	if para[0] < 0.0001:
		check = 1; res = 10**20
	if check == 0:
		m = mids[hga[0] > 0]
		y = np.log(hga[0][hga[0] > 0])
		res = np.sum((y - np.log(10**(para[0]*(m - para[1]))*0.5*mf.erfc((m-para[2])/para[3])))**2)
	return res
	
stpars = minimize(SSst,[0.25,5,20,0.5],method='BFGS') 
gapars = minimize(SSga,[0.45,15,20,0.5],method='BFGS') 

lst = 10**(stpars.x[0]*(m-stpars.x[1])) * 0.5*mf.erfc((m-stpars.x[2])/stpars.x[3])
lga = 10**(gapars.x[0]*(m-gapars.x[1])) * 0.5*mf.erfc((m-gapars.x[2])/gapars.x[3])
plt.plot(m, np.log(lst), color = 'c')
plt.plot(m, np.log(lga), color = 'm')
plt.show()
#plt.savefig("numcounthist.png".format(file))

def MxCfGl(m):
    res = ( 10**(stpars.x[0]*(m-stpars.x[1]))*0.5*mf.erfc((m-stpars.x[2])/stpars.x[3]) ) / ( (10**(stpars.x[0]*(m-stpars.x[1])))*0.5*mf.erfc((m-stpars.x[2])/stpars.x[3]) + 10**(gapars.x[0]*(m-gapars.x[1]))*0.5*mf.erfc((m-gapars.x[2])/gapars.x[3]) )
    #res = ( 10**(stpars.x[0]*(m-stpars.x[1])) ) / ( (10**(stpars.x[0]*(m-stpars.x[1]))) + 10**(gapars.x[0]*(m-gapars.x[1])) )
    return res

plt.plot(m,MxCfGl(m), color = 'k')
plt.plot(mids,(hst[0]/(hst[0]+hga[0]+0.0)), color = 'r')
plt.show()
#plt.savefig("starperhist.png".format(file))
  
 ###################################
 #  Number Counts Plot for Thesis  #
 ###################################
 
 #############################
 #  Assessing the Model Fit  #
 #############################


print clock()-start, "seconds elapsed"
