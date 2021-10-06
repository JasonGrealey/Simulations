print("script beginning")
import pandas as pd
import numpy as np
from sys import argv
import time,random

import seaborn as sns
import argparse
import os
import sys
import statsmodels.api as sm
import csv
from sklearn.metrics import explained_variance_score
from sklearn.linear_model import LinearRegression

'''
input:option - string -  detailing if epistatic interactions should be modelled or not
input:loadstring - string - of what dataframe/file we are loading in for the genotypes
input:runno - integer - (not important for simulation - can be anything) but keeps track of saved files for numerous runs. i.e. If we simulated 2000 times, runno would go from 1...2000
input:sigmabeta - float - standard deviation of the distribution (normal usually) from which we pull the beta values
input:sigmagamma - float - standard deviation of the distribution (normal usually) from which we pull the gamma values which are our epistatic interactions
input:rescaleval - float - a number to which rescale the mean of the beta distribution (only in use when we are having our beta sampled distribution as a function of minor allele frequency)
'''

from numpy import count_nonzero
def t_model(dataframe,loci,locj,indexs):
		'''
		this function takes the dataframe and recalculates the interaction term for two epistatic SNPs at loci and locj in the dataframe
		pass dataframe (episnps)
		pass loci and locj randepis[loci] randepis[locj]
		pass indexs = samples i.e. only numeric calculations
		'''
		print("running t-model of epistatic interactions")
		assert isinstance(indexs,(list, tuple, np.ndarray))
		#print("performing multiplications")
		multiples=np.multiply(dataframe.loc[loci,indexs].values, dataframe.loc[locj,indexs].values)
		epistatic_samples=np.where(multiples >= 1.5,1,0)#(dataframe.loc['interaction'+str(loci)+'_'+str(locj)][indexs].values >= 1.5,1,0)
		epistatic_samples=np.reshape(epistatic_samples,-1)
		print(dataframe.loc[loci,indexs].values[np.where(epistatic_samples==1)])
		print(dataframe.loc[locj,indexs].values[np.where(epistatic_samples==1)])

		print(count_nonzero(epistatic_samples)/float(epistatic_samples.size))
		print("percentage non-zero above")
		vals=pd.Series(epistatic_samples, index=indexs,name=str('interaction_'+str(loci)+'_'+str(locj)))
		dataframe=dataframe.append(vals)
		print("penetrance classes with epistasis - should be in epi_table")
		print(np.unique(np.core.defchararray.add(dataframe.loc[loci,indexs].values[np.where(epistatic_samples==1)].astype(str),dataframe.loc[locj,indexs].values[np.where(epistatic_samples==1)].astype(str))))

		'''
		B = np.where(A > 0.5, 1, 0)#where A > 0.5 return 1 else 0
		for i in range(len(indexs)):
				if dataframe.loc['interaction'+str(loci)+'_'+str(locj)][indexs[i]]>=1.5:
						dataframe.loc['interaction'+str(loci)+'_'+str(locj)][indexs[i]]=1
		'''
		return dataframe
def RR_model(dataframe,loci,locj,indexs):
		'''
		this function takes the dataframe and recalculates the interaction term for two epistatic SNPs at loci and locj in the dataframe
		pass dataframe (episnps)
		pass loci and locj randepis[loci] randepis[locj]
		pass indexs = samples i.e. only numeric calculations
		'''
		print("running RR-model of epistatic interactions")
		assert isinstance(indexs,(list, tuple, np.ndarray))
		#print("performing multiplications")
		multiples=np.multiply(dataframe.loc[loci,indexs].values, dataframe.loc[locj,indexs].values)
		epistatic_samples=np.where(multiples >= 3.5,1,0)#(dataframe.loc['interaction'+str(loci)+'_'+str(locj)][indexs].values >= 1.5,1,0)
		epistatic_samples=np.reshape(epistatic_samples,-1)
		print(dataframe.loc[loci,indexs].values[np.where(epistatic_samples==1)])
		print(dataframe.loc[locj,indexs].values[np.where(epistatic_samples==1)])

		print(count_nonzero(epistatic_samples)/float(epistatic_samples.size))
		print("percentage non-zero above")
		print("penetrance classes with epistasis - should be in epi_table")
		print(np.unique(np.core.defchararray.add(dataframe.loc[loci,indexs].values[np.where(epistatic_samples==1)].astype(str),dataframe.loc[locj,indexs].values[np.where(epistatic_samples==1)].astype(str))))

		vals=pd.Series(epistatic_samples, index=indexs,name=str('interaction_'+str(loci)+'_'+str(locj)))
		dataframe=dataframe.append(vals)
		'''
		B = np.where(A > 0.5, 1, 0)#where A > 0.5 return 1 else 0
		for i in range(len(indexs)):
				if dataframe.loc['interaction'+str(loci)+'_'+str(locj)][indexs[i]]>=1.5:
						dataframe.loc['interaction'+str(loci)+'_'+str(locj)][indexs[i]]=1
		'''
		return dataframe
def m16_model(dataframe,loci,locj,indexs):
		'''
		this function takes the dataframe and recalculates the interaction term for two epistatic SNPs at loci and locj in the dataframe
		pass dataframe (episnps)
		pass loci and locj randepis[loci] randepis[locj]
		pass indexs = samples i.e. only numeric calculations
		double het model
		'''
		print("running m16-model of epistatic interactions")
		assert isinstance(indexs,(list, tuple, np.ndarray))
		#print("performing multiplications")
		multiples=np.multiply(dataframe.loc[loci,indexs].values, dataframe.loc[locj,indexs].values)
		#print(multiples)
		epistatic_samples=np.where(multiples == 1,1,0)#(dataframe.loc['interaction'+str(loci)+'_'+str(locj)][indexs].values >= 1.5,1,0)
		print(epistatic_samples)
		#sys.exit(0)
		#epistatic_samples=np.where(multiples > 1.5, 0, 
		# (np.where(consumption_energy < 0.5, 0, 1)))
		#print(dataframe.loc[loci,indexs].values[np.where(epistatic_samples==1)])
		#print(dataframe.loc[locj,indexs].values[np.where(epistatic_samples==1)])
		#sys.exit(0)
		#epistatic_samples=np.where(dists >= r) and np.where(dists <= r+dr)
		epistatic_samples=np.reshape(epistatic_samples,-1)
		print("penetrance classes with epistasis - should be in epi_table")
		print(np.unique(np.core.defchararray.add(dataframe.loc[loci,indexs].values[np.where(epistatic_samples==1)].astype(str),dataframe.loc[locj,indexs].values[np.where(epistatic_samples==1)].astype(str))))

		print(count_nonzero(epistatic_samples)/float(epistatic_samples.size))
		print("percentage non-zero above")
		vals=pd.Series(epistatic_samples, index=indexs,name=str('interaction_'+str(loci)+'_'+str(locj)))
		print(vals)
		print(dataframe)
		#sys.exit(0)
		#dataframe=pd.concat([dataframe,vals])
		dataframe=dataframe.append(vals)
		print(dataframe)
		#sys.exit(0)
		'''
		B = np.where(A > 0.5, 1, 0)#where A > 0.5 return 1 else 0
		for i in range(len(indexs)):
				if dataframe.loc['interaction'+str(loci)+'_'+str(locj)][indexs[i]]>=1.5:
						dataframe.loc['interaction'+str(loci)+'_'+str(locj)][indexs[i]]=1
		'''
		return dataframe
def XOR_model(dataframe,loci,locj,indexs):
		'''
		this function takes the dataframe and recalculates the interaction term for two epistatic SNPs at loci and locj in the dataframe with an XOR two locus model
		pass dataframe (episnps)
		pass loci and locj randepis[loci] randepis[locj]
		pass indexs = samples i.e. only numeric calculations
		double het model
		'''
		print("running XOR-model of epistatic interactions")
		assert isinstance(indexs,(list, tuple, np.ndarray))
		vals=dataframe.loc[loci,indexs].values.astype(int)
		vals2=dataframe.loc[locj,indexs].values.astype(int)
		vals=vals.astype(str)
		vals2=vals2.astype(str)
		multiples=np.core.defchararray.add(vals,vals2)
		#print(multiples)
		epi_table=['02','12','20','21']
		#sys.exit(0)	
		#multiples=np.multiply(dataframe.loc[loci,indexs].values, dataframe.loc[locj,indexs].values)
		epistatic_samples=np.where(np.isin(multiples, epi_table),1,0)#(dataframe.loc['interaction'+str(loci)+'_'+str(locj)][indexs].values >= 1.5,1,0)
		print("penetrance classes with epistasis - should be in epi_table")
		print(np.unique(np.core.defchararray.add(dataframe.loc[loci,indexs].values[np.where(epistatic_samples==1)].astype(str),dataframe.loc[locj,indexs].values[np.where(epistatic_samples==1)].astype(str))))
		epistatic_samples=np.reshape(epistatic_samples,-1)
		print(count_nonzero(epistatic_samples)/float(epistatic_samples.size))
		print("percentage non-zero above")
		vals=pd.Series(epistatic_samples, index=indexs,name=str('interaction_'+str(loci)+'_'+str(locj)))
		dataframe=dataframe.append(vals)
		#print(dataframe)
		return dataframe

def DD_model(dataframe,loci,locj,indexs):
		'''
		this function takes the dataframe and recalculates the interaction term for two epistatic SNPs at loci and locj in the dataframe
		pass dataframe (episnps)
		pass loci and locj randepis[loci] randepis[locj]
		pass indexs = samples i.e. only numeric calculations
		'''
		print("running DD-model of epistatic interactions")
		#assert isinstance(loci, int)#checking types
		assert isinstance(indexs,(list, tuple, np.ndarray))
		#assert isinstance(locj, int)
		dataframe.loc['interaction'+str(loci)+'_'+str(locj)]=pd.Series(dataframe.loc[loci].values*dataframe.loc[locj],index=indexs)#hadamard product of values
		for i in range(len(indexs)):#rescaling values to required depth
				if dataframe.loc['interaction'+str(loci)+'_'+str(locj)][indexs[i]]>=0.5:
						dataframe.loc['interaction'+str(loci)+'_'+str(locj)][indexs[i]]=1
		return dataframe

def environ(standarderror,length):#returnsvalues for environment variable
		return np.random.normal(scale=standarderror,size=length)
def save_chunks(dataframe,count,savename):
	assert isinstance(savename,str)
	name =[x for x in globals() if globals()[x] is dataframe][0]
	if count ==0:#saving file name with initial loaded chunk
		print("saving first chunk {} for {}".format(count,name))
		dataframe.to_csv(name+".csv")
	else:#appending to csv without headers
		print("appending chunk {} for {}".format(count,name))
		dataframe.to_csv(name+".csv",mode='a',headers=False)
def load_maf_snp(dataframe,neff,maf_constraint,nsnps):
	
	#loads dataframe name at given rows
	randints=np.random.choice(nsnps, size=nsnps, replace=False)	
	df=pd.DataFrame()
	df_snps=pd.DataFrame()
	cols=pd.read_csv(dataframe,nrows=1,index_col=0)
	print(cols)
	print("loading for {}".format(dataframe))
	assert isinstance (dataframe,str)
	assert isinstance (maf_constraint,float)
	count=0
	while len(df) < neff:
		randnumber=randints[count]#np.random.randint(0,nsnps)
		#print(randnumber)
		count=count+1
		print("reading snp {} out of {}".format(len(df),neff))

		if len(df)<1:
			#first load
			df=pd.read_csv(dataframe,skiprows=randnumber,nrows=1,index_col=0 )
			df.columns=cols.columns
			maf=np.sum(df[samples],axis=1)/(2*num_samples)
			print(maf)
			#df_snps['snp']=maf.index.values
			#print(float(df.loc[: ,'maf']) )
			if maf.values<maf_constraint:#float(df.loc[: , 'maf'])<maf_constraint:
				print("SNP maf too low")
				df=pd.DataFrame()
				df_snps=pd.DataFrame()
			#count=count+1
		else:
			a=pd.read_csv(dataframe,skiprows=randnumber,nrows=1,index_col=0)
			a.columns=cols.columns
			maf=np.sum(a[samples],axis=1)/(2*num_samples)
			print(maf)
			a.index=maf.index
			#count=count+1
			if maf.values<maf_constraint:#float(a.loc[:,'maf'])<maf_constraint:
				print("SNP maf too low")
				continue
			#the correct one
			#df=df.append(a,ignore_index=True)
			#df_snps.append(maf.index.value)
			df=df.append(a)
			#=df.append(a)
			#df=pd.concat([df,a],ignore_index=True)#df.append(a)
			#print(df)
	print("SNPs loaded with maf greater than {}".format(maf_constraint))
	return df#,df_snps
def load_maf_snp_hdf(dataframe,neff,maf_constraint,nsnps):
	
	#loads dataframe name at given rows and maf
	randints=np.random.choice(nsnps, size=nsnps, replace=False)	
	#nrows=len(rows)
	df=pd.DataFrame()
	df_snps=pd.DataFrame()
	#df_len=dataframe.shape[0]
	cols=pd.read_hdf(dataframe,'data',start=0,stop=1,index_col=0)
	print(cols)
	
	#print(cols)
	print("loading for {}".format(dataframe))
	count=0
	assert isinstance (dataframe,str)
	assert isinstance (maf_constraint,float)
	while len(df) < neff:
		#randnumber=np.random.randint(0,nsnps)
		#print(randnumber)
		print("reading snp {} out of {}".format(len(df)+1,neff))
		randnumber=randints[count]#np.random.randint(0,nsnps)
		#print(randnumber)
		print(randnumber)
		count=count+1
		if len(df)<1:
			#first load
			df=pd.read_hdf(dataframe,'data',start=randnumber,stop=randnumber+1,index_col=0 )
			df.columns=cols.columns
			maf=np.sum(df[samples],axis=1)/(2*num_samples)
			print(maf)
			df['maf']=maf.values
			#df_snps['snp']=maf.index.values
			#print(float(df.loc[: ,'maf']) )
			if maf.values<maf_constraint:#float(df.loc[: , 'maf'])<maf_constraint:
				print("SNP maf too low")
				df=pd.DataFrame()
				#df_snps=pd.DataFrame()
			elif df.index.values in effect_snps.index.values:
				print("repreat index")
				df=pd.DataFrame()
		else:
			a=pd.read_hdf(dataframe,'data',start=randnumber,stop=randnumber+1,index_col=0)
			a.columns=cols.columns
			maf=np.sum(a[samples],axis=1)/(2*num_samples)
			print(maf)
			a['maf']=maf.values
			a.index=maf.index
			if maf.values<maf_constraint:#float(a.loc[:,'maf'])<maf_constraint:
				print("SNP maf too low")
				continue
			elif a.index.values in effect_snps.index.values:
				print("repeat index")
				continue
			#the correct one
			#df=df.append(a,ignore_index=True)
			#df_snps.append(maf.index.value)
			df=df.append(a)
			#=df.append(a)
			#df=pd.concat([df,a],ignore_index=True)#df.append(a)
			#print(df)
	print("SNPs loaded with maf greater than {}".format(maf_constraint))
	return df#


def load_any_snp(dataframe,neff,nsnps):

	#	loads dataframe name at given rows

	df=pd.DataFrame()
	cols=pd.read_csv(dataframe,nrows=1,index_col=0)
	print(cols)
	assert isinstance (dataframe,str)
	print("loading for {}".format(dataframe))
	while len(df) < neff:#load desired number of SNPs
		randnumber=np.random.randint(0,nsnps)#randomly choose SNP
		print("reading snp {} out of {}".format(len(df),neff))
		if len(df)<1:#load SNPs and create dataframe
			df=pd.read_csv(dataframe,skiprows=randnumber,nrows=1,index_col=0 )
			df.columns=cols.columns
			#df=pd.DataFrame()
		else:#append snp to dataframe
			a=pd.read_csv(dataframe,skiprows=randnumber,nrows=1,index_col=0)
			a.columns=cols.columns
			#df=df.append(a,ignore_index=True)
			df=df.append(a)
	print("SNPs loaded")
	return df

def load_any_snp_hdf(dataframe,neff,nsnps):
	
	#loads dataframe name at given rows
	
	df=pd.DataFrame()
	cols=pd.read_hdf(dataframe,'data',start=0,stop=1,index_col=0)
	print(cols)
	assert isinstance (dataframe,str)
	print("loading for {}".format(dataframe))
	while len(df) < neff:#load desired number of SNPs
		randnumber=np.random.randint(0,nsnps)#randomly choose SNP
		print("reading snp {} out of {}".format(len(df),neff))
		if len(df)<1:#load SNPs and create dataframe
			df=pd.read_hdf(dataframe,'data',start=randnumber,stop=randnumber+1,index_col=0 )
			df.columns=cols.columns
			maf=np.sum(df[samples],axis=1)/(2*num_samples)
			print(maf)
			print("doing maf here")
			df['maf']=maf.values
			#df=pd.DataFrame()
		else:#append snp to dataframe
			a=pd.read_hdf(dataframe,'data',start=randnumber,stop=randnumber+1,index_col=0)
			a.columns=cols.columns
			#df=df.append(a,ignore_index=True)
			maf=np.sum(a[samples],axis=1)/(2*num_samples)
			print(maf)
			a['maf']=maf.values
			df=df.append(a)
	print("SNPs loaded")
	return df
def load_rows_small(dataframe,neff,maf_constraint,nsnps,ncols):
		
	#loads dataframe name at given rows
	
	#nrows=len(rows)
	df=pd.DataFrame()
	#df_len=dataframe.shape[0]
	cols=pd.read_csv(dataframe,nrows=1,index_col=0,usecols=ncols)
	#print(cols)
	print(cols)
	assert isinstance (dataframe,str)
	assert isinstance (maf_constraint,float)
	while len(df) < neff:
		randnumber=np.random.randint(0,nsnps)
		print(randnumber)
		print("reading snp {} out of {}".format(len(df),neff))

		if len(df)<1:
			#first load
			df=pd.read_csv(dataframe,skiprows=randnumber,nrows=1,index_col=0,usecols=ncols )
			df.columns=cols.columns
			maf=np.sum(df)/(2*num_samples)
			print(maf)
			#print(float(df.loc[: ,'maf']) )
			if maf<maf_constraint:#float(df.loc[: , 'maf'])<maf_constraint:
				print("SNP maf too low")
				df=pd.DataFrame()
		else:
			a=pd.read_csv(dataframe,skiprows=randnumber,nrows=1,index_col=0,usecols=ncols)
			a.columns=cols.columns
			maf=np.sum(df)/(2*num_samples)
			print(maf)

			if maf<maf_constraint:#float(a.loc[:,'maf'])<maf_constraint:
				print("SNP maf too low")
				continue
			df=df.append(a,ignore_index=True)
			#df=pd.concat([df,a],ignore_index=True)#df.append(a)
			#print(df)


	return df

def main():
		from datetime import date
		today = date.today()
		print("Today's date:", today)
		print(str(sys.argv))
		
		currentDirectory = os.getcwd()
	
		print(currentDirectory)
		start=time.time()#TIMING RUN
		print(str(sys.argv))
		my_path = str(os.getcwd()) # Figures out the absolute path for you in case your working directory moves around.
		parser = argparse.ArgumentParser()
		parser.add_argument('--option', action='store', default='epistatic',
					help='Interaction types')
		parser.add_argument('--beta_norm', action='store', default=0.0,
					help='Mean for baseline  beta distribution')
		parser.add_argument('--sigbeta', action='store', default=0.1,
					help='Standard Deviation for beta distribution')
		parser.add_argument('--siggamma', action='store', default=0.1,
					help='Standard Deviation for gamma (interaction term)')
		parser.add_argument('--loadstr', action='store', default='ten_k_genotypes',
					help='Name of file to load')
		parser.add_argument('--runno', action='store', default=1,
					help='Number of run in script repeats')
		parser.add_argument('--signoise', action='store', default=1.0,
					help='Standard Deviation for Noise term')
		parser.add_argument('--rescale', action='store', default=1,
					help='Recaling MAF dependent distribution')
		parser.add_argument('--sigmaenviro', action='store', default=0.3,
					help='Standard Deviation for Environmental effect')
		parser.add_argument('--betaenviro', action='store', default=5,
					help='scale of beta vs mubaseline for Environmental effect')
		parser.add_argument('--NSNPs', action='store', default=100,
					help='Number of effect SNPs')
		parser.add_argument('--outdir', action='store', default='test',
					help='directory for output')
		parser.add_argument('--NEPIs', action='store', default=20,
					help='Number of epistatic SNPs')
		parser.add_argument('--heritability', action='store', default=0.5,
					help='Genetic heritability of phenotype')
		parser.add_argument('--epistasis_levels', action='store', default=0.2,
					help='proportion of variance is epistatic')
		parser.add_argument('--if_epi', action='store', default=1,
					help='if epistatic')
		parser.add_argument('--seed',action='store',default=10,
					help="Number to set seed")
		parser.add_argument('--epi_model', action='store', default="m16",
					help='epistatic two locus model, pick T, RR, or m16 (double het) is default')
		parser.add_argument('--if_main', action='store', default=True,
					help='If false, then no main only effects are exhibited')
		parser.add_argument('--version', action='version', version='%(prog)s 1.0')
		results = parser.parse_args()
		paramdict=vars(results)#dictionary of all parsed arguments 

		print('option for run is  = {} \nbeta for environment is {} \nsigma beta = {} \nsigma gamma = {}\nsigma noise = {}\nrun number = {} \ninput file = {}'.format(results.option,results.betaenviro,results.sigbeta,results.siggamma,results.signoise,results.runno,results.loadstr))
		sigmanoise=results.signoise#STANDARD ERROR FOR THE NOISE
		option=str(results.option)
		sigmabeta=float(results.sigbeta)
		beta_baseline=float(results.beta_norm)
		sigmagamma=float(results.siggamma)
		loadstring=str(results.loadstr)
		runno=int(results.runno)
		rescaleval=float(results.rescale)
		sigmaenvironment=float(results.sigmaenviro)
		betaenv=float(results.betaenviro)	  
		betaenv=float(betaenv*beta_baseline)#environmental strength scaled w.r.t baseline beta
		outstring=os.getcwd()+"/"+str(results.outdir)
		herit=float(results.heritability)
		epi_level=float(results.epistasis_levels)
		print(epi_level)
		#sys.exit(0)
		if_main=str(results.if_main)
		if_epi=int(results.if_epi)
		epi_model=str(results.epi_model)
		seed=int(results.seed)
		np.random.seed(seed)
		
		if epi_model not in ["RR","T","XOR","m16"]:
			print("invalid method for two locus epistasis")
			print("must be either m16 (default)\nor T\nor XOR\nor RR\nexiting program")
			sys.exit(0)
		else:
			print("The epistatic model being simulated is the {} model".format(epi_model))
		print(if_main)
		print("Heritability of Phenotype is aiming for: {}".format(herit))
		print(outstring)
		print("saving parameters")
		print(paramdict)
		print(if_main,"if main is being called")
		os.makedirs(outstring, exist_ok=True)
		with open(outstring+'/args.txt', 'w') as fp:#saving parameters
			w=csv.DictWriter(fp,paramdict.keys())
			w.writeheader()
			w.writerow(paramdict)
		numeff=int(results.NSNPs)
		numepi=int(results.NEPIs)
		if numepi%2 != 0 or numepi <=0:
			print("number of epistatic SNPs cannot be negative or odd\nrequested number of epistatic snps {}".format(numepi))
			print("exiting program")
			sys.exit(0)
		scale_factor=float((numeff+numepi)/100)#how many times is numeff SNPs bigger than 100

		print("Reading genotypes")
		import h5py
		chunk_file = h5py.File(loadstring, "r")
		chunk=chunk_file['data']
		print(chunk['/data/axis0'].shape)
		print(chunk['/data/axis1'].shape)
		nrows=chunk['/data/axis1'].shape[0]
		ncols=chunk['/data/axis0'].shape[0]
		col=pd.read_hdf(loadstring,start=0,stop=0,index_col=0)
		print("rows {} \ncolumns {}".format(nrows,ncols))
		global num_samples
		num_samples=int(ncols)

		count=0
		global samples
		samples = col.columns#[1:-1]#ALL OF THE SAMPLE IDS
		print(col.columns)
		print(samples)
		#print(samples)
		#sys.exit(0)
		dfpheno=pd.Series(0,index=samples,dtype=float)#CREATING THE  phenotype lists

		global effect_snps
		effect_snps=load_any_snp_hdf(loadstring,numeff,nrows)
		epi_snps=load_maf_snp_hdf(loadstring,numepi,0.2,nrows)
		SNPs_effect=effect_snps.index
		SNPs_epi=epi_snps.index
		print(effect_snps)
		print(epi_snps)
		print(SNPs_effect)
		print(SNPs_epi)
		
		scalingconst=rescaleval*float(np.mean(effect_snps.maf.values))
		mu_eqn=beta_baseline/scale_factor
		#now we are hitting interactions
		randepis=np.arange(numepi)
		repeats=1#2
		rsqenv=np.zeros(repeats)
		for l in range(repeats):
				startnew=time.time()
				if option == "epistatic":
						print("epstatic interactions")
						epi_snps['beta']=pd.Series(0,index=epi_snps.index,dtype=float)

						for j in range(0,int(numepi),2):#ilooping over each interaction
								if epi_model == "T":
									print("T model")
									#epi_snps=t_model(epi_snps,j,j+1,samples.values)
									epi_snps=t_model(epi_snps,epi_snps.index.values[j],epi_snps.index.values[j+1],samples.values)
								elif epi_model == "RR":
									print("RR Model")
									epi_snps=RR_model(epi_snps,epi_snps.index.values[j],epi_snps.index.values[j+1],samples.values)
								elif epi_model == "XOR":
									print("XOR Model")
									epi_snps=XOR_model(epi_snps,epi_snps.index.values[j],epi_snps.index.values[j+1],samples.values)
								else:
									print("m16 Model")
									epi_snps=m16_model(epi_snps,epi_snps.index.values[j],epi_snps.index.values[j+1],samples.values)
						print(epi_snps)
						#sys.exit(0)
						for i in range(numepi):
							epi_snps['beta'][i]=float(np.random.normal(loc=float(beta_baseline/scale_factor),scale=(sigmabeta/scale_factor),size=1))#epibetatest[i]#float(np.random.normal(loc=sigmabeta,scale=sigmabeta,size=1))
						#print("\n\n\n\n\nHERE are betas for epistatic snps\n\n\n\n\n")
						epielements=np.arange(numepi).astype('int')
						effect_snps['beta']=pd.Series(0,index=effect_snps.index,dtype=float)
						#-----------------------interactions!----------------------------------------------------------------
						if if_main=="True":
							print('main eff snps present')
							effect_snps['beta']=np.random.normal(loc=float(beta_baseline/scale_factor),scale=(sigmabeta/scale_factor),size=numeff)
							#drawn from a slightly positive mu and with a variance that scales w.r.t total number of SNPs
				else:#here we are doing the normal additive SNPS
						print("non epistatic interactions")
						for i in randnums:
								mu=float(scalingconst)/(float(effect_snps['maf'][i]))
								effectsnps['beta'][i]=np.random.normal(loc=mu,scale=sigmabeta)			
				if epi_level==1:
					print("only epistatic contributions")
					print("This option hasn't been implemented")
					sys.exit(0)
					effect_snps['beta']=effect_snps['beta']*0
					epi_snps['beta']=epi_snps['beta']*0
					print("these should be zero")
					print(effect_snps['beta'].values[0:5],epi_snps['beta'].values[0:5])	
				if option == "epistatic":
						#old epistatic effects
						#epi_var=np.random.normal(loc=0.0,scale=sigmagamma,size=(int(numepi/2)))#for each interaction, draw random normal
						epi_var=np.random.normal(loc=float(beta_baseline/scale_factor),scale=(sigmagamma/scale_factor),size=(int(numepi/2)))
						epi_var=epi_var.reshape(int(numepi/2),1)
						if if_epi==0:
							print("no epistasis so setting betas to zero")
							epi_var=epi_var*0
							epi_level=epi_level*0
						if epi_level==0:
							print("no epistasis so setting betas to zero")
							epi_var=epi_var*0
							if_epi=0

						print(epi_var.shape)
						print("calculating interactions")
						interactions=epi_snps.drop(SNPs_epi)
						interactions=interactions[samples]
						print(interactions.shape)
						print("calculating interaction effects")
						epi_val=np.multiply(epi_var,interactions)#part 3 of equation
						print(epi_val.shape)
						print("summing effects")
						epistatic=np.sum(epi_val,axis=0)
						print(epistatic.shape)
						print(epi_snps[samples].shape)
						print(effect_snps.beta.values.T.shape)
						print("calculating main effects")
						print(effect_snps[samples].values.shape)#(2, 100000)
						#print(effect_snps['beta'][samples].values.shape)#(100k,)
						#print(effect_snps['beta'].values.shape)#(2,)
						main_eff=np.matmul(effect_snps.beta.values.T,effect_snps[samples].values)#effect_snps[samples].values,effect_snps['beta'][samples].values)#should be (Nxp x px1)=Nx1
						print(main_eff.shape)
						print("main effect calculated")
						print("shapes betas {}\n samples {}".format(epi_snps.loc[SNPs_epi].beta.values.T,epi_snps.loc[SNPs_epi][samples].values))
						epi_eff=np.matmul(epi_snps.loc[SNPs_epi].beta.values.T,epi_snps.loc[SNPs_epi][samples].values)
						print("epistatic main effects calculated")
						main_genetics=np.add(epi_eff, main_eff)
						var_main=np.var(main_eff)
						var_main_epi=np.var(epi_eff)
						var_interactions=np.var(epistatic)
						sigma_noise=np.sqrt((var_main+var_main_epi +var_interactions)*((1/float(herit)) -1))
						total_genetics=epistatic+main_genetics
						#epistatic + main_genetics
						if if_epi==0:
							epistatic=epistatic*0
							epi_level=0
						if epi_level==0:
							epistatic=epistatic*0
							if_epi=0
						if epi_level==1:
							print("all epistasis, must be no linear main genetics")
							total_genetics=total_genetics-main_genetics
							main_genetics=main_genetics*0
						if if_main==False:
							main_genetics=main_genetics*0	
							print("no main genetics")
						#covering flags for epistasis presence
						environmentaleffect=environ(sigma_noise,len(samples))
						all_effects=effect_snps[samples].append(epi_snps.loc[SNPs_epi][samples])
						#total_genetics=epistatic+main_genetics
						modelepi=sm.OLS(total_genetics,sm.add_constant(all_effects.values.T))
						resultsepi=modelepi.fit()
						#print(resultsepi.summary())
						#print(all_effects)#(effectsnps and epi snps) x samples
						all_genetics_epi=all_effects.append(epi_val)
						print(all_genetics_epi)	
						
						print("need to rescale epistatic contributions here")
						#total_genetics=epistatic+main_genetics
						oldmodelepi=sm.OLS(total_genetics,sm.add_constant(epistatic))#main_genetics))#,epistatic)
						oldresultsepi=oldmodelepi.fit()
						print("multi rsq {} ".format(resultsepi.rsquared))
						print("old rsq {}".format(oldresultsepi.rsquared))
						print(epistatic.values)
						print(main_genetics)
						print(total_genetics)
						#sys.exit(0)
						factor_scale_noise=20
						rsq_epi=resultsepi.rsquared#resultsepi.rsquared
						#rsq_epi=explained_variance_score(total_genetics,sm.add_constant(epistatic))#main_genetics))#(y_true, y_pred) 
						#herit_est=1-rsq_epi
						err=(100*(np.absolute(rsq_epi-epi_level,dtype=np.float32)))/epi_level
						print("currently {}% error from desired heritability\n{} vs actual {}".format(err,rsq_epi,epi_level))
						decay=0
						scale=1
						direction=0
						change=False

						
						while err > 0.5 and if_epi !=0:#5%error and only fix if there are main effects
							if rsq_epi <= epi_level:
								#print(epistatic)
								epistatic=np.multiply(epistatic,factor_scale_noise)
								scale=scale*factor_scale_noise
								#print(epistatic)
								#print(pheno[0:1])
								if direction == 1:
									#print('direction change')
									change=True#now reverse so weaken noise factor
								direction=-1
							else:
								epistatic=np.divide(epistatic,factor_scale_noise)
								#print(epistatic)
								scale=float(scale)/factor_scale_noise
								if direction == -1:
									#print("direction change")
									change = True
								direction=1
								#pheno=pheno+environmentaleffect
								#print(pheno[0:1])
							#print(direction)
							total_genetics=epistatic+main_genetics
							modelepi=sm.OLS(total_genetics,sm.add_constant(epistatic))#main_genetics))#epistatic)
							resultsepi=modelepi.fit()
							print(epistatic.values[0:5])
							#	print(resultsenv.summary())
							#print(resultsepi.rsquared)
							reg = LinearRegression(fit_intercept=True).fit(epistatic.values.reshape(-1, 1),total_genetics)#all_effects.values.T, pheno)
							linr2=reg.score(epistatic.values.reshape(-1, 1),total_genetics)#all_effects.values.T, pheno)#X, y)

							rsq_epi=float(resultsepi.rsquared)#resultsepi.rsquared
							#rsq_epi=1-explained_variance_score(total_genetics,main_genetics)##herit_est=1-rsq_env
							err=(100*(np.absolute(rsq_epi-epi_level,dtype=np.float32)))/epi_level#(100*(rsq_epi-epi_level))/epi_level
							#print("currently {}% error from desired level of epistatic contribution of {}".format(err,epi_level))
							print(rsq_epi,linr2,epi_level)
							decay=decay+1
							#print("old factor {}".format(factor_scale_noise))
							if change == True:
								#print('rescaling scale')
								factor_scale_noise=float(factor_scale_noise*0.9)#(0.5/decay))#(1.0 -float(decay)/100)
								if factor_scale_noise<1:
									factor_scale_noise=1.0001
									change=False
								print("new factor {}".format(factor_scale_noise))

						pheno=np.add(total_genetics,environmentaleffect)
						print(pheno.shape)
						print(pheno)
				else:#nonepistatic phenotypes	
						print("this has not been implemented")
						sys.exit(0)
						for i in range(len(samples)):
								#phenotype = X.Beta + Gaussian Noise(0,1)
								#dfpheno[i]=float(np.dot(effect_snps[samples[i]].values,effect_snps.beta.values)) + float(np.random.normal(loc = 0.0, scale= sigmanoise, size=1))
								pheno[i]=float(np.dot(effect_snps[samples[i]].values,effect_snps.beta.values))+float(np.dot(epi_snps.loc[randepis][samples[i]].values,epi_snps.loc[randepis].beta.values))+environmentaleffect[i]

								#dfpheno2[i]=float(np.dot(effectsnps[samples[i]].values,effectsnps.beta.values)) + float(np.random.normal(loc = 0.0, scale= sigmanoise, size=1))
				#df
				#-----------------------------------SAVING PHENOTYPE-------------------------------------------------------------
				alleff=int(len(epi_snps)+len(effect_snps))
				print("saving phenotype")
				#pickling phenotype
				#print(effectsnps)
				print(outstring)
				print("checking heritability and rescaling")

				#HERE WE WILL INSERT A FITTING FOR THE ENVIRONMENTAL EFFECT
				no_noise=pheno-environmentaleffect
				all_effects=effect_snps[samples].append(epi_snps.loc[SNPs_epi][samples])
				print(epistatic.values)
				print(main_genetics)
				print(total_genetics.values)
				total_genetics=epistatic+main_genetics
				if if_epi==0:#no epistasis
					modelgen=sm.OLS(pheno,sm.add_constant(all_effects.values.T))#total_genetics,all_effects.values.T)
					reg = LinearRegression(fit_intercept=True).fit(all_effects.values.T, pheno)
					linr2=reg.score(all_effects.values.T, pheno)
				else:#including epistasis
					modelgen=sm.OLS(pheno,sm.add_constant(all_genetics_epi.values.T))#all_genetics_epi=
					reg = LinearRegression(fit_intercept=True).fit(all_genetics_epi.values.T, pheno)
					linr2=reg.score(all_genetics_epi.values.T, pheno)#X, y)
				#modelenv=sm.OLS(pheno,environmentaleffect)
				resultsgen=modelgen.fit()
				herit_est=resultsgen.rsquared
				print(herit_est,linr2)
				
				print(resultsgen.summary())#print(resultsenv.summary())
				factor_scale_noise=1.25
				#error between trait heritabiliy and desired heritability
				err=(100*(np.absolute(herit_est-herit,dtype=np.float32)))/herit
				print("currently {}% error from desired heritability\n{} vs actual {}".format(err,herit_est,herit))
				decay=0
				flip=0
				change=False
				while err > 0.5 and if_epi==0:#5%error #no epi
					if herit_est <= herit:
						environmentaleffect=np.divide(environmentaleffect,factor_scale_noise)#scaling noise
						betaenv=float(betaenv)/factor_scale_noise
						pheno=np.add(total_genetics,environmentaleffect)#pheno=pheno+environmentaleffect
						if flip <=0:
							flip =1
							change=True
					else:
						environmentaleffect=np.multiply(environmentaleffect,factor_scale_noise)#scaling noise
						betaenv=float(betaenv)*factor_scale_noise
						pheno=np.add(total_genetics,environmentaleffect)#pheno+environmentaleffect
						if flip >= 0:
							flip=-1
							change = True
						
					total_genetics=epistatic+main_genetics
					modelgen=sm.OLS(pheno,sm.add_constant(all_effects.values.T))
					resultsgen=modelgen.fit()
					herit_est=resultsgen.rsquared
					
					reg = LinearRegression(fit_intercept=True).fit(all_effects.values.T, pheno)
					linr2=reg.score(all_effects.values.T, pheno)#X, y)
					print(herit_est,linr2, herit )
					err=(100*(np.absolute(herit_est-herit,dtype=np.float32)))/herit
					#print("currently {}% error from desired heritability of {}".format(err,herit))
					if change == True:
						print("changing side",factor_scale_noise)
						decay=decay+1
						flip=flip*-1
						factor_scale_noise=factor_scale_noise*(1.0 - (0.5*float(decay))/100)
						if factor_scale_noise <1:
							factor_scale_noise=1.001
						change=False
				while err > 0.5 and if_epi!=0:#5%error with epistasis
					if herit_est <= herit:
						environmentaleffect=np.divide(environmentaleffect,factor_scale_noise)
						betaenv=float(betaenv)/factor_scale_noise
						pheno=np.add(total_genetics,environmentaleffect)#pheno=pheno-environmentaleffect
						if flip <=0:
							flip =1
							change=True
					else:
						environmentaleffect=np.multiply(environmentaleffect,factor_scale_noise)
						betaenv=float(betaenv)*factor_scale_noise
						pheno=np.add(total_genetics,environmentaleffect)#pheno+environmentaleffect
						if flip >= 0:
							flip=-1
							change = True
						
					total_genetics=epistatic+main_genetics
					modelgen=sm.OLS(pheno,sm.add_constant(all_genetics_epi.values.T))
					resultsgen=modelgen.fit()
					herit_est=resultsgen.rsquared
					
					reg = LinearRegression(fit_intercept=True).fit(all_genetics_epi.values.T, pheno)
					linr2=reg.score(all_genetics_epi.values.T, pheno)#X, y)
					print(herit_est,linr2, herit )
					err=(100*(np.absolute(herit_est-herit,dtype=np.float32)))/herit
					#print("currently {}% error from desired heritability of {}".format(err,herit))
					if change == True:
						print("changing side",factor_scale_noise)
						decay=decay+1
						flip=flip*-1
						factor_scale_noise=factor_scale_noise*(1.0 - (0.5*float(decay))/100)
						if factor_scale_noise <1:
							factor_scale_noise=1.001
						change=False

				#
				with open(outstring+"/fitting_summary_repeat_herit_"+str(l)+".csv", 'w') as f:
					f.write(resultsgen.summary().as_csv())
				
				print(epi_snps,epi_snps.loc[interactions.index,'beta'])

				#sys.exit(0)
				print("rescaling")
				print(outstring)
				effect_snps.to_csv(outstring+"/effectsnps"+"numsnp"+str(alleff)+'sigmabeta'+str(sigmabeta)+'sigmagamma'+str(sigmagamma)+'sigmanoise'+str(sigmanoise)+'-'+str(runno)+'repeats'+str(repeats)+'-'+str(l)+".csv")
				#epi_snps.to_csv(outstring+"/episnps"+"numsnp"+str(alleff)+'sigmabeta'+str(sigmabeta)+'sigmagamma'+str(sigmagamma)+'sigmanoise'+str(sigmanoise)+'-'+str(runno)+'repeats'+str(repeats)+'-'+str(l)+".csv")
				epistatic.to_csv(outstring+"/epitatic_contributions"+"numsnp"+str(alleff)+'sigmabeta'+str(sigmabeta)+'sigmagamma'+str(sigmagamma)+'sigmanoise'+str(sigmanoise)+'-'+str(runno)+'repeats'+str(repeats)+'-'+str(l)+".csv")
				np.save(outstring+"/unscaled_epistatic_paired_effect_sizes"+"numsnp"+str(alleff)+'sigmabeta'+str(sigmabeta)+'sigmagamma'+str(sigmagamma),epi_var)
				epi_var=epi_var*scale
				epi_snps.loc[interactions.index,'beta']=epi_var.flatten()
				print(epi_var)
				print(epi_snps)
				print(epistatic)

				np.save(outstring+"/epistatic_paired_effect_sizes"+"numsnp"+str(alleff)+'sigmabeta'+str(sigmabeta)+'sigmagamma'+str(sigmagamma),epi_var)
				epi_snps.to_csv(outstring+"/episnps"+"numsnp"+str(alleff)+'sigmabeta'+str(sigmabeta)+'sigmagamma'+str(sigmagamma)+'sigmanoise'+str(sigmanoise)+'-'+str(runno)+'repeats'+str(repeats)+'-'+str(l)+"_fixed.csv")
				#epi_var.to_csv(outstring+"/epistatic_paired_effect_sizes"+"numsnp"+str(alleff)+'sigmabeta'+str(sigmabeta)+'sigmagamma'+str(sigmagamma)+'sigmanoise'+str(sigmanoise)+'-'+str(runno)+'repeats'+str(repeats)+'-'+str(l)+".csv")
				count=count+1
				print("saving noise and phenotype")
				np.save(outstring+"/environmental_effect_repeat"+str(l),environmentaleffect)#saving all environmental variables
				dfpheno=pd.DataFrame(pheno)
				dfpheno.columns=['phenotype']
				dfpheno.index=samples
				dfpheno.to_csv(outstring+"/phenotype_numeff_"+str(numeff)+"_numepi_"+str(numepi)+"_herit_"+str(herit)+"_epi_"+str(epi_model)+"_"+str(epi_level)+".csv")
				with open(outstring+'/scale_epistatic.txt', 'w') as f:
				#with open('filename.txt', 'w') as f:
					f.write('%f\n' %scale)
					f.write('%f\n' %betaenv)
					f.write('%f\n' %herit_est)
					f.write('%f\n' %rsq_epi)
				#	f.write('%d' % scale)
				#	f.write('%d' % s

if __name__ == '__main__':
		main()
