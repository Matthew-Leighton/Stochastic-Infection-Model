import numpy as np

###Best-fit parameters:
params=[2,79,0.00153,0.01,0.1976,0.7902,0.2010]

def runmodel(t_max,moi,params,s=1,H=1000):
	tlist=[]
	datalist=[]

	###Parameters
	r_max=params[0]
	x_max=params[1]
	Gamma_a=params[2]
	Gamma_ar=params[3]
	Gamma_x_a=params[4]
	Gamma_x_ar=params[5]
	Gamma_r=params[6]

	###Stochasticity Stuff
	#s is the stochasticity coefficient (default value s=1)
	#Account for different values of stochasticity coefficient
	Gamma_ar/=s
	r_max*=s
	x_max*=s

	###Cell Stuff:
	c=1                #Confluency
	R=0                #Total number of ruffles
	#H is the total number of cells (default H=1000)

	###Bacteria Stuff:
	B_tot=moi*H*s      #Total number of bacteria
	B_f=B_tot        #Number of free bacteria

	#For each cell:
	n_a=np.zeros(H)			#Number of bacteria attached
	n_ar=np.zeros(H)		#Number of bacteria attached to ruffles
	n_x=np.zeros(H)			#Number of bacteria internalized
	n_x_r=np.zeros(H)		#Number of bacteria internalized via ruffles
	r=np.zeros(H)			#Number of ruffles


	#Record Initial Conditions:
	t=0
	data=np.zeros((5,H))
	data[0,:]=n_a
	data[1,:]=n_ar
	data[2,:]=n_x
	data[3,:]=n_x_r
	data[4,:]=r

	tlist.append(t)
	datalist.append(data)

	###############Begin Dynamics:

	while t<t_max:

		###Update rates:
		Gamma_A=Gamma_a*c*B_f

		Gamma_AR=Gamma_ar*c*B_f*R/H

		Gamma_X_A=Gamma_x_a*np.sum((n_a)*(1-((n_x)/x_max)))
		Gamma_X_AR=Gamma_x_ar*np.sum((n_ar)*(1-((n_x)/x_max)))

		Gamma_R=Gamma_r*np.sum((1-(r/r_max))*n_a)

		Gamma_tot=Gamma_A+Gamma_AR+Gamma_X_A+Gamma_X_AR+Gamma_R
		if Gamma_tot==0.00000000000:
			print('equilibrium reached')
			break

		###React Things:
		reaction=np.where(np.random.multinomial(1,[Gamma_A/Gamma_tot,Gamma_AR/Gamma_tot,Gamma_X_A/Gamma_tot,Gamma_X_AR/Gamma_tot,Gamma_R/Gamma_tot]))[0][0]

		if reaction==0:       #attach a bacteria (no ruffles)
			#choose a cell 
			cell=np.random.randint(H)
			B_f-=1
			n_a[cell]+=1

		if reaction==1:        #attach a bacteria to a ruffle
			#choose a cell:
			cell=np.where(np.random.multinomial(1,r/R))[0][0]
			B_f-=1
			n_ar[cell]+=1

		elif reaction==2:      #primary-attached bacteria-->Internalized
			#Choose a cell:
			cell=np.where(np.random.multinomial(1,(Gamma_x_a*(n_a)*(1-((n_x)/x_max)))/Gamma_X_A))[0][0]
			n_a[cell]-=1
			n_x[cell]+=1

		elif reaction==3:      #ruffle-attached bacteria-->Internalized
			#Choose a cell:
			cell=np.where(np.random.multinomial(1,(Gamma_x_ar*(n_ar)*(1-((n_x)/x_max)))/Gamma_X_AR))[0][0]
			n_ar[cell]-=1
			n_x[cell]+=1
			n_x_r[cell]+=1

		elif reaction==4:       #form a ruffle
			#choose a cell:
			cell=np.where(np.random.multinomial(1,(1-(r/r_max))*n_a*Gamma_r/Gamma_R))[0][0]
			r[cell]+=1
			n_a[cell]-=1
			n_ar[cell]+=1
			R+=1

		###Calculate Time Step:
		deltaT=np.random.exponential(1/Gamma_tot)
		t+=deltaT
		###Repeat until t_max

		#Update Statistics
		data[0,:]=n_a/s
		data[1,:]=n_ar/s
		data[2,:]=n_x/s
		data[3,:]=n_x_r/s
		data[4,:]=r/s

		tlist.append(t)
		datalist.append(data)
 
	return np.array(tlist),np.array(datalist)
