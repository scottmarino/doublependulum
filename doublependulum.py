#!/usr/bin/env python3

# Differential equation solver for double pendulum using symplectic methods
#
# For derivation of Hamiltonian methods refer to LaTeX document
#

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import threading
import sys
import time

initstr = """
                    Symplectic Double Pendulum Simulator

Please enter the bob masses, pendulum lengths, and starting angles in degrees
as such:

ANGLE1 ANGLE2 M1 M2 L1 L2

Or enter just starting angles (other values defaulted to 1):

ANGLE1 ANGLE2

When you are satisfied with your pendulum, type <Run>.
"""
print(initstr)

def inputthread(): 
	global m1          # Declaring and initializing global variable so the
	m1 = 1.0           # user can communicate with the animation function    
	global m2          # in real time
	m2 = 1.0
	global l1
	l1 = 1.0
	global l2
	l2 = 1.0
	global a0
	a0 = np.pi/2
	global b0
	b0 = np.pi/2
	global a0deg
	a0deg = 90
	global b0deg
	b0deg = 90
	global anistart # Communication variable to end first animation
	anistart = False
	while True:
		instr = input('Pendulum Parameters: ')
		if instr.upper() == 'RUN': # user input not case-specific
			anistart = True
			break
		try:
			instrlist = instr.split(' ')
			nlist = [float(x) for x in instrlist]
			a = len(nlist)
			if a == 6:
				a0deg = nlist[0] # Keeping degree for printing EPS file
				b0deg = nlist[1]
				a0 = (a0deg/360)*2*np.pi # Convert to radians
				b0 = (b0deg/360)*2*np.pi
				m1 = nlist[2]
				m2 = nlist[3]
				l1 = nlist[4]
				l2 = nlist[5]
			if a == 2:
				a0deg = nlist[0]
				b0deg = nlist[1]
				a0 = (a0deg/360)*2*np.pi
				b0 = (b0deg/360)*2*np.pi
			if a != 2 and a != 6:
				raise IndexError # Covers all bases of errors
		except ValueError:
			print('\nInvalid input, try again\n',file=sys.stderr)
		except IndexError:
			print('\nInvalid number of inputs, try again\n',file=sys.stderr)

def Hp(a,b,pa,pb,i): # Derivative of H w/ respect to p
	dif = a[i] - b[i]
	# dH/dpa
	num1 = l2*pa[i]-l1*pb[i]*np.cos(dif)
	den = l1*l2*(m1 + m2*(np.sin(dif)*np.sin(dif)))
	val1 = num1/(l1*den)
	# dH/dpb
	num2 = (m1+m2)*l1*pb[i]-m2*l2*pa[i]*np.cos(dif)
	val2 = num2/(l2*m2*den)
	return val1, val2

def Hq(a,b,pa,pb,i): # Derivative of H w/ respect to q
	dif = a[i+1] - b[i+1] # will have been found beforehand
	ta1 = -1*(m1+m2)*g*l1*np.sin(a[i+1])
	tb1 = -m2*g*l2*np.sin(b[i+1])
	
	num2 = pa[i]*pb[i]*np.sin(dif)
	den = l1*l2*(m1 + m2*(np.sin(dif)*np.sin(dif)))
	t2 = num2/den
	
	num3 = (m2*l2*l2*pa[i] + (m1+m2)*l1*l1*pb[i]*pb[i]-2*m2*l1*l2*pa[i]*pb[i]*np.cos(dif))*np.sin(2*dif)
	t3 = num3/(2*l1*l2*den*den)
	
	val1 = ta1 - t2 +t3
	val2 = tb1 + t2 - t3
	return val1, val2 # Note these are negative, so will be added to our algo


def __init__(): # For input animation
	line.set_data([],[])
	x10 = np.sin(a0)
	y10 = -np.cos(a0)
	x20 = np.sin(b0) + x10
	y20 = -np.cos(b0) + y10
	return line, 

def __init2__(): # For simulation
	line2.set_data([],[])
	time_text.set_text('')
	return line2, time_text
	
def framegen(): # Shut off input animation when user types <run>
	global anistart
	i = 0
	while anistart is False:
		i += 1
		yield i,	


def update1(i): # animate function for animation 1
	x1i = l1*np.sin(a0)
	y1i = -l1*np.cos(a0)

	x2i = l2*np.sin(b0) + x1i
	y2i = -l2*np.cos(b0) + y1i
		
	xgi = [0, x1i, x2i]
	ygi = [0, y1i, y2i]
	
	line.set_data(xgi, ygi)
	axlim = l1 + l2 + 1
	ax.set_xlim(-axlim,axlim)
	ax.set_ylim(-axlim,axlim)
	
	global anistart # closes plot at the same time animation finishes
	if anistart == True:
		plt.close()
	return line, 

def update2(i): # animation 2
	xg = [0, x1[i], x2[i]]
	yg = [0, y1[i], y2[i]]
	
	line2.set_data(xg, yg)
	time_text.set_text(time_temp % (i*dt))
	return line2, time_text

T = 20
dt = 0.02
npoints = 1000 # T and dt chosen so we have 1000 points to plot
g = 9.83
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Main sequence~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if __name__ == '__main__':
	thr = threading.Thread(target = inputthread)
	thr.start()
	
	fig = plt.figure()
	ax = fig.add_subplot(111, autoscale_on=False)
	ax.axis('off')
	
	line, = ax.plot([],[], '-o', linewidth=2, color='black')
	ani = animation.FuncAnimation(fig, update1, frames=framegen, interval=dt*1000, blit=True, init_func=__init__,repeat=False)
	plt.show(block=False) # Need this block function or else next lines cant run
	while anistart is False:
		plt.pause(1)
	
	print('\nCalculating path...')
	plt.close(fig) 
	a = np.zeros(npoints) # Arrays for plot calculation
	b = np.zeros(npoints)
	pa = np.zeros(npoints)
	pb = np.zeros(npoints)
	
	a[0] = a0 # Initializing arrays w/ given parameters
	b[0] = b0
	pa[0] = 0
	pb[0] = 0

# Symplectic Euler's method described in accompanying PDF:
	for i in range(npoints-1):
		Hpa, Hpb = Hp(a,b,pa,pb,i)
		a[i+1] = a[i] + dt*Hpa
		b[i+1] = b[i] + dt*Hpb


		Ha, Hb = Hq(a,b,pa,pb,i)
		pa[i+1] = pa[i] + dt*Ha
		pb[i+1] = pb[i] + dt*Hb
	x1 = l1*np.sin(a) # Convert to Cartesian for graph
	y1 = -l1*np.cos(a)

	x2 = l2*np.sin(b) + x1
	y2 = -l2*np.cos(b) + y1
	
	fig2 = plt.figure()
	ax2 = fig2.add_subplot(111, autoscale_on=False)
	ax2.axis('off') # No labels
	
	axlim = l1 + l2 + 1
	ax2.set_xlim(-axlim,axlim)
	ax2.set_ylim(-axlim,axlim)
	
	line2, = ax2.plot([],[], 'o-', linewidth=2, color='black')
	time_temp = 'Time = %.1f seconds' # Keep time, ensure simulation is correct
	time_text = ax2.text(0.05, 0.9, '', transform=ax2.transAxes)

	print('\nSimulation ready')

	ani2 = animation.FuncAnimation(fig2, update2, np.arange(1,npoints), interval=dt*1000, blit=True, init_func=__init2__)

# One more while loop to get user input on saving files
	while True:
		saveinst = input('\nSave .eps figure of trajectory path and angular motion after animation finishes? <y> or <n>:\n')
		
		f3, (ax3,ax4) = plt.subplots(2,figsize=(5,5)) # Data plots
		suptitlestr = 'First and second initial values: %.1f and %.1f degrees,\n %.1f and %.1f kg, %.1f and %.1f meters' %(a0deg, b0deg, m1, m2, l1, l2)
		f3.suptitle(suptitlestr)
		
		ax3.plot(x2 ,y2, color='tab:orange') # Trace subplot
		ax3.plot([0,x1[0],x2[0]],[0,y1[0],y2[0]], 'o-',linewidth=2, color='black')
		ax3.axis('square')
		ax3.set_xlim(-axlim,axlim)
		ax3.set_ylim(-axlim,axlim)
		
		ax4.plot(a,b,color='tab:red') # Angular change subplot
		ax4.set_xlabel('Angle 1 in radians')
		ax4.set_ylabel('Angle 2 in radians')
		print('\nOk. When you are finished viewing the graphs exit out of the figures to \nend program.\n')
		time.sleep(1) # Give user time to read prompt
		try:
			save = saveinst.upper()
			if save == 'Y':
				cond = [int(a0deg),int(b0deg),int(m1),int(l1),int(l2)]
				condstr = '_'.join([str(i) for i in cond])
				savestr = 'dp_'+condstr+'.eps'
				f3.savefig(savestr,format='eps') # save as eps file
				plt.show()
				break
			if save == 'N':
				plt.show()
				break
			else:
				raise ValueError
		except ValueError:
			print('Invalid input', file=sys.stderr)
		else:
			break
		sys.exit()
