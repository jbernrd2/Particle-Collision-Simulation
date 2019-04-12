"""
Assignment: Homework 6 part 2

Author(s): Jack Bernard

Collaborators:  None

File: Homework_6.2

Date: 3/12/2017

Reference(s):  
Physics 298owl course notes

Physics 298 owl, University of Illinois at Urbana-Champaign, spring 2017
"""
###################################
# Import libraries 
###################################
import numpy as np
import time 
import decimal
decimal.getcontext().prec = 19
###################################
# Define and initialize variables/Arrays
###################################
start_time = time.time()
bins=0
m=decimal.Decimal(5.972*10**24)
earth_radius=decimal.Decimal(6.371*10**6)
G=decimal.Decimal(6.674*10**-11)
tmax=decimal.Decimal(2*24*60*60)
dt=decimal.Decimal(.1)
dtfine=decimal.Decimal(10**-5)
x0=decimal.Decimal(6.471*10**6)                   
x=x0
y=decimal.Decimal(0)
vx=decimal.Decimal(0)
vy=decimal.Decimal(1.08*10**4)
count=0
smax=decimal.Decimal(1e-6)
smin=decimal.Decimal(1e-9)

dxlist=[]
dylist=[]
tlist=[]
dtlist=[]
siglist=[]
rlist=[]
t=decimal.Decimal(0)
sigma=decimal.Decimal(0)
C=decimal.Decimal(0)
circumference=True
count1=decimal.Decimal(0)
###########################################################################
# Set up lists for x,y coordinates
###########################################################################
xlist=[]
ylist=[]
###########################################################################
# Loops through until time t=tmax, uses midpoint integration technique to 
# update position and velocity. Assigns variable sigma, which is the accuracy
# of the approximation. If sigma is too big or too small, it will increase or
# decrease the step size respectively. 
###########################################################################    
while t<tmax:                
    switch=False    
    t+=dt    
    ax=-G*m*x/(x**2+y**2)**(decimal.Decimal(1.5))
    ay=-G*m*y/(x**2+y**2)**(decimal.Decimal(1.5))    
    xlist.append(x)                      
    ylist.append(y)                                          
    xmid,ymid=x+vx*dt*decimal.Decimal(.5),y+vy*dt*decimal.Decimal(.5)                             
    vxmid=vx+ax*dt*decimal.Decimal(.5)                           
    vymid=vy+ay*dt*decimal.Decimal(.5)
    axmid=-G*m*xmid/(xmid**2+ymid**2)**(decimal.Decimal(1.5))                             
    aymid=-G*m*ymid/(xmid**2+ymid**2)**(decimal.Decimal(1.5))      
    dx,dy=vxmid*dt,vymid*dt    
    dvx=axmid*dt    
    dvy=aymid*dt    
    x,y=dx+x,dy+y
    rlist.append(np.sqrt(x**2+y**2))
    vx,vy=dvx+vx,dvy+vy 
    dL=np.sqrt(dx**2+dy**2)
    if circumference==True:    
        C+=dL    
    if x>0 and y>-.1 and y<.1 and t>5000:
        circumference==False
    dxlist.append(dx)
    dylist.append(dy)
    if count>1:
        deltax=2*dxlist[-2]-dxlist[-3]-dxlist[-1]
        deltay=2*dylist[-2]-dylist[-3]-dylist[-1]
        sigma=decimal.Decimal(1/81)*np.sqrt(deltax**2+deltay**2)
        if sigma>smax:
            dt=dt/decimal.Decimal(2)
            switch=True
            dxlist,dylist=[],[]
        if sigma<smin:
            dt=decimal.Decimal(2)*dt
            switch=True
            dxlist,dylist=[],[]
    if y>-25000 and y<25000 and x>0 and t>1000:
        dt=dtfine
        bins+=1
    siglist.append(sigma)
    tlist.append(t)
    dtlist.append(dt)
    if switch==True:
        count=0
        continue
    count+=1
    count1+=1
    
minval=np.min(rlist[100:])
for i in range(len(rlist)):
    if rlist[i]==minval:
        indexmin=i
    
#########################################################################
# Loop finished, calculate time elapsed and print final values
#########################################################################
    
    
rmax,rmin=np.max(rlist),np.min(rlist)
elapsed_time = time.time() - start_time
print("Maximum Orbital Radius:",rmax)
print("Minimum Orbital Radius:",rmin)
print("Circumference:",C)
print("Number of Bins in the Fine-grained Scan:",bins)
print('Index perigee:',indexmin,'  Value:',rmin)
print("Closure miss difference:",np.sqrt((x0-xlist[indexmin])**2+(ylist[indexmin])**2))
print("Elapsed running time = ", elapsed_time)

'''fig=plt.figure()
ax=fig.gca()
ax.set_xlim(-1.2e8,3e7)
ax.set_ylim(-3e7,3e7)
ax.set_aspect('equal')
plt.plot(xlist,ylist,color='c')'''