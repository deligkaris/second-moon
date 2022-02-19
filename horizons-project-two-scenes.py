#####################################
#AUTHOR: CHRISTOS DELIGKARIS
#PREDICTION OF PATH OF OBJECT ?
#####################################

from __future__ import division
from visual import *
from visual.graph import *

###########################
#SCENES

scene3 = display(title='VECTORS', x=700, y=0, width=350, height=350, center=(5,0,0), background=(1,1,1),autoscale=True)
scene2 = display(title='WITH RESPECT TO EARTH', x=350, y=0, width=350, height=350, center=(5,0,0), background=(1,1,1),autoscale=True)
scene1 = display(title='WITH RESPECT TO SUN', x=0, y=0, width=350, height=350, center=(5,0,0), background=(1,1,1))

###########################
#GET POSITIONS AND VELOCITIES OF SUN

scene1.select()
Sun=sphere(pos=vector(0,0,0),radius=6.4e6,color=color.orange)

###########################
#GET POSITIONS AND VELOCITIES OF SECOND MOON

input=open('horizons-results-OtherMoon-wrt-Sun.txt','r')    #open the file with the Horizons system data
L=input.readlines() #read the file
i=0 #i will represent the line number (0 for the first line, 1 for the second line, etc)
imax=len(L) #imax-1 is the number of lines in the input file
pos={}  #an empty list for positions
vel={}  #an empty list for velocities
t=0
deltat=86400

while i<imax:
    if L[i] == "$$SOE\n":   #when we reach the beginning of the data
        i=i+1
        while 1:
            if L[i] == "$$EOE\n": break #stop reading positions and velocities if end of data is reached
            position = L[i+1].split()   #store the position of the object (need to split the string)
            pos[t] = vector(float(position[0]), float(position[1]), float(position[2]))     #convert the strings to floats and then vector
            pos[t] = 1000.*pos[t]   #convert to meters
            velocity = L[i+2].split()   #store the velocity of the object (need to split the string)
            vel[t] = vector(float(velocity[0]), float(velocity[1]), float(velocity[2]))     #convert the strings to floats and then vector
            vel[t] = 1000.*vel[t]   #convert to m/s
            try:
                SecondMoonObs.pos=pos[t]
            except NameError:
                SecondMoonObs = sphere (radius=1.75e6, color=color.black, make_trail=True, pos=pos[t])
            i=i+3       #move to the next position and velocity of the object
            t=t+deltat
    i=i+1

craft=sphere(radius=1e6,color=color.red, make_trail=True, pos=pos[90*86400])

###########################
#GET POSITIONS AND VELOCITIES OF MOON
    
input=open('horizons-results-Moon-wrt-Sun.txt','r') #open the file with the Horizons system data
L=input.readlines() #read the file
i=0 #i will represent the line number (0 for the first line, 1 for the second line, etc)
imax=len(L) #imax-1 is the number of lines in the input file
posMoon={}  #an empty list for positions
velMoon={}  #an empty list for velocities
t=0
deltat=86400

while i<imax:
    if L[i] == "$$SOE\n":   #when we reach the beginning of the data
        i=i+1
        while 1:
            if L[i] == "$$EOE\n": break #stop reading positions and velocities if end of data is reached
            position = L[i+1].split()   #store the position of the object (need to split the string)
            posMoon[t] = vector(float(position[0]), float(position[1]), float(position[2]))     #convert the strings to floats and then vector
            posMoon[t] = 1000.*posMoon[t]   #convert to meters
            velocity = L[i+2].split()   #store the velocity of the object (need to split the string)
            velMoon[t] = vector(float(velocity[0]), float(velocity[1]), float(velocity[2]))     #convert the strings to floats and then vector
            velMoon[t] = 1000.*velMoon[t]   #convert to m/s
            try:
                MoonObs.pos=posMoon[t]
            except NameError:
                MoonObs = sphere (radius=1.75e6, color=color.red, make_trail=True, pos=posMoon[t])
            i=i+3       #move to the next position and velocity of the object
            t=t+deltat
    i=i+1

Moon=sphere(radius=2*1.75e6,color=color.green, pos=posMoon[90*86400])

###########################
#GET POSITIONS AND VELOCITIES OF EARTH

input=open('horizons-results-Earth-wrt-Sun.txt','r')    #open the file with the Horizons system data
L=input.readlines() #read the file
i=0 #i will represent the line number (0 for the first line, 1 for the second line, etc)
imax=len(L) #imax-1 is the number of lines in the input file
posEarth={} #an empty list for positions
velEarth={} #an empty list for velocities
t=0
deltat=86400

while i<imax:
    if L[i] == "$$SOE\n":   #when we reach the beginning of the data
        i=i+1
        while 1:
            if L[i] == "$$EOE\n": break #stop reading positions and velocities if end of data is reached
            position = L[i+1].split()   #store the position of the object (need to split the string)
            posEarth[t] = vector(float(position[0]), float(position[1]), float(position[2]))     #convert the strings to floats and then vector
            posEarth[t] = 1000.*posEarth[t]   #convert to meters
            try: 
                craft2Obs.pos=pos[t]-posEarth[t]
            except NameError:
                scene2.select()
                craft2Obs=sphere(radius=1e6,color=color.black,make_trail=True, pos=pos[t]-posEarth[t])
                scene1.select()
            velocity = L[i+2].split()   #store the velocity of the object (need to split the string)
            velEarth[t] = vector(float(velocity[0]), float(velocity[1]), float(velocity[2]))     #convert the strings to floats and then vector
            velEarth[t] = 1000.*velEarth[t]   #convert to m/s
            try:
                EarthObs.pos=posEarth[t]
            except NameError:
                EarthObs = sphere (radius=1.75e6, color=color.cyan, make_trail=True, pos=posEarth[t])
            i=i+3       #move to the next position and velocity of the object
            t=t+deltat
    i=i+1

Earth=sphere(radius=6.4e6,color=color.blue,pos=posEarth[90*86400])
scene2.select()
craft2=sphere(radius=1e6,color=color.red,make_trail=True, pos=craft.pos-Earth.pos)

###########################
#CREATE FORCE, MOMENTUM ARROWS

scene3.select()
parr_c=arrow(color=color.blue)
fnetarr_c=arrow(color=color.red)#,shaftwidth=2e6)
fnetpararr_c=arrow(color=color.red)
fnetperarr_c=arrow(color=color.red)
fearr_c=arrow(color=color.cyan)#,shaftwidth=2e6)
fmarr_c=arrow(color=color.white)#,shaftwidth=2e6)
fsunarr_c=arrow(color=color.green)
fscale=1e8
pscale=1e1
scene1.select()

###########################
#GRAPHS

eg = gdisplay(title='ERROR',xtitle='Time [Seconds]',ytitle='Error [meters]',x=0,y=350,width=350,height=350)
ErrorGraph = gcurve(color=color.red)

fng = gdisplay(title='MAGNITUDES OF F NET COMPONENTS', xtitle='Time [Seconds]',ytitle='Force Magnitudes [Newtons]',x=350,y=350,width=350,height=350)
FnetGraph = gcurve(color=color.red)
FnetparGraph = gcurve(color=color.blue)
FnetperGraph = gcurve(color=color.green)

fperg = gdisplay(title='F NET PARALLEL', xtitle='Time [Seconds]',ytitle='F NET PARALLEL [Newtons]',x=700,y=350,width=350,height=350)
FnetparPGraph = gcurve(color=color.blue)
PFnetparGraph = gcurve(color=color.red)

fparg = gdisplay(title='F NET PERPENDICULAR', xtitle='Time [Seconds]',ytitle='F NET PERPENDICULAR [Newtons]',x=1000,y=350,width=350,height=350)
FnetperPGraph = gcurve(color=color.blue)
PFnetperGraph = gcurve(color=color.red)

eng = gdisplay(title='ENERGY-WORK',xtitle='Seconds',ytitle='Joules',x=1000,y=0,width=350,height=350)
EkGraph = gcurve(color=color.magenta)
#Ugraph=gcurve(color=color.yellow)
#Etotal=gcurve(color=color.red)
WGraph = gcurve(color=color.cyan)

enrg = gdisplay(title='ENERGY-WORK',xtitle='Meters',ytitle='Joules',x=1000,y=0,width=350,height=350)
EkrGraph = gcurve(color=color.magenta, dot=True)
EprGraph = gcurve(color=color.green, dot=True)

cg = gdisplay(title='ENERGY-CONSERVATION',xtitle='Seconds',ytitle='Joules',x=1000,y=0,width=350,height=350)
DeltaEGraph = gcurve(color=color.red)

###########################
#CALCULATIONS

G=6.67e-11
Earth.m=6e24
craft.m=15e3
Moon.m=7e22
Sun.m=2e30

craft.v=vel[90*86400]
craft.p=craft.m*craft.v
craft.pinitial = craft.p
craft.pmag=mag(craft.p)
craft.phat = craft.p/craft.pmag
Eki=(mag(craft.pinitial)**2)/(2*craft.m)     #system: craft

Moon.v=velMoon[90*86400]
Moon.p = Moon.m*Moon.v

Earth.v=velEarth[90*86400]
Earth.p=Earth.m*Earth.v

###########################
#CALCULATE FORCES ON CRAFT

r_cE=craft.pos-Earth.pos
rmag_cE=sqrt(r_cE.x**2+r_cE.y**2+r_cE.z**2)
if rmag_cE<Earth.radius:
    print("crashed into Earth")   
rhat_cE=r_cE/rmag_cE
Fmag_cE=G*craft.m*Earth.m/(rmag_cE**2)
F_cE=-Fmag_cE*rhat_cE

r_cM=craft.pos-Moon.pos
rmag_cM=sqrt(r_cM.x**2+r_cM.y**2+r_cM.z**2)
if rmag_cM<Moon.radius:
    print("Crashed into Moon")
rhat_cM=r_cM/rmag_cM
Fmag_cM=G*craft.m*Moon.m/(rmag_cM**2)
F_cM=-Fmag_cM*rhat_cM

r_cS=craft.pos-Sun.pos
rmag_cS=mag(r_cS)
if rmag_cS<Sun.radius:
    print("Crashed into Sun")
rhat_cS=r_cS/rmag_cS
Fmag_cS=G*craft.m*Sun.m/(rmag_cS**2)
F_cS=-Fmag_cS*rhat_cS

Fnet=F_cE + F_cM + F_cS

###########################
#CALCULATE FORCES ON MOON

r_ME = Moon.pos-Earth.pos
rmag_ME=mag(r_ME)
rhat_ME=r_ME/rmag_ME
Fmag_ME = G*Earth.m*Moon.m/(rmag_ME**2)
F_ME=-Fmag_ME*rhat_ME

r_Mc=-r_cM
rmag_Mc=rmag_cM
rhat_Mc=-rhat_cM
Fmag_Mc=Fmag_cM
F_Mc=-F_cM

r_MS = Moon.pos-Sun.pos
rmag_MS=mag(r_MS)
rhat_MS=r_MS/rmag_MS
Fmag_MS = G*Moon.m*Sun.m/(rmag_MS**2)
F_MS=-Fmag_MS*rhat_MS

Fnet_M = F_Mc + F_ME + F_MS

###########################
#CALCULATE FORCES ON EARTH

r_ES = Earth.pos-Sun.pos
rmag_ES=mag(r_ES)
rhat_ES=r_ES/rmag_ES
Fmag_ES = G*Earth.m*Sun.m/(rmag_ES**2)
F_ES=-Fmag_ES*rhat_ES

F_EM = - F_ME
F_Ec = - F_cE

Fnet_E = F_Ec + F_EM + F_ES

##############################

#r_cE=craft.pos-Earth.pos
#rmag_cE=sqrt(r_cE.x**2+r_cE.y**2+r_cE.z**2)
#this is needed in order to compare deltaU (system=craft+earth)
#with the work done (system=craft)
#U0 =-G*Earth.m*craft.m/rmag_cE 

#scene.autoscale=1
#scene.center=vector(2e8,0,0)
#scene.center=Earth.pos
#scene.range=10

###########################
# INITIALIZE LOOP

W_oncraft=0

t=90*86400
deltat=864

while t<1.3*365*24*60*60:
    rate(1e3)
    #scene2.autoscale=False

    craft.phatbefore = craft.phat
    
    ###########################
    #UPDATE POSITIONS OF ALL OBJECTS
  
    craft.p=craft.p+Fnet*deltat
    craft.pos=craft.pos+(craft.p/craft.m)*deltat

    Earth.p = Earth.p + Fnet_E*deltat
    Earth.pos = Earth.pos + (Earth.p/Earth.m)*deltat

    Moon.p = Moon.p + Fnet_M*deltat
    Moon.pos = Moon.pos + (Moon.p/Moon.m)*deltat

    craft2.pos=craft.pos-Earth.pos
    
    ###########################
    #CALCULATE FORCES ON CRAFT

    r_cE=craft.pos-Earth.pos
    rmag_cE=sqrt(r_cE.x**2+r_cE.y**2+r_cE.z**2)
    if rmag_cE<Earth.radius:
        print("crashed into Earth")
        break    
    rhat_cE=r_cE/rmag_cE
    Fmag_cE=G*craft.m*Earth.m/(rmag_cE**2)
    F_cE=-Fmag_cE*rhat_cE

    r_cM=craft.pos-Moon.pos
    rmag_cM=sqrt(r_cM.x**2+r_cM.y**2+r_cM.z**2)
    if rmag_cM<Moon.radius:
        print("Crashed into Moon")
        break
    rhat_cM=r_cM/rmag_cM
    Fmag_cM=G*craft.m*Moon.m/(rmag_cM**2)
    F_cM=-Fmag_cM*rhat_cM

    r_cS=craft.pos-Sun.pos
    rmag_cS=mag(r_cS)
    if rmag_cS<Sun.radius:
        print("Crashed into Sun")
        break
    rhat_cS=r_cS/rmag_cS
    Fmag_cS=G*craft.m*Sun.m/(rmag_cS**2)
    F_cS=-Fmag_cS*rhat_cS

    Fnet=F_cE + F_cM + F_cS

    ########################### 
    #CALCULATE FORCES ON MOON

    r_ME = Moon.pos-Earth.pos
    rmag_ME=mag(r_ME)
    rhat_ME=r_ME/rmag_ME
    Fmag_ME = G*Earth.m*Moon.m/(rmag_ME**2)
    F_ME=-Fmag_ME*rhat_ME

    r_Mc=-r_cM
    rmag_Mc=rmag_cM
    rhat_Mc=-rhat_cM
    Fmag_Mc=Fmag_cM
    F_Mc=-F_cM

    r_MS = Moon.pos-Sun.pos
    rmag_MS=mag(r_MS)
    rhat_MS=r_MS/rmag_MS
    Fmag_MS = G*Moon.m*Sun.m/(rmag_MS**2)
    F_MS=-Fmag_MS*rhat_MS
    
    Fnet_M = F_Mc + F_ME + F_MS

    ###########################
    #CALCULATE FORCES ON EARTH

    r_ES = Earth.pos-Sun.pos
    rmag_ES=mag(r_ES)
    rhat_ES=r_ES/rmag_ES
    Fmag_ES = G*Earth.m*Sun.m/(rmag_ES**2)
    F_ES=-Fmag_ES*rhat_ES
    
    F_EM = - F_ME
    F_Ec = - F_cE

    Fnet_E = F_Ec + F_EM + F_ES

    ###########################
    #CALCULATE F NET COMPONENTS, WORK, ENERGIES

    craft.pmag=mag(craft.p)
    craft.phat=craft.p/craft.pmag
    Fnetpar = dot(Fnet,craft.phat)*craft.phat
    Fnetper = Fnet - Fnetpar

    Ek=(craft.pmag**2)/(2*craft.m)     #system: craft

    #U=-G*Earth.m*craft.m/rmag_cE    #system: craft+earth
    #U=0     #system: craft
    #E=K+U   #system: craft+earth

    deltar=(craft.p/craft.m)*deltat
    W_increment=dot(Fnet, deltar)
    W_oncraft=W_oncraft+W_increment

    #E=K+U+W_oncraft     #system:craft

    #Kgraph.plot(pos=(t,K))
    #this U graph will be compared with the work done (system=craft)
    #Ugraph.plot(pos=(t,U-U0))
    #Ugraph.plot(pos=(t,U))
    #Etotal.plot(pos=(t,E))
    #Wgraph.plot(pos=(t,W_oncraft))

    #if abs(t%86400)<10**(-1):
            #print(t,abs(t%86400),mag(craft.pos-pos[round(t)]))
            #ErrorGraph.plot(pos=(t,mag(craft.pos-pos[round(t)])))
            #ErrorGraph.plot(pos=(t,mag(craft2.pos-(pos[t]-posEarth[t]))))
    FnetGraph.plot(pos=(t,mag(Fnet)))
    FnetparGraph.plot(pos=(t,mag(Fnetpar)))
    FnetperGraph.plot(pos=(t,mag(Fnetper)))
    FnetparPGraph.plot(pos=(t,1e7*dot(Fnet,craft.phat)))
    PFnetparGraph.plot(pos=(t,craft.pmag))
    FnetperPGraph.plot(pos=(t,mag(Fnetper)))
    #PFnetperGraph.plot(pos=(t,100*math.acos(dot(craft.p,craft.pinitial)/(mag(craft.p)*mag(craft.pinitial)))))
    PFnetperGraph.plot(pos=(t,1e5*math.acos(dot(craft.phat,craft.phatbefore))/deltat))
    EkGraph.plot(pos=(t,Ek))
    WGraph.plot(pos=(t,W_oncraft))
    DeltaEGraph.plot(pos=(t,Ek-Eki-W_oncraft))
    EkrGraph.plot(pos=(rmag_cE,Ek))
    #print(t,100*math.acos(dot(craft.p,craft.pinitial)/(mag(craft.p)*mag(craft.pinitial))))

    #trail.append(pos=craft.pos)
    #parr_c.pos=craft.pos
    parr_c.axis=craft.p*pscale
    #fnetarr_c.pos=craft.pos
    fnetarr_c.axis=fscale*Fnet
    fnetpararr_c.axis = fscale*Fnetpar
    fnetperarr_c.axis = fscale*Fnetper
    #print(Fnet,craft.p)
    #fearr_c.pos=craft.pos
    fearr_c.axis=fscale*F_cE
    #fmarr_c.pos=craft.pos
    fmarr_c.axis=fscale*F_cM
    fsunarr_c.axis=fscale*F_cS

    #scene.center=Sun.pos
    #print(t,1.3*365*24*60*60)
    t=t+deltat

