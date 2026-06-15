"""
-Gazebo communication functionality preserved for future reuse
-Part lengths measured from axis to axis, eg distance from hip's z axis to knee's xy axis
-Only works for theta is a multiple of pi/2
"""



import tkinter as tk
from math import *
from scipy.optimize import root
from time import time
import os
from matplotlib import *
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon



#____________________CONSTANTS____________________
l1=4.1      #shin length
l2=11.5     #foot length
l0=1.16     #thigh length
bodyr=6.66  #radius of main body (centre of pelvis to hip)
bodylen=sqrt(2)*bodyr   #length of body from hip to closest hip

#phase displacements for each leg
nOffs=[0,0.25,0.5,0.75]

H=l2-l1+1                 #pelvis height from z=0
SW=8                        #distance from foot tip path to hip axis perpendicular to said path [step width]
PL=8                        #length of foot tip path from start to end  [pace length]
PS=0                        #minimum distance to closest foot path  [pace seperation]
theta=-pi/2                 #direction of motion

timeMod=1       #speed of program
st=time()*10    #start time
zchange=0.1     #part of n reserved by raising/lowering foot, time taken from moving lowered foot
raisePhase=0.15 #part of n reserved for moving raised foot

colors=['red','green','blue','magenta'] #color of each leg


#Code that will only run with this variable True labelled [M] (unless otherwise clear)
allowMatplotlib = True

#Code that will only run with this variable True labelled [G] (unless otherwise clear)
allowGazebo = False


#if recieving math domain error for below function, modify this variable
BPinit = 9 


def equation(BP):
    return asin(l2*sin(acos((l1**2+l2**2-BP**2)/(2*l1*l2)))/BP)-pi/2

BP=round(root(lambda BP: equation(BP[0]),BPinit).x[0],3)    #threshhold of angle1 being accurate


def dirIndex(theta,index):      #change given index value for given direction (theta) [direction index]
    return (theta*2/pi+index)%4

#Calculate angles for knee and ankle
def formulaSideAngles(d,z,h):
    try:
        angle2=acos((l1**2+l2**2-(h-z)**2-d**2)/(2*l1*l2))  #knee angle
    
        angle1=asin(l2*sin(angle2)/sqrt((h-z)**2+d**2))+atan(d/(h-z))   #ankle angle
        if sqrt((h-z)**2+d**2)<BP:
            angle1=pi-angle1+2*atan(d/(h-z))

    except Exception as err:        #If failed to publish, error information will be displayed
        # print(err,'in formulaSideAngles(d,z,h)')
        angle1,angle2=-pi,-pi

    return angle1,angle2

#Calculate coords for side views of legs   [M]
def formulaSide(angle1,angle2,h):

    x1=(l2-l1)/5+l1*cos(angle1-pi/2)
    y1=h+l1*sin(angle1-pi/2)
    x2=x1+l2*cos(angle2+angle1+pi/2)
    y2=y1+l2*sin(angle2+angle1+pi/2)
    
    return x1,y1,x2,y2


#Determine next hip angle and foot to hip distance
def formulaBird(n,ymult,Bmod,z):
    
    y=(n*PL-(bodylen-PS)/2)*ymult #foot's y value in hip's reference frame

    new_modulus=sqrt(y**2+SW**2)    #foot to hip distance

    new_bearing=asin(y/new_modulus) #hip angle

    if z==H/4:
        print(n)
        end_y=((bodylen-PS)/2)*ymult
        new_modulus=(sqrt(y**2+SW**2)-sqrt(end_y**2+SW**2))*n+sqrt(end_y**2+SW**2)

    return new_bearing+Bmod, new_modulus-l0

#Draw triangle between 3 feet touching ground, if a foot is raised
def drawContactTriangle(legs,contact_triangle):
    if legs[0].z > 0:
        contact_triangle.set_data(
            [legs[1].foot_coords[0], legs[2].foot_coords[0], legs[3].foot_coords[0],legs[1].foot_coords[0]],
            [legs[1].foot_coords[1], legs[2].foot_coords[1], legs[3].foot_coords[1],legs[1].foot_coords[1]]
        )

    if legs[1].z > 0:
        contact_triangle.set_data(
            [legs[0].foot_coords[0], legs[2].foot_coords[0], legs[3].foot_coords[0],legs[0].foot_coords[0]],
            [legs[0].foot_coords[1], legs[2].foot_coords[1], legs[3].foot_coords[1],legs[0].foot_coords[1]]
        )

    if legs[2].z > 0:
        contact_triangle.set_data(
            [legs[0].foot_coords[0], legs[1].foot_coords[0], legs[3].foot_coords[0],legs[0].foot_coords[0]],
            [legs[0].foot_coords[1], legs[1].foot_coords[1], legs[3].foot_coords[1],legs[0].foot_coords[1]]
        )

    if legs[3].z > 0:
        contact_triangle.set_data(
            [legs[0].foot_coords[0], legs[1].foot_coords[0], legs[2].foot_coords[0],legs[0].foot_coords[0]],
            [legs[0].foot_coords[1], legs[1].foot_coords[1], legs[2].foot_coords[1],legs[0].foot_coords[1]]
        )


class leg():
    def __init__(self,index):

        self.index=index
        self.nOff=nOffs[index]  #Determine phase displacement of leg [n offset]

        self.Bmod=pi/4  #Determine the value to be added to bearing at return [bearing mod]
        self.ymult=1    #Determine multiplier for y value in formulaBird [y multiplier]
        self.invN=False #[Inverse n]

        #change variables based on direction to make robot walk in said direction
        if self.index==dirIndex(theta,1) or self.index==dirIndex(theta,3): self.Bmod-=pi/2
        if self.index==dirIndex(theta,0) or self.index==dirIndex(theta,2): self.ymult=-1
        if self.index==dirIndex(theta,1) or self.index==dirIndex(theta,2): self.invN=True

        self.z=0    #Determine starting foot height

        #   Calculate starting ankle and knee joint angles
        self.a1,self.a2=formulaSideAngles(l1-l2,self.z,H)

        if allowMatplotlib:

            self.foot_coords = [0,0]
            
            self.color=colors[index]

            #   Draw initial BE view leg
            self.line,=ax.plot([bodyr*sin(index*pi/2),(bodyr+l2-l1)*cos(index*pi/2)],
                               [bodyr*cos(index*pi/2),-(bodyr+l2-l1)*sin(index*pi/2)], self.color, linewidth=2)
        
            #   Calculate starting side view leg coords and draw
            x1,y1,x2,y2=formulaSide(self.a1,self.a2,H)
            
            #   Draw initial side view leg
            self.pelvis,= Sax.plot([-5-l0, (l2-l1)/5-l0], [H,H], 'grey', linewidth=2)
            self.thigh, = Sax.plot([(l2-l1)/5-l0, (l2-l1)/5], [H,H], self.color, linewidth=2)
            self.shin, = Sax.plot([(l2-l1)/5, x1], [H, y1], self.color, linewidth=2)
            self.foot, = Sax.plot([x1, x2], [y1, y2], self.color, linewidth=2)

    def update(self):
        self.n=(T/50+self.nOff)%1   

        self.z=0
        
        if self.n<1-raisePhase:
            if self.n<zchange/2:                        #Smooth foot lower
                self.z=H*0.4*(1-self.n/(zchange/2))

            elif self.n<1-raisePhase-zchange/2:                  #Lowered foot
                self.z=0

            else:                                           #Smooth foot raise
                self.z=(self.n-1+raisePhase+zchange/2)*(H/4)/(zchange/2)
                print(self.z)

            self.n=self.n/(1-raisePhase)
                
        else:                                           #Raised foot
            self.z=H*0.4
            self.n=(1-self.n)/raisePhase

        if self.invN:
            self.n=1-self.n

        #   Determine next bearing and hip to foot distance
        self.bearing, self.modulus=formulaBird(self.n,self.ymult,self.Bmod,self.z)
            
        #   Gen side view coords and move side view leg
        self.a1,self.a2=formulaSideAngles(self.modulus,self.z,H)
                  
        if allowMatplotlib:

            #   Calculate next BE view foot coordinates
            self.foot_coords = [(self.modulus+l0) * sin(self.bearing+pi*self.index/2)+bodyr*sin(self.index*pi/2), (self.modulus+l0) * cos(self.bearing+pi*self.index/2)+bodyr*cos(self.index*pi/2)]
            
            #   Calculate next side view leg coordinates
            x1,y1,x2,y2=formulaSide(self.a1,self.a2,H)

            #   Update BE leg view
            self.line.set_data([bodyr*sin(self.index*pi/2),
                                self.foot_coords[0]],
                               [bodyr*cos(self.index*pi/2),
                                self.foot_coords[1]])

            # #   Draw black dot on foot position in BE view
            # ax.plot((self.modulus+l0) * sin(self.bearing+pi*self.index/2)+bodyr*sin(self.index*pi/2),
            # (self.modulus+l0) * cos(self.bearing+pi*self.index/2)+bodyr*cos(self.index*pi/2) , 'o', markersize=3, color='black')

            #   Update side leg view
            self.pelvis.set_data([-8-l0, (l2-l1)/5-l0], [H,H])
            self.thigh, = Sax.plot([(l2-l1)/5-l0, (l2-l1)/5], [H,H], self.color, linewidth=2)
            self.shin.set_data([(l2-l1)/5, x1], [H, y1])
            self.foot.set_data([x1, x2], [y1, y2])

            #   If raised, change color to yellow
            if self.z>0:
                self.line.set_color('yellow')
            else:
                self.line.set_color(self.color)

        if allowGazebo:
            self.publishAngles()

    def publishAngles(self):   #Publish joint angles to gazebo topic so that simulation update joint to said angles [G]
        try:
            os.system('gz topic -p /gazebo/servos/test\ plugin/target_angle -m x:'+str(self.index*3)+',y:'+str(self.bearing-pi/4)+
                      '&& gz topic -p /gazebo/servos/test\ plugin/target_angle -m x:'+str(self.index*3+1)+',y:'+str(self.a1-pi/2)+
                      '&& gz topic -p /gazebo/servos/test\ plugin/target_angle -m x:'+str(self.index*3+2)+',y:'+str(self.a2-pi/4))
            print("published")
        except Exception as err:        #If failed to publish, error information will be displayed
            raise Exception(err)


if allowMatplotlib:
    #___________________________bird's eye figure and axis___________________________
    
    #   BE view axis setup
    fig, ax = plt.subplots()
    ax.set_xlim([-(l1+l2+l0+bodyr), l1+l2+l0+bodyr])
    ax.set_ylim([-(l1+l2+l0+bodyr), l1+l2+l0+bodyr])
    ax.set_aspect('equal')

    #   Draw pelvis
    ax.plot([bodyr,-bodyr],[0,0],'black')
    ax.plot([0,0],[bodyr,-bodyr],'black')
    ax.add_patch(Polygon([[bodyr,0],[0,bodyr],[-bodyr,0],[0,-bodyr]],
                         closed=True,color='grey'))

    #   Hide the axis lines
    ax.spines['left'].set_color('none')
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')
    ax.spines['bottom'].set_color('none')
    ax.tick_params(axis='both', which='both', bottom=False, top=False, labelbottom=False, right=False, left=False, labelleft=False)


    # #___________________________side view figure and axis___________________________

    #   Side view axis setup
    Sfig,Sax=plt.subplots()
    Sax.set_xlim([-l0, (l2+l1)])
    Sax.set_ylim([0, (l2)*1.1])
    Sax.set_aspect('equal')

    #   Side view draw floor
    floor, = Sax.plot([-10,10],[0,0], 'black', linewidth=2)

#init legs
legs=[leg(0),leg(1),leg(2),leg(3)]

#create contact triangle
contact_triangle, = ax.plot([0,0,0,0],[0,0,0,0],'orange')

#main loop
while True:
    T=(time()*10-st)*timeMod%50
    if T > 50:
        T = 0
        st = time()*10
    for i in legs:
        i.update()
    drawContactTriangle(legs, contact_triangle)
    plt.pause(0.00001)
