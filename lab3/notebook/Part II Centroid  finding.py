
# coding: utf-8

## Tracking Asteriod Daphane

# In[58]:
import numpy as np
# get_ipython().magic(u'pylab inline')


# In[59]:

import matplotlib.pyplot as plt, mpld3
# mpld3.enable_notebook()


# In[60]:
import os
os.chdir("../data/2013Data/101613/")
os.getcwd()

# In[61]:

datafile = 'd112.fits'
flat1 = 'd106.fits'


## Data Reduction Bias and Flats

# In[62]:

import astropy.io.fits as pf
import heapq


# In[63]:

def bsub (fits1,plot=False):
    x =pf.getdata(fits1)
    average = x
    n=0
    #Bias Subtraction
    for i in np.arange(100,105):
        dark = 'd{}.fits'.format(i)
        hdulist = pf.open(fits1)
        average = average + hdulist[0].data
        n = n + 1
    xb = average / float(n)
    if (plot):
        plt.title("After Bias Subtraction")
        imp = plt.imshow(xb[600:300:-1],cmap='gray_r',vmin=1000,vmax=2000)#,norm=LogNorm()) 
        plt.colorbar()
        plt.figure()
        plt.title("Before Bias Subtraction")
        imp = plt.imshow(x[600:300:-1],cmap='gray_r',vmin=700,vmax=1500)#,norm=LogNorm()) 
        plt.colorbar()
    return xb


# In[64]:

from matplotlib.colors import LogNorm
xb = bsub(datafile)
flat = pf.getdata(flat1)
flatb = bsub(flat1) # Bias subtract flatb = flatb/np.median(flatb) # normalize
fig = plt.figure()
ax1 = fig.add_subplot(211)
processed =(xb/flatb)* np.median(flatb)
plt.title("Without normalization", fontsize=12)
imp = plt.imshow((xb/flatb)[600:300:-1],cmap='gray_r',vmin=0,vmax=0.1)#norm=LogNorm()) 
plt.colorbar()
ax2 = fig.add_subplot(212)
plt.title("With normalization", fontsize=12)
imp = plt.imshow(processed[600:300:-1],cmap='gray_r',vmin=1000,vmax=2000)
plt.colorbar()
plt.tight_layout()


## Centroid-Finding on this Image 

# Ok this matches  20:28:28.31 -13:34:34.2 on the Aladin Sky Atlas, so this is correctly displayed north up

# In[65]:

a = plt.hist(processed.flatten(), 1000)


# In[66]:

# imp = plt.imshow(processed,cmap='gray_r',vmin=1000,vmax=2000,origin="lower")
# plt.colorbar()
# plt.figure()
# imp = plt.imshow(processed.T,cmap='gray_r',vmin=1000,vmax=2000,origin="lower")
# plt.colorbar()
# plt.figure()
imp = plt.imshow(processed[::-1],cmap='gray_r',vmin=0,vmax=8000,origin="lower")
plt.axis("tight")
plt.colorbar()


# In[67]:

def find_centroid(data,threshold=40000,plot =False):
    #Don't modify original data
    x = np.copy(data)
    x = x[::-1]
    # we actually find that most of this stuff can be cut away by simply eliminating the stuff in overscan column (dark pixels)
    x = x[:,:1024,] #Truncate the overscan columns
    #See if we can use maximum values to automatically determine the threshold?
    # print heapq.nlargest(5,x.flatten())
    #Preprocessing: Masking Bad Columns
#     x[100:,]=0
#     Bad Column masking
    x[:,1023:]=0
    x[1023:,]=0
    x[:,256]=0
    x[255]=0
    x[:,1002]=0
    x[0,:]=0
    #Four Corner Masking
    x[950:,0:60,]=0
    x[0:50,0:50,]=0
    x[0:15,0:18,]=0
    x[1016:,1016:,]=0
    # **** Remember that your data is flipped around when you are plotting!! So you also need to do this when you are finding these values, or else these values are flipped. (This took 3 hours to debug!!)
    coord_lst =  np.array(np.where(x>threshold))
    plt.figure()
    if (plot):
        imp = plt.imshow(x,cmap='gray_r',vmin=0,vmax=8000,origin="lower") 
        plt.plot(coord_lst[1],coord_lst[0],'*',color='r')
        plt.colorbar()
        plt.axis("tight")
    # Finding centroid
    star_data = [] #each star is a group of points, used later for centroid calculation
    star_lst  = np.array([[],[],[]]) #each star is denoted by (x,y,brightness)
    ################################################
    # METHOD : Sum all values above threshold     #
    ################################################
    #Look around a certain radius the coordinate specified by the list
#     r = 3 #3 pixels for now, sum the intensity profile around the star on a 3-by-3 box
    star_data = []
    for i in np.arange(coord_lst.shape[1]):
    #     print a[0][i],a[1][i]
        tmp = []
        xi = coord_lst[0][i]
        yi = coord_lst[1][i]
        for n in np.arange(coord_lst.shape[1]): #looping through other points
            nx = coord_lst[0][n]
            ny = coord_lst[1][n]
    #         print nx, ny
            r_i = np.sqrt((nx-xi)**2+(ny-yi)**2)
            if (r_i<20 and r_i!=0):
    #             print r_i
                tmp.append([nx,ny])#for clustering only, storing x,y locations
        star_data.append(tmp)
    prev_length = 0 #no length zero, enforce first to store
    cleaned_star_data=[]
    for i in star_data:
        if prev_length!=len(i):
            cleaned_star_data.append(i)        
        prev_length  = len(i)
    unique_len = set([len(i) for i in star_data ])
#     assert len(unique_len)==len(cleaned_star_data)
    return cleaned_star_data


# In[68]:

processed[::-1][530,440]


# In[69]:

star_lst = find_centroid(processed,4000,True) #Don't go below 1000


# $x = \frac{\sum_{i} x_i I_i}{\sum_{i} I_i}$

# In[70]:

def centroid_calc(star_locs):
    #given a list of grouped points, compute the centroid pixel position (round up to center pixel (int))
#     print star_locs
    x_lst = np.array(star_locs)[:,0]
    y_lst = np.array(star_locs)[:,1]
    I=[]
    for i in np.arange(len(x_lst)):
        I.append(processed[x_lst[i],y_lst[i]])
    I = np.array(I)
    x_cm = sum(x_lst*I)/sum(I)
    y_cm = sum(y_lst*I)/sum(I)
    #return x_cm, y_cm and intensity of star (summed over star)
    return [x_cm,y_cm,sum(I)]


# In[71]:

centroid_calc(star_lst[0])


# In[72]:

len(star_lst)


# In[73]:

x_lst=[]
y_lst=[]
img_data_stars = [[],[]]


# In[74]:

star_lst= [i for i in star_lst if len(i)>0]
len(star_lst)


# 19:34:43.52 -10:04:01.6

# In[75]:

star_info=[]
img_data_stars = [[],[]]
for i in star_lst:
    centroid = centroid_calc(i)
    star_info.append(centroid)
    imp = plt.imshow(processed[::-1],cmap='gray_r',vmin=0,vmax=8000,origin="lower")
    img_data_stars[0].append(centroid[1])
    img_data_stars[1].append(centroid[0])
plt.plot(img_data_stars[0],img_data_stars[1],'*',color="red") #Note the switch    
plt.axis("tight")
plt.show()
#Visually verfied with above pattern


## Matching with USNO Catalog

# subtract 00:00:30 from RA and 00:12:00 from Dec from the numbers you get from fits header

# In[76]:

# import os
# os.chdir("data/2013Data/101613/")


# In[77]:

import astropy.io.fits as pf 
import matplotlib.pyplot as plt 
import numpy as np
import urllib as url
s1 = pf.open(datafile)
# Read position from the FITS file and convert RA/DEC to degrees # be sure to check that the header data is reliable. If not
# edit the position by hand.
ras = s1[0].header['ra']
des = s1[0].header['dec']
radeg = 15*(float(ras[0:2]) + float(ras[3:5])/60. + float(ras[6:])/3600.)
dsgn = np.sign(float(des[0:3]))
dedeg = float(des[0:3]) + dsgn*float(des[4:6])/60. + dsgn*float(des[7:])/3600.
fovam = 3.0 # size of square search field in arc min
#chdir to where usno.py lies
import os 
os.chdir("../../../")
from usno import usno
# Need to put in values for epoch but in what form???
#EPOCH FOR POCO POSITION IS CURRENT DATE
epoch = s1[0].header['EQUINOXU']
name,rad,ded,rmag = usno(radeg,dedeg,fovam,epoch) #Epoch data taken in year of 2013
w = np.where(rmag < 17.)[0] # select only bright stars r < 15 mag.
plt.plot(rad[w],ded[w],'g.') 
plt.locator_params(axis='x',nbins=4) 
plt.locator_params(axis='y',nbins=4) 
plt.tick_params('x',pad=10)
plt.xlabel('RA [Deg]')
plt.ylabel('Dec [Deg]') 
plt.ticklabel_format(useOffset=False)
# plt.axis('scaled')
# plt.xlim([266.11,266.03]) # reverse the x-axis direction
os.chdir("data/2013Data/101613/")


##### Find in Aladin Sky Atlas 

# 19:34:43.19 -10:04:06.9

# In[78]:

#subtract 00:00:30 from RA
ras


# In[79]:

#subtract 00:12:00 from Dec 
des


# In[80]:

usno_star_lst =  array([rad[w],ded[w]]) #returned star list from USNO 


# In[81]:

fig  = plt.figure()
ax1 = fig.add_subplot(121)
plt.plot(usno_star_lst[0],usno_star_lst[1],'g.') 
plt.locator_params(axis='x',nbins=4) 
plt.locator_params(axis='y',nbins=4) 
plt.tick_params('x',pad=10)
plt.xlabel('RA [Deg]')
plt.ylabel('Dec [Deg]') 
plt.ticklabel_format(useOffset=False)
plt.axis('scaled')
ax2 = fig.add_subplot(122)
plt.plot(img_data_stars[0],img_data_stars[1],'o')
plt.xlabel('X [pixel]')
plt.ylabel('Y [pixel]') 
plt.axis('scaled')
plt.tight_layout()


### Tangent plane transformation

# In[82]:

a = usno_star_lst[0] #alpha : list of usno ra [deg]
d = usno_star_lst[1]  #delta : list of usno dec [deg]
a_0 = radeg #alpha_0 : center of CCD image (ra)
d_0 = dedeg #delta_0 : center of CCD image (dec)

#convert degrees to radians 
d = d*np.pi/180.
d_0 = d_0*np.pi/180.
a = a*np.pi/180.
a_0 = a_0*np.pi/180.
#tangent plane projection (get image position on CCD  (X,Y) )
X_stdcoord = -(cos(d)*sin(a-a_0))/(cos(d_0)*cos(d)*cos(a-a_0)+sin(d)*sin(d_0))
Y_stdcoord = -(sin(d_0)*cos(d)*cos(a-a_0)-cos(d_0)*sin(d))/(cos(d_0)*cos(d)*cos(a-a_0)+sin(d)*sin(d_0))


# In[83]:

X_stdcoord #looks right, CCD should be small.


#### Pixel (X,Y) to celestial coord (x,y) transformation

# In[84]:

#For Nickel (assuming its an "ideal" camera)
f = 16840 #[mm]
p= 0.015 #[mm]
x_0 = 512#center of pixel coord (basically 1024/2)
y_0 = 512
x_pixcoord = f*(X_stdcoord/p)+x_0
y_pixcoord=f*(Y_stdcoord/p)+y_0


# In[85]:

x_pixcoord #looks right, this is now in pixel coordinates!


### Why is USNO data flipped around after tangent plane projection?

# In[86]:

from matplotlib.legend_handler import HandlerLine2D
plt.plot(x_pixcoord,y_pixcoord,'o',label= "USNO")
plt.plot(img_data_stars[0],img_data_stars[1],'*',color="red",label="CCD data")
plt.legend(loc='lower right',prop={'size':12},numpoints=1)
plt.xlabel("x [pixel]",fontsize=12)
plt.ylabel("y [pixel]",fontsize=12)


# In[87]:

plt.plot(x_pixcoord , "*" , label="x")
plt.plot(y_pixcoord,".", color="red", label= "y")
plt.xlabel("x or y [pixel]",fontsize=12)
plt.ylabel("Pixel offset",fontsize=12)
plt.legend(loc='lower right',prop={'size':12},numpoints=1)


## Least Squares

# Do not confuse x_pixcoord with x, x_pixcoord and X_stdcoord is USNO coordinate;
# The goal here is to convert CCD imaging coordinate to RA DEC

# In[89]:

a_x= np.matrix(np.c_[x_pixcoord])


# In[90]:

B = np.matrix(np.c_[(f/p)*X_stdcoord,(f/p)*Y_stdcoord,np.ones_like(Y_stdcoord)])


# In[91]:

c_x = np.linalg.inv(B.T*B)*B.T*a_x


# In[92]:

a_y= np.matrix(np.c_[y_pixcoord])


# In[93]:

c_y = np.linalg.inv(B.T*B)*B.T*a_y


# In[94]:

T = np.matrix([[(f/p)*c_x.item(0),(f/p)*c_x.item(1),c_x.item(2)],[(f/p)*c_y.item(0),(f/p)*c_y.item(1),c_y.item(2)],[0,0,1]])


## Why is f/p different from worksheet

# In[102]:

f/p


# In[103]:

sqrt(np.linalg.det(T))


# Since X = (x,y,1)

# In[104]:

x_pixcoord = np.matrix(np.c_[img_data_stars[0]])
y_pixcoord = np.matrix(np.c_[img_data_stars[1]])


# In[105]:

x = np.matrix((img_data_stars[0],img_data_stars[1],np.ones_like(img_data_stars[0])))


# In[106]:

X = np.linalg.inv(T)*x


# In[107]:

X_std_coord_final = X[0]


# In[108]:

Y_std_coord_final = X[1]


# In[109]:

X_std_coord_final = np.array(X_std_coord_final)
Y_std_coord_final = np.array(Y_std_coord_final)


# Then convert this into celestial coordinates using the tangent plane equation 
# 
# (X,Y)$\rightarrow$ ($\alpha$,$\delta$)
# 
# and keeping the same a_0,d_0

# In[110]:

dedeg


# In[111]:

a_0 = radeg
d_0 = dedeg


# In[112]:

a = np.arctan(-X_std_coord_final/(cos(d_0)-Y_std_coord_final*sin(d_0)))+a_0


# In[113]:

d = np.arcsin((sin(d_0)+Y_std_coord_final*cos(d_0))/((1+np.array(X_std_coord_final)**2+np.array(Y_std_coord_final)**2)**0.5))


# In[114]:

a


# In[115]:

d_0


# In[116]:

d


# In[124]:

fig  = plt.figure()
plt.plot(usno_star_lst[0],usno_star_lst[1],'g.');
plt.locator_params(axis='x',nbins=4) 
plt.locator_params(axis='y',nbins=4) 
plt.tick_params('x',pad=10)
plt.xlabel('RA [Deg]')
plt.ylabel('Dec [Deg]') 
plt.ticklabel_format(useOffset=False)
plt.axis('scaled')
# plt.plot(293.8,-9.9,'o',color="red")
# plt.plot(a,-9.9,'o',color="red");
plt.plot(a,d,'o',color="red");


# In[53]:

#Comparing the x values together first
#How do you know which star matches which, and subtract for residual
# a-usno_star_lst


# In[ ]:




# In[ ]:



