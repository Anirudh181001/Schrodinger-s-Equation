import numpy as np
import scipy.constants as c

def deg_to_rad(deg):
    rad=deg*c.pi/180
    return rad
  
def rad_to_deg(rad):
    deg=rad*180/c.pi
    return deg
def spherical_to_cartesian(r,theta,phi):
    x=round(r*np.sin(theta)*np.cos(phi),5)
    y=round(r*np.sin(theta)*np.sin(phi),5)
    z=round(r*np.cos(theta),5)
    return (x,y,z)

def cartesian_to_spherical(x, y, z):
    r=round(((x**2)+(y**2)+(z**2))**0.5,5)
    if(r==0):
        theta=1.00000
    else:
        theta=round(np.arccos(z/r),5)
        if(y==0):
            phi=0.00000
        elif(x==0):
            if(y>=0):
                phi=round(c.pi/2,5)
            else:
                phi=round(-(c.pi/2),5)
        else:
            phi=round(np.arctan(y/x),5)
    return (r,theta,phi)
def angular_wave_func(m,l,theta,phi):
    if l==0:
        y=(1/(4*c.pi))**0.5*np.exp(phi*0j)
    elif l==1:
        if m==0:
            y=(3/(4*c.pi))**0.5*np.cos(theta)*np.exp(phi*0j)
        elif m==1:
            y=-1*(3/(8*c.pi))**0.5*np.sin(theta)*np.exp(phi*1j)
        else:
            y=1*(3/(8*c.pi))**0.5*np.sin(theta)*np.exp(phi*-1j)
    elif l==2:
        if m==0:
            y=(5/(16*c.pi))**0.5*(3*(np.cos(theta))**2-1)*np.exp(phi*0j)
        elif m==1:
            y=-(15/(8*c.pi))**0.5*np.cos(theta)*np.sin(theta)*np.exp(phi*1j)
        elif m==-1:
            y=(15/(8*c.pi))**0.5*np.cos(theta)*np.sin(theta)*np.exp(-1j*phi)
        elif m==2:
            y=(15/(32*c.pi))**0.5*(np.sin(theta))**2*np.exp(2j*phi)
        else:
            y=(15/(32*c.pi))**0.5*(np.sin(theta))**2*np.exp(-2j*phi)
    elif l==3:
        if m==3:
            y=-(35/(64*c.pi))**0.5*(np.sin(theta))**3*np.exp(3j*phi)
        elif m==-3:
            y=(35/(64*c.pi))**0.5*(np.sin(theta))**3*np.exp(-3j*phi)
        elif m==2:
            y=(105/(32*c.pi))**0.5*np.cos(theta)*(np.sin(theta))**2*np.exp(2j*phi)
        elif m==-2:
            y=(105/(32*c.pi))**0.5*np.cos(theta)*(np.sin(theta))**2*np.exp(-2j*phi)
        elif m==1:
            y=-(21/(64*c.pi))**0.5*np.sin(theta)*(5*(np.cos(theta))**2-1)*np.exp(1j*phi)
        elif m==-1:
            y=-(21/(64*c.pi))**0.5*np.sin(theta)*(5*(np.cos(theta))**2-1)*np.exp(-1j*phi)
        else:
            y=(7/(16*c.pi))**0.5*(5*(np.cos(theta))**3-3*(np.cos(theta)))*np.exp(0j*phi)
            
    return np.around(y,5)
# Angular Wavefunction sympy code
from sympy import *
import math
import numpy as np

l, x, m = symbols('l x m')
init_printing(use_unicode = True)

def legendre_polynomial(l):
# Takes in l as an integer and returns a expression
    expr = (1 / (2**l*math.factorial(l))) * diff((x**2 - 1)**  l, x, l)
    return expr

def associated_legendre_function(m, l, expr, theta):
    second_expr = (1 - x**  2) ** (abs(m) / 2) * diff(expr, x, abs(m))
    return second_expr.subs(x, np.cos(theta))


def angular_wave_func_poly(m, l, theta, phi):
    epsilon = 0
    
    polynomial = legendre_polynomial(l)
    
    function_value = associated_legendre_function(m, l, polynomial, theta)
    
    
    if m > 0:
        epsilon = (-1)**  m
    else:
        epsilon = 1
    
    solution = epsilon * np.sqrt(((2 * l + 1) / (4 * np.pi)) * (math.factorial(l - abs(m)) / math.factorial(l + abs(m))))
    solution *= np.exp(complex(0, m * phi)) * function_value
    return np.around(np.complex(solution), 5)
a=c.physical_constants['Bohr radius'][0]
def radial_wave_func(n,l,r):
    if n==1:
        z=2*(a**(-3/2))*np.exp(-r/a)
    elif n==2:
        if l==0:
            z=(1/(2**0.5))*(a**(-3/2))*(1-r/(2*a))*np.exp((-1*r)/(2*a))
        else:
            z=(1/(24**0.5))*(a**(-3/2))*(r/a)*np.exp((-1*r)/(2*a))
    elif n==3:
        if l==0:
            z=(2/(81*3**0.5))*(a**(-3/2))*(27-18*(r/a)+2*(r/a)**2)*np.exp((-1*r)/(3*a))
        elif l==1:
            z=(8/(27*6**0.5))*(a**(-3/2))*(1-(r/(6*a)))*(r/a)*np.exp((-r)/(3*a))
        else:
            z=(4/(81*(30**0.5)))*(a**(-3/2))*((r/a)**2)*np.exp((-1*r)/(3*a))
    else:
        if l==0:
            z=(1/4)*a**(-3/2)*(1-(3/4)*(r/a)+(1/8)*(r/a)**2-(1/192)*(r/a)**3)*np.exp((-1*r)/(4*a))
        elif l==1:
            z=((5**0.5)/(16*3**0.5))*(a**(-3/2))*(r/a)*(1-(1/4)*(r/a)+(1/80)*(r/a)**2)*np.exp((-1*r)/(4*a))
        elif l==2:
            z=(1/(64*5**0.5))*(a**(-3/2))*((r/a)**2)*(1-(1/12)*(r/a))*np.exp((-1*r)/(4*a))
        else:
            z=(1/(768*35**0.5))*(a**(-3/2))*((r/a)**3)*np.exp((-1*r)/(4*a))
            
    return np.around(z/(a**(-3/2)),5)
def mgrid2d(xstart, xend, xpoints, ystart, yend, ypoints):
    # initialize a list to store the grid points that will be returned
    xr=[xstart]
    yr=[ystart]
    xD=(xend-xstart)/(xpoints-1)
    yD=(yend-ystart)/(ypoints-1)
    for i in range(xpoints - 1):
        xval = xr[i] + xD
        xr.append(xval)
    for i in range(ypoints - 1):
        yval = yr[i] + yD
        yr.append(yval)
    yl = []
    xl = []
    
    for i in xr:
        k = []
        for j in range(len(yr)):
            k.append(i)
        xl.append(k)
            
    for i in range(len(xr)):
        yl.append(yr)
    
    return(xl,yl)
   
# Test:
# assert statement will throw error if the result is wrong
# no output will be produced for correct results

assert np.shape(mgrid2d(0, 4, 5, 0, 4, 5)) == np.shape(np.mgrid[0:4:5j,0:4:5j])
assert np.allclose(mgrid2d(0, 4, 5, 0, 4, 5), np.mgrid[0:4:5j,0:4:5j])

assert np.shape(mgrid2d(0, 5, 15, 0, 4, 10)) == np.shape(np.mgrid[0:5:15j,0:4:10j])
assert np.allclose(mgrid2d(0, 5, 15, 0, 4, 10), np.mgrid[0:5:15j,0:4:10j])
def mgrid3d(xstart, xend, xpoints, 
            ystart, yend, ypoints, 
            zstart, zend, zpoints):
    xr = []
    yr = []
    zr = []
    xval = xstart
    xcount = 0

    xsike = abs(xend- xstart)/(xpoints-1)
    ysike = abs(yend- ystart)/(ypoints-1)
    zsike = abs(zend- zstart)/(zpoints-1)

    while xcount < xpoints:
        yval = ystart
        ycount = 0
        xops = []
        yops = []
        zops = []
    
        while ycount < ypoints:
            zcount = 0
            zval = zstart
            xindie = []
            yindie = []
            zindie = []

            while zcount < zpoints:
                xindie.append(xval)
                yindie.append(yval)
                zindie.append(zval)
                zval += zsike
                zcount += 1

            xops.append(xindie)
            yops.append(yindie)
            zops.append(zindie)
            yval += ysike
            ycount += 1

        xr.append(xops)
        yr.append(yops)
        zr.append(zops)
        xval += xsike
        xcount += 1
        
    return xr, yr, zr

            


# Test:
# assert statement will throw error if the result is wrong
# no output will be produced for correct results

assert np.shape(mgrid3d(0, 4, 5, 0, 4, 5, 0, 4, 5)) == np.shape(np.mgrid[0:4:5j,0:4:5j,0:4:5j])
assert np.allclose(mgrid3d(0, 4, 5, 0, 4, 5, 0, 4, 5), np.mgrid[0:4:5j,0:4:5j,0:4:5j])

assert np.shape(mgrid3d(0, 5, 15, 0, 4, 10, 1, 2, 3)) == np.shape(np.mgrid[0:5:15j,0:4:10j,1:2:3j])
assert np.allclose(mgrid3d(0, 5, 15, 0, 4, 10, 1, 2, 3), np.mgrid[0:5:15j,0:4:10j,1:2:3j])
###Code:
def hydrogen_wave_func(n,l,m,roa,Nx,Ny,Nz):
    vmgrid3d = np.vectorize(mgrid3d)
    vcartesian_to_spherical = np.vectorize(cartesian_to_spherical)
    vangular_wave_func = np.vectorize(angular_wave_func)
    
    xx, yy, zz = np.around(vmgrid3d(-roa, roa, Nx, -roa, roa, Ny, -roa, roa, Nz), 5)

    # Vectorized, r, phi and theta
    r, theta, phi = vcartesian_to_spherical(xx, yy, zz) # phi and theta in radians
                
    y_lookm = None
    if m < 0:
        y_lookm = 1j * (1 / (2 ** 0.5)) * (vangular_wave_func(m, l, theta, phi) + (-(-1) ** m) * vangular_wave_func(-m, l, theta, phi))
    elif m > 0:
        y_lookm = ((1 / (2 ** 0.5))) * (vangular_wave_func(-m, l, theta, phi) + (-1) ** m * vangular_wave_func(m, l, theta, phi))
    else:
        y_lookm = (vangular_wave_func(0, l, theta, phi))
        
    radial_solution = radial_wave_func(n, l, r*a)
    
    real_solution = np.around(np.abs(y_lookm * radial_solution) ** 2, 5)
    
    return xx, yy, zz, real_solution

###Test:
print('Test 1')
x,y,z,mag = hydrogen_wave_func(2 ,1 ,1 ,8 ,2 ,2 ,2)
print('x, y, z:')
print(x, y, z)
print('mag:')
print(mag)

print('\n')
print('Test 2')
x,y,z,mag = hydrogen_wave_func(2 ,1 ,1 ,5 ,3 ,4 ,2)
print('x, y, z:')
print(x, y, z)
print('mag:')
print(mag)

print('\n')
print('Test 3')
x,y,z,mag = hydrogen_wave_func(2 ,0 ,0 ,3 ,5 ,4 ,3)
print('x, y, z:')
print(x, y, z)
print('mag:')
print(mag)
# Code to save the data to a file so that 
# you don't have to keep on computing it:
import numpy as np
import scipy.constants as c
import math

a = c.physical_constants['Bohr radius'][0]

import scipy.constants as c
import numpy as np
def spherical_to_cartesian(r,theta,phi):
    x=round(r*np.sin(theta)*np.cos(phi),5)
    y=round(r*np.sin(theta)*np.sin(phi),5)
    z=round(r*np.cos(theta),5)
    return (x,y,z)

def cartesian_to_spherical(x, y, z):
    r=round(((x**2)+(y**2)+(z**2))**0.5,5)
    if(r==0):
        theta=1.00000
    else:
        theta=round(np.arccos(z/r),5)
        if(y==0):
            phi=0.00000
        elif(x==0):
            if(y>=0):
                phi=round(c.pi/2,5)
            else:
                phi=round(-(c.pi/2),5)
        else:
            phi=round(np.arctan(y/x),5)
    return (r,theta,phi)



def angular_wave_func(m,l,theta,phi):
    if l==0:
        y=(1/(4*c.pi))**0.5*np.exp(phi*0j)
    elif l==1:
        if m==0:
            y=(3/(4*c.pi))**0.5*np.cos(theta)*np.exp(phi*0j)
        elif m==1:
            y=-1*(3/(8*c.pi))**0.5*np.sin(theta)*np.exp(phi*1j)
        else:
            y=1*(3/(8*c.pi))**0.5*np.sin(theta)*np.exp(phi*-1j)
    elif l==2:
        if m==0:
            y=(5/(16*c.pi))**0.5*(3*(np.cos(theta))**2-1)*np.exp(phi*0j)
        elif m==1:
            y=-(15/(8*c.pi))**0.5*np.cos(theta)*np.sin(theta)*np.exp(phi*1j)
        elif m==-1:
            y=(15/(8*c.pi))**0.5*np.cos(theta)*np.sin(theta)*np.exp(-1j*phi)
        elif m==2:
            y=(15/(32*c.pi))**0.5*(np.sin(theta))**2*np.exp(2j*phi)
        else:
            y=(15/(32*c.pi))**0.5*(np.sin(theta))**2*np.exp(-2j*phi)
    elif l==3:
        if m==3:
            y=-(35/(64*c.pi))**0.5*(np.sin(theta))**3*np.exp(3j*phi)
        elif m==-3:
            y=(35/(64*c.pi))**0.5*(np.sin(theta))**3*np.exp(-3j*phi)
        elif m==2:
            y=(105/(32*c.pi))**0.5*np.cos(theta)*(np.sin(theta))**2*np.exp(2j*phi)
        elif m==-2:
            y=(105/(32*c.pi))**0.5*np.cos(theta)*(np.sin(theta))**2*np.exp(-2j*phi)
        elif m==1:
            y=-(21/(64*c.pi))**0.5*np.sin(theta)*(5*(np.cos(theta))**2-1)*np.exp(1j*phi)
        elif m==-1:
            y=-(21/(64*c.pi))**0.5*np.sin(theta)*(5*(np.cos(theta))**2-1)*np.exp(-1j*phi)
        else:
            y=(7/(16*c.pi))**0.5*(5*(np.cos(theta))**3-3*(np.cos(theta)))*np.exp(0j*phi)
            
    return np.around(y,5)

def radial_wave_func(n,l,r):
    if n==1:
        z=2*(a**(-3/2))*np.exp(-r/a)
    elif n==2:
        if l==0:
            z=(1/(2**0.5))*(a**(-3/2))*(1-r/(2*a))*np.exp((-1*r)/(2*a))
        else:
            z=(1/(24**0.5))*(a**(-3/2))*(r/a)*np.exp((-1*r)/(2*a))
    elif n==3:
        if l==0:
            z=(2/(81*3**0.5))*(a**(-3/2))*(27-18*(r/a)+2*(r/a)**2)*np.exp((-1*r)/(3*a))
        elif l==1:
            z=(8/(27*6**0.5))*(a**(-3/2))*(1-(r/(6*a)))*(r/a)*np.exp((-r)/(3*a))
        else:
            z=(4/(81*(30**0.5)))*(a**(-3/2))*((r/a)**2)*np.exp((-1*r)/(3*a))
    else:
        if l==0:
            z=(1/4)*a**(-3/2)*(1-(3/4)*(r/a)+(1/8)*(r/a)**2-(1/192)*(r/a)**3)*np.exp((-1*r)/(4*a))
        elif l==1:
            z=((5**0.5)/(16*3**0.5))*(a**(-3/2))*(r/a)*(1-(1/4)*(r/a)+(1/80)*(r/a)**2)*np.exp((-1*r)/(4*a))
        elif l==2:
            z=(1/(64*5**0.5))*(a**(-3/2))*((r/a)**2)*(1-(1/12)*(r/a))*np.exp((-1*r)/(4*a))
        else:
            z=(1/(768*35**0.5))*(a**(-3/2))*((r/a)**3)*np.exp((-1*r)/(4*a))
            
    return np.around(z/(a**(-3/2)),5)

def mgrid3d(xstart, xend, xpoints, 
            ystart, yend, ypoints, 
            zstart, zend, zpoints):
    xr = []
    yr = []
    zr = []
    xval = xstart
    xcount = 0

    xsike = abs(xend- xstart)/(xpoints-1)
    ysike = abs(yend- ystart)/(ypoints-1)
    zsike = abs(zend- zstart)/(zpoints-1)

    while xcount < xpoints:
        yval = ystart
        ycount = 0
        xops = []
        yops = []
        zops = []
    
        while ycount < ypoints:
            zcount = 0
            zval = zstart
            xindie = []
            yindie = []
            zindie = []

            while zcount < zpoints:
                xindie.append(xval)
                yindie.append(yval)
                zindie.append(zval)
                zval += zsike
                zcount += 1

            xops.append(xindie)
            yops.append(yindie)
            zops.append(zindie)
            yval += ysike
            ycount += 1

        xr.append(xops)
        yr.append(yops)
        zr.append(zops)
        xval += xsike
        xcount += 1
        
    return xr, yr, zr



def hydrogen_wave_func(n,l,m,roa,Nx,Ny,Nz):
    vmgrid3d = np.vectorize(mgrid3d)
    vcartesian_to_spherical = np.vectorize(cartesian_to_spherical)
    vangular_wave_func = np.vectorize(angular_wave_func)
    
    xx, yy, zz = np.around(vmgrid3d(-roa, roa, Nx, -roa, roa, Ny, -roa, roa, Nz), 5)

    # Vectorized, r, phi and theta
    r, theta, phi = vcartesian_to_spherical(xx, yy, zz) # phi and theta in radians
                
    y_lookm = None
    if m < 0:
        y_lookm = 1j * (1 / (2 ** 0.5)) * (vangular_wave_func(m, l, theta, phi) + (-(-1) ** m) * vangular_wave_func(-m, l, theta, phi))
    elif m > 0:
        y_lookm = ((1 / (2 ** 0.5))) * (vangular_wave_func(-m, l, theta, phi) + (-1) ** m * vangular_wave_func(m, l, theta, phi))
    else:
        y_lookm = (vangular_wave_func(0, l, theta, phi))
        
    radial_solution = radial_wave_func(n, l, r*a)
    
    real_solution = np.around(np.abs(y_lookm * radial_solution) ** 2, 5)
    
    return xx, yy, zz, real_solution


print('Test ')
x,y,z,mag=hydrogen_wave_func(4,3,3,40,100,100,100)
print('x, y, z:')
print(x, y, z)
print('mag:')
print(mag)
print (x,y,z,mag)
x.dump('x_test.dat')
y.dump('y_test.dat')
z.dump('z_test.dat')
mag.dump('den_test.dat')
# Mayavi code:

from mayavi import mlab

mu, sigma = 0, 0.1 
x = np.load('x_test.dat')
y = np.load('y_test.dat')
z = np.load('z_test.dat')

density = np.load('den_test.dat')
figure = mlab.figure('DensityPlot')
pts = mlab.contour3d(density,contours=40,opacity=0.4)
mlab.axes()
mlab.show()
###Volume slicer code:
import numpy as np

from traits.api import HasTraits, Instance, Array, \
    on_trait_change
from traitsui.api import View, Item, HGroup, Group

from tvtk.api import tvtk
from tvtk.pyface.scene import Scene

from mayavi import mlab
from mayavi.core.api import PipelineBase, Source
from mayavi.core.ui.api import SceneEditor, MayaviScene, \
                                MlabSceneModel

################################################################################
# Create some data
data = np.load('den_test.dat')

################################################################################
# The object implementing the dialog
class VolumeSlicer(HasTraits):
    # The data to plot
    data = Array()

    # The 4 views displayed
    scene3d = Instance(MlabSceneModel, ())
    scene_x = Instance(MlabSceneModel, ())
    scene_y = Instance(MlabSceneModel, ())
    scene_z = Instance(MlabSceneModel, ())

    # The data source
    data_src3d = Instance(Source)

    # The image plane widgets of the 3D scene
    ipw_3d_x = Instance(PipelineBase)
    ipw_3d_y = Instance(PipelineBase)
    ipw_3d_z = Instance(PipelineBase)

    _axis_names = dict(x=0, y=1, z=2)


    #---------------------------------------------------------------------------
    def __init__(self, **traits):
        super(VolumeSlicer, self).__init__(**traits)
        # Force the creation of the image_plane_widgets:
        self.ipw_3d_x
        self.ipw_3d_y
        self.ipw_3d_z


    #---------------------------------------------------------------------------
    # Default values
    #---------------------------------------------------------------------------
    def _data_src3d_default(self):
        return mlab.pipeline.scalar_field(self.data,
                            figure=self.scene3d.mayavi_scene)

    def make_ipw_3d(self, axis_name):
        ipw = mlab.pipeline.image_plane_widget(self.data_src3d,
                        figure=self.scene3d.mayavi_scene,
                        plane_orientation='%s_axes' % axis_name)
        return ipw

    def _ipw_3d_x_default(self):
        return self.make_ipw_3d('x')

    def _ipw_3d_y_default(self):
        return self.make_ipw_3d('y')

    def _ipw_3d_z_default(self):
        return self.make_ipw_3d('z')


    #---------------------------------------------------------------------------
    # Scene activation callbaks
    #---------------------------------------------------------------------------
    @on_trait_change('scene3d.activated')
    def display_scene3d(self):
        outline = mlab.pipeline.outline(self.data_src3d,
                        figure=self.scene3d.mayavi_scene,
                        )
        self.scene3d.mlab.view(40, 50)
        # Interaction properties can only be changed after the scene
        # has been created, and thus the interactor exists
        for ipw in (self.ipw_3d_x, self.ipw_3d_y, self.ipw_3d_z):
            # Turn the interaction off
            ipw.ipw.interaction = 0
        self.scene3d.scene.background = (0, 0, 0)
        # Keep the view always pointing up
        self.scene3d.scene.interactor.interactor_style = \
                                 tvtk.InteractorStyleTerrain()


    def make_side_view(self, axis_name):
        scene = getattr(self, 'scene_%s' % axis_name)

        # To avoid copying the data, we take a reference to the
        # raw VTK dataset, and pass it on to mlab. Mlab will create
        # a Mayavi source from the VTK without copying it.
        # We have to specify the figure so that the data gets
        # added on the figure we are interested in.
        outline = mlab.pipeline.outline(
                            self.data_src3d.mlab_source.dataset,
                            figure=scene.mayavi_scene,
                            )
        ipw = mlab.pipeline.image_plane_widget(
                            outline,
                            plane_orientation='%s_axes' % axis_name)
        setattr(self, 'ipw_%s' % axis_name, ipw)

        # Synchronize positions between the corresponding image plane
        # widgets on different views.
        ipw.ipw.sync_trait('slice_position',
                            getattr(self, 'ipw_3d_%s'% axis_name).ipw)

        # Make left-clicking create a crosshair
        ipw.ipw.left_button_action = 0
        # Add a callback on the image plane widget interaction to
        # move the others
        def move_view(obj, evt):
            position = obj.GetCurrentCursorPosition()
            for other_axis, axis_number in self._axis_names.items():
                if other_axis == axis_name:
                    continue
                ipw3d = getattr(self, 'ipw_3d_%s' % other_axis)
                ipw3d.ipw.slice_position = position[axis_number]

        ipw.ipw.add_observer('InteractionEvent', move_view)
        ipw.ipw.add_observer('StartInteractionEvent', move_view)

        # Center the image plane widget
        ipw.ipw.slice_position = 0.5*self.data.shape[
                    self._axis_names[axis_name]]

        # Position the view for the scene
        views = dict(x=( 0, 90),
                     y=(90, 90),
                     z=( 0,  0),
                     )
        scene.mlab.view(*views[axis_name])
        # 2D interaction: only pan and zoom
        scene.scene.interactor.interactor_style = \
                                 tvtk.InteractorStyleImage()
        scene.scene.background = (0, 0, 0)


    @on_trait_change('scene_x.activated')
    def display_scene_x(self):
        return self.make_side_view('x')

    @on_trait_change('scene_y.activated')
    def display_scene_y(self):
        return self.make_side_view('y')

    @on_trait_change('scene_z.activated')
    def display_scene_z(self):
        return self.make_side_view('z')


    #---------------------------------------------------------------------------
    # The layout of the dialog created
    #---------------------------------------------------------------------------
    view = View(HGroup(
                  Group(
                       Item('scene_y',
                            editor=SceneEditor(scene_class=Scene),
                            height=250, width=300),
                       Item('scene_z',
                            editor=SceneEditor(scene_class=Scene),
                            height=250, width=300),
                       show_labels=False,
                  ),
                  Group(
                       Item('scene_x',
                            editor=SceneEditor(scene_class=Scene),
                            height=250, width=300),
                       Item('scene3d',
                            editor=SceneEditor(scene_class=MayaviScene),
                            height=250, width=300),
                       show_labels=False,
                  ),
                ),
                resizable=True,
                title='Volume Slicer',
                )


m = VolumeSlicer(data=data)
m.configure_traits()

        

