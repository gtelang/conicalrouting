#!/usr/bin/python2
#-----------------------------------------------------------------------------------------
import networkx as nx
import scipy as sp
import numpy as np
from colorama import Fore, Style
import sys, os, random, time, yaml
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import pyplot as plt, patches
from CGAL import CGAL_Kernel                                                                                                                                                                         
from CGAL.CGAL_Kernel import Point_2              
from CGAL.CGAL_Kernel import Ray_2              
from CGAL.CGAL_Kernel import Vector_2              
from CGAL.CGAL_Kernel import Segment_2              
from CGAL.CGAL_Kernel import Triangle_2 

#-----------------------------------------------------------
# For reproducible runs involving random choices
rseed = 10
np.random.seed(rseed) 


#-----------------------------------------------------------------------------------------
# For Enabling Latex inside plots
plt.rcParams.update({                                                                 
  "text.usetex": True        ,                                                               
  "font.family": "serif"     ,
  "font.serif": ["Palatino"] ,       
 })


#--------------------------------------------------------------------------------------------
# Generic palette of colors to use. all these colors are sequentially distinct. 
# https://matplotlib.org/stable/gallery/color/named_colors.html See this 
# webpage for how these colors look
colorlist = ['red'        , 'green'    ,   'blue', 
             'purple'     , 'salmon'   , 'orange'   , 
             'chartreuse' , 'crimson'  , 'indigo'   , 
             'dodgerblue' , 'darkkhaki', 'limegreen', 
             'teal'       , 'fuchsia'  , 'moccasin' , 
             'firebrick']


#-----------------------------------------------------------------------------------------
# For random colors during plotting
def get_random_color():
    
   r = np.random.rand()
   g = np.random.rand()
   b = np.random.rand()
    
   return (r,g,b)   

#-----------------------------------------------------------------------------------------
def conical_routing(fig, ax, start, target, box, dx):
     """ Give a route from the start to the target using the conical routing scheme
     amidst a box that contain stalactite and stalagmite obstacles. See write-up for
     details on the mathematical details of the scheme. 
     """
     obstsegs     = box.get_obstsegs()
     dx = 0.08
     N  = 10

     spathverts = [box.start] + box.get_obstips() + [box.target]
     conebunches = []
     for i, (left_vert, obstacle_tip, right_vert) in enumerate(zip(spathverts, spathverts[1:], spathverts[2:])):
               print "Generating cone fan at obstacle tip: ",  obstacle_tip
               cones = make_cone_fan( fig, ax, 
                                      obstacle_tip, dx, obstsegs, 
                                      left_vert, 
                                      right_vert, 
                                      wedge_numangles = N, 
                                      color           = colorlist[i])
               conebunches.append(cones)

     # Extract the first arm of each cone 
     armbunches = []
     for conebunch in conebunches:
           armbunch = []
           for cone in conebunch:
              armbunch.append(cone.seg1)
           armbunches.append(armbunch) 

     lastarm = Segment(box.get_obstips()[-1], target) 
     # Render all the armbunches, each bunch in a different color 
     for i, armbunch in enumerate(armbunches):
          for arm in armbunch:
             arm.draw(fig,ax, color=colorlist[i], alpha=0.8)

     
     


#------------------------------------------------------------------------------------------
def make_cone_fan( fig, ax,  obstacle_tip, dx, obstsegs, 
                   left_vert, right_vert,   
                   wedge_numangles, color  ):
     """ Make the cone fan at the given obstacle tip given the bend points of the shortest path to the left 
     and right of the obstacle tip. 
     """
     first_ray       = Ray(obstacle_tip, left_vert-obstacle_tip)
     last_ray        = Ray(obstacle_tip , right_vert-obstacle_tip)
     cones_arms_rays = rock_line(obstacle_tip, first_ray, last_ray, 
                                 wedge_numangles=wedge_numangles)


     midsegs=[]
     for i, raybunch in enumerate(cones_arms_rays):
          if i==0:
               raybunch = raybunch[1:] # disregard the first ray from the first bunch. This omission if fixed after this for loop
          elif i==2:
               raybunch = raybunch[:-1] # disregard the last ray from the last bunch. This omission is fixed after this for loop
               
          for ray in raybunch:
               midsegs.append( Segment( obstacle_tip, ray_shoot(ray, obstsegs)) )

     cones_arms_segs = [Segment(obstacle_tip,  left_vert)] +\
                       midsegs                             +\
                       [Segment(obstacle_tip,  right_vert)]

     cones           = [Cone(obstacle_tip, segl, segr, dx, eps=0.04) 
                        for (segl, segr) in zip(cones_arms_segs, cones_arms_segs[1:])]

     #for cone in cones:
     #     cone.draw(fig,ax, edgecolor=color)
     
     return cones 



#------------------------------------------------------------------------------------------
class Box:
  def __init__(self, xleft, xright, ydown, yup, segs, start=None, target=None):
    self.xleft  = xleft
    self.xright = xright
    self.yup    = yup    
    self.ydown  = ydown
    self.ites   = segs # should be a list of `Segments` 
    self.start  = start
    self.target = target

  def llc (self):
     return np.asarray([self.xleft, self.ydown])
     
  def ulc (self):
     return np.asarray([self.xleft, self.yup])

  def lrc (self):
     return np.asarray([self.xright, self.ydown])
     
  def urc (self):
     return np.asarray([self.xright, self.yup])

  def get_obstsegs (self):
       wall1 = Segment(self.llc(), self.lrc()) # lower wall 
       wall2 = Segment(self.ulc(), self.urc()) # upper wall
       wall3 = Segment(self.llc(), self.ulc()) # left wall
       wall4 = Segment(self.lrc(), self.urc()) # right wall
       obstsegs = list(self.ites) + [wall1, wall2, wall3, wall4] 
       return obstsegs

  def get_obstips(self):
        """ The obstacle tip is assumed to be the 2nd point of each stalactite stalagmite segment
        Should probablty make this into a separae dtaa type. 
        """
        return [ite.q for ite in self.ites]


  def draw(self,fig, ax):
     boxpatch = mpl.patches.Rectangle(  self.llc()                     , 
                                        width  = self.xright - self.xleft       , 
                                        height = self.yup    - self.ydown       ,
                                        angle=0.0,  facecolor='yellow' , 
                                        alpha=0.3, zorder = 0  )
     ax.add_patch(boxpatch)
     
     ax.plot(  [self.xleft  ,self.xright]      , [self.ydown, self.ydown], 'ko-') # bottom edge
     ax.plot(  [self.xleft  ,self.xright]      , [self.yup, self.yup]    , 'ko-') # top edge

     ax.plot(  [self.xleft  ,self.xleft]  , [self.ydown, self.yup], 'ko-')   # left edge
     ax.plot(  [self.xright ,self.xright] , [self.ydown, self.yup], 'ko-')   # right edge

     # draw all the stalactites and stalagmites
     for ite in self.ites:
          ite.draw(fig,ax)

#------------------------------------------------------------------------------------------
def cgalpt (p):
     """ Convert a tuple or numpy array of length 2 into a CGAL point
     """
     return Point_2(p[0],p[1])
#------------------------------------------------------------------------------------------
class Segment :
     def __init__(self, p, q):
           self.p = np.asarray(p)
           self.q = np.asarray(q)
           
     def draw(self,fig,ax, color='black', alpha=0.5):
          ax.plot( [self.p[0],self.q[0]], [self.p[1],self.q[1]], 'o-', color=color, alpha=alpha, lw=0.7 )
          
     def draw_directed(self,fig,ax,color='blue', alpha=0.5):
          """ Draw a directed segment as directed arrow from the member point p to the member point q
          """
          dx = self.q[0]-self.p[0]
          dy = self.q[1]-self.p[1]
          ax.arrow(self.p[0], self.p[1], dx,dy, head_width=0.02, 
                   edgecolor=color, 
                   facecolor=color, 
                   length_includes_head=True, alpha=alpha)

     def tocgal(self):
          p, q = self.p, self.q
          return Segment_2(Point_2(p[0], p[1]), Point_2(q[0], q[1]))
     
     def reverse(self):
          """ This reverse the orientation of a segment, swapping p and q
          """
          return Segment(self.q, self.p)

     def __str__(self):
          return Fore.GREEN  + "Segment (" + str(self.p) + ", " + str(self.q) + ")" + Style.RESET_ALL
#------------------------------------------------------------------------------------------
class Ray :
     def __init__(self, p, dirn):
           self.p    = np.asarray(p)
           self.dirn = np.asarray(dirn)/np.linalg.norm(dirn) # normalize for future convenience

     def tocgal(self):
          return Ray_2 (Point_2(self.p[0], self.p[1]), Vector_2(self.dirn[0] , self.dirn[1] ))

     
     def reverse(self):
          """ Reverses a ray with opposite direction and same base point
          """
          return Ray(self.p, -self.dirn)

     def draw(self,fig,ax, length=0.5, color='purple', alpha=0.5):
          """ Draw an arrow with base at p
          """
          p = self.p
          q = p + length * self.dirn
          dx =q[0]-p[0]
          dy =q[1]-p[1]
          ax.arrow(p[0], p[1], dx,dy, head_width=0.02, 
                   edgecolor=color, 
                   facecolor=color, 
                   length_includes_head=True, alpha=alpha)



     def __str__(self):
          return Fore.YELLOW  + "Ray :: (Base: " + str(self.p) + "),  (Dirn:" + str(self.dirn) + ")" + Style.RESET_ALL


#------------------------------------------------------------------------------------------
class Cone:
    def __init__(self, p, seg1, seg2, dx, eps):
           self.p = np.asarray(p)
           
           assert np.linalg.norm(seg1.p-seg2.p) <= 1e-6, "No cone for you. Bases of Segments  r1 and r2 must match." 
           self.seg1 = seg1
           self.seg2 = seg2
           
           self.stpts_seg1 = seg_steiner_pts(seg1,dx, eps=eps)
           self.stpts_seg2 = seg_steiner_pts(seg2,dx, eps=eps)

    def draw(self, fig, ax, edgecolor, pointcolor='yellow' ):
           self.seg1.draw_directed(fig,ax,color=edgecolor)
           self.seg2.draw_directed(fig,ax,color=edgecolor)

           for pt in self.stpts_seg1:
                draw_point(fig, ax, pt, radius = 0.005, facecolor=pointcolor)
           for pt in self.stpts_seg2:
                draw_point(fig, ax, pt, radius = 0.005, facecolor=pointcolor)

                
           # for p in self.stpts_seg1:
           #      for q in self.stpts_seg2:
           #           seg = Segment(p,q)
           #           seg.draw(fig,ax, color='gray', alpha=0.5)


#------------------------------------------------------------------------------------------
def seg_seg_intersect_p(sega, segb):
     """ This works evern if two segments overlap completely or partially.
     """
     sa     = Segment_2(Point_2(sega.p[0], sega.p[1]), Point_2(sega.q[0], sega.q[1]))
     sb     = Segment_2(Point_2(segb.p[0], segb.p[1]), Point_2(segb.q[0], segb.q[1]))
     object = CGAL_Kernel.intersection(sa, sb)

     intersect_p = False
     if object.is_Segment_2() or object.is_Point_2(): 
          intersect_p = True
     return intersect_p

#------------------------------------------------------------------------------------------

def seg_seg_intersection_point(sega,segb):
     """ Warning: if two segments overlap i.e. lie flush with one another, completely or partially
     the routine returns None, because there is no canonical point to select for the segment segment 
     intersection. 
     """
     sa     = Segment_2(Point_2(sega.p[0], sega.p[1]), Point_2(sega.q[0], sega.q[1]))
     sb     = Segment_2(Point_2(segb.p[0], segb.p[1]), Point_2(segb.q[0], segb.q[1]))
     object = CGAL_Kernel.intersection(sa, sb)
     
     intersection_point = None
     if object.is_Point_2(): 
          pt                 = object.get_Point_2()
          intersection_point = np.asarray([pt.x(), pt.y()])
     elif object.is_Segment_2():
          assert False, "Segment Overlap detected! seg_seg_intersection_point makes sense only for non-overlapping (complete or partial) segments, i.e. they must have either exactly one or no points in common, for the routine to work as intended.,"

     return intersection_point
#------------------------------------------------------------------------------------------

def ray_seg_intersect_p(ray, seg):
     s = Segment_2(Point_2(seg.p[0], seg.p[1]), Point_2(seg.q[0]    , seg.q[1]   ))
     r = Ray_2    (Point_2(ray.p[0], ray.p[1]), Vector_2(ray.dirn[0] , ray.dirn[1] ))
     
     object = CGAL_Kernel.intersection(r, s)

     intersect_p = False
     if object.is_Segment_2() or object.is_Point_2(): 
          intersect_p = True
     return intersect_p

#------------------------------------------------------------------------------------------
def ray_seg_intersection_point(ray, seg):
     """ Warning: if the ray and segment overlap i.e. segment lies flush along the ray. 
     the routine returns the endpoint of the overlapping segment that is closer to the 
     ray's base, since that is where the ray first `interacts` with the segment obstacle
     """
     r = Ray_2    (Point_2(ray.p[0], ray.p[1]), Vector_2(ray.dirn[0] , ray.dirn[1] ))
     s = Segment_2(Point_2(seg.p[0], seg.p[1]), Point_2(seg.q[0]    , seg.q[1]   ))
     
     object = CGAL_Kernel.intersection(r, s)

     intersecn_pt = None
     if object.is_Point_2(): 
          pt = object.get_Point_2()
          intersecn_pt = np.asarray([pt.x(), pt.y()])

     elif object.is_Segment_2():
          iseg = object.get_Segment_2()
          e1   = np.asarray([iseg.source().x(), iseg.source().y()])
          e2   = np.asarray([iseg.target().x(), iseg.target().y()])
               
          d1 = np.linalg.norm(e1-ray.p)
          d2 = np.linalg.norm(e2-ray.p)
               
          if d1 <= d2:
               pt = e1
          else:
               pt = e2
          intersecn_pt = pt

     return intersecn_pt

#------------------------------------------------------------------------------------------
def ray_shoot(ray, obstseg):
     """ For a given ray and obstacle segments find the point on the obstacle segments closest to 
     the base point of the ray. Note that since all ray shooting will be done from inside the box, 
     at least one intersection point with some obstacle is guaranteed. 
     """
     
     #print Fore.YELLOW, "Shooting ray ", ray 
     intersection_points = []
     for seg in obstseg:
          object = CGAL_Kernel.intersection(ray.tocgal(), seg.tocgal())

          intersecn_pt = None
          if object.is_Point_2(): 
               pt = object.get_Point_2()
               intersecn_pt = np.asarray([pt.x(), pt.y()])

          elif object.is_Segment_2():
               iseg = object.get_Segment_2()
               e1   = np.asarray([iseg.source().x(), iseg.source().y()])
               e2   = np.asarray([iseg.target().x(), iseg.target().y()])
               
               d1 = np.linalg.norm(e1-ray.p)
               d2 = np.linalg.norm(e2-ray.p)
               
               if d1 <= d2:
                    pt = e1
               else:
                    pt = e2
               intersecn_pt = pt

          intersection_points.append(intersecn_pt)

     #print "..........Intersection points for ray" , ray, " are--->\n" , intersection_points
     dmin  = np.inf
     ptmin = None 
     for ipt in intersection_points:
          if  (ipt is not None)                 and \
              (cgalpt(ipt) != cgalpt(ray.p))    and \
              np.linalg.norm(ipt-ray.p) < dmin:
                     ptmin = ipt
                     dmin  = np.linalg.norm(ipt-ray.p)

     assert ptmin is not None, "At least one intersection point of a ray with a given set of obstacles must exists as we are ray shootinfg from inside the box/."

     return ptmin

               

#-------------------------------------------------------------------------------------------
def mutually_visible_p(p,q, obstsegs):
     """ Two points p and q are considered mutually visible iff 
     the segment pq lies in the closure of the free space. 
     Here the complement of the free space is represented by the list of 
     obstacle segments, which in this case means the walls of the box 
     and the stalactites and stalagmites. 
     """
     pq       = Segment(p,q).tocgal()
     obstsegs = [seg.tocgal() for seg in obstsegs]
     
     globflag = True
     for seg in obstsegs:
          locflag = False
          object  = CGAL_Kernel.intersection(pq, seg)
          
          # Check if the intersection of two segments is a segment itself. 
          if object.is_Segment_2():
               locflag = True

          # Check if the intersection of two segments is a point, whether it is
          # one of the obstacle segment endpoints. 
          elif object.is_Point_2():
               intpt = object.get_Point_2()
               if intpt == cgalpt(seg.p) or intpt == cgalpt(seg.q): #### Does the == use CGAL routines for comparing points? I made some sanity checks with Python code , and it seems so, but will later have to be sure. 
                    locflag = True

          globflag = globflag and locflag
          
     return globflag

#---------------------------------------------------------------------------------------
def seg_steiner_pts(seg,dx, eps=1e-5):
     """ Generate Steiner points along the given segment 
     with two consecutive points spaced dx distance apart. 
     All points are placed in the interior of the segment
     starting at (seg.p). eps represents the shift 
     wrt the initial point p. if eps = 0, then the first steiner point
     if self.p itself. 
     """
     dirn        = (seg.q-seg.p)
     dirn        = dirn/np.linalg.norm(dirn) # make direction into a unit vector
     steiner_pts = []

     ctr = 1
     while True:
         newpt =  seg.p + (eps + (ctr-1)*dx) * dirn
         if (newpt[0] - seg.p[0]) * (newpt[0] - seg.q[0]) <= 0 and \
            (newpt[1] - seg.p[1]) * (newpt[1] - seg.q[1]) <= 0 :
              steiner_pts.append( newpt )
              ctr   +=  1
         else:
              break # no more steiner points to add!

     return steiner_pts
#------------------------------------------------------------------------------------------
def draw_point(fig,ax, pt, radius = 0.03, facecolor='red', edgecolor='black', zorder=3):
     ptmark = mpl.patches.Circle(pt, 
                              radius    = radius, 
                              facecolor = facecolor, 
                              edgecolor = edgecolor,
                              zorder=3)
     ax.add_patch(ptmark)

#------------------------------------------------------------------------------------------

def draw_curve(fig, ax, pts, marker='o',  color='green', linewidth=1.3, linestyle='dashed'):
     """ Render a curve given a list of points. This is more convenient thant using 
     default plot of matplotlib where we have to separate the x and y coordinates. 
     This function is essentially a wrapper around that plot function/. 
     """
     xs = [pt[0] for pt in pts]
     ys = [pt[1] for pt in pts]

     ax.plot(xs, ys, 
             color     =color, 
             linewidth =linewidth, 
             marker    =marker, 
             linestyle =linestyle)
#---------------------------------------------------------------------------------------
def rock_line(obstacle_tip, first_ray, last_ray, wedge_numangles=10):
     """ Assuming no obstacles present and just two rays both sharing the same base at the 
     obstacle tip. Generate a sequence of rays obtained by rocking a line as mentioned in 
     the writeup. Returns a list of three lists of rays. Each sublist corresponds to 
     a list of rays inside of a wedge in the order yellow,---green---yellow as encountered
     going from the start to the destination. 
     """
     #assert (first_ray.p).tocgal() == (last_ray.p).tocgal(), " Base points of the first_ray and last_ray are not equal"
     wedges = [[ first_ray          , last_ray.reverse()],  # yellow wedge
               [ last_ray.reverse() , first_ray.reverse()],  # green wedge
               [ first_ray.reverse(), last_ray]]   # yellow wedge
     cones_arms_rays = []
     for i, wedge in enumerate(wedges):
          start_dirn          = wedge[0].dirn/np.linalg.norm(wedge[0].dirn)
          end_dirn            = wedge[1].dirn/np.linalg.norm(wedge[1].dirn)
          angle_between_dirns = np.arccos(np.dot(start_dirn, end_dirn))   
          
          if i == 0 or i==1:
               angles = np.linspace(0, angle_between_dirns, wedge_numangles) [:-1]# angle corresponding to last ray is dropped since it is the first ray of the next wedge
          else:
               angles = np.linspace(0, angle_between_dirns, wedge_numangles) # no angle dropped here!

          rays = []
          for angle in angles:
               new_dirn = None
               if np.cross(start_dirn, end_dirn) < 0: # right hand rule makes thumb point downward
                  new_dirn = rot_clock(start_dirn, angle)
               else:                                  # ditto, but  thumb now upward
                  new_dirn = rot_anticlock(start_dirn, angle)
               rays.append(   Ray( obstacle_tip , new_dirn)   )

          cones_arms_rays.append(rays)

          

     return cones_arms_rays
#-------------------------------------------------------------------------------------------
def rot_anticlock(vec, angle):
     """ Rotate a vector anticlockwise by angle theta
     """
     vec    = np.asarray(vec)
     c      = np.cos(angle)
     s      = np.sin(angle)
     rotmat = np.asarray([[c,-s],
                          [s,c]])
     return rotmat.dot(vec)
#-------------------------------------------------------------------------------------------
def rot_clock(vec, angle):
     """ Rotate a vector clockwise by angle theta
     """
     vec    = np.asarray(vec)
     c      = np.cos(angle)
     s      = np.sin(angle)
     rotmat = np.asarray([[c,s],
                          [-s,c]])
     return rotmat.dot(vec)
#-------------------------------------------------------------------------------------------

def example_environment ():
     start  = np.asarray([0,0])
     target = np.asarray([3,1])
     xleft  = start[0]
     xright = target[0]
     ydown  = start[1]
     yup    = target[1]
     #stalagmites (mites go up)                                                     
     mtips  = list(map(np.asarray, [[0.5,0.7], 
                                    [1.2,0.5],  
                                    [1.8,0.7]]))
     mites = []                     
     for i in range(len(mtips)):
         p = np.asarray([mtips[i][0],ydown])
         q = mtips[i]
         mites.append( Segment(p,q) )

     #stalactites (tites go down)                                                     
     ttips  = list(map(np.asarray, [(0.8,0.2), 
                                    (1.5,0.3),  
                                    (2.5,0.2)]))
     tites = []
     for i in range(len(ttips)):
         p = np.asarray([ttips[i][0],yup])
         q = ttips[i]
         tites.append (Segment(p,q))

     ites = []
     for m,t in zip(mites, tites):
          ites.extend([m,t])
     # Create the environment box. 
     box = Box(xleft,xright,ydown,yup, ites, start=start, target=target);
     return box

#---------------------------------------------------------------------------------


def main():
     # Generate one of several possible environments using a dedicated function. 
     # Each environment consists of a start a target and a box that in turn 
     # contains stalactites and stalagmites. 
     box    = example_environment()
     start  = box.start
     target = box.target

     # Set up the Matplotlib canvas. 
     fig, ax = plt.subplots()
     ax.set_aspect(1.0)
     ax.set_title("Conical routing for generating lookahead curves")
     from matplotlib.ticker import (AutoMinorLocator, MultipleLocator)
     ax.xaxis.set_major_locator(MultipleLocator(0.1))
     ax.yaxis.set_major_locator(MultipleLocator(0.1))
     ax.grid(color = 'gray', linestyle = '--', linewidth = 0.3, alpha     = 0.5)

     # Draw the box
     box.draw(fig,ax)

     # Draw the starting and ending points
     draw_point(fig,ax,start)
     draw_point(fig,ax,target)

     # Perform computation to compute lookaheadcurve. since fig, ax
     # have been passedd to the routine, this will enable graphical debugging. 
     # later we can move the graphical plotting part outside the function. 
     lookahead_curve = conical_routing(fig, ax, start, target, box, dx=0.1)
     plt.show()
#-----------------------------------------------------------------------------    
if __name__=="__main__":
     main()
