#!/usr/bin/python2
"""
Solving Lookahead problems via conical routing.
Author : Gaurish Telang
"""
#--------------------------------------------------------------------------------------------------------------
# Python 2 is used here because CGAL seems to work only with Python 2 on my machine. 
#---------------------------------------------------------------------------------------------------------------
# TODO ::: in the case where the last subwedge in the green wedge is very small fuse the arm 
# with the previous or the latter cone. This small angle might happen in some cases, and you need 
# to be careful....for the testing case, you don't need to worry. since all angles are fat. 
# and more or less equal. 
import sys, os, random, time, yaml
from colorama import Fore, Style
# For scientific calculations, vectors matrices, linear algebra
import scipy as sp
import numpy as np
import itertools
import math
# For plotting/rendering etc.
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import pyplot as plt, patches
# For icky geometric computations that need to handled exactly
from CGAL import CGAL_Kernel                                                                                                                                                                         
from CGAL.CGAL_Kernel import Point_2              
from CGAL.CGAL_Kernel import Ray_2              
from CGAL.CGAL_Kernel import Vector_2              
from CGAL.CGAL_Kernel import Segment_2              
from CGAL.CGAL_Kernel import Triangle_2 
# For Graph theoretic stuff
import networkx as nx
# For solving integer, 0-1 and linear programs. 
# See https://coin-or.github.io/pulp/CaseStudies/a_blending_problem.html for a good example. 
# This modelling software is widely used by the operations research community. 
# It has default solvers included but can also interact with CPLEX and many other powerful solvers. 
import pulp

#-------------------------------------------------------------------------------------------------
# Fixed initial seed for reproducible runs involving random choices/randomization/random numbers
rseed = 10
np.random.seed(rseed) 
#-------------------------------------------------------------------------------------------------
# For Enabling Latex inside plots
plt.rcParams.update({                                                                 
  "text.usetex": True        ,                                                               
  "font.family": "serif"     ,
  "font.serif": ["Palatino"] ,       
 })

#--------------------------------------------------------------------------------------------
# Generic all-purpose palette of colors. All these colors are sequentially distinct. 
# https://matplotlib.org/stable/gallery/color/named_colors.html See this 
# webpage for how these colors look. Wayyyy better than using random colors
# which often given horrible selections, sometimes two successively generated colors are very close 
# in appearance.
#--------------------------------------------------------------------------------------------
colorlist = ['red'        , 'green'    , 'blue'     , 
             'purple'     , 'salmon'   , 'orange'   , 
             'chartreuse' , 'crimson'  , 'indigo'   , 
             'dodgerblue' , 'darkkhaki', 'limegreen', 
             'teal'       , 'fuchsia'  , 'moccasin' , 
             'firebrick']

#------------------------------------------------------------------------------------------------
def printlist(l):
     """ Print out a list with each item on its own line. 
     This is a great convenience function while debugging.
     """
     for i, item in enumerate(l):
          print "[", (i+1), "]", item

#-------------------------------------------------------------------------------
def partial_sums(xs, pad_zero_at_start_p=False):
     """ Given xn array [x1, x2, ... xn]
     return the sequence of partial sums 
     [x1, 
     x1+x2, 
     x1+x2+x3,
     ...
     x1+x2....+xn]
     i.e. [sum(xs[:i]) for i in range(1,len(xs)+1)]
     
     if pad_zero_start_p==True then return
     
     [0,
     x1,
     x1+x2,
     ....
     x1+x2_....+xn]

     i.e. [sum(xs[:i]) for i in range(0,len(xs)+1)]
     """
     
     psums       = []

     if pad_zero_at_start_p: 
          psums.append(0)

     sum_so_far  = 0
     for x in xs:
          psums.append(sum_so_far+x)
          sum_so_far += x
     return psums

#--------------------------------------------------------------------------------

def min_length_intprog_lowerbound_via_lp_relaxation(G, start_node_idx, target_node_idx, lookahead_lb, hplanes, makeS_p=False, lpfilename="lookahead_lp_relxation.lp"):
     assert nx.is_directed(G), " G should be a directed graph."
     print (" .......................................................")
     print ("Getting lookahead lower bound via lp relaxation   ")
     print ("........................................................")
     #....................................
     # INITIALIZE PROBLEM
     #....................................
     prob = pulp.LpProblem("lookahead lower bound via lp relxation", pulp.LpMinimize)
     
     #....................................
     # SETUP VARIABLES
     #...................................
     edges      = [ (u,v) for u,v in G.edges()]                                  
     edge_costs = { (u,v) : attr['weight'] for u,v,attr in G.edges(data=True)  } 
     edge_vars  = pulp.LpVariable.dicts("Edges",edges,lowBound=0.0, upBound=1.0, cat='Continuous') 
     
     #....................................
     # SETUP OBJECTIVE
     #...................................
     prob += pulp.lpSum([edge_costs[e]*edge_vars[e] for e in edges]), "Objective: total cost of edges"

     #....................................
     # SETUP CONSTRAINTS
     #...................................
     #----> Type [A] constraints
     for coneidx in sorted(list({a for ((a,_),_) in edges})): 
          print Fore.YELLOW, ".....Generating Type [A] (one-choice-per-cone) constraint for cone ",  coneidx,  Style.RESET_ALL
          prob += pulp.lpSum([edge_vars[(u,v)] for (u,v) in edges if u[0] == coneidx]) == 1 , "Cone " + str(coneidx) + " contributes exactly 1 red edge"

     #-----> Type [B] constraints
     print("......................\n")
     print Fore.YELLOW, ".....Generating Type [B] (start/target edges) constraint for cone 0" ,  Style.RESET_ALL
     numarms = len(list({a for ((a,_),_) in edges})) + 1 
     prob   += pulp.lpSum([edge_vars[(u,v)] for (u,v) in edges if (u==start_node_idx) and v[0] == 1               ]) == 1 , \
               "The red edge in Cone 0  has its tail vertex at the start node " 

     print Fore.YELLOW, ".....Generating Type [B] constraint for cone ",  coneidx,  Style.RESET_ALL
     prob   += pulp.lpSum([edge_vars[(u,v)] for (u,v) in edges if u[0] == numarms-2   and v    == target_node_idx] ) == 1 ,\
               "The red edge in Cone " + str(numarms-2) + "  has its head vertex at the target node " 
 
     #-----> Type [C] constraints
     print("......................\n")
     for nidx in list(G.nodes()):
          if nidx[0] not in [start_node_idx[0], target_node_idx[0]]: #flow constraints are not imposed at nodes belonging to the zeroth or the last arm.(i.e. arms to which the start and target respecitively belong)
               print Fore.YELLOW, ".....Generating Type [C] (flow) constraint for node ",  nidx,  Style.RESET_ALL
               inedges  = list(G.in_edges (nidx))
               outedges = list(G.out_edges(nidx))
               prob     += pulp.lpSum([edge_vars[ie] for ie in inedges])  == pulp.lpSum([edge_vars[oe] for oe in outedges]) , \
                          (" Inflow == Outflow at node " + str(nidx))
               
     #------> Type [D] constraints
     print("......................\n")
     # Remember! An hplane consists of a pair of arms that share the same basepoint but have opposite direction.
     for hplane in hplanes: 
            halfplane_edges = [(u,v) for (u,v) in list(G.edges()) if (u[0] >= hplane[0]) and v[0]<=hplane[1]   ] # a halfplane edge is a graph edge whose tail and head lies in the halfplane. By construction of the cones we can just use arm indices to decide this. 
            prob           += pulp.lpSum( [ edge_costs[hpe]*edge_vars[hpe] for hpe in halfplane_edges] ) >= lookahead_lb    , \
                              ("String in half-plane [" + str(hplane[0]) + '<--->' + str(hplane[1]) + "] should have length >= " +\
                               str(round(lookahead_lb,4)))

     # Problem data is written to an .lp file. This helps debugging. 
     print(".........................................")
     print("Writing model to file:---->" + lpfilename )
     print(".........................................")
     prob.writeLP(lpfilename)

     #.............................
     # SOLVE THE PROGRAM
     #.............................
     prob.solve() # The problem is solved using PuLP's choice of Solver
     print Fore.GREEN, "---> Status:", pulp.LpStatus[prob.status], Style.RESET_ALL # It can be one of "Not Solved", "Infeasible", "Unbounded", "Undefined" or "Optimal". 

     # Each of the variables is printed toa file with its resolved optimum value
     f= open("lowerbound_vector.txt",'w')
     for v in prob.variables():
          val = v.varValue
          f.write(     v.name + "=" + str(val) + "\n")
     f.close()

     non_zero_sol_edges = []
     for eidx in list(sorted(G.edges())): # sorted is important to get the edges as they appeat along the path.
          if edge_vars[eidx].varValue > 1e-5: # 1e-5 is just chosen as a tolerance value to separate the edge value from zero 
              non_zero_sol_edges.append(eidx)
     
     return pulp.value(prob.objective), non_zero_sol_edges

#-------------------------------------------------------------------------------------------------
def min_length_intprog(G, start_node_idx, target_node_idx, lookahead_lb, hplanes, makeS_p=False, lpfilename="lookahead.lp"):
     """
      This function implements the most direct "brute force" attack on solving the lookahead problem, by using the cone based discretization. 
      This will hopefully provide a good base case implementation. 
      NOTE: G should be a directed graph. Without this the loop setting up the type A constraints  will produce garbage. 

      NOTE: Variables of Type (2) and constraints of Type (E) stated below are not used when `makeS_p = False`
      Setting the  makeS_p flag to True ostensibly stops the curve from wiggling around and makes the curve nicer
      to drive along like an `S` where all inflection points are located near or on the shortest path edges joining the points.
      ................................................................
       There are 2 categories of VARIABLES associated with this problem
      ................................................................
      (1) A 0-1 variable associated with each edge of the graph. Each edge of 
          this base graph is considered colored RED 
      (2) A 0-1 variable associated with each pair of directed edges A-B such that head of A = tail of B
          and the associated obstacle tip satisfies certain sidedness constraints with respect to the  
          segment joining the tail of A and head of B. These segments can be considered as additional edges added to the graph G. 
          and are taken to be colored BLUE.  (This `red/blue` terminology is to make it easier to state the constraints.)

      The solution consists of a choice of red edges and blue edges satisfying the following constraints.  
      (Note that the required lookahead path itself is just the union of the red edges so reported by the solver)

      ..................................
      There are 5 groups of CONSTRAINTS:
      ..................................
      (A) Exactly *one* red edge must be chosen from each cone
      (B) The red edge chosen from the first cone  must have its tail at start_node_idx and the red 
          edge from the last cone must have its head at target_node_idx. Both can be imposed via summation 
          constraints at the start node. 
      (C) These red edges must form a *connected* sequence i.e. The head of the edge chosen from cone [i] 
          must be the tail of the edge chosen from cone[i+1] This requires the imposition of an appropriate 
          summation constraints at appropriate nodes. 
      (D) The amount of string in each half plane while rotating the halfplane from one shortest path edge to its successor 
          is >= lookahead_lb. (This is of course the most important constraint of the problem!) It also possibly 
          generates the longest expressions when written out.  However there will only be as many of these 
          constraints as the sum of the number of cones in the first yellow wedge associated with each obstacle tip. 
      (E) Exactly *one* blue edge must be chosen from each `expanded` cone (an expanded cone is the union of 
          cones [i] and [i+1] for each i)

     While 0-1 integer programming might be inefficient, there is probably a good rounding scheme here
     to give an efficient constant factor approximation. Possibly a knapsack style DP also works in giving
     a pseudo-polynomial $\vareps$ approximation scheme? 
     
     While using 0-1 integer restriction helps solve the lookahead problem exactly is possibly impracticaly, 
     solving the LP relaxation provides good lower bounds to test quality of various heuristics. 
     So far the only lower bound we knew was that of the shortest path, which is a stupid lower bound especially if the lookahead is large. 
     This lower bound can be used to evaluate the quality of Mayank's curve for the single segment obstacle
     case, and also for your two parameter method for the same case. 
     """
     
     assert nx.is_directed(G), " G should be a directed graph."
     print (" .......................................................")
     print ("Solving lookahead exactly via 0-1 integer programming   ")
     print ("........................................................")
     
     #....................................
     # INITIALIZE PROBLEM
     #....................................
     prob = pulp.LpProblem("Min Length lookahead via 0-1 integer programming", pulp.LpMinimize)
     
     #....................................
     # SETUP VARIABLES
     #...................................
     edges      = [ (u,v) for u,v in G.edges()]                                  
     edge_costs = { (u,v) : attr['weight'] for u,v,attr in G.edges(data=True)  } 
     edge_vars  = pulp.LpVariable.dicts("Edges",edges,cat='Binary') # A dictionary called 'edge_vars' is created to contain the referenced Variables. the keys are the elements of `edges`
     
     #....................................
     # SETUP OBJECTIVE
     #...................................
     prob += pulp.lpSum([edge_costs[e]*edge_vars[e] for e in edges]), "Objective: total cost of edges"

     #....................................
     # SETUP CONSTRAINTS
     #...................................
     #----> Type [A] constraints
     for coneidx in sorted(list({a for ((a,_),_) in edges})): 
          print Fore.YELLOW, ".....Generating Type [A] (one-choice-per-cone) constraint for cone ",  coneidx,  Style.RESET_ALL
          prob += pulp.lpSum([edge_vars[(u,v)] for (u,v) in edges if u[0] == coneidx]) == 1 , "Cone " + str(coneidx) + " contributes exactly 1 red edge"

     #-----> Type [B] constraints
     print("......................\n")
     print Fore.YELLOW, ".....Generating Type [B] (start/target edges) constraint for cone 0" ,  Style.RESET_ALL
     numarms = len(list({a for ((a,_),_) in edges})) + 1 
     prob   += pulp.lpSum([edge_vars[(u,v)] for (u,v) in edges if (u==start_node_idx) and v[0] == 1               ]) == 1 , \
               "The red edge in Cone 0  has its tail vertex at the start node " 

     print Fore.YELLOW, ".....Generating Type [B] constraint for cone ",  coneidx,  Style.RESET_ALL
     prob   += pulp.lpSum([edge_vars[(u,v)] for (u,v) in edges if u[0] == numarms-2   and v    == target_node_idx] ) == 1 ,\
               "The red edge in Cone " + str(numarms-2) + "  has its head vertex at the target node " 
 
     #-----> Type [C] constraints
     print("......................\n")
     for nidx in list(G.nodes()):
          if nidx[0] not in [start_node_idx[0], target_node_idx[0]]: #flow constraints are not imposed at nodes belonging to the zeroth or the last arm.(i.e. arms to which the start and target respecitively belong)
               print Fore.YELLOW, ".....Generating Type [C] (flow) constraint for node ",  nidx,  Style.RESET_ALL
               inedges  = list(G.in_edges (nidx))
               outedges = list(G.out_edges(nidx))
               prob     += pulp.lpSum([edge_vars[ie] for ie in inedges])  == pulp.lpSum([edge_vars[oe] for oe in outedges]) , \
                          (" Inflow == Outflow at node " + str(nidx))
               
     #------> Type [D] constraints
     print("......................\n")
     # Remember! An hplane consists of a pair of arms that share the same basepoint but have opposite direction.
     for hplane in hplanes: 
            halfplane_edges = [(u,v) for (u,v) in list(G.edges()) if (u[0] >= hplane[0]) and v[0]<=hplane[1]   ] # a halfplane edge is a graph edge whose tail and head lies in the halfplane. By construction of the cones we can just use arm indices to decide this. 
            prob           += pulp.lpSum( [ edge_costs[hpe]*edge_vars[hpe] for hpe in halfplane_edges] ) >= lookahead_lb    , \
                              ("String in half-plane [" + str(hplane[0]) + '<--->' + str(hplane[1]) + "] should have length >= " +\
                               str(round(lookahead_lb,4)))

     # All problem data is written to an .lp file. This helps debugging. 
     print(".........................................")
     print("Writing model to file:---->" + lpfilename )
     print(".........................................")
     prob.writeLP(lpfilename)

     #.............................
     # SOLVE THE PROGRAM
     #.............................
     prob.solve() # The problem is solved using PuLP's choice of Solver
     print Fore.GREEN, "---> Status:", pulp.LpStatus[prob.status], Style.RESET_ALL # It can be one of "Not Solved", "Infeasible", "Unbounded", "Undefined" or "Optimal". 
     

     f= open("intprog_sol_vector.txt",'w')
     for v in prob.variables():
          val = v.varValue
          f.write(     v.name + "=" + str(val) + "\n")
     f.close()

     #.........................................................................
     # Extract the path edges of the lookahead path from the solution vector
     #.........................................................................

     path_edges = []
     for eidx in list(sorted(G.edges())): # sorted is important to get the edges as they appeat along the path.
          if abs(edge_vars[eidx].varValue - 1) < 1e-5: # solver returns 0 and 1 as floats...hence we compare to 1, using abs()<tol method.
              path_edges.append(eidx)

     # simple sanity check on path returned.
     assert path_edges[0][0]   == start_node_idx , " "
     assert path_edges[-1][-1] == target_node_idx, " "
     
     #............................................
     # Get path nodes from edges extracted above
     #............................................
     path_nodes = []
     for (u,_) in path_edges:
          path_nodes.append(u)          # first node of each edge is marked as a path node.
     path_nodes.append(target_node_idx) # This is not appended in the for loop above, so we correct for that here.

     return path_nodes

#-------------------------------------------------------------------------------------------------
def conical_routing(fig, ax, box, dx, eps, W, lookahead_lb, lookaheadmethod=min_length_intprog, draw_mesh_p=True, draw_lp_relax_curve_p=True):
     """ Give a route from the start to the target using the conical routing scheme
     amidst a box that contain stalactite and stalagmite obstacles. See write-up for
     details on the mathematical details of the scheme. 
     lookahead_lb is a lower bound on the amount of lookahead needed for the path. 
     """
     assert W>=2              , "Number of angles in each yellow wedge around an obstacle ti[ ]should be >=2"
     assert lookahead_lb >= 0 , "Lookahead lower bound should be >= 0"

     obstsegs    = box.get_obstsegs()
     spathverts  = [box.start] + box.get_obstips() + [box.target]
     conebunches = []
     shiftvals   = []
     for i, (left_vert, obstacle_tip, right_vert) in enumerate(zip(spathverts, 
                                                                   spathverts[1:], 
                                                                   spathverts[2:])):
               print Fore.MAGENTA, "Generating cone fan at obstacle tip: ", obstacle_tip, Style.RESET_ALL
               cones, shiftval = make_cone_fan( fig, ax, 
                                                obstacle_tip, 
                                                dx              = dx, 
                                                eps             = eps, 
                                                obstsegs        = obstsegs, 
                                                left_vert       = left_vert, 
                                                right_vert      = right_vert, 
                                                wedge_numangles = W, 
                                                color           = colorlist[i])
               conebunches.append(cones)
               shiftvals.append(shiftval)

     #......................................................................................................
     # Extract the first arm of each cone 
     #......................................................................................................
     armbunches = []
     for conebunch in conebunches:
           armbunch = []
           for cone in conebunch:
              armbunch.append(cone.seg1)
           armbunches.append(armbunch) 

     #......................................................................................................
     # Add the arm joining the last obstacle tip to the target. 
     #......................................................................................................
     lastarm = Segment(box.get_obstips()[-1], box.target, dx=dx, eps=eps) 
     armbunches[-1].append(lastarm)
     
     #..............................................................................
     # Get the list of hplanes where the lookahead constraint needs to be imposed.
     #..............................................................................
     hplanes = []
     #offsets = [0] + [len(armbunch) for armbunch in armbunches]
     offsets = partial_sums( [len(armbunch) for armbunch in armbunches], pad_zero_at_start_p=True)

     assert len(shiftvals) == len(armbunches), " "
     assert len(offsets) == len(armbunches) + 1, " "

     for offset, shiftval in zip(offsets, shiftvals):
          for r in range(W+1): 
               hplanes.append([ r + offset  , r + offset +   shiftval ])

     #printlist(hplanes)

     #......................................................................................................
     # Get a flat list of cone arms connecting the start and the target. 
     #......................................................................................................
     flat_list_arms =[arm for armbunch in armbunches for arm in armbunch]

     #print "Total number of wedges are ", len(flat_list_arms) - 1
     #printlist(flat_list_arms)
     #......................................................................................................
     # Set up  complete bipartite directed graphs between every cone and take the union.
     #......................................................................................................
     G = nx.DiGraph() # empty DIRECTED graph. 
     
     #......................................................................................................
     # Add nodes to the graph by iterating over each node
     # in each arm. The edges are set in the next loop. 
     #......................................................................................................
     nodectr     = 0
     nodehash    = {} # hash table mapping (anum,ptnum) to a serial value nodectr. Useful when setting up a linear integer program
     nodehashinv = {} # inverse of hash table nodehash. Useful for extracting solutions. 
     for anum, a in enumerate(flat_list_arms):
          for ptnum, pt in enumerate(a.stpts):
               G.add_node(  (anum, ptnum) , cood=pt  )
               nodehash[(anum,ptnum)] = nodectr
               nodehashinv[nodectr]   = (anum,ptnum) 
               #print "Added node " , (anum, ptnum), " to hash table with id ", nodectr
               nodectr += 1
     #..........................................................................................................................................
     # Add edges to the graph by iterating over each successive cone. This is achieved by iterating over every arm 
     # and its successor arm (which of course together form a cone). Also store for each edge which cone it belongs to
     # A cone is identified by the index number of its arms. This might be useful later while setting up the integer/linear/dynamic programs. 
     #..........................................................................................................................................
     for a1num, (a1, a2) in enumerate(zip(flat_list_arms, flat_list_arms[1:])):
          a2num = a1num+1
          for pnum, p in enumerate(a1.stpts):
               for qnum, q in enumerate(a2.stpts):
                    G.add_edge( (a1num, pnum), (a2num, qnum), 
                                weight           = np.linalg.norm(p-q), 
                                coneid           = a1num              , # id of the cone containing the edge
                                circum_cone_arms = {a1num, a2num}     , # id of the arms of the cones. this information is redundant because of coneid, but included nevertheless for possible future convenience
                                edgeseg          = Segment(p,q))      
                    #print Fore.YELLOW, "Added edge between nodes ----> ", (a1num, pnum), (a2num, qnum) , "  to graph ", Style.RESET_ALL

     if draw_mesh_p:
          #......................................................................................................
          # Draw  the graph edges by rendering the geometrical segment associated with each edge of the graph
          #......................................................................................................
          draw_graph_edges(fig, ax, G,  alpha=0.05, markendpts_p=False)               

          #..............................................
          # Draw all the armbunches
          #..............................................
          draw_armbunches(fig,ax,armbunches,alpha=0.3) # this particular rendering is not absolutely necessary except for possible future debugging. The rendering of the graph itself to make out what the arms are. 


     #................................................................................................................
     # Compute and render the shortest path between start and target and print path length between start and target. 
     # At the moment this is just a preliminary computation / sanity check.  to make sure the graph has been computed correctly. 
     #................................................................................................................
     start_node_idx  = (0                    ,len(flat_list_arms[0].stpts) -1) # start is the last node of first arm
     target_node_idx = (len(flat_list_arms)-1,len(flat_list_arms[-1].stpts)-1) # target is the last node of the last arm
     spathnodeids    = nx.shortest_path(G, source=start_node_idx, target=target_node_idx, weight='weight')
     draw_path(fig,ax, G, spathnodeids,  color='red', linewidth = 2, alpha=1.0)

     #............................................................................
     # Solve the LP relaxation of the problem first, in order to get lower bounds
     #............................................................................

     path_length_lb, non_zero_sol_edges = min_length_intprog_lowerbound_via_lp_relaxation(G, start_node_idx, target_node_idx, lookahead_lb=lookahead_lb, hplanes=hplanes)
     #printlist (non_zero_sol_edges)
     if draw_lp_relax_curve_p:
          draw_edge_set(fig, ax, G, non_zero_sol_edges, color='magenta', alpha=1.0, linewidth=5, markendpts_p=False)
     
     #......................................................................................................
     # Solve the integer programming model here. 
     #......................................................................................................
     start_time      = time.time() 
     lookpathnodeids = lookaheadmethod(G,start_node_idx, target_node_idx, lookahead_lb=lookahead_lb, hplanes=hplanes)
     end_time        = time.time() 
     time_elapsed    = end_time - start_time
     print "Time elapsed for lookaheadmethod function to run : ", round(time_elapsed,4) , " seconds"


     return G,  lookpathnodeids, spathnodeids, {'lookpathlength'    : pathlength(G,lookpathnodeids), \
                                                'spathlength'       : pathlength(G,spathnodeids ),   \
                                                'path_length_lb'    : path_length_lb,                \
                                                'numcones'          : len(flat_list_arms)-1,                      \
                                                'non_zero_sol_edges': non_zero_sol_edges }

#---------------------------------------------------------------------------------------------
def pathlength(G,pathnodeidxs):
      """ Sum of the edge weights of the edges of a path. A path is represented by the list of its nodes. 
      """
      edgewtsum = 0.0
      for n1, n2 in zip(pathnodeidxs, pathnodeidxs[1:]):
              assert G.has_edge(n1,n2),  "Graph does not have edge " + str(n1) + " -- " + str(n2) + "....."
              edgewtsum +=  G.edges[n1, n2]['weight'] 

      return edgewtsum


#--------------------------------------------------------------------------------------------
def draw_armbunches(fig,ax,armbunches,alpha=0.3):
     """ Render all the armbunches, each bunch in a different color 
     """
     for i, armbunch in enumerate(armbunches):
          for arm in armbunch:
             arm.draw(fig,ax, color=colorlist[i], alpha=alpha)
     return

#-------------------------------------------------------------------------------------------
def draw_edge_set(fig, ax, G, edges, color, alpha=0.1, linewidth=2, markendpts_p=False):
     """ Render a set of edges from the graph G by rendering all their associated edge segments. 
     """
     
     for u,v in edges:
          G.edges[u, v]['edgeseg'].draw(fig,ax, color=color, alpha=alpha, linewidth=linewidth, markendpts_p=markendpts_p) 

#--------------------------------------------------------------------------------------------
def draw_graph_edges(fig, ax, G, alpha=0.1, linewidth=2, markendpts_p=False):
     """ Render all the graph edges by means of rendering all their associated edge segments. 
     """
     
     for u,v in G.edges():
          G.edges[u, v]['edgeseg'].draw(fig,ax, alpha=alpha, linewidth=linewidth, markendpts_p=markendpts_p) 

#---------------------------------------------------------------------------------------------

def draw_path(fig,ax, G, pathnodeidxs,  color='red', linewidth = 2, alpha=1.0):
         """ Render a path on canvas given the node ids of the path from start to target.
         """
         for n1, n2 in zip(pathnodeidxs, pathnodeidxs[1:]):
              assert G.has_node(n1) , " Graph does not have node " + str(n1)
              assert G.has_node(n2) , " Graph does not have node " + str(n2)
              p = G.node[n1]['cood'] 
              q = G.node[n2]['cood'] 
              Segment(p,q).draw(fig, ax, color=color, linewidth=linewidth, alpha=1.0, markendpts_p=True)

#------------------------------------------------------------------------------------------
def make_cone_fan( fig, ax,  obstacle_tip, dx, eps, obstsegs, 
                   left_vert, right_vert,   
                   wedge_numangles, color  ):
     """ Make the cone fan at the given obstacle tip given the bend points of the shortest path to the left 
     and right of the obstacle tip. 
     """
     first_ray                 = Ray(obstacle_tip, left_vert-obstacle_tip)
     last_ray                  = Ray(obstacle_tip , right_vert-obstacle_tip)
     cones_arms_rays, shiftval = rock_line(obstacle_tip, first_ray, last_ray, 
                                 wedge_numangles=wedge_numangles)
     midsegs=[]
     for i, raybunch in enumerate(cones_arms_rays):
          if i==0:
               raybunch = raybunch[1:] # disregard the first ray from the first bunch. This omission if fixed after this for loop
          elif i==2:
               raybunch = raybunch[:-1] # ditto
               
          for ray in raybunch:
               midsegs.append( Segment( obstacle_tip, ray_shoot(ray, obstsegs), dx=dx, eps=eps) )

     cones_arms_segs = [Segment(obstacle_tip,  left_vert, dx=dx, eps=eps)] +\
                       midsegs                             +\
                       [Segment(obstacle_tip,  right_vert, dx=dx, eps=eps)]

     cones           = [Cone(obstacle_tip, segl, segr) 
                        for (segl, segr) in zip(cones_arms_segs, cones_arms_segs[1:])]
     return cones, shiftval

#------------------------------------------------------------------------------------------
class Box:
  """ Class representing the environment the car doing lookahead must move in. 
  """
  def __init__(self, xleft, xright, ydown, yup, segs, start=None, target=None):
    self.xleft  = xleft
    self.xright = xright
    self.yup    = yup    
    self.ydown  = ydown
    self.ites   = segs # should be a list of `Segments` 
    self.start  = start
    self.target = target

  # Functions to get the different corners of the box.   
  def llc (self): # lower left corner
     return np.asarray([self.xleft, self.ydown])
     
  def ulc (self): # upper left corner
     return np.asarray([self.xleft, self.yup])

  def lrc (self): # lower right corner
     return np.asarray([self.xright, self.ydown])
     
  def urc (self): # upper right corner
     return np.asarray([self.xright, self.yup])

  def get_obstsegs (self):
       """ The obstacle segments include both the walls of the box *and* the stalactites/stalagmites. 
       """
       wall_lower = Segment(self.llc(), self.lrc()) # lower wall 
       wall_upper = Segment(self.ulc(), self.urc()) # upper wall
       wall_left  = Segment(self.llc(), self.ulc()) # left  wall
       wall_right = Segment(self.lrc(), self.urc()) # right wall
       obstsegs = list(self.ites) + [wall_lower, wall_upper, wall_left, wall_right] 
       return obstsegs

  def get_obstips(self):
        """ The obstacle tip is assumed to be the 2nd point of each stalactite/stalagmite segment
        (Should probablty make this into a separae dtaa type...???)
        """
        return [ite.q for ite in self.ites]

  def draw(self,fig, ax):
     """ Render the box ie. the walls and mites/tites along with other auxiliary informatiuon if needed. 
     """
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
     """ A class representing directed Segments. 
     """
     def __init__(self, p, q,  dx = None, eps=1e-3,):
           self.p     = np.asarray(p)
           self.q     = np.asarray(q)
           self.stpts = []

           if dx is not None:
                self.stpts = seg_steiner_pts(self.p, self.q, dx, eps=eps)

     def mkstpts(self, dx, eps):
           """ Subdivide the segment with steiner points according to given parameters. 
           Any existing subdivision is overwritten.  
           """
           self.stpts = seg_steiner_pts(self.p, self.q, dx, eps=eps)

     def draw(self,fig,ax, color='black', alpha=1.0, linewidth=0.7, drawstpts_p=True, markendpts_p=True):
          """ Draw segment as undirected as just a simple plain segment possibly with the 
          steiner points thrown in and an option to plot or not plot the endpoints of the segment 
          """
          if markendpts_p:
               marker = 'o'
          else:
               marker = ''
          ax.plot( [self.p[0],self.q[0]], [self.p[1],self.q[1]], '-', color=color, alpha=alpha, linewidth=linewidth, marker=marker)
          if drawstpts_p:
                for pt in self.stpts:
                   draw_point(fig, ax, pt, radius = 0.005, facecolor='yellow', edgecolor='yellow', alpha=alpha)

          
     def draw_directed(self,fig,ax,color='blue', alpha=0.5, drawstpts_p=True):
          """ Draw a directed segment as directed arrow from the member point p to the member point q
          """
          dx = self.q[0]-self.p[0]
          dy = self.q[1]-self.p[1]
          ax.arrow(self.p[0], self.p[1], dx,dy, head_width=0.02, 
                   edgecolor=color, 
                   facecolor=color, 
                   length_includes_head=True, alpha=alpha)
          if drawstpts_p:
                for pt in self.stpts:
                   draw_point(fig, ax, pt, radius = 0.005, facecolor='yellow')

     def tocgal(self):
          """ Convert segment to a Segment_2 instance of CGAL. 
          Useful when using geometric predicates. 
          """
          p, q = self.p, self.q
          return Segment_2(Point_2(p[0], p[1]), Point_2(q[0], q[1]))
     
     def reverse(self):
          """ This reverse the orientation of a segment, swapping p and q
          Note that no Steiner point subdivision parameters are passed i.e. dx is None by default.
          """
          return Segment(self.q, self.p)

     def __str__(self):
          return Fore.GREEN  + "Segment (" + str(self.p) + ", " + str(self.q) + ")" + Style.RESET_ALL
#------------------------------------------------------------------------------------------
class Ray :
     """ Class for a ray. Represented by a base point and a vector (which points in the same direction as the ray)
     The vector supplied need not be normalized, since the constructor does it during ....well....the construction.
     """
     def __init__(self, p, dirn):
           self.p    = np.asarray(p)
           self.dirn = np.asarray(dirn)/np.linalg.norm(dirn) # normalize for future convenience

     def tocgal(self):
          """ Convert to the CGAL Ray_2 representation. This is useful for performing intersetion tests.
          """
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
          delx =q[0]-p[0]
          dely =q[1]-p[1]
          ax.arrow(p[0], p[1], delx,dely, head_width=0.02, 
                   edgecolor=color, 
                   facecolor=color, 
                   length_includes_head=True, alpha=alpha)

     def __str__(self):
          return Fore.YELLOW  + "Ray :: (Base: " + str(self.p) + "),  (Dirn:" + str(self.dirn) + ")" + Style.RESET_ALL
#------------------------------------------------------------------------------------------
class Cone:
    def __init__(self, p, seg1, seg2):
           self.p = np.asarray(p)
           
           assert np.linalg.norm(seg1.p-seg2.p) <= 1e-6, "No cone for you. Bases of Segments  r1 and r2 must match." 
           self.seg1 = seg1
           self.seg2 = seg2
           

    def draw(self, fig, ax, edgecolor ):
           self.seg1.draw_directed(fig,ax,color=edgecolor)
           self.seg2.draw_directed(fig,ax,color=edgecolor)
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
     """ WARNING: if two segments overlap i.e. lie flush with one another, completely or partially
     the routine executes an assert stopping the code., because there is no canonical point to select for the segment segment 
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
     """ Check if a ray and Segment intersect. The intersection could be a segment or a point. 
     The routine works in either case. 
     """
     s = Segment_2(Point_2(seg.p[0], seg.p[1]), Point_2(seg.q[0]    , seg.q[1]   ))
     r = Ray_2    (Point_2(ray.p[0], ray.p[1]), Vector_2(ray.dirn[0] , ray.dirn[1] ))
     
     object = CGAL_Kernel.intersection(r, s)

     intersect_p = False
     if object.is_Segment_2() or object.is_Point_2(): 
          intersect_p = True
     return intersect_p
#------------------------------------------------------------------------------------------
def ray_seg_intersection_point(ray, seg):
     """ NOTE!: if the ray and segment overlap i.e. segment lies flush along the ray. 
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
          if  ipt is not None                   and \
              (cgalpt(ipt) != cgalpt(ray.p))    and \
              np.linalg.norm(ipt-ray.p) < dmin:
                     ptmin = ipt
                     dmin  = np.linalg.norm(ipt-ray.p)

     # My favorite assertion statement! This caught a major bug ! 
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
def seg_steiner_pts(p,q,dx, eps=1e-5):
     """ Generate Steiner points along the given segment pq
     with two consecutive points spaced dx distance apart. 
     All points are placed in the interior of the segment
     starting at (seg.p). eps represents the shift 
     wrt the initial point p. if eps = 0, then the first steiner point
     if self.p itself. 

     Also add the point q to the end of the 
     list of the steiner pts. This is necessary so that the endpoint of all rays 
     (which lie on the boundary of the free space) can participate in the routing 
     and also so that the start and target are added to the graph automatically
     as a resulty of this step.
     """
     p = np.asarray(p)
     q = np.asarray(q)
     dirn        = q-p
     dirn        = dirn/np.linalg.norm(dirn) # make direction into a unit vector
     steiner_pts = []

     ctr = 1
     while True:
         newpt =  p + (eps + (ctr-1)*dx) * dirn
         if (newpt[0] - p[0]) * (newpt[0] - q[0]) <= 0 and \
            (newpt[1] - p[1]) * (newpt[1] - q[1]) <= 0 :
              steiner_pts.append( newpt )
              ctr   +=  1
         else:
              break # no more steiner points to add!

     steiner_pts.append(q)
     return steiner_pts
#------------------------------------------------------------------------------------------
def draw_point(fig,ax, pt, radius = 0.03, facecolor='red', edgecolor='black', alpha=1.0, zorder=3):
     ptmark = mpl.patches.Circle(pt, 
                              radius    = radius, 
                              facecolor = facecolor, 
                              edgecolor = edgecolor,
                              alpha    = alpha,  
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
     the writeup. Returns a list of three lists of rays (yes a list of 3 lists!). Each sublist corresponds to 
     a list of rays inside of a wedge in the order yellow,---green---yellow as encountered
     going from the start to the destination. 
     wedge_numangles=10 refers to number of wedge angles in each  yellow wedge. the number of cones in green wedge is adjusted accordingly.
     """
     assert cgalpt(first_ray.p) == cgalpt(last_ray.p), " Base points of the first_ray and last_ray are not equal"
     wedges = [[ first_ray          , last_ray.reverse()],  # yellow wedge
               [ last_ray.reverse() , first_ray.reverse()],  # green wedge
               [ first_ray.reverse(), last_ray]]   # yellow wedge
     cones_arms_rays = []

     dtheta = None
     for i, wedge in enumerate(wedges):
          start_dirn          = wedge[0].dirn/np.linalg.norm(wedge[0].dirn)
          end_dirn            = wedge[1].dirn/np.linalg.norm(wedge[1].dirn)
          angle_between_dirns = np.arccos(np.dot(start_dirn, end_dirn))   
          
          if   i == 0 :   # first yellow wedge
               angles = np.linspace(0, angle_between_dirns, wedge_numangles+1) [:-1]# angle corresponding to last ray is dropped since it is the first ray of the next wedge
               dtheta = angles[1]-angles[0]
               assert (dtheta > 0), "dtheta should be positive"
               shiftval = int(np.ceil(np.pi/dtheta)) # the number of cones in a half-plane.....

          elif i == 1:    # green edge
               assert dtheta is not None, "dtheta should have been set to nonNull in previous iteration"
               angles = mkarr(0, angle_between_dirns, dtheta)[:-1] # angle corresponding to last ray is dropped since it is the first ray of the next wedge

          elif i == 2 :   # second yellow wedge
               angles = np.linspace(0, angle_between_dirns, wedge_numangles+1) # no angle dropped here!
               
          else:
               assert (i in [0,1,2]), " i should be in [0,1,2]"

          rays = []
          for angle in angles:
               new_dirn = None
               if np.cross(start_dirn, end_dirn) < 0: # right hand rule makes thumb point downward
                  new_dirn = rot_clock(start_dirn, angle)
               else:                                  # ditto, but  thumb now upward
                  new_dirn = rot_anticlock(start_dirn, angle)
               rays.append(   Ray( obstacle_tip , new_dirn)   )

          cones_arms_rays.append(rays)

     return cones_arms_rays, shiftval

#--------------------------------------------------------------------------------------------
def mkarr(a,b,dx):
     """ Return [a+i*dx while (a+i*dx < b) ++i] + [b]
     """
     assert (a<b), "a should be l<= than b"
     i=0;
     arr = []
     while (a + i*dx) < b:
          arr.append(a + i*dx)
          i += 1
     arr.append(b)
     return arr

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
def make_box (start, target, xleft, xright, ydown, yup, mtips, ttips):
     """ Make an environment consisting of start, target, positions of the bounding box sides
     and the positions of the tips of the stalactites and stalagmites. 
     The bases of the *ites are created automatically, as they are assumed perpendicular t
     the bounding box sides. 

     Note that the mites and tites should be alternating in position 
     (i.e. no two (t)mites come successively after one another )
     Also tites and mites x coordinates should be given in sorted order. 
     """
     start  = np.asarray(start)
     target = np.asarray(target)
     mtips  = list(map(np.asarray, mtips ))
     ttips  = list(map(np.asarray, ttips ))

     # Stalagmites (mites go up)
     mites = []                     
     for i in range(len(mtips)):
         p = np.asarray([mtips[i][0],ydown])
         q = mtips[i]
         mites.append( Segment(p,q) )

     # Stalactites (tites go down)                                                     
     tites = []
     for i in range(len(ttips)):
         p = np.asarray([ttips[i][0],yup])
         q = ttips[i]
         tites.append (Segment(p,q))

     # interleave the two lists mites and tites.
     ites = [ite for ite in itertools.chain(*itertools.izip_longest(mites, tites)) if ite is not None]

     # Create the environment box with all the information and return
     box = Box(xleft,xright,ydown,yup, ites, start=start, target=target);
     return box

#-----------------------------------------------------------------------------------------------
def make_non_zero_sol_edges_bar_plot(fig, ax, numcones, non_zero_sol_edges):
          bins     = [0 for i in range(numcones)]
          
          for (n1, n2) in non_zero_sol_edges:
               conenum       = n1[0]
               bins[conenum] += 1

          # for marking only integers on the Y axis
          yint = range(min(bins), int(math.ceil(max(bins)))+1)
          mpl.pyplot.yticks(yint)
          
          ax.bar([str(i).zfill(2) for i in range(numcones)], bins, color='maroon', width=0.4)
          ax.grid(True, ls='dotted', axis='y')
          ax.set_title("Number of non zero edges in cones", fontsize=20)
          ax.set_xlabel("Cone id", fontsize=20)
          ax.set_ylabel("Number of non-zero edges in cone", fontsize=20)



#---------------------------------------------------------------------------------
def main():
     dx           = 0.05
     eps          = 0.02
     W            = 5
     lookahead_lb = 0.5

     box = make_box(start=(0,0), target=(1,0), xleft=0, xright=1, ydown=0, yup=1, mtips=list(map(np.asarray, [ [0.2,0.7] ])), ttips=[])
     #box = make_box(start=(0,0), target=(1,1), xleft=0, xright=1, ydown=0, yup=1, mtips=list(map(np.asarray, [ [0.2,0.7] ])), ttips=list(map(np.asarray, [[0.7,0.2]])))
     #box = make_box(start=(0,0), target=(3,0.6), xleft=0, xright=3, ydown=0, yup=1, mtips=list(map(np.asarray, [[0.5,0.7], [1.2,0.5]])), ttips  = list(map(np.asarray, [(0.8,0.2), (2.5,0.2)])) )
     #box = make_box(start=(0,0), target=(3,0.6), xleft=0, xright=3, ydown=0, yup=1, mtips=list(map(np.asarray, [[0.5,0.7], [1.2,0.5],[2,0.5]])), ttips  = list(map(np.asarray, [(0.8,0.2), (1.5,0.2)  , (2.5,0.2)])) )

     # Get the coordinates of the start and target
     start  = box.start
     target = box.target

     # Set up the Matplotlib canvas. 
     fig, ax = plt.subplots()
     ax.set_aspect(1.0)
     ax.set_title("Routing scheme for generating (near optimal) arclength lookahead curves", fontsize=15)
     from matplotlib.ticker import (AutoMinorLocator, MultipleLocator)
     ax.xaxis.set_major_locator(MultipleLocator(0.1))
     ax.yaxis.set_major_locator(MultipleLocator(0.1))
     ax.grid(color = 'gray', linestyle = '--', linewidth = 0.3, alpha     = 0.5)

     # Draw the box
     box.draw(fig,ax)

     # Draw the starting and ending points
     draw_point(fig,ax,start, radius=0.02)
     draw_point(fig,ax,target, radius=0.02)

     # Perform computation to compute lookaheadcurve. since fig, ax
     # have been passedd to the routine, this will enable graphical debugging. 
     # later we can move the graphical plotting part outside the function. 
     G,  lookpathnodeids, spathnodeids, info =  conical_routing(fig, ax, box, 
                                                                dx                    = dx, 
                                                                eps                   = eps, 
                                                                W                     = W, 
                                                                lookahead_lb          = lookahead_lb, 
                                                                draw_mesh_p           = True, 
                                                                draw_lp_relax_curve_p = False)
     #....................................
     # Print summary 
     print Fore.GREEN,    "Shortest path length is ", round(info['spathlength'] ,4), Style.RESET_ALL
     print Fore.YELLOW, "Lookahead path length is " , round(info['lookpathlength'],4), Style.RESET_ALL
     print "......................................................................"
     print "LP relxation lower bound on optimal path length is : " , round(info['spathlength'],4)
     
     # Render the path and set appropriate axes information
     draw_path(fig,ax,G,lookpathnodeids, color='green', linewidth=2, alpha=1.0)
     ax.set_xlabel("Shortest Path (Red) length: "       + str(round(pathlength(G,spathnodeids   ),3))  +\
                    "\nLookahead Path (Green) length :" + str(round(pathlength(G,lookpathnodeids),3)) +\
                    "\nLP relaxation lower bound (Magenta) \underline{weighted} length :" + str(round(info['path_length_lb'],3)),\
                   fontsize=12)
     ax.set_ylabel("Lookahead Lower Bound = " + str(round(lookahead_lb,4)), fontsize=15)

     # bar chart of the non-zero edges
     fig_bar, ax_bar = plt.subplots()
     make_non_zero_sol_edges_bar_plot(fig_bar,ax_bar,info['numcones'],info['non_zero_sol_edges'])
     plt.show()


#-----------------------------------------------------------------------------    
if __name__=="__main__":
     main()
