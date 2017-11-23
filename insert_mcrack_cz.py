# -*- coding: utf-8 -*-
"""
Created on Sat Nov 22 21:13:31 2014

@author: nvcarvalho
"""

import numpy as np

import math as mt

import timeit    


#def findadjels_byedge_inplane(all_elements,edge,indexes):
    
#    elcount = len(all_elements)
#    commonnodes = np.zeros((elcount,4))

#    for i in range(2):
#       nodei = all_elements[:,indexes] == edge[i]   
#       commonnodes[:,i] = nodei.any(axis=1)    
      
#    commonnodesint = commonnodes.astype(int)
#    adj_el = all_elements[commonnodesint.sum(axis=1) == 2,:]
    
#    return adj_el  

#according to subroutine used for visualization
def findadjels_byedge_inplane(all_elements,edge,indexes):
    
    elcount = len(all_elements)
    commonnodes = np.zeros((elcount,4))

#    for i in range(2):
#       nodei = all_elements[:,indexes] == edge[i] +
    ind0 =  all_elements[:,indexes] == edge[0]
    ind1 =  all_elements[:,indexes] == edge[1]
    
    commonnodes = ind0+ind1    
      
      
    edge_ref = []
    commonnodesint = commonnodes.astype(int)
    index = commonnodesint.sum(axis=1) == 2
    adj_el = all_elements[commonnodesint.sum(axis=1) == 2,:]
    if len(adj_el) != 0:
        if np.all(commonnodesint[index][0] == np.array([1,1,0,0])):
            edge_ref = 0
        elif  np.all(commonnodesint[index][0] == np.array([0,1,1,0])):
            edge_ref = 1
        elif  np.all(commonnodesint[index][0] == np.array([0,0,1,1])):  
            edge_ref = 2    
        elif np.all(commonnodesint[index][0] == np.array([1,0,0,1])):
            edge_ref = 3
        
    
    return adj_el, edge_ref 


def findadjels_inplane(all_elements,elnumber):
#find adjacent elements to element elnumber  
    adj_els = np.zeros(5)
    
    el = all_elements[elnumber-1,:]
    
    elcount = len(all_elements)
    commonnodes = np.zeros((elcount,4))

    for i in range(4):
       nodei = all_elements[:,1:5] == el[i+1]   
       commonnodes[:,i] = nodei.any(axis=1)    
      
    commonnodesint = commonnodes.astype(int)
    adj_els = all_elements[commonnodesint.sum(axis=1) == 2,:]
    
    return adj_els


def update_els_tthk(all_elements,all_nodes,elno,mel,tip_val,mig,theta_i):
    #updates elemnts in the throuh-thickness direction

    #if not propagate up only, or down only depending on whether mel=1 and mel = 2 is cracked
    adj_zp = True
    mel_aux = mel
    el = elno.astype(int)-1
    #propagates towards positive z
    while adj_zp:
        if mel_aux == 1:
            if all_elements[el,15] == all_elements[el,16]: #checks if mel2 has thesame orientation as mel1
                all_elements[el,21:25] = all_elements[el,17:21]
                all_elements[el,11] = tip_val
                all_elements[el,65] = theta_i
                if mig == 1:
                    all_elements[el,25] = 10.0
            else:         
                adj_zp = False
                break
        else:
            mel_aux = 1
                
        index_adj_zpos = np.all(all_elements[el,5:9] == all_elements[:,1:5],axis=1)
        #split element (it has the same orientation)
        if index_adj_zpos.sum() != 0:
            if all_elements[index_adj_zpos,10] == 0.0: #guarantees no healing can occur
                all_elements[index_adj_zpos,10] = tip_val #ctip
                all_elements[index_adj_zpos,17:21] = all_elements[el,21:25]#edge_values
                all_elements[index_adj_zpos,65] = theta_i
                el = all_elements[index_adj_zpos,0].astype(int)-1
                if mig == 1:
                    all_elements[el,25] = 10
            else:
                #insert delamination if prior element was failed
                adj_zp = False
        else:
            adj_zp = False
                    
    adj_zn = True
    mel_aux = mel
    el = elno.astype(int)-1
    #propagates towards negative z        
    while adj_zn:

        if mel_aux == 2:
            if all_elements[el,16] == all_elements[el,15]: #checks if mel1 has the same orientation as mel2
                all_elements[el,17:21] = all_elements[el,21:25]      
                all_elements[el,10] = tip_val
                all_elements[el,65] = theta_i
                if mig == 1:
                    all_elements[el,25] = 10
            else:
                adj_zn = False
                break
        else:
            mel_aux = 2
         
         
        index_adj_zneg = np.all(all_elements[el,1:5] == all_elements[:,5:9],axis=1)
        if index_adj_zneg.sum() != 0:       
            if all_elements[index_adj_zneg,11] == 0.0: #guarantees no healing can occur           
                all_elements[index_adj_zneg,11] = tip_val #ctip
                all_elements[index_adj_zneg,21:25] = all_elements[el,17:21]#edge_values
                all_elements[index_adj_zneg,65] = theta_i
                el = all_elements[index_adj_zneg,0][0].astype(int)-1
                if mig == 1:
                    all_elements[el,25] = 10
            else:
                adj_zn = False 
        else:
            adj_zn = False
    
    return all_elements


def updateels(all_elements,all_nodes,elstoupdate,mel,mig,theta_i):
#update elements
    for el in elstoupdate:
        #updates all elements status variables (fstat del, fstat m1, fstat m2)
        mels_ind = all_elements[:,0] == el[0] #selects element to update

#       updates m1                
        if mel == 1:
#   checks if updating variable is different from 0 and (it is smaller then the current damage state (meaning it is more severe) or the element is pristine)
#   guarantees no healing can occur
            if (el[2] != 0 and  all_elements[mels_ind,10]  == 0):
                    all_elements[mels_ind,10] = el[2] #updates status variable
                    all_elements[mels_ind,17:21] = el[4:8] #updates edges
                    all_elements[mels_ind,65] = theta_i
                    # updates of adjacent elements through-thickness
                    #print 'el[0]',el[0] 
                    #if el[0] == 174:
                    #    print 'el[0]',el[0] 
                    all_elements = update_els_tthk(all_elements,all_nodes,el[0],mel,el[2],mig,theta_i)
                    #print 'all_elements[1,11]', all_elements[1,11]
#       updates m2
        if mel == 2:
            if (el[3] != 0 and all_elements[mels_ind,11]  == 0): #checks if it is not already split 
                    all_elements[mels_ind,11] = el[3] #updates status variable
                    all_elements[mels_ind,21:25] = el[4:8] #updates edges
                    all_elements[mels_ind,65] = theta_i
                    #updates adjacent elements through-thickness
                    all_elements = update_els_tthk(all_elements,all_nodes,el[0],mel,el[3],mig,theta_i)
        if mig == 1:
            all_elements[mels_ind,25] = 2.0 #migration info

    return all_elements     

def intlines(xline1,tline1,xline2,tline2):
##intersection between two lines
    
    ##finds the cross product betwen the tangents to the two lines
    vaux_den = np.cross(tline2,tline1)
    
    if np.sqrt(vaux_den.dot(vaux_den)) <= 1E-9:
        
#        intersection = np.array([])
#        t = 0.0
        intersection = np.array([1E9,1E9,1E9])
        t = 1E9
    else:
        xvector = xline1-xline2
        vaux_num = np.cross(tline2,xvector)
        
        if vaux_den[0] != 0.0:
            t = -vaux_num[0]/vaux_den[0]
        elif vaux_den[1] != 0.0:
            t = -vaux_num[1]/vaux_den[1]
        elif vaux_den[2] != 0.0:
            t = -vaux_num[2]/vaux_den[2]
             
        intersection = xline1+tline1*t
    
    return intersection,t

def intlineplane(xline,tline,nplane,xplane):
##line intersection with a plane:   
    #
    cplane = np.dot(nplane,xplane)
    #
    ###project line vector onto plane normal
    nplanevline = np.dot(nplane,tline)
    #
    if nplanevline != 0: 
        t = (np.dot(nplane,xline)-cplane)/nplanevline
        intersection = xline+tline*t
    else:
        intersection = []
  
    return intersection

def mcracklocalori(exlocal,facecentroid,facenormal):
    t = 1
    pointinex = facecentroid + exlocal*t #point along ex axis

#   find the projection of a point along the ex axis in the face plane
    intersection = intlineplane(pointinex,facenormal,facenormal,facecentroid) 
    
# vector corresponding to the projection of ex on the face plane
    exface = (intersection-facecentroid)
    
    return exface

def rotatevectorrodriguez(vect,facenormal,deg):
    vdotnface = np.dot(facenormal,vect)
    vcrossnface = np.cross(facenormal,vect)
    
    vectrot = vect*mt.cos(deg/180.0*mt.pi)+vcrossnface*mt.sin(deg/180.0*mt.pi)+facenormal*vdotnface*(1-mt.cos(deg/180.0*mt.pi))
    
    return vectrot

def normalize_vect(vect):
    #normalizes a vetor
    norm_vect = np.zeros(len(vect))
    length = np.linalg.norm(vect)
    
    norm_vect[:] = vect[:]/length
    
    return norm_vect

def normaltoface(face):
    vect1 = face[0,:] - face[1,:]
    vect2 = face[1,:] - face[2,:]
    facenormal =  normalize_vect(np.cross(vect1,vect2))
    
    return facenormal
    
def facecentroid(facenodes):
    #assumes a quad bi-linear face
    F = np.array([1/4.0,1/4.0,1/4.0,1/4.0])
    facecenter = np.dot(F,facenodes)    
    return facecenter
    
#def findedgesintersect(face,crackx,crackv):
#    #finds intersection between edges of a face and the crack
#    #returns intersection with edges and the values of the ctips
#    edge = np.zeros(4)
#    ctips = np.zeros([2,3])
#    edgetonode = np.array([[1,2],[2,3],[3,4],[4,1]])
#    count = 0 
#    t = 0
#    tol = 1E-9
#    for i in range(4):
#       xline = face[edgetonode[i,0]-1,:]
#       vline = face[edgetonode[i,1]-1,:] - face[edgetonode[i,0]-1,:]
#
#       intersection,t = intlines(xline,vline,crackx,crackv)
#       if len(intersection) != 0:
#           if t >= 0.0-tol and t <= 1.0+tol:
#        #intersection occurs whithin an edge
#               if t < 0.1-tol: 
#        #avoids a corner by redifining crcakv locally in this element
#                   crackv= (xline + vline*0.1) - crackx
#                   intersection = (xline + vline*0.1)
#                   t=0.1
#               elif t > 0.9+tol:
#                    crackv= (xline + vline*0.9) - crackx
#                    intersection = (xline + vline*0.9)
#                    t = 0.9
#        #avoids a corner by redifining crcakv locally in this element 
#               ctips[count,:] = intersection
#               edge[i] = t
#               count = count+1
#            
#    nedges = edgetonode[edge>0,:] 
#            
#    #nodes that correspond to the cracked edges
#    return edge,ctips,nedges

def findedgesintersect(face,crackx,crackv):
    #finds intersection between edges of a face and the crack
    #returns intersection with edges and the values of the ctips
    edge = np.zeros(4)
    ctips = np.zeros([2,3])
    nedges = np.zeros((2,2)).astype(int)
    edgetonode = np.array([[1,2],[2,3],[3,4],[4,1]])
    t = np.zeros(4)
    intersection = np.zeros((4,3))
    tol = 1E-9
    tol_edge = 0.15
    
    index_start_pos = 5
    index_end_pos = 5
    
    #finds the indexes to the start and end crack position
    for i in range(4):
       xline = face[edgetonode[i,0]-1,:]
       vline = face[edgetonode[i,1]-1,:] - face[edgetonode[i,0]-1,:]

       intersection[i,:],t[i] = intlines(xline,vline,crackx,crackv)
       if t[i] >= 0.0-tol and t[i] <= 1.0+tol:
           if np.all(abs(crackx -  intersection[i,:]) < tol):#this was previous coded using an equality, and was not working due to tolerances not being met
               index_start_pos = i                            #scan the coed for similar cases
           else:
               index_end_pos = i
               
    if index_start_pos == 5:  #crack start position is at the edge (e.g. crack initiated at the centroid)       
        #choose a start position based on whether there is an acceptable t
        cond1 = t[:] >= tol_edge
        cond2 = t[:] <= 1-tol_edge
        cond = cond1*cond2
        index_start_pos = np.nonzero(cond)[0]
        if len(index_start_pos) == 0: #no t is acceptable; attribute start nd end positiion at random; correct start
            cond1 = t[:] > 0.0
            cond2 = t[:] < 1.0
            cond = cond1*cond2
            index_start_pos = np.nonzero(cond)[0][0] #choose the start to be the first
            #correct the start position:
            if t[index_start_pos] >= 0.0-tol and t[index_start_pos] <= 1.0+tol:
               #intersection occurs whithin an edge
               xline = face[edgetonode[index_start_pos,0]-1,:]
               vline = face[edgetonode[index_start_pos,1]-1,:] - face[edgetonode[index_start_pos,0]-1,:]
               if t[index_start_pos] < tol_edge-tol: 
                   intersection[index_start_pos,:] = (xline + vline*tol_edge)
                   t[index_start_pos]=tol_edge
#                   print 'corrected', t
                   
               elif t[index_start_pos] > (1.0-tol_edge)+tol:
                    intersection[index_start_pos,:] = (xline + vline*(1.0-tol_edge))
                    t[index_start_pos] = 1-tol_edge
#                    print 'corrected', t
            #udpated the end position based onthe new start position:
        elif len(index_start_pos) >= 1:
            index_start_pos = np.nonzero(cond)[0][0] #choose the first (or only) index start to be crackx
            
            
        crackx = intersection[index_start_pos,:] 
        for i in range(4):
           xline = face[edgetonode[i,0]-1,:]
           vline = face[edgetonode[i,1]-1,:] - face[edgetonode[i,0]-1,:]
           
           intersection[i,:],t[i] = intlines(xline,vline,crackx,crackv)
           if t[i] >= 0.0-tol and t[i] <= 1.0+tol:
               if np.all(crackx ==  intersection[i,:]):
                   index_start_pos = i
               else:
                   index_end_pos = i

#    print 'index_start_pos', index_start_pos
#    print 't[index_end_pos]',t[index_start_pos]
#    print 'index_end_pos', index_end_pos
#    print 't[index_end_pos]',t[index_end_pos]
    
    
    
   #intersection occurs whithin an edge
    xline = face[edgetonode[index_end_pos,0]-1,:]
    vline = face[edgetonode[index_end_pos,1]-1,:] - face[edgetonode[index_end_pos,0]-1,:]   
    if t[index_end_pos] < tol_edge-tol:                              
       intersection[index_end_pos,:] = (xline + vline*tol_edge)
       t[index_end_pos]=tol_edge
       #correction needed"
       
       #recalculates starting point assuming the new end point
       xline = face[edgetonode[index_start_pos,0]-1,:]
       vline = face[edgetonode[index_start_pos,1]-1,:] - face[edgetonode[index_start_pos,0]-1,:] 
       intersectionaux,taux = intlines(xline,vline,intersection[index_end_pos,:],crackv)
       if taux >= tol_edge and taux <= 1-tol_edge:
           intersection[index_start_pos,:] = intersectionaux
           t[index_start_pos] = taux
       
   #intersection,t = intlines(xline,vline,crackx,crackv)      
    elif t[index_end_pos] > (1-tol_edge)+tol:
        #crackv= (xline + vline*0.9) - crackx
        intersection[index_end_pos,:] = (xline + vline*(1-tol_edge))
        t[index_end_pos] = (1-tol_edge)
        #correction needed
        xline = face[edgetonode[index_start_pos,0]-1,:]
        vline = face[edgetonode[index_start_pos,1]-1,:] - face[edgetonode[index_start_pos,0]-1,:] 
        intersectionaux,taux = intlines(xline,vline,intersection[index_end_pos,:],crackv)
        if taux >= tol_edge  and taux <= (1-tol_edge):
            intersection[index_start_pos,:] = intersectionaux
            t[index_start_pos] = taux
            
#avoids a corner by redifining crcakv locally in this element 
    
    ctips[0,:] = intersection[index_start_pos,:]
    ctips[1,:] = intersection[index_end_pos,:]
    vctips = ctips[1,:] - ctips[0,:]    
    
    print 'theta', np.arctan(vctips[0]/vctips[1])/np.pi*180.0    
    
    edge[index_start_pos] = t[index_start_pos] 
    edge[index_end_pos] = t[index_end_pos]        
    nedges[0,:] = edgetonode[index_start_pos,:] 
    nedges[1,:] = edgetonode[index_end_pos,:] 
    
   
    #nodes that correspond to the cracked edges
    return edge,ctips,nedges


    
def initialmcrackel(eli,cracklength,nodes,elements,mel,mig,theta_i):
#determines adjacent split elements given an element eli assumed to be split at the center 
#eli - reference element to be split
#cracklenght - lenght of the crack to be inserted - NOT USED CURENTLY
#nodes - node whithin a range that can be affected by thecrack (presently all_nodes)
#elemens - elements whithin a rane that can be affected by the crack(presently all elements)
#mel - sub elemenet into which the crack will be inserted
#mig - whether the crack results from a migration
#theta_i - delamination direction as an inplane angle

    ########NOTE!!!!!########################################
    #normal to face is going to be assumed to be n=[0,0,1]
    ########################################################

#   Initializations
    tolerance = 1E-9
    
#    elsplit = np.zeros([1,12])
    
    elsplit = np.zeros([1,4+12]) #[el_no, node_i...nodes_i+4,[delam, m1 cracked, m2 pristine],edge_i...edge_i+4,0,0,0,0]first 8 correspond to the nodes
    no_cols = len(elements[0,:])
    elsplit_el_nodes = np.zeros([1,no_cols])   
    
    countcyclemax = 1000000

    #element failed
#    elfnodes = elements[eli-1,1+4*(mel-1):5+4*(mel-1)].astype(int)   ##obtains first four nodes/tracking of mcracks is done by mel = 1 or 2 (lower or upper half ply)

    if mel == 1:
        track_nodes = np.array([1,2,3,4])  ##obtains first four nodes/tracking of mcracks is done by mel = 1 or 2 (lower or upper half ply)
    elif mel == 2:
        track_nodes = np.array([5,6,7,8])
        
    inds_el = np.array(np.hstack([0,track_nodes[:]]))
    
    ##obtains first four nodes/
    elfnodes = elements[eli-1,track_nodes].astype(int)   
 
    # element face
    face = nodes[elfnodes-1,1:4]

    # normal to face 
    #facenormal = normaltoface(face) 
    facenormal = np.array([0.0,0.0,1.0]) #assumed to be [0.0,0.0,1.0]
     

    # face center
    facecenter = facecentroid(face)   
    # local ex
    exlocal = elements[eli-1,12:15]
    exface = mcracklocalori(exlocal,facecenter,facenormal)    
    # local ex rotated acording to the orientation of the matrix

    deg = elements[eli-1,15+(mel-1)]

    #determines vector aligned with the ply orientation
    exfacerot = rotatevectorrodriguez(exface,facenormal,deg)
    
    #defines a parametric line startin at the center of the element and following the vector exfacerot
    crackx = facecenter
    crackv = exfacerot

    #finds the intersection between the line defined above and the edges of the element
    edge,ctips,nedges=findedgesintersect(face,crackx,crackv)
         
    nedgescrk = elfnodes[nedges.flatten()-1]
    nedgescrk.shape=(2,2)
       
    #computes the crack lenght from cracking the eleemen        
    alength = np.linalg.norm(ctips[1,:]-ctips[0,:])
    
    ##Checks if the element has split already (m1 or m2)
    el_prev_split = False
    if mel == 1:
        if elements[elements[:,0]==eli,10].astype(int) == 0:
            #fails element:
            elsplit = np.vstack([elsplit,np.hstack([elements[eli-1,inds_el],np.array([0,3,0]),edge,np.zeros(4)])]) #[1 1 0] = delam, m1 cracked, m2 pristine
        else:
            print 'element has split already'
            el_prev_split = True
    elif mel == 2:
        if  elements[elements[:,0]==eli,11].astype(int) == 0: 
            elsplit = np.vstack([elsplit,np.hstack([elements[eli-1,inds_el],np.array([0,0,3]),edge,np.zeros(4)])]) #[1 0 1] = delam, m1 pristine, m2 cracked        
        else:
            print 'element has split already'
            el_prev_split = True
            
      
    countcycle = 0 #avoid an inifinite loop
     
    ##fails all elements up to alenght
    elementsaux = elements.astype(int)
#    elementsaux = np.delete(elementsaux, elsplit[elsplit[:,2+(mel-1)]==3,0]-1, 0) #removes the element(s) that have been cracked (m1 or m2)    
    elementsaux = np.delete(elementsaux, elsplit[elsplit[:,6+(mel-1)]==3,0]-1, 0) #removes the element(s) that have been cracked (m1 or m2)   
    while (countcycle <= countcyclemax):
        
        eledge = np.zeros([2,5])
        countcycle = countcycle+1
        adjels = np.zeros((2,no_cols))
        ctipsnew = np.zeros([2,3])
        alengthnew = np.zeros(2)
        elfnew = np.zeros([2,5])
        edgenew =  np.zeros([2,4])
        nedgescrknew = np.zeros([2,2])
        
#        inds_el = np.array([0,1+4*(mel-1),2+4*(mel-1),3+4*(mel-1),4+4*(mel-1)])
        inds_el = np.array(np.hstack([0,track_nodes[:]]))
#        for el in elementsaux[:,inds_el]: #find adjacent elements (searches all elements that have not failed)
               
        for edg in range(2): #checks edge1 and edge2 (j=1,j=2)
            edge = nedgescrk[edg,:]
#            adjels = findadjels_byedge_inplane(elementsaux,edge,inds_el[1:])
            adjels,edge_ref_adj = findadjels_byedge_inplane(elementsaux,edge,inds_el[1:])
            split_adj,edge_ref_split  =  findadjels_byedge_inplane(elsplit,edge,inds_el[1:])
            
            if len(adjels) != 0:
            	
            	if len(split_adj):
            		#print 'split_adj', split_adj
            		index_upd_split = elsplit[:,0] == split_adj[0][0]
            	
                el = adjels[0,inds_el]
                
                #element el adjacent to edge = edg
                eledge[edg,:] = el
                #extract face
                face = nodes[eledge[edg,1:].astype(int)-1,1:4]
                
                # normal to face 
                #facenormal = normaltoface(face) 
                facenormal = np.array([0.0,0.0,1.0]) #assumed to be [0.0,0.0,1.0]

                
                # face center
                facecenter = facecentroid(face)      
                # local ex
                exface = mcracklocalori(exlocal,facecenter,facenormal)          
                # local ex rotated according to the orientation of the matrix
                exfacerot = rotatevectorrodriguez(exface,facenormal,deg)              
                # finds intersection betwen the crack and the edges of the face
                edgeaux,ctipsaux,nedgesaux = findedgesintersect(face,ctips[edg,:],exfacerot)
                
                if len(split_adj): #
                    elsplit[index_upd_split,[8+edge_ref_split]] = 1.0-edgeaux[edge_ref_adj]
                    
                nedgescrkaux= el[nedgesaux.flatten()]
                nedgescrkaux.shape=(2,2)                              

                #####remove previous tip from ctipsaux ### (first needs to be selected!!!!!)

                ctipsnew[edg,:] = np.delete(ctipsaux,0,0) 
                #obtains the edges partions position
                edgenew[edg,:] = edgeaux
                #obtains the nodes at the new tip by removing the nodes at the old tip from nedgescrkaux
                
                nedgescrknew[edg,:] = np.delete(nedgescrkaux,0,0)
                #obtains the alength corresponding to the propagation from ctip to ctipnew
                alengthnew[edg] = np.linalg.norm(ctipsnew[edg,:]-ctips[edg,:])
                
                elfnew[edg,:] = adjels[0][inds_el]
        
        #exit condition - no in-plane adjacent element found
        if np.sum(alengthnew) == 0.0: #if alengthnew is equal to zero then
            ##no adjcent elements found; no crack extesion in this increment
            break
        else: #split element
            for i in range(2):
                if alengthnew[i] != 0.0: ##if the crack tip could advance
                    if (len(elsplit[elsplit[:,0]==elfnew[i,0],2]) == 0):  #check if the element has no yet been split
                            #add to the split elements list
                        if mel == 1:
                            if elements[elements[:,0]==elfnew[i,0],10].astype(int) == 0:
                                elsplit = np.vstack([elsplit,np.hstack([elfnew[i,:],np.array([0,3,0]),edgenew[i,:],np.zeros(4)])]) #[0,3,0] = delam, m1 split, m2 pristine
                            else:
                                el_prev_split = True
                        else:
                            if elements[elements[:,0]==elfnew[i,0],11].astype(int) == 0:
                                elsplit = np.vstack([elsplit,np.hstack([elfnew[i,:],np.array([0,0,3]),edgenew[i,:],np.zeros(4)])]) #[0,0,3] = delam, m1 split, m2 split
                            else:
                                el_prev_split = True
                        ctips[i,:] = ctipsnew[i,:]
                        nedgescrk[i,:]=nedgescrknew[i,:]
                        elementsaux = elementsaux[elementsaux[:,0] != elfnew[i,0],:]
        if el_prev_split == True:
            print 'elsplit found in the split path'
            break
    
    if el_prev_split:
        elsplit = np.zeros([1,12])
        print 'elsplit found in the split path'
    if countcycle > 999999:
        print "possible infinite cycle in initialmcrackel"    
        print "otherwise countcyclemax should be increased"
    if all(elsplit[0,:] == 0.0):
        elsplit =  np.delete(elsplit,0,0)
           
    elements = updateels(elements,nodes,elsplit[:,[0,5,6,7,8,9,10,11,12,13,14,15]],mel,mig,theta_i)
    return elements,nodes
        
        
        
    
    
    




    
