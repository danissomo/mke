
import math
from os import pread

import matplotlib
from matplotlib.lines import Line2D
import numpy as np
from enum import Enum
import matplotlib.pyplot as plt
import numpy.linalg as linalg
from numpy.lib.function_base import copy
from primitives import *

class Geom:
    def __init__(self, g, k_mater):
        self.head =Ring()
        self.head.init_geom(Point(g[-1]["cords"]), g[-1]["brdr_type_next"], 
                            Point(g[ 0]["cords"]), g[ 0]["brdr_type_next"], Point(g[1]["cords"]))
        bf = self.head
        for i in range(1, len(g)):
            bf.def_next_geom(Point(g[i]["cords"]), g[i-1]["brdr_type_next"])
            bf = bf.next
        self.head.close()
        self.points = [p["cords"] for p in g] 
        self.geom_in = []
    def __str__(self) -> str:
        res = ""
        bf = self.head
        while True:
            res += "{} {} {}\n".format(bf.prev.point.cords, bf.prevBorder.type._name_, bf.point.cords)
            bf = bf.next
            if bf is self.head:
                break
        return res


    def print(self):
        bf = self.head
        while True:
            print(bf.point.cords, bf.prevBorder.type, bf.nextBorder.type)
            bf = bf.next
            if bf is self.head:
                break
    

    def _circle_cross(self, ci1:Circle, ci2:Circle):
        v = ci1.c.cords - ci2.c.cords
        d = math.sqrt(v[0]*v[0]+ v[1]*v[1])
        a = (ci2.r*ci2.r - ci1.r*ci1.r + d*d)/(2*d)
        # if ci2.r*ci2.r - a*a < 0:
        #     return None
        h = math.sqrt(ci2.r*ci2.r - a*a) 
        p2 = ci2.c.cords + v*a/d
        p3_plus = p2 + np.array([h*v[1]/d, -h*v[0]/d])
        p3_minus = p2 - np.array([h*v[1]/d, -h*v[0]/d])
        return [Point(p3_plus), Point(p3_minus)]


    def show(self, pltshow= True):
        bf = self.head
        plt.axes()
        while True:
            line  = plt.Line2D((bf.point.cords[0], bf.next.point.cords[0]),
                               (bf.point.cords[1], bf.next.point.cords[1]) ,
                               lw =2,color='r', markersize=5, marker='.', markerfacecolor='cyan', markeredgewidth= 0)
            plt.gca().add_line(line)
            bf = bf.next
            if bf is self.head:
                break

        #plt.gca().add_line(plt.Polygon(self.points))
        plt.axis('scaled')
        if pltshow:
            plt.show()

    
    def _precalc_split(self):
        result = []
        bf = self.head
        while True:
            result.append(bf.nextBorder.param_to_kart(0.5))
            bf = bf.next
            if bf == self.head:
                return result


    def split(self, i = 0, not_splited = True):
        bf = self.head
        while not_splited:
            loc = bf.nextBorder.locate_brdrs(bf.prev.nextBorder, bf.next.nextBorder)
            if not loc:
                g = [{"cords" : bf.prev.prev.point.cords, "brdr_type_next": bf.prev.prevBorder.type},
                    {"cords" : bf.prev.point.cords, "brdr_type_next": bf.prevBorder.type}, 
                    {"cords" : bf.point.cords, "brdr_type_next": bf.nextBorder.type}, 
                    {"cords" : bf.next.point.cords, "brdr_type_next": BrdrType.NON}]
                newgeom = Geom(g, 0)
                newgeom.split(i, False)
                #newgeom.show()
                self.head.del_elem(bf.prev)
                self.head.del_elem(bf)
                bf = self.head
                continue
            
            bf= bf.next
            if bf == self.head:
                break

        bf = self.head
        new_points = []
        while True:
            new_points.append(bf.split())
            bf = bf.next.next
            if bf is self.head:
                break
        sum = np.array([0.0]*len(new_points[0].point.cords))
        ln = 0
        c_cross = []
        for point in new_points:
            c1 = Circle()
            c1.from_triangle(point.prev.point, point.prev.prev.point, point.point)
            c2 = Circle()
            c2.from_triangle(point.next.point, point.next.next.point, point.point,)
            newp = self._circle_cross(c1, c2)
            dl  = math.sqrt((point.point.cords - newp[0].cords)[0] *(point.point.cords - newp[0].cords)[0]  +  (point.point.cords - newp[0].cords)[1] *(point.point.cords - newp[0].cords)[1] )
            if dl < math.pow(0.00000001, i):
                newp  = newp[1]
            else:
                newp = newp[0]

            sum += newp.cords
            c_cross.append(newp)          
            
        sum/=len(c_cross)
        #print(len(c_cross))
        # for p in c_cross:
        #     print(p.cords)
        newp = Point(sum)
        k=0
        while k < len(new_points):
            point = new_points[k]
            g = [{"cords" : point.prev.point.cords, "brdr_type_next": point.prev.nextBorder.type}, 
                 {"cords" : point.point.cords, "brdr_type_next": BrdrType.NON}, 
                 {"cords" : newp.cords, "brdr_type_next": BrdrType.NON},
                 {"cords" : point.prev.prev.point.cords, "brdr_type_next": point.prev.prev.nextBorder.type}]
            newgeom = Geom(g, 0)
            recurce_c = 1
            if i < recurce_c: 
                newgeom.split(i+1, False)
            newgeom.show(pltshow=False)
            self.geom_in.append(newgeom)
            k+=1
        
    
    def _cross(self, line:Line):
        bf = self.head
        res = []
        while True:
            t = bf.nextBorder._cross(line)
            if t is not None and t == True:
                res.append(bf)
            bf = bf.next
            if bf == self.head:
                break
        return res
    


class Ring:
    def __init__(self, prev_ring = None) -> None:
        if prev_ring == None:
            self.point = None
            self.next = None
            self.prev= None
            self.nextBorder = None
            self.prevBorder = None
            self.closed = False
        elif self.__class__.__name__  == prev_ring.__class__.__name__:
            self.point =prev_ring.nextBorder.p1
            self.next = None
            self.prev = prev_ring
            self.prevBorder = prev_ring.nextBorder
            self.nextBorder = None
            self.closed = False
        
    def close(self):
        bf = self
        while bf.next != None:
            bf.closed = True
            bf = bf.next
        bf.closed = True
        bf.next = self
        self.prev = bf
        bf.nextBorder = self.prevBorder



    def init_geom(self, p1:Point,type1,  point:Point, type2,  p2:Point):
        self.point = point
        self.prevBorder = Border(p1, point, type1)
        self.nextBorder = Border(p2, point,type2)
    

    def def_next_geom(self, p1:Point, type):
        self.nextBorder = Border(p1, self.point, type)
        self.next = Ring(self)
        


    def split(self,):
        bf = self.next
        new_point = self.nextBorder.param_to_kart(0.5)
        self.nextBorder = Border(new_point, self.point, self.nextBorder.type)
        self.next = Ring(self)
        self.next.next = bf
        self.next.nextBorder = Border(bf.point, new_point, self.nextBorder.type)
        bf.prev = self.next
        bf.prevBorder = self.next.next.nextBorder
        return bf.prev
    

    def del_elem(self, elem):
        bf = self
        while bf!= elem:
            bf = bf.next
            if bf == self:
                return None
        bf.prev.next = bf.next
        bf.prev.nextBorder = Border(bf.next.point, bf.prev.point, BrdrType.NON)
        bf.next.prev = bf.prev
        bf.next.prevBorder = bf.prev.nextBorder

    
    def __str__(self) -> str:
        return "{prev Border %%\nnext Border%%\n point %%}".format(self.prevBorder, self.nextBorder, self.point)
    






            