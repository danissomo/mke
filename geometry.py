
import math
from operator import le
from os import pread

import matplotlib
from matplotlib.lines import Line2D
import numpy as np
from enum import Enum
import matplotlib.pyplot as plt
import numpy.linalg as linalg
from numpy.lib.function_base import copy, select
from primitives import *

class Mesh:
    def __init__(self) -> None:
        self.mesh ={}
        pass
    
    def add_rel(self, key,  value):

        if not key in self.mesh:
            self.mesh[key] = {value}
            return
        self.mesh[key].add(value)
    
    def show(self):
        plt.axes()
        items = self.mesh.items()
        for  pair in items:
            for value in pair[1]:
                line  = plt.Line2D((value[0][0], pair[0][0]),
                                (value[0][1], pair[0][1]) ,
                                lw =2,color='r', markersize=5, marker='.', markerfacecolor=(0.1, 0.2, 0.5), markeredgewidth= 0)
                #plt.gca().add_line(line)
         
        for point_cords in self.mesh.keys():
            t_point = self.temps[self.help_dict[point_cords]]
            high =  max(self.temps)
            c = plt.Circle(point_cords, 0.02, color = (t_point/high,0,0))
            plt.gca().add_patch(c)
        plt.axis('scaled')
        plt.show()
    

    def _dict_to_np(self, keys,  arr):
        help_dict = {}
        for i, key in zip(range(len(keys)), keys):
            help_dict[key] = i
        result = np.array([[0.0]*len(keys)]* len(keys))
        for row in arr.items():
            r_i = row[1].items()
            
            for val in r_i:
                result[help_dict[row[0]]][help_dict[val[0]]] = val[1]
        return result

    def _dict_vec_to_np(self, keys, vec):
        help_dict = {}
        for i, key in zip(range(len(keys)), keys):
            help_dict[key] = i
        result = np.array([0.0]*len(keys))
        for pair in vec.items():
            result[help_dict[pair[0]]] = pair[1]
        self.help_dict = help_dict
        return result


    def calc(self, alpha, temp, q, lambd):
        items = self.mesh.items()
        row = {}
        k_matr = {}
        for pair in items:
            key = pair[0]
            row[key] = 0.0
        for pair in items:
            key = pair[0]
            k_matr[key] = row.copy()
        vec = row.copy()
    
        i_point = 0
        i_btype = 1

        for pair_row in items:
            for_point = pair_row[0]
            values = pair_row[1]
            for contact_point in values:
                if contact_point[i_btype] == BrdrType.CON:
                    vec[for_point] += alpha*100*self.len_pp(for_point, contact_point[i_point])*temp
                elif contact_point[i_btype] == BrdrType.INP:
                    vec[for_point] += q*100*self.len_pp(for_point, contact_point[i_point])
                

        for pair_row in items:
            for_point = pair_row[0]
            values = pair_row[1]
            for contact_point in values:
                if contact_point[i_btype] == BrdrType.CON:
                    k_matr[for_point][for_point] +=alpha*100*self.len_pp(for_point, contact_point[i_point])
                    k_matr[contact_point[i_point]][contact_point[i_point]] += alpha*100*self.len_pp(for_point, contact_point[i_point])
                
                cur_lambd = lambd/(100 * self.len_pp(for_point, contact_point[i_point]))
                k_matr[for_point][for_point] += cur_lambd
                k_matr[contact_point[i_point]][contact_point[i_point]] += cur_lambd
                k_matr[for_point][contact_point[i_point]] -= cur_lambd
                k_matr[contact_point[i_point]][for_point] -= cur_lambd
                
        
        file = open("m_before.txt", "w")
        for row in k_matr.items():
            row_items = row[1].items()
            for val in row_items:
                file.write(val[1].__str__()+",")
            file.write("\n")
        
        k_matr = self._dict_to_np(self.mesh.keys(),k_matr)
        vec= self._dict_vec_to_np(self.mesh.keys(), vec)
        file = open("m.txt", "w")
        for row in k_matr:
            for val in row:
                file.write(val.__str__()+",")
            file.write("\n")
        #print(k_matr[0:][0:10])
        k_matr =  np.array(k_matr)
        m_inv = linalg.inv(k_matr)
        t_answ = np.matmul(m_inv, vec)
        print(t_answ)
        self.temps = t_answ
        return t_answ
    

    def len_pp(self, p1, p2):
        a1 = p1[0] - p2[0]
        a2 = p1[1] - p2[1]
        return math.sqrt(a1*a1 + a2*a2)

class Geom:
    def __init__(self, g, k_mater, mesh = Mesh()):
        self.g = g
        self.head =Ring()
        self.mesh = mesh
        self.head.init_geom(Point(g[-1]["cords"]), g[-1]["brdr_type_next"],
                            Point(g[ 0]["cords"]), g[ 0]["brdr_type_next"], Point(g[1]["cords"]))
        bf = self.head
        for i in range(1, len(g)):
            bf.def_next_geom(Point(g[i]["cords"]), g[i-1]["brdr_type_next"])
            bf = bf.next
        self.head.close()
        self.points = [p["cords"] for p in g]
        self.geom_in = []
        self.is_convex = False
        self.size = len(g)
        self.end_elem = False
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
        if pltshow:
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


    def make_convex_shape(self, rec_deep):#придумать структуру хранения отрезанных элементов
        bf = self.head
        while not self.is_convex:
            loc = bf.nextBorder.locate_brdrs(bf.prev.nextBorder, bf.next.nextBorder)
            if not loc:
                g = []
                back_roller = bf.next
                while True:
                    g.append({"cords" : back_roller.point.cords, "brdr_type_next": back_roller.prevBorder.type})
                    #print(back_roller.point)
                    back_roller = back_roller.prev
                    if not bf.nextBorder.locate_points(bf.prev.point, back_roller.point) and back_roller != bf and back_roller !=bf.next:
                        break
                back_roller = back_roller.next.next
                #print("end back")
                while True:
                    if back_roller == self.head:
                        self.head = self.head.prev
                        
                    self.head.del_elem(back_roller)
                    #print(back_roller.point)
                    if back_roller == bf:
                        self.head.del_elem(back_roller)
                        break
                    back_roller=back_roller.next
                
                g[-1]["brdr_type_next"] = BrdrType.NON
                #print("--")
                newgeom = Geom(g, 0, self.mesh)
                newgeom.is_convex = True
                newgeom.split(0, rec_deep)
                #self.head.del_elem(bf)
                bf = self.head
                continue
                
            bf= bf.next
            if bf == self.head:
                self.is_convex = True
                break


    def split(self, i = 0, rec_deep = 2):
        self.make_convex_shape(rec_deep)

        bf = self.head
        new_points = []
        sum = np.array([0.0]*2)
        while True:
            new_points.append(bf.split())
            sum += new_points[-1].point.cords
            bf = bf.next.next
            if bf is self.head:
                break

        sum/=len(new_points)
        newp = Point(sum)

        for point in new_points:
            g = [{"cords" : point.prev.prev.point.cords, "brdr_type_next": point.prev.prev.nextBorder.type},
                {"cords" : point.prev.point.cords, "brdr_type_next": point.prev.nextBorder.type},
                 {"cords" : point.point.cords, "brdr_type_next": BrdrType.NON},
                 {"cords" : newp.cords, "brdr_type_next": BrdrType.NON}]
            newgeom = Geom(g, 0, self.mesh)

            if i < rec_deep:
                newgeom.split(i+1, rec_deep)
            else: 
                newgeom.end_elem = True
                newgeom.add_to_mesh()
            #newgeom.show(pltshow=False)
            self.geom_in.append(newgeom)


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

    def add_to_mesh(self):
        bf = self.head
        while True:
            self.mesh.add_rel(tuple(bf.point.cords),  (tuple(bf.next.point.cords), bf.nextBorder.type))
            self.mesh.add_rel(tuple(bf.point.cords), (tuple(bf.prev.point.cords), bf.nextBorder.type))
            bf = bf.next
            if bf == self.head:
                break
        pass



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



    def split(self):
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






