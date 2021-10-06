
import math
import numpy as np
import numpy.linalg as linalg
from numpy.lib.function_base import copy
from enum import Enum

class BrdrType(Enum):
    NON = 0
    INS = 1
    CON = 2
    INP = 3


class Point:
    def __init__(self, cords:np.ndarray) -> None:
        self.cords = copy(cords)
    
    def __str__(self) -> str:
        return str(self.cords)
    

class Line:
    def __init__(self, p1:Point, p2:Point) -> None:
        self.p1 = p1
        self.p2 = p2
        self.v = p2.cords -p1.cords
        
    
    def __str__(self) -> str:
        return "Line P1%% P2%% V%%".format(self.p1, self.p2, self.v)


    def param_to_kart(self, t:float)->Point:
        return Point(self.p1.cords + (self.p2.cords - self.p1.cords)*t)
    

    def _len(self):
        v =self.p1.cords - self.p2.cords
        l = 0
        for i in range(len(v)):
            l+=v[i]*v[i]
        return math.sqrt(l)


    def _cross(self, line):
        if (self.v/line.v)[0] - (self.v/line.v)[1] < 0.0001:
            return False 
        if line.v[0] != 0 and self.v[0] != 0:
            k = np.array([line.p2.cords[0], self.p2.cords[0]]) - np.array([line.p2.cords[1]*line.v[1]/line.v[0], self.p2.cords[1]*self.v[1]/self.v[0]]) 
        elif line.v[0] == 0 and self.v[0]!=0: 
            k = np.array([line.p2.cords[0], self.p2.cords[0]]) - np.array([0.0, self.p2.cords[1]*self.v[1]/self.v[0]])
        elif line.v[0] != 0 and self.v[0] == 0:
            k = np.array([line.p2.cords[0], self.p2.cords[0]]) - np.array([line.p2.cords[1]*line.v[1]/line.v[0], 0.0]) 
        else:
            k = np.array([line.p2.cords[0], self.p2.cords[0]])
        
        try:
            m = np.array([[1, -line.v[1]/line.v[0]],[1, -self.v[1]/self.v[0]]])
            inv_m = linalg.inv(m)
        except:
            return False

        result  = np.matmul(inv_m, k)
        if  math.isnan(result[0]) or math.isnan(result[1]):
            return False
        if self.v[0] != 0:
            t1 = (result[0]-self.p2.cords[0])/self.v[0]
        else:
            t1 = (result[1]-self.p2.cords[1])/self.v[1]
        if line.v[0] != 0:
            t2 = (result[0]-line.p2.cords[0])/line.v[0]
        else:
            t2 = (result[1]-line.p2.cords[1])/line.v[1]
        
        
        if math.isinf(t1) or math.isnan(t1) or math.isinf(t2) or math.isnan(t2):
            return False
        return t1 < 1 and t2 < 1 and t1 >0  and t2 > 0


    def y_from_x(self, x):
        if self.p1.cords[0] == self.p2.cords[0]:
            return None
        m = np.array([[self.p1.cords[0], 1],[self.p2.cords[0], 1]])
        m_inv = linalg.inv(m)
        k_p = np.matmul(m_inv, np.array([self.p1.cords[1], self.p2.cords[1]]))
        return k_p[0]*x + k_p[1]

    
    def locate_brdrs(self,  b1, b2):
        if np.array_equal(self.p1.cords , b1.p1.cords):
            p1 = b1.p2
        elif np.array_equal(self.p1.cords, b1.p2.cords):
            p1 = b1.p1
        elif np.array_equal(self.p2.cords, b1.p1.cords):
            p1 = b1.p2
        elif np.array_equal(self.p2.cords, b1.p2.cords):
            p1 = b1.p1
        else:
            return None
        
        if np.array_equal(self.p1.cords , b2.p1.cords):
            p2 = b2.p2
        elif np.array_equal(self.p1.cords, b2.p2.cords):
            p2 = b2.p1
        elif np.array_equal(self.p2.cords, b2.p1.cords):
            p2 = b2.p2
        elif np.array_equal(self.p2.cords, b2.p2.cords):
            p2 = b2.p1
        else:
            return None

        y1 = self.y_from_x(p1.cords[0])
        y2 = self.y_from_x(p2.cords[0])

        if y1 == None or y2 == None:
            return (p1.cords[0] > self.p1.cords[0] and p2.cords[0] > self.p1.cords[0]) or (p1.cords[0] < self.p1.cords[0] and p2.cords[0] < self.p1.cords[0])
        
        return (y1 > p1.cords[1] and y2 > p2.cords[1]) or (y1 < p1.cords[1] and y2 < p2.cords[1]) 
            
    def locate_points(self, p1:Point, p2:Point):
        y1 = self.y_from_x(p1.cords[0])
        y2 = self.y_from_x(p2.cords[0])
        if y1 == None or y2 == None:
            return (p1.cords[0] > self.p1.cords[0] and p2.cords[0] > self.p1.cords[0]) or (p1.cords[0] < self.p1.cords[0] and p2.cords[0] < self.p1.cords[0])
        
        return (y1 > p1.cords[1] and y2 > p2.cords[1]) or (y1 < p1.cords[1] and y2 < p2.cords[1]) 




class Circle:
    def  __init__(self,  c:Point=None, r:float=None):
        self.c = c
        self.r = r
    
    def from_triangle(self, p1:Point, p2:Point, p3:Point):
        v12 = p1.cords - p2.cords
        v23 = p2.cords - p3.cords
        h12 = (p1.cords + p2.cords)/2
        h23 = (p2.cords + p3.cords)/2
        x=0
        y=1
        Y = (-h23[y]*v12[x]*v23[y]-h23[x]*v12[x]*v23[x]+h12[y]*v12[y]*v23[x] + h12[x]*v12[x]*v23[x])/(v12[y]*v23[x]-v12[x]*v23[y])
        if v23[x] != 0:
            X = v23[y]*(-Y+h23[y])/v23[x]+h23[x]
        else:
            X = h23[x]
        self.c = Point([X, Y])
        vr = self.c.cords  - p1.cords
        self.r = math.sqrt(vr[0]*vr[0]+vr[1]*vr[1])
        pass



class Border(Line):
    def __init__(self, p1: Point, p2: Point, type) -> None:
        super().__init__(p1, p2)
        self.type = type

        
    def __str__(self) -> str:
        return "{%%, type: %%}".format(super().__str__(), self.type)  
   