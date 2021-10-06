from geometry import *

def example(A, B, a, b):
     g = [{"cords" : [0,0.0], "brdr_type_next": BrdrType.CON}, 
         {"cords" : [A,0.0], "brdr_type_next": BrdrType.INS}, 
         {"cords" : [A,B], "brdr_type_next": BrdrType.INP}, 
         {"cords" : [A- (A-a)/2,B], "brdr_type_next": BrdrType.CON},
         {"cords" : [A- (A-a)/2,B-b], "brdr_type_next": BrdrType.CON},
         {"cords" : [(A-a)/2,B-b], "brdr_type_next": BrdrType.CON},
         {"cords" : [(A-a)/2,B], "brdr_type_next": BrdrType.INP},
         {"cords" : [0,B], "brdr_type_next": BrdrType.INP}]
     return g

if __name__ == "__main__":
     g = example(15, 10, 5, 4)
     m = Mesh()
     geom = Geom(g,0,m)
     geom.split(rec_deep=4)
     m.calc_temps(10,22, 150, 75)
     geom.show()