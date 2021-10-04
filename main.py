from geometry import *

if __name__ == "__main__":
    g = [{"cords" : [0,0.0], "brdr_type_next": BrdrType.CON}, 
         {"cords" : [0.5,0.0], "brdr_type_next": BrdrType.INS}, 
         {"cords" : [0.5,0.5], "brdr_type_next": BrdrType.INP}, 
         {"cords" : [0.3,0.5], "brdr_type_next": BrdrType.CON},
         {"cords" : [0.3,0.3], "brdr_type_next": BrdrType.CON},
         {"cords" : [0.2,0.3], "brdr_type_next": BrdrType.CON},
         {"cords" : [0.2,0.5], "brdr_type_next": BrdrType.INP},
         {"cords" : [0,0.5], "brdr_type_next": BrdrType.INP}]
    m = Mesh()
    geom = Geom(g,0,m)
    geom.split(rec_deep=2)
    #print(geom)
    #geom.show()
    m.calc(10, 50, 150, 70)
    m.show()