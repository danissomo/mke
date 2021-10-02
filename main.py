from geometry import *

if __name__ == "__main__":
    g = [{"cords" : [0,0.0], "brdr_type_next": BrdrType.GAS}, 
         {"cords" : [0.5,0.0], "brdr_type_next": BrdrType.LIQ}, 
         {"cords" : [0.5,0.5], "brdr_type_next": BrdrType.GAS}, 
         {"cords" : [0.3,0.5], "brdr_type_next": BrdrType.INS},
         {"cords" : [0.3,0.3], "brdr_type_next": BrdrType.INS},
         {"cords" : [0.2,0.3], "brdr_type_next": BrdrType.INS},
         {"cords" : [0.2,0.5], "brdr_type_next": BrdrType.INS},
         {"cords" : [0,0.5], "brdr_type_next": BrdrType.INS}]
    
    geom = Geom(g,0)
    geom.split()
    #print(geom)
    geom.show()