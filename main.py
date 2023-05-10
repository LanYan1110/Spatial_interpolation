
import sys
import math
import csv
import random
import json 
import time

import Interpolations

def main():
    #-- read the needed parameters from the file 'params.json' (must be in same folder)
    try:
        jparams = json.load(open('params.json'))
    except:
        print("ERROR: something is wrong with the params.json file.")
        sys.exit()
    
    #-- store the input 3D points in list
    list_pts_3d = []
    with open(jparams['input-file']) as csvfile:
        r = csv.reader(csvfile, delimiter=' ')
        header = next(r)
        for line in r:
            p = list(map(float, line)) 
            assert(len(p) == 3)
            list_pts_3d.append(p)
    
    #-- interpolations if in the params
    if 'nn' in jparams:
        start_time = time.time()
        print("=== Nearest neighbour interpolation ===")
        Interpolations.nn_interpolation(list_pts_3d, jparams['nn'])
        print("-->%ss" % round(time.time() - start_time, 2))        
    if 'idw' in jparams:
        start_time = time.time()
        print("=== IDW interpolation ===")
        Interpolations.idw_interpolation(list_pts_3d, jparams['idw'])
        print("-->%ss" % round(time.time() - start_time, 2))        
    if 'tin' in jparams:
        start_time = time.time()
        print("=== TIN interpolation ===")
        Interpolations.tin_interpolation(list_pts_3d, jparams['tin'])
        print("-->%ss" % round(time.time() - start_time, 2))        
    if 'laplace' in jparams:
        start_time = time.time()
        print("=== Laplace interpolation ===")
        Interpolations.laplace_interpolation(list_pts_3d, jparams['laplace'])
        print("-->%ss" % round(time.time() - start_time, 2))        

if __name__ == '__main__':
    main()





