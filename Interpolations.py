
import math
import numpy
import scipy.spatial
import startinpy 


def bbox(list_pts_3d, jparams):
    #Define the bounding box of the output raster using giving sample points
    cellsize=jparams['cellsize'] # read cell size in json file
    x_min,x_max=min(item[0] for item in list_pts_3d),max(item[0] for item in list_pts_3d)
    ncols=math.ceil((x_max-x_min)/cellsize)
    x_left=x_min
    y_min,y_max=min(item[1] for item in list_pts_3d),max(item[1] for item in list_pts_3d)
    nrows=math.ceil((y_max-y_min)/cellsize)
    y_left=y_min
    bbox={"x_left":x_left,"y_left": y_left,"ncols":ncols,"nrows":nrows}
    return bbox

def asciiheader(bbox,jparams):
    size=jparams['cellsize']
    with open(jparams['output-file'], 'w') as fout:
        fout.write(f'NCOLS {bbox["ncols"]}\n')
        fout.write(f'NROWS {bbox["nrows"]}\n')
        fout.write(f'XLLCORNER {bbox["x_left"]}\n')
        fout.write(f'YLLCORNER {bbox["y_left"]}\n')
        fout.write(f'CELLSIZE {size}\n')
        fout.write("NODATA_VALUE -9999\n")

def grid_construction(list_pts_3d,bbox,jparams):
    size=jparams['cellsize']
    nrows=bbox['nrows']
    ncols=bbox['ncols']
    list_pts=[]
    for point in list_pts_3d:
        pt=[point[0],point[1]]
        list_pts.append(pt)
    grid_centers=[]
    for cur_row in range(nrows):
        for cur_col in range(ncols):
            x=bbox['x_left']+cur_col*size+size/2
            y=bbox['y_left']+(nrows-cur_row-1)*size+size/2
            grid_center=[x,y]
            grid_centers.append(grid_center)
    return list_pts,grid_centers

def tri_area(p1,p2,p3):
    x1,y1=p1[0],p1[1]
    x2,y2=p2[0],p2[1]
    x3,y3=p3[0],p3[1]
    s=1/2*abs(x1*(y2-y3)+x2*(y3-y1)+x3*(y1-y2))
    return s

def circumcircle(a,b,c):
    ax,ay=a[0],a[1]
    bx,by =b[0],b[1]
    cx,cy=c[0],c[1]
    D=2.0*(ax*(by-cy)+bx*(cy-ay)+cx*(ay-by))
    ux=((ax**2+ay**2)*(by-cy)+(bx**2+by**2)*(cy-ay)+(cx**2+cy**2)*(ay-by))/D
    uy=((ax**2+ay**2)*(cx-bx)+(bx**2+by**2)*(ax-cx)+(cx**2+cy**2)*(bx-ax))/D
    u=[ux,uy]
    return u

def is_out_hull(list_pts,grid_centers):
    c_hull=scipy.spatial.Delaunay(list_pts)
    is_out=scipy.spatial.Delaunay.find_simplex(c_hull,grid_centers)
    return is_out

def nn_interpolation(list_pts_3d, jparams):
    """
    Function that writes the output raster with nearest neighbour interpolation
     
    Input:
        list_pts_3d: the list of the input points (in 3D)
        jparams:     the parameters of the input for "nn"
    Output:
        the output raster file for interpolated digital elevation model
    """  

    size = jparams['cellsize']

    #calculate the bounding box using bbox function
    nn_bbox=bbox(list_pts_3d, jparams)
    nrows=nn_bbox['nrows']
    ncols=nn_bbox['ncols']

    #define ASCII header
    asciiheader(nn_bbox,jparams)

    #construct the grid
    listandgrid=grid_construction(list_pts_3d,nn_bbox,jparams)
    list_pts=listandgrid[0]
    grid_centers=listandgrid[1]
    kd = scipy.spatial.KDTree(list_pts)
    d, i = kd.query(grid_centers, k=1)
    is_out=is_out_hull(list_pts,grid_centers)
    for j in range(len(grid_centers)):
        z=list_pts_3d[i[j]][2]
        if is_out[j]==-1:
            z=-9999
        with open(jparams['output-file'],'a') as fout:
                fout.write(f'{z} ')
        if j!=0 and (j+1)%ncols==0:
            with open(jparams['output-file'],'a') as fout:
                fout.write('\n')
    print("File written to", jparams['output-file'])

def idw_interpolation(list_pts_3d, jparams):
    """

    Function that writes the output raster with IDW
     
    Input:
        list_pts_3d: the list of the input points (in 3D)
        jparams:     the parameters of the input for "idw"
    Output:
        (output file written to disk)
 
    """  
    size=jparams['cellsize']
    idw_bbox=bbox(list_pts_3d, jparams)
    nrows,ncols=idw_bbox['nrows'],idw_bbox['ncols']
    angle=jparams['angle']
    radius1,radius2=jparams['radius1'],jparams['radius2']
    a,b=max(radius1,radius2),min(radius1,radius2)

    #define ASCII header
    asciiheader(idw_bbox,jparams) 
    listandgrid=grid_construction(list_pts_3d,idw_bbox,jparams)
    list_pts=listandgrid[0]
    grid_centers=listandgrid[1]
    is_out = is_out_hull(list_pts, grid_centers)
    
    for index in range(len(grid_centers)):
        g=grid_centers[index]
        x,y=grid_centers[index][0],grid_centers[index][1]
        nearpoints=[]
        weight_heights=[]
        weights=[]
        #for each grid_center, calculate its nearest point set
        for point in list_pts_3d:
            x_p,y_p=point[0],point[1]
            in_ellipse=(math.cos(angle)*(x_p-x)+math.sin(angle)*(y_p-y))**2/(a*a)+(math.sin(angle)*(x_p-x)+math.cos(angle)*(y_p-y))**2/(b*b)
            if in_ellipse<=1:
                nearpoints.append(point)
        if jparams["max_points"]==0:
            nearpoints_in_use=nearpoints
        else:
            if jparams["max_points"]>=len(nearpoints):
                nearpoints_in_use = nearpoints
            else:
                nearpoints_in_use = nearpoints[:jparams["max_points"]]
        for p in nearpoints_in_use:
            if len(nearpoints)<jparams["min_points"]:
                z=-9999
                continue
            if len(nearpoints)==0:
                z=-9999
                continue
            distance=math.sqrt((p[0]-g[0])**2+(p[1]-g[1])**2)
            if distance==0:
                z=p[2]
                continue
            weight=distance**(-jparams["power"])
            weight_height=weight*p[2]
            weight_heights.append(weight_height)
            weights.append(weight)
            z=sum(weight_heights)/sum(weights)
        if is_out[index]==-1:
            z=-9999
        with open(jparams['output-file'],'a') as fout:
                fout.write(f'{z} ')
        if index!=0 and (index+1)% ncols == 0:
            with open(jparams['output-file'], 'a') as fout:
                fout.write('\n')
    print("File written to", jparams['output-file'])

def tin_interpolation(list_pts_3d, jparams):
    """

    Function that writes the output raster with linear in TIN interpolation
     
    Input:
        list_pts_3d: the list of the input points (in 3D)
        jparams:     the parameters of the input for "tin"
    Output:
        (output file written to disk)
 
    """  

    size=jparams['cellsize']
    tin_bbox=bbox(list_pts_3d, jparams)
    nrows=tin_bbox['nrows']
    ncols=tin_bbox['ncols']
    #define ASCII header
    asciiheader(tin_bbox,jparams)
    listandgrid=grid_construction(list_pts_3d,tin_bbox,jparams)
    list_pts=listandgrid[0]
    grid_centers=listandgrid[1]
    is_out = is_out_hull(list_pts, grid_centers)
    #construct the delaunay
    tri= scipy.spatial.Delaunay(list_pts)
    for index in range(len(grid_centers)):
        grid=grid_centers[index]
        j=tri.simplices[tri.find_simplex(grid)]
        p0=list_pts[j[0]]
        p1=list_pts[j[1]]
        p2=list_pts[j[2]]
        area=tri_area(p0,p1,p2)
        w0=tri_area(p1,p2,grid)/area
        w1=tri_area(p0,p2,grid)/area
        w2=tri_area(p0,p1,grid)/area
        z=w0*list_pts_3d[j[0]][2]+w1*list_pts_3d[j[1]][2]+w2*list_pts_3d[j[2]][2]
        if is_out[index]==-1:
            z=-9999
        with open(jparams['output-file'],'a') as fout:
                fout.write(f'{z} ')
        if index!=0 and (index+1)% ncols == 0:
            with open(jparams['output-file'], 'a') as fout:
                fout.write('\n')
    print("File written to", jparams['output-file'])

def laplace_interpolation(list_pts_3d, jparams):
    """
     
    Function that writes the output raster with Laplace interpolation
     
    Input:
        list_pts_3d: the list of the input points (in 3D)
        jparams:     the parameters of the input for "laplace"
    Output:
        (output file written to disk)
 
    """  
    
    size=jparams['cellsize']
    ll_bbox=bbox(list_pts_3d, jparams)
    nrows=ll_bbox['nrows']
    ncols=ll_bbox['ncols']
    #define ASCII header
    asciiheader(ll_bbox,jparams)
    listandgrid=grid_construction(list_pts_3d,ll_bbox,jparams)
    list_pts=listandgrid[0]
    grid_centers=listandgrid[1]
    is_out = is_out_hull(list_pts, grid_centers)
    #construct the delaunay
    t = startinpy.DT()
    t.insert(list_pts_3d)
    for index in range(len(grid_centers)):
        wi_s=[]
        wi_ai_s=[]
        grid=grid_centers[index]
        #insert grid coordinates to construct a new delaunay diagram
        grid_i=t.insert_one_pt(grid[0],grid[1],0)
        #list of incident triangles to the grid center
        itris=t.incident_triangles_to_vertex(grid_i)
        for j in range(len(itris)-1):
            tri_1=itris[j]
            p0,p1,p2=t.get_point(tri_1[0]),t.get_point(tri_1[1]),t.get_point(tri_1[2])
            tri_2=itris[j+1]
            p3,p4,p5=t.get_point(tri_2[0]),t.get_point(tri_2[1]),t.get_point(tri_2[2])
            a=circumcircle(p0,p1,p2)
            b=circumcircle(p3,p4,p5)
            edge=math.sqrt((a[0]-b[0])**2+(a[1]-b[1])**2)
            xpi=math.sqrt((p0[0]-p2[0])**2+(p0[1]-p2[1])**2)
            wi=edge/xpi
            wi_ai=wi*p2[2]
            wi_s.append(wi)
            wi_ai_s.append(wi_ai)
        z=sum(wi_ai_s)/sum(wi_s)
        if is_out[index]==-1:
            z=-9999
        with open(jparams['output-file'],'a') as fout:
                fout.write(f'{z} ')
        if index!=0 and (index+1)% ncols == 0:
            with open(jparams['output-file'], 'a') as fout:
                fout.write('\n')
        r=t.remove(grid_i)



