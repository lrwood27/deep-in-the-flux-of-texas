#Made by Madalyn Oliver

#This script creates a mesh of Texas, including an exterior convex hull that is treated as a void
#Elevation data is overlaid on the mesh and cell materials are assigned based on centroid elevation
#Population density is overlaid on the mesh and source strengths are assigned based on centroid P.D.
#All relevant parameters for calculation are exported in .npz format for other scripts to use
#Includes plot functions for the resulting mesh, elevation map, P.D. map, material map, and source map 


import numpy as np
import geopandas as gpd #used to process Texas map and create boundaries
#^^this uses the shapely module for its geometric objects!
from pathlib import Path #used for reading elevation files
import gmsh #meshing module
import meshio #used to read created mesh, get nodes, and get quad cell coordinates
import rasterio #used to read elevation files, pop density file, transform coordinates,
                    #and grab indices and data values at specific coordinates
from shapely.geometry import Polygon #using this to simplify the Texas polygon further


#modules below used for plotting
import matplotlib.pyplot as plt
from rasterio.plot import plotting_extent
from rasterio.warp import transform
from matplotlib.collections import PolyCollection
from matplotlib.colors import PowerNorm, ListedColormap, BoundaryNorm
from matplotlib.patches import Patch


states = gpd.read_file('cb_2023_us_state_500k.zip') #reading geo data for all states
tx = states.loc[states['NAME']=='Texas'].copy() #getting geo data for texas
tx_meters = tx.to_crs('EPSG:5070') #change coordinate system from degrees to meters


tx_pieces = tx_meters.explode(ignore_index=True) #breaking up the multi polygon into rows of polygons
piece_areas = tx_pieces.geometry.area #calculating the area of each polygon piece
big_piece_id = piece_areas.idxmax() #get index of piece with largest area (mainland)


tx_main = tx_pieces.loc[[big_piece_id]].copy() #table containing just mainland


tx_polygon = tx_main.geometry.iloc[0] #getting the actual polygon
#^While this is GeoPandas, it uses the Shapely module for geometry so some of the functions being used overlap
 # even though the functions are from different modules


#Polygon takes the exterior coordinates to create a single polygon, so any holes that may exist are closed
mainland = Polygon(tx_polygon.exterior) #the .exterior function grabs all boundary points from tx_polygon


#used to check polygon validity or return the largest polygon if multiple exist
def largest_polygon(geom, process):
    
    if geom.is_empty:
        return(f'{process} processing got rid of the whole polygon :(')
        
    if geom.geom_type == 'Polygon': 
        return geom
        
    if geom.geom_type == 'MultiPolygon':
        return max(geom.geoms, key=lambda p: p.area) #get the largest of the polygons
      #^^polygon objects need a specified value for "largest" so we specify the largest area using the key, 
          #p=polygon
    
    raise ValueError(f'Invalid geometry type: {geom.geom_type}') #if for some reason geometry is invalid


bay_radius = 10000. #trying to close the bay, estimated size
ridge_radius = 7000. #trying to smooth ridges, estimated size
simplify_tolerance = 10000. #allowed deviation between original geometry and simplified geometry


#Closing the Houston bay area
process_tx = mainland.buffer(bay_radius).buffer(-bay_radius)

#^equivalent form:
#expand_tx = mainland.buffer(bay_radius) <- expands polygon to fill holes, specifically the bay
#process_tx = expand_tx.buffer(-bay_radius) <- shrinks polygon back to return it to its normal size, but
                                                #with holes filled

#Grab the largest polygon after processing, check if the polygon is still valid
process_tx = largest_polygon(process_tx, 'Bay')


#Resolving the finer ridges
process_tx = process_tx.buffer(ridge_radius).buffer(-ridge_radius) #same process as the bay resolving

#Grab the largest polygon and check its validity again
process_tx = largest_polygon(process_tx, 'Ridge')


#Enclose a new polygon after processing
new_tx = Polygon(process_tx.exterior)


simple_tx = new_tx.simplify(tolerance=simplify_tolerance, preserve_topology=True) #simpler version
#tolerance defines the maximum allowed distance between the original geometry and the simplified lines
#preserve topology prevents self intersections or invalid geometries


#creating a convex hull around texas, not simplifying texas even further
hull = simple_tx.convex_hull
#^separate geometry enclosing texas


#Numbering region IDs + boundary so that cells do not overlap betw. tx & the void
tx_region = 1
void_region = 2
outer_boundary = 3


og = gpd.GeoSeries([tx_polygon], crs=tx_meters.crs,)
simple = gpd.GeoSeries([simple_tx], crs=tx_meters.crs,)
convex = gpd.GeoSeries([hull], crs=tx_meters.crs,)
#geoseries so you can see the differences between the og and simplified bounds, plus the hull shape


# Function to create an "OpenCascade" gmsh surface
# ^ OpenCascade is used for CAD applications, its a geometric modeling kernel
# This function is needed to overlay the hull with texas as separate geometries for meshing
# Specifically we no longer use gmsh.model.geo and we replace it with gmsh.model.occ (occ=opencascade)

def add_polygon_surface(polygon, h):
    
    point_list = [] #list for coordinate points to add to mesh
    line_list = [] #list for lines to make boundary
    
    coords = list(polygon.exterior.coords)[:-1] #new bounds list for each polygon, remove repeat endpoint

    for x, y in coords:
        point_entry = gmsh.model.occ.addPoint(float(x), float(y), 0.0, h) #adding points to mesh
        point_list.append(point_entry)

    for i in range(len(point_list)):
        start = point_list[i]
        end = point_list[(i+1)%len(point_list)] #remainder calc is used to loop back to zero on last point
        line_entry = gmsh.model.occ.addLine(start, end) #connecting lines to points
        line_list.append(line_entry)
        
    boundary = gmsh.model.occ.addCurveLoop(line_list) #taking lines to make a boundary
    surface = gmsh.model.occ.addPlaneSurface([boundary]) #creating region inside boundary

    return surface


#Requested mesh size in meters
h = 30000. 

gmsh.initialize() #meshing time

try:
    gmsh.model.add('tx_with_void') #adding a new mesh model of texas with a void
    
    hull_surface = add_polygon_surface(hull, h) #adding the hull surface with my requested mesh size
    tx_surface = add_polygon_surface(simple_tx, h) #adding tx surface with my requested mesh size

    fragments, input_map = gmsh.model.occ.fragment([(2, hull_surface)], [(2, tx_surface)])
         #^2=2D, this breaks up the polygons into fragments so the hull surface breaks from the tx surface
           # the input map tells us which fragments belong to each polygon, 1st index is hull, 2nd is tx

    hull_surfaces, tx_surfaces = set(), set() #recording surface numbers & subtracting tx pieces from hull
    #^we use sets for this because subtracting them will automatically remove repeat values

    #my indexing below makes more sense a little farther down where i explain entities
    #adding hull surface IDs from the input map into the hull surface set
    for surface in input_map[0]:
        hull_surfaces.add(surface[1])

    #adding texas surface IDs from the input map into the texas surface set
    for surface in input_map[1]:
        tx_surfaces.add(surface[1])

    #finally we get all of our void surfaces in the form of a set
    void_surfaces = hull_surfaces - tx_surfaces

    gmsh.model.occ.synchronize() #transfer all the geometry into gmsh's model data structure

    #2=surface/2D, sorted set of tx surfaces, using the tx_region ID assigned earlier (1)
    gmsh.model.addPhysicalGroup(2, sorted(tx_surfaces), tx_region) 
    gmsh.model.setPhysicalName(2, tx_region, 'Texas') #naming my texas group Texas, 2=2D

    #2=surface/2D, sorted set of void surfaces, using the void_region ID assigned earlier (2)
    gmsh.model.addPhysicalGroup(2, sorted(void_surfaces), void_region) 
    gmsh.model.setPhysicalName(2, void_region, 'Void') #naming my void group Void, 2=2D

    #^In the above addPhysicalGroup lines, first and second arguments together are "model entities"
    # gmsh uses model entities which consist of a dimension and a tag, where the tag is the surface #
    # so together the entity format is (dimension, tag) -> (2, surface number)

    #The entity format will be used below so I wanted to provide some explanation first :)
    
    #We now need to create our boundary group of the whole domain & combine our surface entities
    
    #Making a new combined set of all surfaces:
    all_surfaces = tx_surfaces | void_surfaces #that | operator is union, special thing for sets

    #Putting the surface tags into the entity format explained earlier:
    surface_entities = []
    
    for tag in sorted(all_surfaces):
        surface_entities.append((2, tag)) #appending surface tag in the entity format

    #now we get the combined boundary in the form of boundary entities, so the form is (1, boundary #)
    boundary_entities = gmsh.model.getBoundary(surface_entities, combined=True,
                                              oriented = False, recursive=False)
    #^combined=True means we get outside boundary of combined regions (default is True)
    #^oriented=False means there is no minus on tags if the orientation is reversed (default is True)
    #recursive=False returns 1 less dimension than input (2-1 = 1), True mentions dim=0 (default is False)

    #now we have curve tags to combine into a new group:
    outer_curves = []
    for entity in sorted(boundary_entities):
        outer_curves.append(entity[1])

    #add the outer boundary group to the mesh:
    gmsh.model.addPhysicalGroup(1, outer_curves, outer_boundary) #add group, 1=1D
    gmsh.model.setPhysicalName(1, outer_boundary, "Outer Boundary") #set name

    #setting the requested mesh size, getEntities is built in which gets all the geometric points
      #in the model, 0=0D which means points
    gmsh.model.mesh.setSize(gmsh.model.getEntities(0), h)
    
    #All setNumber options remain the same
    gmsh.option.setNumber('Mesh.ElementOrder', 1) #1=first order, one node per vertex
    gmsh.option.setNumber('Mesh.RecombineAll', 1) #1=all quads, recombine mesh into all quads
    gmsh.option.setNumber('Mesh.SubdivisionAlgorithm', 1) #1=all quads, subdivide mesh into all quads
    
    gmsh.model.mesh.generate(2) #creates the mesh, 2=2D
    gmsh.write('tx_mesh.msh') #write mesh file
    
finally: #using finally closes out the gmsh session even if there's an error
    gmsh.finalize()


#read the quadrilateral mesh
mesh = meshio.read('tx_mesh.msh')
nodes = mesh.points[:, :2] #get nodes
quads = mesh.get_cells_type('quad') #get cell info

#getting physical group labels but only for the quads, not the boundaries
cell_region = mesh.get_cell_data('gmsh:physical', "quad") 


#Boolean array to identify cell region indices which are void regions
void_bool = (cell_region==void_region)


#make sure no triangles snuck in the mesh
if len(mesh.get_cells_type('triangle')) > 0:
    raise ValueError('Mesh generation left non-quadrilateral surface cells.')


#make sure theres quads in the mesh
if len(quads) == 0:
    raise ValueError('No quadrilaterals found in tx_mesh.msh.')


#quad/cell count
print('Quadrilaterals:', len(quads))


#Cell region labels, should be the same length as quads
print('Number of region labels:', len(cell_region))


#Number of cells belonging to texas and the void
print('Texas cells:', np.count_nonzero(cell_region==tx_region))
print('Void cells:', np.count_nonzero(void_bool)) #we already made the void boolean array


#calculating area centroids of the quadrilaterals, in EPSG:5070 meters
#centroids are used to analyze elevation/population density data at a single point
#uses the general polygon centroid formula
centers = np.zeros((len(quads), 2))


for i, quad in enumerate(quads): #each quad has 4 node entries
    points = nodes[quad]
    double_area = 0.0 #calculating area except its double per the formula (aka 2*A)
    x_sum = 0.0 
    y_sum = 0.0

    for j in range(4):
        x0, y0 = points[j]
        x1, y1 = points[(j + 1) % 4] #remainder used to loop back to start point on last entry
        cross = (x0 * y1) - (x1 * y0) #cross = cross product per the formula
        double_area += cross
        x_sum += (x0 + x1) * cross
        y_sum += (y0 + y1) * cross

    #formula uses 1/6*A, this is equivalently expressed as 1/3*2*A
    centers[i] = np.array([x_sum, y_sum]) / (3.0 * double_area)


#should be the same as quads & cell_region lengths
print(len(centers))


#reading all of the elevation files
cell_elevation = np.full(len(quads), np.nan) #no number array for quads to record elevation
dem_files = sorted(Path('Elevation').glob('*.tif')) #sorted list of tif files (glob finds the files)


for file_num, dem_path in enumerate(dem_files, start=1):
    with rasterio.open(dem_path) as dem:
        #longitude and latitude of quad centers, change from meters to the elevation file units
        center_long, center_lat = transform('EPSG:5070', dem.crs, centers[:,0], centers[:,1])
        
        for i in range(len(quads)): #quad loop
            #if i already got elevation data for this quad i skip the loop
            if np.isfinite(cell_elevation[i]):
                continue
            long = center_long[i] #assign current longitude and latitude
            lat = center_lat[i]

            #skip quads that aren't in this elevation data region (geographic coordinates)
            if not (dem.bounds.left <= long < dem.bounds.right 
                    and dem.bounds.bottom < lat <= dem.bounds.top): 
                continue

            row, col = dem.index(long, lat) #long/lat indices
                
            #this next if statement is kind of redundant but it makes it so not every quad
                #is getting its indices called bc we don't have that kind of time
                #not used for population density, only using it bc i'm reading so many files

            #grab entry if current quad is in the elevation data region   
            if 0 <= row < dem.height and 0 <= col < dem.width: 
                sample = next(dem.sample([(long,lat)], indexes=1, masked=True)) #reading elevation
                #^^next is used to grab the next entry of the elevation data since we aren't looping it
                #masked=True makes values with no data get masked, used for next if statement
                
                if not np.ma.is_masked(sample[0]): #checking if there is some hidden value
                    if np.isfinite(float(sample[0])): #making sure we actually have a number
                        cell_elevation[i] = float(sample[0]) #if so, we record elevation
        #print(f'recorded: {file_num}')
        #^uncomment that if you want to see the progress
print('done')


#reading population density file
pop_file = 'usa_pd_2020_1km.tif'
cell_pop = np.full(len(quads), np.nan) #population density per cell array, initialized as nan


with rasterio.open(pop_file) as pop:
    #transforming centroid units into same units of the pop file
    pop_long, pop_lat = transform('EPSG:5070', pop.crs, centers[:,0], centers[:,1])
    
    for i in range(len(quads)):
        long, lat = pop_long[i], pop_lat[i] #current longitude and latitude
        row, col = pop.index(long, lat) #indices for long and lat

        #check if indices are within population region
        #once again less if statements because this is just one file
        if 0 <= row < pop.height and 0 <= col < pop.width:
            sample = next(pop.sample([(long,lat)], indexes=1, masked=True))

            if not np.ma.is_masked(sample[0]):
                if np.isfinite(float(sample[0])):
                    cell_pop[i]=float(sample[0])


#just if you want to see how many cells don't have data
no_pop_dats = ~np.isfinite(cell_pop) #~ returns reverse, so all non-finite entries
print(np.count_nonzero(no_pop_dats))


#We will first assign the material and source IDs as usual, then use the void boolean to
 # reassign the voids with the "void material"


#material defining
cell_mat = np.full(len(quads), -1, dtype=int) #-1 is a placeholder, int array


for i in range(len(quads)):
    elev = cell_elevation[i]

    #material 0: water
    if elev < 100:
        cell_mat[i] = 0
        
    #material 1: graphite
    elif 100 <= elev < 300:
        cell_mat[i] = 1

    #material 2: boron
    elif 300 <= elev < 700:
        cell_mat[i] = 2

    #material 3: LEU
    elif 700 <= elev < 1500:
        cell_mat[i] = 3

    #material 4: weapons-grade plutonium
    else:
        cell_mat[i] = 4


for mat in range(5):
    count = np.count_nonzero(cell_mat==mat)
    print(f'material ID: {mat}, count: {count}')


#source defining
#initial array filled with zeros because areas missing population data will be no source
cell_source = np.zeros(len(quads), dtype=int) 


for i in range(len(quads)):
    pop = cell_pop[i]
    if np.isfinite(pop): #keeps the no data entries as zero
        #source 0: no source
        if pop < 100:
            cell_source[i] = 0
        #source 1: weak source  
        elif 100 <= pop < 1000:
            cell_source[i] = 1
        #source 2: strong source
        else:
            cell_source[i] = 2


for Q in range(3):
    count = np.count_nonzero(cell_source == Q)
    print(f'source ID: {Q}, count: {count}')


#Now we overwrite the cell IDs, material 5 = void, source ID is 0 (no source)
void_mat_ID = 5
cell_mat[void_bool] = void_mat_ID #anywhere that the boolean is true becomes the void material ID
cell_source[void_bool] = 0 #same process as the material ID, but just declare them as no source


#Check what has changed with the addition of the new material ID
for mat in range(6):
    count = np.count_nonzero(cell_mat==mat)
    print(f'material ID: {mat}, count: {count}')


#Check what has changed after assigning voids as no source (0)
for Q in range(3):
    count = np.count_nonzero(cell_source == Q)
    print(f'source ID: {Q}, count: {count}')


#Writing a file with all the needed parameters that can be loaded into our other scripts
#If we imported parameters from this script it would probably run it again and we don't have that time
np.savez("mesh_params.npz", nodes=nodes, quads=quads,
         centers=centers, cell_region=cell_region, void_bool=void_bool,
         cell_mat=cell_mat, cell_source=cell_source, tx_region=tx_region, void_region=void_region)

# .npz is a zip file that stores multiple numpy arrays, you need to use np.load to read them
# below is how you grab all of the parameters:
#with np.load('mesh_params.npz', allow_pickle=False) as dat:
     #^allow_pickle=False just means it won't think any of the files are malicious
#    nodes = dat['nodes']
#    quads = dat['quads']
#    centers = dat['centers']
#    cell_mat = dat['cell_mat']
#    cell_source = dat['cell_source']


# everything below is plot functions


def mesh_map():
    fig, ax = plt.subplots(figsize=(10, 10))

    polygons = nodes[quads]

    facecolors = np.where(
        void_bool,
        "lightgray",
        "lightblue"
    )

    mesh_plot = PolyCollection(
        polygons,
        facecolors=facecolors,
        edgecolors="black",
        linewidths=0.3
    )

    ax.add_collection(mesh_plot)
    ax.autoscale_view()
    ax.set_aspect("equal")

    ax.set_xlabel("Easting (meters, EPSG:5070)")
    ax.set_ylabel("Northing (meters, EPSG:5070)")
    ax.set_title("Texas and surrounding void: quadrilateral mesh")

    ax.legend(
        handles=[
            Patch(
                facecolor="lightblue",
                edgecolor="black",
                label="Texas"
            ),
            Patch(
                facecolor="lightgray",
                edgecolor="black",
                label="Void"
            )
        ],
        loc="best"
    )

    fig.tight_layout()
    plt.show()


#exclude void region
def elev_map():
    polygons = nodes[quads]
    
    fig, ax = plt.subplots(figsize=(10, 10))
    
    image = PolyCollection(
        polygons,
        array=np.ma.masked_where(
        void_bool | ~np.isfinite(cell_elevation),
        cell_elevation
        ),
        edgecolors="none",
        cmap="terrain"
    )
    
    fig.colorbar(
        image,
        ax=ax,
        label="Sampled elevation (meters, NAVD88)"
    )
    
    ax.add_collection(image)
    ax.autoscale_view()
    ax.set_aspect("equal")
    ax.set_xlabel("Easting (meters, EPSG:5070)")
    ax.set_ylabel("Northing (meters, EPSG:5070)")
    ax.set_title("Elevation sampled at Texas quadrilateral centers")
    
    plt.show()


#exclude void region
def pop_map():
    polygons = nodes[quads]
    
    
    fig, ax = plt.subplots(figsize=(10, 10))
    ax.set_facecolor("lightgray")
    
    image = PolyCollection(
        polygons,
        array=np.ma.masked_where(
        void_bool | ~np.isfinite(cell_pop),
        cell_pop
        ),
        edgecolors="none",
        cmap="viridis",
        norm=PowerNorm(gamma=0.5, vmin=0) #power law curve to make colors more variant
    )
    
    fig.colorbar(
        image,
        ax=ax,
        label="Population density (people/km²)"
    )
    
    ax.add_collection(image)
    ax.autoscale_view()
    ax.set_aspect("equal")
    ax.set_xlabel("Easting (meters, EPSG:5070)")
    ax.set_ylabel("Northing (meters, EPSG:5070)")
    ax.set_title("Population density sampled at Texas quadrilateral centers")
    
    plt.show()


#material map
def mat_map():

    polygons = nodes[quads]
    
    material_cmap = ListedColormap([
        "royalblue", #water
        "gray",      #graphite
        "gold",      #boron
        "orange",    #LEU
        "purple",    #Plutonium
        "lightgray"  #void
    ])
    
    material_norm = BoundaryNorm(
        [-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5],
        material_cmap.N
    )
    
    fig, ax = plt.subplots(figsize=(10, 10))
    
    image = PolyCollection(
        polygons,
        array=cell_mat,
        edgecolors="none",
        cmap=material_cmap,
        norm=material_norm
    )
    
    colorbar = fig.colorbar(image, ax=ax, ticks=[0, 1, 2, 3, 4, 5])
    colorbar.ax.set_yticklabels([
        "0: Water",
        "1: Graphite",
        "2: Boron",
        "3: LEU (5%)",
        "4: Weapons-Grade Plutonium",
        "5: Void"
    ])
    
    ax.add_collection(image)
    ax.autoscale_view()
    ax.set_aspect("equal")
    ax.set_xlabel("Easting (meters, EPSG:5070)")
    ax.set_ylabel("Northing (meters, EPSG:5070)")
    ax.set_title("Prescribed material regions")
    
    plt.show()


def source_map():
    polygons = nodes[quads]
    
    source_cmap = ListedColormap([
        "lightgray",
        "orange",
        "red"
    ])
    
    source_norm = BoundaryNorm(
        [-0.5, 0.5, 1.5, 2.5],
        source_cmap.N
    )
    
    fig, ax = plt.subplots(figsize=(10, 10))
    
    image = PolyCollection(
        polygons,
        array=cell_source,
        edgecolors="none",
        cmap=source_cmap,
        norm=source_norm
    )
    
    colorbar = fig.colorbar(image, ax=ax, ticks=[0, 1, 2])
    colorbar.ax.set_yticklabels([
        "0: No source",
        "1: Weak source",
        "2: Strong source"
    ])
    
    ax.add_collection(image)
    # Texas outline
    tx_x, tx_y = simple_tx.exterior.xy

    ax.plot(
    tx_x,
    tx_y,
    color="darkgray",
    linewidth=0.5,
    zorder=10, ls='--'
    )
    ax.autoscale_view()
    ax.set_aspect("equal")
    ax.set_xlabel("Easting (meters, EPSG:5070)")
    ax.set_ylabel("Northing (meters, EPSG:5070)")
    ax.set_title("Prescribed source regions")
    
    plt.show()




