#Made by Madalyn Oliver

#This script creates a mesh of Texas, including an exterior convex hull that is treated as a void
# Mesh refinement levels have been added for error evaluation
# Material and source ID assignments have been moved to a new script


import numpy as np
import geopandas as gpd #used to process Texas map and create boundaries
#^^this uses the shapely module for its geometric objects!
from pathlib import Path #used for reading elevation files
import gmsh #meshing module
from shapely.geometry import Polygon #using this to simplify the Texas polygon further
import matplotlib.pyplot as plt


# Repository paths added because of Luke
# This script is in scripts/preprocessing/Quads_Meshing.py.
PROJECT_ROOT = Path(__file__).resolve().parents[2]

RAW_DATA_DIR = PROJECT_ROOT / "data" / "raw"

MESH_DIR = PROJECT_ROOT / "meshes" / "generated"

BOUNDARY_FILE = (
    RAW_DATA_DIR / "boundaries" / "cb_2023_us_state_500k.zip"
)

# Create output mesh folder if needed
MESH_DIR.mkdir(parents=True, exist_ok=True)


states = gpd.read_file(BOUNDARY_FILE) #reading geo data for all states

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
#Panhandle ID has been added because I'm adding this to the mesh for later BC use
tx_region = 1
void_region = 2
outer_boundary = 3
panhandle = 4


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


#Requested mesh size in meters for course mesh
h = 30000. 

#Amount of refinement levels for adding finer meshes
#Refined meshes should multiply the amount of cells in the previous mesh by 4
max_refinement_lvl = 0

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
    
    #panhandle curve is curve # 42 from outer curves, did some temporary checks and that's the only panhandle curve
    #^referring to the VERY TOP, the flat edge up against Oklahoma where we shoot the beam
    panhandle_curve = [42] #gmsh needs this in list form even though it's just 1 curve

    #getting rid of curve 42 from the outer curves list so the groups don't overlap
    outer_curves.remove(42) 
            
    #add the outer boundary group to the mesh:
    gmsh.model.addPhysicalGroup(1, outer_curves, outer_boundary) #add group, 1=1D
    gmsh.model.setPhysicalName(1, outer_boundary, "Outer Boundary") #set name

    #recording panhandle curve for applying the boundary conditions later on
    gmsh.model.addPhysicalGroup(1, panhandle_curve, panhandle)
    gmsh.model.setPhysicalName(1, panhandle, 'Panhandle Boundary')


    #setting the requested mesh size, getEntities is built in which gets all the geometric points
      #in the model, 0=0D which means points
    gmsh.model.mesh.setSize(gmsh.model.getEntities(0), h)
    
    #All setNumber options remain the same
    gmsh.option.setNumber('Mesh.ElementOrder', 1) #1=first order, one node per vertex
    gmsh.option.setNumber('Mesh.RecombineAll', 1) #1=all quads, recombine mesh into all quads
    gmsh.option.setNumber('Mesh.SubdivisionAlgorithm', 1) #1=all quads, subdivide mesh into all quads
    
    gmsh.model.mesh.generate(2) #creates the mesh, 2=2D

    #course mesh file will be written separately, for more refinements the next loop is used
    gmsh.write(str(MESH_DIR/'tx_mesh_ref0.msh')) 

    #Refinement loop, new .msh files written for each refinement level
    for lvl in range(max_refinement_lvl + 1):

        if lvl > 0:
            #Each time you call this, the mesh gets more refined (hence the loop)
            gmsh.model.mesh.refine() 

        gmsh.write(str(MESH_DIR/f'tx_mesh_ref{lvl}.msh')) #write mesh file
    
    
finally: #using finally closes out the gmsh session even if there's an error
    gmsh.finalize()

