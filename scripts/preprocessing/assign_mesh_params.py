#Made by Madalyn Oliver

#This script takes care of source and material assignments for the current mesh refinement
#Elevation data is overlaid on the mesh and cell materials are assigned based on centroid elevation
#Population density is overlaid on the mesh and source strengths are assigned based on centroid P.D.
#All relevant parameters for calculation are exported in .npz format for other scripts to use
#Includes plot functions for the resulting mesh, elevation map, P.D. map, material map, and source map 


import numpy as np
from pathlib import Path #used for reading elevation files
import meshio #used to read created mesh, get nodes, and get quad cell coordinates
import rasterio #used to read elevation files, pop density file, transform coordinates,
                    #and grab indices and data values at specific coordinates

#modules below used for plotting
import matplotlib.pyplot as plt
from rasterio.plot import plotting_extent
from rasterio.warp import transform
from matplotlib.collections import PolyCollection
from matplotlib.colors import PowerNorm, ListedColormap, BoundaryNorm
from matplotlib.patches import Patch

PROJECT_ROOT = Path(__file__).resolve().parents[2]

RAW_DATA_DIR = PROJECT_ROOT / "data" / "raw"
ELEVATION_DIR = PROJECT_ROOT / "Elevation"

MESH_DIR = PROJECT_ROOT / "meshes" / "generated"
PROCESSED_DIR = PROJECT_ROOT / "data" / "processed"

POP_FILE = RAW_DATA_DIR / "population" / "usa_pd_2020_1km.tif"

PROCESSED_DIR.mkdir(parents=True, exist_ok=True)

level = 0 #current mesh refinement level, manual input
MESH_FILE = MESH_DIR / f"tx_mesh_ref{level}.msh"
MESH_PARAMS_FILE = PROCESSED_DIR / f"mesh_params_ref{level}.npz"

#read the quadrilateral mesh
mesh = meshio.read(MESH_FILE)
nodes = mesh.points[:, :2] #get nodes
quads = mesh.get_cells_type('quad') #get cell info, contains 4 node indices per line
lines = mesh.get_cells_type('line') #get line info, contains two node indices per line

#Reading the region IDs/tags from the mesh, tx_region = 1, void_region = 2, panhandle = 4
void_region = int(mesh.field_data["Void"][0])
tx_region = int(mesh.field_data["Texas"][0])
panhandle = int(mesh.field_data["Panhandle Boundary"][0])

#getting physical group labels for the quads
cell_region = mesh.get_cell_data('gmsh:physical', "quad") 

#also getting physical group labels for the lines to eventually grab the panhandle node/edge indices
line_tags = mesh.get_cell_data('gmsh:physical', 'line')

#Boolean array to identify cell region indices which are void regions
void_bool = (cell_region==void_region)

#Boolean array to identify line indices which are the pandhandle edges
panhandle_edges = lines[line_tags==panhandle]
#^this gets added to the .npz file so that we know the top edge node pairs


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
dem_files = sorted(ELEVATION_DIR.glob('*.tif')) #sorted list of tif files (glob finds the files)


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
cell_pop = np.full(len(quads), np.nan) #population density per cell array, initialized as nan


with rasterio.open(POP_FILE) as pop:
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

print('Before adding void regions:')
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



#Now we overwrite the cell IDs, material 5 = void, source ID is 0 (no source)
void_mat_ID = 5
cell_mat[void_bool] = void_mat_ID #anywhere that the boolean is true becomes the void material ID
cell_source[void_bool] = 0 #same process as the material ID, but just declare them as no source

print('After adding void regions:')
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
np.savez(MESH_PARAMS_FILE, nodes=nodes, quads=quads,
         centers=centers, cell_region=cell_region, void_bool=void_bool,
         cell_mat=cell_mat, cell_source=cell_source, tx_region=tx_region, 
         void_region=void_region, panhandle_edges=panhandle_edges)

# .npz is a zip file that stores multiple numpy arrays, you need to use np.load to read them
# below is how you grab all of the parameters, some have been added but you get the format:
#with np.load(MESH_PARAMS_FILE, allow_pickle=False) as dat:
     #^allow_pickle=False just means it won't think any of the files are malicious
#    nodes = dat['nodes']
#    quads = dat['quads']
#    centers = dat['centers']
#    cell_mat = dat['cell_mat']
#    cell_source = dat['cell_source']


#Plots

#view current mesh
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
