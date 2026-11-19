#Made by Madalyn Oliver

#This script calculates the face length and outward normal vector components for all cell faces, as well as cell area
   #^Face lengths are scaled to down to centimeters by treating every 100,000 m as 1 cm
#It then records which cells own a face, allowing us to know each cell's neighboring cells
#If a cell does not have 4 neighbors, "-1" will appear as a neighbor, allowing us to identify it as a boundary cell
#Assigns the boundary conditions for all of the boundary faces, 0 for vacuum, 1 for incoming current (panhandle top), default -1 if its an interior face
#All relevant parameters are then returned in a .npz file


import numpy as np
from pathlib import Path
from collections import defaultdict #this speeds up my dictionary creation for identifying boundary cells


#directories
PROJECT_ROOT = Path(__file__).resolve().parents[2]
PROCESSED_DIR = PROJECT_ROOT / "data" / "processed"

level = 0 #once again, manual input based on which refinement level you have for your current mesh
MESH_PARAMS_FILE = PROCESSED_DIR / f"mesh_params_ref{level}.npz"

#loading the parameters from Quads_Meshing.py
with np.load(MESH_PARAMS_FILE, allow_pickle=False) as dat:
    #^allow_pickle=False just means it won't think any of the files are malicious
    nodes = dat['nodes']
    quads = dat['quads']
    centers = dat['centers']
    cell_mat = dat['cell_mat']
    cell_source = dat['cell_source']
    cell_region = dat['cell_region']
    panhandle_edges = dat['panhandle_edges']
    #these next ones aren't used but I just had them around
    void_bool = dat['void_bool']
    tx_region = dat['tx_region'].item()
    void_region = dat['void_region'].item()

#Defining variables, arrays, and lists
m_to_cm_scaler = 100000. #every 100,000 meters will become 1 cm
n_cells = len(quads)
normals = np.zeros((n_cells, 4, 2)) #one row per cell, 4 faces, 2 coordinates, 0=x, 1=y
face_lengths = np.zeros((n_cells, 4)) #one row per cell, 4 faces
face_owners = defaultdict(list) #my auto dictionary, default value is an empty list
#^this basically makes it so i don't need to create a key every time, i can just declare a key
  # and it automatically creates it
cell_neighbors = np.full((n_cells, 4), -1, dtype=int)
 # ^ recording each cell's 4 neighbors, if -1 is leftover after iteration then it is a boundary cell
boundary_cells = [] #recording boundary cells as well in case it is needed
transport_nodes = nodes/m_to_cm_scaler #applying the scale to the node coordinates
cell_areas = np.zeros(n_cells) #getting area for each cell


#face length, normal vector, and cell area calculations + assigning cell face ownerships

#in this I use both the nodes and quads parameters. 
    # the first quad index is the cell number, the second index is the face number
    # ^ together, they hold a node label which is a single integer
    # we use this to find our node, which holds the x and y coordinates of the node
    # we do this for the next cell on the node by adding 1 on the quad face number index
    # ^ now that we have two adjacent node coordinates, we can calculate the face length
    # using "signed area" it was determined that nodes are ordered clockwise
    # in the CW orientation the outward normal is: n_x = -delta_y / face_length, n_y = delta_x / face_length

for i in range(n_cells): #iterating over each cell

    area_dblsum = 0. #area summation for each cell, it is double because we are doing the centroid cross product

    for f in range(4): #iterating over each face

        #grabbing current face node labels
        node_a = quads[i,f]
        node_b = quads[i, (f+1)%4] #using the remainder to loop back to node 0

        #getting node coordinates
        start = transport_nodes[node_a]
        stop = transport_nodes[node_b]

        cross = (start[0]*stop[1]) - (stop[0]*start[1]) #cross product which gives double area (centroid area calculation)
        area_dblsum += cross

        delta_coords = stop-start #calculates [delta_x, delta_y]

        face_length = np.sqrt((delta_coords[0]**2)+(delta_coords[1]**2)) #face length calculation

        n_x = -delta_coords[1]/face_length #outward normal x-component calculation
        n_y = delta_coords[0]/face_length #outward normal y-component calculation

        face_lengths[i, f] = face_length #recording current face length for current cell in the iteration
        normals[i, f] = [n_x, n_y] #recording outward normal vector components for the current face

        face_key = tuple(sorted((node_a, node_b))) #sorted tuple containing a face's node labels
        #^need this to be sorted because a face may have a node pair (22, 47) and another could have (47, 22)
         # ^without the sorting, these same node pairs will register as different faces & not count as neighbors

        face_owners[face_key].append((i, f)) #now it is recorded that cell i, face f owns node pair (a, b)

    #Finally, recording area of the cell, dividing by 2 because it is double area
    cell_areas[i] = 0.5 * np.abs(area_dblsum) #some cross products can be negative so we take the absolute value


# Neighbor assignment

for face_key, owners in face_owners.items():

    #If a face is interior it will have two owners, I record the ownership
     #If a face is on the boundary it has 1 owner, leaving it as -1 indicates this face is on the boundary
    if len(owners) == 2:
        cell_a, face_a = owners[0]
        cell_b, face_b = owners[1]

    #this assignment lets us know which cell is neighboring our current cell, and on which current cell's face
    cell_neighbors[cell_a, face_a] = cell_b
    cell_neighbors[cell_b, face_b] = cell_a

#Recording boundary cell IDs
for i, cell in enumerate(cell_neighbors):
    #if -1 is a cell neighbor, we know it is a boundary cell
    if -1 in cell:
        boundary_cells.append(i)

boundary_cells = np.asarray(boundary_cells, dtype=int)

#getting the boundary conditions for all cells, 0 = vacuum, 1 = incoming current (panhandle top), -1 = interior face, no BC
face_bc = np.full((n_cells, 4), -1, dtype=int)

#default all boundary faces as the vacuum BC ID, 0. remember all -1 entries in cell neighbors indicate a boundary face
#this tells you both the boundary condition and which faces are on the boundary
face_bc[cell_neighbors==-1]=0


#going through the panhandle edges and finding the face keys in face owners which match up
for edge in panhandle_edges:
    face_key = tuple(sorted(edge)) #same process as earlier, face keys hold the node label pairs, now we are finding them for the panhandle edges
    owner = face_owners.get(face_key) #returns [(cell ID, local face ID)]
    cell_ID = owner[0][0] #get cell ID, first index of face_bc
    face_ID = owner[0][1] #get face ID, whatever face index of the remaining 4 of face_bc
    face_bc[cell_ID, face_ID] = 1 #rewrite it as the incoming current BC


#Saving a .npz file with all the parameters necessary for calculations
np.savez(PROCESSED_DIR / f"calculation_params_ref{level}.npz", cell_mat=cell_mat, cell_source=cell_source,
         quads=quads, normals=normals, face_lengths=face_lengths, cell_neighbors=cell_neighbors,
         boundary_cells=boundary_cells, cell_areas=cell_areas, face_bc=face_bc)


