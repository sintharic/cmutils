"""
-----------------------
  STL file generation
-----------------------

This is a simple "brute force" conversion resulting in large files.
Each rectangle in the rastered data is simply converted to two triangles.

The resulting stl file therefore includes:

- vertices: x-, y-, and z-coordinates of all points in the original surface  
- faces: triangles represented by the indices of their corner points  
- vectors: normals of triangle faces  

You can reduce the size of the file afterward in MeshLab using
"Filters -> Simplification: Quadratic Edge Collapse Decimation",
where "Percentage reduction" is the target mesh size relative to the original.

------------
  Examples
------------

Generate a blurred uniform random surface profile and save it to an stl file 
with an added border around it:

.. code-block:: python

    import numpy as np
    import stlutils
    import scipy.ndimage as ndi
    
    data = np.random.default_rng(1234).uniform(0, 1, (128,128))
    data = ndi.gaussian_filter(data, 2)
    stlutils.convertArray(data, 'random.stl', Lx=1., Ly=1., border=8)

Save the same array to a config file and convert it without the border around it:

.. code-block:: python

    import cmutils as cm
    
    cm.dumpConfig(data, 'random.dat')
    stlutils.convertFile('random.dat', 'random2.stl')

The resulting stl file (with border) looks like this:

.. image:: snapshot.png
    :alt: The generated stl file viewed in MeshLab.
    :align: center




---------------------
  API documentation
---------------------

"""



import numpy as np
from stl import mesh
from os import path

# global variables
nx = ny = 0



def from_file(inpath:str, flip:bool=True):
  """
  Convert 2D height topography to a 3D vertex list assuming uniform lattice spacing.

  Parameters
  ----------
  inpath : str
    Filepath to a contMech config file containing nx*ny points.
  flip : bool, optional
    Whether or not to flip the topography upside-down. Default is True.

  Returns
  -------
  vertex_list : np.ndarray
    3D vertex positions as an array of shape (nx*ny, 3).

  Warning
  -------
  These files assume that the surface is periodically repeatable!

  Warning
  -------
  flip=True is default since these files store the surface upside down!
  """

  # read config file and convert it to 3D vertex list 
  global nx,ny
  print("Reading config file. CAUTION: This assumes periodic boundaries!")

  # read file
  fid = open(inpath,"r"); line1 = fid.readline().split(); fid.close()
  nx, ny = [int(line1[0][1:]), int(line1[1])]
  vertex_list = np.loadtxt(inpath, usecols=[0,1,2], dtype=np.double)

  # print surface stats
  #minX = vertex_list[:,0].min();   maxX = vertex_list[:,0].max()
  #minY = vertex_list[:,1].min();   maxY = vertex_list[:,1].max()
  #lengthX = maxX - minX
  #lengthY = maxY - minY
  #heightZ = vertex_list[:,2].max() - vertex_list[:,2].min()
  #print(lengthX,"x",lengthY,(nx,ny),"height:",heightZ)

  # update nx,ny because of periodic repetition of first line after last line
  nx += 1; ny += 1

  # turn upside down
  if flip: vertex_list[:,2] = vertex_list[:,2].max() - vertex_list[:,2]

  #print("without foundation:",vertex_list.shape) #DEBUG
  return(vertex_list)



def from_array(array:np.ndarray, Lx:float=1, Ly:float=1, flip:bool=False):
  """
  Convert 2D height topography to a 3D vertex list assuming uniform lattice spacing.

  Parameters
  ----------
  array : np.ndarray
    Array of shape (nx,ny) containing z coordinates.
  Lx : float, optional
    Physical dimension of the topography in x direction. Default is 1.
  Ly : float, optional
    Physical dimension of the topography in y direction. Default is 1.
  flip : bool, optional
    Whether or not to flip the topography upside-down. Default is False.

  Returns
  -------
  vertex_list : np.ndarray
    3D vertex positions as an array of shape (nx*ny, 3).
  """

  global nx,ny

  # fill z data
  nx,ny = array.shape 
  vertex_list = np.zeros((nx*ny,3))
  vertex_list[:,2] = array.flatten()
  
  # fill x and y data
  dx = Lx/(nx-1)
  y = np.linspace(0,Ly,ny)
  for ix in range(nx):
    vertex_list[ix*ny : (ix+1)*ny, 0] = ix*dx
    vertex_list[ix*ny : (ix+1)*ny, 1] = y

  # turn upside down
  if flip: vertex_list[:,2] = vertex_list[:,2].max() - vertex_list[:,2]

  #print("Lx:",Lx,"\tnx:",nx,"\tdx:",dx) #DEBUG
  return(vertex_list)
  


def add_border(array:np.ndarray, border:int, flip:bool):
  """
  Add border around 2D height array assuming uniform lattice spacing.

  Parameters
  ----------
  array : np.ndarray
    Array of shape (nx,ny) containing z coordinates.
  border : int
    Width of the border in 'pixels'.
  flip : bool
    Whether or not to flip the topography upside-down.

  Returns
  -------
  bordered_array : np.ndarray
    Updated 2D array of shape (nx+2*border, ny+2*border).
  """

  global nx,ny

  nx,ny = array.shape

  if flip: border_val = array.max()
  else: border_val = array.min()

  bordered_array = border_val*np.ones((nx+2*border, ny+2*border), dtype=np.double)
  bordered_array[border:-border,border:-border] = array.reshape((nx,ny))

  nx += 2*border; ny += 2*border 

  return(bordered_array)



def add_foundation(vertex_list:np.ndarray):
  """
  Add to vertex list the vertices representing the bottom foundation of the 3D model.

  Parameters
  ----------
  vertex_list : np.ndarray
    Array of shape (nx*ny, 3).

  Returns
  -------
  new_vertices : np.ndarray
    Updated vertex list of shape ((nx+2)*(ny+2), 3).
  """

  global nx,ny
  print("Adding foundation...")

  new_vertices = np.zeros(((nx+2)*(ny+2),3))

  # add border to z data
  heightZ = vertex_list[:,2].max() - vertex_list[:,2].min()
  data = -0.15*heightZ*np.ones((nx+2,ny+2),dtype=np.double) #adjust height prefactor
  data[1:nx+1,1:ny+1] = vertex_list[:,2].reshape((nx,ny))
  new_vertices[:,2] = data.flatten()

  # add border to x data
  data = np.zeros((nx+2,ny+2))
  data[1:nx+1,1:ny+1] = vertex_list[:,0].reshape((nx,ny))
  data[0,:] = vertex_list[:,0].min();
  data[-1,:] = vertex_list[:,0].max()
  data[:,0] = data[:,1]; 
  data[:,-1] = data[:,-2]
  new_vertices[:,0] = data.flatten()

  # add border to y data
  data = np.zeros((nx+2,ny+2))
  data[1:nx+1,1:ny+1] = vertex_list[:,1].reshape((nx,ny))
  data[:,0] = vertex_list[:,1].min();
  data[:,-1] = vertex_list[:,1].max()
  data[0,:] = data[1,:]; 
  data[-1,:] = data[-2,:]
  new_vertices[:,1] = data.flatten()

  # update nx,ny
  nx += 2; ny +=2

  #print("with foundation:",new_vertices.shape) #DEBUG
  return(new_vertices)



def create_faces(foundation:bool=True):
  """
  Calculate the 2*(nx-1)*(ny-1) triangles contained in a nx*ny surface.

  Parameters
  ----------
  foundation : bool, optional
    Whether or not to add 2 triangles representing the bottom of the 3D model. Default is True.

  Returns
  -------
  faces_list : np.ndarray
    Array of shape (2*(nx-1)*(ny-1) + foundation*2, 3) containing the indices of the 3 vertices of each triangle.
  """

  global nx,ny
  print("Generating triangles...")

  nxM1 = nx-1; nyM1 = ny-1

  #TODO: can this be done more efficiently using np.einsum()?
  faces_list = np.zeros((0,3),dtype=np.uint32)
  for ix in range(nxM1):
    
    # always going through triangle corners counter-clockwise
    a = np.zeros((nyM1,3),dtype=np.uint32)
    a[:,0] = np.arange(nyM1) + ix*ny # top left
    a[:,1] = np.arange(nyM1) + (ix+1)*ny + 1 # bottom right
    a[:,2] = np.arange(nyM1) + ix*ny + 1 # top right
    faces_list = np.vstack([faces_list,a])
    a[:,2] = np.arange(nyM1) + (ix+1)*ny # bottom left
    faces_list = np.vstack([faces_list,a[:,[0,2,1]]])

  # add triangles representing bottom of foundation
  if foundation:
    faces_list = np.vstack([faces_list, np.array([0,nyM1,nx*ny-1],dtype=np.uint32)])
    faces_list = np.vstack([faces_list, np.array([0,nx*ny-1,(nx-1)*ny],dtype=np.uint32)])

  #print("faces_list:",faces_list.shape)#DEBUG
  return(faces_list)



def save_mesh(vertex_list:np.ndarray, faces_list:np.ndarray, outpath:str):
  """
  Create mesh from vertices and faces and save it to an stl file.

  Parameters
  ----------
  vertex_list : np.ndarray
    Array of shape (nx*ny, 3).
  faces_list : np.ndarray
    Array of 3-tuples of indices, where each of those 3-tuples forms a triangle in vertex_list.
  outpath : str
    Filepath to the stl output file.
  """

  print("Generating vectors...")

  result = mesh.Mesh(np.zeros(faces_list.shape[0], dtype=mesh.Mesh.dtype))
  for i, f in enumerate(faces_list):
    for j in range(3):
      result.vectors[i][j] = vertex_list[f[j],:]

  print("Saving",outpath,"...")
  result.save(outpath)



def convertFile(inpath:str, outpath:str="", norm=1, flip=True, foundation=True):
  """ 
  Create stl file from config file.

  Parameters
  ----------
  inpath : str
    Filepath to a contMech config file containing nx*ny points.
  outpath : str, optional
    Filepath to the stl output file. Default is "".
  norm : float, optional
    Factor, by which to multiply coordinates. Default is 1.
  flip : bool, optional
    Whether or not to flip the topography upside-down. Default is True.
  foundation : bool, optional
    Whether or not to add 2 triangles representing the bottom of the 3D model. Default is True.
  """

  if outpath=="": outpath = path.splitext(inpath)[0] + ".stl"
  vertices = from_file(inpath, flip)
  vertices = norm*vertices
  if foundation: vertices = add_foundation(vertices)
  faces = create_faces(foundation)
  save_mesh(vertices, faces, outpath)



def convertArray(array:np.ndarray, outpath:str, Lx:float=1, Ly:float=0, flip:bool=False, foundation:bool=True, border:int=0):
  """ 
  Create stl file from 2D numpy array.

  Parameters
  ----------
  array : np.ndarray
    Array of shape (nx,ny) containing z coordinates.
  outpath : str
    Filepath to the stl output file.
  Lx : float, optional
    Physical dimension of the topography in x direction. Default is 1.
  Ly : float, optional
    Physical dimension of the topography in y direction. Default is 0.
  flip : bool, optional
    Whether or not to flip the topography upside-down. Default is False.
  foundation : bool, optional
    Whether or not to add 2 triangles representing the bottom of the 3D model. Default is True.
  border : int, optional
    Width of the border in 'pixels'. Default is 0.
  """

  if Ly==0: Ly=Lx
  if border > 0: 
    array = add_border(array, border, flip)
    Lx *= (nx-1)/(nx-2*border-1)
    Ly *= (ny-1)/(ny-2*border-1)
  vertices = from_array(array,Lx,Ly,flip)
  if foundation: vertices = add_foundation(vertices)
  faces = create_faces(foundation)
  save_mesh(vertices, faces, outpath)
