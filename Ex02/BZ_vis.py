import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

plt.style.use('tableau-colorblind10')

def get_brillouin_zone_3d(cell):
    """
    Generate the Brillouin Zone of a given cell. The BZ is the Wigner-Seitz cell
    of the reciprocal lattice, which can be constructed by Voronoi decomposition
    to the reciprocal lattice.  A Voronoi diagram is a subdivision of the space
    into the nearest neighborhoods of a given set of points. 

    https://en.wikipedia.org/wiki/Wigner%E2%80%93Seitz_cell
    https://docs.scipy.org/doc/scipy/reference/tutorial/spatial.html#voronoi-diagrams
    """

    cell = np.asarray(cell, dtype=float)
    assert cell.shape == (3, 3)

    px, py, pz = np.tensordot(cell, np.mgrid[-1:2, -1:2, -1:2], axes=[0, 0])
    points = np.c_[px.ravel(), py.ravel(), pz.ravel()]

    from scipy.spatial import Voronoi
    vor = Voronoi(points)

    bz_facets = []
    bz_ridges = []
    bz_vertices = []

    for pid, rid in zip(vor.ridge_points, vor.ridge_vertices):
        if(pid[0] == 13 or pid[1] == 13):
            bz_ridges.append(vor.vertices[np.r_[rid, [rid[0]]]])
            bz_facets.append(vor.vertices[rid])
            bz_vertices += rid

    bz_vertices = list(set(bz_vertices))

    return vor.vertices[bz_vertices], bz_ridges, bz_facets

# fcc as used by QE
cell = 0.5*np.array([[-1, 0, 1],
                     [ 0, 1, 1],
                     [-1, 1, 0]])

# basis vectors of the reciprocal lattice
icell = np.linalg.inv(cell).T

# normalising basis vectors
b1, b2, b3 = np.linalg.norm(icell, axis=1)

# getting the info on the BZ
v, e, f = get_brillouin_zone_3d(icell)

# get vectors to symmetry points
G = np.dot(np.array([0    ,0    ,0    ]),icell)
L = np.dot(np.array([?.???,?.???,?.???]),icell)
X = np.dot(np.array([?.???,?.???,?.???]),icell)
U = np.dot(np.array([?.???,?.???,?.???]),icell)

# combine into a path-array
path = np.array([L,G,X,U,G])

# Now make a picture of the BZ including the k-path
fig = plt.figure()
ax = plt.subplot(111, projection='3d')
for xx in e:
    ax.plot(xx[:, 0], xx[:, 1], xx[:, 2], color='C00', lw=1.0)
ax.plot(path[:,0],path[:,1],path[:,2],'C03',ls='--')
for i in range(len(path)):
    ax.text(*path[i],np.array(['L','$\\Gamma$','X','U','$\\Gamma$'])[i])
ax.set_aspect('equal', adjustable='box')
ax.set_axis_off()
#with this line, you rotate the BZ (elev and azim) and the angle it is viewed from (roll)
#ax.view_init(elev=60, azim=70, roll=15)
fig.savefig('BZ_path_fcc.png',bbox_inches='tight',dpi=350)
plt.show()




