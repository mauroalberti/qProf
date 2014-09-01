
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

  
def view_surface( geosurface ):
           
    X, Y, Z = geosurface
        
    # tripcolor plot.
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1, projection='3d')
    ax.plot_trisurf( X, Y, Z, cmap=cm.jet, linewidth=0.1 )
    ax.autoscale(enable=True, axis='both', tight=True)
    plt.show()