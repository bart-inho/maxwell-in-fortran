import numpy as np
import matplotlib.pyplot as plt

nt = 4000

n_figures = 4000*200

# Create a function that generate a live heat map of the field  
def LiveImshow(array):
    # function that generates a live heat map of the field
    print('Generating a live heat map for the start ...')
    plt.ion()
    fig = plt.figure(figsize=((32, 8)))
    ax = fig.add_subplot(111)
    ax.set_ylim((0, 40))
    ax.set_xlim((0, 40))
    im = ax.imshow(array, cmap='jet')
    plt.title('Ex filed')
    plt.show()
    n=0
    for i in np.arange(0,np.size(array, 1)-1 , np.size(array, 0)):
        im.set_data(array[:, n:n+40])
        fig.canvas.draw()
        fig.canvas.flush_events()
        plt.pause(0.1)
        n+=200
    print('Live heat map successfully finished !!')

hello = np.loadtxt('Hz.txt')
# hello = hello[:, 2000:4000] # Crop into the region of interest
LiveImshow(hello) # Live plot