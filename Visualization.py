import numpy as np
import matplotlib.pyplot as plt

def LivePlot(array):
    # function that generates a live plot of the field
    print('Generating a live plot for the start ...')
    plt.ion()
    fig = plt.figure(figsize=((32, 8)))
    ax = fig.add_subplot(111)
    ax.set_ylim((-1,1))
    ax.set_xlim((0, len(array[0, :])))
    line, = ax.plot(array[0, :])
    plt.title('Ex filed')
    plt.show()
    for i in np.arange(0, len(array)):
        line.set_ydata(array[i, :])
        fig.canvas.draw()
        fig.canvas.flush_events()
        plt.pause(0.0001)
    print('Live plot successfully finished !!')

hello = np.loadtxt('Ex.txt')
# hello = hello[:, 2000:4000] # Crop into the region of interest
LivePlot(hello) # Live plot