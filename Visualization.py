import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import os

def update(frame, lines, dots, comparison_data, fig, axs):
    artists = []
    # Save the first 30 frames and then every 10th frame
    if frame < 30 or frame % 10 == 0:
        # Check and create the "frames" directory if it doesn't exist
        if not os.path.exists('frames'):
            os.makedirs('frames')
        # Save the figure for the current frame
        fig.savefig(f'frames/frame_{frame:04d}.png')

    # Update all lines and dots for the main plot and comparison plots
    for i, (line, dot, data) in enumerate(zip(lines, dots, comparison_data)):
        line.set_ydata(data[frame, :])
        artists.append(line)
        dot.set_data(500, data[frame, 500])  # Ensure dot is always updated
        artists.append(dot)
    return artists

def LivePlot(control, constructive, destructive):
    print('Generating a live plot for the start ...')
    fig, axs = plt.subplots(3, 1, figsize=(8, 6))  # Three subplots
    fields = [control, constructive, destructive]
    colors = ['black', 'blue', 'red']  # Control, Constructive, Destructive
    labels = ['Control', 'Constructive', 'Destructive']
    lines = []
    dots = []

    # Vertical line specifications
    vline_specs = [(None, None), (485, 'green'), (475, 'black')]

    for i, (ax, field, color, label, (vline_x, vline_color)) in enumerate(zip(axs, fields, colors, labels, vline_specs)):
        ax.set_ylim((-1, 1))
        ax.set_xlim((0, len(field[0, :])))
        line, = ax.plot(field[0, :], color=color, label=label)
        lines.append(line)
        dot, = ax.plot(500, field[0, 500], 'ro')  # Add a red dot on each plot
        dots.append(dot)
        ax.title.set_text(f'{label} Field')
        ax.legend(loc="upper right")
        if vline_x:  # Add vertical line if specified
            ax.axvline(x=vline_x, color=vline_color, linestyle='--')
        if i == 2:  # Only for the last subplot
            ax.set_xlabel("n-points")

    plt.subplots_adjust(left=0.1, bottom=0.1, right=0.9, top=0.9, wspace=0.4, hspace=0.5)
    # General y-label for the whole plot
    fig.text(0.04, 0.5, 'Field Strength []', va='center', rotation='vertical')

    # Creating the animation
    anim = FuncAnimation(fig, update, frames=np.arange(0, len(control)),
                         fargs=(lines, dots, fields, fig, axs), blit=True)

    # Saving the animation
    anim.save('live_plot_comparison.mp4', fps=10, extra_args=['-vcodec', 'libx264'])

    print('Live plot successfully finished and saved as video!!')

# Assuming your data loading code here

# Load your data here, assuming the paths are properly set
control = np.loadtxt('control/Ex.txt')[:, 2000:4000]
constructive = np.loadtxt('constructive/Ex.txt')[:, 2000:4000]
destructive = np.loadtxt('destructive/Ex.txt')[:, 2000:4000]

# Now call the function like this
LivePlot(control, constructive, destructive)
