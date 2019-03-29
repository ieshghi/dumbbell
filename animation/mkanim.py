from matplotlib import animation
import matplotlib.pyplot as plt
from simdat import *

x1,x2,y1,y2 = get_traj("/home/data/ie355/Documents/code/dumbbell/data/parab/sample_traj/t01t1parab")
n = size(x1)
xpot,ypot = drawpot(0,50)

fig = plt.figure()
ax = plt.axes(xlim=(0,50),ylim=(0,1))
line, = ax.plot([],[],lw=2)

def init():
    line.set_data(xpot,ypot)
    return ax,

def update(frame):
    line.set_data(x1[frame],y1[frame])
    line.set_data(x2[frame],y2[frame])
    return line,

anim = animation.FuncAnimation(fig, update, init_func=init,frames=200, interval=20, blit=True)
anim.save('basic_animation.mp4', fps=30, extra_args=['-vcodec', 'libx264'])
plt.show()
