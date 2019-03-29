from matplotlib import animation
import matplotlib.pyplot as plt
from simdat import *

x1,y1,x2,y2 = get_traj("t01t1parab")
n = size(x1)
xpot,ypot = drawpot(0,50)
print(xpot)
print(ypot)

fig = plt.figure()
ax = plt.axes(xlim=(0,50),ylim=(0,1))
line, = ax.plot([],[],lw=1)
pt1, = ax.plot([],[],"r.",markersize=10)
pt2, = ax.plot([],[],"b.",markersize=10)

def init():
    line.set_data(xpot,ypot)
    pt1.set_data(x1[0],y1[0])
    pt2.set_data(x2[0],y2[0])
    return line,pt1,pt2

def update(frame):
    pt1.set_data(x1[frame],y1[frame])
    print(x1[frame])
    pt2.set_data(x2[frame],y2[frame])
    print(x2[frame])
    return line,pt1,pt2

anim = animation.FuncAnimation(fig, update, init_func=init,frames=199, interval=20, blit=True)
anim.save('basic_animation.mp4', fps=30, extra_args=['-vcodec', 'libx264'])
plt.show()
