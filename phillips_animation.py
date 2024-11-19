import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import matplotlib.animation as animation

# Following https://climate-cms.org/posts/2019-09-03-python-animation.html

d = xr.open_dataset('phillips_model.nc')

fig, axes = plt.subplots()

# For pressure
# var = d['ps']
# p = var[0][::-1].plot(vmin=-70,vmax=-10)
t500 = d['t500']
ps = d['ps']
time = d['time']
p_T = t500[0][::-1].plot(vmin=-30, vmax=30, add_colorbar=False, cmap='coolwarm')
p_ps = ps[0][::-1].plot.contour(add_colorbar=False,levels=np.arange(-65,-25,7))
# Wind
# var = d['u']
# p = var[0,0].plot(vmin=-10,vmax=60)
axes.invert_yaxis()

def animate(i):
    global p_ps
    p_T.set_array(t500[i][::-1].values.flatten())
    # https://scipython.com/blog/animated-contour-plots-with-matplotlib/
    for coll in p_ps.collections:
        coll.remove()
    p_ps = ps[i][::-1].plot.contour(add_colorbar=False,levels=np.arange(-65,-25,7))
    axes.set_title(f"Day {131+i/24:.2f}")
    axes.invert_yaxis()
    return p_ps


ani = animation.FuncAnimation(fig, animate, frames=len(ps),
                              interval=25) # , blit=True, init_func=init)

plt.show()
