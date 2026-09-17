import numpy as np
import matplotlib.pyplot as plt
import h5py
import matplotlib.animation as animation
import sys

# Load result
with h5py.File(sys.argv[1] if len(sys.argv) > 1 else "output.h5", "r") as f:
    x = f["tallies/mesh_flux/grid/x"][:]
    x_mid = 0.5 * (x[:-1] + x[1:])
    y = f["tallies/mesh_flux/grid/y"][:]
    y_mid = 0.5 * (y[:-1] + y[1:])
    t = f["tallies/mesh_flux/grid/time"][:]
    t_mid = 0.5 * (t[:-1] + t[1:])
    X, Y = np.meshgrid(y, x)

    phi = f["tallies/mesh_flux/flux/mean"][:]
    phi_sd = f["tallies/mesh_flux/flux/sdev"][:]

    phi_total = f["tallies/density/density/mean"][:]
    phi_total_sd = f["tallies/density/density/sdev"][:]
    capture = f["tallies/detector/capture/mean"][:]
    capture_sd = f["tallies/detector/capture/sdev"][:]
    t_detector = f["tallies/detector/grid/time"][:]
    t_detector_mid = 0.5 * (t_detector[:-1] + t_detector[1:])

# Animate result
fig, ax = plt.subplots(
    1,
    3,
    figsize=(12, 4),
    gridspec_kw={"width_ratios": [1.0, 2, 1.0]},
    layout="constrained",
)
#
cax = ax[1].pcolormesh(X, Y, phi[0], vmin=phi[0].min(), vmax=phi[0].max())
ax[1].set_aspect("equal", "box")
ax[1].set_xlabel("$y$ [cm]")
ax[1].set_ylabel("$x$ [cm]")
#
ax[0].plot(t_mid, phi_total)
ax[0].set_xlabel("$t$ [s]")
ax[0].set_ylabel("Neutron density")
ax[0].set_yscale("log")
ax[0].plot(t_mid, phi_total, "b-")
ax[0].fill_between(
    t_mid, phi_total - phi_total_sd, phi_total + phi_total_sd, alpha=0.2, color="b"
)
ax[0].grid()
ax[0].set_box_aspect(1)
(line,) = ax[0].plot([], [], "ok", fillstyle="none")

# Cell tallies are integrated over each time bin; divide by its width for rates.
capture_rate = capture / np.diff(t_detector)
capture_rate_sd = capture_sd / np.diff(t_detector)
(detector_line,) = ax[2].plot(t_detector_mid, capture_rate)
ax[2].fill_between(
    t_detector_mid,
    capture_rate - capture_rate_sd,
    capture_rate + capture_rate_sd,
    alpha=0.2,
    color=detector_line.get_color(),
)
ax[2].set_xlabel("$t$ [s]")
ax[2].set_ylabel("Detector capture rate [s$^{-1}$]")
ax[2].set_yscale("log")
ax[2].grid()
ax[2].set_box_aspect(1)


#
def animate(i):
    n = np.zeros_like(t_mid)
    n[i] = phi_total[i]
    line.set_data(t_mid, n)
    cax.set_array(phi[i])
    cax.set_clim(phi[i].min(), phi[i].max())


#
K = len(t) - 1
anim = animation.FuncAnimation(fig, animate, frames=K)
plt.show()
