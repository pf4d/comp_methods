import numpy             as np
import matplotlib        as mpl
import matplotlib.pyplot as plt

# use the LaTeX fonts for everything :
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['text.usetex'] = True

data_dir = "../data/"

pi = np.pi

# exact solution :
def u_e(x):
	val = - 0.5*x**2 + x + 19 * np.cos(5 * pi * x) / (100 * pi) * (1 - x) \
	      + 19 / (500 * pi * pi) * np.sin(5 * pi * x) \
	      + (100 * pi - 19) / (100 * pi)
	return val

# create a figure with twin x-coord axes :
fig = plt.figure(figsize=(6,3.5))
gs  = mpl.gridspec.GridSpec(nrows=2, ncols=1,
                            hspace=0.1, height_ratios=[2,1],
                            left=0.1, right=0.95, top=0.95, bottom=0.15)
ax1 = plt.subplot(gs[0,0])
ax2 = plt.subplot(gs[1,0])

ls_a = [':', '-.', '--', '-']
n_a  = np.loadtxt(data_dir + "c/n.txt").astype('int')
n_a  = np.array([10, 40, 160, 640])

x_e = np.linspace(0, 1, 1000)
ax2.plot(x_e, u_e(x_e), 'k-')

# load the data :
for n, ls in zip(n_a, ls_a):

	x = np.loadtxt(data_dir + "c/x_%i.txt" % n)
	u = np.loadtxt(data_dir + "c/u_%i.txt" % n)

	# plot the data :
	ax1.semilogy(x, abs(u - u_e(x)), color='k', ls=ls, label=r"$n = %i$" % n)

# format the plot to look nice and save it as a pdf :
ax1.set_ylabel(r"$|u(x) - u_e(x)|$")
ax1.set_xlabel(r"$x$")
ax1.xaxis.set_ticklabels([])
ax1.legend(loc="lower right", prop={'size':8})
ax1.grid()

ax2.set_ylabel(r"$u_e(x)$")
ax2.set_xlabel(r"$x$")
ax2.grid()

plt.savefig("../images/prob_2_part_c.pdf")
plt.close()



