import numpy             as np
import matplotlib        as mpl
import matplotlib.pyplot as plt

# use the LaTeX fonts for everything :
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['text.usetex'] = True

data_dir = "../data/"

# create a figure with twin x-coord axes :
fig = plt.figure(figsize=(6,3.5))
gs  = mpl.gridspec.GridSpec(nrows=2, ncols=1,
                            hspace=0.1, height_ratios=[1,1],
                            left=0.1, right=0.95, top=0.95, bottom=0.15)
ax1 = plt.subplot(gs[0,0])
ax2 = plt.subplot(gs[1,0])

simulations = ['a', 'b', 'c']
symbols     = ['.', 'x', '+']

for sim,sym in zip(simulations,symbols):

	n_a   = np.loadtxt(data_dir + sim + "/n.txt").astype('int')
	itr_a = np.loadtxt(data_dir + sim + "/i.txt").astype('int')
	tme_a = np.loadtxt(data_dir + sim + "/t.txt")
	err_a = np.loadtxt(data_dir + sim + "/e.txt")

	ax1.semilogy(err_a, 'k' + sym)

	if sim == 'a':
		h = 1 / n_a
		y = h**2
		ax1.plot(y, c='r', label="$h = 1/n^2$")

	ax2.semilogy(itr_a, 'k' + sym, label="simulation (" + sim + ")")

ax2.set_xticklabels(np.append(0,n_a))

ax2.set_xlabel(r"$n$")
ax2.set_ylabel(r"iteration count")
ax2.legend(loc="lower right", prop={'size':8})
ax2.grid()

# format the plot to look nice and save it as a pdf :
ax1.set_ylabel(r"$\Vert \underline{u} - \underline{u}_e \Vert_2 / n$")
ax1.set_xlabel(r"$x$")
ax1.xaxis.set_ticklabels([])
ax1.legend(loc="center left", prop={'size':8})
ax1.grid()

plt.savefig("../images/prob_2_convergence.pdf")
plt.close()



