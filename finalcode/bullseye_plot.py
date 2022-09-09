"""
=======================
Left ventricle bullseye
=======================

This example demonstrates how to create the 17 segment model for the left
ventricle recommended by the American Heart Association (AHA).
"""
#import xlrd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

plt.rcParams["font.family"] = "times new roman"
plt.rcParams.update({'font.size': 12})

def bullseye_plot(ax, data, seg_bold=None, cmap=None, norm=None):
    """
    Bullseye representation for the left ventricle.

    Parameters
    ----------
    ax : axes
    data : list of int and float
        The intensity values for each of the 17 segments
    seg_bold : list of int, optional
        A list with the segments to highlight
    cmap : ColorMap or None, optional
        Optional argument to set the desired colormap
    norm : Normalize or None, optional
        Optional argument to normalize data into the [0.0, 1.0] range


    Notes
    -----
    This function create the 17 segment model for the left ventricle according
    to the American Heart Association (AHA) [1]_

    References
    ----------
    .. [1] M. D. Cerqueira, N. J. Weissman, V. Dilsizian, A. K. Jacobs,
        S. Kaul, W. K. Laskey, D. J. Pennell, J. A. Rumberger, T. Ryan,
        and M. S. Verani, "Standardized myocardial segmentation and
        nomenclature for tomographic imaging of the heart",
        Circulation, vol. 105, no. 4, pp. 539-542, 2002.
    """
    if seg_bold is None:
        seg_bold = []

    linewidth = 2
    data = np.array(data).ravel()

    if cmap is None:
        cmap = plt.cm.viridis

    if norm is None:
        norm = mpl.colors.Normalize(vmin=data.min(), vmax=data.max())

    theta = np.linspace(0, 2 * np.pi, 768)
    r = np.linspace(0.2, 1.0, 4)

    # Create the bound for the segment 17
    for i in range(r.shape[0]):
        ax.plot(theta, np.repeat(r[i], theta.shape), '-k', lw=linewidth)

    # Create the bounds for the segments 1-12
    for i in range(6):
        theta_i = np.deg2rad(i * 60)
        ax.plot([theta_i, theta_i], [r[1], 1], '-k', lw=linewidth)

    # Create the bounds for the segments 13-16
    for i in range(4):
        theta_i = np.deg2rad(i * 90 - 45)
        ax.plot([theta_i, theta_i], [r[0], r[1]], '-k', lw=linewidth)

    # Fill the segments 1-6
    r0 = r[2:4]
    r0 = np.repeat(r0[:, np.newaxis], 128, axis=1).T
    for i in range(6):
        # First segment start at 60 degrees
        theta0 = theta[i * 128:i * 128 + 128] + np.deg2rad(60)
        theta0 = np.repeat(theta0[:, np.newaxis], 2, axis=1)
        z = np.ones((128, 2)) * data[i]
        ax.pcolormesh(theta0, r0, z, cmap=cmap, norm=norm)
        if i + 1 in seg_bold:
            ax.plot(theta0, r0, '-k', lw=linewidth + 2)
            ax.plot(theta0[0], [r[2], r[3]], '-k', lw=linewidth + 1)
            ax.plot(theta0[-1], [r[2], r[3]], '-k', lw=linewidth + 1)

    # Fill the segments 7-12
    r0 = r[1:3]
    r0 = np.repeat(r0[:, np.newaxis], 128, axis=1).T
    for i in range(6):
        # First segment start at 60 degrees
        theta0 = theta[i * 128:i * 128 + 128] + np.deg2rad(60)
        theta0 = np.repeat(theta0[:, np.newaxis], 2, axis=1)
        z = np.ones((128, 2)) * data[i + 6]
        ax.pcolormesh(theta0, r0, z, cmap=cmap, norm=norm)
        if i + 7 in seg_bold:
            ax.plot(theta0, r0, '-k', lw=linewidth + 2)
            ax.plot(theta0[0], [r[1], r[2]], '-k', lw=linewidth + 1)
            ax.plot(theta0[-1], [r[1], r[2]], '-k', lw=linewidth + 1)

    # Fill the segments 13-16
    r0 = r[0:2]
    r0 = np.repeat(r0[:, np.newaxis], 192, axis=1).T
    for i in range(4):
        # First segment start at 45 degrees
        theta0 = theta[i * 192:i * 192 + 192] + np.deg2rad(45)
        theta0 = np.repeat(theta0[:, np.newaxis], 2, axis=1)
        z = np.ones((192, 2)) * data[i + 12]
        ax.pcolormesh(theta0, r0, z, cmap=cmap, norm=norm)
        if i + 13 in seg_bold:
            ax.plot(theta0, r0, '-k', lw=linewidth + 2)
            ax.plot(theta0[0], [r[0], r[1]], '-k', lw=linewidth + 1)
            ax.plot(theta0[-1], [r[0], r[1]], '-k', lw=linewidth + 1)

    # Fill the segments 17
    if data.size == 17:
        r0 = np.array([0, r[0]])
        r0 = np.repeat(r0[:, np.newaxis], theta.size, axis=1).T
        theta0 = np.repeat(theta[:, np.newaxis], 2, axis=1)
        z = np.ones((theta.size, 2)) * data[16]
        ax.pcolormesh(theta0, r0, z, cmap=cmap, norm=norm)
        if 17 in seg_bold:
            ax.plot(theta0, r0, '-k', lw=linewidth + 2)

    ax.set_ylim([0, 1])
    ax.set_yticklabels([])
    ax.set_xticklabels([])

def tick_mat(data):

    MIN = round(np.amin(data),4)
    MAX = round(np.amax(data),4)
    h = (MAX-MIN)/3
    m1 = round(MIN+h, 4)
    m2 = round(MAX-h, 4)
    #m0 = round(MIN-(h/2),3)
    #m3 = round(MAX+(h/2),3)
    m =np.array([MIN, m1, m2, MAX])

    return m


def tick_mat2(data1, data2):

    MIN = round(np.amin(data1),3)
    MAX = round(np.amax(data2),3)
    h = (MAX-MIN)/3
    m1 = round(MIN+h, 4)
    m2 = round(MAX-h, 4)

    m =np.array([MIN, m1, m2, MAX])

    return m




def LVbullseyemulti(data, filename, Stimulus, DirectorY, maxm, minm):

	plt.rcParams["font.family"] = "times new roman"
	plt.rcParams.update({'font.size': 14})
	plt.rcParams.update({'font.weight': 'bold'})
	min1 = minm #np.min(data)
	max1 = maxm #np.max(data)


	m2 = tick_mat2(min1, max1)

	# Make a figure and axes with dimensions as desired.
	fig, ax = plt.subplots(figsize=(12, 8), nrows=1, ncols=3,
                       subplot_kw=dict(projection='polar'))
	fig.canvas.set_window_title(Stimulus)

	# Create the axis for the colorbars

	axl4 = fig.add_axes([0.12, 0.25, 0.78, 0.06])


	# Set the colormap and norm to correspond to the data for which
	# the colorbar will be used.

	cmap = mpl.cm.viridis # mpl.cm.plasma#mpl.cm.cool#mpl.cm.inferno# #
	norm2 = mpl.colors.Normalize(vmin=round(np.amin(min1),3), vmax=round(np.amax(max1),3))

	cb3 = mpl.colorbar.ColorbarBase(axl4, cmap=cmap, norm=norm2, ticks = m2, orientation='horizontal')
	cb3.ax.tick_params(labelsize=14)
		

	bullseye_plot(ax[0], data[0][:], cmap=cmap, norm=norm2)
	ax[0].set_title('Control', fontweight= 'bold')

	bullseye_plot(ax[1], data[1][:], cmap=cmap, norm=norm2)
	ax[1].set_title('Non-obstructive', fontweight ='bold')


	bullseye_plot(ax[2], data[2][:], cmap=cmap, norm=norm2)
	ax[2].set_title('Obstructive', fontweight = 'bold')


	plt.savefig(DirectorY + filename +'_1.png')

	#plt.show()

	plt.close()



def LVbullseyemulti6(data, filename, Stimulus, DirectorY, maxm, minm):


	plt.rcParams["font.family"] = "times new roman"
	plt.rcParams.update({'font.size': 14})
	plt.rcParams.update({'font.weight': 'bold'})
	#data = data/abs(np.median(data))
	min1 = minm #np.min(data)
	max1 = maxm #np.max(data)

	m2 = tick_mat2(min1, max1)

	# Make a figure and axes with dimensions as desired.
	fig, ax = plt.subplots(figsize=(12, 8), nrows=2, ncols=3,
                       subplot_kw=dict(projection='polar'))
	fig.canvas.set_window_title(Stimulus)

	# Create the axis for the colorbars


	axl4 = fig.add_axes([0.12, 0.05, 0.78, 0.05])

	# Set the colormap and norm to correspond to the data for which
	# the colorbar will be used.

	cmap = mpl.cm.viridis # mpl.cm.plasma#mpl.cm.cool#mpl.cm.inferno# #
	norm2 = mpl.colors.Normalize(vmin=round(np.amin(min1),3), vmax=round(np.amax(max1),3))

	cb3 = mpl.colorbar.ColorbarBase(axl4, cmap=cmap, norm=norm2, ticks = m2, orientation='horizontal')
	cb3.ax.tick_params(labelsize=20)
		

	bullseye_plot(ax[0][0], data[0][:], cmap=cmap, norm=norm2)
	ax[0][0].set_title(r'$\kappa $ = 0', fontweight = 'bold')
	#ax[0][0].set_title('m = 0', fontweight = 'bold')

	bullseye_plot(ax[0][1], data[1][:], cmap=cmap, norm=norm2)
	ax[0][1].set_title(r'$\kappa $ = 0.07', fontweight = 'bold')
	#ax[0][1].set_title('m = 1', fontweight = 'bold')


	bullseye_plot(ax[0][2], data[2][:], cmap=cmap, norm=norm2)
	ax[0][2].set_title(r'$\kappa $= 0.1', fontweight = 'bold')
	#ax[0][2].set_title('m = 2', fontweight = 'bold')

	bullseye_plot(ax[1][0], data[3][:], cmap=cmap, norm=norm2)
	ax[1][0].set_title(r'$\kappa $ = 0.14', fontweight = 'bold')


	bullseye_plot(ax[1][1], data[4][:], cmap=cmap, norm=norm2)
	ax[1][1].set_title(r'$\kappa $= 0.18', fontweight = 'bold')


	bullseye_plot(ax[1][2], data[5][:], cmap=cmap, norm=norm2)
	ax[1][2].set_title(r'$\kappa $ = 0.22', fontweight = 'bold')

	plt.savefig(DirectorY+ filename +'_1.png')

	#plt.show()

	plt.close()


def LVbullseyemulti5(data, filename, Stimulus, DirectorY, maxm, minm):


	plt.rcParams["font.family"] = "times new roman"
	plt.rcParams.update({'font.size': 14})
	plt.rcParams.update({'font.weight': 'bold'})
	#data = data/abs(np.median(data))
	min1 = minm #np.min(data)
	max1 = maxm #np.max(data)

	m2 = tick_mat2(min1, max1)

	# Make a figure and axes with dimensions as desired.
	fig, ax = plt.subplots(figsize=(12, 8), nrows=2, ncols=3,
                       subplot_kw=dict(projection='polar'))
	fig.canvas.set_window_title(Stimulus)

	# Create the axis for the colorbars


	axl4 = fig.add_axes([0.12, 0.05, 0.78, 0.05])

	# Set the colormap and norm to correspond to the data for which
	# the colorbar will be used.

	cmap = mpl.cm.viridis # mpl.cm.plasma#mpl.cm.cool#mpl.cm.inferno# #
	norm2 = mpl.colors.Normalize(vmin=round(np.amin(min1),3), vmax=round(np.amax(max1),3))

	cb3 = mpl.colorbar.ColorbarBase(axl4, cmap=cmap, norm=norm2, ticks = m2, orientation='horizontal')
	cb3.ax.tick_params(labelsize=20)
		

	bullseye_plot(ax[0][0], data[0][:], cmap=cmap, norm=norm2)
	ax[0][0].set_title(r'$\kappa $ = 0', fontweight = 'bold')
	#ax[0][0].set_title('m = 0', fontweight = 'bold')

	bullseye_plot(ax[0][1], data[1][:], cmap=cmap, norm=norm2)
	ax[0][1].set_title(r'$\kappa $ = 0.07', fontweight = 'bold')
	#ax[0][1].set_title('m = 1', fontweight = 'bold')


	bullseye_plot(ax[0][2], data[2][:], cmap=cmap, norm=norm2)
	ax[0][2].set_title(r'$\kappa $= 0.1', fontweight = 'bold')
	#ax[0][2].set_title('m = 2', fontweight = 'bold')

	bullseye_plot(ax[1][0], data[3][:], cmap=cmap, norm=norm2)
	ax[1][0].set_title(r'$\kappa $ = 0.14', fontweight = 'bold')


	bullseye_plot(ax[1][1], data[4][:], cmap=cmap, norm=norm2)
	ax[1][1].set_title(r'$\kappa $= 0.18', fontweight = 'bold')



	plt.savefig(DirectorY+ filename +'_1.png')

	#plt.show()

	plt.close()
