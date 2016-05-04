"""Methods for plotting and assessing convergence of a CACTUS run."""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm


def new_figax(ax=None):
    """Create a figure and axis handle.

    Parameters
    ----------
    ax : None, optional
        Description
    """
    if not ax:
        fig = plt.figure()
        ax = fig.add_subplot(111)

    return ax


def plot_rotor_vs_time(run, qty_name, ax=None, rev=False, **plot_kwargs):
    """Plot rotor data (time and rev-averaged) against simulation time.

    Parameters
    ----------
    run : CactusRun
        CactusRun instance.
    qty_name : str
        Name of column from rev and time data to plot.
    ax : axis handle
        Axis handle on which to plot.
    **plot_kwargs : dict
        Keyword arguments for plot command.
    """
    ax = new_figax(ax)

    if rev:
        xtime = run.time_to_rev(run.time_data['Normalized Time (-)'])
        xrev  = run.rev_data['Rev']
        xlabel = 'Rev'
    else:
        xtime = run.time_data['Normalized Time (-)']
        xrev  = run.rev_to_time(run.rev_data['Rev'])
        xlabel = 'Normalized Time (-)'

    # plot time series
    line, = ax.plot(xtime, run.time_data[qty_name], alpha=0.5, **plot_kwargs)

    # plot rev-averaged series
    ax.plot(xrev,
            run.rev_data[qty_name].values,
            color=line.get_color(),
            linewidth=2, marker='o')

    ax.set_xlabel(xlabel)
    ax.set_ylabel(qty_name)
    ax.grid(True)
    return ax


def plot_blade_1rev_avg(run, qty_name, ax=None, blade_num=0, **plot_kwargs):
    """Plot blade data averaged over the final revolution (last nti steps).

    Parameters
    ----------
    run : CactusRun
        CactusRun instance.
    qty_name : str
        Name of column from rev and time data to plot.
    ax : axis handle
        Axis handle on which to plot.
    blade_num : int, optional
        Description
    **plot_kwargs : dict
        Keyword arguments for plot command.
    """
    ax = new_figax(ax)
    nti = run.input.namelist['configinputs']['nti']
    avg_data = run.bladeelem_data.data_time_average(qty_name,
                                                    range(-nti, 0), blade_num)

    ax.plot(run.geom.blades[blade_num]['r_over_R_elem'],
            avg_data, **plot_kwargs)

    ax.set_xlabel('r/R')
    ax.set_ylabel(qty_name)
    ax.grid(True)
    return ax


def plot_blade_inst(run, qty_name, time_index,
                    ax=None,
                    blade_num=0,
                    **plot_kwargs):
    """Plot instantaneous blade data at a specified timestep.

    Parameters
    ----------
    run : CactusRun
        CactusRun instance.
    qty_name : str
        Name of column from rev and time data to plot.
    time_index : int
        Index of timestep at which to plot instantaneous blade data.
    ax : axis handle
        Axis handle on which to plot.
    blade_num : int, optional
        Blade number to plot (default is zero).
    **plot_kwargs : dict
        Keyword arguments for plot command.
    """
    ax = new_figax(ax)
    time, dfs = run.bladeelem_data.data_at_time_index(time_index)

    ax.plot(run.geom.blades[blade_num]['r_over_R_elem'],
            dfs[blade_num][qty_name],
            **plot_kwargs)

    ax.set_xlabel('r/R')
    ax.set_ylabel(qty_name)
    ax.grid(True)
    return ax


def plot_radial_kymograph(run, qty_name,
                          timesteps=None,
                          type='absolute',
                          blade_num=0,
                          contour_levels=None,
                          **pcolor_kwargs):
    """Plot blade quantity of a CactusRun instance on a radial kymograph.

    This is likely only useful for HAWT turbines.

    Parameters
    ----------
    run : CactusRun
        CactusRun instance.
    qty_name : str
        Name of column from rev and time data to plot.
    timesteps : list, optional
        List of timesteps (int) from which to generate the radial kymograph. If
        no timesteps are specified, defaults to all the timesteps in the final
        revolution.
    type : str, optional
        Type of kymograph -- `absolute` or `delta`. If absolute, the kymograph
        will be colored by the absolute value of the blade quantity in each
        location. If `delta`, the azimuthal average of each blade element is
        computed, and the delta between the azimuthal average and the
        instantaneous value will be shown. Default is `absolute`.
    blade_num : int, optional
        Blade number to plot (default is zero).
    contour_levels : None, optional
        Description
    **pcolor_kwargs : dict
        Keyword arguments for the `pcolor` command used to create plot.
    """
    if not timesteps:
        timesteps = range(-1, -run.nti - 2, -1)

    # get the number of elements
    nelem = run.geom.blades[0]['NElem'][0]

    # allocate storage
    Rs     = np.zeros([nelem, len(timesteps)])
    Thetas = np.zeros([nelem, len(timesteps)])
    Values = np.zeros([nelem, len(timesteps)])

    # draw a polar contour plot colored by the variable value
    fig, ax = plt.subplots(subplot_kw=dict(projection='polar'))
    ax.set_theta_zero_location('N')
    ax.set_xticks(np.pi / 180. * np.array([0, 90, 180, 270]))
    ax.set_rlabel_position(0)

    # if type is delta, we are plotting the instantaneous-avg, so compute avg
    # over the selected timesteps
    if type == 'delta':
        avg = run.bladeelem_data.data_time_average(qty_name,
                                                   timesteps,
                                                   blade_num)

    for theta_index, t_index in enumerate(timesteps):
        time, bladedata = run.bladeelem_data.data_at_time_index(t_index)
        if type == 'absolute':
            Values[:, theta_index] = bladedata[blade_num][qty_name]
        elif type == 'delta':
            Values[:, theta_index] = bladedata[blade_num][qty_name] - avg

        # store the (R,theta) locations of the elements
        Thetas[:, theta_index] = bladedata[blade_num]['Theta (rad)']
        Rs[:, theta_index]     = run.geom.blades[0]['r_over_R_elem']

    if contour_levels is None:
        im = ax.pcolor(Thetas, Rs, Values, cmap=cm.coolwarm, **pcolor_kwargs)

    vallim = [Values.min(), Values.max()]

    # display colorbar
    fig.colorbar(im)

    return fig, ax, im, vallim
