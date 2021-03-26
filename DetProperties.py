"""
Shock and Detonation Toolbox Demo Program

Calculates the CJ speed using the Minimum Wave Speed Method and
then finds the equilibrium state of the gas behind a shock wave
traveling at the CJ speed.

################################################################################
Theory, numerical methods and applications are described in the following report:

SDToolbox - Numerical Tools for Shock and Detonation Wave Modeling,
Explosion Dynamics Laboratory, Contributors: S. Browne, J. Ziegler,
N. Bitter, B. Schmidt, J. Lawson and J. E. Shepherd, GALCIT
Technical Report FM2018.001 Revised January 2021.

Please cite this report and the website if you use these routines.

Please refer to LICENCE.txt or the above report for copyright and disclaimers.

http://shepherd.caltech.edu/EDL/PublicResources/sdt/

################################################################################
Updated September 2018
Tested with:
    Python 3.5 and 3.6, Cantera 2.3 and 2.4
Under these operating systems:
    Windows 8.1, Windows 10, Linux (Debian 9)
"""
import cantera as ct
from sdtoolbox.postshock import CJspeed
from sdtoolbox.postshock import PostShock_eq
from sdtoolbox.thermo import soundspeed_eq, soundspeed_fr
import numpy as np
import matplotlib.pyplot as plt

if __name__ == "__main__":
    inittemp = np.array([300, 300, 600])
    initpress = np.array([100, 500, 100])
    compositions = np.linspace(2.2, 9.2, 30)
    # make arrays to hold detonations characteristics
    property_array = np.zeros_like(compositions)
    detionation_speed = np.stack(
        (property_array, property_array, property_array))
    detonation_speed = np.stack(
        (property_array, property_array, property_array))
    detonation_mach = np.stack(
        (property_array, property_array, property_array))
    product_temp = np.stack((property_array, property_array, property_array))
    product_press = np.stack((property_array, property_array, property_array))
    product_specificheat = np.stack(
        (property_array, property_array, property_array))
    calc_detspeed = np.stack((property_array, property_array, property_array))

    for IC, _ in enumerate(inittemp):
        for idx, propane_comp in enumerate(compositions):
            # Initial state specification:
            # P1 = Initial Pressure
            # T1 = Initial Temperature
            # U = Shock Speed
            # q = Initial Composition
            # mech = Cantera mechanism File name

            # calculate composition based on percentage of propane
            HC_comp = propane_comp
            O_comp = (100 - HC_comp) * 0.21
            N2_comp = (100 - HC_comp) * 0.79
            # initial pressures and temperatures
            P1 = initpress[IC] * 1000
            T1 = inittemp[IC]
            q = {'C3H8': HC_comp, 'O2': O_comp, 'N2': N2_comp}
            mech = 'gri30.cti'
            # make solution and initialize to starting conditions, store density
            gas_initial = ct.Solution(mech)
            gas_initial.TPX = T1, P1, q
            rho_1 = gas_initial.density
            csound = soundspeed_eq(gas_initial)

            # compute CJ speed, print, and add to speed array
            [cj_speed, R2, plot_data] = CJspeed(
                P1, T1, q, mech, fullOutput=True)
            # compute equilibrium CJ state parameters
            gas = PostShock_eq(cj_speed, P1, T1, q, mech)
            ae = soundspeed_eq(gas)
            af = soundspeed_fr(gas)
            rho_2 = gas.density
            gammae = ae**2 * rho_2 / gas.P
            gammaf = af**2 * rho_2 / gas.P
            w2 = cj_speed * rho_1 / rho_2
            u2 = cj_speed - w2
            D = cj_speed# u2 * rho_2 / rho_1
            Dcalc = (1.29 + 1)/1.29 * np.power(1.29*8314*2700/28, 0.5)
            mach = D / csound
            print('\nCJ computation for ' + mech + ' with composition:')
            print(q)
            print('Initial conditions: P1 = %.3e Pa & T1 = %.2f K' % (P1, T1))
            print('CJ Speed   %.1f m/s' % cj_speed)
            print('CJ State')
            print('   Pressure   %.3e Pa' % gas.P)
            print('   Temperature  %.1f K' % gas.T)
            print('   Density  %.3f kg/m3' % gas.density)
            print('   Entropy  %.3e J/K' % gas.entropy_mass)
            print('   w2 (wave frame) %.1f m/s' % w2)
            print('   u2 (lab frame) %.1f m/s' % u2)
            print('   c2 frozen %.1f m/s' % af)
            print('   c2 equilbrium %.1f m/s' % ae)
            print('   gamma2 frozen %.3f ' % gammaf)
            print('   gamma2 equilbrium %.3f ' % gammae)
            print('   Speed of sound %.3f ' % csound)
            print('   Det Speed %.3f ' % D)
            print('   Calc 2700K Det Speed %.3f ' % Dcalc)
            print('   Mach Number %.3f ' % mach)
            print(f'   Density 1:{rho_1}, Density 2:{rho_2}')

            product_temp[IC, idx] = gas.T
            detonation_speed[IC, idx] = D
            product_press[IC, idx] = gas.P
            detonation_mach[IC, idx] = mach
            product_specificheat[IC, idx] = gammaf

    plt.plot(compositions, detonation_speed[0, :], 'r-', label='T$_0$=300 K, P$_0$=100 kPa')
    plt.plot(compositions, detonation_speed[1, :], 'g-', label='T$_0$=300 K, P$_0$=500 kPa')
    plt.plot(compositions, detonation_speed[2, :], 'b-', label='T$_0$=600 K, P$_0$=100 kPa')
    plt.legend()
    plt.xlabel('C$_3$H$_8$ Initial Composition (%)')
    plt.ylabel('Detonation Speed (m/s)')
    plt.title('C3H8-Air Detonation Speed')
    figname = 'Figures/DetSpeed.png'
    plt.savefig(figname)
    plt.show()

    plt.plot(compositions, product_temp[0, :], 'r-', label='T$_0$=300 K, P$_0$=100 kPa')
    plt.plot(compositions, product_temp[1, :], 'g-', label='T$_0$=300 K, P$_0$=500 kPa')
    plt.plot(compositions, product_temp[2, :], 'b-', label='T$_0$=600 K, P$_0$=100 kPa')
    plt.legend()
    plt.xlabel('C$_3$H$_8$ Initial Composition (%)')
    plt.ylabel('Product Temperature (K)')
    plt.title('C3H8-Air Detonation Product T')
    figname = 'Figures/DetTemp.png'
    plt.savefig(figname)
    plt.show()

    plt.plot(compositions, product_press[0, :] / 101325, 'r-', label='T$_0$=300 K, P$_0$=100 kPa')
    plt.plot(compositions, product_press[1, :] / 101325, 'g-', label='T$_0$=300 K, P$_0$=500 kPa')
    plt.plot(compositions, product_press[2, :] / 101325, 'b-', label='T$_0$=600 K, P$_0$=100 kPa')
    plt.legend()
    plt.xlabel('C$_3$H$_8$ Initial Composition (%)')
    plt.ylabel('Product Pressure (atm)')
    plt.title('C3H8-Air Detonation Product P')
    figname = 'Figures/DetPress.png'
    plt.savefig(figname)
    plt.show()

    plt.plot(compositions, detonation_mach[0, :], 'r-', label='T$_0$=300 K, P$_0$=100 kPa')
    plt.plot(compositions, detonation_mach[1, :], 'g-', label='T$_0$=300 K, P$_0$=500 kPa')
    plt.plot(compositions, detonation_mach[2, :], 'b-', label='T$_0$=600 K, P$_0$=100 kPa')
    plt.legend()
    plt.xlabel('C$_3$H$_8$ Initial Composition (%)')
    plt.ylabel('Detonation Mach Number (unitless)')
    plt.title('C3H8-Air Detonation Mach Number')
    figname = 'Figures/DetMach.png'
    plt.savefig(figname)
    plt.show()

    plt.plot(compositions, product_specificheat[0, :], 'r-', label='T$_0$=300 K, P$_0$=100 kPa')
    plt.plot(compositions, product_specificheat[1, :], 'g-', label='T$_0$=300 K, P$_0$=500 kPa')
    plt.plot(compositions, product_specificheat[2, :], 'b-', label='T$_0$=600 K, P$_0$=100 kPa')
    plt.legend()
    plt.xlabel('C$_3$H$_8$ Initial Composition (%)')
    plt.ylabel('Specific Heat Ratio (unitless)')
    plt.title('C3H8-Air Specific Heat Ratio ($\gamma$)')
    figname = 'Figures/DetSpecHeat.png'
    plt.savefig(figname)
    plt.show()
