"""
This file experiments with the postshock detonation delay for the Propane-Air CJ-detonation limits for different stoichiometric conditions`
"""
import cantera as ct
import numpy as np
from sdtoolbox.postshock import CJspeed
from sdtoolbox.postshock import PostShock_eq, PostShock_fr
from sdtoolbox.thermo import soundspeed_eq, soundspeed_fr

if __name__ == "__main__":
    # describe composition of mixture
    inittemp, initpress = (300,1)
    # calculate composition based on percentage of propane
    HC_comp = 0.092
    O_comp = (1 - HC_comp) * 0.21
    N2_comp = (1 - HC_comp) * 0.79
    q = {'C3H8': HC_comp, 'O2': O_comp, 'N2': N2_comp}
    mech = 'gri30.cti'
    
    
    print('\nCJ computation for ' + mech + ' with composition:')
    print(q)
    print('Initial conditions: P1 = %.3e atm & T1 = %.2f K' % (initpress, inittemp))
    # generate gas solution and set initial conditions
    preshockgas = ct.Solution(mech)
    preshockgas.TPX = inittemp, initpress * ct.one_atm, q

    # calculate CJ speed for gas mixture
    cj_speed = CJspeed(initpress*ct.one_atm, inittemp, q, mech, fullOutput=False)

    # calculate postshock gas condition
    postshockgas = PostShock_fr(cj_speed, initpress*ct.one_atm, inittemp, q, mech)
    # save and print initial shock temp and pressure
    shockedpress = postshockgas.P/101325
    shockedtemp = postshockgas.T
    print(f'    PostShock Pressure: {postshockgas.P/101325} atm')
    print(f'    PostShock Temperature: {postshockgas.T} K')

    # Resolution: The PFR will be simulated by 'n_steps' time steps
    n_steps = 300000
    
    # create a new reactor
    r1 = ct.IdealGasConstPressureReactor(postshockgas)
    # create a reactor network for performing time integration
    sim1 = ct.ReactorNet([r1])
    
    # approximate a time step to achieve a similar resolution as in the next method
    
    dt = 0.25 / n_steps
    det_time = 0
    # define timesteps
    timesteps = (np.arange(n_steps) + 1) * dt
    for t_i in timesteps:
        # perform time integration
        sim1.advance(t_i)
        # store current time
        det_time = t_i
        # compute velocity and transform into space
        # print(r1.thermo.state.T[0])
        if(r1.thermo.state.T[0] > shockedtemp + 150):
            break

    print(f'    PostShock Autoignition Delay: {det_time} s')
