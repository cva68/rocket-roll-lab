"""
    ENEL321 Control Systems Rocket Lab

    Take a given transfer function for rocket roll angle and input disturbance
    with PID control, plot the simulated step response against experimental
    step response, and calculate the time domain specifications for both.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import TransferFunction, lsim

# Plot paramaters. Let CSV_FILE = None to produce only experimental values
CSV_FILE = "data/1.9 0 0.2.csv"
MIN_T = 0
MAX_T = 5.01
STEP_SIZE = 0.01

# Control Paramaters
K_P = 1.9
K_I = 0
K_D = 0.2

# System Paramaters
BETA = 80
ALPHA = 4
D = 400
INPUT_ANGLE = 90


def simulate_step_response(time_values):
    """
        Simulate the step response of the system
    """
    # Transfer functions
    denominator = [1, (ALPHA + BETA * K_D), BETA * K_P, BETA * K_I]
    numerator_r = [BETA * K_P, BETA * K_I]
    numerator_d = [D, 0]
    G_r = TransferFunction(numerator_r, denominator)
    G_d = TransferFunction(numerator_d, denominator)

    # Inputs
    U_r = INPUT_ANGLE * np.ones_like(time_values)  # 90 deg
    U_d = np.ones_like(time_values)  # No scaling

    # System response with lsim
    _t, y_r, _ = lsim(G_r, U=U_r, T=time_values)
    _t, y_d, _ = lsim(G_d, U=U_d, T=time_values)

    # Superposition
    y = y_r + y_d

    return y


def load_step_from_csv(csv_file):
    """
        Load a single step response from a given CSV file
    """
    csv_values = np.loadtxt(csv_file, delimiter=",", dtype=float, skiprows=1)

    # Find first rising edge
    first_0 = np.where(csv_values[:, 1] == 0)[0][0]
    first_90_after_0 = np.where(csv_values[:, 1][first_0:] == 90)[0][0]
    start = first_0 + first_90_after_0

    # Find first falling edge
    first_0_after_90 = np.where(csv_values[:, 1][start:] == 0)[0][0]
    end = start + first_0_after_90

    # Extract this one segment of data
    segment = csv_values[start:end]

    # Make time start at 0
    min_csv_time_value = segment[0, 0]
    segment[:, 0] -= min_csv_time_value

    return segment


def fit_csv_to_simulation(data, simulation_time):
    """
        Scale the time values of the experimental data to match
        that of the simulation data
    """
    # Find the time values associated with this section of data
    max_csv_time_value = data[-1, 0]

    # Scale the time to match the simulation time
    scale_value = simulation_time / max_csv_time_value
    time_values = data[:, 0] * scale_value

    # Return data from this time period
    angle_values = data[:, 2]
    u_in = data[:, 3]

    return time_values, angle_values, u_in


def plot_sim_vs_real(sim_t, sim_y, exp_t, exp_y):
    # Generate Plot
    plt.figure(figsize=(10, 5))
    plt.plot(sim_t, sim_y, label='Simulated output', linewidth=2)
    plt.plot(sim_t,
             90 * np.ones_like(sim_t),
             'k:',
             label='Reference (90°)',
             linewidth=2)

    if exp_t is not None:
        plt.plot(exp_t, exp_y, label='Experimental output', linewidth=2)
    plt.xlabel('Time (s)')
    plt.ylabel('Roll angle (deg)')
    plt.legend()
    plt.grid(True)
    plt.title(f'System Response with Kp={K_P}, Ki={K_I}, Kd={K_D}')

    # Save plot
    filename = f"plots/{K_P}-{K_I}-{K_D}-sim.png" if exp_t is None \
        else f"plots/{K_P}-{K_I}-{K_D}-exp.png"
    plt.savefig(filename)

    # Display the plot
    plt.show()


def find_time_domain_params(time_values, angle_values):
    """
        From the angle values of a CSV file, find time domain paramaters
    """
    rel_max = max(angle_values)
    print(rel_max)

    # Damping, exclusively from sim values
    damping = (ALPHA + BETA * K_D) / (2 * (BETA*K_P) ** 0.5)

    # Percent overshoot
    M_p_per = 100 * (rel_max / INPUT_ANGLE - 1)

    # Steady-state error, from last 10% of data
    e_ss = np.mean(angle_values[-int(0.1*len(angle_values)):]) - INPUT_ANGLE

    # Get closest recorded values to 10% and 90% of reference angle
    low_thresh = INPUT_ANGLE * 0.1
    high_thresh = INPUT_ANGLE * 0.9

    # First crossing of 10%
    idx_start = np.argmax(angle_values >= low_thresh)
    # First crossing of 90% *after* 10% has been crossed
    idx_end = idx_start + np.argmax(angle_values[idx_start:] >= high_thresh)

    t_r = time_values[idx_end] - time_values[idx_start]

    return damping, M_p_per, e_ss, t_r


def main():
    if CSV_FILE is not None:
        # Load CSV
        data = load_step_from_csv(CSV_FILE)
        exp_t, exp_y, _u = fit_csv_to_simulation(data, MAX_T)

        # Get time domain values from experimental data
        exp_damping, exp_M_p_per, exp_e_ss, exp_t_r = \
            find_time_domain_params(data[:, 0], data[:, 2])
    else:
        exp_M_p_per, exp_e_ss, exp_t_r = 0, 0, 0
        exp_t, exp_y = None, None

    # Simulate response
    sim_t = np.arange(MIN_T, MAX_T, STEP_SIZE)
    sim_y = simulate_step_response(sim_t)

    # Get time domain values from simulation data
    sim_damping, sim_M_p_per, sim_e_ss, sim_t_r = \
        find_time_domain_params(sim_t, sim_y)

    # Display nice table
    table_str = []

    table_str.append(f"Results for Kp={K_P}, Ki={K_I}, Kd={K_D}:\n")
    table_str.append("Val\t\tExp\t\tSim\n")
    table_str.append(f"ζ\t\t-\t\t{sim_damping:.1f}\n")
    table_str.append(f"Mp%\t\t{exp_M_p_per:.1f} %\t\t{sim_M_p_per:.1f} %\n")
    table_str.append(f"e_ss\t\t{exp_e_ss:.1f} deg\t\t{sim_e_ss:.1f} deg\n")
    table_str.append(f"t_r\t\t{exp_t_r:.2f} s\t\t{sim_t_r:.2f} s\n")

    print(''.join(table_str))
    with open(f"specs/{K_P}-{K_I}-{K_D}-sim.txt", 'w') as output_spec:
        output_spec.writelines(table_str)

    # Plot data
    plot_sim_vs_real(sim_t, sim_y, exp_t, exp_y)


if __name__ == "__main__":
    main()
