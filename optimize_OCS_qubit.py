#!/usr/bin/env python
import matplotlib
import pandas as pd
from qiskit_metal.qlibrary.qubits.transmon_cross import TransmonCross
import qiskit_metal as metal
from qiskit_metal import designs, draw
from qiskit_metal import MetalGUI, Dict, Headings
from qiskit_metal.analyses.quantization import EPRanalysis
import pyEPR as epr
from pyEPR.calcs import Convert
from datetime import datetime
from LC_JJ import *
from tqdm import tqdm


def launch_Metal_GUI():
    design = designs.DesignPlanar()

    design._chips['main']['size']['size_x'] = '5mm'
    design._chips['main']['size']['size_y'] = '5mm'

    design.variables['cpw_width'] = '10 um'
    design.variables['cpw_gap'] = '6 um'

    gui = MetalGUI(design)
    return gui, design



def create_OCS_qubit(cross_length, cross_width, cross_gap):
    ########## OCS QUBUIT ########
    options_d = dict(
        cross_width = '{}um'.format(cross_width),
        cross_length = '{}um'.format(cross_length),
        cross_gap = '{}um'.format(cross_gap),
        chip='main',
        connection_pads=dict(
            ground_pin  = dict(connector_location = '0', connector_type = '0', claw_width = '19.1um', claw_length = '35.9um',
                            ground_spacing = '5um')
        )
    )

    q_d = TransmonCross(design, 'Q1', options = dict(
        pos_x='-1250um', pos_y='-1mm', orientation = '90', **options_d))

    gui.rebuild()
    return q_d



def get_all_parameters_of_interest(hfss, eig_qb, Lj, Cj, pass_num=10):
    hfss.clean_active_design()
    hfss.render_design(['Q1'], [])


    # Analysis properties
    setup = hfss.pinfo.setup
    setup.passes = pass_num
    print(f"""
    Number of eigenmodes to find             = {setup.n_modes}
    Number of simulation passes              = {setup.passes}
    Convergence freq max delta percent diff  = {setup.delta_f}
    """)

    pinfo = hfss.pinfo
    pinfo.design.set_variable('Lj', '{} nH'.format(Lj))
    pinfo.design.set_variable('Cj', '{} fF'.format(Cj))

    setup.analyze()

    df1 = eig_qb.get_frequencies()

    qubit_freq  = df1["Freq. (GHz)"][0]
    print(f"qubit frequency = {qubit_freq} GHz")


    pinfo = hfss.pinfo
    pinfo.junctions['jj'] = {'Lj_variable': 'Lj', 'rect': 'JJ_rect_Lj_Q1_rect_jj', 
                                'line': 'JJ_Lj_Q1_rect_jj_',  'Cj_variable': 'Cj'}
    pinfo.validate_junction_info() # Check that valid names of variables and objects have been supplied
    pinfo.dissipative['dielectrics_bulk'] = ['main'] # Dissipative elements: specify

    eprd = epr.DistributedAnalysis(pinfo)
    eprd.do_EPR_analysis()

    epra = epr.QuantumAnalysis(eprd.data_filename)
    sim_info = epra.analyze_all_variations(cos_trunc = 8, fock_trunc = 7)

    df = pd.DataFrame(sim_info, columns=sim_info.keys())

    alpha = -df["0"]["chi_O1"].values[0][0]
    print(f"Anharmonicity is {alpha} MHz")

    w03 = 3*(qubit_freq) + 2*alpha*1e-3
    print(f"w_03 is {w03} GHz")

    E_c = -alpha
    print(f"E_c is {E_c} MHz")


    E_j = Convert.Ej_from_Lj(10, 'nH', "GHz")
    print(f"E_j is {E_j} GHz")

    ratio = E_j / E_c*1e3
    print(f"E_j / E_c = {ratio}")

    return qubit_freq, ratio, alpha


def keep_record(data, fname):
    df = pd.DataFrame([data])
    df.to_csv(fname)


def run_simulation(hfss, eig_qb, Lj, Cj, cross_length, cross_width, cross_gap, pass_num):
    qubit_freq, ratio, alpha = get_all_parameters_of_interest(hfss, eig_qb, Lj, Cj, pass_num)
    w03 = 3*(qubit_freq) + 2*alpha*1e-3
    data = {"cross_length (um)": cross_length, "cross_width (um)": cross_width, "cross_gap (um)": cross_gap,
            "Lj (nH)": Lj, "Cj (fF)": Cj, "qubit_freq (GHz)": qubit_freq, "ratio": ratio, "alpha (MHz)": alpha, "w03 (GHz)": w03}

    fname = "data/simulation_info_{}.csv".format(datetime.today().strftime('%Y-%m-%d-%H-%M-%S'))

    keep_record(data, fname)
    print(f"\n\n{fname} has been created.\n\n")
    print(30*"=")


def sweep_Lj(Lj_range, hfss, eig_qb, cross_length, cross_width, cross_gap, pass_num=10, target_qubit_frequency=3, sim_offset=1.47):
    target_qubit_frequency = target_qubit_frequency - sim_offset #GHz
    for Lj in Lj_range:
        Cj = round(get_Cj_from_Lj(Lj, target_qubit_frequency),2) #fF
        run_simulation(hfss, eig_qb, Lj, Cj, cross_length, cross_width, cross_gap, pass_num)

def sweep_cross_length(p_range, hfss, eig_qb, Lj, cross_width, cross_gap, pass_num=10,target_qubit_frequency=3, sim_offset=1.47):
    target_qubit_frequency = target_qubit_frequency - sim_offset #GHz
    for cross_length in p_range:
        Cj = round(get_Cj_from_Lj(Lj, target_qubit_frequency),2) #fF
        run_simulation(hfss, eig_qb, Lj, Cj, cross_length, cross_width, cross_gap, pass_num)

def sweep_cross_width(p_range, hfss, eig_qb, Lj, cross_length, cross_gap, pass_num,target_qubit_frequency=3, sim_offset=1.47):
    target_qubit_frequency = target_qubit_frequency - sim_offset #GHz
    for cross_length in p_range:
        Cj = round(get_Cj_from_Lj(Lj, target_qubit_frequency),2) #fF
        run_simulation(hfss, eig_qb, Lj, Cj, cross_length, cross_width, cross_gap, pass_num)

def sweep_cross_gap(p_range, hfss, eig_qb, Lj, cross_length, cross_width, pass_num=10,target_qubit_frequency=3, sim_offset=1.47):
    target_qubit_frequency = target_qubit_frequency - sim_offset #GHz
    for cross_length in p_range:
        Cj = round(get_Cj_from_Lj(Lj, target_qubit_frequency),2) #fF
        run_simulation(hfss, eig_qb, Lj, Cj, cross_length, cross_width, cross_gap, pass_num)




if __name__ == "__main__":
    matplotlib.use("Agg")
    gui, design = launch_Metal_GUI()
    #Allow running the same cell here multiple times to overwrite changes
    design.overwrite_enabled = True

    cross_length, cross_width, cross_gap = 150, 25, 40 #um

    q_ocs = create_OCS_qubit(cross_length, cross_width, cross_gap)

    pass_num = 15
    eig_qb = EPRanalysis(design, "hfss")
    hfss = eig_qb.sim.renderer
    hfss.start()
    hfss.activate_ansys_design("ocs_optimize", 'eigenmode')  # use new_ansys_design() to force creation of a blank design

    """
    # L_j sweep
    Lj_range = np.arange(10,101,10)
    sweep_Lj(Lj_range, hfss, eig_qb, cross_length, cross_width, cross_gap, pass_num,target_qubit_frequency=3, sim_offset=1.47)
    """

    # Sim Once
    target_qubit_frequency = 3 #GHz
    sim_offset = 1.47 #GHz
    target_qubit_frequency = target_qubit_frequency - sim_offset #GHz
    Lj = 40 #nH
    Cj = round(get_Cj_from_Lj(Lj, target_qubit_frequency),2) #fF
    run_simulation(hfss, eig_qb, Lj, Cj, cross_length, cross_width, cross_gap, pass_num)

    """
    # 3 D sweep

    target_qubit_frequency = 3 #GHz
    sim_offset = 1.47 #GHz
    target_qubit_frequency = target_qubit_frequency - sim_offset #GHz
    Lj = 20 #(nH) from Lj Sweep
    Cj = round(get_Cj_from_Lj(Lj, target_qubit_frequency),2) #fF

    cross_lengths = np.arange(100,221,25)
    cross_widths = np.arange(10,51,10)
    cross_gaps = np.arange(10,51,10)

    for cross_length in tqdm(cross_lengths):
        for cross_gap in cross_gaps:
            for cross_width in cross_widths:
                run_simulation(hfss, eig_qb, Lj, Cj, cross_length, cross_width, cross_gap, pass_num)
    """


