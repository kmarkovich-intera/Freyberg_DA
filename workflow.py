# script to run paired simple complex analysis of MF6 Freyberg model to compare batch and sequential DA

import os
import sys
import shutil
import platform
import string
import numpy as np
import pandas as pd
import platform
import pyemu
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
import matplotlib as mpl
import flopy
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.patches as mpatches

plt.rcParams.update({'font.size': 8})

port = 4021

# set path to pestpp executables
bin_path = os.path.join("bin")
if "linux" in platform.platform().lower():
    bin_path = os.path.join(bin_path, "linux")
elif "macos" in platform.platform().lower():
    bin_path = os.path.join(bin_path, "mac")
else:
    bin_path = os.path.join(bin_path, "win")

exe = ""
if "windows" in platform.platform().lower():
    exe = ".exe"
da_path = os.path.join(bin_path, "pestpp-da" + exe)
ies_path = os.path.join(bin_path, "pestpp-ies" + exe)

keep = ['arrobs_head_k:0_i:13_j:10', 'arrobs_head_k:2_i:2_j:9', 'arrobs_head_k:2_i:33_j:7', 'sfr_usecol:gage_1']
keep_labels = ["gw_1","gw_2","gw_3","sw_1"]
keep_units = ["$m$","$m$","$m$","$\\frac{m^3}{d}$"]
keep_dict = {k:l for k,l in zip(keep,keep_labels)}
keep_dict2 = {k:l for k,l in zip(keep,keep_units)}

# forecast = ["sfr_usecol:tailwater","sfr_usecol:headwater","arrobs_head_k:0_i:9_j:1"]
# forecast_labels = ["tailwater sw-gw exchg","headwater sw-gw exchg","gw forecast"]
# forecast_dict = {k:l for k,l in zip(forecast,forecast_labels)}
# forecast_units = ["$\\frac{m^3}{d}$","$\\frac{m^3}{d}$","$m$"]

# forecast = ["sfr_usecol:tailwater","sfr_usecol:headwater","arrobs_head_k:0_i:9_j:1",
#             "cum_mass_usecol:wel_out","arrobs_conc_k:2_i:26_j:11",]
# forecast_labels = ["tailwater sw-gw exchg","headwater sw-gw exchg","gw forecast","cumulative well mass removed",
#                   "gw conc"]
forecast = ["sfr_usecol:tailwater","sfr_usecol:headwater","arrobs_head_k:0_i:9_j:1"]
forecast_labels = ["tailwater sw-gw exchg","headwater sw-gw exchg","gw forecast"]
forecast_dict = {k:l for k,l in zip(forecast,forecast_labels)}
forecast_units = ["$\\frac{m^3}{d}$","$\\frac{m^3}{d}$","$m$"]
forecast_dict2 = {k:l for k,l in zip(forecast,forecast_units)}

keep_sngl_lyr = ['arrobs_head_k:0_i:13_j:10', 'arrobs_head_k:0_i:2_j:9', 'arrobs_head_k:0_i:33_j:7', 'sfr_usecol:gage_1']
sngl_lyr_dct = {k:l for k,l in zip(keep_sngl_lyr,keep)}
m_lrc = (1,25,5)

def clean_master_dirs():
    for i in range(100):
        ies_d = os.path.join('simple_master2_ies_{0}'.format(i))
        da_d = os.path.join('simple_master2_da_{0}'.format(i))

        os.chdir(ies_d)
        try:
            os.remove('prior.jcb')
        except:
            print('no jcb to remove')
        try:
            os.remove('ies_prior.jcb')
        except:
            print('no jcb to remove')

        os.chdir('..')
        os.chdir(da_d)

        try:
            os.remove('da_prior.jcb')
        except:
            print('no jcb to remove')

        os.chdir('..')


def compare_mf6_freyberg(num_workers=10,num_reals=100,num_replicates=100,use_sim_states=True,
                         run_ies=True,run_da=True,adj_init_states=True,seq_noptmax=1):
    complex_dir = os.path.join('daily_model_files_master_prior')
    bat_dir = os.path.join('monthly_model_files_template')
    seq_dir = "seq_" + bat_dir
    for ireal in range(num_replicates):

        ies_t_d = map_complex_to_simple_bat(complex_dir,bat_dir,ireal)

        # run batch and sequential simple models
        # ies stuff
        ies_pst = pyemu.Pst(os.path.join(ies_t_d, "freyberg.pst"))


        ies_pst.pestpp_options["ies_par_en"] = "prior.jcb"

        # set pestpp options for batch da
        ies_pst.pestpp_options.pop("ies_num_reals", None)
        ies_pst.pestpp_options.pop("da_num_reals", None)
        ies_pst.pestpp_options["ies_no_noise"] = False
        ies_pst.pestpp_options["ies_verbose_level"] = 1
        ies_pst.pestpp_options.pop("ies_localizer", None)
        ies_pst.pestpp_options["ies_autoadaloc"] = False
        ies_pst.pestpp_options["ies_save_lambda_en"] = False
        ies_pst.pestpp_options["ies_drop_conflicts"] = False
        ies_pst.pestpp_options["ies_num_reals"] = num_reals
        ies_pst.pestpp_options["ies_use_mda"] = False
        ies_pst.control_data.noptmax = 3

        #mark the future pars fixed
        par = ies_pst.parameter_data
        rch_par = par.loc[par.parnme.str.contains("direct_const_rch_recharge"),:].copy()
        rch_par.loc[:,"kper"] = rch_par.parnme.apply(lambda x: int(x.split('_')[4]) - 1)
        rch_par = rch_par.loc[rch_par.kper>=13,"parnme"]
        par.loc[rch_par,"partrans"] = "fixed"

        twel_par = par.loc[par.parnme.str.contains("twel_mlt"), :].copy()
        twel_par.loc[:, "kper"] = twel_par.parnme.apply(lambda x: int(x.split('_')[2]))
        twel_par = twel_par.loc[twel_par.kper>=13,"parnme"]
        par.loc[twel_par,"partrans"] = "fixed"


        ies_pst.write(os.path.join(ies_t_d, "freyberg.pst"), version=2)


        # run da          

        if run_da:
            da_t_d = map_simple_bat_to_seq(ies_t_d, seq_dir)

            # prep that prior ensemble for da
            da_pst = pyemu.Pst(os.path.join(da_t_d, "freyberg.pst"))

            # set pestpp options for sequential da
            da_pst.pestpp_options.pop("da_num_reals", None)
            da_pst.pestpp_options.pop("ies_num_reals", None)
            da_pst.pestpp_options["ies_no_noise"] = False
            da_pst.pestpp_options["ies_verbose_level"] = 1
            da_pst.pestpp_options.pop("ies_localizer", None)
            da_pst.pestpp_options["ies_autoadaloc"] = False
            da_pst.pestpp_options["ies_save_lambda_en"] = False
            da_pst.pestpp_options["ies_drop_conflicts"] = False
            da_pst.pestpp_options["ies_num_reals"] = num_reals
            da_pst.pestpp_options["ies_use_mda"] = False
            da_pst.pestpp_options["da_use_simulated_states"] = use_sim_states
            if not adj_init_states:
                par = da_pst.parameter_data
                istate_pars = par.loc[par.parnme.str.startswith("direct_head"),"parnme"]
                par.loc[istate_pars,"partrans"] = "fixed"
            da_pst.control_data.noptmax = seq_noptmax
            da_pst.write(os.path.join(da_t_d, "freyberg.pst"), version=2)

            m_da_dir = da_t_d.replace("template", "master")
            pyemu.os_utils.start_workers(da_t_d, 'pestpp-da', "freyberg.pst", port=port,
                                         num_workers=num_workers, master_dir=m_da_dir, verbose=True)

            shutil.rmtree(da_t_d)

        # run ies
        m_ies_dir = ies_t_d.replace("template","master")

        if run_ies:
            pyemu.os_utils.start_workers(ies_t_d, 'pestpp-ies', "freyberg.pst", port=port,
                                      num_workers=num_workers, master_dir=m_ies_dir, verbose=True)

        shutil.rmtree(ies_t_d)



def run_complex_prior_mc(c_d,num_workers=4):
    m_c_d = c_d.replace("template", "master_prior")
    pyemu.os_utils.start_workers(c_d, "pestpp-ies", "freyberg.pst", num_workers=num_workers, worker_root=".",
                                 master_dir=m_c_d)
    return m_c_d

def run_bat_prior_mc(b_d,num_workers=4):
    m_b_d = b_d.replace("template", "master_prior")
    pyemu.os_utils.start_workers(b_d, "pestpp-ies", "freyberg.pst", num_workers=num_workers, worker_root=".",
                                 master_dir=m_b_d)
    return m_b_d


def plot_phi_seq_bat():
    seq_phi_master = []
    bat_phi_master = []
    bat_phi_master = pd.DataFrame(bat_phi_master)

    for i in range(100):
        seq_dir = os.path.join('simple_master2_da_{0}'.format(i))
        bat_dir = os.path.join('simple_master2_ies_{0}'.format(i))

        bat_phi = pd.read_csv(os.path.join(bat_dir, 'freyberg6_run_ies.phi.actual.csv'))
        bat_phi_master = bat_phi_master.append(bat_phi.iloc[3, :])
        # print(bat_phi.iloc[3,:])

        seq_phi = pd.read_csv(os.path.join(seq_dir, 'freyberg6_run_da.global.phi.actual.csv'))
        seq_cyc_mean = 0.
        for i in range(25):
            seq_cyc = seq_phi.loc[seq_phi.cycle == i, :]
            seq_cyc = seq_cyc.loc[seq_cyc.iteration == 3,]
            # print(seq_cyc)
            try:
                seq_cyc_mean += float(seq_cyc['mean'])
            except:
                print('no data assimilated')
        seq_phi_master.append(seq_cyc_mean)

    plt.hist(bat_phi_master.iloc[:, 2] - seq_phi_master, label='BAT-SEQ', alpha=.3)
    plt.hist([bat_phi_master.iloc[:, 2], seq_phi_master], label=['BAT', 'SEQ'])
    plt.legend(loc='upper right')
    plt.ylabel('Freqeuncy')
    plt.xlabel('phi')
    # plt.title("mean phi for batch data assimilation of simple model across 100 complex reals")
    # plt.show()
    plt.savefig('phi_hists.pdf')
    plt.close()


def s_plot():
    complex_dir = os.path.join('complex_master')
    complex_obs = pd.read_csv(os.path.join(complex_dir, 'new_complex_obs.csv'))
    complex_gw_all = []
    complex_tail_all = []
    complex_head_all = []
    simple_bat_gw_all = []
    simple_bat_tail_all = []
    simple_bat_head_all = []
    simple_seq_gw_all = []
    simple_seq_tail_all = []
    simple_seq_head_all = []

    for i in range(100):
        print('working on realization ', i)
        complex_tail = complex_obs.loc[i, 'tailwater_time:20170131']
        complex_head = complex_obs.loc[i, 'headwater_time:20171031']
        complex_gw = complex_obs.loc[i, 'trgw_0_29_5_time:20171031']
        complex_gw_all.append(complex_gw)
        complex_head_all.append(complex_head)
        complex_tail_all.append(complex_tail)

        simple_bat_dir = os.path.join('simple_master2_ies_{0}'.format(i))
        simple_bat = pd.read_csv(os.path.join(simple_bat_dir, 'freyberg6_run_ies.3.obs.csv'))
        simple_bat_gw = simple_bat.loc[:, 'trgw_0_9_1_20171031']
        simple_bat_tail = simple_bat.loc[:, 'tailwater_20170131']
        simple_bat_head = simple_bat.loc[:, 'headwater_20171031']
        simple_bat_tail_all.append(simple_bat_tail)
        simple_bat_head_all.append(simple_bat_head)
        simple_bat_gw_all.append(simple_bat_gw)

        simple_seq_dir = os.path.join('simple_master2_da_{0}'.format(i))
        simple_seq = pd.read_csv(os.path.join(simple_seq_dir, 'freyberg6_run_da.23.0.obs.csv'))
        simple_seq_gw = simple_seq.loc[:, 'head_00_009_001']
        simple_seq_head = simple_seq.loc[:, 'headwater']
        simple_seq_head_all.append(simple_seq_head)
        simple_seq_gw_all.append(simple_seq_gw)

        simple_seq = pd.read_csv(os.path.join(simple_seq_dir, 'freyberg6_run_da.14.0.obs.csv'))
        simple_seq_tail = simple_seq.loc[:, 'tailwater']
        simple_seq_tail_all.append(simple_seq_tail)

    fig, ax = plt.subplots(1, 1)
    for ye, xe in zip(complex_gw_all, simple_bat_gw_all):
        ax.scatter(xe, [ye] * len(xe), color='blue', s=1, label='BAT')
    for ze, le in zip(complex_gw_all, simple_seq_gw_all):
        ax.scatter(le, [ze] * len(le), color='orange', s=1, label='SEQ')
    ax.set_xlim(33.5, 38.5)
    ax.set_ylim(33.5, 38.5)
    mn = min(ax.get_xlim()[0], ax.get_ylim()[0])
    mx = max(ax.get_xlim()[1], ax.get_ylim()[1])
    ax.plot([mn, mx], [mn, mx], 'k')
    ax.set_title('GW_3 Forecast')
    ax.set_xlabel('Simple Forecast (ft)')
    ax.set_ylabel('Complex Forecast (ft)')
    plt.savefig('gw_3_forecast.pdf')
    plt.close(fig)
    # plt.show()

    fig, ax = plt.subplots(1, 1)
    for ye, xe in zip(complex_head_all, simple_bat_head_all):
        ax.scatter(xe, [ye] * len(xe), color='blue', s=1, label='BAT')
    for ze, le in zip(complex_head_all, simple_seq_head_all):
        ax.scatter(le, [ze] * len(le), color='orange', s=1, label='SEQ')
    ax.set_xlim(-1500, 500)
    ax.set_ylim(-1500, 500)
    mn = min(ax.get_xlim()[0], ax.get_ylim()[0])
    mx = max(ax.get_xlim()[1], ax.get_ylim()[1])
    ax.plot([mn, mx], [mn, mx], 'k')
    ax.set_title('Headwater Forecast')
    ax.set_xlabel('Simple Forecast (ft)')
    ax.set_ylabel('Complex Forecast (ft)')
    plt.savefig('headwater_forecast.pdf')
    plt.close(fig)

    fig, ax = plt.subplots(1, 1)
    for ye, xe in zip(complex_tail_all, simple_bat_tail_all):
        ax.scatter(xe, [ye] * len(xe), color='blue', s=1, label='BAT')
    for ze, le in zip(complex_tail_all, simple_seq_tail_all):
        ax.scatter(le, [ze] * len(le), color='orange', s=1, label='SEQ')
    ax.set_xlim(-1750, 0)
    ax.set_ylim(-1750, 0)
    mn = min(ax.get_xlim()[0], ax.get_ylim()[0])
    mx = max(ax.get_xlim()[1], ax.get_ylim()[1])
    ax.plot([mn, mx], [mn, mx], 'k')
    ax.set_title('Tailwater Forecast')
    ax.set_xlabel('Simple Forecast (ft)')
    ax.set_ylabel('Complex Forecast (ft)')
    plt.savefig('tailwater_forecast.pdf')
    plt.close(fig)


def plots_obs_v_sim():
    noptmax = 3
    start_date = pd.to_datetime('20151231', format='%Y%m%d') - pd.to_timedelta(10000,unit="d")
    c_d = 'daily_model_files_master_prior'

    cpst = pyemu.Pst(os.path.join(c_d,"freyberg.pst"))
    # load in obs ensemble
    oe_f = pd.read_csv(os.path.join(c_d, "freyberg.0.obs.csv"), index_col=0)
    oe_f = oe_f.T

    # process obs
    cobs = cpst.observation_data
    cobs.loc[:,"time"] = cobs.time.apply(float)

    with PdfPages(os.path.join("obs_v_sim.pdf")) as pdf:
        for i in range(100):
            da_m_d = os.path.join("seq_monthly_model_files_master_{0}".format(i))
            ies_m_d = os.path.join("monthly_model_files_master_{0}".format(i))
            ies_case = "freyberg"
            da_case = "freyberg"
            if not os.path.exists(da_m_d) or not os.path.exists(ies_m_d):
                break

            ies_pst = pyemu.Pst(os.path.join(ies_m_d, ies_case + ".pst"))
            ies_obs = ies_pst.observation_data
            ies_obs = ies_obs.loc[ies_obs.weight > 0, :]
            ies_obs.loc[:,"time"] = ies_obs.time.apply(float)

            ies_pr_oe = pd.read_csv(os.path.join(ies_m_d, ies_case + ".0.obs.csv"))
            ies_pt_oe = pd.read_csv(os.path.join(ies_m_d, ies_case + ".{0}.obs.csv".format(ies_pst.control_data.noptmax)))

            da_pst = pyemu.Pst(os.path.join(da_m_d, da_case + ".pst"))
            da_obs = da_pst.observation_data.loc[da_pst.nnz_obs_names, :].copy()


            da_pr_dict = {}
            da_pt_dict = {}
            times = ies_obs.time.unique()
            times.sort()
            for cycle,time in enumerate(times):
                print(cycle)
                da_pr_oe = pd.read_csv(os.path.join(da_m_d, da_case + ".{0}.0.obs.csv".format(cycle)))
                da_pr_dict[cycle] = da_pr_oe
                pt_file = os.path.join(da_m_d, da_case + ".{0}.{1}.obs.csv".format(cycle, da_pst.control_data.noptmax))
                if (os.path.exists(pt_file)):
                    da_pt_oe = pd.read_csv(pt_file)
                    da_pt_dict[cycle] = da_pt_oe
                else:
                    print("missing posterior", cycle)
            ies_og_uvals = ies_obs.obgnme.unique()
            ies_og_uvals.sort()

            for og in ies_og_uvals:
                ies_obs_og = ies_obs.loc[ies_obs.obgnme == og, :].copy()
                ies_obs_og.sort_values(by="time", inplace=True)
                da_obs_og = da_obs.loc[da_obs.obgnme == og, :]
                cmplx_obs_og = cobs.loc[cobs.obgnme == og, :].copy()
                cmplx_obs_og.sort_values(by="time", inplace=True)
                assert ies_obs_og.shape[0] > 0
                assert da_obs_og.shape[0] > 0
                assert cmplx_obs_og.shape[0] > 0
                #tms = pd.date_range(start='2015-12-31', periods=731, freq='D')
                #dts = ies_obs_og.datetime.values

                def make_plot(axes):
                    ax = axes[0]
                    ax.set_title("batch DA, complex real {0}, observation location: ".format(i) + og, loc="left")
                    [ax.plot(ies_obs_og.time, ies_pr_oe.loc[idx, ies_obs_og.obsnme], "0.5", alpha=0.5, lw=0.25) for idx in
                     ies_pr_oe.index]
                    [ax.plot(ies_obs_og.time, ies_pt_oe.loc[idx, ies_obs_og.obsnme], "b", alpha=0.5, lw=0.5) for idx in
                     ies_pt_oe.index]
                    ax.plot(ies_obs_og.time, ies_pr_oe.loc[ies_pr_oe.index[0], ies_obs_og.obsnme], "0.5", alpha=0.5,
                            lw=0.1, label="prior real")
                    ax.plot(ies_obs_og.time, ies_pt_oe.loc[ies_pt_oe.index[0], ies_obs_og.obsnme], "b", alpha=0.5,
                            lw=0.1, label="post real")
                    # ax.plot(dts, ies_obs_og.obsval, "r", label="truth")
                    ax.plot(cmplx_obs_og.time, cmplx_obs_og.iloc[:, i], 'r', label="truth")
                    ies_obs_nz = ies_obs_og.loc[ies_obs_og.weight > 0, :]

                    ax.scatter(ies_obs_nz.time, ies_obs_nz.obsval, marker="^", color="r", s=50,
                               zorder=10, label="obs")
                    ax.legend(loc="upper right")
                    ax = axes[1]

                    ax.set_title(
                        "sequential DA, complex real {0}, observation location: ".format(i) + og,
                        loc="left")
                    # ax.plot(dts, ies_obs_og.obsval, "r", label="truth")
                    ax.plot(cmplx_obs_og.time, cmplx_obs_og.iloc[:, i], 'r', label="truth")
                    ax.scatter(ies_obs_nz.datetime.values, ies_obs_nz.obsval, marker="^", color="r",
                               s=50, zorder=10, label="obs")

                    post_labeled = False

                    ax.set_ylim(axes[0].get_ylim())
                    ax.legend(loc="upper right")

                    # prior

                fig, axes = plt.subplots(2, 1, figsize=(8, 8))
                make_plot(axes)
                cycles = list(da_pr_dict.keys())
                cycles.sort()
                for cycle in cycles:
                    da_pr_oe = da_pr_dict[cycle]
                    ax = axes[1]
                    ax.scatter([times[cycle] for _ in range(da_pr_oe.shape[0])],
                               da_pr_oe.loc[:, da_obs_og.obsnme[0]].values,
                               marker=".", color="0.5")
                    ax.set_ylim(axes[0].get_ylim())

                    if cycle in da_pt_dict:
                        da_pt_oe = da_pt_dict[cycle]
                        ax.scatter([times[cycle] for _ in range(da_pt_oe.shape[0])],
                                   da_pt_oe.loc[:, da_obs_og.obsnme[0]].values,
                                   marker=".", color="b")
                        ax.set_ylim(axes[0].get_ylim())
                        # plt.tight_layout()
                pdf.savefig()
                plt.close(fig)
                # plt.show()


def invest():
    m_d = "daily_model_files_master_prior"
    pst = pyemu.Pst(os.path.join(m_d,"freyberg.pst"))
    pe = pyemu.ParameterEnsemble.from_binary(pst=pst,filename=os.path.join(m_d,"prior.jcb"))
    print(pe.shape,pst.npar)
    pst.parameter_data.loc[pe.columns.values,"parval1"] = pe._df.loc["0",:]
    pst.control_data.noptmax = 0
    pst.write(os.path.join(m_d,"test.pst"))
    #pyemu.os_utils.run("pestpp-ies test.pst",cwd=m_d)
    oe = pd.read_csv(os.path.join(m_d,"freyberg.0.obs.csv"),index_col=0)
    print(oe.loc[:,oe.columns.map(lambda x: "cum_mass_usecol:wel" in x)].mean().values)


def test_extract_state_obs(t_d):
    cwd = os.getcwd()
    os.chdir(t_d)
    fnames = extract_state_obs()
    os.chdir(cwd)
    return fnames


def extract_hds_state_obs():
    hds = flopy.utils.HeadFile('freyberg6_freyberg.hds')
    arr = hds.get_data()
    fnames = []
    for k, a in enumerate(arr):
        fname = 'heads_' + str(k) + '.dat'
        np.savetxt(fname, a, fmt='%15.6E')
        fnames.append(fname)
    return fnames

def test_extract_hds_state_obs(t_d):
    cwd = os.getcwd()
    os.chdir(t_d)
    fnames = extract_hds_state_obs()
    os.chdir(cwd)
    return fnames

def extract_HK():
    HK_lay1 = np.loadtxt('freyberg6.npf_k_layer1.txt')
    try:
        HK_lay2 = np.loadtxt('freyberg6.npf_k_layer2.txt')
    except:
        print('no layer 2')
    try:
        HK_lay3 = np.loadtxt('freyberg6.npf_k_layer3.txt')
    except:
        print('no layer 3')
    try:
        HK_lay1 = (HK_lay1 + HK_lay2 + HK_lay3)/3
    except:
        print('just a single layer model in a multilayer world')

    np.savetxt('HK.txt', HK_lay1, fmt='%15.6E')

def test_extract_HK(t_d):
    cwd = os.getcwd()
    os.chdir(t_d)
    fnames = extract_HK()
    os.chdir(cwd)
    return fnames

def extract_conc_state_obs():
    ucn = flopy.utils.HeadFile('freyberg6_trns.ucn', precision="double", text="CONCENTRATION")
    arr = ucn.get_data()
    fnames = []
    for k, a in enumerate(arr):
        fname = 'conc_' + str(k) + '.dat'
        np.savetxt(fname, a, fmt='%15.6E')
        fnames.append(fname)
    return fnames

def test_extract_conc_state_obs(t_d):
    cwd = os.getcwd()
    os.chdir(t_d)
    fnames = extract_conc_state_obs()
    os.chdir(cwd)
    return fnames

def setup_interface(org_ws, num_reals=10):
    """copied from auto_pest.py

    """
    np.random.seed(123456)

    # run mf6
    tmp_ws = org_ws + "_temp"
    if os.path.exists(tmp_ws):
        shutil.rmtree(tmp_ws)
    shutil.copytree(org_ws, tmp_ws)
    #to make sure we get a consistent version of pyemu...
    shutil.copytree("pyemu",os.path.join(tmp_ws,"pyemu"))
    pyemu.os_utils.run("mf6", cwd=tmp_ws)

    # load the mf6 model with flopy to get the spatial reference
    sim = flopy.mf6.MFSimulation.load(sim_ws=tmp_ws)
    m = sim.get_model("freyberg6")

    # work out the spatial rediscretization factor
    redis_fac = m.dis.nrow.data / 40

    # where the pest interface will be constructed
    template_ws = org_ws.replace("_newstress","").replace("_1lyr","").replace("_trnsprt","").replace("_muted_rch","") + "_template"

    # instantiate PstFrom object
    pf = pyemu.utils.PstFrom(original_d=tmp_ws, new_d=template_ws,
                             remove_existing=True,
                             longnames=True, spatial_reference=m.modelgrid,
                             zero_based=False, start_datetime="1-1-2018")

    file_names = test_process_list_files(template_ws)
    for file_name in file_names:
        df = pd.read_csv(os.path.join(template_ws,file_name),index_col=0,delim_whitespace=True)
        names = df.columns.tolist()
        pf.add_observations(file_name, index_cols=["time"], use_cols=names, ofile_skip=0,
                            obsgp=file_name.split(".")[0],
                            prefix=file_name.split(".")[0], ofile_sep=" ")

    pf.add_py_function("workflow.py", "process_list_files()", is_pre_cmd=False)

    # inc,cum = test_process_lst_file(template_ws)
    # df = pf.add_observations("inc.csv",index_cols=["time"],use_cols=inc.columns.tolist(),ofile_sep=",",prefix="inc",obsgp="inc")
    # df = pf.add_observations("cum.csv", index_cols=["time"], use_cols=cum.columns.tolist(), ofile_sep=",", prefix="cum",
    #                          obsgp="cum")
    # pf.add_py_function("workflow.py","process_lst_file()",is_pre_cmd=False)


    # add observations from the sfr observation output file
    df = pd.read_csv(os.path.join(template_ws, "sfr.csv"), index_col=0)
    pf.add_observations("sfr.csv", insfile="sfr.csv.ins", index_cols="time",
                        use_cols=list(df.columns.values),
                        prefix="sfr")

    # add observations from the sft observation output file
    df = pd.read_csv(os.path.join(template_ws, "sft.csv"), index_col=0)
    pf.add_observations("sft.csv", insfile="sft.csv.ins", index_cols="time",
                        use_cols=list(df.columns.values),
                        prefix="sft")

    # add observations for the heads observation output file
    df = pd.read_csv(os.path.join(template_ws, "heads.csv"), index_col=0)
    pf.add_observations("heads.csv", insfile="heads.csv.ins",
                        index_cols="time", use_cols=list(df.columns.values),
                        prefix="hds")

    # add observations for the concs observation output file
    df = pd.read_csv(os.path.join(template_ws, "concs.csv"), index_col=0)
    pf.add_observations("concs.csv", insfile="concs.csv.ins",
                        index_cols="time", use_cols=list(df.columns.values),
                        prefix="cnc")

    # add observations for the simulated K values
    #note that for the complex model, this is the avg of the 3 layers
    pf.add_py_function("workflow.py","extract_HK()", is_pre_cmd=False)
    test_extract_HK(template_ws)
    prefix = "HK_k:1"
    pf.add_observations('HK.txt', prefix=prefix, obsgp=prefix)

    # add observations for simulated hds states
    pf.add_py_function("workflow.py", "extract_hds_state_obs()", is_pre_cmd=False)
    fnames = test_extract_hds_state_obs(template_ws)
    for k, fname in enumerate(fnames):
        prefix = "head_k:{0}".format(k)
        pf.add_observations(fname, prefix=prefix, obsgp=prefix)

    # add observations for simulated conc states
    pf.add_py_function("workflow.py", "extract_conc_state_obs()", is_pre_cmd=False)
    fnames = test_extract_conc_state_obs(template_ws)
    for k, fname in enumerate(fnames):
        prefix = "conc_k:{0}".format(k)
        pf.add_observations(fname, prefix=prefix, obsgp=prefix)


    # write a really simple instruction file to read the MODPATH end point file
    # out_file = "freyberg6_mp_forward.mpend"
    # ins_file = out_file + ".ins"
    # with open(os.path.join(pst_helper.new_model_ws, ins_file), 'w') as f:
    #     f.write("pif ~\n")
    #     f.write("l7 w w w !part_status! w w !part_time!\n")
    #
    # #add observations for particle endpoints
    # pf.add_py_function("workflow.py", "postprocess_mp7_ep_data()", is_pre_cmd=False)

    # the geostruct object for grid-scale parameters
    grid_v = pyemu.geostats.ExpVario(contribution=1.0, a=500)
    grid_gs = pyemu.geostats.GeoStruct(variograms=grid_v)

    # the geostruct object for pilot-point-scale parameters
    pp_v = pyemu.geostats.ExpVario(contribution=1.0, a=2000)
    pp_gs = pyemu.geostats.GeoStruct(variograms=pp_v)

    # the geostruct for recharge grid-scale parameters
    rch_v = pyemu.geostats.ExpVario(contribution=1.0, a=1000)
    rch_gs = pyemu.geostats.GeoStruct(variograms=rch_v)

    # the geostruct for temporal correlation
    #temporal_gs = pyemu.geostats.GeoStruct(variograms=pyemu.geostats.ExpVario(contribution=1.0, a=30))

    # import flopy as part of the forward run process
    pf.extra_py_imports.append('flopy')

    # use the idomain array for masking parameter locations
    try:
        ib = m.dis.idomain[0].array
    except:
        ib = m.dis.idomain[0]

    # define a dict that contains file name tags and lower/upper bound information
    tags = {"npf_k_": [0.2, 5.], "npf_k33_": [.2, 5], "sto_ss": [.5, 2], "sto_sy": [.8, 1.2],
            "rch_recharge": [.5, 1.5]}
    dts = pd.to_datetime("1-1-2018") + \
          pd.to_timedelta(np.cumsum(sim.tdis.perioddata.array["perlen"]), unit="d")

    # loop over each tag, bound info pair
    for tag, bnd in tags.items():
        lb, ub = bnd[0], bnd[1]
        # find all array based files that have the tag in the name
        arr_files = [f for f in os.listdir(template_ws) if tag in f and f.endswith(".txt")]

        if len(arr_files) == 0:
            print("warning: no array files found for ", tag)
            continue

        # make sure each array file in nrow X ncol dimensions (not wrapped)
        # for arr_file in arr_files:
        #     print(ib.shape)
        #     arr = np.loadtxt(os.path.join(template_ws, arr_file)).reshape(ib.shape)
        #     print(arr)
        #     np.savetxt(os.path.join(template_ws, arr_file), arr, fmt="%15.6E")

        # if this is the recharge tag
        if "rch" in tag:
            # add one set of grid-scale parameters for all files
            pf.add_parameters(filenames=arr_files, par_type="grid", par_name_base="rch_gr",
                              pargp="rch_gr", zone_array=ib, upper_bound=ub, lower_bound=lb,
                              geostruct=rch_gs)

            # add one constant parameter for each array, and assign it a datetime
            # so we can work out the temporal correlation
            for arr_file in arr_files:
                arr = np.loadtxt(os.path.join(template_ws, arr_file))
                print(arr_file,arr.mean(),arr.std())
                uub = arr.mean() * ub
                llb = arr.mean() * lb
                if "daily" in org_ws.lower():
                    uub *= 5
                    llb /= 5
                kper = int(arr_file.split('.')[1].split('_')[-1]) - 1
                pf.add_parameters(filenames=arr_file, par_type="constant", par_name_base=arr_file.split('.')[1] + "_cn",
                                  pargp="rch_const", zone_array=ib, upper_bound=uub, lower_bound=llb,
                                  par_style="direct")


        # otherwise...
        else:
            # for each array add both grid-scale and pilot-point scale parameters
            for arr_file in arr_files:
                pf.add_parameters(filenames=arr_file, par_type="grid", par_name_base=arr_file.split('.')[1] + "_gr",
                                  pargp=arr_file.split('.')[1] + "_gr", zone_array=ib, upper_bound=ub, lower_bound=lb,
                                  geostruct=grid_gs)
                pf.add_parameters(filenames=arr_file, par_type="constant", par_name_base=arr_file.split('.')[1] + "_cn",
                                  pargp=arr_file.split('.')[1] + "_cn", zone_array=ib, upper_bound=ub, lower_bound=lb,
                                  geostruct=grid_gs)
                pf.add_parameters(filenames=arr_file, par_type="pilotpoints",
                                  par_name_base=arr_file.split('.')[1] + "_pp",
                                  pargp=arr_file.split('.')[1] + "_pp", zone_array=ib, upper_bound=ub, lower_bound=lb,
                                  pp_space=5, geostruct=pp_gs)

    # add direct pars for the ic strt values
    tag = "ic_strt"
    arr_files = [f for f in os.listdir(template_ws) if tag in f and f.endswith(".txt")]
    for arr_file in arr_files:
        # make sure each array file in nrow X ncol dimensions (not wrapped)
        # arr = np.loadtxt(os.path.join(template_ws, arr_file)).reshape(ib.shape)
        arr = np.loadtxt(os.path.join(template_ws, arr_file))
        # np.savetxt(os.path.join(template_ws, arr_file), arr, fmt="%15.6E")
        k = int(arr_file.split('.')[1][-1]) - 1
        prefix = "head_k:{0}".format(k)
        zn_arr = np.ones_like(arr, dtype=int)
        zn_arr[arr < 0] = 0
        zn_arr[arr > 1000] = 0
        pf.add_parameters(arr_file, par_type="grid", par_style="direct", pargp=prefix, par_name_base=prefix,
                          transform="none",
                          lower_bound=-10000, upper_bound=10000, zone_array=zn_arr)

        # add direct pars for the conc strt values
    tag = "conc"
    arr_files = [f for f in os.listdir(template_ws) if tag in f and f.endswith(".txt")]
    for arr_file in arr_files:
        print(arr_file)
        # make sure each array file in nrow X ncol dimensions (not wrapped)
        # arr = np.loadtxt(os.path.join(template_ws, arr_file)).reshape(ib.shape)
        arr = np.loadtxt(os.path.join(template_ws, arr_file))
        # np.savetxt(os.path.join(template_ws, arr_file), arr, fmt="%15.6E")
        k = int(arr_file.split('.')[1][-1]) - 1
        prefix = "conc_k:{0}".format(k)
        # zn_arr = np.ones_like(arr, dtype=int)
        # zn_arr[arr < 0] = 0
        # zn_arr[arr > 1000] = 0
        pf.add_parameters(arr_file, par_type="grid", par_style="direct", pargp=prefix, par_name_base=prefix,
                                transform="none",
                                lower_bound=-10000, upper_bound=10000, zone_array=zn_arr)

    # get all the list-type files associated with the wel package
    list_files = [f for f in os.listdir(org_ws) if "freyberg6.wel_stress_period_data_" in f and f.endswith(".txt")]
    # for each wel-package list-type file
    for list_file in list_files:
        kper = int(list_file.split(".")[1].split('_')[-1]) - 1
        # add spatially constant, but temporally correlated parameter
        pf.add_parameters(filenames=list_file, par_type="constant", par_name_base="twel_mlt_{0}".format(kper),
                          pargp="twel_mlt".format(kper), index_cols=[0, 1, 2], use_cols=[3],
                          upper_bound=5, lower_bound=0.2, datetime=dts[kper])

    lb,ub = .2,5
    if "daily" in org_ws.lower():
        lb, ub = .1, 10
    # add temporally indep, but spatially correlated grid-scale parameters, one per well
    pf.add_parameters(filenames=list_files, par_type="grid", par_name_base="wel_grid",
                      pargp="wel_grid", index_cols=[0, 1, 2], use_cols=[3],
                      upper_bound=ub, lower_bound=lb)

    # add grid-scale parameters for SFR reach conductance.  Use layer, row, col and reach
    # number in the parameter names
    pf.add_parameters(filenames="freyberg6.sfr_packagedata.txt", par_name_base="sfr_rhk",
                      pargp="sfr_rhk", index_cols=[0, 1, 2, 3], use_cols=[9], upper_bound=20.,
                      lower_bound=0.05,
                      par_type="grid")

    # add model run command
    pf.mod_sys_cmds.append("mf6")

    # pf.add_py_function("workflow.py","remove_ats()",is_pre_cmd=True)

    # add modpath run command
    # pf.mod_sys_cmds.append("mp7 freyberg6_mp_forward.mpsim")

    # build pest control file
    pst = pf.build_pst('freyberg.pst')
    par = pst.parameter_data
    searchfor = ['head_k', 'conc_k']
    strt_pars = par.loc[par.parnme.str.contains('|'.join(searchfor)), "parnme"]

    # set the first stress period to no pumping
    first_wpar = "twel_mlt_0_inst:0_usecol:3"
    assert first_wpar in set(pst.par_names)
    pf.pst.parameter_data.loc[first_wpar,"partrans"] = "fixed"
    pf.pst.parameter_data.loc[first_wpar, "parval1"] = 5.0e-10
    pf.pst.parameter_data.loc[first_wpar, "parlbnd"] = 1.0e-10
    pf.pst.parameter_data.loc[first_wpar, "parubnd"] = 1.0e-9

    # fix the new stress well if present
    new_wpar = par.loc[par.parnme.apply(lambda x: "wel_grid" in x and "idx0:0" in x),"parnme"]
    m_lrc = (1,25,5)
    if "daily" in template_ws:
        m_lrc = (m_lrc[0],m_lrc[1]*3,m_lrc[2]*3)
    if new_wpar.shape[0] > 0:
        new_wpar = "wel_grid_inst:0_usecol:3_idx0:{0}_idx1:{1}_idx2:{2}".format(m_lrc[0]-1,m_lrc[1]-1,m_lrc[2]-1)
        assert new_wpar in set(pf.pst.par_names),new_wpar
    pf.pst.parameter_data.loc[new_wpar, "partrans"] = "fixed"
    pf.pst.parameter_data.loc[new_wpar, "parval1"] = 1.0
    pf.pst.parameter_data.loc[new_wpar, "parlbnd"] = 0.999999
    pf.pst.parameter_data.loc[new_wpar, "parubnd"] = 1.000001

    # draw from the prior and save the ensemble in binary format
    pe = pf.draw(num_reals, use_specsim=True)
    # replace the ic strt pars with the control file values

    print(strt_pars)
    for idx in pe.index:
        pe.loc[idx, strt_pars] = par.loc[strt_pars, "parval1"]
    par.loc[strt_pars, "parchglim"] = "relative"
    pe.to_binary(os.path.join(template_ws, "prior.jcb"))

    # set some algorithmic controls
    pst.control_data.noptmax = 0
    #pst.pestpp_options["additional_ins_delimiters"] = ","

    # ident the obs-par state linkage
    obs = pst.observation_data
    state_obs = obs.loc[obs.obsnme.str.contains("arrobs_head"), :].copy()
    state_par = par.loc[par.parnme.str.contains('d_head'), :].copy()
    for v in ["k", "i", "j"]:
        state_par.loc[:, v] = state_par.loc[:, v].apply(int)
        state_obs.loc[:, v] = state_obs.loc[:, v].apply(int)
    state_par_dict = {"{0}_{1}_{2}".format(k, i, j): n for k, i, j, n in
                      zip(state_par.k, state_par.i, state_par.j, state_par.parnme)}
    obs.loc[:, "state_par_link"] = np.nan
    state_obs.loc[:, "kij"] = state_obs.apply(lambda x: "{0}_{1}_{2}".format(x.k, x.i, x.j), axis=1)

    obs.loc[state_obs.obsnme, "state_par_link"] = state_obs.apply(lambda x: state_par_dict.get((x.kij), np.nan), axis=1)
    print(obs.state_par_link.dropna().shape)

    # ident the obs-par state linkage
    obs = pst.observation_data
    state_obs = obs.loc[obs.obsnme.str.contains("arrobs_conc"), :].copy()
    state_par = par.loc[par.parnme.str.contains('d_conc'), :].copy()
    for v in ["k", "i", "j"]:
        state_par.loc[:, v] = state_par.loc[:, v].apply(int)
        state_obs.loc[:, v] = state_obs.loc[:, v].apply(int)
    state_par_dict = {"{0}_{1}_{2}".format(k, i, j): n for k, i, j, n in
                      zip(state_par.k, state_par.i, state_par.j, state_par.parnme)}
    state_obs.loc[:, "kij"] = state_obs.apply(lambda x: "{0}_{1}_{2}".format(x.k, x.i, x.j), axis=1)

    obs.loc[state_obs.obsnme, "state_par_link"] = state_obs.apply(lambda x: state_par_dict.get((x.kij), np.nan), axis=1)
    print(obs.state_par_link.dropna().shape)

    # write the control file
    pst.write(os.path.join(pf.new_d, "freyberg.pst"),version=2)

    # run with noptmax = 0
    pyemu.os_utils.run("{0} freyberg.pst".format(
        os.path.join("pestpp-ies")), cwd=pf.new_d)

    # make sure it ran
    res_file = os.path.join(pf.new_d, "freyberg.base.rei")
    assert os.path.exists(res_file), res_file
    pst.set_res(res_file)
    print(pst.phi)

    pst.control_data.noptmax = -1

    # define what file has the prior parameter ensemble
    pst.pestpp_options["ies_par_en"] = "prior.jcb"

    # write the updated pest control file
    pst.write(os.path.join(pf.new_d, "freyberg.pst"),version=2)
    return template_ws

def remove_ats():
    ats_flag = np.loadtxt("ats.txt")
    print(ats_flag)
    if ats_flag == 0:
        with open("freyberg6.tdis", 'w') as f:
            f.write("BEGIN Options\n  TIME_UNITS  days\nEND Options\n")
            f.write("BEGIN Dimensions\n  NPER  1\nEND Dimensions\n")
            f.write("BEGIN PERIODDATA\n31.00000000  1       1.00000000\nEND PERIODDATA\n")

        # write a tdis template file
        with open("freyberg6.tdis.tpl", 'w') as f:
            f.write("ptf  ~\n")
            f.write("BEGIN Options\n  TIME_UNITS  days\nEND Options\n")
            f.write("BEGIN Dimensions\n  NPER  1\nEND Dimensions\n")
            f.write("BEGIN PERIODDATA\n~  perlen  ~  1       1.00000000\nEND PERIODDATA\n")


def process_lst_file():
    lst = flopy.utils.Mf6ListBudget("freyberg6.lst")
    inc,cum = lst.get_dataframes(start_datetime=None,diff=True)
    inc.columns = inc.columns.map(str.lower)
    inc.index.name = "time"
    cum.columns = cum.columns.map(str.lower)
    cum.index.name = "time"
    inc.to_csv("inc.csv")
    cum.to_csv("cum.csv")
    return inc,cum

def process_list_files():
    """process the gwf and gwt list files into a format that is actually useful....
        took this lil number from benchmark.py in zp1 to include transport list file components
    """
    class Mf6TListBudget(flopy.utils.mflistfile.ListBudget):
        """"""

        def set_budget_key(self):
            self.budgetkey = "MASS BUDGET FOR ENTIRE MODEL"
            return

    ucn = flopy.utils.HeadFile(os.path.join("freyberg6_trns.ucn"), precision="double", text="CONCENTRATION")
    lst = flopy.utils.Mf6ListBudget("freyberg6.lst")
    fdf,cdf = lst.get_dataframes(diff=True)
    fdf.index = ucn.times
    cdf.index = ucn.times
    cdf.index.name = "time"
    fdf.index.name = "time"
    cdf.fillna(0.0,inplace=True)
    fdf.fillna(0.0,inplace=True)
    fdf.to_csv("inc_flow.dat",sep=" ")
    fdf.to_csv("cum_flow.dat", sep=" ")

    lst = Mf6TListBudget("freyberg6_trns.lst")
    fdf,cdf = lst.get_dataframes(diff=False)
    cdf.index = ucn.times
    fdf.index = ucn.times
    cdf.index.name = "time"
    fdf.index.name = "time"
    cdf.fillna(0.0, inplace=True)
    fdf.fillna(0.0, inplace=True)
    cdf.to_csv("cum_mass.dat", sep=" ")
    fdf.to_csv("inc_mass.dat", sep=" ")
    return ["inc_flow.dat","cum_flow.dat","cum_mass.dat","inc_mass.dat"]

def test_process_lst_file(d):
    cwd = os.getcwd()
    os.chdir(d)
    inc,cum = process_lst_file()
    os.chdir(cwd)
    return inc,cum

def test_process_list_files(d):
    cwd = os.getcwd()
    os.chdir(d)
    inc_flow,cum_flow,cum_mass,inc_mass = process_list_files()
    os.chdir(cwd)
    return inc_flow,cum_flow,cum_mass,inc_mass


def monthly_ies_to_da(org_d,include_est_states=False):
    """convert the batch monthly model to sequential formulation"""

    # load the existing multi-stress period model
    org_sim = flopy.mf6.MFSimulation.load(sim_ws=org_d)

    # setup the destination dir
    t_d = "seq_" + org_d
    if os.path.exists(t_d):
        shutil.rmtree(t_d)
    shutil.copytree(org_d, t_d)

    #first modify the multiplication key file
    mm_df = pd.read_csv(os.path.join(t_d, "mult2model_info.csv"))
    # remove all but the first stress period for time-varying recharge pars
    drop_rch = mm_df.loc[mm_df.model_file.apply(lambda x: "rch_" in x and not "_1." in x), :].index
    mm_df = mm_df.drop(drop_rch)
    # remove all but the first stress period for time-varying well pars
    drop_wel = mm_df.loc[mm_df.model_file.apply(lambda x: ".wel_" in x and not "_1." in x), :].index
    mm_df = mm_df.drop(drop_wel)
    mm_df.to_csv(os.path.join(t_d, "mult2model_info.csv"))

    # modify the tdis
    with open(os.path.join(t_d, "freyberg6.tdis"), 'w') as f:
        # f.write("BEGIN Options\n  TIME_UNITS  days\n  ATS6 FILEIN freyberg6.ats\nEND Options\n")
        f.write("BEGIN Options\n  TIME_UNITS  days\nEND Options\n")
        f.write("BEGIN Dimensions\n  NPER  1\nEND Dimensions\n")
        f.write("BEGIN PERIODDATA\n31.00000000  1       1.00000000\nEND PERIODDATA\n")

    # write a tdis template file
    with open(os.path.join(t_d, "freyberg6.tdis.tpl"), 'w') as f:
        f.write("ptf  ~\n")
        # f.write("BEGIN Options\n  TIME_UNITS  days\n  ATS6 FILEIN freyberg6.ats\nEND Options\n")
        f.write("BEGIN Options\n  TIME_UNITS  days\nEND Options\n")
        f.write("BEGIN Dimensions\n  NPER  1\nEND Dimensions\n")
        f.write("BEGIN PERIODDATA\n~  perlen  ~  1       1.00000000\nEND PERIODDATA\n")
    new_tpl, new_in = [os.path.join(t_d, "freyberg6.tdis.tpl")], [os.path.join(t_d, "freyberg6.tdis")]
    new_tpl_cycle = [-1]

    # make sure it runs...
    pyemu.os_utils.run("mf6", cwd=t_d)

    # load the existing control file
    pst = pyemu.Pst(os.path.join(t_d, "freyberg.pst"))

    # add the new tdis template file
    for tpl, inf, cy in zip(new_tpl, new_in, new_tpl_cycle):
        print(tpl)
        df = pst.add_parameters(tpl, inf, pst_path=".")
        pst.parameter_data.loc[df.parnme, "cycle"] = cy
        pst.parameter_data.loc[df.parnme, "partrans"] = "fixed"

        # write a CNC template file
    with open(os.path.join(t_d, "freyberg6_trns.cnc.tpl"), 'w') as f:
        f.write("ptf  ~\n")
        # f.write("BEGIN Options\n  TIME_UNITS  days\n  ATS6 FILEIN freyberg6.ats\nEND Options\n")
        f.write("BEGIN Options\n  SAVE_FLOWS\nEND Options\n")
        f.write("BEGIN Dimensions\n  MAXBOUND  1\nEND Dimensions\n")
        f.write("BEGIN PERIOD 1\n  1 15 3 ~  concval  ~\nEND PERIOD\n")
    cnc_tpl = os.path.join(t_d, "freyberg6_trns.cnc.tpl")
    cnc_in = os.path.join(t_d, "freyberg6_trns.cnc")
    cnc_tpl_cycle = [0, 1000.0, 1000.0, 1000.0, 1000.0, 1000.0, 1000.0, 1000.0, 1000.0, 1000.0,
                     1000.0, 1000.0, 1000.0, 1000.0, 1000.0, 1000.0, 1000.0, 1000.0, 1000.0,
                     1000.0, 1000.0, 1000.0, 1000.0, 1000.0, 1000.0]

    df = pst.add_parameters(cnc_tpl, cnc_in, pst_path=".")
    pst.parameter_data.loc[df.parnme, "cycle"] = -1
    pst.parameter_data.loc[df.parnme, "partrans"] = "fixed"

    # #add dummy parameter for ats
    # with open(os.path.join(t_d, "ats.txt.tpl"), 'w') as f:
    #     f.write("ptf  ~\n")
    #     f.write("~ats_flag~")
    # ats_tpl = os.path.join(t_d, "ats.txt.tpl")
    # ats_in = os.path.join(t_d, "ats.txt")
    # ats_tpl_cycle = [1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
    #
    # df = pst.add_parameters(ats_tpl, ats_in, pst_path=".")
    # pst.parameter_data.loc[df.parnme, "cycle"] = -1
    # pst.parameter_data.loc[df.parnme, "partrans"] = "fixed"

    # write par cycle table for the perlen par (fixed)
    pers = org_sim.tdis.perioddata.array["perlen"]
    pdf = pd.DataFrame(index=['perlen'], columns=np.arange(org_sim.tdis.nper.data))
    pdf.loc['perlen', :] = pers
    # pdf.loc['ats_flag',:] = ats_tpl_cycle
    pdf.loc['concval',:] = cnc_tpl_cycle
    pdf.to_csv(os.path.join(t_d, "par_cycle_table.csv"))
    pst.pestpp_options["da_parameter_cycle_table"] = "par_cycle_table.csv"

    # save this for later!
    org_obs = pst.observation_data.copy()

    # now drop all existing heads obs package observations since
    # those will be replaced by the state obs
    fname = 'heads.csv'
    fname_ins = fname + ".ins"
    pst.drop_observations(os.path.join(t_d, fname_ins), '.')

    # now drop all existing heads concs package observations since
    # those will be replaced by the state obs
    fname = 'concs.csv'
    fname_ins = fname + ".ins"
    pst.drop_observations(os.path.join(t_d, fname_ins), '.')

    # fix the sfr and list instruction files - only need one output time
    fname = 'sfr.csv'
    ins_names = ["sfr.csv.ins","sft.csv.ins","cum_mass.dat.ins","cum_flow.dat.ins","inc_mass.dat.ins","inc_flow.dat.ins"]
    # ins_names = ["sfr.csv.ins","inc.csv.ins","cum.csv.ins"]
    for ins_name in ins_names:
        pst.drop_observations(os.path.join(t_d, ins_name), '.')

        ins_lines = open(os.path.join(t_d, ins_name), 'r').readlines()
        with open(os.path.join(t_d, ins_name), 'w') as f:
            for line in ins_lines[:3]:
                f.write(line)
        pst.add_observations(os.path.join(t_d, ins_name), pst_path=".")


    if include_est_states:
        # first add gw level sim state pars - just copy the "direct head" par template files
        ic_tpl_files = [f for f in os.listdir(t_d) if ".ic_" in f and f.endswith(".txt.tpl")]
        for ic_tpl_file in ic_tpl_files:
            print(ic_tpl_file)
            lines = open(os.path.join(t_d,ic_tpl_file)).readlines()
            new_tpl = ic_tpl_file.replace(".txt.tpl",".est.txt.tpl")
            with open(os.path.join(t_d,new_tpl),'w') as f:
                for line in lines:
                    f.write(line.replace("d_head","est_d_head"))
            df = pst.add_parameters(os.path.join(t_d,new_tpl),pst_path=".")
            pst.parameter_data.loc[df.parnme,"parval1"] = 40
            pst.parameter_data.loc[df.parnme, "parubnd"] = 60
            pst.parameter_data.loc[df.parnme, "parlbnd"] = 20

        # first add conc sim state pars - just copy the "direct head" par template files
        ic_tpl_files = [f for f in os.listdir(t_d) if ".icc_" in f and f.endswith(".txt.tpl")]
        for ic_tpl_file in ic_tpl_files:
            print(ic_tpl_file)
            lines = open(os.path.join(t_d,ic_tpl_file)).readlines()
            new_tpl = ic_tpl_file.replace(".txt.tpl",".est.txt.tpl")
            with open(os.path.join(t_d,new_tpl),'w') as f:
                for line in lines:
                    f.write(line.replace("d_conc","est_d_conc"))
            df = pst.add_parameters(os.path.join(t_d,new_tpl),pst_path=".")
            pst.parameter_data.loc[df.parnme,"parval1"] = 40
            pst.parameter_data.loc[df.parnme, "parubnd"] = 1000
            pst.parameter_data.loc[df.parnme, "parlbnd"] = 0

        # now add double fake pars for the forecasts just so they are getting
        # estimated one-step-ahead values
        obs = pst.observation_data
        fnames = obs.loc[obs.obsnme.apply(lambda x: "arrobs" not in x),"obsnme"].tolist()
        tpl_file = "double_state_forecast.dat.tpl"
        with open(os.path.join(t_d,tpl_file),'w') as f:
            f.write("ptf ~\n")
            for fname in fnames:
                f.write("{0}  ~   ipar_{0}    ~  ~   est_ipar_{0}    ~\n".format(fname))
        df = pst.add_parameters(os.path.join(t_d,tpl_file),pst_path=".")
        pst.parameter_data.loc[df.parnme, "parval1"] = 0
        pst.parameter_data.loc[df.parnme, "parubnd"] = 1.0e+10
        pst.parameter_data.loc[df.parnme, "parlbnd"] = -1.0e+10
        pst.parameter_data.loc[df.parnme, "partrans"] = "none"
        pst.parameter_data.loc[df.parnme, "parchglim"] = "relative"

        # now set the state_par_link for all these
        # both par data and obs data
        pst.parameter_data.loc[:,"state_par_link"] = np.nan
        par = pst.parameter_data
        est_pars = par.loc[par.parnme.str.startswith("est_"),"parnme"]
        pst.parameter_data.loc[est_pars,"state_par_link"] = est_pars.apply(lambda x: x.replace("est_",""))


        sim_pars = par.loc[par.parnme.str.startswith("ipar_"),"parnme"]
        obs_equiv = sim_pars.apply(lambda x: x.replace("ipar_",""))
        pst.observation_data.loc[obs_equiv.values,"state_par_link"] = sim_pars.values

    # now set cycles
    pst.observation_data.loc[:, 'cycle'] = -1

    # need to set cycle vals and reset the model_file attr
    # for each cycle-specific template files (rch and wel)
    pst.model_input_data.loc[:, "cycle"] = -1
    for i in range(len(pst.model_input_data)):
        # for time-varying well pars
        if pst.model_input_data.iloc[i, 0].startswith('twel_mlt'):
            cy = int(pst.model_input_data.iloc[i, 0].split('_')[2])
            # set the cycle
            pst.model_input_data.iloc[i, 2] = cy
            # point all template files to write the same input file (the first stress period file)
            pst.model_input_data.iloc[i, 1] = pst.model_input_data.iloc[i, 1].replace(str(cy), "0")
        # spatially constant but temporally varying
        elif 'rch_recharge' in pst.model_input_data.iloc[i, 0] and "cn" in pst.model_input_data.iloc[i, 0]:
            cy = int(pst.model_input_data.iloc[i, 0].split('_')[2])
            # set the cycle (one based)
            pst.model_input_data.iloc[i, 2] = cy - 1
            # point all template files to write the same input file (the first stress period file)
            pst.model_input_data.iloc[i, 1] = pst.model_input_data.iloc[i, 1].replace(str(cy), "1")
        # spatially varying but temporal constant
        elif 'rch_recharge' in pst.model_input_data.iloc[i, 0] and pst.model_input_data.iloc[i, 0].endswith(".txt.tpl"):
            cy = int(pst.model_input_data.iloc[i, 0].split('.')[1].split("_")[-1])
            # set the cycle value (one based)
            pst.model_input_data.iloc[i, 2] = cy - 1
            # point all template files to write the same input file (the first stress period file)
            pst.model_input_data.iloc[i, 1] = pst.model_input_data.iloc[i, 1].replace(str(cy), "1")

    # read all instruction files every cycle
    pst.model_output_data.loc[:, "cycle"] = -1

    # need to set the cycle value for all pars - static properties and
    # multi-stress period broadcast forcing pars
    # should get a cycle value of -1.
    pst.parameter_data.loc[:, "cycle"] = -1

    for pname in pst.par_names:
        if pname.startswith('twel_mlt'):
            cy = int(pname.split('_')[2])
            pst.parameter_data.loc[pname, "cycle"] = cy
        elif 'rch_recharge' in pname and "cn" in pname:
            cy = int(pname.split('_')[4])
            pst.parameter_data.loc[pname, "cycle"] = cy - 1

    pst.control_data.noptmax = 3
    pst.pestpp_options["da_verbose_level"] = 3
    pst.write(os.path.join(t_d, "freyberg.pst"), version=2)

    # fill any missing pars in the prior ensemble with ctl file values (esp the perlen par)
    pe = pyemu.ParameterEnsemble.from_binary(pst=pst, filename=os.path.join(t_d, "prior.jcb"))
    pepars = set(pe.columns.tolist())
    pstpars = set(pst.par_names)
    par = pst.parameter_data
    d = pepars.symmetric_difference(pstpars)
    missing = par.loc[par.parnme.apply(lambda x: x in d), "parnme"]
    mvals = par.parval1.loc[missing]
    pe.loc[:, missing] = np.NaN
    for idx in pe.index:
        pe.loc[idx, missing] = mvals
    pe.to_binary(os.path.join(t_d, "prior.jcb"))

    # make a test run
    pst.control_data.noptmax = 0
    pst.write(os.path.join(t_d, "test.pst"), version=2)
    pyemu.os_utils.run("pestpp-da test.pst", cwd=t_d)
    return t_d


def run_batch_seq_prior_monte_carlo(b_d,s_d):
    """run prior monte carlo for the batch and seq monthly models
    """
    pst = pyemu.Pst(os.path.join(s_d,"freyberg.pst"))
    pst.control_data.noptmax = -1
    pst.write(os.path.join(s_d,"freyberg.pst"),version=2)
    m_s_d = s_d.replace("template", "master_prior")
    pyemu.os_utils.start_workers(s_d, "pestpp-da", "freyberg.pst", num_workers=40,
                                 master_dir=m_s_d)

    pst = pyemu.Pst(os.path.join(b_d, "freyberg.pst"))
    pst.control_data.noptmax = -1
    pst.write(os.path.join(b_d, "freyberg.pst"),version=2)
    m_b_d = b_d.replace("template", "master_prior")
    pyemu.os_utils.start_workers(b_d, "pestpp-ies", "freyberg.pst", num_workers=40,
                                 master_dir=m_b_d)

    return m_b_d,m_s_d



def plot_prior_mc_all():
    """plot the prior monte carlo results for daily, monthly batch and monthly sequential

    """

    ognames = keep.copy()
    ognames.extend(forecast)
    ognames.sort()
    label_dict = keep_dict
    label_dict.update(forecast_dict)
    unit_dict = {n: u for n, u in zip(forecast, forecast_units)}
    unit_dict.update({n: u for n, u in zip(keep, keep_units)})
    print(unit_dict)


    c_m_d = "daily_model_files_master_prior"
    s_b_m_d = "monthly_model_files_master_prior"
    s_s_m_d = "seq_monthly_model_files_master_prior"

    top = np.loadtxt(os.path.join(s_b_m_d,"freyberg6.dis_top.txt"))


    c_pst = pyemu.Pst(os.path.join(c_m_d, "freyberg.pst"))
    obs = c_pst.observation_data
    cobs = obs.loc[obs.obsnme.str.startswith("hds_usecol:arrobs_head_"),:]
    cobs.loc[:,"time"] = cobs.time.apply(float)

    s_b_pst = pyemu.Pst(os.path.join(s_b_m_d, "freyberg.pst"))
    c_oe = pd.read_csv(os.path.join(c_m_d, "freyberg.0.obs.csv"), index_col=0)
    log = c_oe.columns.map(lambda x: "mass" in x)
    c_oe.loc[:,log] = c_oe.loc[:,log].apply(np.log10)
    s_b_oe = pd.read_csv(os.path.join(s_b_m_d, "freyberg.0.obs.csv"), index_col=0)
    log = s_b_oe.columns.map(lambda x: "mass" in x)
    s_b_oe.loc[:, log] = s_b_oe.loc[:, log].apply(np.log10)

    s_s_pst = pyemu.Pst(os.path.join(s_s_m_d,"freyberg.pst"))
    seq_oe_files = [f for f in os.listdir(s_s_m_d) if f.endswith(".oe.csv") and "global" in f and f.startswith("freyberg")]
    s_s_oe_dict = {int(f.split(".")[2]):pd.read_csv(os.path.join(s_s_m_d,f),index_col=0) for f in seq_oe_files}
    for key,df in s_s_oe_dict.items():
        log = df.columns.map(lambda x: "mass" in x)
        df.loc[:, log] = df.loc[:, log].apply(np.log10)
        s_s_oe_dict[key] = df
    #oct = pd.read_csv(os.path.join(s_s_m_d,"obs_cycle_tbl.csv"),index_col=0)

    #ognames = list(cobs.obgnme.unique())
    #ognames.sort()
    #ognames.append("sfr_usecol:gage_1")


    with PdfPages("prior_lst_budget_inc.pdf") as pdf:
        c_lst_obs = c_pst.observation_data.copy()
        c_lst_obs = c_lst_obs.loc[c_lst_obs.obsnme.apply(lambda x: x.startswith("inc")),:]
        c_lst_obs.loc[:,"time"] = c_lst_obs.time.apply(float)
        b_lst_obs = s_b_pst.observation_data.copy()
        b_lst_obs = b_lst_obs.loc[b_lst_obs.obsnme.apply(lambda x: x.startswith("inc")), :]
        b_lst_obs.loc[:, "time"] = b_lst_obs.time.apply(float)
        lst_grps = c_lst_obs.obgnme.unique().tolist()
        lst_grps.sort()
        for lst_grp in lst_grps:
            print(lst_grp)
            bgobs = b_lst_obs.loc[b_lst_obs.obgnme==lst_grp,:].copy()
            fig, ax = plt.subplots(1, 1, figsize=(8, 8))
            inc = bgobs.loc[bgobs.obsnme.str.startswith("inc"),:].copy()
            inc.sort_values(by="time",inplace=True)
            [ax.plot(inc.time,s_b_oe.loc[idx,inc.obsnme],"0.5", alpha=0.1) for idx in s_b_oe.index]
            cgobs = c_lst_obs.loc[c_lst_obs.obgnme == lst_grp, :].copy()
            inc = cgobs.loc[cgobs.obsnme.str.startswith("inc"), :].copy()
            inc.sort_values(by="time", inplace=True)
            [ax.plot(inc.time, c_oe.loc[idx, inc.obsnme], "b--",alpha=0.1) for idx in c_oe.index]
            ax.set_title(lst_grp)
            pdf.savefig()
            plt.close(fig)

    with PdfPages("prior_lst_budget_cum.pdf") as pdf:
        c_lst_obs = c_pst.observation_data.copy()
        c_lst_obs = c_lst_obs.loc[c_lst_obs.obsnme.apply(lambda x: x.startswith("cum")),:]
        c_lst_obs.loc[:,"time"] = c_lst_obs.time.apply(float)
        b_lst_obs = s_b_pst.observation_data.copy()
        b_lst_obs = b_lst_obs.loc[b_lst_obs.obsnme.apply(lambda x: x.startswith("cum")), :]
        b_lst_obs.loc[:, "time"] = b_lst_obs.time.apply(float)
        lst_grps = c_lst_obs.obgnme.unique().tolist()
        lst_grps.sort()
        for lst_grp in lst_grps:
            print(lst_grp)
            bgobs = b_lst_obs.loc[b_lst_obs.obgnme==lst_grp,:].copy()
            fig, ax = plt.subplots(1, 1, figsize=(8, 8))
            inc = bgobs.loc[bgobs.obsnme.str.startswith("cum"),:].copy()
            inc.sort_values(by="time",inplace=True)
            [ax.plot(inc.time,s_b_oe.loc[idx,inc.obsnme],"0.5", alpha=0.1) for idx in s_b_oe.index]
            cgobs = c_lst_obs.loc[c_lst_obs.obgnme == lst_grp, :].copy()
            inc = cgobs.loc[cgobs.obsnme.str.startswith("cum"), :].copy()
            inc.sort_values(by="time", inplace=True)
            [ax.plot(inc.time, c_oe.loc[idx, inc.obsnme], "b--",alpha=0.1) for idx in c_oe.index]
            ax.set_title(lst_grp)
            pdf.savefig()
            plt.close(fig)


    with PdfPages("prior_obs_v_sim.pdf") as pdf:
        fig, axes = plt.subplots(len(ognames), 1, figsize=(8, 10))
        for iax,ogname in enumerate(ognames):
            ax = axes[iax]
            cgobs = c_pst.observation_data.loc[c_pst.observation_data.obsnme.str.contains(ogname),:].copy()
            sgobs = s_b_pst.observation_data.loc[s_b_pst.observation_data.obsnme.str.contains(ogname),:].copy()
            k0ogname = ogname
            if sgobs.shape[0] == 0:
                k0ogname = ogname.replace("k:2","k:0")
                sgobs = s_b_pst.observation_data.loc[s_b_pst.observation_data.obsnme.str.contains(k0ogname), :].copy()
            sgobs.loc[:,"time"] = sgobs.time.apply(float)
            cgobs.loc[:, "time"] = cgobs.time.apply(float)

            sgobs.sort_values(by="time", inplace=True)
            cgobs.sort_values(by="time",inplace=True)




            seq_name = k0ogname
            if "arrobs" not in k0ogname:
                seq_name = k0ogname + "_time:10000.0"
            print(seq_name)
            for itime,time in enumerate(sgobs.time):
                if itime in s_s_oe_dict:
                    oe = s_s_oe_dict[itime]
                    #print(oe.loc[:,seq_name])
                    ax.scatter([time for _ in range(oe.shape[0])], oe.loc[:, seq_name], marker=".",color="0.5",alpha=0.5)

            [ax.plot(cgobs.time, c_oe.loc[idx, cgobs.obsnme], "b", lw=0.01, alpha=0.5) for idx in c_oe.index]

            [ax.plot(sgobs.time, s_b_oe.loc[idx, sgobs.obsnme], "0.5", lw=0.01, alpha=0.5) for idx in s_b_oe.index]


            if "arrobs" in ogname:
                i = sgobs.i.apply(int)[0]
                j = sgobs.j.apply(int)[0]
                t = top[i, j]
                ax.plot(ax.get_xlim(),[t,t],"k--",lw=3)
            ax.set_title("{0}) {1}".format(string.ascii_uppercase[iax], label_dict[ogname]),loc="left")
            ax.set_ylabel(unit_dict[ogname])

            #if "gage" not in ogname:
            #    ax.set_ylim(30,ax.get_ylim()[1])
        for ax in axes[:-1]:
            ax.set_xticklabels([])
        axes[-1].set_xlabel("time ($days$)")
        plt.tight_layout()
        pdf.savefig()
        plt.close(fig)






def map_complex_to_simple_bat(c_d,b_d,real_idx):
    """map the daily model outputs to the monthly model observations

    Args:
        c_d (str): the daily (complex) model prior monte carlo master dir
        b_d (str): the monthly (simple) batch model tmeplate dir
        real_idx (int): the complex model prior monte carlo obs ensemble index position
            to use as the observations for the monthly model

    """

    # load the daily model prior MC obs ensemble and get the simulated values
    # coe = pd.read_csv(os.path.join(c_d,"freyberg.0.obs.csv"),index_col=0)
    coe = pd.read_csv(os.path.join(c_d, "freyberg.0.obs.csv"))
    coe.set_index('real_name')
    cvals = coe.loc[real_idx,:]

    # load the batch monthly model control file
    bpst = pyemu.Pst(os.path.join(b_d,"freyberg.pst"))
    obs = bpst.observation_data
    is_1_lay = True
    if True in [True if "k:2" in o else False for o in bpst.obs_names]:
        is_1_lay = False

    #need to do some trickery for complex gw obs in layer 3 to pretend to be in layer 1
    #cvals.index = [item.replace('hds_usecol:arrobs_head_k:0_i:2_j:9', 'hds_usecol:arrobs_head_k:2_i:2_j:9') for item in cvals.index]
    #cvals.index = [item.replace('hds_usecol:arrobs_head_k:0_i:33_j:7', 'hds_usecol:arrobs_head_k:2_i:33_j:7') for item in cvals.index]
    #cvals.index = [item.replace('hds_usecol:arrobs_head_k:2_i:2_j:9', 'hds_usecol:arrobs_head_k:0_i:2_j:9') for item in
    #               cvals.index]
    #cvals.index = [item.replace('hds_usecol:arrobs_head_k:2_i:33_j:7', 'hds_usecol:arrobs_head_k:0_i:33_j:7') for item
    #               in cvals.index]
    # assign all common observations
    if is_1_lay:
        #idx = cvals.index.copy()
        for k in keep:
            if "k:2" not in k:
                continue
            d = cvals.index.map(lambda x: k not in x)
            cvals = cvals.loc[d]
            kk = k.replace("k:2","k:0")
            cvals.index = cvals.index.map(lambda x: x.replace("k:2","k:0") if kk in x else x)

        for k in forecast:
            if "k:2" not in k:
                continue
            d = cvals.index.map(lambda x: k not in x)
            cvals = cvals.loc[d]
            kk = k.replace("k:2", "k:0")
            cvals.index = cvals.index.map(lambda x: x.replace("k:2", "k:0") if kk in x else x)
    print(cvals)
    d = set(bpst.obs_names) - set(cvals.index.tolist())
    print(d)
    obs.loc[:,"obsval"] = cvals.loc[bpst.obs_names]
    assert obs.obsval.shape[0] == obs.obsval.dropna().shape[0]
    assert cvals.loc[bpst.obs_names].shape[0] == cvals.loc[bpst.obs_names].dropna().shape[0]

    # set some weights - only the first 12 months (not counting spin up time)
    # just picked weights that seemed to balance-ish the one realization I looked at...
    obs.loc[:,"weight"] = 0.0
    for k in keep.copy():
        if is_1_lay:
            k = k.replace("k:2","k:0")
        kobs = obs.loc[obs.obsnme.str.contains(k),:].copy()
        kobs.loc[:,"time"] = kobs.time.apply(float)
        kobs.sort_values(by="time",inplace=True)
        if "gage" in k:
            obs.loc[kobs.obsnme[1:13], "weight"] = 0.005
        else:
            obs.loc[kobs.obsnme[1:13],"weight"] = 2.0

    assert bpst.nnz_obs == 48,bpst.nnz_obs

    # setup a template dir for this complex model realization
    t_d = b_d + "_{0}".format(real_idx)
    if os.path.exists(t_d):
        shutil.rmtree(t_d)
    shutil.copytree(b_d,t_d)
    bpst.write(os.path.join(t_d,"freyberg.pst"),version=2)
    return t_d

def map_simple_bat_to_seq(b_d,s_d):
    """map obs vals and weights from monthly batch to sequential

    Args:
        b_d (str): batch template dir with obs vals and weights set
        s_d (str): generic sequential template dir


    """
    org_sim = flopy.mf6.MFSimulation.load(sim_ws=b_d)

    # load the batch control file
    bpst = pyemu.Pst(os.path.join(b_d, "freyberg.pst"))
    bobs = bpst.observation_data
    # work on the sequential obs names - yuck
    seq_names = [o.split("_time")[0].replace("hds_usecol:","") for o in bpst.nnz_obs_names]
    seq_names = [n+"_time:10000.0" if "sfr" in n else n for n in seq_names]
    seq_names = set(seq_names)
    assert len(seq_names) == len(keep)

    # load the sequential control file
    spst = pyemu.Pst(os.path.join(s_d, "freyberg.pst"))
    spst.observation_data.loc[:,"weight"] = 0.0
    snames = set(spst.obs_names)
    # assert len(seq_names - snames) == 0
    seq_names = list(seq_names)
    seq_names.sort()

    sobs = spst.observation_data.loc[seq_names]

    # build the weight and obs cycle tables
    odf = pd.DataFrame(columns=np.arange(org_sim.tdis.nper.data),index=seq_names)
    wdf = pd.DataFrame(columns=np.arange(org_sim.tdis.nper.data), index=seq_names)
    wdf.loc[:,:] = 0.0
    for seq_name in seq_names:
        print(seq_name)
        # the batch obs info for this sequential obs name
        bsobs = bobs.loc[bobs.obsnme.str.contains(seq_name.replace("_time:10000.0","")),:].copy()
        # sort by float time
        bsobs.loc[:,"time"] = bsobs.time.apply(float)
        bsobs.sort_values(by="time",inplace=True)

        # set the values for the 2nd thru 13th stress period/cycle
        odf.loc[seq_name,np.arange(1,13)] = bsobs.obsval.iloc[1:13].values
        wdf.loc[seq_name, np.arange(1,13)] = bsobs.weight.iloc[1:13].values
        spst.observation_data.loc[seq_name,"weight"] = 1.0 #just a generic weight, overridden by weight table

    # some sanity checks
    print(odf.iloc[0,1])
    print(bobs.loc[bpst.nnz_obs_names[1],"obsval"])
    assert odf.iloc[0,1] == bobs.loc[bpst.nnz_obs_names[0],"obsval"]
    assert wdf.iloc[0, 1] == bobs.loc[bpst.nnz_obs_names[0], "weight"]

    # setup the seq template for this complex realization
    t_d = "seq_" + b_d
    if os.path.exists(t_d):
        shutil.rmtree(t_d)
    shutil.copytree(s_d,t_d)
    # save and point to the cycle tables
    odf.to_csv(os.path.join(t_d,"obs_cycle_tbl.csv"))
    wdf.to_csv(os.path.join(t_d,"weight_cycle_tbl.csv"))
    spst.pestpp_options["da_weight_cycle_table"] = "weight_cycle_tbl.csv"
    spst.pestpp_options["da_observation_cycle_table"] = "obs_cycle_tbl.csv"

    # test run with noptmax = 0
    spst.control_data.noptmax = 0
    spst.pestpp_options["da_verbose_level"] = 3
    spst.write(os.path.join(t_d,"freyberg.pst"),version=2)
    pyemu.os_utils.run("pestpp-da freyberg.pst",cwd=t_d)
    return t_d


def plot_obs_v_sim2(subdir=".",post_iter=None):
    """plot the results for daily, monthly batch and monthly sequential

    """
    c_m_d = "daily_model_files_master_prior"
    c_pst = pyemu.Pst(os.path.join(c_m_d, "freyberg.pst"))
    cobs = c_pst.observation_data
    #cobs = obs.loc[obs.obsnme.str.startswith("hds_usecol:arrobs_head_"), :]
    cobs.loc[:, "time"] = cobs.time.apply(float)
    c_oe = pd.read_csv(os.path.join(c_m_d, "freyberg.0.obs.csv"), index_col=0)
    cw_cols = c_oe.columns.map(lambda x: "mass" in x)
    c_oe.loc[:,cw_cols] = c_oe.loc[:,cw_cols].apply(np.log10)
    pname = os.path.join(subdir,"obs_v_sim.pdf")
    if post_iter is not None:
        pname = os.path.join(subdir,"obs_v_sim_postier_{0}.pdf".format(post_iter))

    with PdfPages(pname) as pdf:
        for ireal in range(100):
            s_b_m_d = os.path.join(subdir,"monthly_model_files_master_{0}".format(ireal))
            s_s_m_d = os.path.join(subdir,"seq_monthly_model_files_master_{0}".format(ireal))
            if not os.path.exists(s_s_m_d) or not os.path.exists(s_b_m_d):
                break
            try:
                s_b_pst = pyemu.Pst(os.path.join(s_b_m_d, "freyberg.pst"))
                s_b_oe_pr = pd.read_csv(os.path.join(s_b_m_d, "freyberg.0.obs.csv"), index_col=0)
                sw_cols = s_b_oe_pr.columns.map(lambda x: "mass" in x)
                s_b_oe_pr.loc[:,sw_cols] = s_b_oe_pr.loc[:,sw_cols].apply(np.log10)
                bpost_iter = s_b_pst.control_data.noptmax
                if post_iter is not None:
                    bpost_iter = post_iter
                s_b_oe_pt = pd.read_csv(os.path.join(s_b_m_d, "freyberg.{0}.obs.csv".format(bpost_iter)),
                                        index_col=0)
                s_b_oe_pt.loc[:, sw_cols] = s_b_oe_pt.loc[:, sw_cols].apply(np.log10)
                s_s_pst = pyemu.Pst(os.path.join(s_s_m_d,"freyberg.pst"))
                seq_oe_files_pr = [f for f in os.listdir(s_s_m_d) if f.endswith("0.obs.csv") and f.startswith("freyberg")]
                spost_iter = s_s_pst.control_data.noptmax
                if post_iter is not None:
                    spost_iter = post_iter
                seq_oe_files_pt = [f for f in os.listdir(s_s_m_d) if
                                   f.endswith("{0}.obs.csv".format(spost_iter)) and f.startswith("freyberg")]

                s_s_oe_dict_pr = {int(f.split(".")[1]):pd.read_csv(os.path.join(s_s_m_d,f),index_col=0) for f in seq_oe_files_pr}
                s_s_oe_dict_pt = {int(f.split(".")[1]): pd.read_csv(os.path.join(s_s_m_d, f), index_col=0) for f in
                                  seq_oe_files_pt}
                for key,df in s_s_oe_dict_pr.items():
                    log_cols = df.columns.map(lambda x: "mass" in x)
                    df.loc[:,log_cols] = df.loc[:,log_cols].apply(np.log10)
                    s_s_oe_dict_pr[key] = df
                for key, df in s_s_oe_dict_pt.items():
                    log_cols = df.columns.map(lambda x: "mass" in x)
                    df.loc[:, log_cols] = df.loc[:, log_cols].apply(np.log10)
                    s_s_oe_dict_pt[key] = df

            except:
                break


            ognames = keep.copy()
            ognames.extend(forecast)
            label_dict = keep_dict.copy()
            label_dict.update(forecast_dict)


            is_1_lay = True
            if True in [True if "k:2" in o else False for o in s_b_pst.obs_names]:
                is_1_lay = False

            for ogname in ognames:
                #if "cum" not in ogname:
                #    continue
                k0name = ogname
                if is_1_lay:
                    k0ogname = ogname.replace("k:2","k:0")
                fig, axes = plt.subplots(2, 1, figsize=(8, 8))
                cgobs = cobs.loc[cobs.obsnme.str.contains(ogname),:].copy()
                sgobs = s_b_pst.observation_data.loc[s_b_pst.observation_data.obsnme.str.contains(k0ogname),:].copy()
                if cgobs.shape[0] == 0:
                   raise Exception("complex empty "+ogname)
                if sgobs.shape[0] == 0:
                    raise Exception("batch empty "+k0ogname+","+ogname)

                sgobs.loc[:,"time"] = sgobs.time.apply(float)
                cgobs.loc[:, "time"] = cgobs.time.apply(float)
                sgnzobs = sgobs.loc[sgobs.weight > 0, :].copy()

                sgobs.sort_values(by="time", inplace=True)
                cgobs.sort_values(by="time",inplace=True)


                ax = axes[0]
                ax.set_title("A) batch formulation {0}, replicate {1}".format(label_dict[ogname],c_oe.index[ireal]),loc="left")


                [ax.plot(sgobs.time,s_b_oe_pr.loc[idx,sgobs.obsnme],"0.5",lw=0.01,alpha=0.5) for idx in s_b_oe_pr.index]
                [ax.plot(sgobs.time, s_b_oe_pt.loc[idx, sgobs.obsnme], "b", lw=0.01, alpha=0.5) for idx in
                 s_b_oe_pt.index]
                ax.plot(cgobs.time, c_oe.loc[c_oe.index[ireal], cgobs.obsnme], "r", lw=2.0, alpha=0.85)
                ax.scatter(sgnzobs.time, sgnzobs.obsval, marker="^", color="r")
                ax = axes[1]

                seq_name = k0ogname
                if "arrobs" not in k0ogname:
                    seq_name = k0ogname + "_time:10000.0"
                print(ireal,seq_name)
                for itime,time in enumerate(sgobs.time):
                    #itime += 1


                    if itime in s_s_oe_dict_pr:
                        oe = s_s_oe_dict_pr[itime]
                        #print(oe.loc[:,seq_name])
                        ax.scatter([time for _ in range(oe.shape[0])], oe.loc[:, seq_name], marker=".",color="0.5",alpha=0.5)
                    if itime in s_s_oe_dict_pt:
                        oe = s_s_oe_dict_pt[itime]
                        # print(oe.loc[:,seq_name])
                        ax.scatter([time for _ in range(oe.shape[0])], oe.loc[:, seq_name], marker=".", color="b",
                                   alpha=0.5)
                ax.plot(cgobs.time, c_oe.loc[c_oe.index[ireal], cgobs.obsnme], "r", lw=2.0, alpha=0.85)
                ax.scatter(sgnzobs.time, sgnzobs.obsval, marker="^", color="r")
                ax.set_title("B) sequential formulation {0}, replicate {1}".format(label_dict[ogname], c_oe.index[ireal]),loc="left")
                #if "gage" not in ogname:
                #    ax.set_ylim(30,ax.get_ylim()[1])
                mn = 1.0e+10
                mx = -1.0e+10
                for ax in axes.flatten():
                    mn = min(ax.get_ylim()[0],mn)
                    mx = max(ax.get_ylim()[1],mx)
                for ax in axes.flatten():
                    ax.set_ylim(mn,mx)
                plt.tight_layout()
                pdf.savefig()
                plt.close(fig)


def plot_domain():
    model_ws = "monthly_model_files_master_0"
    sim = flopy.mf6.MFSimulation.load(sim_ws=model_ws)
    m = sim.get_model("freyberg6")
    wel_data = m.wel.stress_period_data.array[0]
    sfr_data = m.sfr.packagedata.array
    ghb_data = m.ghb.stress_period_data.array[0]
    ib = m.dis.idomain.array[0, :, :]
    ib_cmap = plt.get_cmap("Greys_r")
    ib_cmap.set_bad(alpha=0.0)

    pst = pyemu.Pst(os.path.join(model_ws,"freyberg.pst"))

    fig,ax = plt.subplots(1,1,figsize=(6,5))

    ib = np.ma.masked_where(ib!=0,ib)

    ax.imshow(ib,cmap=ib_cmap,extent=m.modelgrid.extent)
    for ii,cid in enumerate(wel_data.cellid):
        i,j = cid[1],cid[2]
        verts = m.modelgrid.get_cell_vertices(i, j)
        if ii == wel_data.shape[0] - 1:
            p = Polygon(verts, facecolor='b',label="extraction well")
        else:
            p = Polygon(verts, facecolor='b')
        ax.add_patch(p)
    for ii,cid in enumerate(ghb_data.cellid):
        i,j = cid[1],cid[2]
        verts = m.modelgrid.get_cell_vertices(i, j)
        if ii == ghb_data.shape[0] - 1:
            p = Polygon(verts, facecolor='m',label="GHB cell")
        else:
            p = Polygon(verts, facecolor='m')
        ax.add_patch(p)
    for cid in sfr_data.cellid:
        i, j = cid[1], cid[2]
        verts = m.modelgrid.get_cell_vertices(i, j)
        if i < 20:
            c = "g"
        else:
            c = "c"
        if i == 19:
            p = Polygon(verts, facecolor=c, label="SFR headwater reaches")
        elif i == 39:
            p = Polygon(verts, facecolor=c, label="SFR tailwater reaches")
        else:
            p = Polygon(verts, facecolor=c)
        ax.add_patch(p)
    #ax.text(m.modelgrid.xcellcenters[10,j],m.modelgrid.ycellcenters[10,j],"headwater sw-gw exchange forecast reaches",rotation=90,ha="center",va="center")
    #spnspecs.add_text(ax=ax,x=m.modelgrid.xcellcenters[10,j],y=m.modelgrid.ycellcenters[10,j],text="headwater sw-gw exchange forecast reaches",
    #                  rotation=90,ha="center",va="center",bold=False, italic=False,transform=False,bbox={"facecolor":"none","edgecolor":"none"})
    #ax.text(m.modelgrid.xcellcenters[30, j], m.modelgrid.ycellcenters[30, j], "tailwater sw-gw exchange forecast reaches", rotation=90, ha="center",
    #        va="center")
    #spnspecs.add_text(ax=ax, x=m.modelgrid.xcellcenters[30, j], y=m.modelgrid.ycellcenters[30, j],
    #                  text="tailwater sw-gw exchange forecast reaches", rotation=90, ha="center", va="center")

    x = m.modelgrid.xcellcenters[i, j]
    y = m.modelgrid.ycellcenters[i, j]
    ax.scatter([x],[y],marker="^",c='r',s=100,zorder=10,label="point observation/forecast location")
    ax.text(x + 150, y+150, "sw_1".format(1), zorder=11, bbox=dict(facecolor='w', alpha=1,edgecolor="none",pad=1))
    ylim = ax.get_ylim()

    nz_obs = pst.observation_data.loc[pst.nnz_obs_names,:].copy()
    nz_obs = nz_obs.loc[nz_obs.obsnme.str.contains("arrobs"),:]

    nz_obs.loc[:, "i"] = nz_obs.i.apply(int)
    nz_obs.loc[:, "j"] = nz_obs.j.apply(int)
    for ii,g in enumerate(nz_obs.obgnme.unique()):
        i = nz_obs.loc[nz_obs.obgnme == g, "i"][0]
        j = nz_obs.loc[nz_obs.obgnme == g, "j"][0]
        x = m.modelgrid.xcellcenters[i, j]
        y = m.modelgrid.ycellcenters[i, j]
        ax.scatter([x], [y], marker="^", c='r', s=100, zorder=10)
        ax.text(x + 150, y+150,"gw_{0}".format(ii+1),zorder=11, bbox=dict(facecolor='w', alpha=1,edgecolor="none",pad=1))

    gw_fore = [f for f in forecast if "arrobs" in f][0]
    i = int(gw_fore.split('_')[3].split(":")[1])
    j = int(gw_fore.split('_')[4].split(":")[1])
    x = m.modelgrid.xcellcenters[i, j]
    y = m.modelgrid.ycellcenters[i, j]
    ax.scatter([x], [y], marker="^", c='r', s=100, zorder=10)
    ax.text(x + 150, y + 150, forecast_dict[gw_fore], zorder=11,
            bbox=dict(facecolor='w', alpha=1, edgecolor="none", pad=1))

    top = m.dis.top.array.reshape(ib.shape)
    top[top<0] = np.NaN
    cb = ax.imshow(top,extent=m.modelgrid.extent,cmap="bone")
    cb = plt.colorbar(cb,pad=0.01)
    cb.set_label("top $m$")
    ax.set_xlabel("x $m$")
    ax.set_ylabel("y $m$")

    ax.set_ylim(0,ylim[1])
    plt.tight_layout()
    plt.savefig("domain.pdf")
    plt.close("all")


def plot_s_vs_s(summarize=False, subdir=".", post_iter=None):
    include_est_states = False
    ognames = keep
    ognames.extend(forecast)
    label_dict = keep_dict
    label_dict.update(forecast_dict)

    # first rip thru all the dirs and load...
    s_b_dict = {}
    s_s_dict = {}
    s_s_est_dict = {}
    print("loading results...")

    for ireal in range(100):
        s_b_m_d = os.path.join(subdir,"monthly_model_files_master_{0}".format(ireal))
        s_s_m_d = os.path.join(subdir,"seq_monthly_model_files_master_{0}".format(ireal))

        if not os.path.exists(s_s_m_d) or not os.path.exists(s_s_m_d):
            break
        try:
            s_b_pst = pyemu.Pst(os.path.join(s_b_m_d, "freyberg.pst"))
            obs = s_b_pst.observation_data.loc[s_b_pst.nnz_obs_names, :]
            obs = obs.loc[obs.obsnme.str.contains("sfr_usecol:gage_1")]
            if obs.obsval.max() > 8000:
                continue
            obs = s_b_pst.observation_data
            obs = obs.loc[obs.obsnme.str.contains("hds_usecol:arrobs_head_k:0_i:33_j:7"),:]
            if obs.obsval.min() < 30:
                continue
            s_b_oe_pr = pd.read_csv(os.path.join(s_b_m_d, "freyberg.0.obs.csv"), index_col=0)
            log_cols = s_b_oe_pr.columns.map(lambda x: "mass" in x or "cnc" in x)
            s_b_oe_pr.loc[:,log_cols] = s_b_oe_pr.loc[:,log_cols].apply(np.log10)

            bpost_iter = s_b_pst.control_data.noptmax
            if post_iter is not None:
                bpost_iter = post_iter
            s_b_oe_pt = pd.read_csv(os.path.join(s_b_m_d, "freyberg.{0}.obs.csv".format(bpost_iter)),
                                    index_col=0)
            log_cols = s_b_oe_pt.columns.map(lambda x: "mass" in x or "cnc" in x)
            s_b_oe_pt.loc[:, log_cols] = s_b_oe_pt.loc[:, log_cols].apply(np.log10)


            s_s_pst = pyemu.Pst(os.path.join(s_s_m_d, "freyberg.pst"))
            seq_oe_files_pr = [f for f in os.listdir(s_s_m_d) if f.endswith("0.obs.csv") and f.startswith("freyberg")]
            spost_iter = s_s_pst.control_data.noptmax
            if post_iter is not None:
                spost_iter = post_iter
            seq_oe_files_pt = [f for f in os.listdir(s_s_m_d) if
                               f.endswith("{0}.obs.csv".format(spost_iter)) and f.startswith(
                                   "freyberg")]

            s_s_oe_dict_pr = {int(f.split(".")[1]): pd.read_csv(os.path.join(s_s_m_d, f), index_col=0) for f in
                              seq_oe_files_pr}
            s_s_oe_dict_pt = {int(f.split(".")[1]): pd.read_csv(os.path.join(s_s_m_d, f), index_col=0) for f in
                              seq_oe_files_pt}

            for key,df in s_s_oe_dict_pr.items():
                log_cols = df.columns.map(lambda x: "mass" in x or "cnc" in x)


                df.loc[:,log_cols] = df.loc[:,log_cols].apply(np.log10)
                df = df.replace(-np.Inf, np.nan)
                s_s_oe_dict_pr[key] = df

            for key, df in s_s_oe_dict_pt.items():
                log_cols = df.columns.map(lambda x: "mass" in x or "cnc" in x)
                df.loc[:, log_cols] = df.loc[:, log_cols].apply(np.log10)
                df = df.replace(-np.Inf, np.nan)
                s_s_oe_dict_pt[key] = df

            s_b_dict[ireal] = [s_b_pst,s_b_oe_pr,s_b_oe_pt]
            s_s_dict[ireal] = [s_s_pst,s_s_oe_dict_pr,s_s_oe_dict_pt]

            if include_est_states:
                # key these one cycle ahead since the posterior est states for this cycle are equiv to the prior sim states
                # of next cycle
                s_s_pe_dict_pt = {int(f.split(".")[1]) + 1: pd.read_csv(os.path.join(s_s_m_d, f.replace(".obs.",".par.")),
                                                                    index_col=0) for f in seq_oe_files_pt}
                s_s_est_dict[ireal] = s_s_pe_dict_pt
            print(ireal)
        except:
            break

    obs = s_s_pst.observation_data
    #sobs_to_sipar = obs.loc[pd.notna(obs.state_par_link),"state_par_link"].to_dict()
    #par = s_s_pst.parameter_data
    #sfpar = par.loc[pd.notna(par.state_par_link),:]
    #sipar_to_sfpar = {si:sf for si,sf in zip(sfpar.state_par_link,sfpar.parnme)}


    if len(s_b_dict) == 0:
        raise Exception()

    ireals = list(s_s_dict.keys())
    ireals.sort()
    sbobs_org = s_b_pst.observation_data
    print("plotting")
    size,lw=3,0.5
    pname = os.path.join(subdir,"s_vs_s.pdf")
    if post_iter is not None:
        pname = os.path.join(subdir,"s_vs_s_postiter_{0}.pdf".format(post_iter))

    is_1_lay = True
    if True in [True if "k:2" in o else False for o in s_b_pst.obs_names]:
        is_1_lay = False


    with PdfPages(pname) as pdf:
        for ogname in ognames:
            #if "mass" not in ogname:
            #    continue
            k0ogname = ogname
            if is_1_lay:
                k0ogname = ogname.replace("k:2","k:0")
            sgobs = sbobs_org.loc[sbobs_org.obsnme.str.contains(k0ogname),:].copy()
            sgobs = sgobs.loc[sgobs.obsnme.str.contains("_time"),:]
            sgobs.loc[:, "time"] = sgobs.time.apply(float)
            sgobs.sort_values(by="time", inplace=True)
            figall,axesall = plt.subplots(2, 2, figsize=(8, 8))
            for itime,oname in enumerate(sgobs.obsnme):
                if itime == 0:
                    continue
                #if itime != 3:
                #    continue
                fig, axes = plt.subplots(2, 2, figsize=(8, 6))
                for ireal in ireals:
                    s_b_pst,s_b_oe_pr,s_b_oe_pt = s_b_dict[ireal]
                    sbobs = s_b_pst.observation_data
                    sgobs = sbobs.loc[sbobs.obsnme.str.contains(k0ogname), :].copy()

                    s_s_pst,s_s_oe_dict_pr,s_s_oe_dict_pt = s_s_dict[ireal]

                    cval = sgobs.loc[oname,"obsval"].copy()
                    if "conc" in oname:
                        if np.abs(cval) > 1.0e+10:
                            continue
                    if "mass" in oname or "cnc" in oname:
                        cval = np.log10(cval)
                        if ~np.isfinite(cval):
                            continue
                    weight = sgobs.loc[oname,"weight"]
                    print(ireal,oname,cval)

                    if summarize:
                        mn = s_b_oe_pr.loc[:, oname].mean()
                        lq = s_b_oe_pr.loc[:, oname].quantile(0.05)
                        uq = s_b_oe_pr.loc[:, oname].quantile(0.95)
                        axes[0, 0].scatter(mn, cval,
                                           marker="o", color="0.5", alpha=0.5,s=size)
                        axes[0, 0].plot([lq,uq], [cval,cval],
                                           color="0.5", alpha=0.5,lw=lw)

                        axesall[0, 0].scatter(mn, cval,
                                           marker="o", color="0.5", alpha=0.5,s=size)
                        axesall[0, 0].plot([lq, uq], [cval, cval],
                                        color="0.5", alpha=0.5, lw=lw)

                        mn = s_b_oe_pt.loc[:, oname].mean()
                        lq = s_b_oe_pt.loc[:, oname].quantile(0.05)
                        uq = s_b_oe_pt.loc[:, oname].quantile(0.95)
                        axes[1, 0].scatter(mn, cval,
                                           marker="o", color="b", alpha=0.5,s=size)
                        axes[1, 0].plot([lq, uq], [cval, cval],
                                        color="b", alpha=0.5, lw=lw)

                        axesall[1, 0].scatter(mn, cval,
                                              marker="o", color="b", alpha=0.5,s=size,zorder=10)
                        axesall[1, 0].plot([lq, uq], [cval, cval],
                                           color="b", alpha=0.5, lw=lw,zorder=10)


                    else:
                        axes[0,0].scatter(s_b_oe_pr.loc[:, oname],[cval for _ in range(s_b_oe_pr.shape[0])],marker="o",color="0.5",alpha=0.5,s=size)
                        axes[1,0].scatter(s_b_oe_pt.loc[:, oname], [cval for _ in range(s_b_oe_pt.shape[0])], marker="o", color="b",
                                   alpha=0.5,s=size)
                        axesall[0,0].scatter(s_b_oe_pr.loc[:, oname], [cval for _ in range(s_b_oe_pr.shape[0])], marker="o",
                                   color="0.5", alpha=0.5,s=size)
                        axesall[1,0].scatter(s_b_oe_pt.loc[:, oname], [cval for _ in range(s_b_oe_pt.shape[0])], marker="o",
                                   color="b",alpha=0.5,s=size,zorder=10)

                    seq_name = k0ogname
                    if "arrobs" not in k0ogname:
                        seq_name = k0ogname + "_time:10000.0"

                    if itime in s_s_oe_dict_pr:
                        oe = s_s_oe_dict_pr[itime]
                        if summarize:
                            mn = oe.loc[:, seq_name].dropna().mean()
                            lq = oe.loc[:, seq_name].dropna().quantile(0.05)
                            uq = oe.loc[:, seq_name].dropna().quantile(0.95)
                            if itime == 3:
                                if ~np.isfinite(mn):
                                    print("seq",itime,seq_name,cval,mn,lq,uq)
                            axes[0, 1].scatter(mn, cval,
                                               marker="o", color="0.5", alpha=0.5,s=size)
                            axes[0, 1].plot([lq, uq], [cval, cval],
                                            color="0.5", alpha=0.5, lw=lw)

                            axesall[0, 1].scatter(mn, cval,
                                                  marker="o", color="0.5", alpha=0.5,s=size)
                            axesall[0, 1].plot([lq, uq], [cval, cval],
                                               color="0.5", alpha=0.5, lw=lw)

                        else:
                            axes[0,1].scatter(oe.loc[:, seq_name],[cval for _ in range(oe.shape[0])], marker="o", color="0.5",
                                       alpha=0.5,s=size)
                            axesall[0,1].scatter(oe.loc[:, seq_name], [cval for _ in range(oe.shape[0])], marker="o", color="0.5",
                                       alpha=0.5,s=size)

                    # the estimated states...
                    # s_s_pe_dict_pt = s_s_est_dict.get(ireal,None)
                    # if s_s_pe_dict_pt is not None and itime in s_s_pe_dict_pt:
                    #     pe = s_s_pe_dict_pt[itime]
                    #     pname = sobs_to_sipar[seq_name]
                    #     simseq_name = sipar_to_sfpar[pname]
                    #     if summarize:
                    #         mn = pe.loc[:, simseq_name].mean()
                    #         lq = pe.loc[:, simseq_name].quantile(0.05)
                    #         uq = pe.loc[:, simseq_name].quantile(0.95)
                    #         axes[0, 1].scatter(mn, cval,
                    #                            marker="*", color="0.5", alpha=0.5,s=size*2)
                    #         axes[0, 1].plot([lq, uq], [cval, cval],
                    #                         color="0.5", alpha=0.5, lw=lw,dashes=(1,1))
                    #
                    #         axesall[0, 1].scatter(mn, cval,
                    #                               marker="*", color="0.5", alpha=0.5,s=size*2)
                    #         axesall[0, 1].plot([lq, uq], [cval, cval],
                    #                            color="0.5", alpha=0.5, lw=lw,dashes=(1,1))
                    #
                    #     else:
                    #         axes[0,1].scatter(pe.loc[:, simseq_name],[cval for _ in range(pe.shape[0])], marker="*", color="0.5",
                    #                    alpha=0.5,s=size*2)
                    #         axesall[0,1].scatter(pe.loc[:, simseq_name], [cval for _ in range(pe.shape[0])], marker="*", color="0.5",
                    #                    alpha=0.5,s=size*2)

                    if itime in s_s_oe_dict_pt:
                        oe = s_s_oe_dict_pt[itime]
                        if summarize:
                            mn = oe.loc[:, seq_name].dropna().mean()
                            lq = oe.loc[:, seq_name].dropna().quantile(0.05)
                            uq = oe.loc[:, seq_name].dropna().quantile(0.95)
                            axes[1, 1].scatter(mn, cval,
                                               marker="o", color="b", alpha=0.5,s=size)
                            axes[1, 1].plot([lq, uq], [cval, cval],
                                            color="b", alpha=0.5, lw=lw)

                            axesall[1, 1].scatter(mn, cval,
                                                  marker="o", color="b", alpha=0.5,s=size)
                            axesall[1, 1].plot([lq, uq], [cval, cval],
                                               color="b", alpha=0.5, lw=lw)
                        else:

                            axes[1,1].scatter(oe.loc[:, seq_name],[cval for _ in range(oe.shape[0])], marker="o", color="b",
                                       alpha=0.5,s=size,zorder=10)
                            axesall[1,1].scatter(oe.loc[:, seq_name], [cval for _ in range(oe.shape[0])], marker="o", color="b",
                                       alpha=0.5,s=size,zorder=10)
                    elif itime in s_s_oe_dict_pr:
                        oe = s_s_oe_dict_pr[itime]
                        if summarize:
                            mn = oe.loc[:, seq_name].mean()
                            lq = oe.loc[:, seq_name].quantile(0.05)
                            uq = oe.loc[:, seq_name].quantile(0.95)
                            axes[1, 1].scatter(mn, cval,
                                               marker="o", color="0.5", alpha=0.5,s=size)
                            axes[1, 1].plot([lq, uq], [cval, cval],
                                            color="0.5", alpha=0.5, lw=lw)

                            axesall[1, 1].scatter(mn, cval,
                                                  marker="o", color="0.5", alpha=0.5,s=size)
                            axesall[1, 1].plot([lq, uq], [cval, cval],
                                               color="0.5", alpha=0.5, lw=lw)
                        else:

                            axes[1,1].scatter(oe.loc[:, seq_name],[cval for _ in range(oe.shape[0])], marker="o", color="0.5",
                                       alpha=0.5,s=size)
                            axesall[1,1].scatter(oe.loc[:, seq_name], [cval for _ in range(oe.shape[0])], marker="o", color="0.5",
                                       alpha=0.5,s=size)

                axes[0,0].set_title("A) batch prior",loc="left")
                axes[1, 0].set_title("C) batch posterior", loc="left")
                axes[0,1].set_title("B) sequential prior", loc="left")
                axes[1, 1].set_title("D) sequential posterior", loc="left")
                fig.suptitle(oname+", weight:{0}, cycle:{1}".format(weight,itime))
                mn = 1.0e+20
                mx = -1.0e+20
                for ax in axes.flatten():
                    if ax.get_xlim()[0] == 0:
                        continue
                    mn = min(mn,ax.get_ylim()[0],ax.get_xlim()[0])
                    mx = max(mx,ax.get_ylim()[1],ax.get_xlim()[1])
                for ax in axes.flatten():
                    ax.plot([mn,mx],[mn,mx],"k--")
                    ax.set_xlim(mn,mx)
                    ax.set_ylim(mn,mx)
                    ax.set_aspect("equal")
                    ax.set_ylabel("complex")
                    ax.set_xlabel("simple")
                plt.tight_layout()
                pdf.savefig(fig)
                plt.close(fig)

            axes = axesall
            axes[0, 0].set_title("A) batch prior", loc="left")
            axes[1, 0].set_title("C) batch posterior", loc="left")
            axes[0, 1].set_title("B) sequential prior", loc="left")
            axes[1, 1].set_title("D) sequential posterior", loc="left")
            figall.suptitle(ogname + " all times")
            mn = 1.0e+20
            mx = -1.0e+20
            for ax in axes.flatten():
                if ax.get_xlim()[0] == 0:
                    continue
                mn = min(mn, ax.get_ylim()[0], ax.get_xlim()[0])
                mx = max(mx, ax.get_ylim()[1], ax.get_xlim()[1])
            for ax in axes.flatten():
                ax.plot([mn, mx], [mn, mx], "k--")
                ax.set_xlim(mn, mx)
                ax.set_ylim(mn, mx)
                ax.set_aspect("equal")
                ax.set_ylabel("complex")
                ax.set_xlabel("simple")
            plt.tight_layout()
            pdf.savefig(figall)
            plt.close(figall)

def sync_phase(s_d = "monthly_model_files_org"):
    c_d = "daily_model_files_trnsprt_org"


    t_c_d = c_d.replace("_org","")
    t_s_d = s_d.replace("_org","")

    if os.path.exists(t_c_d):
        shutil.rmtree(t_c_d)
    shutil.copytree(c_d,t_c_d)
    if os.path.exists(t_s_d):
        shutil.rmtree(t_s_d)
    shutil.copytree(s_d, t_s_d)

    for d in [t_c_d,t_s_d]:
        rch_file_list = [f for f in os.listdir(d) if "rch_recharge" in f]
        rch_files = {}
        for f in rch_file_list:
            arr = np.loadtxt(os.path.join(d,f)) * 0.75
            np.savetxt(os.path.join(d,f),arr,fmt="%15.6E")
            rch_files[int(f.split(".")[1].split("_")[-1])] = arr
        keys = list(rch_files.keys())
        keys.sort()
        # calc the first sp recharge as the mean of the others
        new_first = np.zeros_like(rch_files[keys[0]])
        for key in keys[1:]:
            new_first += rch_files[key]
        new_first /= float(len(rch_files)-1)
        np.savetxt(os.path.join(d, "freyberg6.rch_recharge_1.txt"), new_first, fmt="%15.6E")
        rch_files[1] = new_first
        if "monthly" in d:
            s_rch_files = rch_files

        # process the wel flux files
        wel_files = [f for f in os.listdir(d) if ".wel_stress_period_data" in f]
        wel_files = {int(f.split(".")[1].split("_")[-1]):pd.read_csv(os.path.join(d,f),header=None,names=["l","r","c","flux"]) for f in wel_files}
        # no pumping in the first sp
        # now doing this in setup_pst
        #wel_files[1].loc[:,"flux"] = 0.0
        for sp,df in wel_files.items():
            df.loc[:, "aux"] = 0
            # uniform base pumping otherwise
            if sp != 1:
                df.loc[:,"flux"] = -300.0
            df.to_csv(os.path.join(d,"freyberg6.wel_stress_period_data_{0}.txt".format(sp)),index=False,header=False,sep=" ")



    c_sim = flopy.mf6.MFSimulation.load(sim_ws=t_c_d)
    cdts = pd.to_datetime("1-1-1900") + pd.to_timedelta(np.cumsum(c_sim.tdis.perioddata.array["perlen"]),unit="d")
    cdts_dict = {dt:i for i,dt in enumerate(cdts)}
    s_sim = flopy.mf6.MFSimulation.load(sim_ws=t_s_d)
    m = c_sim.get_model("freyberg6")

    nrow,ncol = m.dis.nrow.data,m.dis.ncol.data

    # get the mean simple model rch rate for each complex model stress period

    sdts = pd.to_datetime("1-1-1900") + pd.to_timedelta(np.cumsum(s_sim.tdis.perioddata.array["perlen"]),unit='d')
    dts,vals = [],[]
    for iend in range(1,len(sdts)):
        istart = iend - 1
        end = sdts[iend]
        start = sdts[istart]
        for dt,i in cdts_dict.items():
            if dt >=start and dt < end:
                dts.append(dt)
                vals.append(s_rch_files[iend].mean())
                #arr = np.zeros((nrow,ncol)) + s_rch_files[iend].mean()
                #np.savetxt(os.path.join(t_c_d, "freyberg6.rch_recharge_{0}.txt".format(i+1)), arr, fmt='%15.6E')

    #smooth that blocky shit!
    df = pd.DataFrame({"org":vals},index=dts)
    df.loc[:,"roll_month"] = df.rolling(60,center=True,min_periods=1).mean()
    # save the complex model daily rech arrays
    for i,val in enumerate(df.roll_month.values):
        arr = np.zeros((nrow, ncol)) + val
        np.savetxt(os.path.join(t_c_d,"freyberg6.rch_recharge_{0}.txt".format(i+1)),arr,fmt='%15.6E')

    # run em
    pyemu.os_utils.run("mf6", cwd=t_c_d)
    pyemu.os_utils.run("mf6", cwd=t_s_d)

    c_lst = flopy.utils.Mf6ListBudget(os.path.join(t_c_d,"freyberg6.lst"))
    s_lst = flopy.utils.Mf6ListBudget(os.path.join(t_s_d, "freyberg6.lst"))

    c_dfs = c_lst.get_dataframes(diff=True)
    s_dfs = s_lst.get_dataframes(diff=True)
    c_df = c_dfs[0]
    c_df_cum = c_dfs[1]
    s_df = s_dfs[0]
    s_df_cum = s_dfs[1]

    #s_df.index = s_df.index + pd.to_timedelta(15,unit="d")
    with PdfPages("lst_compare.pdf") as pdf:
        for col in c_df.columns:
            fig,ax = plt.subplots(1,1,figsize=(8,8))
            ax.plot(c_df.index.values,c_df.loc[:,col],color="c",label="complex")
            ax.plot(s_df.index.values, s_df.loc[:, col], color="m",label="simple")
            axt = plt.twinx(ax)
            axt.plot(c_df_cum.index.values,c_df_cum.loc[:,col],"c--",label="complex")
            axt.plot(s_df_cum.index.values, s_df_cum.loc[:, col],"m--",label="simple")
            ax.set_title(col)
            pdf.savefig()
            plt.close(fig)
    return t_c_d,t_s_d


def clean_results(subdir="."):
    clean_dir = os.path.join(subdir,"clean")
    if os.path.exists(clean_dir):
        shutil.rmtree(clean_dir)
    os.makedirs(clean_dir)
    #print(os.listdir(subdir))
    b_ds = [d for d in os.listdir(subdir) if os.path.isdir(os.path.join(subdir,d)) and d.startswith("monthly_model_files_master")]
    s_ds = [d for d in os.listdir(subdir) if os.path.isdir(os.path.join(subdir,d)) and d.startswith("seq_monthly_model_files_master")]

    b_tags = set(["obs_data.csv",".obs.csv",".par.csv",".rec"])
    def contains(fname):
        return True if True in [True if tag in fname else False for tag in b_tags] else False

    for b_d in s_ds:
        org_b_d = os.path.join(subdir,b_d)
        new_b_d = os.path.join(clean_dir,b_d)
        os.mkdir(new_b_d)
        files = [f.lower() for f in os.listdir(org_b_d)]
        keep = [f for f in files if contains(f)]
        for k in keep:
            shutil.copy2(os.path.join(org_b_d,k),os.path.join(new_b_d,k))
        print(b_d)

    for b_d in b_ds:
        org_b_d = os.path.join(subdir,b_d)
        new_b_d = os.path.join(clean_dir,b_d)
        os.mkdir(new_b_d)
        files = [f.lower() for f in os.listdir(org_b_d)]
        keep = [f for f in files if contains(f)]
        for k in keep:
            shutil.copy2(os.path.join(org_b_d,k),os.path.join(new_b_d,k))
        print(b_d)



def add_new_stress(m_d_org = "monthly_model_files"):

    m_lrc = (1,25,5)
    m_start_sp = 12
    d_start_sp = m_start_sp * 30
    d_lrc = (m_lrc[0],m_lrc[1]*3,m_lrc[2]*3)
    new_flux = -550
    d_d_org = "daily_model_files_trnsprt"
    d_d_new = d_d_org+"_newstress"
    if os.path.exists(d_d_new):
        shutil.rmtree(d_d_new)
    shutil.copytree(d_d_org,d_d_new)

    wel_files = [f for f in os.listdir(d_d_new) if ".wel_stress_period" in f and f.endswith(".txt")]
    for wel_file in wel_files:
        sp = int(wel_file.split(".")[1].split('_')[-1])
        df = pd.read_csv(os.path.join(d_d_new,wel_file),header=None,names=["l","r","c","flux","aux"],delim_whitespace=True)
        df.loc[6,["l","r","c"]] = [d_lrc[0],d_lrc[1],d_lrc[2]]
        if sp < d_start_sp:
            df.loc[6,"flux"] = 0.0
            df.loc[6, "aux"] = 0.0
        else:
            df.loc[6, "flux"] = new_flux
            df.loc[6, "aux"] = 0.0
        df = df.astype({"l": int, "r": int, "c": int})
        print(df.dtypes)
        df.to_csv(os.path.join(d_d_new,wel_file),header=False,index=False,sep=" ")
    wel_file = os.path.join(d_d_new,"freyberg6.wel")
    lines = open(wel_file,'r').readlines()
    with open(wel_file,'w') as f:
        for line in lines:
            if "maxbound" in line.lower():
                line = "   maxbound  7\n"
            f.write(line)
    pyemu.os_utils.run("mf6",cwd=d_d_new)


    # m_d_org = "monthly_model_files"

    m_d_new = m_d_org+"_newstress"
    if os.path.exists(m_d_new):
        shutil.rmtree(m_d_new)
    shutil.copytree(m_d_org,m_d_new)

    wel_files = [f for f in os.listdir(m_d_new) if ".wel_stress_period" in f and f.endswith(".txt")]
    for wel_file in wel_files:
        sp = int(wel_file.split(".")[1].split('_')[-1])
        df = pd.read_csv(os.path.join(m_d_new,wel_file),header=None,names=["l","r","c","flux","aux"],delim_whitespace=True)
        df.loc[6,["l","r","c"]] = [m_lrc[0],m_lrc[1],m_lrc[2]]
        if sp < m_start_sp:
            df.loc[6,"flux"] = 0.0
            df.loc[6, "aux"] = 0.0
        else:
            df.loc[6, "flux"] = new_flux
            df.loc[6, "aux"] = 0.0
        df = df.astype({"l": int, "r": int, "c": int})
        print(df.dtypes)
        df.to_csv(os.path.join(m_d_new,wel_file),header=False,index=False,sep=" ")
    wel_file = os.path.join(m_d_new,"freyberg6.wel")
    lines = open(wel_file,'r').readlines()
    with open(wel_file,'w') as f:
        for line in lines:
            if "maxbound" in line.lower():
                line = "   maxbound  7\n"
            f.write(line)
    pyemu.os_utils.run("mf6",cwd=m_d_new)


def make_muted_recharge(s_d = 'monthly_model_files_1lyr_org',c_d = 'daily_model_files_org'):
    s_d_new = s_d + '_muted_rch'
    if os.path.exists(s_d_new):
        shutil.rmtree(s_d_new)
    shutil.copytree(s_d, s_d_new)

    #make df of daily model recharge by daily timestep
    c_sim = flopy.mf6.MFSimulation.load(sim_ws = c_d)
    sdts_cmplx = pd.to_datetime("1-1-1900") + pd.to_timedelta(np.cumsum(c_sim.tdis.perioddata.array["perlen"]),unit='d')

    sim = flopy.mf6.MFSimulation.load(sim_ws = s_d)
    sdts_smpl = pd.to_datetime("1-1-1900") + pd.to_timedelta(np.cumsum(sim.tdis.perioddata.array["perlen"]),unit='d')

    m=0
    rch = pd.DataFrame(columns=['rch', 'simday'])
    rch_files = [f for f in os.listdir(c_d) if ".rch_recharge" in f and f.endswith(".txt")]
    for rch_file in rch_files:
        day = int(rch_file.split(".")[1].split('_')[-1])
        df = pd.read_csv(os.path.join(c_d,rch_file),header=None,delim_whitespace=True)
        rch.loc[m,'rch'] = df.iloc[0,0]
        rch.loc[m, 'simday'] = sdts_cmplx[day - 1]
        m+=1
    rch = rch.sort_values(by = 'simday').reset_index()

    rch_df = pd.DataFrame(columns = ['rch_avg'])
    for sp in range(len(sdts_smpl)-1):
        rch_avg = 0
        for day in range(1,len(rch)):
            if rch.simday[day] > sdts_smpl[sp] and rch.simday[day] <= sdts_smpl[sp + 1]:
                rch_avg += rch.rch[day]
        rch_df.loc[sp,'rch_avg'] = rch_avg

    rch_df = pd.concat([rch_df, rch_df])
    rch_df = rch_df.reset_index(drop='True')

    smooth_data = rch_df.rolling(window=9).mean()
    smooth_data = smooth_data.dropna().reset_index(drop = 'True')
    smooth_data = smooth_data[8:32].reset_index(drop = 'True')
    rch_df = rch_df[0:24]
    plt.plot(smooth_data, '*')
    plt.plot(rch_df)
    plt.show()

    rch_df = smooth_data
    rch_df.loc[-1] = rch_df.mean()
    rch_df.index = rch_df.index + 1  # shifting index
    rch_df.sort_index(inplace=True)

    rch_files = [f for f in os.listdir(s_d) if ".rch_recharge" in f and f.endswith(".txt")]
    for rch_file in rch_files:
        sp = int(rch_file.split(".")[1].split('_')[-1])
        df = pd.read_csv(os.path.join(s_d,rch_file),header=None,delim_whitespace=True)
        for row in range(40):
            for col in range(20):
                df[row,col] = rch_df.iloc[sp-1]
        df.to_csv(os.path.join(s_d_new, rch_file), header=False, index=False, sep=" ")


def reduce_simple_pars(t_d):
    pst = pyemu.Pst(os.path.join(t_d,"freyberg.pst"))
    pe = pyemu.ParameterEnsemble.from_binary(pst=pst,filename=os.path.join(t_d,"prior.jcb"))
    par = pst.parameter_data

    # # fix all but constant hks
    # hk_pars = par.loc[par.apply(lambda x: "npf_k" in x.parnme and "_cn" not in x.parnme and x.partrans != "fixed",axis=1),"parnme"]
    # par.loc[hk_pars,"partrans"] = "fixed"
    # pe.loc[:,hk_pars] = 1.0

    # fix all adjustable transient wel pars
    twel_par = par.loc[par.apply(lambda x: x.parnme.startswith("twel_mlt") and x.partrans != "fixed",axis=1),"parnme"]
    par.loc[twel_par,"partrans"] = "fixed"
    pe.loc[:,twel_par] = 1.0

    # fix all adjustable grid wel pars
    grwel_par = par.loc[par.parnme.str.startswith("wel_grid"),"parnme"]
    pe.loc[:,grwel_par] = 1.0
    par.loc[grwel_par,"partrans"] = "fixed"

    #grrch_par = par.loc[par.parnme.str.startswith("m_rch_gr"),"parnme"]
    #par.loc[grrch_par,"partrans"] = "fixed"
    #pe.loc[:,grrch_par] = 1.0

    #crch_par = par.loc[par.parnme.str.startswith("d_const_rch"),:].copy()
    #crch_par.loc[:,"sp"] = crch_par.parnme.apply(lambda x: int(x.split('_')[4]))
    #tie_crch_par = crch_par.parnme
    #val = crch_par.loc[tie_crch_par,"parval1"].mean()
    #pe.loc[:,tie_crch_par] = val
    #par.loc[tie_crch_par,"partrans"] = "fixed"
    #par.loc[tie_crch_par, "partied"] = "d_const_rch_recharge_2_cn_inst:0"

    pe.to_binary(os.path.join(t_d, "prior_reduced.jcb"))
    pst.pestpp_options["ies_par_en"] = "prior_reduced.jcb"
    pst.control_data.noptmax = -2
    pst.write(os.path.join(t_d,"test.pst"),version=2)
    pyemu.os_utils.run("pestpp-ies test.pst",cwd=t_d)

    pst.control_data.noptmax = -1
    pst.write(os.path.join(t_d,"freyberg.pst"),version=2)

def reduce_to_layer_pars(t_d):
    pst = pyemu.Pst(os.path.join(t_d,"freyberg.pst"))
    pe = pyemu.ParameterEnsemble.from_binary(pst=pst,filename=os.path.join(t_d,"prior.jcb"))
    par = pst.parameter_data
    # fix all grid and pilot point pars
    gr_par = par.loc[par.apply(lambda x: "gr" in x.parnme and x.partrans != "fixed",axis=1),"parnme"]
    par.loc[gr_par,"partrans"] = "fixed"
    pe.loc[:,gr_par] = 1.0

    pp_par = par.loc[par.apply(lambda x: "pp" in x.parnme and x.partrans != "fixed",axis=1),"parnme"]
    par.loc[pp_par,"partrans"] = "fixed"
    pe.loc[:,pp_par] = 1.0

    # fix all adjustable grid wel pars
    grwel_par = par.loc[par.parnme.str.startswith("wel_grid"),"parnme"]
    pe.loc[:,grwel_par] = 1.0
    par.loc[grwel_par,"partrans"] = "fixed"

    pe.to_binary(os.path.join(t_d, "prior_reduced.jcb"))
    pst.pestpp_options["ies_par_en"] = "prior_reduced.jcb"
    pst.control_data.noptmax = -2
    pst.write(os.path.join(t_d,"test.pst"),version=2)
    pyemu.os_utils.run("pestpp-ies test.pst",cwd=t_d)

    pst.control_data.noptmax = -1
    pst.write(os.path.join(t_d,"freyberg.pst"),version=2)

def make_prop_histograms(subdir="."):
    c_m_d = "daily_model_files_master_prior"
    c_pst = pyemu.Pst(os.path.join(c_m_d, "freyberg.pst"))
    cobs = c_pst.observation_data
    cobs = cobs.loc[cobs.obsnme.str.startswith("arrobs_hk"), :]
    c_oe = pd.read_csv(os.path.join(c_m_d, "freyberg.0.obs.csv"), index_col=0)
    c_hk = c_oe.loc[:,cobs.obsnme]
    c_hk_means = c_hk.mean(axis=1)
    c_hk_sds = c_hk.std(axis=1)

    bat_hk_means = []
    bat_hk_sds = []
    seq_hk_means = []
    seq_hk_sds = []

    for ireal in range(6):
        s_b_m_d = os.path.join(subdir, "monthly_model_files_master_{0}".format(ireal))
        s_s_m_d = os.path.join(subdir, "seq_monthly_model_files_master_{0}".format(ireal))

        # try:
        s_b_pst = pyemu.Pst(os.path.join(s_b_m_d, "freyberg.pst"))
        sobs = s_b_pst.observation_data
        sobs = sobs.loc[sobs.obsnme.str.startswith("arrobs_hk"), :]
        s_oe = pd.read_csv(os.path.join(s_b_m_d, "freyberg.3.obs.csv"), index_col=0)
        s_hk = s_oe.loc[:, sobs.obsnme]
        s_hk_mn = s_hk.stack().mean()
        s_hk_sd = s_hk.stack().std()
        bat_hk_means.append(s_hk_mn)
        bat_hk_sds.append(s_hk_sd)

        s_s_pst = pyemu.Pst(os.path.join(s_s_m_d, "freyberg.pst"))
        sobs = s_s_pst.observation_data
        sobs = sobs.loc[sobs.obsnme.str.startswith("arrobs_hk"), :]
        s_oe = pd.read_csv(os.path.join(s_s_m_d, "freyberg.24.0.obs.csv"), index_col=0)
        s_hk = s_oe.loc[:, sobs.obsnme]
        s_hk_mn = s_hk.stack().mean()
        s_hk_sd = s_hk.stack().std()
        seq_hk_means.append(s_hk_mn)
        seq_hk_sds.append(s_hk_sd)

    with PdfPages('hk_hists.pdf') as pdf:
        fig, axes = plt.subplots(2, 2, figsize=(8, 10))
        axes[0, 0].hist(c_hk_means, color='red',alpha =0.5)
        axes[0, 0].hist(seq_hk_means,color='blue',alpha =0.5)
        axes[1, 0].hist(c_hk_sds, color='red',alpha =0.5)
        axes[1, 0].hist(seq_hk_sds,color='blue',alpha =0.5)
        axes[0, 1].hist(c_hk_means, color='red',alpha =0.5)
        axes[0, 1].hist(bat_hk_means,color='blue',alpha =0.5)
        axes[1, 1].hist(c_hk_sds, color='red',alpha =0.5)
        axes[1, 1].hist(bat_hk_sds,color='blue',alpha =0.5)
        axes[0, 0].set_title("A) sequential HK first moment", loc="left")
        axes[1, 0].set_title("C) sequential HK second moment", loc="left")
        axes[0, 1].set_title("B) batch HK first moment", loc="left")
        axes[1, 1].set_title("D) batch HK second moment", loc="left")

        for ax in axes.flatten():
            ax.set_xlim(0, 100)
            ax.set_ylim(0, 20)
            ax.set_ylabel("frequency")
            ax.set_xlabel("Mean HK (ft/d)")
        pdf.savefig()
        plt.close(fig)


def plot_prior_mc_pub(subdir="."):
    """plot the prior monte carlo results for daily, monthly batch and monthly sequential

    """



    ognames = keep.copy()
    ognames = [ognames[0],ognames[-1]]
    #ognames.extend(forecast)
    ognames.sort()
    label_dict = keep_dict
    label_dict.update(forecast_dict)
    unit_dict = {n: u for n, u in zip(forecast, forecast_units)}
    unit_dict.update({n: u for n, u in zip(keep, keep_units)})
    print(unit_dict)


    c_m_d = os.path.join(subdir,"daily_model_files_master_prior")
    s_b_m_d = os.path.join(subdir,"monthly_model_files_master_prior")
    s_s_m_d = os.path.join(subdir,"seq_monthly_model_files_master_prior")

    top = np.loadtxt(os.path.join(s_b_m_d,"freyberg6.dis_top.txt"))


    c_pst = pyemu.Pst(os.path.join(c_m_d, "freyberg.pst"))
    obs = c_pst.observation_data
    cobs = obs.loc[obs.obsnme.str.startswith("hds_usecol:arrobs_head_"),:]
    cobs.loc[:,"time"] = cobs.time.apply(float)

    s_b_pst = pyemu.Pst(os.path.join(s_b_m_d, "freyberg.pst"))
    c_oe = pd.read_csv(os.path.join(c_m_d, "freyberg.0.obs.csv"), index_col=0)
    log = c_oe.columns.map(lambda x: "mass" in x)
    c_oe.loc[:,log] = c_oe.loc[:,log].apply(np.log10)
    s_b_oe = pd.read_csv(os.path.join(s_b_m_d, "freyberg.0.obs.csv"), index_col=0)
    log = s_b_oe.columns.map(lambda x: "mass" in x)
    s_b_oe.loc[:, log] = s_b_oe.loc[:, log].apply(np.log10)

    s_s_pst = pyemu.Pst(os.path.join(s_s_m_d,"freyberg.pst"))
    seq_oe_files = [f for f in os.listdir(s_s_m_d) if f.endswith(".oe.csv") and "global" in f and f.startswith("freyberg")]
    s_s_oe_dict = {int(f.split(".")[2]):pd.read_csv(os.path.join(s_s_m_d,f),index_col=0) for f in seq_oe_files}
    for key,df in s_s_oe_dict.items():
        log = df.columns.map(lambda x: "mass" in x)
        df.loc[:, log] = df.loc[:, log].apply(np.log10)
        s_s_oe_dict[key] = df

    pname = "prior_obs_v_sim_pub.pdf"
    if subdir != ".":
        pname = os.path.join(subdir,"prior_obs_v_sim_pub_"+subdir+".pdf")

    with PdfPages(pname) as pdf:
        fig, axes = plt.subplots(len(ognames), 1, figsize=(7, 5))
        for iax,ogname in enumerate(ognames):
            ax = axes[iax]
            cgobs = c_pst.observation_data.loc[c_pst.observation_data.obsnme.str.contains(ogname),:].copy()
            sgobs = s_b_pst.observation_data.loc[s_b_pst.observation_data.obsnme.str.contains(ogname),:].copy()
            k0ogname = ogname
            if sgobs.shape[0] == 0:
                k0ogname = ogname.replace("k:2","k:0")
                sgobs = s_b_pst.observation_data.loc[s_b_pst.observation_data.obsnme.str.contains(k0ogname), :].copy()
            sgobs.loc[:,"time"] = sgobs.time.apply(float)
            cgobs.loc[:, "time"] = cgobs.time.apply(float)

            sgobs.sort_values(by="time", inplace=True)
            cgobs.sort_values(by="time",inplace=True)

            seq_name = k0ogname
            if "arrobs" not in k0ogname:
                seq_name = k0ogname + "_time:10000.0"
            print(seq_name)
            for itime,time in enumerate(sgobs.time):
                if itime in s_s_oe_dict:
                    oe = s_s_oe_dict[itime]
                    #print(oe.loc[:,seq_name])
                    ax.scatter([time for _ in range(oe.shape[0])], oe.loc[:, seq_name], marker=".",color="0.5",alpha=0.5)

            [ax.plot(cgobs.time, c_oe.loc[idx, cgobs.obsnme], "r", lw=0.01, alpha=0.5) for idx in c_oe.index]

            [ax.plot(sgobs.time, s_b_oe.loc[idx, sgobs.obsnme], "0.5", lw=0.01, alpha=0.5) for idx in s_b_oe.index]


            # if "arrobs" in ogname:
            #     i = sgobs.i.apply(int)[0]
            #     j = sgobs.j.apply(int)[0]
            #     t = top[i, j]
            #     ax.plot(ax.get_xlim(),[t,t],"k--",lw=3)
            ax.set_title("{0}) {1}".format(string.ascii_uppercase[iax], label_dict[ogname]),loc="left")
            if "arrobs" in ogname:
                ax.set_ylim(30,ax.get_ylim()[1])
            ax.set_ylabel(unit_dict[ogname])

            #if "gage" not in ogname:
            #    ax.set_ylim(30,ax.get_ylim()[1])
        #for ax in axes[:-1]:
        #    ax.set_xticklabels([])
            ax.set_xlabel("time ($days$)")
        plt.tight_layout()
        pdf.savefig()
        plt.close(fig)

def plot_s_vs_s_pub(summarize=False, subdir=".", post_iter=None):
    include_est_states = False
    ognames = keep.copy()
    ognames.extend(forecast)
    label_dict = keep_dict
    label_dict.update(forecast_dict)

    # first rip thru all the dirs and load...
    s_b_dict = {}
    s_s_dict = {}
    s_s_est_dict = {}
    print("loading results...")

    for ireal in range(50):
        s_b_m_d = os.path.join(subdir,"monthly_model_files_master_{0}".format(ireal))
        s_s_m_d = os.path.join(subdir,"seq_monthly_model_files_master_{0}".format(ireal))

        if not os.path.exists(s_s_m_d) or not os.path.exists(s_s_m_d):
            break
        try:
            s_b_pst = pyemu.Pst(os.path.join(s_b_m_d, "freyberg.pst"))
            obs = s_b_pst.observation_data.loc[s_b_pst.nnz_obs_names,:]
            obs = obs.loc[obs.obsnme.str.contains("sfr_usecol:gage_1")]
            if obs.obsval.max() > 8000:
                continue

            s_b_oe_pr = pd.read_csv(os.path.join(s_b_m_d, "freyberg.0.obs.csv"), index_col=0)
            log_cols = s_b_oe_pr.columns.map(lambda x: "mass" in x or "cnc" in x)
            s_b_oe_pr.loc[:,log_cols] = s_b_oe_pr.loc[:,log_cols].apply(np.log10)

            bpost_iter = s_b_pst.control_data.noptmax
            if post_iter is not None:
                bpost_iter = post_iter
            s_b_oe_pt = pd.read_csv(os.path.join(s_b_m_d, "freyberg.{0}.obs.csv".format(bpost_iter)),
                                    index_col=0)
            log_cols = s_b_oe_pt.columns.map(lambda x: "mass" in x or "cnc" in x)
            s_b_oe_pt.loc[:, log_cols] = s_b_oe_pt.loc[:, log_cols].apply(np.log10)


            s_s_pst = pyemu.Pst(os.path.join(s_s_m_d, "freyberg.pst"))
            seq_oe_files_pr = [f for f in os.listdir(s_s_m_d) if f.endswith("0.obs.csv") and f.startswith("freyberg")]
            spost_iter = s_s_pst.control_data.noptmax
            if post_iter is not None:
                spost_iter = post_iter
            seq_oe_files_pt = [f for f in os.listdir(s_s_m_d) if
                               f.endswith("{0}.obs.csv".format(spost_iter)) and f.startswith(
                                   "freyberg")]

            s_s_oe_dict_pr = {int(f.split(".")[1]): pd.read_csv(os.path.join(s_s_m_d, f), index_col=0) for f in
                              seq_oe_files_pr}
            s_s_oe_dict_pt = {int(f.split(".")[1]): pd.read_csv(os.path.join(s_s_m_d, f), index_col=0) for f in
                              seq_oe_files_pt}

            for key,df in s_s_oe_dict_pr.items():
                log_cols = df.columns.map(lambda x: "mass" in x or "cnc" in x)


                df.loc[:,log_cols] = df.loc[:,log_cols].apply(np.log10)
                df = df.replace(-np.Inf, np.nan)
                s_s_oe_dict_pr[key] = df

            for key, df in s_s_oe_dict_pt.items():
                log_cols = df.columns.map(lambda x: "mass" in x or "cnc" in x)
                df.loc[:, log_cols] = df.loc[:, log_cols].apply(np.log10)
                df = df.replace(-np.Inf, np.nan)
                s_s_oe_dict_pt[key] = df

            s_b_dict[ireal] = [s_b_pst,s_b_oe_pr,s_b_oe_pt]
            s_s_dict[ireal] = [s_s_pst,s_s_oe_dict_pr,s_s_oe_dict_pt]

            if include_est_states:
                # key these one cycle ahead since the posterior est states for this cycle are equiv to the prior sim states
                # of next cycle
                s_s_pe_dict_pt = {int(f.split(".")[1]) + 1: pd.read_csv(os.path.join(s_s_m_d, f.replace(".obs.",".par.")),
                                                                    index_col=0) for f in seq_oe_files_pt}
                s_s_est_dict[ireal] = s_s_pe_dict_pt
            print(ireal)
        except:
            break

    obs = s_s_pst.observation_data
    #sobs_to_sipar = obs.loc[pd.notna(obs.state_par_link),"state_par_link"].to_dict()
    #par = s_s_pst.parameter_data
    #sfpar = par.loc[pd.notna(par.state_par_link),:]
    #sipar_to_sfpar = {si:sf for si,sf in zip(sfpar.state_par_link,sfpar.parnme)}


    if len(s_b_dict) == 0:
        raise Exception()

    ireals = list(s_s_dict.keys())
    ireals.sort()
    sbobs_org = s_b_pst.observation_data
    print("plotting")
    size,lw=3,0.5
    pname = os.path.join(subdir,"s_vs_s_pub.pdf")
    if subdir != ".":
        pname = pname.replace(".pdf","_"+subdir+".pdf")
    #if post_iter is not None:
    #    pname = os.path.join(subdir,"s_vs_s_postiter_{0}.pdf".format(post_iter))

    is_1_lay = True
    if True in [True if "k:2" in o else False for o in s_b_pst.obs_names]:
        is_1_lay = False

    with PdfPages(pname) as pdf:
        figall_keep, axesall_keep = plt.subplots(len(keep), 2, figsize=(6, 8))
        ax_count = 0
        for ikeep,ogname in enumerate(keep):
            #if "mass" not in ogname:
            #    continue
            k0ogname = ogname
            if is_1_lay:
                k0ogname = ogname.replace("k:2","k:0")
            sgobs = sbobs_org.loc[sbobs_org.obsnme.str.contains(k0ogname),:].copy()
            sgobs = sgobs.loc[sgobs.obsnme.str.contains("_time"),:]
            sgobs.loc[:, "time"] = sgobs.time.apply(float)
            #sgobs = sgobs.loc[sgobs.time.apply(lambda x: x > 10000 and x < 10366),:]
            sgobs.sort_values(by="time", inplace=True)
            axesall_keep[ikeep,0].set_title("{0}) batch {1}".format(string.ascii_uppercase[ax_count],
                                                                    label_dict[ogname]),loc="left")
            ax_count += 1
            axesall_keep[ikeep, 1].set_title(
                "{0}) sequential {1}".format(string.ascii_uppercase[ax_count], label_dict[ogname]),loc="left")
            ax_count += 1
            for itime,oname in enumerate(sgobs.obsnme):
                if itime < 1 or itime > 12:
                #if itime != 12:
                    continue

                for ireal in ireals:
                    s_b_pst,s_b_oe_pr,s_b_oe_pt = s_b_dict[ireal]
                    sbobs = s_b_pst.observation_data
                    sgobs = sbobs.loc[sbobs.obsnme.str.contains(k0ogname), :].copy()

                    s_s_pst,s_s_oe_dict_pr,s_s_oe_dict_pt = s_s_dict[ireal]

                    cval = sgobs.loc[oname,"obsval"].copy()
                    if "conc" in oname:
                        if np.abs(cval) > 1.0e+10:
                            continue
                    if "mass" in oname or "cnc" in oname:
                        cval = np.log10(cval)
                        if ~np.isfinite(cval):
                            continue
                    weight = sgobs.loc[oname,"weight"]
                    print(ireal,oname,cval)

                    if summarize:
                        mn = s_b_oe_pr.loc[:, oname].mean()
                        lq = s_b_oe_pr.loc[:, oname].quantile(0.05)
                        uq = s_b_oe_pr.loc[:, oname].quantile(0.95)
                        #axes[0, 0].scatter(mn, cval,
                        #                   marker="o", color="0.5", alpha=0.5,s=size)
                        #axes[0, 0].plot([lq,uq], [cval,cval],
                        #                   color="0.5", alpha=0.5,lw=lw)

                        #axesall_keep[ikeep, 0].scatter(mn, cval,
                        #                   marker="o", color="0.5", alpha=0.5,s=size)
                        #axesall_keep[ikeep, 0].plot([lq, uq], [cval, cval],
                        #                color="0.5", alpha=0.5, lw=lw)

                        mn = s_b_oe_pt.loc[:, oname].mean()
                        lq = s_b_oe_pt.loc[:, oname].quantile(0.05)
                        uq = s_b_oe_pt.loc[:, oname].quantile(0.95)
                        #axes[1, 0].scatter(mn, cval,
                        #                   marker="o", color="b", alpha=0.5,s=size)
                        #axes[1, 0].plot([lq, uq], [cval, cval],
                        #                color="b", alpha=0.5, lw=lw)

                        axesall_keep[ikeep, 0].scatter(mn, cval,
                                              marker="o", color="b", alpha=0.5,s=size,zorder=10)
                        axesall_keep[ikeep, 0].plot([lq, uq], [cval, cval],
                                           color="b", alpha=0.5, lw=lw,zorder=10)


                    else:
                        #axes[0,0].scatter(s_b_oe_pr.loc[:, oname],[cval for _ in range(s_b_oe_pr.shape[0])],marker="o",color="0.5",alpha=0.5,s=size)
                        #axes[1,0].scatter(s_b_oe_pt.loc[:, oname], [cval for _ in range(s_b_oe_pt.shape[0])], marker="o", color="b",
                        #           alpha=0.5,s=size)
                        axesall_keep[ikeep,0].scatter(s_b_oe_pr.loc[:, oname], [cval for _ in range(s_b_oe_pr.shape[0])], marker="o",
                                   color="0.5", alpha=0.5,s=size)
                        axesall_keep[ikeep,0].scatter(s_b_oe_pt.loc[:, oname], [cval for _ in range(s_b_oe_pt.shape[0])], marker="o",
                                   color="b",alpha=0.5,s=size,zorder=10)

                    seq_name = k0ogname
                    if "arrobs" not in k0ogname:
                        seq_name = k0ogname + "_time:10000.0"

                    if itime in s_s_oe_dict_pr:
                        oe = s_s_oe_dict_pr[itime]
                        if summarize:
                            mn = oe.loc[:, seq_name].dropna().mean()
                            lq = oe.loc[:, seq_name].dropna().quantile(0.05)
                            uq = oe.loc[:, seq_name].dropna().quantile(0.95)
                            if itime == 3:
                                if ~np.isfinite(mn):
                                    print("seq",itime,seq_name,cval,mn,lq,uq)
                            #axes[0, 1].scatter(mn, cval,
                            #                   marker="o", color="0.5", alpha=0.5,s=size)
                            #axes[0, 1].plot([lq, uq], [cval, cval],
                            #                color="0.5", alpha=0.5, lw=lw)

                            #axesall_keep[ikeep, 1].scatter(mn, cval,
                            #                      marker="o", color="0.5", alpha=0.5,s=size)
                            #axesall_keep[ikeep, 1].plot([lq, uq], [cval, cval],
                            #                   color="0.5", alpha=0.5, lw=lw)

                        else:
                            #axes[0,1].scatter(oe.loc[:, seq_name],[cval for _ in range(oe.shape[0])], marker="o", color="0.5",
                            #           alpha=0.5,s=size)
                            axesall_keep[ikeep,1].scatter(oe.loc[:, seq_name], [cval for _ in range(oe.shape[0])], marker="o", color="0.5",
                                       alpha=0.5,s=size)


                    if itime in s_s_oe_dict_pt:
                        oe = s_s_oe_dict_pt[itime]
                        if summarize:
                            mn = oe.loc[:, seq_name].dropna().mean()
                            lq = oe.loc[:, seq_name].dropna().quantile(0.05)
                            uq = oe.loc[:, seq_name].dropna().quantile(0.95)
                            #axes[1, 1].scatter(mn, cval,
                            #                   marker="o", color="b", alpha=0.5,s=size)
                            #axes[1, 1].plot([lq, uq], [cval, cval],
                            #                color="b", alpha=0.5, lw=lw)

                            axesall_keep[ikeep, 1].scatter(mn, cval,
                                                  marker="o", color="b", alpha=0.5,s=size)
                            axesall_keep[ikeep, 1].plot([lq, uq], [cval, cval],
                                               color="b", alpha=0.5, lw=lw)
                        else:

                            axesall_keep[ikeep,1].scatter(oe.loc[:, seq_name],[cval for _ in range(oe.shape[0])], marker="o", color="b",
                                       alpha=0.5,s=size,zorder=10)
                            axesall_keep[ikeep,1].scatter(oe.loc[:, seq_name], [cval for _ in range(oe.shape[0])], marker="o", color="b",
                                       alpha=0.5,s=size,zorder=10)
                    # elif itime in s_s_oe_dict_pr:
                    #     oe = s_s_oe_dict_pr[itime]
                    #     if summarize:
                    #         mn = oe.loc[:, seq_name].mean()
                    #         lq = oe.loc[:, seq_name].quantile(0.05)
                    #         uq = oe.loc[:, seq_name].quantile(0.95)
                    #         #axes[1, 1].scatter(mn, cval,
                    #         #                   marker="o", color="0.5", alpha=0.5,s=size)
                    #         #axes[1, 1].plot([lq, uq], [cval, cval],
                    #         #                color="0.5", alpha=0.5, lw=lw)
                    #
                    #         axesall_keep[ikeep, 1].scatter(mn, cval,
                    #                               marker="o", color="0.5", alpha=0.5,s=size)
                    #         axesall_keep[ikeep, 1].plot([lq, uq], [cval, cval],
                    #                            color="0.5", alpha=0.5, lw=lw)
                    #     else:
                    #
                    #         #axes[1,1].scatter(oe.loc[:, seq_name],[cval for _ in range(oe.shape[0])], marker="o", color="0.5",
                    #         #          alpha=0.5,s=size)
                    #         axesall_keep[ikeep,1].scatter(oe.loc[:, seq_name], [cval for _ in range(oe.shape[0])], marker="o", color="0.5",
                    #                    alpha=0.5,s=size)


            mn = min(axesall_keep[ikeep,0].get_xlim()[0],axesall_keep[ikeep,0].get_ylim()[0],
                     axesall_keep[ikeep,1].get_xlim()[0],axesall_keep[ikeep,1].get_ylim()[0])
            mx = max(axesall_keep[ikeep, 0].get_xlim()[1], axesall_keep[ikeep, 0].get_ylim()[1],
                     axesall_keep[ikeep, 1].get_xlim()[1], axesall_keep[ikeep, 1].get_ylim()[1])
            axesall_keep[ikeep,0].plot([mn,mx],[mn,mx],"k--")
            axesall_keep[ikeep, 1].plot([mn, mx], [mn, mx], "k--")

            axesall_keep[ikeep, 0].set_xlim(mn, mx)
            axesall_keep[ikeep, 0].set_ylim(mn, mx)
            axesall_keep[ikeep, 1].set_xlim(mn, mx)
            axesall_keep[ikeep, 1].set_ylim(mn, mx)
            axesall_keep[ikeep,0].set_aspect("equal")
            axesall_keep[ikeep,0].set_ylabel("complex")
            axesall_keep[ikeep,0].set_xlabel("simple")
            axesall_keep[ikeep, 1].set_aspect("equal")
            axesall_keep[ikeep, 1].set_ylabel("complex")
            axesall_keep[ikeep, 1].set_xlabel("simple")


        plt.tight_layout()
        pdf.savefig(figall_keep)
        plt.close(figall_keep)

        figall_keep, axesall_keep = plt.subplots(len(forecast), 2, figsize=(5, 8))
        ax_count = 0
        for ikeep, ogname in enumerate(forecast):
            # if "mass" not in ogname:
            #    continue
            k0ogname = ogname
            if is_1_lay:
                k0ogname = ogname.replace("k:2", "k:0")
            sgobs = sbobs_org.loc[sbobs_org.obsnme.str.contains(k0ogname), :].copy()
            sgobs = sgobs.loc[sgobs.obsnme.str.contains("_time"), :]
            sgobs.loc[:, "time"] = sgobs.time.apply(float)
            # sgobs = sgobs.loc[sgobs.time.apply(lambda x: x > 10000 and x < 10366),:]
            sgobs.sort_values(by="time", inplace=True)
            axesall_keep[ikeep, 0].set_title(
                "{0}) batch {1}".format(string.ascii_uppercase[ax_count], label_dict[ogname]),loc="left")
            ax_count += 1
            axesall_keep[ikeep, 1].set_title(
                "{0}) sequential {1}".format(string.ascii_uppercase[ax_count], label_dict[ogname]),loc="left")
            ax_count += 1
            for itime, oname in enumerate(sgobs.obsnme):
                if itime < 12 or itime > 16:
                    continue

                for ireal in ireals:
                    s_b_pst, s_b_oe_pr, s_b_oe_pt = s_b_dict[ireal]
                    sbobs = s_b_pst.observation_data
                    sgobs = sbobs.loc[sbobs.obsnme.str.contains(k0ogname), :].copy()

                    s_s_pst, s_s_oe_dict_pr, s_s_oe_dict_pt = s_s_dict[ireal]

                    cval = sgobs.loc[oname, "obsval"].copy()
                    if "conc" in oname:
                        if np.abs(cval) > 1.0e+10:
                            continue
                    if "mass" in oname or "cnc" in oname:
                        cval = np.log10(cval)
                        if ~np.isfinite(cval):
                            continue
                    weight = sgobs.loc[oname, "weight"]
                    print(ireal, oname, cval)

                    if summarize:
                        mn = s_b_oe_pr.loc[:, oname].mean()
                        lq = s_b_oe_pr.loc[:, oname].quantile(0.05)
                        uq = s_b_oe_pr.loc[:, oname].quantile(0.95)
                        # axes[0, 0].scatter(mn, cval,
                        #                   marker="o", color="0.5", alpha=0.5,s=size)
                        # axes[0, 0].plot([lq,uq], [cval,cval],
                        #                   color="0.5", alpha=0.5,lw=lw)

                        # axesall_keep[ikeep, 0].scatter(mn, cval,
                        #                   marker="o", color="0.5", alpha=0.5,s=size)
                        # axesall_keep[ikeep, 0].plot([lq, uq], [cval, cval],
                        #                color="0.5", alpha=0.5, lw=lw)

                        mn = s_b_oe_pt.loc[:, oname].mean()
                        lq = s_b_oe_pt.loc[:, oname].quantile(0.05)
                        uq = s_b_oe_pt.loc[:, oname].quantile(0.95)
                        # axes[1, 0].scatter(mn, cval,
                        #                   marker="o", color="b", alpha=0.5,s=size)
                        # axes[1, 0].plot([lq, uq], [cval, cval],
                        #                color="b", alpha=0.5, lw=lw)

                        axesall_keep[ikeep, 0].scatter(mn, cval,
                                                       marker="o", color="b", alpha=0.5, s=size, zorder=10)
                        axesall_keep[ikeep, 0].plot([lq, uq], [cval, cval],
                                                    color="b", alpha=0.5, lw=lw, zorder=10)


                    else:
                        # axes[0,0].scatter(s_b_oe_pr.loc[:, oname],[cval for _ in range(s_b_oe_pr.shape[0])],marker="o",color="0.5",alpha=0.5,s=size)
                        # axes[1,0].scatter(s_b_oe_pt.loc[:, oname], [cval for _ in range(s_b_oe_pt.shape[0])], marker="o", color="b",
                        #           alpha=0.5,s=size)
                        axesall_keep[ikeep, 0].scatter(s_b_oe_pr.loc[:, oname],
                                                       [cval for _ in range(s_b_oe_pr.shape[0])], marker="o",
                                                       color="0.5", alpha=0.5, s=size)
                        axesall_keep[ikeep, 0].scatter(s_b_oe_pt.loc[:, oname],
                                                       [cval for _ in range(s_b_oe_pt.shape[0])], marker="o",
                                                       color="b", alpha=0.5, s=size, zorder=10)

                    seq_name = k0ogname
                    if "arrobs" not in k0ogname:
                        seq_name = k0ogname + "_time:10000.0"

                    if itime in s_s_oe_dict_pr:
                        oe = s_s_oe_dict_pr[itime]
                        if summarize:
                            mn = oe.loc[:, seq_name].dropna().mean()
                            lq = oe.loc[:, seq_name].dropna().quantile(0.05)
                            uq = oe.loc[:, seq_name].dropna().quantile(0.95)
                            if itime == 3:
                                if ~np.isfinite(mn):
                                    print("seq", itime, seq_name, cval, mn, lq, uq)

                            axesall_keep[ikeep, 1].scatter(mn, cval,
                                                  marker="o", color="0.5", alpha=0.5,s=size)
                            axesall_keep[ikeep, 1].plot([lq, uq], [cval, cval],
                                               color="0.5", alpha=0.5, lw=lw)

                        else:
                            # axes[0,1].scatter(oe.loc[:, seq_name],[cval for _ in range(oe.shape[0])], marker="o", color="0.5",
                            #           alpha=0.5,s=size)
                            axesall_keep[ikeep, 1].scatter(oe.loc[:, seq_name], [cval for _ in range(oe.shape[0])],
                                                           marker="o", color="0.5",
                                                           alpha=0.5, s=size)

                    if itime in s_s_oe_dict_pt:
                        oe = s_s_oe_dict_pt[itime]
                        if summarize:
                            mn = oe.loc[:, seq_name].dropna().mean()
                            lq = oe.loc[:, seq_name].dropna().quantile(0.05)
                            uq = oe.loc[:, seq_name].dropna().quantile(0.95)
                            # axes[1, 1].scatter(mn, cval,
                            #                   marker="o", color="b", alpha=0.5,s=size)
                            # axes[1, 1].plot([lq, uq], [cval, cval],
                            #                color="b", alpha=0.5, lw=lw)

                            #axesall_keep[ikeep, 1].scatter(mn, cval,
                            #                               marker="o", color="b", alpha=0.5, s=size)
                            #axesall_keep[ikeep, 1].plot([lq, uq], [cval, cval],
                            #                            color="b", alpha=0.5, lw=lw)
                        else:

                            axesall_keep[ikeep, 1].scatter(oe.loc[:, seq_name], [cval for _ in range(oe.shape[0])],
                                                        marker="o", color="b",
                                                        alpha=0.5, s=size, zorder=10)
                            axesall_keep[ikeep, 1].scatter(oe.loc[:, seq_name], [cval for _ in range(oe.shape[0])],
                                                           marker="o", color="b",
                                                           alpha=0.5, s=size, zorder=10)
                    # elif itime in s_s_oe_dict_pr:
                    #     oe = s_s_oe_dict_pr[itime]
                    #     if summarize:
                    #         mn = oe.loc[:, seq_name].mean()
                    #         lq = oe.loc[:, seq_name].quantile(0.05)
                    #         uq = oe.loc[:, seq_name].quantile(0.95)
                    #         # axes[1, 1].scatter(mn, cval,
                    #         #                   marker="o", color="0.5", alpha=0.5,s=size)
                    #         # axes[1, 1].plot([lq, uq], [cval, cval],
                    #         #                color="0.5", alpha=0.5, lw=lw)
                    #
                    #         axesall_keep[ikeep, 1].scatter(mn, cval,
                    #                                        marker="o", color="0.5", alpha=0.5, s=size)
                    #         axesall_keep[ikeep, 1].plot([lq, uq], [cval, cval],
                    #                                     color="0.5", alpha=0.5, lw=lw)
                    #     else:
                    #
                    #         # axes[1,1].scatter(oe.loc[:, seq_name],[cval for _ in range(oe.shape[0])], marker="o", color="0.5",
                    #         #          alpha=0.5,s=size)
                    #         axesall_keep[ikeep, 1].scatter(oe.loc[:, seq_name], [cval for _ in range(oe.shape[0])],
                    #                                        marker="o", color="0.5",
                    #                                        alpha=0.5, s=size)

            mn = min(axesall_keep[ikeep, 0].get_xlim()[0], axesall_keep[ikeep, 0].get_ylim()[0],
                     axesall_keep[ikeep, 1].get_xlim()[0], axesall_keep[ikeep, 1].get_ylim()[0])
            mx = max(axesall_keep[ikeep, 0].get_xlim()[1], axesall_keep[ikeep, 0].get_ylim()[1],
                     axesall_keep[ikeep, 1].get_xlim()[1], axesall_keep[ikeep, 1].get_ylim()[1])
            axesall_keep[ikeep, 0].plot([mn, mx], [mn, mx], "k--")
            axesall_keep[ikeep, 1].plot([mn, mx], [mn, mx], "k--")
            #print(ogname,mn,mx)
            axesall_keep[ikeep, 0].set_xlim(mn, mx)
            axesall_keep[ikeep, 0].set_ylim(mn, mx)
            axesall_keep[ikeep, 1].set_xlim(mn, mx)
            axesall_keep[ikeep, 1].set_ylim(mn, mx)
            axesall_keep[ikeep, 0].set_aspect("equal")
            axesall_keep[ikeep, 0].set_ylabel("complex")
            axesall_keep[ikeep, 0].set_xlabel("simple")
            axesall_keep[ikeep, 1].set_aspect("equal")
            axesall_keep[ikeep, 1].set_ylabel("complex")
            axesall_keep[ikeep, 1].set_xlabel("simple")

        plt.tight_layout()
        pdf.savefig(figall_keep)
        plt.close(figall_keep)

        cycles = np.arange(12,16,dtype=int)
        norm = mpl.colors.Normalize(vmin=0,vmax=len(cycles))
        cmap = plt.get_cmap("jet")
        figall_keep, axesall_keep = plt.subplots(1, 2, figsize=(8, 4))
        ax_count = 0
        for ikeep, ogname in enumerate(forecast[1:]):
            # if "mass" not in ogname:
            #    continue
            k0ogname = ogname
            if is_1_lay:
                k0ogname = ogname.replace("k:2", "k:0")
            sgobs = sbobs_org.loc[sbobs_org.obsnme.str.contains(k0ogname), :].copy()
            sgobs = sgobs.loc[sgobs.obsnme.str.contains("_time"), :]
            sgobs.loc[:, "time"] = sgobs.time.apply(float)
            # sgobs = sgobs.loc[sgobs.time.apply(lambda x: x > 10000 and x < 10366),:]
            sgobs.sort_values(by="time", inplace=True)
            axesall_keep[ikeep].set_title(
                "{0}) {1}".format(string.ascii_uppercase[ax_count], label_dict[ogname].replace("\n"," ")),loc="left")
            ax_count += 1
            axesall_keep[ikeep].set_title(
                "{0}) {1}".format(string.ascii_uppercase[ax_count], label_dict[ogname].replace("\n"," ")),loc="left")
            ax_count += 1
            iitime = 0
            for itime, oname in enumerate(sgobs.obsnme):
                if itime not in cycles:
                    continue
                color = cmap(norm(iitime))
                iitime += 1
                print(itime,oname)
                for ireal in ireals:
                    s_b_pst, s_b_oe_pr, s_b_oe_pt = s_b_dict[ireal]
                    sbobs = s_b_pst.observation_data
                    sgobs = sbobs.loc[sbobs.obsnme.str.contains(k0ogname), :].copy()

                    s_s_pst, s_s_oe_dict_pr, s_s_oe_dict_pt = s_s_dict[ireal]

                    cval = sgobs.loc[oname, "obsval"].copy()
                    if "conc" in oname:
                        if np.abs(cval) > 1.0e+10:
                            continue
                    if "mass" in oname or "cnc" in oname:
                        cval = np.log10(cval)
                        if ~np.isfinite(cval):
                            continue

                    seq_name = k0ogname
                    if "arrobs" not in k0ogname:
                        seq_name = k0ogname + "_time:10000.0"

                    if itime in s_s_oe_dict_pr:
                        oe = s_s_oe_dict_pr[itime]
                        print(itime,ogname,color)
                        if summarize:
                            mn = oe.loc[:, seq_name].dropna().mean()
                            lq = oe.loc[:, seq_name].dropna().quantile(0.05)
                            uq = oe.loc[:, seq_name].dropna().quantile(0.95)

                            axesall_keep[ikeep].scatter(mn, cval,
                                                  marker="o", color=color, alpha=0.5,s=size)
                            axesall_keep[ikeep].plot([lq, uq], [cval, cval],
                                               color=color, alpha=0.5, lw=lw)

                        #else:
                            # axes[0,1].scatter(oe.loc[:, seq_name],[cval for _ in range(oe.shape[0])], marker="o", color="0.5",
                            #           alpha=0.5,s=size)
                            #axesall_keep[ikeep, 1].scatter(oe.loc[:, seq_name], [cval for _ in range(oe.shape[0])],
                            #                               marker="o", color=color,
                            #                               alpha=0.5, s=size)

                    # if int(itime) in s_s_oe_dict_pt:
                    #     oe = s_s_oe_dict_pt[itime]
                    #     print(ogname,itime, color)
                    #     if summarize:
                    #         mn = oe.loc[:, seq_name].dropna().mean()
                    #         lq = oe.loc[:, seq_name].dropna().quantile(0.05)
                    #         uq = oe.loc[:, seq_name].dropna().quantile(0.95)
                    #         # axes[1, 1].scatter(mn, cval,
                    #         #                   marker="o", color="b", alpha=0.5,s=size)
                    #         # axes[1, 1].plot([lq, uq], [cval, cval],
                    #         #                color="b", alpha=0.5, lw=lw)
                    #
                    #         axesall_keep[ikeep].scatter(mn, cval,
                    #                                       marker="o", color=color, alpha=0.5, s=size)
                    #         axesall_keep[ikeep].plot([lq, uq], [cval, cval],
                    #                                    color=color, alpha=0.5, lw=lw)
                    #     else:
                    #
                    #         axesall_keep[ikeep].scatter(oe.loc[:, seq_name], [cval for _ in range(oe.shape[0])],
                    #                                     marker="o", color=color,
                    #                                     alpha=0.5, s=size, zorder=10)
                    #         axesall_keep[ikeep].scatter(oe.loc[:, seq_name], [cval for _ in range(oe.shape[0])],
                    #                                        marker="o", color=color,
                    #                                        alpha=0.5, s=size, zorder=10)



            ax = axesall_keep[ikeep]
            mn = min(ax.get_xlim()[0], ax.get_ylim()[0])
            mx = max(ax.get_xlim()[1], ax.get_ylim()[1])
            ax.plot([mn, mx], [mn, mx], "k--")

            ax.set_xlim(mn, mx)
            ax.set_ylim(mn, mx)

            axesall_keep[ikeep].set_aspect("equal")
            axesall_keep[ikeep].set_ylabel("complex")
            axesall_keep[ikeep].set_xlabel("simple")
            #axesall_keep[ikeep, 1].set_aspect("equal")
            #axesall_keep[ikeep, 1].set_ylabel("complex")
            #axesall_keep[ikeep, 1].set_xlabel("simple")

        plt.tight_layout()
        pdf.savefig(figall_keep)
        plt.close(figall_keep)










def plot_s_vs_s_pub_2(summarize=False, subdir=".", post_iter=None):
    include_est_states = False
    ognames = keep.copy()
    ognames.extend(forecast)
    label_dict = keep_dict
    label_dict.update(forecast_dict)

    # first rip thru all the dirs and load...
    s_b_dict = {}
    s_s_dict = {}
    s_s_est_dict = {}
    print("loading results...")

    for ireal in range(50):
        s_b_m_d = os.path.join(subdir,"monthly_model_files_master_{0}".format(ireal))
        s_s_m_d = os.path.join(subdir,"seq_monthly_model_files_master_{0}".format(ireal))

        if not os.path.exists(s_s_m_d) or not os.path.exists(s_s_m_d):
            break
        try:
            s_b_pst = pyemu.Pst(os.path.join(s_b_m_d, "freyberg.pst"))
            obs = s_b_pst.observation_data.loc[s_b_pst.nnz_obs_names,:]
            obs = obs.loc[obs.obsnme.str.contains("sfr_usecol:gage_1")]
            if obs.obsval.max() > 8000:
                continue
            obs = s_b_pst.observation_data
            obs = obs.loc[obs.obsnme.str.contains("hds_usecol:arrobs_head_k:0_i:33_j:7"), :]
            if obs.obsval.min() < 30:
                continue
            s_b_oe_pr = pd.read_csv(os.path.join(s_b_m_d, "freyberg.0.obs.csv"), index_col=0)
            log_cols = s_b_oe_pr.columns.map(lambda x: "mass" in x or "cnc" in x)
            s_b_oe_pr.loc[:,log_cols] = s_b_oe_pr.loc[:,log_cols].apply(np.log10)

            bpost_iter = s_b_pst.control_data.noptmax
            if post_iter is not None:
                bpost_iter = post_iter
            s_b_oe_pt = pd.read_csv(os.path.join(s_b_m_d, "freyberg.{0}.obs.csv".format(bpost_iter)),
                                    index_col=0)
            log_cols = s_b_oe_pt.columns.map(lambda x: "mass" in x or "cnc" in x)
            s_b_oe_pt.loc[:, log_cols] = s_b_oe_pt.loc[:, log_cols].apply(np.log10)


            s_s_pst = pyemu.Pst(os.path.join(s_s_m_d, "freyberg.pst"))
            seq_oe_files_pr = [f for f in os.listdir(s_s_m_d) if f.endswith("0.obs.csv") and f.startswith("freyberg")]
            spost_iter = s_s_pst.control_data.noptmax
            if post_iter is not None:
                spost_iter = post_iter
            seq_oe_files_pt = [f for f in os.listdir(s_s_m_d) if
                               f.endswith("{0}.obs.csv".format(spost_iter)) and f.startswith(
                                   "freyberg")]

            s_s_oe_dict_pr = {int(f.split(".")[1]): pd.read_csv(os.path.join(s_s_m_d, f), index_col=0) for f in
                              seq_oe_files_pr}
            s_s_oe_dict_pt = {int(f.split(".")[1]): pd.read_csv(os.path.join(s_s_m_d, f), index_col=0) for f in
                              seq_oe_files_pt}

            for key,df in s_s_oe_dict_pr.items():
                log_cols = df.columns.map(lambda x: "mass" in x or "cnc" in x)


                df.loc[:,log_cols] = df.loc[:,log_cols].apply(np.log10)
                df = df.replace(-np.Inf, np.nan)
                s_s_oe_dict_pr[key] = df

            for key, df in s_s_oe_dict_pt.items():
                log_cols = df.columns.map(lambda x: "mass" in x or "cnc" in x)
                df.loc[:, log_cols] = df.loc[:, log_cols].apply(np.log10)
                df = df.replace(-np.Inf, np.nan)
                s_s_oe_dict_pt[key] = df

            s_b_dict[ireal] = [s_b_pst,s_b_oe_pr,s_b_oe_pt]
            s_s_dict[ireal] = [s_s_pst,s_s_oe_dict_pr,s_s_oe_dict_pt]

            if include_est_states:
                # key these one cycle ahead since the posterior est states for this cycle are equiv to the prior sim states
                # of next cycle
                s_s_pe_dict_pt = {int(f.split(".")[1]) + 1: pd.read_csv(os.path.join(s_s_m_d, f.replace(".obs.",".par.")),
                                                                    index_col=0) for f in seq_oe_files_pt}
                s_s_est_dict[ireal] = s_s_pe_dict_pt
            print(ireal)
        except:
            break

    obs = s_s_pst.observation_data
    #sobs_to_sipar = obs.loc[pd.notna(obs.state_par_link),"state_par_link"].to_dict()
    #par = s_s_pst.parameter_data
    #sfpar = par.loc[pd.notna(par.state_par_link),:]
    #sipar_to_sfpar = {si:sf for si,sf in zip(sfpar.state_par_link,sfpar.parnme)}


    if len(s_b_dict) == 0:
        raise Exception()

    ireals = list(s_s_dict.keys())
    ireals.sort()
    sbobs_org = s_b_pst.observation_data
    print("plotting")
    size,lw=3,0.5
    pname = os.path.join(subdir,"s_vs_s_pub.pdf")
    if subdir != ".":
        pname = pname.replace(".pdf","_"+subdir+".pdf")
    #if post_iter is not None:
    #    pname = os.path.join(subdir,"s_vs_s_postiter_{0}.pdf".format(post_iter))

    is_1_lay = True
    if True in [True if "k:2" in o else False for o in s_b_pst.obs_names]:
        is_1_lay = False
    labels = ["surface-water flux ($\\frac{m^3}{d}$)","groundwater level ($m$)","groundwater level ($m$)","SW-GW flux ($\\frac{m^3}{d}$)"]
    sites = [keep[-1],keep[0],forecast[2],forecast[1]]
    names = ["sw_1","gw_1",
             "gw forecast","headwater forecast\n   "]
    with PdfPages(pname) as pdf:
        for ikeep, ogname in enumerate(sites):
            # if "mass" not in ogname:
            #    continue

            lab = labels[ikeep]
            k0ogname = ogname
            if is_1_lay:
                k0ogname = ogname.replace("k:2", "k:0")
            sgobs = sbobs_org.loc[sbobs_org.obsnme.str.contains(k0ogname), :].copy()
            sgobs = sgobs.loc[sgobs.obsnme.str.contains("_time"), :]
            sgobs.loc[:, "time"] = sgobs.time.apply(float)
            # sgobs = sgobs.loc[sgobs.time.apply(lambda x: x > 10000 and x < 10366),:]
            sgobs.sort_values(by="time", inplace=True)
            fig,axes = plt.subplots(2,2,figsize=(6,6))

            for itime, oname in enumerate(sgobs.obsnme):
                if itime != 12:
                    continue
                print(itime, oname)
                for ireal in ireals:
                    s_b_pst, s_b_oe_pr, s_b_oe_pt = s_b_dict[ireal]
                    sbobs = s_b_pst.observation_data
                    sgobs = sbobs.loc[sbobs.obsnme.str.contains(k0ogname), :].copy()

                    s_s_pst, s_s_oe_dict_pr, s_s_oe_dict_pt = s_s_dict[ireal]

                    cval = sgobs.loc[oname, "obsval"].copy()

                    mn = s_b_oe_pr.loc[:, oname].mean()
                    lq = s_b_oe_pr.loc[:, oname].quantile(0.05)
                    uq = s_b_oe_pr.loc[:, oname].quantile(0.95)
                    #axes[0, 0].scatter(mn, cval,
                    #                   marker="o", color="0.5", alpha=0.5,s=size)
                    #axes[0, 0].plot([lq,uq], [cval,cval],
                    #                   color="0.5", alpha=0.5,lw=lw)

                    mn = s_b_oe_pt.loc[:, oname].mean()
                    lq = s_b_oe_pt.loc[:, oname].quantile(0.05)
                    uq = s_b_oe_pt.loc[:, oname].quantile(0.95)
                    axes[0, 0].scatter(mn, cval,
                                       marker="o", color="b", alpha=0.5, s=size)
                    axes[0, 0].plot([lq, uq], [cval, cval],
                                    color="b", alpha=0.5, lw=lw)

                    seq_name = k0ogname
                    if "arrobs" not in k0ogname:
                        seq_name = k0ogname + "_time:10000.0"

                    if itime in s_s_oe_dict_pr:
                        oe = s_s_oe_dict_pr[itime]

                        mn = oe.loc[:, seq_name].dropna().mean()
                        lq = oe.loc[:, seq_name].dropna().quantile(0.05)
                        uq = oe.loc[:, seq_name].dropna().quantile(0.95)
                        axes[0, 1].scatter(mn, cval,
                                           marker="o", color="0.5", alpha=0.5, s=size)
                        axes[0, 1].plot([lq, uq], [cval, cval],
                                        color="0.5", alpha=0.5, lw=lw)

                        oe = s_s_oe_dict_pt[itime]

                        mn = oe.loc[:, seq_name].dropna().mean()
                        lq = oe.loc[:, seq_name].dropna().quantile(0.05)
                        uq = oe.loc[:, seq_name].dropna().quantile(0.95)
                        axes[0, 1].scatter(mn, cval,
                                           marker="o", color="b", alpha=0.5, s=size)
                        axes[0, 1].plot([lq, uq], [cval, cval],
                                        color="b", alpha=0.35, lw=lw)
            for itime, oname in enumerate(sgobs.obsnme):
                if itime != 24:
                    continue
                print(itime, oname)
                for ireal in ireals:
                    s_b_pst, s_b_oe_pr, s_b_oe_pt = s_b_dict[ireal]
                    sbobs = s_b_pst.observation_data
                    sgobs = sbobs.loc[sbobs.obsnme.str.contains(k0ogname), :].copy()

                    s_s_pst, s_s_oe_dict_pr, s_s_oe_dict_pt = s_s_dict[ireal]

                    cval = sgobs.loc[oname, "obsval"].copy()

                    mn = s_b_oe_pr.loc[:, oname].mean()
                    lq = s_b_oe_pr.loc[:, oname].quantile(0.05)
                    uq = s_b_oe_pr.loc[:, oname].quantile(0.95)
                    #axes[1, 0].scatter(mn, cval,
                    #                   marker="o", color="0.5", alpha=0.5,s=size)
                    #axes[1, 0].plot([lq,uq], [cval,cval],
                    #                   color="0.5", alpha=0.5,lw=lw)

                    mn = s_b_oe_pt.loc[:, oname].mean()
                    lq = s_b_oe_pt.loc[:, oname].quantile(0.05)
                    uq = s_b_oe_pt.loc[:, oname].quantile(0.95)
                    axes[1, 0].scatter(mn, cval,
                                     marker="o", color="b", alpha=0.5, s=size)
                    axes[1, 0].plot([lq, uq], [cval, cval],
                                    color="b", alpha=0.5, lw=lw)

                    seq_name = k0ogname
                    if "arrobs" not in k0ogname:
                        seq_name = k0ogname + "_time:10000.0"

                    if itime in s_s_oe_dict_pr:
                        oe = s_s_oe_dict_pr[itime]

                        mn = oe.loc[:, seq_name].dropna().mean()
                        lq = oe.loc[:, seq_name].dropna().quantile(0.05)
                        uq = oe.loc[:, seq_name].dropna().quantile(0.95)
                        axes[1, 1].scatter(mn, cval,
                                           marker="o", color="0.5", alpha=0.5, s=size)
                        axes[1, 1].plot([lq, uq], [cval, cval],
                                        color="0.5", alpha=0.5, lw=lw)

                        #oe = s_s_oe_dict_pt[itime]

                        # mn = oe.loc[:, seq_name].dropna().mean()
                        # lq = oe.loc[:, seq_name].dropna().quantile(0.05)
                        # uq = oe.loc[:, seq_name].dropna().quantile(0.95)
                        # axes[1, 1].scatter(mn, cval,
                        #                    marker="o", color="b", alpha=0.5, s=size)
                        # axes[1, 1].plot([lq, uq], [cval, cval],
                        #                 color="b", alpha=0.5, lw=lw)



            mn = 1.0e+10
            mx = -12.0e+10
            for ax in axes.flatten():
                mn = min(mn,ax.get_xlim()[0], ax.get_ylim()[0])
                mx = max(mx,ax.get_xlim()[1], ax.get_ylim()[1])

            for ax in axes.flatten():
                ax.plot([mn, mx], [mn, mx], "k--")
                ax.set_xlim(mn,mx)
                ax.set_ylim(mn,mx)
                ax.set_xlabel("simple {0}".format(lab))
                ax.set_ylabel("complex {0}".format(lab))
            axes[0,0].set_title("A) {0} batch cycle 13".format(names[ikeep]),loc="left")
            axes[0, 1].set_title("B) {0} sequential cycle 13".format(names[ikeep]),loc="left")
            axes[1, 0].set_title("C) {0} batch cycle 25".format(names[ikeep]),loc="left")
            axes[1, 1].set_title("D) {0} sequential cycle 25".format(names[ikeep]),loc="left")
            #fig.suptitle(names[ikeep])
            plt.tight_layout()
            pdf.savefig(fig)
            plt.close(fig)

def plot_obs_v_sim3(subdir=".",post_iter=None):
    """plot the results for daily, monthly batch and monthly sequential

     """
    c_m_d = "daily_model_files_master_prior"
    c_pst = pyemu.Pst(os.path.join(c_m_d, "freyberg.pst"))
    cobs = c_pst.observation_data
    # cobs = obs.loc[obs.obsnme.str.startswith("hds_usecol:arrobs_head_"), :]
    cobs.loc[:, "time"] = cobs.time.apply(float)
    c_oe = pd.read_csv(os.path.join(c_m_d, "freyberg.0.obs.csv"), index_col=0)
    cw_cols = c_oe.columns.map(lambda x: "mass" in x)
    c_oe.loc[:, cw_cols] = c_oe.loc[:, cw_cols].apply(np.log10)
    pname = os.path.join('.', "HESS_2022-170_supporting_information.pdf")
    if post_iter is not None:
        pname = os.path.join('.', "obs_v_sim_postier_{0}.pdf".format(post_iter))

    m = 1
    pp = PdfPages(pname)
    for ireal in range(20):
        s_b_m_d = os.path.join('.', "monthly_model_files_master_{0}".format(ireal))
        s_s_m_d = os.path.join('.', "seq_monthly_model_files_master_{0}".format(ireal))
        if not os.path.exists(s_s_m_d) or not os.path.exists(s_b_m_d):
            break
        try:
            s_b_pst = pyemu.Pst(os.path.join(s_b_m_d, "freyberg.pst"))
            s_b_oe_pr = pd.read_csv(os.path.join(s_b_m_d, "freyberg.0.obs.csv"), index_col=0)
            sw_cols = s_b_oe_pr.columns.map(lambda x: "mass" in x)
            s_b_oe_pr.loc[:, sw_cols] = s_b_oe_pr.loc[:, sw_cols].apply(np.log10)
            bpost_iter = s_b_pst.control_data.noptmax
            if post_iter is not None:
                bpost_iter = post_iter
            s_b_oe_pt = pd.read_csv(os.path.join(s_b_m_d, "freyberg.{0}.obs.csv".format(bpost_iter)),
                                    index_col=0)
            s_b_oe_pt.loc[:, sw_cols] = s_b_oe_pt.loc[:, sw_cols].apply(np.log10)
            s_s_pst = pyemu.Pst(os.path.join(s_s_m_d, "freyberg.pst"))
            seq_oe_files_pr = [f for f in os.listdir(s_s_m_d) if
                               f.endswith("0.obs.csv") and f.startswith("freyberg")]
            spost_iter = s_s_pst.control_data.noptmax
            if post_iter is not None:
                spost_iter = post_iter
            seq_oe_files_pt = [f for f in os.listdir(s_s_m_d) if
                               f.endswith("{0}.obs.csv".format(spost_iter)) and f.startswith("freyberg")]

            s_s_oe_dict_pr = {int(f.split(".")[1]): pd.read_csv(os.path.join(s_s_m_d, f), index_col=0) for f in
                              seq_oe_files_pr}
            s_s_oe_dict_pt = {int(f.split(".")[1]): pd.read_csv(os.path.join(s_s_m_d, f), index_col=0) for f in
                              seq_oe_files_pt}
            for key, df in s_s_oe_dict_pr.items():
                log_cols = df.columns.map(lambda x: "mass" in x)
                df.loc[:, log_cols] = df.loc[:, log_cols].apply(np.log10)
                s_s_oe_dict_pr[key] = df
            for key, df in s_s_oe_dict_pt.items():
                log_cols = df.columns.map(lambda x: "mass" in x)
                df.loc[:, log_cols] = df.loc[:, log_cols].apply(np.log10)
                s_s_oe_dict_pt[key] = df

        except:
            break

        ognames = keep.copy()
        ognames.extend(forecast)
        label_dict = keep_dict.copy()
        label_dict.update(forecast_dict)
        units_dict = keep_dict2.copy()
        units_dict.update(forecast_dict2)

        is_1_lay = True
        if True in [True if "k:2" in o else False for o in s_b_pst.obs_names]:
            is_1_lay = False

        for ogname in ognames:
            # if "cum" not in ogname:
            #    continue
            k0name = ogname
            if is_1_lay:
                k0ogname = ogname.replace("k:2", "k:0")
            fig, axes = plt.subplots(2, 1, figsize=(8, 8))
            cgobs = cobs.loc[cobs.obsnme.str.contains(ogname), :].copy()
            sgobs = s_b_pst.observation_data.loc[s_b_pst.observation_data.obsnme.str.contains(k0ogname), :].copy()
            if cgobs.shape[0] == 0:
                raise Exception("complex empty " + ogname)
            if sgobs.shape[0] == 0:
                raise Exception("batch empty " + k0ogname + "," + ogname)

            sgobs.loc[:, "time"] = sgobs.time.apply(float)
            cgobs.loc[:, "time"] = cgobs.time.apply(float)
            sgnzobs = sgobs.loc[sgobs.weight > 0, :].copy()

            sgobs.sort_values(by="time", inplace=True)
            cgobs.sort_values(by="time", inplace=True)

            ax = axes[0]
            ax.set_title("Figure S{0}A) batch formulation {1}, replicate {2}, coarse scenario".format(m, label_dict[ogname], c_oe.index[ireal]),
                         loc="left")
            ax.set_xlabel("Simulation Time (days)")
            ax.set_ylabel("{0}".format(units_dict[ogname]))
            obs = mpatches.Patch(color='red', label='Replicate Simulated Quantity')
            pr = mpatches.Patch(color='grey', label='Prior Realization Simulated Quantity')
            pt = mpatches.Patch(color='blue', label='Posterior Realization Simulated Quantity')
            plt.legend(handles=[obs, pr, pt])
            [ax.plot(sgobs.time, s_b_oe_pr.loc[idx, sgobs.obsnme], "0.5", lw=0.01, alpha=0.5) for idx in
             s_b_oe_pr.index]
            [ax.plot(sgobs.time, s_b_oe_pt.loc[idx, sgobs.obsnme], "b", lw=0.01, alpha=0.5) for idx in
             s_b_oe_pt.index]
            ax.plot(cgobs.time, c_oe.loc[c_oe.index[ireal], cgobs.obsnme], "r", lw=2.0, alpha=0.85)
            ax.scatter(sgnzobs.time, sgnzobs.obsval, marker="^", color="r")
            ax = axes[1]

            seq_name = k0ogname
            if "arrobs" not in k0ogname:
                seq_name = k0ogname + "_time:10000.0"
            print(ireal, seq_name)
            for itime, time in enumerate(sgobs.time):
                # itime += 1

                if itime in s_s_oe_dict_pr:
                    oe = s_s_oe_dict_pr[itime]
                    # print(oe.loc[:,seq_name])
                    ax.scatter([time for _ in range(oe.shape[0])], oe.loc[:, seq_name], marker=".", color="0.5",
                               alpha=0.5)
                if itime in s_s_oe_dict_pt:
                    oe = s_s_oe_dict_pt[itime]
                    # print(oe.loc[:,seq_name])
                    ax.scatter([time for _ in range(oe.shape[0])], oe.loc[:, seq_name], marker=".", color="b",
                               alpha=0.5)
            ax.plot(cgobs.time, c_oe.loc[c_oe.index[ireal], cgobs.obsnme], "r", lw=2.0, alpha=0.85)
            ax.scatter(sgnzobs.time, sgnzobs.obsval, marker="^", color="r")
            ax.set_title(
                "Figure S{0}B) sequential formulation {1}, replicate {2}, coarse scenario".format(m, label_dict[ogname], c_oe.index[ireal]),
                loc="left")
            ax.set_xlabel("Simulation Time (days)")
            ax.set_ylabel("{0}".format(units_dict[ogname]))
            obs = mpatches.Patch(color='red', label='Replicate Simulated Quantity')
            pr = mpatches.Patch(color='grey', label='Prior Realization Simulated Quantity')
            pt = mpatches.Patch(color='blue', label='Posterior Realization Simulated Quantity')
            plt.legend(handles=[obs, pr, pt])
            # if "gage" not in ogname:
            #    ax.set_ylim(30,ax.get_ylim()[1])
            mn = 1.0e+10
            mx = -1.0e+10
            for ax in axes.flatten():
                mn = min(ax.get_ylim()[0], mn)
                mx = max(ax.get_ylim()[1], mx)
            for ax in axes.flatten():
                ax.set_ylim(mn, mx)
            plt.tight_layout()
            plt.savefig(pp, format='pdf')
            # pp.savefig()
            m+=1

    c_m_d = "daily_model_files_master_prior"
    c_pst = pyemu.Pst(os.path.join(c_m_d, "freyberg.pst"))
    cobs = c_pst.observation_data
    # cobs = obs.loc[obs.obsnme.str.startswith("hds_usecol:arrobs_head_"), :]
    cobs.loc[:, "time"] = cobs.time.apply(float)
    c_oe = pd.read_csv(os.path.join(c_m_d, "freyberg.0.obs.csv"), index_col=0)
    cw_cols = c_oe.columns.map(lambda x: "mass" in x)
    c_oe.loc[:, cw_cols] = c_oe.loc[:, cw_cols].apply(np.log10)

    if post_iter is not None:
        pname = os.path.join(subdir, "obs_v_sim_postier_{0}.pdf".format(post_iter))

    # pp = PdfPages("HESS_2022-170_supporting_information2.pdf")
    for ireal in range(20):
        s_b_m_d = os.path.join(subdir, "monthly_model_files_master_{0}".format(ireal))
        s_s_m_d = os.path.join(subdir, "seq_monthly_model_files_master_{0}".format(ireal))
        if not os.path.exists(s_s_m_d) or not os.path.exists(s_b_m_d):
            break
        try:
            s_b_pst = pyemu.Pst(os.path.join(s_b_m_d, "freyberg.pst"))
            s_b_oe_pr = pd.read_csv(os.path.join(s_b_m_d, "freyberg.0.obs.csv"), index_col=0)
            sw_cols = s_b_oe_pr.columns.map(lambda x: "mass" in x)
            s_b_oe_pr.loc[:, sw_cols] = s_b_oe_pr.loc[:, sw_cols].apply(np.log10)
            bpost_iter = s_b_pst.control_data.noptmax
            if post_iter is not None:
                bpost_iter = post_iter
            s_b_oe_pt = pd.read_csv(os.path.join(s_b_m_d, "freyberg.{0}.obs.csv".format(bpost_iter)),
                                    index_col=0)
            s_b_oe_pt.loc[:, sw_cols] = s_b_oe_pt.loc[:, sw_cols].apply(np.log10)
            s_s_pst = pyemu.Pst(os.path.join(s_s_m_d, "freyberg.pst"))
            seq_oe_files_pr = [f for f in os.listdir(s_s_m_d) if
                               f.endswith("0.obs.csv") and f.startswith("freyberg")]
            spost_iter = s_s_pst.control_data.noptmax
            if post_iter is not None:
                spost_iter = post_iter
            seq_oe_files_pt = [f for f in os.listdir(s_s_m_d) if
                               f.endswith("{0}.obs.csv".format(spost_iter)) and f.startswith("freyberg")]

            s_s_oe_dict_pr = {int(f.split(".")[1]): pd.read_csv(os.path.join(s_s_m_d, f), index_col=0) for f in
                              seq_oe_files_pr}
            s_s_oe_dict_pt = {int(f.split(".")[1]): pd.read_csv(os.path.join(s_s_m_d, f), index_col=0) for f in
                              seq_oe_files_pt}
            for key, df in s_s_oe_dict_pr.items():
                log_cols = df.columns.map(lambda x: "mass" in x)
                df.loc[:, log_cols] = df.loc[:, log_cols].apply(np.log10)
                s_s_oe_dict_pr[key] = df
            for key, df in s_s_oe_dict_pt.items():
                log_cols = df.columns.map(lambda x: "mass" in x)
                df.loc[:, log_cols] = df.loc[:, log_cols].apply(np.log10)
                s_s_oe_dict_pt[key] = df

        except:
            break

        ognames = keep.copy()
        ognames.extend(forecast)
        label_dict = keep_dict.copy()
        label_dict.update(forecast_dict)

        is_1_lay = True
        if True in [True if "k:2" in o else False for o in s_b_pst.obs_names]:
            is_1_lay = False

        for ogname in ognames:
            # if "cum" not in ogname:
            #    continue
            k0name = ogname
            if is_1_lay:
                k0ogname = ogname.replace("k:2", "k:0")
            fig, axes = plt.subplots(2, 1, figsize=(8, 8))
            cgobs = cobs.loc[cobs.obsnme.str.contains(ogname), :].copy()
            sgobs = s_b_pst.observation_data.loc[s_b_pst.observation_data.obsnme.str.contains(k0ogname), :].copy()
            if cgobs.shape[0] == 0:
                raise Exception("complex empty " + ogname)
            if sgobs.shape[0] == 0:
                raise Exception("batch empty " + k0ogname + "," + ogname)

            sgobs.loc[:, "time"] = sgobs.time.apply(float)
            cgobs.loc[:, "time"] = cgobs.time.apply(float)
            sgnzobs = sgobs.loc[sgobs.weight > 0, :].copy()

            sgobs.sort_values(by="time", inplace=True)
            cgobs.sort_values(by="time", inplace=True)

            ax = axes[0]
            ax.set_title("Figure S{0}A) batch formulation {1}, replicate {2}, fixed well scenario".format(m, label_dict[ogname], c_oe.index[ireal]),
                         loc="left")
            ax.set_xlabel("Simulation Time (days)")
            ax.set_ylabel("{0}".format(units_dict[ogname]))
            obs = mpatches.Patch(color='red', label='Replicate Simulated Quantity')
            pr = mpatches.Patch(color='grey', label='Prior Realization Simulated Quantity')
            pt = mpatches.Patch(color='blue', label='Posterior Realization Simulated Quantity')
            plt.legend(handles=[obs, pr, pt])
            [ax.plot(sgobs.time, s_b_oe_pr.loc[idx, sgobs.obsnme], "0.5", lw=0.01, alpha=0.5) for idx in
             s_b_oe_pr.index]
            [ax.plot(sgobs.time, s_b_oe_pt.loc[idx, sgobs.obsnme], "b", lw=0.01, alpha=0.5) for idx in
             s_b_oe_pt.index]
            ax.plot(cgobs.time, c_oe.loc[c_oe.index[ireal], cgobs.obsnme], "r", lw=2.0, alpha=0.85)
            ax.scatter(sgnzobs.time, sgnzobs.obsval, marker="^", color="r")
            ax = axes[1]

            seq_name = k0ogname
            if "arrobs" not in k0ogname:
                seq_name = k0ogname + "_time:10000.0"
            print(ireal, seq_name)
            for itime, time in enumerate(sgobs.time):
                # itime += 1

                if itime in s_s_oe_dict_pr:
                    oe = s_s_oe_dict_pr[itime]
                    # print(oe.loc[:,seq_name])
                    ax.scatter([time for _ in range(oe.shape[0])], oe.loc[:, seq_name], marker=".", color="0.5",
                               alpha=0.5)
                if itime in s_s_oe_dict_pt:
                    oe = s_s_oe_dict_pt[itime]
                    # print(oe.loc[:,seq_name])
                    ax.scatter([time for _ in range(oe.shape[0])], oe.loc[:, seq_name], marker=".", color="b",
                               alpha=0.5)
            ax.plot(cgobs.time, c_oe.loc[c_oe.index[ireal], cgobs.obsnme], "r", lw=2.0, alpha=0.85)
            ax.scatter(sgnzobs.time, sgnzobs.obsval, marker="^", color="r")
            ax.set_title(
                "Figure S{0}B) sequential formulation {1}, replicate {2}, fixed well scenario".format(m, label_dict[ogname], c_oe.index[ireal]),
                loc="left")
            ax.set_xlabel("Simulation Time (days)")
            ax.set_ylabel("{0}".format(units_dict[ogname]))
            obs = mpatches.Patch(color='red', label='Replicate Simulated Quantity')
            pr = mpatches.Patch(color='grey', label='Prior Realization Simulated Quantity')
            pt = mpatches.Patch(color='blue', label='Posterior Realization Simulated Quantity')
            plt.legend(handles=[obs, pr, pt])
            # if "gage" not in ogname:
            #    ax.set_ylim(30,ax.get_ylim()[1])
            mn = 1.0e+10
            mx = -1.0e+10
            for ax in axes.flatten():
                mn = min(ax.get_ylim()[0], mn)
                mx = max(ax.get_ylim()[1], mx)
            for ax in axes.flatten():
                ax.set_ylim(mn, mx)
            plt.tight_layout()
            plt.savefig(pp, format='pdf')
            # pp.savefig()
            m+=1

    pp.close()

def plot_obs_v_sim_pub(subdir=".",post_iter=None):
    """plot the results for daily, monthly batch and monthly sequential

     """
    c_m_d = "daily_model_files_master_prior"
    c_pst = pyemu.Pst(os.path.join(c_m_d, "freyberg.pst"))
    cobs = c_pst.observation_data
    # cobs = obs.loc[obs.obsnme.str.startswith("hds_usecol:arrobs_head_"), :]
    cobs.loc[:, "time"] = cobs.time.apply(float)
    c_oe = pd.read_csv(os.path.join(c_m_d, "freyberg.0.obs.csv"), index_col=0)
    cw_cols = c_oe.columns.map(lambda x: "mass" in x)
    c_oe.loc[:, cw_cols] = c_oe.loc[:, cw_cols].apply(np.log10)
    pname = os.path.join('.', "obs_v_sim_pub.pdf")
    if post_iter is not None:
        pname = os.path.join('.', "obs_v_sim_postier_{0}.pdf".format(post_iter))

    m = 1
    pp = PdfPages(pname)
    for ireal in range(12,13):
        s_b_m_d = os.path.join('.', "monthly_model_files_master_{0}".format(ireal))
        s_s_m_d = os.path.join('.', "seq_monthly_model_files_master_{0}".format(ireal))
        if not os.path.exists(s_s_m_d) or not os.path.exists(s_b_m_d):
            break
        try:
            s_b_pst = pyemu.Pst(os.path.join(s_b_m_d, "freyberg.pst"))
            s_b_oe_pr = pd.read_csv(os.path.join(s_b_m_d, "freyberg.0.obs.csv"), index_col=0)
            sw_cols = s_b_oe_pr.columns.map(lambda x: "mass" in x)
            s_b_oe_pr.loc[:, sw_cols] = s_b_oe_pr.loc[:, sw_cols].apply(np.log10)
            bpost_iter = s_b_pst.control_data.noptmax
            if post_iter is not None:
                bpost_iter = post_iter
            s_b_oe_pt = pd.read_csv(os.path.join(s_b_m_d, "freyberg.{0}.obs.csv".format(bpost_iter)),
                                    index_col=0)
            s_b_oe_pt.loc[:, sw_cols] = s_b_oe_pt.loc[:, sw_cols].apply(np.log10)
            s_s_pst = pyemu.Pst(os.path.join(s_s_m_d, "freyberg.pst"))
            seq_oe_files_pr = [f for f in os.listdir(s_s_m_d) if
                               f.endswith("0.obs.csv") and f.startswith("freyberg")]
            spost_iter = s_s_pst.control_data.noptmax
            if post_iter is not None:
                spost_iter = post_iter
            seq_oe_files_pt = [f for f in os.listdir(s_s_m_d) if
                               f.endswith("{0}.obs.csv".format(spost_iter)) and f.startswith("freyberg")]

            s_s_oe_dict_pr = {int(f.split(".")[1]): pd.read_csv(os.path.join(s_s_m_d, f), index_col=0) for f in
                              seq_oe_files_pr}
            s_s_oe_dict_pt = {int(f.split(".")[1]): pd.read_csv(os.path.join(s_s_m_d, f), index_col=0) for f in
                              seq_oe_files_pt}
            for key, df in s_s_oe_dict_pr.items():
                log_cols = df.columns.map(lambda x: "mass" in x)
                df.loc[:, log_cols] = df.loc[:, log_cols].apply(np.log10)
                s_s_oe_dict_pr[key] = df
            for key, df in s_s_oe_dict_pt.items():
                log_cols = df.columns.map(lambda x: "mass" in x)
                df.loc[:, log_cols] = df.loc[:, log_cols].apply(np.log10)
                s_s_oe_dict_pt[key] = df

        except:
            break

        ognames = keep.copy()
        ognames.extend(forecast)
        label_dict = keep_dict.copy()
        label_dict.update(forecast_dict)
        units_dict = keep_dict2.copy()
        units_dict.update(forecast_dict2)

        is_1_lay = True
        if True in [True if "k:2" in o else False for o in s_b_pst.obs_names]:
            is_1_lay = False

        for ogname in ognames:
            print('plotting', ogname)
            # if "cum" not in ogname:
            #    continue
            k0name = ogname
            if is_1_lay:
                k0ogname = ogname.replace("k:2", "k:0")
            fig, axes = plt.subplots(2, 1, figsize=(8, 8))
            cgobs = cobs.loc[cobs.obsnme.str.contains(ogname), :].copy()
            sgobs = s_b_pst.observation_data.loc[s_b_pst.observation_data.obsnme.str.contains(k0ogname), :].copy()
            if cgobs.shape[0] == 0:
                raise Exception("complex empty " + ogname)
            if sgobs.shape[0] == 0:
                raise Exception("batch empty " + k0ogname + "," + ogname)

            sgobs.loc[:, "time"] = sgobs.time.apply(float)
            cgobs.loc[:, "time"] = cgobs.time.apply(float)
            sgnzobs = sgobs.loc[sgobs.weight > 0, :].copy()

            sgobs.sort_values(by="time", inplace=True)
            cgobs.sort_values(by="time", inplace=True)

            ax = axes[0]
            ax.set_title("A) batch formulation {0}, replicate {1}, coarse scenario".format(label_dict[ogname], c_oe.index[ireal]),
                         loc="left")
            ax.set_xlabel("Simulation Time (days)")
            ax.set_ylabel("{0}".format(units_dict[ogname]))
            obs = mpatches.Patch(color='red', label='Replicate Simulated Quantity')
            pr = mpatches.Patch(color='grey', label='Prior Realization Simulated Quantity')
            pt = mpatches.Patch(color='blue', label='Posterior Realization Simulated Quantity')
            plt.legend(handles=[obs, pr, pt])
            [ax.plot(sgobs.time, s_b_oe_pr.loc[idx, sgobs.obsnme], "0.5", lw=0.01, alpha=0.5) for idx in
             s_b_oe_pr.index]
            [ax.plot(sgobs.time, s_b_oe_pt.loc[idx, sgobs.obsnme], "b", lw=0.01, alpha=0.5) for idx in
             s_b_oe_pt.index]
            ax.plot(cgobs.time, c_oe.loc[c_oe.index[ireal], cgobs.obsnme], "r", lw=2.0, alpha=0.85)
            ax.scatter(sgnzobs.time, sgnzobs.obsval, marker="^", color="r")
            ax = axes[1]

            seq_name = k0ogname
            if "arrobs" not in k0ogname:
                seq_name = k0ogname + "_time:10000.0"
            print(ireal, seq_name)
            for itime, time in enumerate(sgobs.time):
                # itime += 1

                if itime in s_s_oe_dict_pr:
                    oe = s_s_oe_dict_pr[itime]
                    # print(oe.loc[:,seq_name])
                    ax.scatter([time for _ in range(oe.shape[0])], oe.loc[:, seq_name], marker=".", color="0.5",
                               alpha=0.5)
                if itime in s_s_oe_dict_pt:
                    oe = s_s_oe_dict_pt[itime]
                    # print(oe.loc[:,seq_name])
                    ax.scatter([time for _ in range(oe.shape[0])], oe.loc[:, seq_name], marker=".", color="b",
                               alpha=0.5)
            ax.plot(cgobs.time, c_oe.loc[c_oe.index[ireal], cgobs.obsnme], "r", lw=2.0, alpha=0.85)
            ax.scatter(sgnzobs.time, sgnzobs.obsval, marker="^", color="r")
            ax.set_title(
                "B) sequential formulation {0}, replicate {1}, coarse scenario".format(label_dict[ogname], c_oe.index[ireal]),
                loc="left")
            ax.set_xlabel("Simulation Time (days)")
            ax.set_ylabel("{0}".format(units_dict[ogname]))
            obs = mpatches.Patch(color='red', label='Replicate Simulated Quantity')
            pr = mpatches.Patch(color='grey', label='Prior Realization Simulated Quantity')
            pt = mpatches.Patch(color='blue', label='Posterior Realization Simulated Quantity')
            plt.legend(handles=[obs, pr, pt])
            # if "gage" not in ogname:
            #    ax.set_ylim(30,ax.get_ylim()[1])
            mn = 1.0e+10
            mx = -1.0e+10
            for ax in axes.flatten():
                mn = min(ax.get_ylim()[0], mn)
                mx = max(ax.get_ylim()[1], mx)
            for ax in axes.flatten():
                ax.set_ylim(mn, mx)
            plt.tight_layout()
            plt.savefig(pp, format='pdf')
            # pp.savefig()
        m+=1
    # pp.close()

    c_m_d = "daily_model_files_master_prior"
    c_pst = pyemu.Pst(os.path.join(c_m_d, "freyberg.pst"))
    cobs = c_pst.observation_data
    # cobs = obs.loc[obs.obsnme.str.startswith("hds_usecol:arrobs_head_"), :]
    cobs.loc[:, "time"] = cobs.time.apply(float)
    c_oe = pd.read_csv(os.path.join(c_m_d, "freyberg.0.obs.csv"), index_col=0)
    cw_cols = c_oe.columns.map(lambda x: "mass" in x)
    c_oe.loc[:, cw_cols] = c_oe.loc[:, cw_cols].apply(np.log10)

    if post_iter is not None:
        pname = os.path.join(subdir, "obs_v_sim_postier_{0}.pdf".format(post_iter))

    # pp = PdfPages("HESS_2022-170_supporting_information2.pdf")
    for ireal in range(12,13):
        s_b_m_d = os.path.join(subdir, "monthly_model_files_master_{0}".format(ireal))
        s_s_m_d = os.path.join(subdir, "seq_monthly_model_files_master_{0}".format(ireal))
        if not os.path.exists(s_s_m_d) or not os.path.exists(s_b_m_d):
            break
        try:
            s_b_pst = pyemu.Pst(os.path.join(s_b_m_d, "freyberg.pst"))
            s_b_oe_pr = pd.read_csv(os.path.join(s_b_m_d, "freyberg.0.obs.csv"), index_col=0)
            sw_cols = s_b_oe_pr.columns.map(lambda x: "mass" in x)
            s_b_oe_pr.loc[:, sw_cols] = s_b_oe_pr.loc[:, sw_cols].apply(np.log10)
            bpost_iter = s_b_pst.control_data.noptmax
            if post_iter is not None:
                bpost_iter = post_iter
            s_b_oe_pt = pd.read_csv(os.path.join(s_b_m_d, "freyberg.{0}.obs.csv".format(bpost_iter)),
                                    index_col=0)
            s_b_oe_pt.loc[:, sw_cols] = s_b_oe_pt.loc[:, sw_cols].apply(np.log10)
            s_s_pst = pyemu.Pst(os.path.join(s_s_m_d, "freyberg.pst"))
            seq_oe_files_pr = [f for f in os.listdir(s_s_m_d) if
                               f.endswith("0.obs.csv") and f.startswith("freyberg")]
            spost_iter = s_s_pst.control_data.noptmax
            if post_iter is not None:
                spost_iter = post_iter
            seq_oe_files_pt = [f for f in os.listdir(s_s_m_d) if
                               f.endswith("{0}.obs.csv".format(spost_iter)) and f.startswith("freyberg")]

            s_s_oe_dict_pr = {int(f.split(".")[1]): pd.read_csv(os.path.join(s_s_m_d, f), index_col=0) for f in
                              seq_oe_files_pr}
            s_s_oe_dict_pt = {int(f.split(".")[1]): pd.read_csv(os.path.join(s_s_m_d, f), index_col=0) for f in
                              seq_oe_files_pt}
            for key, df in s_s_oe_dict_pr.items():
                log_cols = df.columns.map(lambda x: "mass" in x)
                df.loc[:, log_cols] = df.loc[:, log_cols].apply(np.log10)
                s_s_oe_dict_pr[key] = df
            for key, df in s_s_oe_dict_pt.items():
                log_cols = df.columns.map(lambda x: "mass" in x)
                df.loc[:, log_cols] = df.loc[:, log_cols].apply(np.log10)
                s_s_oe_dict_pt[key] = df

        except:
            break

        ognames = keep.copy()
        ognames.extend(forecast)
        label_dict = keep_dict.copy()
        label_dict.update(forecast_dict)

        is_1_lay = True
        if True in [True if "k:2" in o else False for o in s_b_pst.obs_names]:
            is_1_lay = False

        for ogname in ognames:
            # if "cum" not in ogname:
            #    continue
            k0name = ogname
            if is_1_lay:
                k0ogname = ogname.replace("k:2", "k:0")
            fig, axes = plt.subplots(2, 1, figsize=(8, 8))
            cgobs = cobs.loc[cobs.obsnme.str.contains(ogname), :].copy()
            sgobs = s_b_pst.observation_data.loc[s_b_pst.observation_data.obsnme.str.contains(k0ogname), :].copy()
            if cgobs.shape[0] == 0:
                raise Exception("complex empty " + ogname)
            if sgobs.shape[0] == 0:
                raise Exception("batch empty " + k0ogname + "," + ogname)

            sgobs.loc[:, "time"] = sgobs.time.apply(float)
            cgobs.loc[:, "time"] = cgobs.time.apply(float)
            sgnzobs = sgobs.loc[sgobs.weight > 0, :].copy()

            sgobs.sort_values(by="time", inplace=True)
            cgobs.sort_values(by="time", inplace=True)

            ax = axes[0]
            ax.set_title("A) batch formulation {1}, replicate {2}, fixed well scenario".format(m, label_dict[ogname], c_oe.index[ireal]),
                         loc="left")
            ax.set_xlabel("Simulation Time (days)")
            ax.set_ylabel("{0}".format(units_dict[ogname]))
            obs = mpatches.Patch(color='red', label='Replicate Simulated Quantity')
            pr = mpatches.Patch(color='grey', label='Prior Realization Simulated Quantity')
            pt = mpatches.Patch(color='blue', label='Posterior Realization Simulated Quantity')
            plt.legend(handles=[obs, pr, pt])
            [ax.plot(sgobs.time, s_b_oe_pr.loc[idx, sgobs.obsnme], "0.5", lw=0.01, alpha=0.5) for idx in
             s_b_oe_pr.index]
            [ax.plot(sgobs.time, s_b_oe_pt.loc[idx, sgobs.obsnme], "b", lw=0.01, alpha=0.5) for idx in
             s_b_oe_pt.index]
            ax.plot(cgobs.time, c_oe.loc[c_oe.index[ireal], cgobs.obsnme], "r", lw=2.0, alpha=0.85)
            ax.scatter(sgnzobs.time, sgnzobs.obsval, marker="^", color="r")
            ax = axes[1]

            seq_name = k0ogname
            if "arrobs" not in k0ogname:
                seq_name = k0ogname + "_time:10000.0"
            print(ireal, seq_name)
            for itime, time in enumerate(sgobs.time):
                # itime += 1

                if itime in s_s_oe_dict_pr:
                    oe = s_s_oe_dict_pr[itime]
                    # print(oe.loc[:,seq_name])
                    ax.scatter([time for _ in range(oe.shape[0])], oe.loc[:, seq_name], marker=".", color="0.5",
                               alpha=0.5)
                if itime in s_s_oe_dict_pt:
                    oe = s_s_oe_dict_pt[itime]
                    # print(oe.loc[:,seq_name])
                    ax.scatter([time for _ in range(oe.shape[0])], oe.loc[:, seq_name], marker=".", color="b",
                               alpha=0.5)
            ax.plot(cgobs.time, c_oe.loc[c_oe.index[ireal], cgobs.obsnme], "r", lw=2.0, alpha=0.85)
            ax.scatter(sgnzobs.time, sgnzobs.obsval, marker="^", color="r")
            ax.set_title(
                "B) sequential formulation {1}, replicate {2}, fixed well scenario".format(m, label_dict[ogname], c_oe.index[ireal]),
                loc="left")
            ax.set_xlabel("Simulation Time (days)")
            ax.set_ylabel("{0}".format(units_dict[ogname]))
            obs = mpatches.Patch(color='red', label='Replicate Simulated Quantity')
            pr = mpatches.Patch(color='grey', label='Prior Realization Simulated Quantity')
            pt = mpatches.Patch(color='blue', label='Posterior Realization Simulated Quantity')
            plt.legend(handles=[obs, pr, pt])
            # if "gage" not in ogname:
            #    ax.set_ylim(30,ax.get_ylim()[1])
            mn = 1.0e+10
            mx = -1.0e+10
            for ax in axes.flatten():
                mn = min(ax.get_ylim()[0], mn)
                mx = max(ax.get_ylim()[1], mx)
            for ax in axes.flatten():
                ax.set_ylim(mn, mx)
            plt.tight_layout()
            plt.savefig(pp, format='pdf')
            # pp.savefig()
        m+=1

    pp.close()

if __name__ == "__main__":


    #### MAIN WORKFLOW ####
    #coarse scenario
    sync_phase(s_d = "monthly_model_files_1lyr_trnsprt_org")
    add_new_stress(m_d_org = "monthly_model_files_1lyr_trnsprt")
    c_d = setup_interface("daily_model_files_trnsprt_newstress",num_reals=50)
    b_d = setup_interface("monthly_model_files_1lyr_trnsprt_newstress",num_reals=50)
    s_d = monthly_ies_to_da(b_d,include_est_states=False)

    m_c_d = run_complex_prior_mc(c_d,num_workers=10)

    b_d = map_complex_to_simple_bat("daily_model_files_master_prior",b_d,0)
    s_d = map_simple_bat_to_seq(b_d,"seq_monthly_model_files_template")

    compare_mf6_freyberg(num_workers=10, num_replicates=50,num_reals=50,use_sim_states=True,
                       run_ies=True,run_da=True,adj_init_states=True)

    #fixed well scenario
    sync_phase(s_d = "monthly_model_files_1lyr_trnsprt_org")
    add_new_stress(m_d_org = "monthly_model_files_1lyr_trnsprt")
    b_d = setup_interface("monthly_model_files_1lyr_trnsprt_newstress",num_reals=50)
    reduce_simple_pars(b_d)
    s_d = monthly_ies_to_da(b_d,include_est_states=False)

    b_d = map_complex_to_simple_bat("daily_model_files_master_prior",b_d,0)
    s_d = map_simple_bat_to_seq(b_d,"seq_monthly_model_files_template")

    compare_mf6_freyberg(num_workers=10, num_replicates=50,num_reals=50,use_sim_states=True,
                       run_ies=True,run_da=True,adj_init_states=True)


    # plotting
    plot_domain()
    plot_obs_v_sim_pub(subdir="missing_wel_pars")
    plot_obs_v_sim3(subdir="missing_wel_pars")
    plot_s_vs_s_pub_2(summarize=True)
    plot_s_vs_s_pub_2(summarize=True,subdir="missing_wel_pars")

    exit()

    #other fxns
    # m_b_d, m_s_d = run_batch_seq_prior_monte_carlo(b_d, s_d)
    #plot_prior_mc_all()
    #plot_prior_mc_pub()
    #plot_prior_mc_pub(subdir="missing_wel_pars")
    #invest()
    #clean_results("naive_50reals_eststates")
    #make_prop_histograms()
    #plot_s_vs_s_pub(summarize=True)
    #plot_s_vs_s(summarize=True)#,subdir="missing_wel_pars")