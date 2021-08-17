# script to run paired simple complex analysis of MF6 Freyberg model to compare batch and sequential DA

import os
import sys
import shutil
import platform
import numpy as np
import pandas as pd
import platform
import pyemu
import matplotlib.pyplot as plt
import flopy
from matplotlib.backends.backend_pdf import PdfPages

plt.rcParams.update({'font.size': 12})

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

keep = ['arrobs_head_k:0_i:22_j:15', 'arrobs_head_k:2_i:2_j:9', 'arrobs_head_k:2_i:33_j:7', 'sfr_usecol:gage_1']

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


def compare_mf6_freyberg(num_workers=10):
    for ireal in range(100):
        complex_dir = os.path.join('daily_model_files_master_prior')
        bat_dir = os.path.join('monthly_model_files_template')
        seq_dir = os.path.join('seq_monthly_model_files_template')

        ies_t_d = map_complex_to_simple_bat(complex_dir,bat_dir,ireal)
        da_t_d = map_simple_bat_to_seq(ies_t_d,seq_dir)

        # run batch and sequential simple models
        # ies stuff
        ies_pst = pyemu.Pst(os.path.join(ies_t_d, "freyberg.pst"))

        # prep that prior ensemble for da
        da_pst = pyemu.Pst(os.path.join(da_t_d, "freyberg.pst"))

        ies_pst.pestpp_options["ies_par_en"] = "prior.jcb"
        # ies_pe = pyemu.ParameterEnsemble.from_binary(pst=ies_pst, filename=os.path.join(ies_t_d, "prior.jcb"))
        # da_pe = pyemu.ParameterEnsemble.from_gaussian_draw(pst=da_pst,
        #                                                    cov=pyemu.Cov.from_parameter_data(da_pst),
        #                                                    num_reals=ies_pe.shape[0])
        # da_pe.index = ies_pe.index
        # d = set(da_pe.columns.tolist()).symmetric_difference(set(ies_pe.columns.tolist()))
        # print(d)
        # da_pe.loc[:, ies_pe.columns] = ies_pe.values
        # da_pe.to_binary(os.path.join(da_t_d, "da_prior.jcb"))
        # da_pst.pestpp_options["ies_par_en"] = "da_prior.jcb"

        # set pestpp options for batch da
        ies_pst.pestpp_options.pop("ies_num_reals", None)
        ies_pst.pestpp_options.pop("da_num_reals", None)
        ies_pst.pestpp_options["ies_no_noise"] = False
        ies_pst.pestpp_options["ies_verbose_level"] = 1
        ies_pst.pestpp_options.pop("ies_localizer", None)
        ies_pst.pestpp_options["ies_autoadaloc"] = False
        ies_pst.pestpp_options["ies_save_lambda_en"] = False
        ies_pst.pestpp_options["ies_drop_conflicts"] = False
        ies_pst.pestpp_options["ies_num_reals"] = 100
        ies_pst.pestpp_options["ies_use_mda"] = False
        ies_pst.control_data.noptmax = 3
        ies_pst.write(os.path.join(ies_t_d, "freyberg.pst"), version=2)

        # set pestpp options for sequential da
        da_pst.pestpp_options.pop("da_num_reals", None)
        da_pst.pestpp_options.pop("ies_num_reals", None)
        da_pst.pestpp_options["ies_no_noise"] = False
        da_pst.pestpp_options["ies_verbose_level"] = 1
        da_pst.pestpp_options.pop("ies_localizer", None)
        da_pst.pestpp_options["ies_autoadaloc"] = False
        da_pst.pestpp_options["ies_save_lambda_en"] = False
        da_pst.pestpp_options["ies_drop_conflicts"] = False
        da_pst.pestpp_options["ies_num_reals"] = 100
        da_pst.pestpp_options["ies_use_mda"] = False
        da_pst.control_data.noptmax = 3
        da_pst.write(os.path.join(da_t_d, "freyberg.pst"), version=2)

        # run da          
        m_da_dir = da_t_d.replace("template","master")

        pyemu.os_utils.start_workers(da_t_d, 'pestpp-da', "freyberg.pst", port=port,
                                     num_workers=num_workers, master_dir=m_da_dir, verbose=True)

        shutil.rmtree(da_t_d)

        # run ies  
        m_ies_dir = ies_t_d.replace("template","master")

        pyemu.os_utils.start_workers(ies_t_d, 'pestpp-ies', "freyberg.pst", port=port,
                                     num_workers=4, master_dir=m_ies_dir, verbose=True)

        shutil.rmtree(ies_t_d)
        break


def run_complex_prior_mc(c_t):
    pyemu.os_utils.start_workers(c_t, "pestpp-ies", "freyberg.pst", num_workers=10, worker_root=".",
                                 master_dir=c_t.replace("template", "master_prior"))


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
    start_date = pd.to_datetime('20151231', format='%Y%m%d')
    c_d = 'complex_master'
    redis_fac = 3.

    # load in obs ensemble
    oe_f = pd.read_csv(os.path.join(c_d, "freyberg.0.obs.csv"), index_col=0)
    oe_f = oe_f.T

    # process obs
    hds_f = oe_f.loc[oe_f.index.to_series().apply(lambda x: x.startswith("hds")), :].copy()
    hds_f.loc[:, "k"] = hds_f.index.to_series().apply(lambda x: int(x.split('_')[2]))
    hds_f.loc[:, "i"] = hds_f.index.to_series().apply(lambda x: int(x.split('_')[3]))
    hds_f.loc[:, "j"] = hds_f.index.to_series().apply(lambda x: int(x.split('_')[4]))
    hds_f.loc[:, "int_time"] = hds_f.index.to_series().apply(lambda x: float(x.split('_')[-1].split(':')[1]))
    time = hds_f.loc[:, "int_time"]
    time = pd.to_timedelta(time.values - 1, unit='D')
    hds_f.loc[:, "org_time"] = start_date + time.values
    hds_f.loc[:, "org_time"] = hds_f.loc[:, "org_time"].apply(lambda x: x.strftime('%Y%m%d'))
    hds_f.loc[:, "org_i"] = (hds_f.i / redis_fac).apply(np.int)
    hds_f.loc[:, "org_j"] = (hds_f.j / redis_fac).apply(np.int)
    hds_f.loc[:, "org_obgnme"] = hds_f.apply(lambda x: "trgw_{0}_{1}_{2}".format(x.k, x.org_i, x.org_j),
                                             axis=1)

    sfr_f = oe_f.loc[oe_f.index.to_series().apply(lambda x: x.startswith("sfr")), :].copy()
    sfr_f.loc[:, "int_time"] = sfr_f.index.to_series().apply(lambda x: float(x.split('_')[-1].split(':')[1]))
    type = sfr_f.index.to_series().apply(lambda x: x.split(':')[1].split('_')[0])
    type = pd.DataFrame(type)
    # type = type.replace('gage', 'gage_1')
    sfr_f.loc[:, "type"] = type.values
    time = sfr_f.loc[:, "int_time"]
    time = pd.to_timedelta(time.values - 1, unit='D')
    sfr_f.loc[:, "org_time"] = start_date + time.values
    sfr_f.loc[:, "org_time"] = sfr_f.loc[:, "org_time"].apply(lambda x: x.strftime('%Y%m%d'))
    sfr_f.loc[:, "org_obgnme"] = sfr_f.apply(lambda x: "{0}".format(x.type), axis=1)

    frames = [hds_f, sfr_f]

    complex_obs = pd.concat(frames)

    with PdfPages(os.path.join("obs_v_sim.pdf")) as pdf:
        for i in range(100):
            da_m_d = os.path.join("simple_master2_da_{0}".format(i))
            ies_m_d = os.path.join("simple_master2_ies_{0}".format(i))
            ies_case = "freyberg6_run_ies"
            da_case = "freyberg6_run_da"

            ies_pst = pyemu.Pst(os.path.join(ies_m_d, ies_case + ".pst"))
            ies_obs = ies_pst.observation_data  # .loc[ies_pst.nnz_obs_names,:]
            ies_obs = ies_obs.loc[ies_obs.obgnme.apply(lambda x: x in ies_pst.nnz_obs_groups), :]
            ies_obs.loc[:, "datetime"] = pd.to_datetime(ies_obs.obsnme.apply(lambda x: x.split('_')[-1]),
                                                        format='%Y%m%d')

            ies_pr_oe = pd.read_csv(os.path.join(ies_m_d, ies_case + ".0.obs.csv"))
            ies_pt_oe = pd.read_csv(os.path.join(ies_m_d, ies_case + ".{0}.obs.csv".format(noptmax)))

            da_pst = pyemu.Pst(os.path.join(da_m_d, da_case + ".pst"))
            da_obs = da_pst.observation_data.loc[da_pst.nnz_obs_names, :].copy()
            da_obs.loc[da_obs.obsnme.str.contains("gage"), "org_obgnme"] = "gage"

            num_cycles = 25
            da_pr_dict = {}
            da_pt_dict = {}
            for cycle in range(num_cycles):
                print(cycle)
                da_pr_oe = pd.read_csv(os.path.join(da_m_d, da_case + ".{0}.0.obs.csv".format(cycle)))
                da_pr_dict[cycle] = da_pr_oe
                pt_file = os.path.join(da_m_d, da_case + ".{0}.{1}.obs.csv".format(cycle, noptmax))
                if (os.path.exists(pt_file)):
                    da_pt_oe = pd.read_csv(pt_file)
                    da_pt_dict[cycle] = da_pt_oe
                else:
                    print("missing posterior", cycle)
            ies_og_uvals = ies_obs.obgnme.unique()
            ies_og_uvals.sort()

            for og in ies_og_uvals:
                ies_obs_og = ies_obs.loc[ies_obs.obgnme == og, :].copy()
                ies_obs_og.sort_values(by="datetime", inplace=True)
                da_obs_og = da_obs.loc[da_obs.org_obgnme == og, :]
                cmplx_obs_og = complex_obs.loc[complex_obs.org_obgnme == og, :].copy()
                cmplx_obs_og.sort_values(by="int_time", inplace=True)
                print(cmplx_obs_og.int_time)
                tms = pd.date_range(start='2015-12-31', periods=731, freq='D')
                dts = ies_obs_og.datetime.values

                def make_plot(axes):
                    ax = axes[0]
                    ax.set_title("batch DA, complex real {0}, observation location: ".format(i) + og, loc="left")
                    [ax.plot(dts, ies_pr_oe.loc[idx, ies_obs_og.obsnme], "0.5", alpha=0.5, lw=0.25) for idx in
                     ies_pr_oe.index]
                    [ax.plot(dts, ies_pt_oe.loc[idx, ies_obs_og.obsnme], "b", alpha=0.5, lw=0.5) for idx in
                     ies_pt_oe.index]
                    ax.plot(dts, ies_pr_oe.loc[ies_pr_oe.index[0], ies_obs_og.obsnme], "0.5", alpha=0.5,
                            lw=0.1, label="prior real")
                    ax.plot(dts, ies_pt_oe.loc[ies_pt_oe.index[0], ies_obs_og.obsnme], "b", alpha=0.5,
                            lw=0.1, label="post real")
                    # ax.plot(dts, ies_obs_og.obsval, "r", label="truth")
                    ax.plot(tms, cmplx_obs_og.iloc[:, i], 'r', label="truth")
                    ies_obs_nz = ies_obs_og.loc[ies_obs_og.weight > 0, :]

                    ax.scatter(ies_obs_nz.datetime.values, ies_obs_nz.obsval, marker="^", color="r", s=50,
                               zorder=10, label="obs")
                    ax.legend(loc="upper right")
                    ax = axes[1]

                    ax.set_title(
                        "sequential DA, complex real {0}, observation location: ".format(i) + da_obs_og.obsnme.values[
                            0],
                        loc="left")
                    # ax.plot(dts, ies_obs_og.obsval, "r", label="truth")
                    ax.plot(tms, cmplx_obs_og.iloc[:, i], 'r', label="truth")
                    ax.scatter(ies_obs_nz.datetime.values, ies_obs_nz.obsval, marker="^", color="r",
                               s=50, zorder=10, label="obs")

                    post_labeled = False
                    # for ccycle in range(cycle):
                    #     da_pr_oe = da_pr_dict[ccycle]
                    #     label = None
                    #     if ccycle == 0:
                    #         label = "prior real"
                    #     ax.scatter([dts[ccycle] for _ in range(da_pr_oe.shape[0])],
                    #                da_pr_oe.loc[:, da_obs_og.obsnme[0]].values,
                    #                marker=".", color="0.5", label=label)
                    #
                    #     if ccycle in da_pt_dict:
                    #         da_pt_oe = da_pt_dict[ccycle]
                    #         label = None
                    #         if not post_labeled:
                    #             label = "post real"
                    #             post_labeled = True
                    #         ax.scatter([dts[ccycle] for _ in range(da_pt_oe.shape[0])],
                    #                    da_pt_oe.loc[:, da_obs_og.obsnme[0]].values,
                    #                    marker=".", color="b", label=label)
                    ax.set_ylim(axes[0].get_ylim())
                    ax.legend(loc="upper right")

                    # prior

                fig, axes = plt.subplots(2, 1, figsize=(8, 8))
                make_plot(axes)
                for cycle in range(num_cycles):
                    da_pr_oe = da_pr_dict[cycle]
                    ax = axes[1]
                    ax.scatter([dts[cycle] for _ in range(da_pr_oe.shape[0])],
                               da_pr_oe.loc[:, da_obs_og.obsnme[0]].values,
                               marker=".", color="0.5")
                    ax.set_ylim(axes[0].get_ylim())

                # plt.tight_layout()
                # pdf.savefig()
                # plt.close(fig)

                # posterior
                # fig, axes = plt.subplots(2, 1, figsize=(8, 8))
                # make_plot(axes)
                for cycle in range(num_cycles):
                    da_pr_oe = da_pr_dict[cycle]
                    ax.scatter([dts[cycle] for _ in range(da_pr_oe.shape[0])],
                               da_pr_oe.loc[:, da_obs_og.obsnme[0]].values,
                               marker=".", color="0.5")
                    if cycle in da_pt_dict:
                        da_pt_oe = da_pt_dict[cycle]
                        ax.scatter([dts[cycle] for _ in range(da_pt_oe.shape[0])],
                                   da_pt_oe.loc[:, da_obs_og.obsnme[0]].values,
                                   marker=".", color="b")
                        ax.set_ylim(axes[0].get_ylim())
                        # plt.tight_layout()
                pdf.savefig()
                plt.close(fig)
                # plt.show()


def invest():
    csim = flopy.mf6.MFSimulation.load(sim_ws="complex_template")
    cm = csim.get_model("freyberg6")
    ssim = flopy.mf6.MFSimulation.load(sim_ws="simple_template_ies")
    sm = ssim.get_model("freyberg6")
    ctotim = np.cumsum(csim.tdis.perioddata.array["perlen"])
    stotim = np.cumsum(ssim.tdis.perioddata.array["perlen"])

    # for kper in range(csim.tdis.nper):
    #     carr = cm.wel.stress_period_data.array[kper]
    # print(carr)
    # return

    carr = cm.rch.recharge.array
    carr = carr.mean(axis=(1, 2, 3))
    sarr = sm.rch.recharge.array
    sarr = sarr.mean(axis=(1, 2, 3))
    print(sarr.shape, stotim.shape)
    print(carr.shape, ctotim.shape)

    fig, ax = plt.subplots(1, 1, figsize=(6, 6))
    ax.plot(stotim[1:], sarr[1:])
    ax.plot(ctotim[1:], carr[1:])
    plt.show()


def test_extract_state_obs(t_d):
    cwd = os.getcwd()
    os.chdir(t_d)
    fnames = extract_state_obs()
    os.chdir(cwd)
    return fnames


def extract_state_obs():
    hds = flopy.utils.HeadFile('freyberg6_freyberg.hds')
    arr = hds.get_data()
    fnames = []
    for k, a in enumerate(arr):
        fname = 'heads_' + str(k) + '.dat'
        np.savetxt(fname, a, fmt='%15.6E')
        fnames.append(fname)
    return fnames


def setup_interface(org_ws, num_reals=100):
    """copied from auto_pest.py

    """
    np.random.seed(123456)

    # run mf6
    tmp_ws = org_ws + "_temp"
    if os.path.exists(tmp_ws):
        shutil.rmtree(tmp_ws)
    shutil.copytree(org_ws, tmp_ws)
    pyemu.os_utils.run("mf6", cwd=tmp_ws)

    # load the mf6 model with flopy to get the spatial reference
    sim = flopy.mf6.MFSimulation.load(sim_ws=tmp_ws)
    m = sim.get_model("freyberg6")

    # work out the spatial rediscretization factor
    redis_fac = m.dis.nrow.data / 40

    # where the pest interface will be constructed
    template_ws = org_ws + "_template"

    # instantiate PstFrom object
    pf = pyemu.utils.PstFrom(original_d=tmp_ws, new_d=template_ws,
                             remove_existing=True,
                             longnames=True, spatial_reference=m.modelgrid,
                             zero_based=False, start_datetime="1-1-2018")

    # add observations from the sfr observation output file
    df = pd.read_csv(os.path.join(template_ws, "sfr.csv"), index_col=0)
    pf.add_observations("sfr.csv", insfile="sfr.csv.ins", index_cols="time",
                        use_cols=list(df.columns.values),
                        prefix="sfr")

    # add observations for the heads observation output file
    df = pd.read_csv(os.path.join(template_ws, "heads.csv"), index_col=0)
    pf.add_observations("heads.csv", insfile="heads.csv.ins",
                        index_cols="time", use_cols=list(df.columns.values),
                        prefix="hds")

    # add observations for simulated states
    pf.add_py_function("workflow.py", "extract_state_obs()", is_pre_cmd=False)
    fnames = test_extract_state_obs(template_ws)
    for k, fname in enumerate(fnames):
        prefix = "head_k:{0}".format(k)
        pf.add_observations(fname, prefix=prefix, obsgp=prefix)

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
    temporal_gs = pyemu.geostats.GeoStruct(variograms=pyemu.geostats.ExpVario(contribution=1.0, a=60))

    # import flopy as part of the forward run process
    pf.extra_py_imports.append('flopy')

    # use the idomain array for masking parameter locations
    ib = m.dis.idomain[0].array

    # define a dict that contains file name tags and lower/upper bound information
    tags = {"npf_k_": [0.1, 10.], "npf_k33_": [.1, 10], "sto_ss": [.1, 10], "sto_sy": [.9, 1.1],
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
        for arr_file in arr_files:
            arr = np.loadtxt(os.path.join(template_ws, arr_file)).reshape(ib.shape)
            np.savetxt(os.path.join(template_ws, arr_file), arr, fmt="%15.6E")

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
                kper = int(arr_file.split('.')[1].split('_')[-1]) - 1
                pf.add_parameters(filenames=arr_file, par_type="constant", par_name_base=arr_file.split('.')[1] + "_cn",
                                  pargp="rch_const", zone_array=ib, upper_bound=uub, lower_bound=llb,
                                  geostruct=temporal_gs,
                                  datetime=dts[kper],
                                  par_style="direct")


        # otherwise...
        else:
            # for each array add both grid-scale and pilot-point scale parameters
            for arr_file in arr_files:
                pf.add_parameters(filenames=arr_file, par_type="grid", par_name_base=arr_file.split('.')[1] + "_gr",
                                  pargp=arr_file.split('.')[1] + "_gr", zone_array=ib, upper_bound=ub, lower_bound=lb,
                                  geostruct=grid_gs)
                pf.add_parameters(filenames=arr_file, par_type="pilotpoints",
                                  par_name_base=arr_file.split('.')[1] + "_pp",
                                  pargp=arr_file.split('.')[1] + "_pp", zone_array=ib, upper_bound=ub, lower_bound=lb,
                                  pp_space=int(5 * redis_fac), geostruct=pp_gs)

    # add direct pars for the ic strt values
    tag = "ic_strt"
    arr_files = [f for f in os.listdir(template_ws) if tag in f and f.endswith(".txt")]
    for arr_file in arr_files:
        # make sure each array file in nrow X ncol dimensions (not wrapped)
        arr = np.loadtxt(os.path.join(template_ws, arr_file)).reshape(ib.shape)
        np.savetxt(os.path.join(template_ws, arr_file), arr, fmt="%15.6E")
        k = int(arr_file.split('.')[1][-1]) - 1
        prefix = "head_k:{0}".format(k)
        zn_arr = np.ones_like(arr, dtype=int)
        zn_arr[arr < 0] = 0
        zn_arr[arr > 1000] = 0
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
                          upper_bound=1.5, lower_bound=0.5, datetime=dts[kper], geostruct=temporal_gs)

    # add temporally indep, but spatially correlated grid-scale parameters, one per well
    pf.add_parameters(filenames=list_files, par_type="grid", par_name_base="wel_grid",
                      pargp="wel_grid", index_cols=[0, 1, 2], use_cols=[3],
                      upper_bound=1.5, lower_bound=0.5)

    # add grid-scale parameters for SFR reach conductance.  Use layer, row, col and reach
    # number in the parameter names
    pf.add_parameters(filenames="freyberg6.sfr_packagedata.txt", par_name_base="sfr_rhk",
                      pargp="sfr_rhk", index_cols=[0, 1, 2, 3], use_cols=[9], upper_bound=10.,
                      lower_bound=0.1,
                      par_type="grid")

    # add model run command
    pf.mod_sys_cmds.append("mf6")

    # build pest control file
    pst = pf.build_pst('freyberg.pst')
    par = pst.parameter_data
    strt_pars = par.loc[par.parnme.str.contains("head_k"), "parnme"]

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
    state_par = par.loc[par.parnme.str.contains("direct_head"), :].copy()
    for v in ["k", "i", "j"]:
        state_par.loc[:, v] = state_par.loc[:, v].apply(int)
        state_obs.loc[:, v] = state_obs.loc[:, v].apply(int)
    state_par_dict = {"{0}_{1}_{2}".format(k, i, j): n for k, i, j, n in
                      zip(state_par.k, state_par.i, state_par.j, state_par.parnme)}
    obs.loc[:, "state_par_link"] = np.nan
    state_obs.loc[:, "kij"] = state_obs.apply(lambda x: "{0}_{1}_{2}".format(x.k, x.i, x.j), axis=1)

    # for kij,n in state_par_dict.items():
    #    if kij not in state_obs.kij:
    #        print(kij,n)
    obs.loc[state_obs.obsnme, "state_par_link"] = state_obs.apply(lambda x: state_par_dict.get((x.kij), np.nan), axis=1)
    print(obs.state_par_link.dropna().shape)

    keep = ['arrobs_head_k:0_i:22_j:15', 'arrobs_head_k:2_i:2_j:9', 'arrobs_head_k:2_i:33_j:7', 'gage_1']

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


def monthly_ies_to_da(org_d):
    """convert the batch monthly model to sequential formulation"""

    # load the existing 25-stress period model
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
        f.write("BEGIN Options\n  TIME_UNITS  days\nEND Options\n")
        f.write("BEGIN Dimensions\n  NPER  1\nEND Dimensions\n")
        f.write("BEGIN PERIODDATA\n31.00000000  1       1.00000000\nEND PERIODDATA\n")
    # make sure it runs
    pyemu.os_utils.run("mf6", cwd=t_d)

    # write a tdis template file
    with open(os.path.join(t_d, "freyberg6.tdis.tpl"), 'w') as f:
        f.write("ptf  ~\n")
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
        df = pst.add_parameters(tpl, inf, pst_path=".")
        pst.parameter_data.loc[df.parnme, "cycle"] = cy
        pst.parameter_data.loc[df.parnme, "partrans"] = "fixed"

    # write par cycle table for the perlen par (fixed)
    pers = org_sim.tdis.perioddata.array["perlen"]
    pdf = pd.DataFrame(index=['perlen'], columns=np.arange(25))
    pdf.loc['perlen', :] = pers
    pdf.to_csv(os.path.join(t_d, "par_cycle_table.csv"))
    pst.pestpp_options["da_parameter_cycle_table"] = "par_cycle_table.csv"

    # save this for later!
    org_obs = pst.observation_data.copy()

    # now drop all existing heads obs package observations since
    # those will be replaced by the state obs
    fname = 'heads.csv'
    fname_ins = fname + ".ins"
    pst.drop_observations(os.path.join(t_d, fname_ins), '.')

    # fix the sfr instruction file - only need one output time
    fname = 'sfr.csv'
    fname_ins = fname + ".ins"
    pst.drop_observations(os.path.join(t_d, fname_ins), '.')

    ins_lines = open(os.path.join(t_d, fname_ins), 'r').readlines()
    with open(os.path.join(t_d, fname_ins), 'w') as f:
        for line in ins_lines[:3]:
            f.write(line)

    new_ins_files = [os.path.join(t_d, fname_ins)]
    new_out = [os.path.join(t_d, "sfr.csv")]
    new_ins_cycle = [-1]

    # add sfr obs
    for ins_file in new_ins_files:
        pst.add_observations(ins_file, pst_path=".")

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
            pst.model_input_data.iloc[i, 1] = pst.model_input_data.iloc[i, 1].replace(str(cy), "1")
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


def run_batch_seq_prior_monte_carlo():
    """run prior monte carlo for the batch and seq monthly models
    """
    t_d = "seq_monthly_model_files_template"
    pst = pyemu.Pst(os.path.join(t_d,"freyberg.pst"))
    pst.control_data.noptmax = -1
    pst.write(os.path.join(t_d,"freyberg.pst"),version=2)
    pyemu.os_utils.start_workers(t_d, "pestpp-da", "freyberg.pst", num_workers=10,
                                 master_dir=t_d.replace("template", "master_prior"))

    t_d = "monthly_model_files_template"
    pst = pyemu.Pst(os.path.join(t_d, "freyberg.pst"))
    pst.control_data.noptmax = -1
    pst.write(os.path.join(t_d, "freyberg.pst"),version=2)
    pyemu.os_utils.start_workers(t_d, "pestpp-ies", "freyberg.pst", num_workers=10,
                                 master_dir=t_d.replace("template", "master_prior"))



def plot_prior_mc():
    """plot the prior monte carlo results for daily, monthly batch and monthly sequential

    """
    c_m_d = "daily_model_files_master_prior"
    s_b_m_d = "monthly_model_files_master_prior"
    s_s_m_d = "seq_monthly_model_files_master_prior"

    c_pst = pyemu.Pst(os.path.join(c_m_d, "freyberg.pst"))
    obs = c_pst.observation_data
    cobs = obs.loc[obs.obsnme.str.startswith("hds_usecol:arrobs_head_"),:]
    cobs.loc[:,"time"] = cobs.time.apply(float)

    s_b_pst = pyemu.Pst(os.path.join(s_b_m_d, "freyberg.pst"))
    c_oe = pd.read_csv(os.path.join(c_m_d, "freyberg.0.obs.csv"), index_col=0)
    s_b_oe = pd.read_csv(os.path.join(s_b_m_d, "freyberg.0.obs.csv"), index_col=0)

    s_s_pst = pyemu.Pst(os.path.join(s_s_m_d,"freyberg.pst"))
    seq_oe_files = [f for f in os.listdir(s_s_m_d) if f.endswith(".oe.csv") and "global" in f and f.startswith("freyberg")]
    s_s_oe_dict = {int(f.split(".")[2]):pd.read_csv(os.path.join(s_s_m_d,f),index_col=0) for f in seq_oe_files}
    #oct = pd.read_csv(os.path.join(s_s_m_d,"obs_cycle_tbl.csv"),index_col=0)

    ognames = list(cobs.obgnme.unique())
    ognames.sort()
    ognames.append("sfr_usecol:gage_1")
    with PdfPages("prior_obs_v_sim.pdf") as pdf:

        for ogname in ognames:
            cgobs = cobs.loc[cobs.obgnme==ogname,:].copy()
            sgobs = s_b_pst.observation_data.loc[s_b_pst.observation_data.obgnme==ogname,:].copy()
            sgobs.loc[:,"time"] = sgobs.time.apply(float)
            cgobs.loc[:, "time"] = cgobs.time.apply(float)

            sgobs.sort_values(by="time", inplace=True)
            cgobs.sort_values(by="time",inplace=True)

            fig,ax = plt.subplots(1,1,figsize=(8,8))

            [ax.plot(cgobs.time, c_oe.loc[idx, cgobs.obsnme], "b", lw=0.01,alpha=0.5) for idx in c_oe.index]

            [ax.plot(sgobs.time,s_b_oe.loc[idx,sgobs.obsnme],"0.5",lw=0.01,alpha=0.5) for idx in s_b_oe.index]

            if "hds" in ogname:
                seq_name = ogname.replace("hds_usecol:","")
            else:
                seq_name = ogname + "_time:10000.0"
            for itime,time in enumerate(sgobs.time):
                oe = s_s_oe_dict[itime]
                #print(oe.loc[:,seq_name])
                ax.scatter([time for _ in range(oe.shape[0])], oe.loc[:, seq_name], marker=".",color="0.5",alpha=0.5)

            ax.set_title(ogname)
            if "gage" not in ogname:
                ax.set_ylim(30,ax.get_ylim()[1])
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
    coe = pd.read_csv(os.path.join(c_d,"freyberg.0.obs.csv"),index_col=0)
    cvals = coe.loc[real_idx,:]

    # load the batch monthly model control file
    bpst = pyemu.Pst(os.path.join(b_d,"freyberg.pst"))
    obs = bpst.observation_data
    # assign all common observations
    obs.loc[:,"obsval"] = cvals.loc[bpst.obs_names]

    # set some weights - only the first 12 months (not counting spin up time)
    # just picked weights that seemed to balance-ish the one realization I looked at...
    obs.loc[:,"weight"] = 0.0
    for k in keep:
        kobs = obs.loc[obs.obsnme.str.contains(k),:].copy()
        kobs.loc[:,"time"] = kobs.time.apply(float)
        kobs.sort_values(by="time",inplace=True)
        if "gage" in k:
            obs.loc[kobs.obsnme[1:13], "weight"] = 0.005
        else:
            obs.loc[kobs.obsnme[1:13],"weight"] = 2.0

    assert bpst.nnz_obs == 48
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

    # load the batch control file
    bpst = pyemu.Pst(os.path.join(b_d, "freyberg.pst"))
    bobs = bpst.observation_data

    # work on the sequential obs names - yuck
    seq_names = [o.split("_time")[0].replace("hds_usecol:","") if "hds_usecol" in o else o for o in bpst.nnz_obs_names ]
    seq_names = set([o for o in seq_names if "arrobs" in o or "time:10000.0" in o ])

    # load the sequential control file
    spst = pyemu.Pst(os.path.join(s_d, "freyberg.pst"))
    spst.observation_data.loc[:,"weight"] = 0.0
    snames = set(spst.obs_names)
    assert len(seq_names - snames) == 0
    seq_names = list(seq_names)
    seq_names.sort()

    sobs = spst.observation_data.loc[seq_names]

    # build the weight and obs cycle tables
    odf = pd.DataFrame(columns = np.arange(25),index=seq_names)
    wdf = pd.DataFrame(columns=np.arange(25), index=seq_names)
    wdf.loc[:,:] = 0.0
    for seq_name in seq_names:

        # the batch obs info for this sequential obs name
        bsobs = bobs.loc[bobs.obsnme.str.contains(seq_name),:].copy()
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


if __name__ == "__main__":

    #setup_interface("monthly_model_files")
    #monthly_ies_to_da("monthly_model_files_template")
    b_d = map_complex_to_simple_bat("daily_model_files_master_prior","monthly_model_files_template",1)
    s_d = map_simple_bat_to_seq(b_d,"seq_monthly_model_files_template")
    #process_complex_target_output('daily_model_files_master_prior','monthly_model_files_template','seq_monthly_model_files_template',1 )
    #run_batch_seq_prior_monte_carlo()
    #setup_interface("daily_model_files")
    #run_complex_prior_mc('daily_model_files_template')
    #plot_prior_mc()

    #compare_mf6_freyberg()

    exit()

    # invest()
    # exit()

    # BOOLEANS TO SELECT CODE BLOCKS BELOW
    prep_complex_model = False  # do this once before running paired simple/complex analysis
    run_prior_mc = False
    run_simple_complex = False
    clean_dirs = False
    plot_s_vs_s = False
    plot_phis = False
    plot_phi_diffs = False
    plot_obs_sim = False

    if prep_complex_model:
        prep_complex_prior_mc()

    if run_prior_mc:
        run_complex_prior_mc('complex_template')

    if run_simple_complex:
        compare_mf6_freyberg()

    if plot_phis:
        plot_phi_seq_bat()

    if plot_phi_diffs:
        plot_phi_diff_seq_bat()

    if plot_s_vs_s:
        s_plot()

    if plot_obs_sim:
        plots_obs_v_sim()

    if clean_dirs:
        clean_master_dirs()
