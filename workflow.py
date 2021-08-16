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


def mod_tdis_sto(org_t_d, t_d):
    tdis_file = "freyberg6.tdis"
    lines = open(os.path.join(org_t_d, tdis_file), 'r').readlines()

    with open(os.path.join(t_d, tdis_file), 'w') as f:
        iline = 0
        while True:
            line = lines[iline]
            if "begin period" in line.lower():
                lines[iline + 1] = "100000  1   1.0\n"
                print(lines[iline + 1])
            print(line)
            f.write(line)
            iline += 1
            if iline >= len(lines):
                break
    sto_file = "freyberg6.sto"
    lines = open(os.path.join(org_t_d, sto_file), 'r').readlines()

    with open(os.path.join(t_d, sto_file), 'w') as f:

        for line in lines:
            f.write(line)
            if line.lower().startswith("end griddata"):
                break
        f.write("\nbegin period 1\n  transient\nend period 1\n")


def da_prep_4_mf6_freyberg_seq(sync_state_names=True):
    t_d = "simple_template_da"
    if os.path.exists(t_d):
        shutil.rmtree(t_d)
    shutil.copytree("monthly_template", t_d)
    for f in os.listdir(t_d):
        for tag in ["ies", "opt", "glm"]:
            if tag in f.lower():
                os.remove(os.path.join(t_d, f))

    # first modify the tdis
    with open(os.path.join(t_d, "freyberg6.tdis"), 'w') as f:
        f.write("BEGIN Options\n  TIME_UNITS  days\nEND Options\n")
        f.write("BEGIN Dimensions\n  NPER  1\nEND Dimensions\n")
        f.write("BEGIN PERIODDATA\n31.00000000  1       1.00000000\nEND PERIODDATA\n")
    # make sure it runs
    pyemu.os_utils.run("mf6", cwd=t_d)

    # write a tdis template file - could possibly keep all 25 stress periods to
    # simulate a 2-year-ahead forecast...

    with open(os.path.join(t_d, "freyberg6.tdis.tpl"), 'w') as f:
        f.write("ptf  ~\n")
        f.write("BEGIN Options\n  TIME_UNITS  days\nEND Options\n")
        f.write("BEGIN Dimensions\n  NPER  1\nEND Dimensions\n")
        f.write("BEGIN PERIODDATA\n~  perlen  ~  1       1.00000000\nEND PERIODDATA\n")
    new_tpl, new_in = [os.path.join(t_d, "freyberg6.tdis.tpl")], [os.path.join(t_d, "freyberg6.tdis")]
    new_tpl_cycle = [-1]

    # split out the head, sfr and list instruction files into multiple instruction file
    # eventually, want to move to da par and obs cycle tables for heads and gage_1 obs
    # lines = open(os.path.join(t_d,"heads.csv.ins"),'r').readlines()[2:]
    # new_ins,new_out,new_ins_cycle = [],[],[]
    # #print(lines)
    # for icycle, line in enumerate(lines):
    #     ins_name = os.path.join(t_d,"heads_{0}.csv.ins".format(icycle))
    #     with open(ins_name,'w') as f:
    #         f.write("pif ~\nl1\n")
    #         f.write(line)
    #     new_ins.append(ins_name)
    #     new_out.append(os.path.join(t_d,"heads.csv"))
    #     new_ins_cycle.append(icycle)
    # remove_ins = ["heads.csv.ins"]
    hds = pd.read('heads.csv')
    lay = hds.columns.to_series().apply(lambda x: x.split('_')[1])
    row = hds.columns.to_series().apply(lambda x: x.split('_')[2])
    col = hds.columns.to_series().apply(lambda x: x.split('_')[3])
    fname = 'heads.csv'
    fname_ins = fname + ".ins"
    with open(fname_ins, 'w') as f:
        f.write("pif ~\n")
        for i in range(hds.shape[0]):
            f.write("l1 \n")
            f.write("l1")
            oname = "head_k:{0}_i:{1}_j:{2}".format(lay[i], row[i], col[i])
            f.write("w !{0}! ".format(oname))
        new_ins.append(ins_name)
        new_out.append(os.path.join(t_d, "heads.csv"))
        new_ins_cycle.append(i)
    remove_ins.append("heads.csv.ins")

    lines = open(os.path.join(t_d, "sfr.csv.ins"), 'r').readlines()[2:]
    # print(lines)
    for icycle, line in enumerate(lines):
        ins_name = os.path.join(t_d, "sfr_{0}.csv.ins".format(icycle))
        with open(ins_name, 'w') as f:
            f.write("pif ~\nl1\n")
            f.write(line)
        new_ins.append(ins_name)
        new_out.append(os.path.join(t_d, "sfr.csv"))
        new_ins_cycle.append(icycle)
    remove_ins.append("sfr.csv.ins")

    lines = open(os.path.join(t_d, "freyberg6.lst.ins"), 'r').readlines()[1:]
    icycle = 0
    tag_line = lines[0]
    for s in range(0, len(lines), 13):
        ins_name = os.path.join(t_d, "freyberg6_{0}.lst.ins".format(icycle))
        with open(os.path.join(t_d, "freyberg6_{0}.lst.ins".format(icycle)), 'w') as f:
            f.write("pif ~\n")
            f.write(tag_line)
            for line in lines[s + 1:s + 13]:
                f.write(line)
        new_ins.append(ins_name)
        new_out.append(os.path.join(t_d, "freyberg6.lst"))
        new_ins_cycle.append(icycle)
        icycle += 1
    remove_ins.append("freyberg6.lst.ins")

    # modify the ic file
    k = 0
    with open(os.path.join(t_d, "freyberg6.ic"), 'r') as f:
        while True:
            line = f.readline()
            if line == "":
                break
            if line.lower().strip().startswith("internal"):
                arr_lines = []
                while True:
                    line = f.readline()
                    if line == "":
                        raise Exception
                    if line.lower().strip().startswith("end"):
                        break
                    if line.lower().strip().startswith("internal"):
                        arr = np.array(arr_lines, dtype=np.float)
                        np.savetxt(os.path.join(t_d, "heads_{0}.dat_in".format(k)), arr, fmt="%15.6E")
                        k += 1
                        arr_lines = []
                    else:
                        arr_lines.append(line.strip().split())

        arr = np.array(arr_lines, dtype=np.float)
        np.savetxt(os.path.join(t_d, "heads_{0}.dat_in".format(k)), arr, fmt="%15.6E")
    with open(os.path.join(t_d, "freyberg6.ic"), 'w') as f:
        f.write("begin griddata\nstrt layered\n")
        for k in range(3):
            f.write("open/close 'heads_{0}.dat_in' FACTOR 1.0\n".format(k))
        f.write("end griddata\n")

    # write a python script to extract final heads and save to files
    with open(os.path.join(t_d, "forward_run.py"), 'w') as f:
        f.write("import numpy as np\nimport flopy\nimport pyemu\n")
        f.write("pyemu.os_utils.run('mf6')\n")
        f.write("hds = flopy.utils.HeadFile('freyberg6_freyberg.hds')\n")
        f.write("arr = hds.get_data()\n")
        f.write("for k,a in enumerate(arr):\n")
        f.write("    np.savetxt('heads_'+str(k)+'.dat',a,fmt='%15.6E')\n")

    # dont run it so we can harvest the ic values in the arrays for setting the parval1 values
    pyemu.os_utils.run("python forward_run.py", cwd=t_d)

    # now write ins and tpl file for these
    ic_parvals = {}
    obs_to_par_map = dict()
    for k in range(3):
        fname = os.path.join(t_d, "heads_{0}.dat".format(k))
        assert os.path.exists(fname), fname
        arr = np.loadtxt(fname)
        fname_ins = fname + "_out.ins"
        fname_tpl = fname + "_in.tpl"
        in_arr = np.loadtxt(fname_tpl.replace(".tpl", ""))
        ft = open(fname_tpl, 'w')
        with open(fname_ins, 'w') as f:
            f.write("pif ~\n")
            ft.write("ptf ~\n")
            for i in range(arr.shape[0]):
                f.write("l1 ")
                for j in range(arr.shape[1]):
                    if np.abs(arr[i, j]) > 100 or np.abs(in_arr[i, j]) > 100:
                        f.write(" !dum! ")
                        ft.write(" 40 ")
                    else:
                        oname = "head_{0:02d}_{1:03d}_{2:03d}".format(k, i, j)
                        f.write(" !{0}! ".format(oname))
                        if sync_state_names:
                            ft.write(" ~  {0} ~ ".format(oname))
                            ic_parvals[oname] = in_arr[i, j]
                        else:
                            pname = "p" + oname
                            ft.write(" ~  {0} ~ ".format(pname))
                            obs_to_par_map[oname] = pname
                            ic_parvals[pname] = in_arr[i, j]

                f.write("\n")
                ft.write("\n")
        ft.close()
        new_tpl.append(fname_tpl)
        new_in.append(fname_tpl.replace(".tpl", ""))
        new_tpl_cycle.append(-1)
        new_ins.append(fname_ins)
        new_out.append(os.path.join(t_d, "heads_{0}.dat".format(k)))
        new_ins_cycle.append(-1)

        i = pyemu.pst_utils.InstructionFile(fname_ins)
        df = i.read_output_file(fname)
        # print(df)

    # split out the wel and rch tpl files into cycle files

    # lines = []
    # with open(os.path.join(t_d,"freyberg6.wel.tpl"),'r') as f:
    #     for i in range(19):
    #         lines.append(f.readline())
    # print(lines)
    #
    for icycle in range(25):
        tpl_file = os.path.join(t_d, "wel_grid_{0}_inst0_grid.csv.tpl".format(icycle))
        # with open(tpl_file,'w') as f:
        #     for line in lines:
        #         new_line = line.replace("_0","_{0}".format(icycle))
        #         f.write(new_line)

        new_tpl.append(tpl_file)
        new_in.append(os.path.join(t_d, "freyberg6.wel"))
        new_tpl_cycle.append(icycle)

    for icycle in range(25):
        tpl_file = os.path.join(t_d, "twel_mlt_{0}_inst0_constant.csv.tpl".format(icycle))
        # with open(tpl_file,'w') as f:
        #     for line in lines:
        #         new_line = line.replace("_0","_{0}".format(icycle))
        #         f.write(new_line)

        new_tpl.append(tpl_file)
        new_in.append(os.path.join(t_d, "freyberg6.wel"))
        new_tpl_cycle.append(icycle)
    # remove_tpl = ["freyberg6.wel.tpl"]
    #
    # lines = []
    # with open(os.path.join(t_d, "freyberg6.rch.tpl"), 'r') as f:
    #     for i in range(11):
    #         lines.append(f.readline())
    # print(lines)
    #
    for icycle in range(25):
        tpl_file = os.path.join(t_d, "rch_recharge_{0}_cn_inst0_constant.csv.tpl".format(icycle))
        #     with open(tpl_file,'w') as f:
        #         for line in lines:
        #             new_line = line.replace("_0","_{0}".format(icycle))
        #             f.write(new_line)
        #
        new_tpl.append(tpl_file)
        new_in.append(os.path.join(t_d, "freyberg6.rch"))
        new_tpl_cycle.append(icycle)
    # remove_tpl.append('freyberg6.rch.tpl')

    # now for the fun part: modify the pst
    pst = pyemu.Pst(os.path.join(t_d, "freyberg.pst"))
    pst.parameter_data = pst.parameter_data.drop(
        ['extra', 'inst', 'i', 'j', 'x', 'y', 'zone', 'usecol', 'idx0', 'idx1', 'idx2', 'idx3'], axis=1)
    pst.model_input_data.loc[:, "cycle"] = -1
    pst.parameter_data.loc[:, "cycle"] = -1

    for i in range(9766):
        if pst.parameter_data.iloc[i, 0].startswith('wel_grid'):
            cy = int(pst.parameter_data.iloc[i, 0].split('_')[2])
            pst.parameter_data.iloc[i, 10] = cy
        elif pst.parameter_data.iloc[i, 0].startswith('twel_mlt'):
            pst.parameter_data.iloc[i, 10] = cy
        elif pst.parameter_data.iloc[i, 0].startswith('multiplier_const_rch'):
            cy = int(pst.parameter_data.iloc[i, 0].split('_')[4])
            pst.parameter_data.iloc[i, 10] = cy

    # print(pst.npar_adj,pst.nnz_obs)
    #
    # # swap out obs info
    # dropped_dfs = []
    # # for ins in remove_ins:
    # #     dropped_dfs.append(pst.drop_observations(os.path.join(t_d, ins), '.'))
    # for insf, outf, cy in zip(new_ins, new_out, new_ins_cycle):
    #     df = pst.add_observations(insf,outf, pst_path=".")
    #     pst.observation_data.loc[:, "cycle"] = cy
    #     pst.model_output_data.loc[os.path.join(".",os.path.split(insf)[1]),"cycle"] = cy
    # pst.observation_data.loc[:,"weight"] = 0.0
    # for df in dropped_dfs:
    #     for c in ["obsval","weight"]:
    #         pst.observation_data.loc[df.obsnme, c] = df.loc[:, c]
    #
    # # swap out par info
    # dropped_dfs = []
    # # for tpl in remove_tpl:
    # #     dropped_dfs.append(pst.drop_parameters(os.path.join(t_d,tpl),'.'))
    # for tplf, inf, cy in zip(new_tpl,new_in,new_tpl_cycle):
    #     df = pst.add_parameters(tplf,inf,pst_path=".")
    #     pst.parameter_data.loc[df.parnme,"cycle"] = cy
    #     pst.model_input_data.loc[os.path.join(".",os.path.split(tplf)[1]),"cycle"] = cy
    # for df in dropped_dfs:
    #     for c in ["parval1","parubnd","parlbnd","pargp"]:
    #         pst.parameter_data.loc[df.parnme,c] = df.loc[:,c]

    # set the state parameter info
    for p, v in ic_parvals.items():
        pst.parameter_data.loc[p, "parval1"] = v
        pst.parameter_data.loc[p, "parlbnd"] = v * 0.5
        pst.parameter_data.loc[p, "parubnd"] = v * 1.5
        pst.parameter_data.loc[p, "pargp"] = "head_state"
        pst.parameter_data.loc[p, "partrans"] = "none"
        pst.parameter_data.loc[p, "scale"] = 1.0
        pst.parameter_data.loc[p, "offset"] = 0.0
        pst.parameter_data.loc[p, "dercom"] = 1
        pst.parameter_data.loc[p, "cycle"] = -1
        pst.parameter_data.loc[p, "parchglim"] = "factor"
        pst.parameter_data.loc[p, "parnme"] = p

    # add all head locs to obs data
    pst.observation_data = pst.observation_data.drop(['extra', 'usecol', 'time'], axis=1)
    for p, v in ic_parvals.items():
        pst.observation_data.loc[p, "obsval"] = v
        pst.observation_data.loc[p, "obsnme"] = 'o' + p
        pst.observation_data.loc[p, "weight"] = 0.
        pst.observation_data.loc[p, "obgnme"] = 'obgnme'
        pst.observation_data.loc[p, "state_par_link"] = p

    # set up per length parameter
    pst.parameter_data.loc["perlen", "partrans"] = "fixed"
    pst.parameter_data.loc["perlen", "parval1"] = 1.0
    pst.parameter_data.loc["perlen", "parlbnd"] = 1e-10
    pst.parameter_data.loc["perlen", "parubnd"] = 11000000000.0
    pst.parameter_data.loc["perlen", "pargp"] = "pargp"
    pst.parameter_data.loc["perlen", "parchglim"] = "factor"
    pst.parameter_data.loc["perlen", "scale"] = 1.0
    pst.parameter_data.loc["perlen", "offset"] = 0
    pst.parameter_data.loc["perlen", "dercom"] = 1
    pst.parameter_data.loc["perlen", "cycle"] = -1
    pst.parameter_data.loc["perlen", "parnme"] = "perlen"

    pst.control_data.noptmax = 3
    pst.model_command = "python forward_run.py"
    pst.pestpp_options.pop("ies_par_en")
    pst.pestpp_options["ies_num_reals"] = 50
    pst.pestpp_options["da_num_reals"] = 50
    # if not sync_state_names:
    #     pst.observation_data.loc[:,"state_par_link"] = np.NaN
    #     obs = pst.observation_data
    #     obs.loc[:,"state_par_link"] = obs.obsnme.apply(lambda x: obs_to_par_map.get(x,np.NaN))
    pst.write(os.path.join(t_d, "freyberg6_run_da1.pst"), version=2)
    return pst


def da_prep_4_mf6_freyberg_seq_tbl():
    t_d = "simple_template_da"
    if not os.path.exists(t_d):
        da_prep_4_mf6_freyberg_seq()
    pst = pyemu.Pst(os.path.join(t_d, "freyberg6_run_da1.pst"))
    pdf = pd.DataFrame({"perlen": 31}, index=np.arange(25))
    pdf.T.to_csv(os.path.join(t_d, "par_cycle_tbl.csv"))
    pst.pestpp_options["da_parameter_cycle_table"] = "par_cycle_tbl.csv"

    # mod sfr_out
    sfr_ins_file = os.path.join(t_d, "sfr.csv.ins")
    with open(sfr_ins_file, 'w') as f:
        f.write("pif ~\n")
        f.write("l1\n")
        f.write("l1 ~,~ !headwater!  ~,~ !tailwater!  ~,~ !gage_1!\n")

    # and lst budget
    lines = open(os.path.join(t_d, "freyberg6_0.lst.ins"), 'r').readlines()
    lst_ins_file = os.path.join(t_d, "freyberg6.lst.ins")
    with open(lst_ins_file, 'w') as f:
        for line in lines:
            f.write(line.replace("_20151231", ""))

    obs = pst.observation_data.copy()
    tr_obs = obs.loc[obs.obsnme.str.startswith("trgw"), :].copy()
    tr_obs.loc[tr_obs.obsnme, "datetime"] = pd.to_datetime(tr_obs.obsnme.apply(lambda x: x.split('_')[-1]))
    tr_obs.loc[tr_obs.obsnme, "k"] = tr_obs.obsnme.apply(lambda x: np.int(x.split('_')[1]))
    tr_obs.loc[tr_obs.obsnme, "i"] = tr_obs.obsnme.apply(lambda x: np.int(x.split('_')[2]))
    tr_obs.loc[tr_obs.obsnme, "j"] = tr_obs.obsnme.apply(lambda x: np.int(x.split('_')[3]))
    tr_obs.loc[tr_obs.obsnme, "obgnme"] = tr_obs.obsnme.apply(lambda x: "_".join(x.split("_")[:-1]))

    head_obs = obs.loc[obs.obsnme.str.startswith("head_"), :].copy()
    head_obs.loc[head_obs.obsnme, "k"] = head_obs.obsnme.apply(lambda x: np.int(x.split('_')[1]))
    head_obs.loc[head_obs.obsnme, "i"] = head_obs.obsnme.apply(lambda x: np.int(x.split('_')[2]))
    head_obs.loc[head_obs.obsnme, "j"] = head_obs.obsnme.apply(lambda x: np.int(x.split('_')[3]))

    # print(pst.nobs)
    for ins_file in pst.model_output_data.pest_file.copy():
        if "heads_" in ins_file and ins_file.endswith("dat_out.ins"):
            continue
        pst.drop_observations(os.path.join(t_d, ins_file), pst_path=".")
    for ins_file in [sfr_ins_file, lst_ins_file]:
        pst.add_observations(ins_file, pst_path=".")

    # work out which heads are obs site
    obs_heads = {}
    odf_names = []
    pst.observation_data.loc[:, "org_obgnme"] = np.NaN
    pst.observation_data.loc[:, "weight"] = 0.0
    for og in tr_obs.obgnme.unique():
        site_obs = tr_obs.loc[tr_obs.obgnme == og, :]
        site_obs.sort_values(by="datetime", inplace=True)
        head_name = "head_{0:02d}_{1:03d}_{2:03d}".format(site_obs.k[0], site_obs.i[0], site_obs.j[0])
        for i, oname in enumerate(site_obs.obsnme):
            obs_heads[oname] = (head_name, i)
        # assign max weight in the control file since some have zero weight and
        # we are covering weights in the weight table
        # print(og,site_obs.weight.max())
        pst.observation_data.loc[head_name, "weight"] = site_obs.weight.max()
        pst.observation_data.loc[head_name, "org_obgnme"] = og
        odf_names.append(head_name)

    odf_names.append("gage_1")

    odf = pd.DataFrame(columns=odf_names, index=np.arange(25))
    wdf = pd.DataFrame(columns=odf_names, index=np.arange(25))
    for tr_name, (head_name, icycle) in obs_heads.items():
        odf.loc[icycle, head_name] = obs.loc[tr_name, "obsval"]
        wdf.loc[icycle, head_name] = obs.loc[tr_name, "weight"]

    g_obs = obs.loc[obs.obsnme.str.startswith("gage_1"), :].copy()
    # give these obs the max weight since some have zero weight
    pst.observation_data.loc["gage_1", "weight"] = g_obs.weight.max()
    g_obs.sort_index(inplace=True)
    for i, name in enumerate(g_obs.obsnme):
        odf.loc[i, "gage_1"] = g_obs.loc[name, "obsval"]
        wdf.loc[i, "gage_1"] = g_obs.loc[name, "weight"]

    # now drop any entries that have zero weight across all cycles
    print(pst.nnz_obs_names)
    odf = odf.loc[:, pst.nnz_obs_names]
    wdf = wdf.loc[:, pst.nnz_obs_names]

    odf.T.to_csv(os.path.join(t_d, "obs_cycle_tbl.csv"))
    pst.pestpp_options["da_observation_cycle_table"] = "obs_cycle_tbl.csv"
    wdf.T.to_csv(os.path.join(t_d, "weight_cycle_tbl.csv"))
    pst.pestpp_options["da_weight_cycle_table"] = "weight_cycle_tbl.csv"

    pst.observation_data.loc[:, "cycle"] = -1
    pst.model_output_data.loc[:, "cycle"] = -1

    pst.write(os.path.join(t_d, "freyberg6_run_da2.pst"), version=2)
    return pst


def process_complex_target_output(c_d, s_d, d_d, real):
    # process for BAT
    redis_fac = 3

    start_date = pd.to_datetime('20151231', format='%Y%m%d')

    # load in obs ensemble
    oe_f = pd.read_csv(os.path.join(c_d, "freyberg.0.obs.csv"), index_col=0)
    oe_f = oe_f.T

    # process obs
    hds_f = oe_f.loc[oe_f.index.to_series().apply(lambda x: x.startswith("hds")), :].copy()
    hds_f.loc[:, "k"] = hds_f.index.to_series().apply(lambda x: int(x.split('_')[2]))
    hds_f.loc[:, "i"] = hds_f.index.to_series().apply(lambda x: int(x.split('_')[3]))
    hds_f.loc[:, "j"] = hds_f.index.to_series().apply(lambda x: int(x.split('_')[4]))
    hds_f.loc[:, "time"] = hds_f.index.to_series().apply(lambda x: float(x.split('_')[-1].split(':')[1]))
    time = hds_f.loc[:, "time"]
    time = pd.to_timedelta(time.values - 1, unit='D')
    hds_f.loc[:, "org_time"] = start_date + time.values
    hds_f.loc[:, "org_time"] = hds_f.loc[:, "org_time"].apply(lambda x: x.strftime('%Y%m%d'))
    hds_f.loc[:, "org_i"] = (hds_f.i / redis_fac).apply(np.int)
    hds_f.loc[:, "org_j"] = (hds_f.j / redis_fac).apply(np.int)
    hds_f.loc[:, "org_obgnme"] = hds_f.apply(lambda x: "trgw_{0}_{1}_{2}_{3}".format(x.k, x.org_i, x.org_j, x.org_time),
                                             axis=1)

    sfr_f = oe_f.loc[oe_f.index.to_series().apply(lambda x: x.startswith("sfr")), :].copy()
    sfr_f.loc[:, "time"] = sfr_f.index.to_series().apply(lambda x: float(x.split('_')[-1].split(':')[1]))
    type = sfr_f.index.to_series().apply(lambda x: x.split(':')[1].split('_')[0])
    type = pd.DataFrame(type)
    type = type.replace('gage', 'gage_1')
    sfr_f.loc[:, "type"] = type.values
    time = sfr_f.loc[:, "time"]
    time = pd.to_timedelta(time.values - 1, unit='D')
    sfr_f.loc[:, "org_time"] = start_date + time.values
    sfr_f.loc[:, "org_time"] = sfr_f.loc[:, "org_time"].apply(lambda x: x.strftime('%Y%m%d'))
    sfr_f.loc[:, "org_obgnme"] = sfr_f.apply(lambda x: "{0}_{1}".format(x.type, x.org_time), axis=1)

    pst = pyemu.Pst(os.path.join(s_d, 'freyberg.pst'))
    obs_s = pst.observation_data

    for j in range(len(hds_f)):
        cv = str(hds_f.org_obgnme[j])
        for i, sv in enumerate(obs_s.obsnme):
            if sv in cv:
                obs_s.obsval[i] = hds_f.iloc[j, real]
                break

    for j in range(len(sfr_f)):
        cv = str(sfr_f.org_obgnme[j])
        for i, sv in enumerate(obs_s.obsnme):
            if sv in cv:
                obs_s.obsval[i] = sfr_f.iloc[j, real]
                break

    # write pst files to template ies dir (for current real)
    m_ies_dir = os.path.join('simple_template_ies_{0}'.format(real))
    if os.path.exists(m_ies_dir):
        shutil.rmtree(m_ies_dir)
    shutil.copytree(s_d, m_ies_dir)
    pst.write(os.path.join(m_ies_dir, "freyberg.pst"), version=2)

    # process for SEQ
    redis_fac = 3
    pst = pyemu.Pst(os.path.join(d_d, "freyberg.pst"))

    start_date = pd.to_datetime('20151231', format='%Y%m%d')
    dates = pd.date_range(start='2015-12-31', periods=25, freq='M').strftime("%Y%m%d")

    # load in obs cycle table
    oe_f = pd.read_csv(os.path.join(c_d, "freyberg.0.obs.csv"), index_col=0)
    oe_f = oe_f.T
    obs_d = pd.read_csv(os.path.join(d_d, 'obs_cycle_tbl.csv'))
    obs_d.loc["date", 1:] = dates

    hds_f = oe_f.loc[oe_f.index.to_series().apply(lambda x: x.startswith("hds")), :].copy()
    hds_f.loc[:, "k"] = hds_f.index.to_series().apply(lambda x: int(x.split('_')[2]))
    hds_f.loc[:, "i"] = hds_f.index.to_series().apply(lambda x: int(x.split('_')[3]))
    hds_f.loc[:, "j"] = hds_f.index.to_series().apply(lambda x: int(x.split('_')[4]))
    hds_f.loc[:, "time"] = hds_f.index.to_series().apply(lambda x: float(x.split('_')[-1].split(':')[1]))
    time = hds_f.loc[:, "time"]
    time = pd.to_timedelta(time.values - 1, unit='D')
    hds_f.loc[:, "org_time"] = start_date + time.values
    hds_f.loc[:, "org_time"] = hds_f.loc[:, "org_time"].apply(lambda x: x.strftime('%Y%m%d'))
    hds_f.loc[:, "org_i"] = (hds_f.i / redis_fac).apply(np.int)
    hds_f.loc[:, "org_j"] = (hds_f.j / redis_fac).apply(np.int)
    hds_f.loc[:, "org_obgnme"] = hds_f.apply(lambda x: "head_{:02d}_{:03d}_{:03d}".format(x.k, x.org_i, x.org_j),
                                             axis=1)

    sfr_f = oe_f.loc[oe_f.index.to_series().apply(lambda x: x.startswith("sfr")), :].copy()
    sfr_f.loc[:, "time"] = sfr_f.index.to_series().apply(lambda x: float(x.split('_')[-1].split(':')[1]))
    type = sfr_f.index.to_series().apply(lambda x: x.split(':')[1].split('_')[0])
    type = pd.DataFrame(type)
    type = type.replace('gage', 'gage_1')
    sfr_f.loc[:, "type"] = type.values
    time = sfr_f.loc[:, "time"]
    time = pd.to_timedelta(time.values - 1, unit='D')
    sfr_f.loc[:, "org_time"] = start_date + time.values
    sfr_f.loc[:, "org_time"] = sfr_f.loc[:, "org_time"].apply(lambda x: x.strftime('%Y%m%d'))
    sfr_f.loc[:, "org_obgnme"] = sfr_f.apply(lambda x: "{0}".format(x.type), axis=1)

    for i in range(25):
        for j, (cv, ct) in enumerate(zip(hds_f.org_obgnme, hds_f.org_time)):
            for k in range(2):
                if obs_d.iloc[k, 0] in cv and obs_d.iloc[3, i] == ct:
                    obs_d.iloc[k, i] = hds_f.iloc[j, real]

    for i in range(25):
        for j, (cv, ct) in enumerate(zip(hds_f.org_obgnme, hds_f.org_time)):
            if obs_d.iloc[2, 0] in cv and obs_d.iloc[3, i] == ct:
                obs_d.iloc[2, i] = sfr_f.iloc[j, real]
    obs_d = obs_d.drop(['date'], axis=0)

    # write pst files and obs cycle table to master da dir (for current real)
    m_da_dir = os.path.join('simple_template_da_{0}'.format(real))
    if os.path.exists(m_da_dir):
        shutil.rmtree(m_ies_dir)
    shutil.copytree(d_d, m_da_dir)
    pst.write(os.path.join(m_da_dir, "freyberg.pst"), version=2)
    obs_d.to_csv(os.path.join(m_da_dir, 'obs_cycle_tbl.csv'), index=False)


def balance_weights(ireal):
    ies_dir = os.path.join('simple_template_ies_{0}'.format(ireal))
    da_dir = os.path.join('simple_template_da_{0}'.format(ireal))
    ies_file = 'freyberg6_run_ies.pst'

    # run pestpp ies to get phi
    pyemu.os_utils.run("pestpp-ies.exe {0}".format(ies_file), cwd=ies_dir)

    pst = pyemu.Pst(os.path.join(ies_dir, ies_file))
    obs = pst.observation_data

    # balance weights based on phi components
    pst.res.loc[:, 'residual'] = pst.res.loc[:, 'residual'] * pst.res.loc[:, 'weight']
    phi_comps = pst.res.groupby(["group"]).sum()
    phi_comp = phi_comps.loc[:, 'residual']
    phi_grps = phi_comps.index
    obs_phi_dict = dict(zip(phi_grps, phi_comp))
    pst._adjust_weights_by_phi_components(obs_phi_dict, original_ceiling=False)
    pst.write(os.path.join(ies_dir, ies_file), version=2)

    # modify da weight cycle table based on new weights
    start_date = pd.to_datetime('20151231', format='%Y%m%d')
    dates = pd.date_range(start='2015-12-31', periods=25, freq='M').strftime("%Y%m%d")

    da_wt = pd.read_csv(os.path.join(da_dir, 'weight_cycle_tbl.csv'))

    hds_f = obs.loc[obs.index.to_series().apply(lambda x: x.startswith("trgw")), :].copy()
    hds_f.loc[:, "time"] = hds_f.obsnme.apply(lambda x: x.split('_')[-1])
    hds_f.loc[:, "k"] = hds_f.obsnme.apply(lambda x: int(x.split('_')[1]))
    hds_f.loc[:, "i"] = hds_f.obsnme.apply(lambda x: int(x.split('_')[2]))
    hds_f.loc[:, "j"] = hds_f.obsnme.apply(lambda x: int(x.split('_')[3]))
    hds_f.loc[:, "org_obgnme"] = hds_f.apply(lambda x: "head_{:02d}_{:03d}_{:03d}".format(x.k, x.i, x.j),
                                             axis=1)
    sfr_f = obs.loc[obs.index.to_series().apply(lambda x: x.startswith("gage_1")), :].copy()
    sfr_f.loc[:, "time"] = sfr_f.obsnme.apply(lambda x: x.split('_')[-1])
    dates = pd.DataFrame(dates)
    dates.loc[:, 'cycle'] = range(25)
    cycle_dict = dict(zip(dates.iloc[:, 0], dates.loc[:, 'cycle']))
    hds_f.loc[:, 'time'] = hds_f.loc[:, 'time'].replace(cycle_dict, regex=True)
    sfr_f.loc[:, 'time'] = sfr_f.loc[:, 'time'].replace(cycle_dict, regex=True)

    for i in range(3):
        for j in range(24):
            for k in range(650):
                m = j + 2
                if da_wt.iloc[i, 0] in hds_f.iloc[k, 9] and hds_f.iloc[k, 5] == m:
                    da_wt.iloc[i, m] = hds_f.iloc[k, 2]

    for j in range(24):
        for k in range(25):
            m = j + 2
            if sfr_f.iloc[k, 5] == m:
                da_wt.iloc[2, m] = sfr_f.iloc[k, 2]

    da_wt.to_csv(os.path.join(da_dir, 'weight_cycle_tbl.csv'), index=False)


def compare_mf6_freyberg():
    for ireal in range(100):
        complex_dir = os.path.join('complex_master')
        # ies_dir = os.path.join('simple_template_ies')
        # da_dir = os.path.join('simple_template_da')
        ies_dir = os.path.join('monthly_template')
        da_dir = os.path.join('seq_monthly_template')

        process_complex_target_output(complex_dir, ies_dir, da_dir, ireal)

        balance_weights(ireal)
        return
        # run batch and sequential simple models
        # ies stuff
        ies_t_d = os.path.join('simple_template_ies_{0}'.format(ireal))
        ies_pst = pyemu.Pst(os.path.join(ies_t_d, "freyberg6_run_ies.pst"))

        # prep that prior ensemble for da
        da_t_d = os.path.join('simple_template_da_{0}'.format(ireal))
        da_pst = pyemu.Pst(os.path.join(da_t_d, "freyberg6_run_da.pst"))
        par = da_pst.parameter_data
        par.loc[par.parnme.str.contains("welflx"), "scale"] = -1.0

        ies_pst.pestpp_options["ies_par_en"] = "prior.jcb"
        ies_pe = pyemu.ParameterEnsemble.from_binary(pst=ies_pst, filename=os.path.join(ies_t_d, "prior.jcb"))
        da_pe = pyemu.ParameterEnsemble.from_gaussian_draw(pst=da_pst,
                                                           cov=pyemu.Cov.from_parameter_data(da_pst),
                                                           num_reals=ies_pe.shape[0])
        da_pe.index = ies_pe.index
        d = set(da_pe.columns.tolist()).symmetric_difference(set(ies_pe.columns.tolist()))
        print(d)
        da_pe.loc[:, ies_pe.columns] = ies_pe.values
        da_pe.to_binary(os.path.join(da_t_d, "da_prior.jcb"))
        da_pst.pestpp_options["ies_par_en"] = "da_prior.jcb"

        # set pestpp options for batch da
        ies_pst.pestpp_options.pop("ies_num_reals", None)
        ies_pst.pestpp_options.pop("da_num_reals", None)
        ies_pst.pestpp_options["ies_no_noise"] = False
        ies_pst.pestpp_options["ies_verbose_level"] = 1
        ies_pst.pestpp_options.pop("ies_localizer", None)
        ies_pst.pestpp_options["ies_autoadaloc"] = False
        ies_pst.pestpp_options["ies_save_lambda_en"] = False
        ies_pst.pestpp_options["ies_drop_conflicts"] = False
        ies_pst.pestpp_options["ies_num_reals"] = 50
        ies_pst.pestpp_options["ies_use_mda"] = False
        ies_pst.control_data.noptmax = 3
        ies_pst.write(os.path.join(ies_t_d, "freyberg6_run_ies.pst"), version=2)

        # set pestpp options for sequential da
        da_pst.pestpp_options.pop("da_num_reals", None)
        da_pst.pestpp_options.pop("ies_num_reals", None)
        da_pst.pestpp_options["ies_no_noise"] = False
        da_pst.pestpp_options["ies_verbose_level"] = 1
        da_pst.pestpp_options.pop("ies_localizer", None)
        da_pst.pestpp_options["ies_autoadaloc"] = False
        da_pst.pestpp_options["ies_save_lambda_en"] = False
        da_pst.pestpp_options["ies_drop_conflicts"] = False
        da_pst.pestpp_options["ies_num_reals"] = 50
        da_pst.pestpp_options["ies_use_mda"] = False
        da_pst.control_data.noptmax = 3
        da_pst.write(os.path.join(da_t_d, "freyberg6_run_da.pst"), version=2)

        # run da          
        m_da_dir = os.path.join('simple_master2_da_{0}'.format(ireal))

        pyemu.os_utils.start_workers(da_t_d, 'pestpp-da.exe', "freyberg6_run_da.pst", port=port,
                                     num_workers=4, master_dir=m_da_dir, verbose=True)

        shutil.rmtree(da_t_d)

        # run ies  
        m_ies_dir = os.path.join('simple_master2_ies_{0}'.format(ireal))

        pyemu.os_utils.start_workers(ies_t_d, 'pestpp-ies.exe', "freyberg6_run_ies.pst", port=port,
                                     num_workers=4, master_dir=m_ies_dir, verbose=True)

        shutil.rmtree(ies_t_d)


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


def monthly_ies_to_da_old(org_d="monthly_template"):
    org_sim = flopy.mf6.MFSimulation.load(sim_ws=org_d)
    org_totim = np.cumsum(org_sim.tdis.perioddata.array["perlen"])

    t_d = "seq_" + org_d
    if os.path.exists(t_d):
        shutil.rmtree(t_d)
    shutil.copytree(org_d, t_d)

    mm_df = pd.read_csv(os.path.join(t_d, "mult2model_info.csv"))
    print(mm_df.shape)
    drop_rch = mm_df.loc[mm_df.model_file.apply(lambda x: "rch_" in x and not "_1." in x), :].index
    mm_df = mm_df.drop(drop_rch)
    print(mm_df.shape)
    drop_wel = mm_df.loc[mm_df.model_file.apply(lambda x: ".wel_" in x and not "_1." in x), :].index
    mm_df = mm_df.drop(drop_wel)
    print(mm_df.shape)

    mm_df.to_csv(os.path.join(t_d, "mult2model_info.csv"))

    # first modify the tdis
    with open(os.path.join(t_d, "freyberg6.tdis"), 'w') as f:
        f.write("BEGIN Options\n  TIME_UNITS  days\nEND Options\n")
        f.write("BEGIN Dimensions\n  NPER  1\nEND Dimensions\n")
        f.write("BEGIN PERIODDATA\n31.00000000  1       1.00000000\nEND PERIODDATA\n")
    # make sure it runs
    pyemu.os_utils.run("mf6", cwd=t_d)

    # write a tdis template file - could possibly keep all 25 stress periods to
    # simulate a 2-year-ahead forecast...
    with open(os.path.join(t_d, "freyberg6.tdis.tpl"), 'w') as f:
        f.write("ptf  ~\n")
        f.write("BEGIN Options\n  TIME_UNITS  days\nEND Options\n")
        f.write("BEGIN Dimensions\n  NPER  1\nEND Dimensions\n")
        f.write("BEGIN PERIODDATA\n~  perlen  ~  1       1.00000000\nEND PERIODDATA\n")
    new_tpl, new_in = [os.path.join(t_d, "freyberg6.tdis.tpl")], [os.path.join(t_d, "freyberg6.tdis")]
    new_tpl_cycle = [-1]

    pyemu.os_utils.run("mf6", cwd=t_d)

    pst = pyemu.Pst(os.path.join(t_d, "freyberg.pst"))
    for tpl, inf, cy in zip(new_tpl, new_in, new_tpl_cycle):
        df = pst.add_parameters(tpl, inf, pst_path=".")
        pst.parameter_data.loc[df.parnme, "cycle"] = cy
        pst.parameter_data.loc[df.parnme, "partrans"] = "fixed"

    # write par cycle table
    pers = org_sim.tdis.perioddata.array["perlen"]
    # pers[0] = 1000
    pdf = pd.DataFrame(index=['perlen'], columns=np.arange(25))
    pdf.loc['perlen', :] = pers
    pdf.to_csv(os.path.join(t_d, "par_cycle_table.csv"))
    pst.pestpp_options["da_parameter_cycle_table"] = "par_cycle_table.csv"

    # add initial condition parameters (all cells)
    ic_files = [f for f in os.listdir(t_d) if "ic_strt" in f.lower() and f.lower().endswith(".txt")]
    print(ic_files)

    for ic_file in ic_files:
        k = int(ic_file.split(".")[1][-1]) - 1
        # ib = org_sim.get_model("freyberg6").dis.idomain[k].array
        ib = np.loadtxt(os.path.join(t_d, 'freyberg6.dis_idomain_layer{0}.txt'.format(k + 1)))
        arr = np.loadtxt(os.path.join(t_d, ic_file))
        tpl_file = os.path.join(t_d, ic_file + ".tpl")
        ic_val_dict = {}
        with open(tpl_file, 'w') as f:
            f.write("ptf ~\n")
            for i in range(arr.shape[0]):
                for j in range(arr.shape[1]):
                    if ib[i, j] > 0 and arr[i, j] > 0:
                        pname = "head_k:{0}_i:{1}_j:{2}".format(k, i, j)
                        f.write(" ~ {0} ~ ".format(pname))
                        ic_val_dict[pname] = arr[i, j]
                    else:
                        f.write(" -9999999 ".format(pname))
                f.write("\n")
        df = pst.add_parameters(tpl_file, pst_path=".")
        for pname, parval1 in ic_val_dict.items():
            pst.parameter_data.loc[pname, "parval1"] = parval1
            pst.parameter_data.loc[pname, "partrans"] = "none"
            pst.parameter_data.loc[pname, "parchglim"] = "relative"

            pst.parameter_data.loc[pname, "parlbnd"] = -1000000
            pst.parameter_data.loc[pname, "parubnd"] = 10000000

    # and add final simulated water level observations (all cells - will need a new post processor)
    # and link these two via the observation data state_par_link attr.
    frun = os.path.join(t_d, "forward_run.py")
    frun_lines = open(frun).readlines()
    frun_lines.append("import flopy\n")
    frun_lines.append("hds = flopy.utils.HeadFile('freyberg6_freyberg.hds')\n")
    frun_lines.append("arr = hds.get_data()\n")
    frun_lines.append("for k,a in enumerate(arr):\n")
    frun_lines.append("    np.savetxt('sim_head_{0}.dat'.format(k),a,fmt='%15.6E')\n")
    with open(frun, 'w') as f:
        for line in frun_lines:
            f.write(line)
    hds = flopy.utils.HeadFile(os.path.join(t_d, "freyberg6_freyberg.hds"))
    arr = hds.get_data()
    new_ins_files = []
    new_out = []
    new_ins_cycle = []

    for k in range(arr.shape[0]):
        out_file = os.path.join(t_d, "sim_head_{0}.dat".format(k))
        np.savetxt(out_file, arr[k, :, :], fmt="%15.6E")
        with open(out_file + ".ins", 'w') as f:
            f.write("pif ~\n")
            for i in range(arr.shape[1]):
                f.write("l1 ")
                for j in range(arr.shape[2]):
                    oname = "head_k:{0}_i:{1}_j:{2}".format(k, i, j)
                    f.write(" !{0}! ".format(oname))
                f.write("\n")
        # just to check that the ins file is valid...
        new_ins_files.append(out_file + ".ins")
        i = pyemu.pst_utils.InstructionFile(out_file + ".ins")
        df = i.read_output_file(out_file)

    # save this for later!
    org_obs = pst.observation_data.copy()

    # now drop all existing heads obs since those will be replaced by the state obs
    fname = 'heads.csv'
    fname_ins = fname + ".ins"
    pst.drop_observations(os.path.join(t_d, fname_ins), '.')

    sfr = pd.read_csv(os.path.join(t_d, 'sfr.csv'))
    sfr = sfr.drop(['time'], axis=1)
    nms = sfr.columns.to_series()
    fname = 'sfr.csv'
    fname_ins = fname + ".ins"
    pst.drop_observations(os.path.join(t_d, fname_ins), '.')
    with open(os.path.join(t_d, fname_ins), 'w') as f:
        f.write("pif ~\n")
        f.write("l1 \n")
        f.write("l1")
        for i in range(sfr.shape[1]):
            oname = "{0}".format(nms[i])
            f.write(" ~,~ !{0}! ".format(oname))
    new_ins_files.append(os.path.join(t_d, fname_ins))
    new_out.append(os.path.join(t_d, "sfr.csv"))
    new_ins_cycle.append(-1)

    # add sfr obs
    for ins_file in new_ins_files:
        pst.add_observations(ins_file, pst_path=".")

    pst.observation_data.loc[:, 'cycle'] = -1
    tr_obs = org_obs.loc[org_obs.obsnme.str.contains("hds_usecol:trgw"), :].copy()
    tr_obs.loc[tr_obs.obsnme, "time"] = tr_obs.obsnme.apply(lambda x: x.split(':')[-1])
    tr_obs.loc[tr_obs.obsnme, "k"] = tr_obs.obsnme.apply(lambda x: np.int(x.split('_')[2]))
    tr_obs.loc[tr_obs.obsnme, "i"] = tr_obs.obsnme.apply(lambda x: np.int(x.split('_')[3]))
    tr_obs.loc[tr_obs.obsnme, "j"] = tr_obs.obsnme.apply(lambda x: np.int(x.split('_')[4]))
    tr_obs.loc[tr_obs.obsnme, "obgnme"] = tr_obs.obsnme.apply(lambda x: "_".join(x.split("_")[:-1]))

    head_obs = pst.observation_data.loc[pst.observation_data.obsnme.str.startswith("head_"), :].copy()
    head_obs.loc[head_obs.obsnme, "k"] = head_obs.obsnme.apply(lambda x: np.int(x.split('_')[1].split(':')[1]))
    head_obs.loc[head_obs.obsnme, "i"] = head_obs.obsnme.apply(lambda x: np.int(x.split('_')[2].split(':')[1]))
    head_obs.loc[head_obs.obsnme, "j"] = head_obs.obsnme.apply(lambda x: np.int(x.split('_')[3].split(':')[1]))

    obs_heads = {}
    odf_names = []
    pst.observation_data.loc[:, "org_obgnme"] = np.NaN
    pst.observation_data.loc[:, "weight"] = 0.0
    pst.observation_data.loc["gage_1", "weight"] = 1.0

    for og in tr_obs.obgnme.unique():
        site_obs = tr_obs.loc[tr_obs.obgnme == og, :]
        site_obs.sort_values(by="time", inplace=True)
        head_name = "head_k:{0}_i:{1}_j:{2}".format(site_obs.k[0], site_obs.i[0], site_obs.j[0])
        for i, oname in enumerate(site_obs.obsnme):
            obs_heads[oname] = (head_name, i)
        # assign max weight in the control file since some have zero weight and
        # we are covering weights in the weight table
        pst.observation_data.loc[head_name, "weight"] = site_obs.weight.max()
        pst.observation_data.loc[head_name, "org_obgnme"] = og
        odf_names.append(head_name)
    odf_names.append("gage_1")

    odf = pd.DataFrame(columns=odf_names, index=np.arange(25))
    wdf = pd.DataFrame(columns=odf_names, index=np.arange(25))
    for tr_name, (head_name, icycle) in obs_heads.items():
        if icycle > 12:
            odf.loc[icycle, head_name] = org_obs.loc[tr_name, "obsval"]
            wdf.loc[icycle, head_name] = 0
        else:
            odf.loc[icycle, head_name] = org_obs.loc[tr_name, "obsval"]
            wdf.loc[icycle, head_name] = org_obs.loc[tr_name, "weight"]

    g_obs = org_obs.loc[org_obs.obsnme.str.startswith("sfr_usecol:gage_1"), :].copy()
    # give these obs the max weight since some have zero weight
    # pst.observation_data.loc["gage_1", "weight"] = g_obs.weight.max()
    g_obs.sort_index(inplace=True)
    for i, name in enumerate(g_obs.obsnme):
        if i > 12:
            odf.loc[i, "gage_1"] = g_obs.loc[name, "obsval"]
            wdf.loc[i, "gage_1"] = 0
        else:
            odf.loc[i, "gage_1"] = g_obs.loc[name, "obsval"]
            wdf.loc[i, "gage_1"] = g_obs.loc[name, "weight"]

    odf.T.to_csv(os.path.join(t_d, "obs_cycle_tbl.csv"))
    pst.pestpp_options["da_observation_cycle_table"] = "obs_cycle_tbl.csv"
    wdf.T.to_csv(os.path.join(t_d, "weight_cycle_tbl.csv"))
    pst.pestpp_options["da_weight_cycle_table"] = "weight_cycle_tbl.csv"

    # need to set cycle vals and reset the model_file attr for each cycle-specific template files (rch and wel)
    pst.model_input_data.loc[:, "cycle"] = -1
    for i in range(len(pst.model_input_data)):
        if pst.model_input_data.iloc[i, 0].startswith('wel_grid'):
            cy = int(pst.model_input_data.iloc[i, 0].split('_')[2])
            pst.model_input_data.iloc[i, 2] = cy - 1
            pst.model_input_data.iloc[i, 1] = pst.model_input_data.iloc[i, 1].replace(str(cy), "1")
        if pst.model_input_data.iloc[i, 0].startswith('twel_mlt'):
            cy = int(pst.model_input_data.iloc[i, 0].split('_')[2])
            pst.model_input_data.iloc[i, 2] = cy - 1
            pst.model_input_data.iloc[i, 1] = pst.model_input_data.iloc[i, 1].replace(str(cy), "1")
        elif 'rch_recharge' in pst.model_input_data.iloc[i, 0] and "cn" in pst.model_input_data.iloc[i, 0]:
            cy = int(pst.model_input_data.iloc[i, 0].split('_')[2])
            pst.model_input_data.iloc[i, 2] = cy - 1
            pst.model_input_data.iloc[i, 1] = pst.model_input_data.iloc[i, 1].replace(str(cy), "1")

    pst.model_output_data.loc[:, "cycle"] = -1

    # need to set the cycle value for all pars - static properties and multi-stress period broadcast forcing pars
    # should get a cycle value of -1.
    pst.parameter_data.loc[:, "cycle"] = -1

    for i in range(len(pst.parameter_data)):
        if pst.parameter_data.iloc[i, 0].startswith('wel_grid'):
            cy = int(pst.parameter_data.iloc[i, 0].split('_')[2])
            pst.parameter_data.iloc[i, 22] = cy - 1
        elif pst.parameter_data.iloc[i, 0].startswith('twel_mlt'):
            cy = int(pst.parameter_data.iloc[i, 0].split('_')[2])
            pst.parameter_data.iloc[i, 22] = cy - 1
        elif 'rch_recharge' in pst.parameter_data.iloc[i, 0] and "cn" in pst.parameter_data.iloc[i, 0]:
            cy = int(pst.parameter_data.iloc[i, 0].split('_')[4])
            pst.parameter_data.iloc[i, 22] = cy - 1

    # pst.observation_data.loc[:, "state_par_link"] = ''
    # print(pst.observation_data.iloc[2429,:])
    # for i in range(len(pst.observation_data)):
    #    if pst.observation_data.iloc[i, 0].startswith('head_'):
    #        pst.observation_data.iloc[i,9] = pst.observation_data.iloc[i,0]

    pst.control_data.noptmax = 3
    # # pst.pestpp_options["ies_num_reals"] = 3
    pst.pestpp_options["da_num_reals"] = 50
    # # if not sync_state_names:
    # #     pst.observation_data.loc[:,"state_par_link"] = np.NaN
    # #     obs = pst.observation_data
    # #     obs.loc[:,"state_par_link"] = obs.obsnme.apply(lambda x: obs_to_par_map.get(x,np.NaN))
    print(pst.nnz_obs_names)
    pst.write(os.path.join(t_d, "freyberg6_run_da.pst"), version=2)
    # return pst

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

    pst.pestpp_options["da_num_reals"] = 100
    pst.control_data.noptmax = 3
    pst.write(os.path.join(t_d, "test.pst"), version=2)
    # pyemu.os_utils.run("pestpp-da test.pst",cwd=t_d)
    return
    pst.pestpp_options["da_num_reals"] = 100
    pst.write(os.path.join(t_d, "test.pst"), version=2)
    pyemu.os_utils.start_workers(t_d, "pestpp-da", "test.pst", num_workers=10, master_dir=t_d + "prior_test")


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
                #pf.add_parameters(filenames=arr_file, par_type="constant", par_name_base=arr_file.split('.')[1] + "_cn",
                #                  pargp="rch_const", zone_array=ib, upper_bound=ub, lower_bound=lb,
                #                  geostruct=temporal_gs,
                #                  datetime=dts[kper])
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
    pst.pestpp_options["additional_ins_delimiters"] = ","

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

    # write the control file
    pst.write(os.path.join(pf.new_d, "freyberg.pst"))

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
    pst.write(os.path.join(pf.new_d, "freyberg.pst"))


def monthly_ies_to_da(org_d):
    org_sim = flopy.mf6.MFSimulation.load(sim_ws=org_d)
    org_totim = np.cumsum(org_sim.tdis.perioddata.array["perlen"])

    t_d = "seq_" + org_d
    if os.path.exists(t_d):
        shutil.rmtree(t_d)
    shutil.copytree(org_d, t_d)

    mm_df = pd.read_csv(os.path.join(t_d, "mult2model_info.csv"))
    print(mm_df.shape)
    drop_rch = mm_df.loc[mm_df.model_file.apply(lambda x: "rch_" in x and not "_1." in x), :].index
    mm_df = mm_df.drop(drop_rch)
    print(mm_df.shape)
    drop_wel = mm_df.loc[mm_df.model_file.apply(lambda x: ".wel_" in x and not "_1." in x), :].index
    mm_df = mm_df.drop(drop_wel)
    print(mm_df.shape)
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

    pyemu.os_utils.run("mf6", cwd=t_d)

    pst = pyemu.Pst(os.path.join(t_d, "freyberg.pst"))
    for tpl, inf, cy in zip(new_tpl, new_in, new_tpl_cycle):
        df = pst.add_parameters(tpl, inf, pst_path=".")
        pst.parameter_data.loc[df.parnme, "cycle"] = cy
        pst.parameter_data.loc[df.parnme, "partrans"] = "fixed"

    # write par cycle table
    pers = org_sim.tdis.perioddata.array["perlen"]
    # pers[0] = 1000
    pdf = pd.DataFrame(index=['perlen'], columns=np.arange(25))
    pdf.loc['perlen', :] = pers
    pdf.to_csv(os.path.join(t_d, "par_cycle_table.csv"))
    pst.pestpp_options["da_parameter_cycle_table"] = "par_cycle_table.csv"

    # save this for later!
    org_obs = pst.observation_data.copy()

    # now drop all existing heads obs since those will be replaced by the state obs
    fname = 'heads.csv'
    fname_ins = fname + ".ins"
    pst.drop_observations(os.path.join(t_d, fname_ins), '.')

    # fix the sfr obs
    sfr = pd.read_csv(os.path.join(t_d, 'sfr.csv'))
    sfr = sfr.drop(['time'], axis=1)
    nms = sfr.columns.to_series()
    fname = 'sfr.csv'
    fname_ins = fname + ".ins"
    pst.drop_observations(os.path.join(t_d, fname_ins), '.')
    with open(os.path.join(t_d, fname_ins), 'w') as f:
        f.write("pif ~\n")
        f.write("l1 \n")
        f.write("l1")
        for i in range(sfr.shape[1]):
            oname = "{0}".format(nms[i])
            f.write(" ~,~ !{0}! ".format(oname))
    new_ins_files = [os.path.join(t_d, fname_ins)]
    new_out = [os.path.join(t_d, "sfr.csv")]
    new_ins_cycle = [-1]

    # add sfr obs
    for ins_file in new_ins_files:
        pst.add_observations(ins_file, pst_path=".")

    pst.observation_data.loc[:, 'cycle'] = -1
    tr_obs = org_obs.loc[org_obs.obsnme.str.contains("hds_usecol:trgw"), :].copy()
    tr_obs.loc[tr_obs.obsnme, "time"] = tr_obs.obsnme.apply(lambda x: x.split(':')[-1])
    tr_obs.loc[tr_obs.obsnme, "k"] = tr_obs.obsnme.apply(lambda x: np.int(x.split('_')[2]))
    tr_obs.loc[tr_obs.obsnme, "i"] = tr_obs.obsnme.apply(lambda x: np.int(x.split('_')[3]))
    tr_obs.loc[tr_obs.obsnme, "j"] = tr_obs.obsnme.apply(lambda x: np.int(x.split('_')[4]))
    tr_obs.loc[tr_obs.obsnme, "obgnme"] = tr_obs.obsnme.apply(lambda x: "_".join(x.split("_")[:-1]))

    head_obs = pst.observation_data.loc[pst.observation_data.obsnme.str.startswith("arrobs_head_"), :].copy()
    # head_obs.loc[head_obs.obsnme, "k"] = head_obs.obsnme.apply(lambda x: np.int(x.split('_')[1].split(':')[1]))
    # head_obs.loc[head_obs.obsnme, "i"] = head_obs.obsnme.apply(lambda x: np.int(x.split('_')[2].split(':')[1]))
    # head_obs.loc[head_obs.obsnme, "j"] = head_obs.obsnme.apply(lambda x: np.int(x.split('_')[3].split(':')[1]))
    for v in ["k", "i", "j"]:
        head_obs.loc[:, v] = head_obs.loc[:, v].apply(int)

    obs_heads = {}
    odf_names = []
    pst.observation_data.loc[:, "org_obgnme"] = np.NaN
    pst.observation_data.loc[:, "weight"] = 0.0
    pst.observation_data.loc["gage_1", "weight"] = 1.0

    for og in tr_obs.obgnme.unique():
        site_obs = tr_obs.loc[tr_obs.obgnme == og, :]
        site_obs.sort_values(by="time", inplace=True)
        head_name = "arrobs_head_k:{0}_i:{1}_j:{2}".format(site_obs.k[0], site_obs.i[0], site_obs.j[0])
        for i, oname in enumerate(site_obs.obsnme):
            obs_heads[oname] = (head_name, i)
        # assign max weight in the control file since some have zero weight and
        # we are covering weights in the weight table
        pst.observation_data.loc[head_name, "weight"] = site_obs.weight.max()
        pst.observation_data.loc[head_name, "org_obgnme"] = og
        odf_names.append(head_name)
    odf_names.append("gage_1")

    odf = pd.DataFrame(columns=odf_names, index=np.arange(25))
    wdf = pd.DataFrame(columns=odf_names, index=np.arange(25))
    for tr_name, (head_name, icycle) in obs_heads.items():
        odf.loc[icycle, head_name] = org_obs.loc[tr_name, "obsval"]
        if icycle > 12:
            wdf.loc[icycle, head_name] = 0
        else:
            wdf.loc[icycle, head_name] = org_obs.loc[tr_name, "weight"]

    g_obs = org_obs.loc[org_obs.obsnme.str.startswith("sfr_usecol:gage_1"), :].copy()
    # give these obs the max weight since some have zero weight
    # pst.observation_data.loc["gage_1", "weight"] = g_obs.weight.max()
    g_obs.sort_index(inplace=True)
    for i, name in enumerate(g_obs.obsnme):
        odf.loc[i, "gage_1"] = g_obs.loc[name, "obsval"]
        if i > 12:
            wdf.loc[i, "gage_1"] = 0
        else:
            wdf.loc[i, "gage_1"] = g_obs.loc[name, "weight"]

    odf.T.to_csv(os.path.join(t_d, "obs_cycle_tbl.csv"))
    pst.pestpp_options["da_observation_cycle_table"] = "obs_cycle_tbl.csv"
    wdf.T.to_csv(os.path.join(t_d, "weight_cycle_tbl.csv"))
    pst.pestpp_options["da_weight_cycle_table"] = "weight_cycle_tbl.csv"

    # need to set cycle vals and reset the model_file attr for each cycle-specific template files (rch and wel)
    pst.model_input_data.loc[:, "cycle"] = -1
    for i in range(len(pst.model_input_data)):
        # if pst.model_input_data.iloc[i,0].startswith('wel_grid'):
        #     cy = int(pst.model_input_data.iloc[i,0].split('_')[2])
        #     pst.model_input_data.iloc[i, 2] = cy
        #     pst.model_input_data.iloc[i, 1] = pst.model_input_data.iloc[i, 1].replace(str(cy), "1")
        if pst.model_input_data.iloc[i, 0].startswith('twel_mlt'):
            cy = int(pst.model_input_data.iloc[i, 0].split('_')[2])
            pst.model_input_data.iloc[i, 2] = cy
            pst.model_input_data.iloc[i, 1] = pst.model_input_data.iloc[i, 1].replace(str(cy), "1")
        elif 'rch_recharge' in pst.model_input_data.iloc[i, 0] and "cn" in pst.model_input_data.iloc[i, 0]:
            cy = int(pst.model_input_data.iloc[i, 0].split('_')[2])
            pst.model_input_data.iloc[i, 2] = cy - 1
            pst.model_input_data.iloc[i, 1] = pst.model_input_data.iloc[i, 1].replace(str(cy), "1")
        elif 'rch_recharge' in pst.model_input_data.iloc[i, 0] and pst.model_input_data.iloc[i, 0].endswith(".txt.tpl"):
            cy = int(pst.model_input_data.iloc[i, 0].split('.')[1].split("_")[-1])
            pst.model_input_data.iloc[i, 2] = cy - 1
            pst.model_input_data.iloc[i, 1] = pst.model_input_data.iloc[i, 1].replace(str(cy), "1")

    pst.model_output_data.loc[:, "cycle"] = -1

    # need to set the cycle value for all pars - static properties and multi-stress period broadcast forcing pars
    # should get a cycle value of -1.
    pst.parameter_data.loc[:, "cycle"] = -1

    for pname in pst.par_names:
        # if pst.parameter_data.iloc[i,0].startswith('wel_grid'):
        #     cy = int(pst.parameter_data.iloc[i,0].split('_')[2])
        #     pst.parameter_data.iloc[i, 22] = cy
        if pname.startswith('twel_mlt'):
            cy = int(pname.split('_')[2])
            pst.parameter_data.loc[pname, "cycle"] = cy
        elif 'rch_recharge' in pname and "cn" in pname:
            cy = int(pname.split('_')[4])
            pst.parameter_data.loc[pname, "cycle"] = cy - 1



    pst.control_data.noptmax = 3
    pst.write(os.path.join(t_d, "freyberg.pst"), version=2)
    # return pst

    # fill any missing pars with ctl file values (esp the perlen par)
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

    pst.control_data.noptmax = 0
    pst.write(os.path.join(t_d, "test.pst"), version=2)
    pyemu.os_utils.run("pestpp-da test.pst", cwd=t_d)


def run_batch_seq_prior_monte_carlo():
    t_d = "seq_monthly_model_files_template"
    pst = pyemu.Pst(os.path.join(t_d,"freyberg.pst"))
    pst.control_data.noptmax = -1
    pst.write(os.path.join(t_d,"freyberg.pst"))
    pyemu.os_utils.start_workers(t_d, "pestpp-da", "freyberg.pst", num_workers=10,
                                 master_dir=t_d.replace("template", "master_prior"))

    t_d = "monthly_model_files_template"
    pst = pyemu.Pst(os.path.join(t_d, "freyberg.pst"))
    pst.control_data.noptmax = -1
    pst.write(os.path.join(t_d, "freyberg.pst"))
    pyemu.os_utils.start_workers(t_d, "pestpp-ies", "freyberg.pst", num_workers=10,
                                 master_dir=t_d.replace("template", "master_prior"))


def plot_prior_mc():
    c_m_d = "daily_model_files_master_prior"
    s_b_m_d = "monthly_model_files_master_prior"
    s_s_m_d = "seq_monthly_model_files_master_prior"

    c_pst = pyemu.Pst(os.path.join(c_m_d, "freyberg.pst"))
    obs = c_pst.observation_data
    ctr_obs = obs.loc[obs.obsnme.str.contains("trgw"), :].copy()
    ctr_obs.loc[:, "k"] = ctr_obs.obsnme.apply(lambda x: int(x.split("_")[2]))
    ctr_obs.loc[:, "i"] = ctr_obs.obsnme.apply(lambda x: int(x.split("_")[3]))
    ctr_obs.loc[:, "j"] = ctr_obs.obsnme.apply(lambda x: int(x.split("_")[4]))
    ctr_obs.loc[:, "time"] = ctr_obs.time.apply(float)
    print(ctr_obs.obgnme.unique())


    s_b_pst = pyemu.Pst(os.path.join(s_b_m_d, "freyberg.pst"))
    obs = s_b_pst.observation_data
    tr_obs = obs.loc[obs.obsnme.str.contains("trgw"),:].copy()
    tr_obs.loc[:,"k"] = tr_obs.obsnme.apply(lambda x: int(x.split("_")[2]))
    tr_obs.loc[:, "i"] = tr_obs.obsnme.apply(lambda x: int(x.split("_")[3]))
    tr_obs.loc[:, "j"] = tr_obs.obsnme.apply(lambda x: int(x.split("_")[4]))
    tr_obs.loc[:,"time"] = tr_obs.time.apply(float)

    c_oe = pd.read_csv(os.path.join(c_m_d, "freyberg.0.obs.csv"), index_col=0)
    s_b_oe = pd.read_csv(os.path.join(s_b_m_d, "freyberg.0.obs.csv"), index_col=0)

    s_s_pst = pyemu.Pst(os.path.join(s_s_m_d,"freyberg.pst"))
    seq_oe_files = [f for f in os.listdir(s_s_m_d) if f.endswith(".oe.csv") and "global" in f and f.startswith("freyberg")]
    s_s_oe_dict = {int(f.split(".")[2]):pd.read_csv(os.path.join(s_s_m_d,f),index_col=0) for f in seq_oe_files}
    #oct = pd.read_csv(os.path.join(s_s_m_d,"obs_cycle_tbl.csv"),index_col=0)

    ognames = list(tr_obs.obgnme.unique())
    ognames.sort()
    ognames.append("sfr_usecol:gage_1")
    with PdfPages("prior_obs_v_sim.pdf") as pdf:

        for ogname in ognames:
            sbobs = tr_obs.loc[tr_obs.obgnme == ogname, :].copy()

            if "gage" in ogname:
                seq_name ="gage_1"
                cogname = ogname
                sbobs = s_b_pst.observation_data.loc[s_b_pst.observation_data.obgnme==ogname,:].copy()
                sbobs.loc[:,"time"] = sbobs.time.apply(float)
                cobs = c_pst.observation_data.loc[c_pst.observation_data.obgnme == ogname, :].copy()
                cobs.loc[:, "time"] = cobs.time.apply(float)

            else:
                sbobs = tr_obs.loc[tr_obs.obgnme == ogname, :].copy()
                sbobs.sort_values(by="time", inplace=True)
                k,i,j = sbobs.k.iloc[0],sbobs.i.iloc[0],sbobs.j.iloc[0]
                seq_name = "arrobs_head_k:{0}_i:{1}_j:{2}".format(k,i,j)
                cogname = "hds_usecol:trgw_{0}_{1}_{2}".format(k,2+(i*3),2+(j*3))
                print(ogname,cogname)
                cobs = ctr_obs.loc[ctr_obs.obgnme==cogname,:].copy()
            sbobs.sort_values(by="time", inplace=True)
            cobs.sort_values(by="time",inplace=True)


            fig,ax = plt.subplots(1,1,figsize=(8,8))

            [ax.plot(cobs.time, c_oe.loc[idx, cobs.obsnme], "b", lw=0.01,alpha=0.5) for idx in c_oe.index]

            [ax.plot(sbobs.time,s_b_oe.loc[idx,sbobs.obsnme],"0.5",lw=0.01,alpha=0.5) for idx in s_b_oe.index]
            for itime,time in enumerate(sbobs.time):
                oe = s_s_oe_dict[itime]
                #print(oe.loc[:,seq_name])
                ax.scatter([time for _ in range(oe.shape[0])], oe.loc[:, seq_name], marker=".",color="0.5",alpha=0.5)


            ax.set_title(ogname)
            if "gage" not in ogname:
                ax.set_ylim(30,ax.get_ylim()[1])
            pdf.savefig()
            plt.close(fig)









if __name__ == "__main__":

    setup_interface("monthly_model_files")
    monthly_ies_to_da("monthly_model_files_template")
    run_batch_seq_prior_monte_carlo()
    setup_interface("daily_model_files")
    run_complex_prior_mc('daily_model_files_template')
    plot_prior_mc()
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
