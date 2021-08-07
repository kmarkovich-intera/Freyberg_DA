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

port = 4021

#set path to pestpp executables
bin_path = os.path.join("bin")
if "linux" in platform.platform().lower():
    bin_path = os.path.join(bin_path,"linux")
elif "macos" in platform.platform().lower():
    bin_path = os.path.join(bin_path,"mac")
else:
    bin_path = os.path.join(bin_path,"win")
    
exe = ""
if "windows" in platform.platform().lower():
    exe = ".exe"
da_path = os.path.join(bin_path, "pestpp-da" + exe)
ies_path = os.path.join(bin_path, "pestpp-ies" + exe)


def mod_tdis_sto(org_t_d,t_d):
    tdis_file = "freyberg6.tdis"
    lines = open(os.path.join(org_t_d,tdis_file),'r').readlines()
 
    with open(os.path.join(t_d,tdis_file),'w') as f:
        iline = 0
        while True:
            line = lines[iline]
            if "begin period" in line.lower():
                lines[iline+1] = "100000  1   1.0\n"
                print(lines[iline+1])
            print(line)
            f.write(line)
            iline += 1
            if iline >= len(lines):
                break
    sto_file = "freyberg6.sto"
    lines = open(os.path.join(org_t_d,sto_file),'r').readlines()
 
    with open(os.path.join(t_d,sto_file),'w') as f:
 
         for line in lines:
             f.write(line)
             if line.lower().startswith("end griddata"):
                 break
         f.write("\nbegin period 1\n  transient\nend period 1\n")

def da_prep_4_mf6_freyberg_seq(sync_state_names=True):
    t_d = os.path.join("mf6_freyberg","template_seq")
    if os.path.exists(t_d):
        shutil.rmtree(t_d)
    shutil.copytree(os.path.join("mf6_freyberg","template"),t_d)
    for f in os.listdir(t_d):
        for tag in ["ies","opt","glm"]:
            if tag in f.lower():
                os.remove(os.path.join(t_d,f))

    #first modify the tdis
    with open(os.path.join(t_d,"freyberg6.tdis"),'w') as f:
        f.write("BEGIN Options\n  TIME_UNITS  days\nEND Options\n")
        f.write("BEGIN Dimensions\n  NPER  1\nEND Dimensions\n")
        f.write("BEGIN PERIODDATA\n31.00000000  1       1.00000000\nEND PERIODDATA\n")
    #make sure it runs
    pyemu.os_utils.run("mf6",cwd=t_d)

    # write a tdis template file - could possibly keep all 25 stress periods to
    # simulate a 2-year-ahead forecast...

    with open(os.path.join(t_d,"freyberg6.tdis.tpl"),'w') as f:
        f.write("ptf  ~\n")
        f.write("BEGIN Options\n  TIME_UNITS  days\nEND Options\n")
        f.write("BEGIN Dimensions\n  NPER  1\nEND Dimensions\n")
        f.write("BEGIN PERIODDATA\n~  perlen  ~  1       1.00000000\nEND PERIODDATA\n")
    new_tpl,new_in = [os.path.join(t_d,"freyberg6.tdis.tpl")],[os.path.join(t_d,"freyberg6.tdis")]
    new_tpl_cycle = [-1]

    # split out the head, sfr and list instruction files into multiple instruction file
    #eventually, want to move to da par and obs cycle tables for heads and gage_1 obs
    lines = open(os.path.join(t_d,"heads.csv.ins"),'r').readlines()[2:]
    new_ins,new_out,new_ins_cycle = [],[],[]
    #print(lines)
    for icycle, line in enumerate(lines):
        ins_name = os.path.join(t_d,"heads_{0}.csv.ins".format(icycle))
        with open(ins_name,'w') as f:
            f.write("pif ~\nl1\n")
            f.write(line)
        new_ins.append(ins_name)
        new_out.append(os.path.join(t_d,"heads.csv"))
        new_ins_cycle.append(icycle)
    remove_ins = ["heads.csv.ins"]

    lines = open(os.path.join(t_d,"sfr.csv.ins"),'r').readlines()[2:]
    #print(lines)
    for icycle, line in enumerate(lines):
        ins_name = os.path.join(t_d,"sfr_{0}.csv.ins".format(icycle))
        with open(ins_name,'w') as f:
            f.write("pif ~\nl1\n")
            f.write(line)
        new_ins.append(ins_name)
        new_out.append(os.path.join(t_d,"sfr.csv"))
        new_ins_cycle.append(icycle)
    remove_ins.append("sfr.csv.ins")

    lines = open(os.path.join(t_d,"freyberg6.lst.ins"),'r').readlines()[1:]
    icycle = 0
    tag_line = lines[0]
    for s in range(0,len(lines),13):
        ins_name = os.path.join(t_d,"freyberg6_{0}.lst.ins".format(icycle))
        with open(os.path.join(t_d,"freyberg6_{0}.lst.ins".format(icycle)),'w') as f:
            f.write("pif ~\n")
            f.write(tag_line)
            for line in lines[s+1:s+13]:
                f.write(line)
        new_ins.append(ins_name)
        new_out.append(os.path.join(t_d,"freyberg6.lst"))
        new_ins_cycle.append(icycle)
        icycle += 1
    remove_ins.append("freyberg6.lst.ins")

    # modify the ic file
    k = 0
    with open(os.path.join(t_d,"freyberg6.ic"),'r') as f:
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
                        arr = np.array(arr_lines,dtype=np.float)
                        np.savetxt(os.path.join(t_d,"heads_{0}.dat_in".format(k)),arr,fmt="%15.6E")
                        k += 1
                        arr_lines = []
                    else:
                        arr_lines.append(line.strip().split())

        arr = np.array(arr_lines, dtype=np.float)
        np.savetxt(os.path.join(t_d, "heads_{0}.dat_in".format(k)), arr, fmt="%15.6E")
    with open(os.path.join(t_d,"freyberg6.ic"),'w') as f:
        f.write("begin griddata\nstrt layered\n")
        for k in range(3):
            f.write("open/close 'heads_{0}.dat_in' FACTOR 1.0\n".format(k))
        f.write("end griddata\n")


    # write a python script to extract final heads and save to files
    with open(os.path.join(t_d,"forward_run.py"),'w') as f:
        f.write("import numpy as np\nimport flopy\nimport pyemu\n")
        f.write("pyemu.os_utils.run('mf6')\n")
        f.write("hds = flopy.utils.HeadFile('freyberg6_freyberg.hds')\n")
        f.write("arr = hds.get_data()\n")
        f.write("for k,a in enumerate(arr):\n")
        f.write("    np.savetxt('heads_'+str(k)+'.dat',a,fmt='%15.6E')\n")

    # dont run it so we can harvest the ic values in the arrays for setting the parval1 values
    pyemu.os_utils.run("python forward_run.py",cwd=t_d)

    # now write ins and tpl file for these
    ic_parvals = {}
    obs_to_par_map = dict()
    for k in range(3):
        fname = os.path.join(t_d,"heads_{0}.dat".format(k))
        assert os.path.exists(fname),fname
        arr = np.loadtxt(fname)
        fname_ins = fname + "_out.ins"
        fname_tpl = fname + "_in.tpl"
        in_arr = np.loadtxt(fname_tpl.replace(".tpl",""))
        ft = open(fname_tpl,'w')
        with open(fname_ins,'w') as f:
            f.write("pif ~\n")
            ft.write("ptf ~\n")
            for i in range(arr.shape[0]):
                f.write("l1 ")
                for j in range(arr.shape[1]):
                    if np.abs(arr[i,j]) > 100 or np.abs(in_arr[i,j]) > 100:
                        f.write(" !dum! ")
                        ft.write(" 40 ")
                    else:
                        oname = "head_{0:02d}_{1:03d}_{2:03d}".format(k,i,j)
                        f.write(" !{0}! ".format(oname))
                        if sync_state_names:
                            ft.write(" ~  {0} ~ ".format(oname))
                            ic_parvals[oname] = in_arr[i, j]
                        else:
                            pname = "p"+oname
                            ft.write(" ~  {0} ~ ".format(pname))
                            obs_to_par_map[oname] = pname
                            ic_parvals[pname] = in_arr[i, j]

                f.write("\n")
                ft.write("\n")
        ft.close()
        new_tpl.append(fname_tpl)
        new_in.append(fname_tpl.replace(".tpl",""))
        new_tpl_cycle.append(-1)
        new_ins.append(fname_ins)
        new_out.append(os.path.join(t_d,"heads_{0}.dat".format(k)))
        new_ins_cycle.append(-1)

        i = pyemu.pst_utils.InstructionFile(fname_ins)
        df = i.read_output_file(fname)
        #print(df)

    # split out the wel and rch tpl files into cycle files
    lines = []
    with open(os.path.join(t_d,"freyberg6.wel.tpl"),'r') as f:
        for i in range(19):
            lines.append(f.readline())
    print(lines)

    for icycle in range(25):
        tpl_file = os.path.join(t_d,"freyberg6.wel_{0}.tpl".format(icycle))
        with open(tpl_file,'w') as f:
            for line in lines:
                new_line = line.replace("_0","_{0}".format(icycle))
                f.write(new_line)

        new_tpl.append(tpl_file)
        new_in.append(os.path.join(t_d,"freyberg6.wel"))
        new_tpl_cycle.append(icycle)
    remove_tpl = ["freyberg6.wel.tpl"]

    lines = []
    with open(os.path.join(t_d, "freyberg6.rch.tpl"), 'r') as f:
        for i in range(11):
            lines.append(f.readline())
    print(lines)

    for icycle in range(25):
        tpl_file = os.path.join(t_d,"freyberg6.rch_{0}.tpl".format(icycle))
        with open(tpl_file,'w') as f:
            for line in lines:
                new_line = line.replace("_0","_{0}".format(icycle))
                f.write(new_line)

        new_tpl.append(tpl_file)
        new_in.append(os.path.join(t_d,"freyberg6.rch"))
        new_tpl_cycle.append(icycle)
    remove_tpl.append('freyberg6.rch.tpl')

    # now for the fun part: modify the pst
    pst = pyemu.Pst(os.path.join(t_d,"freyberg6_run.pst"))

    print(pst.npar_adj,pst.nnz_obs)

    # swap out obs info
    dropped_dfs = []
    for ins in remove_ins:
        dropped_dfs.append(pst.drop_observations(os.path.join(t_d, ins), '.'))
    for insf, outf, cy in zip(new_ins, new_out, new_ins_cycle):
        df = pst.add_observations(insf,outf, pst_path=".")
        pst.observation_data.loc[df.obsnme, "cycle"] = cy
        pst.model_output_data.loc[os.path.join(".",os.path.split(insf)[1]),"cycle"] = cy
    pst.observation_data.loc[:,"weight"] = 0.0
    for df in dropped_dfs:
        for c in ["obsval","weight"]:
            pst.observation_data.loc[df.obsnme, c] = df.loc[:, c]

    # swap out par info
    dropped_dfs = []
    for tpl in remove_tpl:
        dropped_dfs.append(pst.drop_parameters(os.path.join(t_d,tpl),'.'))
    pst.parameter_data.loc[:,"cycle"] = -1
    pst.model_input_data.loc[:,"cycle"] = -1
    for tplf, inf, cy in zip(new_tpl,new_in,new_tpl_cycle):
        df = pst.add_parameters(tplf,inf,pst_path=".")
        pst.parameter_data.loc[df.parnme,"cycle"] = cy
        pst.model_input_data.loc[os.path.join(".",os.path.split(tplf)[1]),"cycle"] = cy
    for df in dropped_dfs:
        for c in ["parval1","parubnd","parlbnd","pargp"]:
            pst.parameter_data.loc[df.parnme,c] = df.loc[:,c]

    #set the state parameter info
    for p,v in ic_parvals.items():
        pst.parameter_data.loc[p,"parval1"] = v
        pst.parameter_data.loc[p, "parlbnd"] = v * 0.5
        pst.parameter_data.loc[p, "parubnd"] = v * 1.5
        pst.parameter_data.loc[p,"pargp"] = "head_state"
        pst.parameter_data.loc[p,"partrans"] = "none"

    pst.control_data.noptmax = 2
    pst.model_command = "python forward_run.py"
    pst.pestpp_options.pop("ies_par_en")
    pst.parameter_data.loc["perlen","partrans"] = "fixed"
    pst.pestpp_options["ies_num_reals"] = 5
    pst.pestpp_options["da_num_reals"] = 5
    if not sync_state_names:
        pst.observation_data.loc[:,"state_par_link"] = np.NaN
        obs = pst.observation_data
        obs.loc[:,"state_par_link"] = obs.obsnme.apply(lambda x: obs_to_par_map.get(x,np.NaN))
    pst.write(os.path.join(t_d,"freyberg6_run_da1.pst"),version=2)
    return pst
    
def process_complex_target_output(c_d, s_d, d_d, real):
    
    # process for BAT
    redis_fac = 3
    
    start_date = pd.to_datetime('20151231', format = '%Y%m%d')
    
    # load in obs ensemble
    oe_f = pd.read_csv(os.path.join(c_d, "freyberg.0.obs.csv"), index_col=0)
    oe_f = oe_f.T
    
    #process obs
    hds_f = oe_f.loc[oe_f.index.to_series().apply(lambda x: x.startswith("hds")), :].copy()
    hds_f.loc[:, "k"] = hds_f.index.to_series().apply(lambda x: int(x.split('_')[2]))
    hds_f.loc[:, "i"] = hds_f.index.to_series().apply(lambda x: int(x.split('_')[3]))
    hds_f.loc[:, "j"] = hds_f.index.to_series().apply(lambda x: int(x.split('_')[4]))
    hds_f.loc[:, "time"] = hds_f.index.to_series().apply(lambda x: float(x.split('_')[-1].split(':')[1]))
    time = hds_f.loc[:, "time"]
    time = pd.to_timedelta(time.values-1, unit='D')
    hds_f.loc[:, "org_time"] = start_date + time.values
    hds_f.loc[:, "org_time"] = hds_f.loc[:, "org_time"].apply(lambda x: x.strftime('%Y%m%d'))
    hds_f.loc[:, "org_i"] = (hds_f.i / redis_fac).apply(np.int)
    hds_f.loc[:, "org_j"] = (hds_f.j / redis_fac).apply(np.int)
    hds_f.loc[:, "org_obgnme"] = hds_f.apply(lambda x: "trgw_{0}_{1}_{2}_{3}".format(x.k, x.org_i, x.org_j, x.org_time), axis=1)
    
    sfr_f = oe_f.loc[oe_f.index.to_series().apply(lambda x: x.startswith("sfr")), :].copy()
    sfr_f.loc[:, "time"] = sfr_f.index.to_series().apply(lambda x: float(x.split('_')[-1].split(':')[1]))
    type = sfr_f.index.to_series().apply(lambda x: x.split(':')[1].split('_')[0])
    type = pd.DataFrame(type)
    type = type.replace('gage', 'gage_1')
    sfr_f.loc[:, "type"] = type.values
    time = sfr_f.loc[:, "time"]
    time = pd.to_timedelta(time.values-1, unit='D')
    sfr_f.loc[:, "org_time"] = start_date + time.values
    sfr_f.loc[:, "org_time"] = sfr_f.loc[:, "org_time"].apply(lambda x: x.strftime('%Y%m%d'))
    sfr_f.loc[:, "org_obgnme"] = sfr_f.apply(lambda x: "{0}_{1}".format(x.type, x.org_time), axis=1)
    
    pst = pyemu.Pst(os.path.join(s_d, 'freyberg6_run_ies.pst'))
    obs_s = pst.observation_data
    
    for j in range(len(hds_f)):
        for i in range(len(obs_s)):
            if obs_s.obsnme[i] in hds_f.org_obgnme[j]:
                obs_s.obsval[i] = hds_f.iloc[j, real]

    for j in range(len(sfr_f)):
        for i in range(len(obs_s)):
            if obs_s.obsnme[i] in sfr_f.org_obgnme[j]:
                obs_s.obsval[i] = sfr_f.iloc[j, real]
                
    #write pst files to template ies dir (for current real)
    m_ies_dir = os.path.join('simple_template_ies_{0}'.format(real)) 
    shutil.copytree(s_d,m_ies_dir)         
    pst.write(os.path.join(m_ies_dir,"freyberg6_run_ies.pst"),version=2)
    
    # process for SEQ
    redis_fac = 3
    pst = pyemu.Pst(os.path.join(d_d, "freyberg6_run_da2.pst"))
    
    start_date = pd.to_datetime('20151231', format='%Y%m%d')
    dates = pd.date_range(start='2015-12-31', periods=25,freq='M').strftime("%Y%m%d")
    
    # load in obs cycle table
    oe_f = pd.read_csv(os.path.join(c_d, "freyberg.0.obs.csv"), index_col=0)
    oe_f = oe_f.T
    obs_d = pd.read_csv(os.path.join(d_d, 'obs_cycle_tbl.csv'))
    obs_d.loc["date",1:] = dates
    
    hds_f = oe_f.loc[oe_f.index.to_series().apply(lambda x: x.startswith("hds")), :].copy()
    hds_f.loc[:, "k"] = hds_f.index.to_series().apply(lambda x: int(x.split('_')[2]))
    hds_f.loc[:, "i"] = hds_f.index.to_series().apply(lambda x: int(x.split('_')[3]))
    hds_f.loc[:, "j"] = hds_f.index.to_series().apply(lambda x: int(x.split('_')[4]))
    hds_f.loc[:, "time"] = hds_f.index.to_series().apply(lambda x: float(x.split('_')[-1].split(':')[1]))
    time = hds_f.loc[:, "time"]
    time = pd.to_timedelta(time.values-1, unit='D')
    hds_f.loc[:, "org_time"] = start_date + time.values
    hds_f.loc[:, "org_time"] = hds_f.loc[:, "org_time"].apply(lambda x: x.strftime('%Y%m%d'))
    hds_f.loc[:, "org_i"] = (hds_f.i / redis_fac).apply(np.int)
    hds_f.loc[:, "org_j"] = (hds_f.j / redis_fac).apply(np.int)
    hds_f.loc[:, "org_obgnme"] = hds_f.apply(lambda x: "head_{:02d}_{:03d}_{:03d}".format(x.k, x.org_i, x.org_j), axis=1)
    
    sfr_f = oe_f.loc[oe_f.index.to_series().apply(lambda x: x.startswith("sfr")), :].copy()
    sfr_f.loc[:, "time"] = sfr_f.index.to_series().apply(lambda x: float(x.split('_')[-1].split(':')[1]))
    type = sfr_f.index.to_series().apply(lambda x: x.split(':')[1].split('_')[0])
    type = pd.DataFrame(type)
    type = type.replace('gage', 'gage_1')
    sfr_f.loc[:, "type"] = type.values
    time = sfr_f.loc[:, "time"]
    time = pd.to_timedelta(time.values-1, unit='D')
    sfr_f.loc[:, "org_time"] = start_date + time.values
    sfr_f.loc[:, "org_time"] = sfr_f.loc[:, "org_time"].apply(lambda x: x.strftime('%Y%m%d'))
    sfr_f.loc[:, "org_obgnme"] = sfr_f.apply(lambda x: "{0}".format(x.type), axis=1)
    
    for i in range(25):
        for j in range(len(hds_f)):
            for k in range(2):
                if obs_d.iloc[k,0] in hds_f.org_obgnme[j] and obs_d.iloc[3,i] == hds_f.org_time[j]:
                    obs_d.iloc[k,i] = hds_f.iloc[j, real]
    
    for i in range(25):
        for j in range(len(sfr_f)):
                if obs_d.iloc[2,0] in sfr_f.org_obgnme[j] and obs_d.iloc[3,i] == sfr_f.org_time[j]:
                    obs_d.iloc[2,i] = sfr_f.iloc[j, real]
    obs_d = obs_d.drop(['date'], axis = 0)
    
    #write pst files and obs cycle table to master da dir (for current real)
    m_da_dir = os.path.join('simple_template_da_{0}'.format(real)) 
    shutil.copytree(d_d,m_da_dir)         
    pst.write(os.path.join(m_da_dir,"freyberg6_run_da.pst"),version=2)
    obs_d.to_csv(os.path.join(m_da_dir,'obs_cycle_tbl.csv'), index=False)

    
    
def balance_weights(ireal):

    ies_dir = os.path.join('simple_template_ies_{0}'.format(ireal))    
    da_dir =  os.path.join('simple_template_da_{0}'.format(ireal))   
    ies_file = 'freyberg6_run_ies.pst'
    
    #run pestpp ies to get phi
    pyemu.os_utils.run("pestpp-ies.exe {0}".format(ies_file),cwd=ies_dir)
    
    pst = pyemu.Pst(os.path.join(ies_dir, ies_file))
    obs = pst.observation_data
    
    #balance weights based on phi components
    pst.res.loc[:,'residual'] = pst.res.loc[:,'residual'] * pst.res.loc[:,'weight']
    phi_comps = pst.res.groupby(["group"]).sum()
    phi_comp = phi_comps.loc[:,'residual']
    phi_grps = phi_comps.index
    obs_phi_dict = dict(zip(phi_grps, phi_comp))
    pst._adjust_weights_by_phi_components(obs_phi_dict, original_ceiling = False)
    pst.write(os.path.join(ies_dir, ies_file), version = 2)
    
    #modify da weight cycle table based on new weights 
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
    dates.loc[:,'cycle']=range(25)
    cycle_dict = dict(zip(dates.iloc[:,0],dates.loc[:,'cycle']))
    hds_f.loc[:,'time'] = hds_f.loc[:,'time'].replace(cycle_dict, regex=True)
    sfr_f.loc[:,'time'] = sfr_f.loc[:,'time'].replace(cycle_dict, regex=True)
    
    for i in range(3):
        for j in range(24):
            for k in range(650):
                m = j+2
                if da_wt.iloc[i,0] in hds_f.iloc[k,9] and hds_f.iloc[k,5] == m:
                    da_wt.iloc[i,m] = hds_f.iloc[k, 2]
    
    for j in range(24):
        for k in range(25):
            m = j+2
            if sfr_f.iloc[k,5] == m:
                da_wt.iloc[2,m] = sfr_f.iloc[k, 2]
    
    da_wt.to_csv(os.path.join(da_dir, 'weight_cycle_tbl.csv'), index=False)

     
def compare_mf6_freyberg():
    for ireal in range(100):
        complex_dir = os.path.join('complex_master')
        ies_dir = os.path.join('simple_template_ies')
        da_dir = os.path.join('simple_template_da')
        
        process_complex_target_output(complex_dir, ies_dir, da_dir, ireal)
        
        balance_weights(ireal)
        
        #run batch and sequential simple models 
        #ies stuff 
        ies_t_d = os.path.join('simple_template_ies_{0}'.format(ireal))
        ies_pst = pyemu.Pst(os.path.join(ies_t_d,"freyberg6_run_ies.pst"))
        
        # prep that prior ensemble for da
        da_t_d = os.path.join('simple_template_da_{0}'.format(ireal))
        da_pst = pyemu.Pst(os.path.join(da_t_d,"freyberg6_run_da.pst"))
        par = da_pst.parameter_data
        par.loc[par.parnme.str.contains("welflx"),"scale"] = -1.0
    
        ies_pst.pestpp_options["ies_par_en"] = "prior.jcb"
        ies_pe = pyemu.ParameterEnsemble.from_binary(pst=ies_pst,filename=os.path.join(ies_t_d,"prior.jcb"))
        da_pe = pyemu.ParameterEnsemble.from_gaussian_draw(pst=da_pst,
            cov=pyemu.Cov.from_parameter_data(da_pst),num_reals=ies_pe.shape[0])
        da_pe.index = ies_pe.index
        d = set(da_pe.columns.tolist()).symmetric_difference(set(ies_pe.columns.tolist()))
        print(d)
        da_pe.loc[:,ies_pe.columns] = ies_pe.values
        da_pe.to_binary(os.path.join(da_t_d,"da_prior.jcb"))
        da_pst.pestpp_options["ies_par_en"] = "da_prior.jcb"
        
        # set pestpp options for batch da
        ies_pst.pestpp_options.pop("ies_num_reals",None)
        ies_pst.pestpp_options.pop("da_num_reals",None)
        ies_pst.pestpp_options["ies_no_noise"] = False
        ies_pst.pestpp_options["ies_verbose_level"] = 1
        ies_pst.pestpp_options.pop("ies_localizer",None)
        ies_pst.pestpp_options["ies_autoadaloc"] = False
        ies_pst.pestpp_options["ies_save_lambda_en"] = False
        ies_pst.pestpp_options["ies_drop_conflicts"] = False
        ies_pst.pestpp_options["ies_num_reals"] = 50
        ies_pst.pestpp_options["ies_use_mda"] = False
        ies_pst.control_data.noptmax = 3
        ies_pst.write(os.path.join(ies_t_d,"freyberg6_run_ies.pst"),version=2)
        
        # set pestpp options for sequential da
        da_pst.pestpp_options.pop("da_num_reals",None)
        da_pst.pestpp_options.pop("ies_num_reals",None)
        da_pst.pestpp_options["ies_no_noise"] = False
        da_pst.pestpp_options["ies_verbose_level"] = 1
        da_pst.pestpp_options.pop("ies_localizer", None)
        da_pst.pestpp_options["ies_autoadaloc"] = False
        da_pst.pestpp_options["ies_save_lambda_en"] = False
        da_pst.pestpp_options["ies_drop_conflicts"] = False
        da_pst.pestpp_options["ies_num_reals"] = 50
        da_pst.pestpp_options["ies_use_mda"] = False
        da_pst.control_data.noptmax = 3
        da_pst.write(os.path.join(da_t_d,"freyberg6_run_da.pst"),version=2)
    
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
    pyemu.os_utils.start_workers(c_t,"pestpp-ies","freyberg.pst",num_workers=4,worker_root=".",
                                 master_dir=c_t.replace("template","master"))
    
def plot_phi_seq_bat():
    seq_phi_master = []
    bat_phi_master = []
    bat_phi_master = pd.DataFrame(bat_phi_master)

    for i in range(100):
        seq_dir = os.path.join('simple_master_da_{0}'.format(i))
        bat_dir = os.path.join('simple_master_ies_{0}'.format(i))

        bat_phi = pd.read_csv(os.path.join(bat_dir, 'freyberg6_run_ies.phi.actual.csv'))
        bat_phi_master = bat_phi_master.append(bat_phi.iloc[3,:])
        # print(bat_phi.iloc[3,:])

        seq_phi = pd.read_csv(os.path.join(seq_dir, 'freyberg6_run_da.global.phi.actual.csv'))
        seq_phi = seq_phi.loc[seq_phi.iteration.max(),'mean'].sum()
        seq_phi_master.append(seq_phi)

        # print(seq_phi)

    plt.hist([bat_phi_master.iloc[:,3], seq_phi_master],  label=['BAT', 'SEQ'])
    plt.legend(loc='upper right')
    plt.title("mean phi for batch data assimilation of simple model across 100 complex reals")
    plt.show()

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
        complex_tail = complex_obs.loc[i,'tailwater_time:20170131']
        complex_head = complex_obs.loc[i,'headwater_time:20171031']
        complex_gw = complex_obs.loc[i,'trgw_0_29_5_time:20171031']
        complex_gw_all.append(complex_gw)
        complex_head_all.append(complex_head)
        complex_tail_all.append(complex_tail)

        simple_bat_dir = os.path.join('simple_master_ies_{0}'.format(i))
        simple_bat = pd.read_csv(os.path.join(simple_bat_dir, 'freyberg6_run_ies.3.obs.csv'))
        simple_bat_gw = simple_bat.loc[:,'trgw_0_9_1_20171031']
        simple_bat_tail = simple_bat.loc[:,'tailwater_20170131']
        simple_bat_head = simple_bat.loc[:,'headwater_20171031']
        simple_bat_tail_all.append(simple_bat_tail)
        simple_bat_head_all.append(simple_bat_head)
        simple_bat_gw_all.append(simple_bat_gw)

        simple_seq_dir = os.path.join('simple_master_da_{0}'.format(i))
        simple_seq = pd.read_csv(os.path.join(simple_seq_dir, 'freyberg6_run_da.23.0.obs.csv'))
        simple_seq_gw = simple_seq.loc[:,'head_00_009_001']
        simple_seq_head = simple_seq.loc[:,'headwater']
        simple_seq_head_all.append(simple_seq_head)
        simple_seq_gw_all.append(simple_seq_gw)

        simple_seq = pd.read_csv(os.path.join(simple_seq_dir, 'freyberg6_run_da.14.0.obs.csv'))
        simple_seq_tail = simple_seq.loc[:,'tailwater']
        simple_seq_tail_all.append(simple_seq_tail)

    for ye, xe in zip(complex_gw_all, simple_bat_gw_all):
        plt.scatter(xe, [ye] * len(xe), color='blue', s=1, label='BAT')
    for ze, le in zip(complex_gw_all, simple_seq_gw_all):
        plt.scatter(le, [ze] * len(le), color='orange', s = 1,label = 'SEQ')
    plt.xlim(33.5, 38.5)
    plt.ylim(33.5, 38.5)
    plt.title('GW_3 Forecast')
    plt.xlabel('Simple Forecast (ft)')
    plt.ylabel('Complex Forecast (ft)')
    plt.legend(loc = 'upper right')
    plt.tight_layout()
    plt.savefig('gw_3_forecast.pdf')
    plt.close()
    # plt.show()

    for ye, xe in zip(complex_head_all, simple_bat_head_all):
        plt.scatter(xe, [ye] * len(xe), color='blue', s=1, label='BAT')
    for ze, le in zip(complex_head_all, simple_seq_head_all):
        plt.scatter(le, [ze] * len(le), color='orange', s = 1,label = 'SEQ')
    plt.xlim(-1500,500)
    plt.ylim(-1500,500)
    plt.title('Headwater Forecast')
    plt.xlabel('Simple Forecast (ft)')
    plt.ylabel('Complex Forecast (ft)')
    plt.legend(loc = 'upper right')
    plt.tight_layout()
    plt.savefig('headwater_forecast.pdf')
    plt.close()
    # plt.show()

    for ye, xe in zip(complex_tail_all, simple_bat_tail_all):
        plt.scatter(xe, [ye] * len(xe), color='blue', s=1, label='BAT')
    for ze, le in zip(complex_tail_all, simple_seq_tail_all):
        plt.scatter(le, [ze] * len(le), color='orange', s = 1,label = 'SEQ')
    plt.xlim(-2500,0)
    plt.ylim(-2500,0)
    plt.title('Tailwater Forecast')
    plt.xlabel('Simple Forecast (ft)')
    plt.ylabel('Complex Forecast (ft)')
    plt.legend(loc = 'upper right')
    plt.tight_layout()
    plt.savefig('tailwater_forecast.pdf')
    plt.close()
    # plt.show()

if __name__ == "__main__":
    
    # BOOLEANS TO SELECT CODE BLOCKS BELOW (currently a wishlist)
    prep_complex_model = False #do this once before running paired simple/complex analysis
    run_prior_mc = False
    run_simple_complex = False
    plot_s_vs_s = True
    plot_phis = False
    plot_phi_diffs = False
    
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
        


