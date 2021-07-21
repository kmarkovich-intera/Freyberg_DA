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


bin_path = os.path.join("test_bin")
if "linux" in platform.platform().lower():
    bin_path = os.path.join(bin_path,"linux")
elif "darwin" in platform.platform().lower():
    bin_path = os.path.join(bin_path,"mac")
else:
    bin_path = os.path.join(bin_path,"win")

bin_path = os.path.abspath("test_bin")
os.environ["PATH"] += os.pathsep + bin_path


bin_path = os.path.join("..","..","..","bin")
exe = ""
if "windows" in platform.platform().lower():
    exe = ".exe"
exe_path = os.path.join(bin_path, "pestpp-da" + exe)

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

    # mod the sto to make sp 1 transient - or should this be a pre-processor so that cycle 0 is 
    # ss and the rest are transient?

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
    
def process_complex_target_output(c_d, s_d, real):

    sim = flopy.mf6.MFSimulation.load(sim_ws=c_d)
    m = sim.get_model("freyberg6")
    redis_fac = m.dis.nrow.data / 40
    pst = pyemu.Pst(os.path.join(c_d,"freyberg.pst"))
    
    start_date = pd.to_datetime('20151231', format = '%Y%m%d')
    
    # load in obs ensemble
    oe_f = pd.read_csv(os.path.join(c_d, "freyberg.0.obs.csv"), index_col=0)
    oe_f = oe_f.T
    obs_f = pst.observation_data
    
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
                
    #write pst files to master ies dir
    m_ies_dir = os.path.join('simple_master_ies_{0}'.format(real))          
    pst.write(os.path.join(m_ies_dir,"freyberg6_run_ies.pst"),version=2)
     
def compare_mf6_freyberg():
    
    for ireal in range(100):
        
        process_complex_target_output(complex_dir, simple_dir, ireal)
        
        #run batch and sequential simple models 
        # prep that prior ensemble for da
        da_test_d = "mf6_freyberg"
        da_t_d = os.path.join(da_test_d, "template_seq_native")
        org_da_t_d = da_t_d
        da_t_d = da_t_d + "_compare"
        if os.path.exists(da_t_d):
            shutil.rmtree(da_t_d)
        shutil.copytree(org_da_t_d,da_t_d)
        da_pst = pyemu.Pst(os.path.join(da_t_d,"freyberg6_run_da2.pst"))
        par = da_pst.parameter_data
        par.loc[par.parnme.str.contains("welflx"),"scale"] = -1.0
        par.loc["perlen","parval1"] = 100000
        #par.loc[~par.parnme.str.contains("head"),"partrans"] = "fixed"
          
        m_ies_dir = os.path.join('simple_master_ies_{0}'.format(real)) 
        org_dir
        if os.path.exists(ies_t_d):
            shutil.rmtree(ies_t_d)
        shutil.copytree(org_ies_t_d,ies_t_d)
        ies_pst = pyemu.Pst(os.path.join(ies_t_d,"freyberg6_run_ies.pst"))
    
        mod_tdis_sto(org_ies_t_d,ies_t_d)
        
        pyemu.os_utils.run("mf6",cwd=ies_t_d)
        pyemu.os_utils.run("mf6",cwd=da_t_d)
    
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
    
        # load in complex model target outputs as obs
        ies_obs = ies_pst.observation_data
        da_obs = da_pst.observation_data
            
        #ies_obs.loc[[o for o in ies_pst.nnz_obs_names if not "gage" in o],"weight"] = 50
        #ies_obs.loc[[o for o in ies_pst.nnz_obs_names if "gage" in o],"weight"] = 1.0 / (ies_obs.loc[[o for o in ies_pst.nnz_obs_names if "gage" in o],"obsval"] * 0.15)
    
        # check the phi contribs
        # ies_pst.control_data.noptmax = 0
        # ies_pst.write(os.path.join(ies_t_d,"reweight.pst"))
        # pyemu.os_utils.run("{0} reweight.pst".format(exe_path.replace("-da","-ies")),cwd=ies_t_d)
        # ies_pst.set_res(os.path.join(ies_t_d,"reweight.base.rei"))
        # print({n:p for n,p in ies_pst.phi_components.items() if n in ies_pst.nnz_obs_groups})
        # ies_pst.plot(kind="phi_pie")
        # import matplotlib.pyplot as plt
        # plt.show()
        # return
    
        tr_obs.loc[:,"k"] = tr_obs.obsnme.apply(lambda x: int(x.split('_')[1]))
        tr_obs.loc[:, "i"] = tr_obs.obsnme.apply(lambda x: int(x.split('_')[2]))
        tr_obs.loc[:, "j"] = tr_obs.obsnme.apply(lambda x: int(x.split('_')[3]))
        tr_obs.loc[:, "kij"] = tr_obs.apply(lambda x: (x.k,x.i,x.j),axis=1)
        ukij = set(tr_obs.kij.unique().tolist())
        da_obs = da_pst.observation_data
        hd_obs = da_obs.loc[da_obs.obsnme.str.startswith("head_"),:].copy()
        print(hd_obs.obsnme)
        hd_obs.loc[:,"k"]= hd_obs.obsnme.apply(lambda x: int(x.split('_')[1]))
        hd_obs.loc[:, "i"] = hd_obs.obsnme.apply(lambda x: int(x.split('_')[2]))
        hd_obs.loc[:, "j"] = hd_obs.obsnme.apply(lambda x: int(x.split('_')[3]))
        hd_obs.loc[:,"kij"] = hd_obs.apply(lambda x: (x.k,x.i,x.j),axis=1)
        hd_obs.loc[:,"org_obgnme"] = hd_obs.apply(lambda x: "trgw_{0}_{1}_{2}".format(x.k,x.i,x.j),axis=1)
    
        include = hd_obs.loc[hd_obs.kij.apply(lambda x: x in ukij),"obsnme"].tolist()
        #include.append("gage_1")
        otbl = pd.DataFrame(columns=include,index=np.arange(25)).T
        wtbl = pd.DataFrame(columns=include,index=np.arange(25)).T
        for uog in hd_obs.loc[include,"org_obgnme"].unique():
            oname = hd_obs.loc[hd_obs.org_obgnme==uog,"obsnme"].values[0]
            og_obs = ies_obs.loc[ies_obs.obgnme==uog,:]
            og_obs.sort_index(inplace=True)
            otbl.loc[oname,:] = og_obs.obsval.values
            wtbl.loc[oname,:] = og_obs.weight.values
    
        otbl.loc["gage_1",:] = ies_obs.loc[ies_obs.obsnme.str.contains("gage"),"obsval"].values
        wtbl.loc["gage_1",:] = ies_obs.loc[ies_obs.obsnme.str.contains("gage"), "weight"].values
        otbl.to_csv(os.path.join(da_t_d,"obs_cycle_tbl.csv"))
        wtbl.to_csv(os.path.join(da_t_d, "weight_cycle_tbl.csv"))
        da_pst.observation_data.loc[otbl.index.values,"weight"] = 1.0
        
        # set pestpp options for batch da
        ies_pst.pestpp_options.pop("ies_num_reals",None)
        ies_pst.pestpp_options.pop("da_num_reals",None)
        ies_pst.pestpp_options["ies_no_noise"] = False
        ies_pst.pestpp_options["ies_verbose_level"] = 1
        ies_pst.pestpp_options.pop("ies_localizer",None)
        ies_pst.pestpp_options["ies_autoadaloc"] = False
        ies_pst.pestpp_options["ies_save_lambda_en"] = False
        ies_pst.pestpp_options["ies_drop_conflicts"] = True
        ies_pst.pestpp_options["ies_num_reals"] = 50
        
        # set pestpp options for sequential da
        da_pst.pestpp_options.pop("da_num_reals",None)
        da_pst.pestpp_options.pop("ies_num_reals",None)
        da_pst.pestpp_options["ies_no_noise"] = False
        da_pst.pestpp_options["ies_verbose_level"] = 1
        da_pst.pestpp_options.pop("ies_localizer", None)
        da_pst.pestpp_options["ies_autoadaloc"] = False
        da_pst.pestpp_options["ies_save_lambda_en"] = False
        da_pst.pestpp_options["ies_drop_conflicts"] = True
        da_pst.pestpp_options["ies_num_reals"] = 50
        #da_pst.pestpp_options["da_stop_cycle"] = 1
    
        ies_pst.control_data.noptmax = 3
        da_pst.control_data.noptmax = 3
    
        # run da                     
        da_pst.pestpp_options["ies_use_mda"] = False
        da_pst.write(os.path.join(da_t_d,"freyberg6_run_da.pst"),version=2)
        da_m_d_glm = os.path.join(da_test_d, "master_da_glm_{0}".format(i))
        pyemu.os_utils.start_workers(da_t_d, exe_path.replace("-ies","-da"), "freyberg6_run_da.pst",
                                    num_workers=20, worker_root=da_test_d, port=port,
                                    master_dir=da_m_d_glm, verbose=True)
    
        # run ies    
        ies_pst.pestpp_options["ies_use_mda"] = False
        ies_pst.write(os.path.join(ies_t_d,"freyberg6_run_ies.pst"),version=2)
        ies_m_d_glm = os.path.join(ies_test_d, "master_ies_glm_{0}".format(i))
        pyemu.os_utils.start_workers(ies_t_d, exe_path.replace("-da","-ies"), "freyberg6_run_ies.pst",
                                    num_workers=20, worker_root=ies_test_d, port=port,
                                    master_dir=ies_m_d_glm, verbose=True)
                                
def run_complex_prior_mc(c_t):
    pyemu.os_utils.start_workers(c_t,"pestpp-ies","freyberg.pst",num_workers=4,worker_root=".",
                                 master_dir=c_t.replace("template","master"))
    
    

if __name__ == "__main__":
    
    # BOOLEANS TO SELECT CODE BLOCKS BELOW (currently a wishlist)
    prep_complex_model = False #do this once before running paired simple/complex analysis
    run_prior_mc = False
    run_simple_complex = False
    plot_s_vs_s = False


    
    if prep_complex_model:
        prep_complex_prior_mc()
        
    if run_prior_mc:
        run_complex_prior_mc('complex_template')
        
    if run_simple_complex:
        compare_mf6_freyberg()
        
        
#     if plot_s_vs_s:
        


