import os.path
import os
import glob
import pandas as pd
import numpy as np

## Global variables
BASE_SEARCH_ORDER = [
    "/ghds/groups/lunar/flowcells/processed/lunar1v1.0_rc6/",
    "/ghds/groups/lunar/flowcells/processed/lunar1v1.0_rc5/",
    "/ghds/groups/lunar/flowcells/processed/lunar1v1.0_rc4/",
    "/ghds/groups/lunar/flowcells/processed/lunar1v1.0_rc3/",
    "/ghds/groups/lunar/flowcells/processed/lunar1v1.0_rc2/",
    "/ghds/groups/lunar/flowcells/processed/lunar1v1.0_rc1/",
    "/ghds/flowcentral/",
    "/ghds/groups/lunar/flowcells/processed/bip38_33/",
    "/ghds/groups/lunar/flowcells/processed/bip38_29/",
    "/ghds/groups/lunar/flowcells/processed/bip38_25/",
    "/ghds/groups/lunar/flowcells/processed/bip38/",
    "/mnt/flowcentral/"]

OUTPUT_PATH = "/home/ykong/notebook/qcdashboard/output_test"

## FILE PROCESSING FUNCTION

## 1. Give a flowcell, find the most recent output folder
def get_fc_path(fc):
    for path in BASE_SEARCH_ORDER:
        fc_name = glob.glob(path + fc + "*")
        if fc_name:
            return fc_name[0]

    ## raise error message if the flowcell not found
    raise FileNotFoundError("The flowcell {} can not be found in all searched dirs".format(fc))


## 2. Given one or more flowcells, 0 or more sample names, return a df with corresponding path
def process_fc_sample(fcs, samples=""):
    sample_df = pd.DataFrame()

    ## check the input is one fc string or multiple fcs in a list
    if isinstance(fcs, str):
        fcs = [fcs]
    ## find the most recent lunar output path for each input fc
    for fc in fcs:
        tmp_dir = get_fc_path(fc)

        # list all the run sample ids on this fc
        all_samples = glob.glob(tmp_dir + "/*.mbdqc_report.hdr.tsv")
        all_samples = list(map(lambda x: x.split("/")[-1].split(".")[0], all_samples))
        sub_sample_df = pd.DataFrame(
            {"fc": fc, "path": tmp_dir, "run_sample_id": all_samples}
        )
        sample_df = sample_df.append(sub_sample_df)

    ## filter samples to keep only input samples
    ## if no sample name supplied, return all samples on the fc
    if samples:
        sample_df = sample_df[sample_df["run_sample_id"].isin(samples)]

    ## check if there is any input samples not on the input flowcells
    sample_not_on = set(samples) - set(sample_df["run_sample_id"])
    if sample_not_on:
        raise FileNotFoundError(
            "Sample: {} can not be found in any of the provided flowcells".format(
                sample_not_on
            )
        )

    ## check if output is empty
    if sample_df.shape[0] < 1:
        raise FileNotFoundError(
            "None of the flowcells provided can be found in all searched dirs"
        )
    else:
        sample_df = sample_df.reset_index().drop("index", axis=1)
        return sample_df

## 3. Given fcs, sample_ids, generate df with reportable MBD QC metrics
def get_mbd_qc_report(fcs, samples="", subpanel="backbone",write=True):

    ## get the df containing paths to fcs and samples
    sample_df = process_fc_sample(fcs, samples)

    mbd_qc_report_df = pd.DataFrame()

    for index, row in sample_df.iterrows():
        qc_report = pd.read_csv(os.path.join(row["path"], row["run_sample_id"] + ".mbdqc_report.hdr.tsv"), sep="\t")

        ## 1. 0 CG in hyper frac
        cg_zero_hyper = qc_report[(qc_report["subpanel"] == subpanel) & \
                                  (qc_report["mbdqc_metric"] == "0cg_molecule_hyper_pct")]["value"].iloc[0]

        ## 2. unmethylated in hyper frac
        unmeth_hyper = qc_report[(qc_report["subpanel"] == subpanel) & \
                                 (qc_report["mbdqc_metric"] == "unmeth_hyper_pct")]["value"].iloc[0]

        ## 3. methy 1/2
        methy_1_2 = qc_report[(qc_report["subpanel"] == subpanel) & \
                              (qc_report["mbdqc_metric"] == "methyl_1_2")]["value"].iloc[0]

        ## 4. MBD fraction distribution
        hyper_frac = qc_report[(qc_report["subpanel"] == subpanel) & \
                               (qc_report["mbdqc_metric"] == "hyper_molecule_pct")]["value"].iloc[0]
        hypo_frac = qc_report[(qc_report["subpanel"] == subpanel) & \
                              (qc_report["mbdqc_metric"] == "hypo_molecule_pct")]["value"].iloc[0]
        residual_frac = qc_report[(qc_report["subpanel"] == subpanel) & \
                                  (qc_report["mbdqc_metric"] == "residual_molecule_pct")]["value"].iloc[0]
        mixed_frac = qc_report[(qc_report["subpanel"] == subpanel) & \
                               (qc_report["mbdqc_metric"] == "mixed_molecule_pct")]["value"].iloc[0]

        ## 5. methy hypo pct
        methy_hypo = qc_report[(qc_report["subpanel"] == subpanel) & \
                               (qc_report["mbdqc_metric"] == "meth_hypo_pct")]["value"].iloc[0]

        ## 6. total molecules
        total_mol = qc_report[(qc_report["subpanel"] == subpanel) & \
                              (qc_report["mbdqc_metric"] == "total_mole")]["value"].iloc[0]
        total_mol_ex_mix = qc_report[(qc_report["subpanel"] == subpanel) & \
                                     (qc_report["mbdqc_metric"] == "total_mole_exclude_mixed")]["value"].iloc[0]

        sub_mbd_qc_report = pd.DataFrame({"fc": row["fc"],
                                          "run_sample_id": row["run_sample_id"],
                                          "path": row["path"],
                                          "cg_zero_hyper": cg_zero_hyper,
                                          "unmeth_hyper": unmeth_hyper,
                                          "methy_1_2": methy_1_2,
                                          "methy_hypo": methy_hypo,
                                          "hyper_frac": hyper_frac,
                                          "hypo_frac": hypo_frac,
                                          "residual_frac": residual_frac,
                                          "mixed_frac": mixed_frac,
                                          "total_mol": total_mol,
                                          "total_mol_ex_mix": total_mol_ex_mix}, index=[0])

        mbd_qc_report_df = mbd_qc_report_df.append(sub_mbd_qc_report)
        mbd_qc_report_df = mbd_qc_report_df.reset_index().drop("index", axis=1)

        if write==True:
            for fc in sample_df["fc"].unique():
                folder = os.path.join(OUTPUT_PATH, sample_df[sample_df["fc"] == fc]["path"].unique()[0].split("/")[-1])
                if not os.path.exists(folder):
                    os.makedirs(folder)
                mbd_qc_report_df.to_csv(os.path.join(folder, "mbd_qc_default_report.hdr.tsv"),
                                              sep="\t",index=False)
    return mbd_qc_report_df

## 4. Given fcs, sample_ids, generate df for MBD binding curves, and length distribution for control probes
def get_binding_curve_length_dist(fcs, samples="", write=True):
    ## get the df containing paths to fcs and samples
    sample_df = process_fc_sample(fcs, samples)

    binding_curve_df = pd.DataFrame()
    insert_len_dist_df = pd.DataFrame()

    for index, row in sample_df.iterrows():

        binding = pd.read_csv(os.path.join(row["path"], row["run_sample_id"] + ".mbdqc_binding_curve.hdr.tsv"),
                              sep="\t")

        # calculate binding curves
        # filter by insert size
        binding_curve = binding[binding["insert_size"].between(120, 240)]
        total_count_at_cg_num = binding_curve.drop("insert_size", axis=1).groupby(
            ["control_type", "num_cg"]).sum().reset_index()
        total_count_at_cg_num.columns = ["control_type", "num_cg", "total_molecule_count"]
        partition_count = binding_curve.drop("insert_size", axis=1).groupby(
            ["control_type", "partition", "num_cg"]).sum().reset_index()

        # merge
        sub_binding_curve = pd.merge(partition_count, total_count_at_cg_num, on=["control_type", "num_cg"], how="outer")
        sub_binding_curve["mole_frac"] = (sub_binding_curve["molecule_count"] + 1) / \
                                         (sub_binding_curve["total_molecule_count"] + 3)
        sub_binding_curve["fc"] = row["fc"]
        sub_binding_curve["run_sample_id"] = row["run_sample_id"]
        sub_binding_curve = sub_binding_curve.sort_values(by="num_cg")
        binding_curve_df = binding_curve_df.append(sub_binding_curve, ignore_index=True)

        # calculate length distribution for control probes
        # drop unused cols
        length = binding.drop(["num_cg", "partition"], axis=1)

        # bin insert size into 10bp bins, up to 1000bp
        bins = np.arange(0, 1000, 10)
        length["binned"] = pd.cut(length["insert_size"],bins)
        length["insert_size"] = length["binned"].apply(lambda x: float(x.left))

        length_control = length.groupby(["control_type", "insert_size"]).sum().reset_index()
        control_total_mol = length.drop("insert_size", axis=1).groupby(["control_type"]).sum().reset_index()
        control_total_mol.columns = ["control_type", "total_mol"]

        sub_insert_len_dist_df = pd.merge(length_control, control_total_mol, on="control_type", how="outer")
        sub_insert_len_dist_df["insert_len_frac"] = sub_insert_len_dist_df["molecule_count"] / sub_insert_len_dist_df[
            "total_mol"]

        sub_insert_len_dist_df["fc"] = row["fc"]
        sub_insert_len_dist_df["run_sample_id"] = row["run_sample_id"]
        sub_insert_len_dist_df = sub_insert_len_dist_df.sort_values(by="insert_size")
        insert_len_dist_df = insert_len_dist_df.append(sub_insert_len_dist_df)
    insert_len_dist_df = insert_len_dist_df[insert_len_dist_df["molecule_count"] > 0]
    if write == True:
        for fc in sample_df["fc"].unique():
            folder = os.path.join(OUTPUT_PATH,sample_df[sample_df["fc"] == fc]["path"].unique()[0].split("/")[-1])
            if not os.path.exists(folder):
                os.makedirs(folder)
            binding_curve_df.to_csv(os.path.join(folder, "mbd_binding_curve.hdr.tsv"),
                                          sep="\t", index=False)
            insert_len_dist_df.to_csv(os.path.join(folder, "mbd_length_dist_control.hdr.tsv"),
                                          sep="\t", index=False)


    return binding_curve_df,insert_len_dist_df

## 5. Given binding curve, generate signal to noise ratio
def get_value(x):
    if len(x) == 0:
        return 0
    else:
        return int(x)

def ratio(binding_curve, partition, cg):
    control_hyper = (get_value(binding_curve[(binding_curve["control_type"] == "hyper") & \
                                             (binding_curve["partition"] == partition) & \
                                             (binding_curve["num_cg"] == cg)]["molecule_count"]) + 1) / \
                    (get_value(binding_curve[(binding_curve["control_type"] == "hyper") & \
                                             (binding_curve["num_cg"] == cg)]["total_molecule_count"].unique()) + 3)

    control_hypo = (get_value(binding_curve[(binding_curve["control_type"] == "hypo") & \
                                            (binding_curve["partition"] == partition) & \
                                            (binding_curve["num_cg"] == cg)]["molecule_count"]) + 1) / \
                   (get_value(binding_curve[(binding_curve["control_type"] == "hypo") & \
                                            (binding_curve["num_cg"] == cg)]["total_molecule_count"].unique()) + 3)

    ratio = np.log(control_hyper / control_hypo)
    return (ratio)


def get_sig_noise_radio_df(binding_curve_df,write=True):
    samples = binding_curve_df["run_sample_id"].unique()

    sig_noise_ratio_df = pd.DataFrame()

    for sample in samples:
        sample_binding_curve = binding_curve_df[binding_curve_df["run_sample_id"] == sample]

        for cg in sample_binding_curve["num_cg"].unique():
            hyper_sig_ratio = ratio(sample_binding_curve, "hyper", cg)
            residual_sig_ratio = ratio(sample_binding_curve, "residual", cg)
            hypo_sig_ratio = ratio(sample_binding_curve, "hypo", cg)

            sub_sig_noise_ratio = pd.DataFrame({"fc": sample_binding_curve["fc"].iloc[0],
                                                "run_sample_id": sample,
                                                "num_cg": cg,
                                                "sig_to_noise_ratio": max(hyper_sig_ratio, residual_sig_ratio),
                                                "hypo_sig_ratio": hypo_sig_ratio},
                                               index=[0])
            sig_noise_ratio_df = sig_noise_ratio_df.append(sub_sig_noise_ratio)
    sig_noise_ratio_df = sig_noise_ratio_df.sort_values(by="num_cg")

    if write==True:
        for fc in binding_curve_df["fc"].unique():
            fc_folder = get_fc_path(fc)
            folder = os.path.join(OUTPUT_PATH, fc_folder.split("/")[-1])
            if not os.path.exists(folder):
                os.makedirs(folder)
            sig_noise_ratio_df.to_csv(os.path.join(folder, "mbd_qc_sig_noise_ratio.hdr.tsv"),
                                      sep="\t",index=False)
    return sig_noise_ratio_df

## 6. Given fcs, samples_ids, generate diversities
def get_diversity(fcs, samples="", subpanel="backbone",write=True):
    ## get the df containing paths to fcs and samples
    sample_df = process_fc_sample(fcs, samples)

    diversity_df = pd.DataFrame()

    for index, row in sample_df.iterrows():
        sub_diversity = pd.read_csv(os.path.join(row["path"], row["run_sample_id"] + ".subpanel_qc_stats.hdr.tsv"),
                                    sep="\t")

        sub_diversity_df = pd.DataFrame({"fc": row["fc"],
                                         "run_sample_id": row["run_sample_id"],
                                         "diversity": sub_diversity[(sub_diversity["subpanel"] == subpanel) & \
                                                                    (sub_diversity["metric"] == "diversity")]["value"]})
        diversity_df = diversity_df.append(sub_diversity_df, ignore_index=True)

    if write == True:
        for fc in sample_df["fc"].unique():
            folder = os.path.join(OUTPUT_PATH,sample_df[sample_df["fc"] == fc]["path"].unique()[0].split("/")[-1])
            if not os.path.exists(folder):
                os.makedirs(folder)
            diversity_df.to_csv(os.path.join(folder, "mbd_diversity.hdr.tsv"),
                                sep="\t", index=False)
    return (diversity_df)

## 7. Given fcs, samples_ids, generate length distribution for all the mols in each partiton
def get_frag_length(fcs, samples="",write=True):
    sample_df = process_fc_sample(fcs, samples)

    frag_len_df = pd.DataFrame()

    # bin insert size into 10bp bins, up to 1000bp
    bins = np.arange(0, 1000, 10)

    for index, row in sample_df.iterrows():
        ## 1. hyper partition
        hyper = pd.read_csv(os.path.join(row["path"], row["run_sample_id"] + ".mbd_hyper_molecules.tsv"),
                            sep="\t",
                            low_memory=False)



        # calc insert size
        hyper["insert_size"] = hyper["end"] - hyper["start"]
        hyper["binned"] = pd.cut(hyper["insert_size"], bins)
        hyper_len = hyper.drop(["start", "end", "num_reads", "mean_mapq", "cg_count","insert_size"], axis=1).groupby(
            "binned").sum().reset_index()
        hyper_len["frac"] = hyper_len["num_molecules"].apply(lambda x: x / hyper_len["num_molecules"].sum())
        hyper_len["partition"] = "hyper"
        hyper_len["fc"] = row["fc"]
        hyper_len["run_sample_id"] = row["run_sample_id"]

        ## 2. hypo partition
        hypo = pd.read_csv(os.path.join(row["path"], row["run_sample_id"] + ".mbd_hypo_molecules.tsv"),
                           sep="\t",
                           low_memory=False)
        # calc insert size
        hypo["insert_size"] = hypo["end"] - hypo["start"]
        hypo["binned"] = pd.cut(hypo["insert_size"], bins)
        hypo_len = hypo.drop(["start", "end", "num_reads", "mean_mapq", "cg_count","insert_size"], axis=1).groupby(
            "binned").sum().reset_index()
        hypo_len["frac"] = hypo_len["num_molecules"].apply(lambda x: x / hypo_len["num_molecules"].sum())
        hypo_len["partition"] = "hypo"
        hypo_len["fc"] = row["fc"]
        hypo_len["run_sample_id"] = row["run_sample_id"]

        ## 3. residual partition
        residual = pd.read_csv(os.path.join(row["path"], row["run_sample_id"] + ".mbd_residual_molecules.tsv"),
                               sep="\t",
                               low_memory=False)
        # calc insert size
        residual["insert_size"] = residual["end"] - residual["start"]
        residual["binned"] = pd.cut(residual["insert_size"], bins)
        residual_len = residual.drop(["start", "end", "num_reads", "mean_mapq", "cg_count","insert_size"], axis=1).groupby(
            "binned").sum().reset_index()
        residual_len["frac"] = residual_len["num_molecules"].apply(lambda x: x / residual_len["num_molecules"].sum())
        residual_len["partition"] = "residual"
        residual_len["fc"] = row["fc"]
        residual_len["run_sample_id"] = row["run_sample_id"]


        ## concat all partitions
        sub_frag_len = pd.concat([hyper_len, hypo_len, residual_len],
                                 ignore_index=True)
        sub_frag_len["insert_size"] = sub_frag_len["binned"].apply(lambda x: float(x.left))
        sub_frag_len = sub_frag_len.sort_values(by="insert_size")
        frag_len_df = frag_len_df.append(sub_frag_len)

    frag_len_df = frag_len_df[frag_len_df["num_molecules"] > 0]
    if write == True:
        for fc in sample_df["fc"].unique():
            folder = os.path.join(OUTPUT_PATH,sample_df[sample_df["fc"] == fc]["path"].unique()[0].split("/")[-1])
            if not os.path.exists(folder):
                os.makedirs(folder)
            frag_len_df.to_csv(os.path.join(folder, "mbd_fragment_length.hdr.tsv"),
                               sep="\t", index=False)
    return (frag_len_df)

## 8. Given fc and samples, generate lambda spike-in result
def get_lambda_spike(fcs, samples="",write=True):
    ## get the df containing paths to fcs and samples
    sample_df = process_fc_sample(fcs, samples)

    lambda_spike_df = pd.DataFrame()

    for index, row in sample_df.iterrows():
        lambda_spike = pd.read_csv(
            os.path.join(row["path"], row["run_sample_id"] + ".mbdqc_lambda_spike_summary.hdr.tsv"),
            sep="\t")

        sub_lambda = pd.DataFrame({"fc": row["fc"],
                                   "run_sample_id": row["run_sample_id"],
                                   "partition": lambda_spike["partition"],
                                   "spike_id": lambda_spike["spike_id"].apply(lambda x: x.split("_")[-1]),
                                   "molecule_count": lambda_spike["molecule_count"]})
        total_lambda_count = lambda_spike.drop("read_count", axis=1).groupby("spike_id").sum().reset_index()
        total_lambda_count["spike_id"] = total_lambda_count["spike_id"].apply(lambda x: x.split("_")[-1])
        total_lambda_count.columns = ["spike_id", "total_molecule_count"]

        sub_lambda = pd.merge(sub_lambda, total_lambda_count, on="spike_id", how="outer")
        sub_lambda["percentage"] = sub_lambda["molecule_count"] / sub_lambda["total_molecule_count"]

        lambda_spike_df = lambda_spike_df.append(sub_lambda, ignore_index=True)
        lambda_spike_df = lambda_spike_df.sort_values(by="spike_id")
        lambda_spike_df["percentage"] = lambda_spike_df["percentage"].apply(lambda x: None if x == 0 else x)
    if write == True:
        for fc in sample_df["fc"].unique():
            folder = os.path.join(OUTPUT_PATH,sample_df[sample_df["fc"] == fc]["path"].unique()[0].split("/")[-1])
            if not os.path.exists(folder):
                os.makedirs(folder)
            lambda_spike_df.to_csv(os.path.join(folder, "mbd_lambda_spike_in.hdr.tsv"),
                               sep="\t", index=False)
    return lambda_spike_df

## 9. Given fc and samples, generate molecules distribution given number of cgs in each partition
def get_mol_dist_at_cg(fcs, samples="",write=True):
    sample_df = process_fc_sample(fcs, samples)

    mol_dist_at_cg_df = pd.DataFrame()

    for index, row in sample_df.iterrows():
        ## 1. hyper partition
        hyper = pd.read_csv(os.path.join(row["path"], row["run_sample_id"] + ".mbd_hyper_molecules.tsv"),
                            sep="\t",
                            low_memory=False)
        # we only consider cg_count up to 40
        hyper = hyper[hyper["cg_count"] <= 40]
        # calc num of mols at each cg
        hyper_dist = hyper.drop(["chrom", "start", "end", "num_reads", "mean_mapq"], axis=1).groupby(
            "cg_count").sum().reset_index()
        hyper_dist["partition"] = "hyper"
        hyper_dist["fc"] = row["fc"]
        hyper_dist["run_sample_id"] = row["run_sample_id"]

        ## 2. hypo partition
        hypo = pd.read_csv(os.path.join(row["path"], row["run_sample_id"] + ".mbd_hypo_molecules.tsv"),
                           sep="\t",
                           low_memory=False)
        # we only consider cg_count up to 40
        hypo = hypo[hypo["cg_count"] <= 40]
        # calc num of mols at each cg
        hypo_dist = hypo.drop(["chrom", "start", "end", "num_reads", "mean_mapq"], axis=1).groupby(
            "cg_count").sum().reset_index()
        hypo_dist["partition"] = "hypo"
        hypo_dist["fc"] = row["fc"]
        hypo_dist["run_sample_id"] = row["run_sample_id"]

        ## 3. residual partition
        residual = pd.read_csv(os.path.join(row["path"], row["run_sample_id"] + ".mbd_residual_molecules.tsv"),
                               sep="\t",
                               low_memory=False)
        # we only consider cg_count up to 40
        residual = residual[residual["cg_count"] <= 40]
        # calc insert size
        residual_dist = residual.drop(["chrom", "start", "end", "num_reads", "mean_mapq"], axis=1).groupby(
            "cg_count").sum().reset_index()
        residual_dist["partition"] = "residual"
        residual_dist["fc"] = row["fc"]
        residual_dist["run_sample_id"] = row["run_sample_id"]

        ## concat all partitions
        sub_mols_dist = pd.concat([hyper_dist, hypo_dist, residual_dist],
                                  ignore_index=True)
        sub_mols_dist["fraction"] = sub_mols_dist["num_molecules"].apply(
            lambda x: x / sub_mols_dist["num_molecules"].sum())

        mol_dist_at_cg_df = mol_dist_at_cg_df.append(sub_mols_dist)

    mol_dist_at_cg_df = mol_dist_at_cg_df.sort_values(by="cg_count")
    if write == True:
        for fc in sample_df["fc"].unique():
            folder = os.path.join(OUTPUT_PATH,sample_df[sample_df["fc"] == fc]["path"].unique()[0].split("/")[-1])
            if not os.path.exists(folder):
                os.makedirs(folder)
            mol_dist_at_cg_df.to_csv(os.path.join(folder, "mbd_molecule_dist_at_cg.tsv"),
                                     sep="\t", index=False)

    return (mol_dist_at_cg_df)


(binding_curve, control_length_dist) = get_binding_curve_length_dist(fcs="180921_A00441_0096_AHG37KDMXX")
get_mbd_qc_report(fcs="180921_A00441_0096_AHG37KDMXX")
get_sig_noise_radio_df(binding_curve)
get_diversity(fcs="180921_A00441_0096_AHG37KDMXX")
get_frag_length(fcs="180921_A00441_0096_AHG37KDMXX")
get_lambda_spike(fcs="180921_A00441_0096_AHG37KDMXX")
#get_mol_dist_at_cg(fcs="180921_A00441_0096_AHG37KDMXX")