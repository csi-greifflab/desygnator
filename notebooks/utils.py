import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import os
import sys
import genmodel

hex_colors = ["#ffffff", "#3dcccc", "#80aaff", "#3939e6", "#a179f2", "#cc2936", "#ff9580", "#B75194", "#3423A6", "#ABD2FA"]

colors = [mpl.colors.hex2color(c) for c in hex_colors]
coef = 0.7
colors = colors + [tuple(np.array(c) * 0.7) for c in colors[1:]]


sample_type_to_color = {"synthetic": colors[1],
                        "data": colors[2],
                        "technical": colors[3],
                        "biological": colors[4],
                        "cross-tissue": colors[7],
                        "twins": colors[5],
                        "twin": colors[5],
                        "unrelated": colors[6]}
sample_type_to_label = {"synthetic": "synthetic replicates",
                        "data": "data replicates",
                        "technical": "technical replicates",
                        "biological": "biological replicates",
                        "cross-tissue": "cross-cell population replicates",
                        "twins": "twin subjects",
                        "twin": "twin subjects",
                        "unrelated": "unrelated subjects"}
sample_sizes = [1000, 3000, 10000, 30000]

def export_legend(legend, filename):
    fig  = legend.figure
    fig.canvas.draw()
    bbox  = legend.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    legend.get_frame().set_linewidth(0.0)
    fig.savefig(filename, dpi="figure", bbox_inches=bbox, format="png")
    
    
def transform_model(model):
    res = {}
    eps = 1e-9

    model_v = model.marginals[0]["v_choice"]
    model_j = model.marginals[0]["j_choice"]
    model_d = model.marginals[0]["d_gene"]
    repeated_v = np.tile(model_v, (model_j.shape[-1], model_d.shape[-1], 1))
    repeated_v = np.transpose(repeated_v, (2, 0, 1))
    repeated_j = np.tile(model_j, (model_d.shape[-1], 1, 1))
    repeated_j = np.transpose(repeated_j, (1, 2, 0))
    repeated_v.shape, repeated_j.shape, model_d.shape
    v_d_j = model_d * repeated_v * repeated_j + eps
    res["v_d_j"] = v_d_j

    model_delv = model.marginals[0]["v_3_del"]
    repeated_v = np.tile(model_v, (model_delv.shape[-1], 1))
    repeated_v = np.transpose(repeated_v, (1, 0))
    v_delv = model_delv * repeated_v + eps
    res["v_delv"] = v_delv

    model_delj = model.marginals[0]["j_5_del"]
    marginal_j = v_d_j.sum(axis=(0, 2))
    repeated_j = np.tile(marginal_j, (model_delj.shape[-1], 1))
    repeated_j = np.transpose(repeated_j, (1, 0))
    repeated_j.shape, model_delj.shape
    j_delj = model_delj * repeated_j + eps
    res["j_delj"] = j_delj

    model_deld5 = model.marginals[0]["d_5_del"]
    model_deld3 = model.marginals[0]["d_3_del"]
    marginal_d = v_d_j.sum(axis=(0, 1))
    repeated_d = np.tile(marginal_d, (model_deld5.shape[-1], model_deld3.shape[-1], 1))
    repeated_d = np.transpose(repeated_d, (2, 0, 1))
    repeated_deld5 = np.tile(model_deld5, (model_deld3.shape[-1], 1, 1))
    repeated_deld5 = np.transpose(repeated_deld5, (1, 2, 0))
    d_deld = model_deld3 * repeated_deld5 * repeated_d + eps
    res["d_deld"] = d_deld

    model_vd_ins = model.marginals[0]["vd_ins"] + eps
    model_vd_dinucl = model.marginals[0]["vd_dinucl"] + eps
    model_dj_ins = model.marginals[0]["dj_ins"] + eps
    model_dj_dinucl = model.marginals[0]["dj_dinucl"] + eps

    res["vd_ins"] = model_vd_ins
    res["dj_ins"] = model_dj_ins

    return res

def merge_vdj_alleles(model):
    v_names = model.events[0].get_realization_vector()
    d_names = model.events[1].get_realization_vector()
    j_names = model.events[2].get_realization_vector()
    v_gene_to_alleles = {}
    for i, allele in enumerate(v_names):
        gene_base_name, allele_id = allele.split("*")
        v_gene_to_alleles.setdefault(gene_base_name, [])
        v_gene_to_alleles[gene_base_name].append(i)
    j_gene_to_alleles = {}
    for i, allele in enumerate(j_names):
        gene_base_name, allele_id = allele.split("*")
        j_gene_to_alleles.setdefault(gene_base_name, [])
        j_gene_to_alleles[gene_base_name].append(i)
    d_gene_to_alleles = {}
    for i, allele in enumerate(d_names):
        gene_base_name, allele_id = allele.split("*")
        d_gene_to_alleles.setdefault(gene_base_name, [])
        d_gene_to_alleles[gene_base_name].append(i)
        
    
    transformed_model = transform_model(model)
    sorted_v_names = sorted(v_gene_to_alleles.keys())
    sorted_d_names = sorted(d_gene_to_alleles.keys())
    sorted_j_names = sorted(j_gene_to_alleles.keys())
    
    merged_transformed_model = {}
    
    merged_transformed_model["v_delv"] = np.array([transformed_model["v_delv"][v_gene_to_alleles[gene]].sum(axis=0) for gene in sorted_v_names])
    merged_transformed_model["v_d_j"] = np.array([transformed_model["v_d_j"][v_gene_to_alleles[gene]].sum(axis=0) for gene in sorted_v_names])
    
    merged_transformed_model["j_delj"] = np.array([transformed_model["j_delj"][j_gene_to_alleles[gene]].sum(axis=0) for gene in sorted_j_names])
    temp_vdj_matrix = merged_transformed_model["v_d_j"].transpose(1, 0, 2)
    merged_transformed_model["v_d_j"] = np.array([temp_vdj_matrix[j_gene_to_alleles[gene]].sum(axis=0) for gene in sorted_j_names]).transpose(1, 0, 2)
    merged_transformed_model["d_deld"] = np.array([transformed_model["d_deld"][d_gene_to_alleles[gene]].sum(axis=0) for gene in sorted_d_names])
    temp_vdj_matrix = merged_transformed_model["v_d_j"].transpose(2, 1, 0)
    merged_transformed_model["v_d_j"] = np.array([temp_vdj_matrix[d_gene_to_alleles[gene]].sum(axis=0) for gene in sorted_d_names]).transpose(2, 1, 0)
    
    for key, value in transformed_model.items():
        merged_transformed_model.setdefault(key, value)
    return merged_transformed_model


def full_js_distance(transformed_model_1, transformed_model_2, exclude=[]):
    res = 0.0
    for key in transformed_model_1.keys():
        if key in exclude:
            continue
        p = transformed_model_1[key]
        q = transformed_model_2[key]
        m = (p + q) / 2
        (p * np.log(p / q)).sum()
        res += (p * np.log(p / m)).sum() + (q * np.log(q / m)).sum()
    return (res / 2) ** 0.5

def full_shannon_entropy(transformed_model, exclude=[]):
    res = 0.0
    for key in transformed_model.keys():
        if key in exclude:
            continue
        p = transformed_model[key]
        res += -(p * np.log(p)).sum()
    return res


def draw_distances(statistics, sample_types, normalize, output_filename, ylim=(2**-3, 2**1)):
    plt.figure(figsize=(8, 5))
    fontsize = 30
    xs = [1000, 3000, 10000, 30000]
    for sample_type_id, sample_type in enumerate(sample_types):
        ys_main_median = []
        ys_main_min = []
        ys_main_max = []
        for x in xs:
            main_distances = statistics[x][sample_type]
            y_median = np.median(main_distances)
            y_min = np.min(main_distances)
            y_max = np.max(main_distances)
            if normalize:
                normalizer = np.mean(statistics[x]["synthetic"])
            else:
                normalizer = 1.0
            ys_main_median.append(y_median / normalizer) 
            ys_main_min.append(y_min / normalizer)
            ys_main_max.append(y_max / normalizer)

        markersize_coef = 0.5
        label = sample_type_to_label[sample_type]
        plt.plot(xs, ys_main_median, color=sample_type_to_color[sample_type], marker='o', markersize=fontsize * markersize_coef, label=f"{label}")
        plt.fill_between(xs, ys_main_min, ys_main_max, color=sample_type_to_color[sample_type], linestyle="--", alpha=0.1)

    plt.xlabel("# sequencing reads", fontsize=fontsize)
    plt.ylabel(f"{'Normalized' if normalize else 'Explicit'} JSD", fontsize=fontsize)

    plt.xscale("log")

    ax = plt.gca()
    ax.set_xticks(xs)
    ax.get_xaxis().set_major_formatter(mpl.ticker.ScalarFormatter())
    
    

    plt.ylim(ylim)

    plt.setp(plt.gca().get_xticklabels(), fontsize=fontsize * 0.8)
    plt.setp(plt.gca().get_yticklabels(), fontsize=fontsize * 0.8)
    plt.minorticks_off()


    plt.tight_layout()
    
    plt.savefig(output_filename, format="png")
    plt.show()

def get_distance_stats_test(TEST_DATA_DIR, n_sequences, names, merge=True, exclude=[]):
    run = 1
    models_a1 = []
    models_a2 = []
    other_models = {name: [[], []] for name in names}
    n_replicates = 30
    for replicate_id in range(n_replicates):
        prefix = os.path.join(TEST_DATA_DIR, f"A1_{replicate_id}_{n_sequences}_run/main_inference/final")
        model_1 = genmodel.GenModel(model_parms_file=prefix + "_parms.txt",
                                    marginals_file=prefix + "_marginals.txt")
        models_a1.append(model_1)
        prefix = os.path.join(TEST_DATA_DIR, f"A2_{replicate_id}_{n_sequences}_run/main_inference/final")
        model_2 = genmodel.GenModel(model_parms_file=prefix + "_parms.txt",
                                    marginals_file=prefix + "_marginals.txt")
        models_a2.append(model_2)
    for name in names:
        for replicate_id in range(n_replicates):
            prefix = os.path.join(TEST_DATA_DIR, f"{name}_A1_{replicate_id}_{n_sequences}_run/main_inference/final")
            model_1 = genmodel.GenModel(model_parms_file=prefix + "_parms.txt",
                                        marginals_file=prefix + "_marginals.txt")
            other_models[name][0].append(model_1)
            prefix = os.path.join(TEST_DATA_DIR, f"{name}_B1_{replicate_id}_{n_sequences}_run/main_inference/final")
            model_2 = genmodel.GenModel(model_parms_file=prefix + "_parms.txt",
                                        marginals_file=prefix + "_marginals.txt")
            other_models[name][1].append(model_2)
    
    transformed_models_a1 = [merge_vdj_alleles(m) if merge else transform_model(m) for m in models_a1]
    transformed_models_a2 = [merge_vdj_alleles(m) if merge else transform_model(m) for m in models_a2]
    transformed_other_models = {name: [[merge_vdj_alleles(m) if merge else transform_model(m) for m in p[0]],
                                       [merge_vdj_alleles(m) if merge else transform_model(m) for m in p[1]]] \
                                for name, p in other_models.items()}

    main_distances = [full_js_distance(transformed_models_a1[i],
                                       transformed_models_a2[i],
                                       exclude=exclude)\
                       for i in range(n_replicates)]
    other_distances = {name : [full_js_distance(p[0][i],
                                                p[1][i],
                                                exclude=exclude) \
                               for i in range(n_replicates)]\
                        for name, p in transformed_other_models.items()}
    return main_distances, other_distances

def draw_pvalues(pvalues, names, output_filename):
    lim_low = 1e-36
    lim_high = 1e2
    plt.figure(figsize=(8, 5))
    fontsize = 30
    xs = sample_sizes
    for size_id, sample_type in enumerate(names):
        pvalues_main = []
        for x in xs:
            p = pvalues[x][sample_type]
            pvalues_main.append(p)
        pvalues_main = np.clip(pvalues_main, a_min=lim_low, a_max=lim_high)
        if sample_type == "twin":
            pvalues_main *= 10 # to be observable in the picture
        markersize_coef = 0.5
        label = sample_type_to_label[sample_type]
        plt.plot(xs, pvalues_main, color=sample_type_to_color[sample_type], marker='o', markersize=fontsize * markersize_coef, label=f"{label}")
        

    significance_threshold = 0.01
    plt.axhline(significance_threshold, color="grey", linestyle="--")


    plt.xlabel("# sequencing reads", fontsize=fontsize)
    plt.ylabel("Adjusted p-value", fontsize=fontsize)

    plt.yscale("log")
    plt.xscale("log")

    ax = plt.gca()
    ax.set_xticks(sample_sizes)
    ax.get_xaxis().set_major_formatter(mpl.ticker.ScalarFormatter())

    tick_locs, tick_labels = plt.yticks()

    tick_locs = [10**i for i in range(-36, 4, 5)]
    tick_labels = ['$\\mathdefault{10^{' + str(i) + '}}$' for i in range(-36, 4, 5)]
    tick_labels[0] = '$\\mathdefault{<10^{' + str(-36) + '}}$'
    plt.yticks(tick_locs, tick_labels)

    plt.setp(plt.gca().get_xticklabels(), fontsize=fontsize * 0.8)
    plt.setp(plt.gca().get_yticklabels(), fontsize=fontsize * 0.7)
    plt.minorticks_off()

    plt.ylim((lim_low / 30, lim_high * 100))

    plt.tight_layout()
    plt.savefig(output_filename, format="png")
    plt.show()