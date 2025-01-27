import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn import linear_model
from scipy import stats
import tangram as tg
import timeit
from datetime import timedelta
import multiprocessing as mp
import anndata as ad

def validate_extrapolation(adata_sc, adata_st, r2_threshold = 0.8, spearman_threshold = 0.4, threads = 1, **kwargs):
    print("Validating Extrapolation...")
    overlap_genes = list(set(adata_sc.var_names.values).intersection(set(adata_st.var_names.values)))
    adata_st.uns["r2_threshold"] = r2_threshold
    adata_st.uns["spearman_threshold"] = spearman_threshold
    scores = []
    scores_cross = []
    start = timeit.default_timer()
    share_sc_X = mp.RawArray('d', len(adata_sc[:, overlap_genes].X.toarray().flatten()))
    share_st_X = mp.RawArray('d', len(adata_st[:, overlap_genes].X.toarray().flatten()))
    scX = adata_sc[:, overlap_genes].X.toarray()
    stX = adata_st[:, overlap_genes].X.toarray()
    share_sc_X_arr = np.frombuffer(share_sc_X).reshape(scX.shape)
    share_st_X_arr = np.frombuffer(share_st_X).reshape(stX.shape)
    np.copyto(share_sc_X_arr, scX)
    np.copyto(share_st_X_arr, stX)
    scores = []
    scores_cross = []
    alphas = []
    for i in range(10):
        s = []
        sc = []
        for gene in overlap_genes:
            X = adata_sc[:, list(set(overlap_genes) - set([gene]))].X.toarray()
            y = adata_sc[:, [gene]].X.toarray()
            X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2)
            clf = linear_model.LassoCV(cv=3, n_jobs=-1)
            clf.fit(X_train, y_train)
            # clf = Lasso(alpha = clf.alpha_, **kwargs)
            # # clf = linear_model.ElasticNet(random_state=0)
            # clf.fit(X_train, y_train)
            s.append(clf.score(X_test, y_test))
            alphas.append(clf.alpha_)
            clf.fit(X, y)
            y_pred = clf.predict(adata_st[:, list(set(overlap_genes) - set([gene]))].X.toarray())
            sc.append(stats.spearmanr(y_pred, adata_st[:, [gene]].X.toarray()).statistic)
            scores.append(s)
            scores_cross.append(sc)
    scores = np.array(scores)
    scores_cross = np.array(scores_cross)
    adata_st.uns["extrapolate_validation_alphas"] = alphas
    adata_st.uns["extrapolate_validation"] = pd.DataFrame(zip(*[overlap_genes, np.mean(scores, axis=0), np.mean(scores_cross, axis=0)]), columns=["genes", "scRNAseq", "ST"]) #scRNAseq is R2 score and ST is spearman correlation
    print("Time to validate extrapolation", timedelta(seconds=timeit.default_timer() - start))
    print(f'Number of genes above scRNA-seq threshold {adata_st.uns["extrapolate_validation"][adata_st.uns["extrapolate_validation"]["scRNAseq"] >= r2_threshold].shape[0]}')
    print(f'Number of genes above scRNA-seq threshold and R2 threshold {adata_st.uns["extrapolate_validation"][(adata_st.uns["extrapolate_validation"]["ST"] >= spearman_threshold) & (adata_st.uns["extrapolate_validation"]["scRNAseq"] >= r2_threshold)].shape[0]}')
    return adata_st

def extrapolate(adata_sc, adata_st, alpha=0.1, **kwargs):
    print("Extrapolating...")
    # adata_st.uns["r2_threshold"] = 0.8
    # adata_st.uns["spearman_threshold"] = 0.4
    overlap_genes = list(set(adata_sc.var_names.values).intersection(set(adata_st.var_names.values)))
    relevant_genes = list(set(adata_sc.var_names.values).difference(set(adata_st.var_names.values)))
    print("Number of Genes to extrapolate for :", len(relevant_genes))
    scores = []
    df_sc_extrapolated = adata_st.to_df().copy()
    start = timeit.default_timer()
    for gene in relevant_genes:
        X = adata_sc[:, overlap_genes].X.toarray()
        y = adata_sc[:, [gene]].X.toarray()
        X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=0)
        clf = linear_model.LassoCV(cv=3, n_jobs=-1)
        clf.fit(X_train, y_train)
        # clf = Lasso(alpha = clf.alpha_, **kwargs)
        # clf.fit(X_train, y_train)
        scores.append(clf.score(X_test, y_test))
        if scores[-1] >= adata_st.uns["r2_threshold"] and scores[-1] < 1:
            clf.fit(X, y)
            y_pred = clf.predict(adata_st[:, overlap_genes].X.toarray())
            df_sc_extrapolated.loc[:, f'{gene}'] = y_pred
    print("Time to extrapolate", timedelta(seconds=timeit.default_timer() - start))
    score_df = pd.DataFrame(zip(*[relevant_genes, scores]), columns=["genes", "scRNAseq"])
    print(f'Number of genes above scRNA-seq threshold {score_df[(score_df["R2 score"] < 1) & (score_df["R2 score"] > adata_st.uns["r2_threshold"])].shape[0]}')
    validation_ratio = adata_st.uns["extrapolate_validation"][(adata_st.uns["extrapolate_validation"]["ST"] >= adata_st.uns["spearman_threshold"]) & (adata_st.uns["extrapolate_validation"]["scRNAseq"] >= adata_st.uns["r2_threshold"])].shape[0] / adata_st.uns["extrapolate_validation"][adata_st.uns["extrapolate_validation"]["scRNAseq"] >= adata_st.uns["r2_threshold"]].shape[0]
    print(f'Based on validation thresholds, we can expect {int(validation_ratio * score_df[(score_df["R2 score"] < 1) & (score_df["R2 score"] > adata_st.uns["r2_threshold"])].shape[0])} genes to have spearman correlation > {adata_st.uns["spearman_threshold"]}')
    adata_st.uns["extrapolation_results"] = score_df
    adata_st.uns["X_extrapolated"] = df_sc_extrapolated
    return adata_st


def _parallelize_training(args):
    gene, n_genes, g_position_sc, g_position_st, overlap_genes, r2 = args[0], args[1], args[2], args[3], args[4], args[5]
    sc_X_local = np.frombuffer(sc_X).reshape(sc_X_shape)
    other_genes_sc = list(range(n_genes))
    other_genes_sc.pop(g_position_sc)
    X = sc_X_local[:, other_genes_sc]
    y = sc_X_local[:, g_position_sc]
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2)
    clf = linear_model.LassoCV(cv=3, n_jobs=-1)
    clf.fit(X_train, y_train)
    score = clf.score(X_test, y_test)
    clf.fit(X, y)
    if score > r2 and score < 1:
        other_genes_st = list(range(overlap_genes))
        other_genes_st.pop(g_position_st)
        st_X_local = np.frombuffer(st_X).reshape(st_X_shape)
        y_pred = clf.predict(st_X_local[:, other_genes_st])
        return gene, score, y_pred
    else:
        return None

def _initializer_func(X1, X1_shape, X2, X2_shape):
    global sc_X,  sc_X_shape
    sc_X = X1
    sc_X_shape = X1_shape
    global st_X, st_X_shape
    st_X = X2
    st_X_shape = X2_shape

def extrapolate_pll(adata_sc, adata_st, alpha=0.1, **kwargs):
    print("Extrapolating...")
    # adata_st.uns["r2_threshold"] = 0.8
    # adata_st.uns["spearman_threshold"] = 0.4
    overlap_genes = list(set(adata_sc.var_names.values).intersection(set(adata_st.var_names.values)))
    relevant_genes = list(set(adata_sc.var_names.values).difference(set(adata_st.var_names.values)))
    print("Number of Genes to extrapolate for :", len(relevant_genes))
    scores = []
    df_sc_extrapolated = adata_st.to_df().copy()
    start = timeit.default_timer()
    to_parallelize = []
    share_sc_X = mp.RawArray('d', len(adata_sc.X.toarray().flatten()))
    share_st_X = mp.RawArray('d', len(adata_st[:, overlap_genes].X.toarray().flatten()))
    scX = adata_sc.X.toarray()
    stX = adata_st[:, overlap_genes].X.toarray()
    share_sc_X_arr = np.frombuffer(share_sc_X).reshape(scX.shape)
    share_st_X_arr = np.frombuffer(share_st_X).reshape(stX.shape)
    np.copyto(share_sc_X_arr, scX)
    np.copyto(share_st_X_arr, stX)
    for gene in relevant_genes:
            to_parallelize.append([gene, len(adata_sc.var_names.values), list(adata_sc.var_names.values).index(gene), list(adata_st[:, overlap_genes].var_names.values).index(gene), len(overlap_genes), adata_st.uns["r2_threshold"]])
    
    with mp.Pool(processes=128, initializer=_initializer_func, initargs=(share_sc_X, scX.shape, share_st_X, stX.shape), maxtasksperchild = 1) as pool:
        results = pool.map(_parallelize_training, to_parallelize)

    results_dict = {}
    for i in results:
        if i[0] not in results_dict:
            results_dict[i[0]] = [[i[1]], [i[2]]]
        else:
            results_dict[i[0]][0].append(i[1])
            results_dict[i[0]][1].append(i[2])
    scores = []
    scores_cross = []
    for gene, in results_dict.keys():
        scores.append(np.nanmean(results_dict[gene][0]))
        scores_cross.append(np.nanmean(results_dict[gene][1]))
        df_sc_extrapolated.loc[:, f'{gene}'] = y_pred
    print("Time to extrapolate", timedelta(seconds=timeit.default_timer() - start))
    score_df = pd.DataFrame(zip(*[relevant_genes, scores]), columns=["genes", "scRNAseq"])
    print(f'Number of genes above scRNA-seq threshold {score_df[(score_df["R2 score"] < 1) & (score_df["R2 score"] > adata_st.uns["r2_threshold"])].shape[0]}')
    validation_ratio = adata_st.uns["extrapolate_validation"][(adata_st.uns["extrapolate_validation"]["ST"] >= adata_st.uns["spearman_threshold"]) & (adata_st.uns["extrapolate_validation"]["scRNAseq"] >= adata_st.uns["r2_threshold"])].shape[0] / adata_st.uns["extrapolate_validation"][adata_st.uns["extrapolate_validation"]["scRNAseq"] >= adata_st.uns["r2_threshold"]].shape[0]
    print(f'Based on validation thresholds, we can expect {int(validation_ratio * score_df[(score_df["R2 score"] < 1) & (score_df["R2 score"] > adata_st.uns["r2_threshold"])].shape[0])} genes to have spearman correlation > {adata_st.uns["spearman_threshold"]}')
    adata_st.uns["extrapolation_results"] = score_df
    adata_st.uns["X_extrapolated"] = df_sc_extrapolated
    return adata_st

def run_tangram(adata_sc, adata_st, device="cpu"):
    '''
    Run Tangram for extrapolation of genes from scRNA-seq to ST.
    Arguments
    ----------
    adata_sc : AnnData
        Anndata object containing scRNA-seq data.
    adata_st : AnnData
        Anndata object containing spatial transcriptomics data.
    device : str
        Device to run Tangram on. Either 'cpu' or 'cuda'.
    '''
    print("Running Tangram...")
    start = timeit.default_timer()
    tg.pp_adatas(adata_sc, adata_st, genes=None)
    adata_map = tg.map_cells_to_space(adata_sc, adata_st, device=device)
    adata_ge = tg.project_genes(adata_map, adata_sc)
    for excluded_gene in list(set(adata_st.var_names.values).difference(set(adata_sc.var_names.values))):
        temp = adata_ge.to_df()
        temp.loc[:, excluded_gene] = adata_st[:, excluded_gene].X.toarray()
        adata = ad.AnnData(X=temp)
        adata.obs_names = temp.index
        adata.obs.index = temp.index
        adata.var_names = temp.columns
        adata_ge = adata.copy()
    adata_ge.raw = adata_st.copy()
    for i in adata_st.obs.keys():
        adata_ge.obs[i] = adata_st.obs[i].copy()
    for i in adata_st.obsm.keys():
        adata_ge.obsm[i] = adata_st.obsm[i].copy()
    for i in adata_st.uns.keys():
        adata_ge.uns[i] = adata_st.uns[i].copy()
    print("Time to run Tangram", timedelta(seconds=timeit.default_timer() - start))
    return adata_ge

# def validate_extrapolation(adata_sc, adata_st, r2_threshold = 0.8, spearman_threshold = 0.4, threads = 1, **kwargs):
#     print("Validating Extrapolation...")
#     overlap_genes = list(set(adata_sc.var_names.values).intersection(set(adata_st.var_names.values)))
#     print(overlap_genes, len(overlap_genes))
#     adata_st.uns["r2_threshold"] = r2_threshold
#     adata_st.uns["spearman_threshold"] = spearman_threshold
#     scores = []
#     scores_cross = []
#     start = timeit.default_timer()
#     share_sc_X = mp.RawArray('d', len(adata_sc[:, overlap_genes].X.toarray().flatten()))
#     share_st_X = mp.RawArray('d', len(adata_st[:, overlap_genes].X.toarray().flatten()))
#     scX = adata_sc[:, overlap_genes].X.toarray()
#     stX = adata_st[:, overlap_genes].X.toarray()
#     share_sc_X_arr = np.frombuffer(share_sc_X).reshape(scX.shape)
#     share_st_X_arr = np.frombuffer(share_st_X).reshape(stX.shape)
#     np.copyto(share_sc_X_arr, scX)
#     np.copyto(share_st_X_arr, stX)
#     to_parallelize = []
#     for i in range(10):
#         for gene in overlap_genes:
#             to_parallelize.append([gene, len(overlap_genes), list(adata_sc[:, overlap_genes].var_names.values).index(gene), list(adata_st[:, overlap_genes].var_names.values).index(gene)])
#     with mp.Pool(processes=threads, initializer=_initializer_func, initargs=(share_sc_X, scX.shape, share_st_X, stX.shape), maxtasksperchild = 1) as pool:
#         results = pool.map(_parallelize_training, to_parallelize)
#     del scX, stX
    
#     results_dict = {}
#     for i in results:
#         if i[0] not in results_dict:
#             results_dict[i[0]] = [[i[1]], [i[2]]]
#         else:
#             results_dict[i[0]][0].append(i[1])
#             results_dict[i[0]][1].append(i[2])
#     scores = []
#     scores_cross = []
#     for gene, in results_dict.keys():
#         scores.append(np.nanmean(results_dict[gene][0]))
#         scores_cross.append(np.nanmean(results_dict[gene][1]))

#         # X = adata_sc[:, list(set(overlap_genes) - set([gene]))].X.toarray()
#         # y = adata_sc[:, [gene]].X.toarray()
#         # X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2)
#         # clf = linear_model.LassoCV(cv=3, n_jobs=-1)
#         # clf.fit(X_train, y_train)
#         # # clf = Lasso(alpha = clf.alpha_, **kwargs)
#         # # # clf = linear_model.ElasticNet(random_state=0)
#         # # clf.fit(X_train, y_train)
#         # s.append(clf.score(X_test, y_test))
#         # clf.fit(X, y)
#         # y_pred = clf.predict(adata_st[:, list(set(overlap_genes) - set([gene]))].X.toarray())
#         # sc.append(stats.spearmanr(y_pred, adata_st[:, [gene]].X.toarray()).statistic)
#         # scores.append(s)
#         # scores_cross.append(sc)
#     scores = np.array(scores)
#     scores_cross = np.array(scores_cross)
#     adata_st.uns["extrapolate_validation"] = pd.DataFrame(zip(*[overlap_genes, scores, scores_cross]), columns=["genes", "scRNAseq", "ST"]) #scRNAseq is R2 score and ST is spearman correlation
#     print("Time to validate extrapolation", timedelta(seconds=timeit.default_timer() - start))
#     print(f'Number of genes above scRNA-seq threshold {adata_st.uns["extrapolate_validation"][adata_st.uns["extrapolate_validation"]["scRNAseq"] >= r2_threshold].shape[0]}')
#     print(f'Number of genes above scRNA-seq threshold and R2 threshold {adata_st.uns["extrapolate_validation"][(adata_st.uns["extrapolate_validation"]["ST"] >= spearman_threshold) & (adata_st.uns["extrapolate_validation"]["scRNAseq"] >= r2_threshold)].shape[0]}')
#     return adata_st