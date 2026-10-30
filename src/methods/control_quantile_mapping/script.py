import sys

import anndata as ad
import numpy as np

## VIASH START
par = {
    "input": "resources_test/task_cyto_batch_integration/mouse_spleen_flow_cytometry_subset/censored_split1.h5ad",
    "output": "resources_test/task_cyto_batch_integration/mouse_spleen_flow_cytometry_subset/output.h5ad",
    "n_quantiles": 201,
    "n_clusters": 0,
    "shrinkage": 200,
    "soft": 0,
    "interp": "linear",
    "target": "mid",
    "allcells_weight": 0.0,
    "seed": 0,
}
meta = {"name": "control_quantile_mapping", "resources_dir": "src/utils"}
## VIASH END

# Control-anchored quantile mapping.
# The control samples (is_control > 0) are aliquots of the same donor measured in every batch, so for each marker the
# difference between their distributions across batches is purely technical. For every batch b and marker j we learn
# the monotone map F_ref^-1(F_b(x)) that sends the control cells of batch b onto the reference distribution of the same
# control donor (reference = the donor's cells pooled over batches with equal weight per batch, i.e. a "mid" target),
# average the maps over control donors and apply it to all cells of batch b.
# Optional (n_clusters > 0): the same maps are learned per cell population (k-means on lineage markers after the
# global correction) and shrunk towards the global map with weight n / (n + shrinkage).

print("Reading input", flush=True)
adata = ad.read_h5ad(par["input"])
X = np.asarray(adata.layers["preprocessed"], dtype=np.float64)
batch = adata.obs["batch"].astype(str).to_numpy()
ctrl = adata.obs["is_control"].astype(int).to_numpy()
to_corr = adata.var["to_correct"].to_numpy().astype(bool)
lineage = (adata.var["marker_type"].astype(str).to_numpy() == "lineage") & to_corr
batches = np.unique(batch)
donors = [d for d in np.unique(ctrl) if d > 0]
cols = np.where(to_corr)[0]
qs = np.linspace(0, 1, par["n_quantiles"])[1:-1]
print(f"{adata.n_obs} cells, batches {list(batches)}, control donors {donors}, {len(cols)} markers", flush=True)


def qmaps(Xc, bc, dc, donor_list=None):
    """For each batch: (source quantiles, target quantiles) [n_q, n_markers], averaged over control donors.
    Returns dict batch -> (src, tgt, n_cells)."""
    out = {}
    for b in batches:
        srcs, tgts, n = [], [], 0
        for d in (donors if donor_list is None else donor_list):
            ref_parts = [Xc[(dc == d) & (bc == bb)] for bb in batches]
            if any(len(p) < 20 for p in ref_parts):
                continue
            # reference quantiles: equal-weight mixture over batches = mean of per-batch quantile functions
            # (quantile averaging = Wasserstein barycenter in 1D)
            if par.get("target", "mid") == "goal":
                # reference = the control donor in the goal batch ("1"; first batch if absent)
                gb = "1" if "1" in list(batches) else batches[0]
                ref_q = np.quantile(ref_parts[list(batches).index(gb)], qs, axis=0)
            else:
                ref_q = np.mean([np.quantile(p, qs, axis=0) for p in ref_parts], axis=0)
            src_q = np.quantile(ref_parts[list(batches).index(b)], qs, axis=0)
            srcs.append(src_q)
            tgts.append(ref_q)
            n += len(ref_parts[list(batches).index(b)])
        if srcs:
            out[b] = (np.mean(srcs, axis=0), np.mean(tgts, axis=0), n)
    return out


def apply_map(x, src, tgt):
    """Monotone piecewise-linear map; outside the quantile range shift by the end offsets."""
    if par.get("interp", "linear") == "pchip":
        from scipy.interpolate import PchipInterpolator

        # thin knots so that consecutive source quantiles differ (pchip needs strictly increasing x)
        keep = np.concatenate([[True], np.diff(src) > 1e-6])
        y = PchipInterpolator(src[keep], tgt[keep], extrapolate=False)(x)
    else:
        y = np.interp(x, src, tgt)
    lo, hi = x < src[0], x > src[-1]
    y[lo] = x[lo] + (tgt[0] - src[0])
    y[hi] = x[hi] + (tgt[-1] - src[-1])
    return y


def correct(Xin, maps, rows_of_batch, w=None, base=None):
    Y = Xin.copy()
    for b, (src, tgt, _) in maps.items():
        r = rows_of_batch[b]
        for k, j in enumerate(cols):
            # make source quantiles strictly increasing for interp
            s = np.maximum.accumulate(src[:, k] + np.arange(len(src)) * 1e-9)
            y = apply_map(Xin[r, j], s, tgt[:, k])
            if w is None:
                Y[r, j] = y
            else:
                Y[r, j] = w[r] * y + (1 - w[r]) * base[r, j]
    return Y


rows = {b: np.where(batch == b)[0] for b in batches}
isc = ctrl > 0
print("Global control quantile maps", flush=True)
gmaps = qmaps(X[isc][:, cols], batch[isc], ctrl[isc])
Y = correct(X, gmaps, rows)

if par["n_clusters"] and par["n_clusters"] > 0:
    from sklearn.cluster import MiniBatchKMeans

    print("Clustering on globally corrected lineage markers", flush=True)
    L = Y[:, lineage]
    mu, sd = L.mean(0), L.std(0) + 1e-9
    km = MiniBatchKMeans(par["n_clusters"], random_state=par["seed"], batch_size=4096, n_init=3)
    rng = np.random.default_rng(par["seed"])
    fit_idx = rng.choice(len(L), min(len(L), 200000), replace=False)
    km.fit((L[fit_idx] - mu) / sd)
    Ls = (L - mu) / sd
    lab = km.predict(Ls)
    soft = par.get("soft", 0)
    if soft:
        # soft assignment to the 3 nearest centroids, Gaussian weights with bandwidth = soft x median
        # distance to the nearest centroid (smooth transitions between populations)
        C = km.cluster_centers_
        R_idx = np.zeros((len(Ls), 3), np.int32)
        R_w = np.zeros((len(Ls), 3), np.float32)
        for s0 in range(0, len(Ls), 200000):
            D = ((Ls[s0:s0 + 200000, None, :] - C[None]) ** 2).sum(-1)
            ix = np.argsort(D, 1)[:, :3]
            R_idx[s0:s0 + 200000] = ix
            R_w[s0:s0 + 200000] = np.take_along_axis(D, ix, 1)
        tau2 = soft ** 2 * np.median(R_w[:, 0])
        R_w = np.exp(-(R_w - R_w[:, :1]) / (2 * tau2))
        R_w /= R_w.sum(1, keepdims=True)
    else:
        R_idx = lab[:, None].astype(np.int32)
        R_w = np.ones((len(Ls), 1), np.float32)
    Z = np.zeros_like(Y)
    Z[:, ~to_corr] = Y[:, ~to_corr]
    acc_w = np.zeros(len(Y))
    for c in range(par["n_clusters"]):
        cm = (lab == c) & isc
        cmaps = qmaps(X[cm][:, cols], batch[cm], ctrl[cm]) if cm.sum() >= 50 else {}
        hit = (R_idx == c)
        rows_c = np.where(hit.any(1))[0]
        if len(rows_c) == 0:
            continue
        wc = (R_w * hit).sum(1)[rows_c]
        amaps = {}
        if par.get("allcells_weight", 0.0) > 0:
            ma = lab == c
            amaps = qmaps(X[ma][:, cols], batch[ma], np.ones(ma.sum(), int), donor_list=[1])
        for b in batches:
            r = rows_c[batch[rows_c] == b]
            if len(r) == 0:
                continue
            w_r = wc[batch[rows_c] == b]
            Yc = Y[np.ix_(r, cols)].copy()
            if b in cmaps:
                src, tgt, n = cmaps[b]
                wgt = n / (n + par["shrinkage"])
                for k, j in enumerate(cols):
                    s = np.maximum.accumulate(src[:, k] + np.arange(len(src)) * 1e-9)
                    Yc[:, k] = wgt * apply_map(X[r, j], s, tgt[:, k]) + (1 - wgt) * Yc[:, k]
            alpha = par.get("allcells_weight", 0.0)
            if alpha > 0 and b in amaps:
                # all cells of the population (controls and non-controls) as one pseudo-donor
                src, tgt, n = amaps[b]
                for k, j in enumerate(cols):
                    s = np.maximum.accumulate(src[:, k] + np.arange(len(src)) * 1e-9)
                    Yc[:, k] = (1 - alpha) * Yc[:, k] + alpha * apply_map(X[r, j], s, tgt[:, k])
            Z[np.ix_(r, cols)] += w_r[:, None] * Yc
    Y = Z

print("Write output", flush=True)
out = ad.AnnData(
    obs=adata.obs[[]],
    var=adata.var[[]],
    layers={"integrated": Y.astype(np.float32) if adata.layers["preprocessed"].dtype == np.float32 else Y},
    uns={
        "dataset_id": adata.uns["dataset_id"],
        "method_id": meta["name"],
        "parameters": {k: par.get(k) for k in ("n_quantiles", "n_clusters", "shrinkage", "soft", "interp", "target", "allcells_weight", "seed")},
    },
)
out.write_h5ad(par["output"], compression="gzip")
print("Written", out.shape, flush=True)
