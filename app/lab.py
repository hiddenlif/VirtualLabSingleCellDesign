"""Synthetic single-cell RNA study: donors, batches, libraries and sequencing depth.

A 40-gene panel is measured in three cell types. Genes G01-G05, G06-G10, G11-G15 are markers
of type_1, type_2, type_3 and are not treatment-responsive. G16-G40 are the 25 test genes.
A treatment (case vs control) changes a few test genes in type_3 only.

Per cell i of type t in a library of donor d placed in batch b
    mean_g = base[t,g] * exp(donor[d,t,g]) * exp(batch[b,g]) * capture[b] * size_i * 2**(lfc[t,g] if case)
    true molecules M_ig ~ NB(mean_g, dispersion 0.15), plus ambient molecules for a contaminated library
    sequencing at cumulative depth r reads per cell observes each molecule with q = 1 - exp(-r / 8000).
More reads on a library add observations of the same cells and do not add cells or donors.
"""

import base64

import numpy as np
from harness import LabError

N_GENES = 40
GENES = [f"G{i + 1:02d}" for i in range(N_GENES)]
TYPES = ["type_1", "type_2", "type_3"]
MARKERS = {"type_1": GENES[0:5], "type_2": GENES[5:10], "type_3": GENES[10:15]}
TEST_GENES = GENES[15:]
N_BATCHES = 4
BATCH_CAPACITY = 4
BUDGET = 130.0
MIN_CELLS, MAX_CELLS = 200, 3000
MIN_READS, MAX_READS, MAX_DEPTH = 1000, 40000, 60000
D0 = 8000.0
POOL = {"case": 12, "control": 12}
SIGMA_DONOR, SIGMA_BATCH, DISPERSION = 0.15, 0.25, 0.15
PILOT = [
    ("P1", "control", 1),
    ("P2", "control", 1),
    ("P3", "case", 2),
    ("P4", "case", 2),
]


def make_case(seed):
    r = np.random.default_rng(seed)
    mu = r.normal(np.log(15.0), 0.7, size=N_GENES)
    base = np.exp(mu[None, :] + r.normal(0, 0.4, size=(3, N_GENES)))
    for t in range(3):
        base[t, :15] = 4.0
        base[t, 5 * t : 5 * t + 5] = 14.0
    lfc = np.zeros((3, N_GENES))
    genes = r.choice(np.arange(15, N_GENES), size=8, replace=False)
    lfc[2, genes] = r.choice([-1.0, 1.0], size=8) * r.uniform(0.7, 1.2, size=8)
    donors = {}
    ids = [p[0] for p in PILOT] + [f"D{i}" for i in range(1, 25)]
    for did in ids:
        donors[did] = {
            "effect": r.normal(0, SIGMA_DONOR, size=(3, N_GENES)).round(4).tolist(),
            "proportions": r.dirichlet(60.0 * np.array([0.62, 0.30, 0.08]))
            .round(4)
            .tolist(),
        }
    return {
        "seed": int(seed),
        "base": base.round(4).tolist(),
        "lfc": lfc.round(4).tolist(),
        "batch_effect": r.normal(0, SIGMA_BATCH, size=(N_BATCHES, N_GENES))
        .round(4)
        .tolist(),
        "capture": np.exp(r.normal(0, 0.15, size=N_BATCHES)).round(4).tolist(),
        "donors": donors,
    }


def encode_counts(a):
    a = np.ascontiguousarray(a, dtype="uint16")
    return {
        "shape": list(a.shape),
        "dtype": "uint16",
        "b64": base64.b64encode(a.tobytes()).decode(),
    }


class Lab:
    def __init__(self, case, budget=BUDGET):
        self.case = case
        self.budget_total = float(budget)
        self.spent = 0.0
        self.rng_n = 0
        self.libs = {}
        self.donor_cond = {}
        self.pool_used = {"case": 0, "control": 0}
        self.log = []
        for i, (did, cond, batch) in enumerate(PILOT):
            self.donor_cond[did] = cond
            self._make_library(did, batch, 500)
            self._sequence(f"L{i + 1}", 4000, charge=False)

    def _rng(self):
        self.rng_n += 1
        return np.random.default_rng([self.case["seed"], 4, self.rng_n])

    def _make_library(self, donor, batch, cells):
        rng = self._rng()
        c = self.case
        d = c["donors"][donor]
        cond = self.donor_cond[donor]
        lid = f"L{len(self.libs) + 1}"
        types = rng.choice(
            3, size=cells, p=np.array(d["proportions"]) / np.sum(d["proportions"])
        )
        base = np.array(c["base"])[types]
        logmu = (
            np.log(base)
            + np.array(d["effect"])[types]
            + np.array(c["batch_effect"])[batch - 1][None, :]
        )
        logmu += np.log(c["capture"][batch - 1]) + 0.25 * rng.standard_normal(
            (cells, 1)
        )
        if cond == "case":
            logmu += np.log(2.0) * np.array(c["lfc"])[types]
        mu = np.exp(logmu)
        ambient = (
            float(rng.uniform(0.0, 0.12))
            if rng.random() > 0.15
            else float(rng.uniform(0.3, 0.5))
        )
        M = rng.poisson(rng.gamma(1.0 / DISPERSION, mu * DISPERSION))
        profile = np.array(c["base"]).mean(axis=0)
        M = M + rng.poisson(ambient * profile[None, :] * np.ones((cells, 1)))
        self.libs[lid] = {
            "donor": donor,
            "condition": cond,
            "batch": batch,
            "cells": cells,
            "types": types,
            "M": M,
            "depth": 0.0,
            "obs": np.zeros((cells, N_GENES), dtype=np.int64),
            "ambient": ambient,
        }
        return lid

    def _sequence(self, lid, reads, charge=True):
        L = self.libs[lid]
        q1 = 1.0 - np.exp(-L["depth"] / D0)
        L["depth"] += reads
        q2 = 1.0 - np.exp(-L["depth"] / D0)
        rng = self._rng()
        remaining = L["M"] - L["obs"]
        L["obs"] = L["obs"] + rng.binomial(remaining, (q2 - q1) / max(1.0 - q1, 1e-12))
        if charge:
            self.spent += L["cells"] * reads / 5e6

    def public_config(self):
        return {
            "genes": GENES,
            "marker_genes": MARKERS,
            "test_genes": TEST_GENES,
            "cell_types": TYPES,
            "batches": N_BATCHES,
            "batch_capacity": BATCH_CAPACITY,
            "budget": self.budget_total,
            "cost": {
                "recruit_donor": "6 + cells/500",
                "new_library_same_donor": "3 + cells/500",
                "sequence": "cells * reads_per_cell / 5e6",
            },
            "cells_range": [MIN_CELLS, MAX_CELLS],
            "reads_range": [MIN_READS, MAX_READS],
            "max_depth": MAX_DEPTH,
            "donors_available": {k: POOL[k] - self.pool_used[k] for k in POOL},
            "pilot_libraries": ["L1", "L2", "L3", "L4"],
        }

    def status(self):
        return {
            "budget_left": round(self.budget_total - self.spent, 4),
            "donors_available": {k: POOL[k] - self.pool_used[k] for k in POOL},
            "batch_load": {
                b: sum(1 for l in self.libs.values() if l["batch"] == b)
                for b in range(1, N_BATCHES + 1)
            },
            "libraries": {
                k: {
                    "donor": l["donor"],
                    "condition": l["condition"],
                    "batch": l["batch"],
                    "cells": l["cells"],
                    "reads_per_cell": l["depth"],
                }
                for k, l in self.libs.items()
            },
        }

    def _check_cost(self, cost):
        if self.spent + cost > self.budget_total + 1e-9:
            raise LabError(
                f"budget exhausted (need {cost:.2f}, have {self.budget_total - self.spent:.2f})"
            )

    def _check_batch(self, batch):
        if (
            not isinstance(batch, int)
            or isinstance(batch, bool)
            or not 1 <= batch <= N_BATCHES
        ):
            raise LabError(f"batch must be an integer from 1 to {N_BATCHES}")
        if sum(1 for l in self.libs.values() if l["batch"] == batch) >= BATCH_CAPACITY:
            raise LabError(f"batch {batch} is full")

    def _check_cells(self, n):
        if (
            not isinstance(n, int)
            or isinstance(n, bool)
            or not MIN_CELLS <= n <= MAX_CELLS
        ):
            raise LabError(f"cells must be an integer in [{MIN_CELLS}, {MAX_CELLS}]")

    def handle(self, action, msg):
        if action == "status":
            return self.status()
        if action == "counts":
            lid = msg.get("library")
            if lid not in self.libs:
                raise LabError(f"unknown library {lid!r}")
            L = self.libs[lid]
            if L["depth"] <= 0:
                raise LabError(f"library {lid} has not been sequenced")
            qc_est = float(
                np.clip(
                    L["ambient"]
                    + np.random.default_rng(
                        [self.case["seed"], 5, int(lid[1:])]
                    ).normal(0, 0.03),
                    0,
                    1,
                )
            )
            return {
                "library": lid,
                "donor": L["donor"],
                "condition": L["condition"],
                "batch": L["batch"],
                "reads_per_cell": L["depth"],
                "counts": encode_counts(L["obs"]),
                "qc": {
                    "median_umi_per_cell": float(np.median(L["obs"].sum(axis=1))),
                    "ambient_fraction_estimate": round(qc_est, 3),
                },
            }
        if action == "recruit":
            cond, batch, cells = (
                msg.get("condition"),
                msg.get("batch"),
                msg.get("cells"),
            )
            if cond not in POOL:
                raise LabError("condition must be case or control")
            self._check_cells(cells)
            self._check_batch(batch)
            if self.pool_used[cond] >= POOL[cond]:
                raise LabError(f"no more {cond} donors available")
            cost = 6.0 + cells / 500.0
            self._check_cost(cost)
            did = f"D{self.pool_used['case'] + self.pool_used['control'] + 1}"
            self.pool_used[cond] += 1
            self.donor_cond[did] = cond
            self.spent += cost
            lid = self._make_library(did, batch, cells)
            self.log.append(
                {
                    "op": "recruit",
                    "donor": did,
                    "condition": cond,
                    "batch": batch,
                    "cells": cells,
                }
            )
            return {
                "donor": did,
                "library": lid,
                "cost": cost,
                "budget_left": round(self.budget_total - self.spent, 4),
            }
        if action == "new_library":
            did, batch, cells = msg.get("donor"), msg.get("batch"), msg.get("cells")
            if did not in self.donor_cond:
                raise LabError(f"unknown donor {did!r}")
            self._check_cells(cells)
            self._check_batch(batch)
            cost = 3.0 + cells / 500.0
            self._check_cost(cost)
            self.spent += cost
            lid = self._make_library(did, batch, cells)
            self.log.append(
                {"op": "new_library", "donor": did, "batch": batch, "cells": cells}
            )
            return {
                "library": lid,
                "cost": cost,
                "budget_left": round(self.budget_total - self.spent, 4),
            }
        if action == "sequence":
            lid, reads = msg.get("library"), msg.get("reads_per_cell")
            if lid not in self.libs:
                raise LabError(f"unknown library {lid!r}")
            if (
                not isinstance(reads, int)
                or isinstance(reads, bool)
                or not MIN_READS <= reads <= MAX_READS
            ):
                raise LabError(
                    f"reads_per_cell must be an integer in [{MIN_READS}, {MAX_READS}]"
                )
            L = self.libs[lid]
            if L["depth"] + reads > MAX_DEPTH:
                raise LabError(f"library depth limit {MAX_DEPTH} reached")
            cost = L["cells"] * reads / 5e6
            self._check_cost(cost)
            self._sequence(lid, reads)
            self.log.append({"op": "sequence", "library": lid, "reads_per_cell": reads})
            return {
                "cost": cost,
                "reads_per_cell": L["depth"],
                "budget_left": round(self.budget_total - self.spent, 4),
            }
        raise LabError(f"unknown action {action!r}")


FOLLOW_UP_STEPS = [
    "recruit_donors",
    "new_libraries",
    "resequence",
    "no_further_acquisition",
]
RULES = {"min_true_called": 5, "max_false_calls": 2}
# Interval coverage (true_covered, null_coverage) is reported but not graded: with the effect range
# disclosed, a fixed interval can meet any coverage rule, and width caps penalised legitimate methods.


def evaluate(lab, result):
    lfc = np.array(lab.case["lfc"])
    out = {
        "spent": round(lab.spent, 3),
        "donors_sequenced": {},
        "batches_used": sorted({l["batch"] for l in lab.libs.values()}),
    }
    for cond in ("case", "control"):
        out["donors_sequenced"][cond] = len(
            {
                l["donor"]
                for l in lab.libs.values()
                if l["condition"] == cond and l["depth"] > 0
            }
        )
    eff = (
        result.get("effects")
        if isinstance(result, dict) and isinstance(result.get("effects"), dict)
        else {}
    )
    tp = fp = cov_true = n_true = n_null_rep = null_cov = n_reported = 0
    missing = 0
    sig_spanning_zero = 0
    for ti, t in enumerate(TYPES):
        for gi in range(15, N_GENES):
            g = GENES[gi]
            truth = lfc[ti, gi]
            try:
                e = eff[t][g]
                est, lo, hi, sig = (
                    float(e["log2fc"]),
                    float(e["ci"][0]),
                    float(e["ci"][1]),
                    bool(e["significant"]),
                )
                if not (np.all(np.isfinite([est, lo, hi])) and lo <= hi):
                    raise ValueError("log2fc and interval must be finite with lo <= hi")
                n_reported += 1
                if sig and lo <= 0.0 <= hi:
                    sig_spanning_zero += 1
            except Exception:
                missing += 1
                est = lo = hi = None
                sig = False
            if truth != 0:
                n_true += 1
                if sig and est is not None and np.sign(est) == np.sign(truth):
                    tp += 1
                if lo is not None and lo <= truth <= hi:
                    cov_true += 1
            else:
                if sig:
                    fp += 1
                if lo is not None:
                    n_null_rep += 1
                    null_cov += lo <= 0.0 <= hi
    led = result.get("ledger") if isinstance(result, dict) else None
    try:
        out["ledger_spent_matches"] = abs(float(led["spent"]) - lab.spent) <= 0.01
        out["ledger_libraries_match"] = sorted(led["libraries"]) == sorted(
            k for k, l in lab.libs.items() if l["depth"] > 0
        )
    except Exception:
        out["ledger_spent_matches"] = out["ledger_libraries_match"] = False
    fu = result.get("follow_up") if isinstance(result, dict) else None
    out["follow_up_valid"] = (
        isinstance(fu, dict)
        and fu.get("next_step") in FOLLOW_UP_STEPS
        and fu.get("cell_type") in TYPES + [None]
        and isinstance(fu.get("rationale"), str)
        and len(fu["rationale"]) > 0
    )
    out.update(
        true_effects=n_true,
        true_called_correct_sign=tp,
        false_calls=fp,
        true_covered=cov_true,
        null_coverage=(null_cov / n_null_rep) if n_null_rep else 0.0,
        reported=n_reported,
        missing_tests=missing,
        significant_with_interval_spanning_0=sig_spanning_zero,
    )
    return out


def check_rules(m):
    why = []
    if not (m["ledger_spent_matches"] and m["ledger_libraries_match"]):
        why.append("ledger does not match the session record")
    if not m["follow_up_valid"]:
        why.append("follow_up missing or malformed")
    if m["missing_tests"]:
        why.append(
            f"{m['missing_tests']} of 75 tests not reported (or not finite with lo <= hi)"
        )
    if m["significant_with_interval_spanning_0"]:
        why.append(
            f"{m['significant_with_interval_spanning_0']} tests called significant have an interval containing 0"
        )
    if m["true_called_correct_sign"] < RULES["min_true_called"]:
        why.append(
            f"only {m['true_called_correct_sign']} of {m['true_effects']} true effects called with the right sign"
        )
    if m["false_calls"] > RULES["max_false_calls"]:
        why.append(f"{m['false_calls']} false calls")
    return not why, why
