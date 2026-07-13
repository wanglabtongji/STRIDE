#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Integration smoke tests for STRIDE.

Runs the STRIDE command-line tool end-to-end on the demo data bundled in this
directory and validates the outputs. Intended to be run after code changes to
catch regressions -- in particular the empty ``*_spot_celltype_frac.txt`` bug
fixed in v1.0.4 (a zero-sum topic column used to propagate NaN and produce an
empty deconvolution result).

Two cases are available:

* ``embryo`` (default, ~1-2 min): mouse embryo GEO-seq data. Plain-text input,
  automatic marker selection. Exercises ``deconvolve``.
* ``heart`` (~5 min): developing human heart data. 10X HDF5 input, custom
  marker list. Exercises ``deconvolve`` and ``map``, and checks that the mapped
  cell counts stay consistent with the deconvolved fractions.

Usage::

    python test/run_tests.py                 # fast case (mouse embryo)
    python test/run_tests.py --case heart    # human heart (deconvolve + map)
    python test/run_tests.py --case all      # both
    python test/run_tests.py --keep          # keep the output files
    python test/run_tests.py --stride /path/to/STRIDE   # override CLI location

Exit code is 0 if every check passes, 1 otherwise -- suitable for CI.
"""

import argparse
import glob
import os
import shutil
import subprocess
import sys
import tempfile

import numpy as np
import pandas as pd

TEST_DIR = os.path.dirname(os.path.abspath(__file__))

# ANSI colours (disabled when output is not a terminal).
_TTY = sys.stdout.isatty()
GREEN = "\033[32m" if _TTY else ""
RED = "\033[31m" if _TTY else ""
BOLD = "\033[1m" if _TTY else ""
RESET = "\033[0m" if _TTY else ""


class Checker(object):
    """Collects named pass/fail checks and prints them as they run."""

    def __init__(self):
        self.failures = 0
        self.total = 0

    def check(self, name, ok, detail=""):
        self.total += 1
        if not ok:
            self.failures += 1
        mark = GREEN + "PASS" + RESET if ok else RED + "FAIL" + RESET
        line = "  [%s] %s" % (mark, name)
        if detail:
            line += " -- %s" % detail
        print(line)
        return ok


def run_stride(stride, args, log_path):
    """Run a STRIDE subcommand, streaming stdout/stderr to a log file.

    Returns the process return code. The (verbose) STRIDE/gensim output is kept
    in ``log_path`` so the test output stays readable; the log tail is shown
    only when the command fails.
    """
    cmd = [stride] + args
    print("  $ " + " ".join(cmd))
    with open(log_path, "w") as log:
        proc = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT)
    if proc.returncode != 0:
        print(RED + "  command failed (exit %d); last log lines:" % proc.returncode + RESET)
        with open(log_path) as log:
            for line in log.readlines()[-20:]:
                print("    " + line.rstrip())
    return proc.returncode


def validate_frac(chk, frac_path, prefix):
    """Validate a ``*_spot_celltype_frac.txt`` deconvolution result.

    This is the core regression guard: the file must exist, be non-empty,
    contain no NaN, have every spot assigned, and (because we run with
    ``--normalize``) have per-spot fractions that sum to 1.
    """
    ok = chk.check("%s: frac file exists and non-empty" % prefix,
                   os.path.isfile(frac_path) and os.path.getsize(frac_path) > 0,
                   frac_path)
    if not ok:
        return None

    df = pd.read_csv(frac_path, sep="\t", index_col=0)
    chk.check("%s: frac has spots and celltypes" % prefix,
              df.shape[0] > 0 and df.shape[1] > 0,
              "shape = %d spots x %d celltypes" % (df.shape[0], df.shape[1]))
    chk.check("%s: no NaN in fractions" % prefix,
              not bool(df.isna().any().any()))
    n_empty = int((df.sum(axis=1) == 0).sum())
    chk.check("%s: no empty (all-zero) spots" % prefix,
              n_empty == 0, "%d empty spots" % n_empty)
    row_sums = df.sum(axis=1)
    chk.check("%s: per-spot fractions sum to 1" % prefix,
              bool(np.allclose(row_sums.values, 1.0, atol=1e-3)),
              "row-sum range %.4f .. %.4f" % (row_sums.min(), row_sums.max()))
    chk.check("%s: fractions within [0, 1]" % prefix,
              bool(df.values.min() >= -1e-9 and df.values.max() <= 1.0 + 1e-9),
              "value range %.4f .. %.4f" % (df.values.min(), df.values.max()))
    return df


def case_deconvolve(chk, stride, outdir, prefix, sc_count, sc_celltype,
                    st_count, gene_use=None, ntopics=None):
    """Run ``STRIDE deconvolve`` and validate the fraction output."""
    print(BOLD + "\n== deconvolve: %s ==" % prefix + RESET)
    args = ["deconvolve",
            "--sc-count", sc_count,
            "--sc-celltype", sc_celltype,
            "--st-count", st_count,
            "--outdir", outdir,
            "--outprefix", prefix,
            "--normalize"]
    if gene_use:
        args += ["--gene-use", gene_use]
    if ntopics:
        args += ["--ntopics"] + [str(n) for n in ntopics]

    log_path = os.path.join(outdir, "%s_deconvolve.log" % prefix)
    rc = run_stride(stride, args, log_path)
    if not chk.check("%s: deconvolve exits 0" % prefix, rc == 0):
        return
    validate_frac(chk, os.path.join(outdir, "%s_spot_celltype_frac.txt" % prefix), prefix)


def case_map(chk, stride, outdir, prefix, sc_celltype):
    """Run ``STRIDE map`` on the deconvolve outputs and validate consistency.

    Depends on the deconvolve case for the same ``prefix`` having run first.
    """
    print(BOLD + "\n== map: %s ==" % prefix + RESET)
    frac_path = os.path.join(outdir, "%s_spot_celltype_frac.txt" % prefix)
    topic_files = glob.glob(os.path.join(outdir, "%s_topic_spot_mat_*.txt" % prefix))
    model_dir = os.path.join(outdir, "model")

    ok = chk.check("%s: deconvolve outputs present for map" % prefix,
                   os.path.isfile(frac_path) and len(topic_files) == 1 and os.path.isdir(model_dir),
                   "topic files found: %d" % len(topic_files))
    if not ok:
        return

    args = ["map",
            "--topic-spot-mat", topic_files[0],
            "--sc-celltype", sc_celltype,
            "--spot-celltype-frac", frac_path,
            "--model-dir", model_dir,
            "--outdir", outdir,
            "--outprefix", prefix]
    log_path = os.path.join(outdir, "%s_map.log" % prefix)
    rc = run_stride(stride, args, log_path)
    if not chk.check("%s: map exits 0" % prefix, rc == 0):
        return

    counts_path = os.path.join(outdir, "%s_spot_mapping_celltype_counts.txt" % prefix)
    sim_pattern = glob.glob(os.path.join(outdir, "%s_spot_mapping_similar_*_cell.txt" % prefix))
    chk.check("%s: mapping count + similar-cell files exist" % prefix,
              os.path.isfile(counts_path) and len(sim_pattern) == 1)
    if not os.path.isfile(counts_path):
        return

    counts = pd.read_csv(counts_path, sep="\t", index_col=0)
    chk.check("%s: no NaN in mapping counts" % prefix,
              not bool(counts.isna().any().any()))

    # The v1.0.3/v1.0.4 fix makes mapping consistent with the deconvolved
    # fractions: per-spot mapped cell proportions should track the fractions.
    frac = pd.read_csv(frac_path, sep="\t", index_col=0)
    common = frac.index.intersection(counts.index)
    if chk.check("%s: mapping and frac share spots" % prefix, len(common) > 0,
                 "%d common spots" % len(common)):
        f = frac.loc[common]
        totals = counts.loc[common].sum(axis=1).replace(0, np.nan)
        c = counts.loc[common].div(totals, axis=0).fillna(0)
        r = np.corrcoef(f.values.ravel(), c.values.ravel())[0, 1]
        chk.check("%s: mapping tracks fractions (corr > 0.5)" % prefix,
                  r > 0.5, "corr = %.3f" % r)


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--case", choices=["embryo", "heart", "all"], default="embryo",
                        help="Which test case(s) to run. DEFAULT: embryo (fast).")
    parser.add_argument("--stride", default=None,
                        help="Path to the STRIDE executable. DEFAULT: found on PATH.")
    parser.add_argument("--outdir", default=None,
                        help="Directory for test outputs. DEFAULT: a temp dir (removed on exit).")
    parser.add_argument("--keep", action="store_true",
                        help="Keep the output files instead of deleting them.")
    parser.add_argument("--ntopics", type=int, nargs="+", default=None,
                        help="Restrict the topic-number search (speeds up the run).")
    args = parser.parse_args()

    stride = args.stride or shutil.which("STRIDE")
    if not stride or not os.path.exists(stride):
        print(RED + "ERROR: STRIDE executable not found. Install it (pip install stridespatial) "
              "or pass --stride /path/to/STRIDE." + RESET)
        return 2
    print("Using STRIDE: %s" % stride)
    try:
        ver = subprocess.check_output([stride, "--version"], stderr=subprocess.STDOUT).decode().strip()
        print("STRIDE version: %s" % ver)
    except Exception as exc:  # pragma: no cover - informational only
        print("Could not read STRIDE version: %s" % exc)

    if args.outdir:
        outdir = os.path.abspath(args.outdir)
        os.makedirs(outdir, exist_ok=True)
        cleanup = False
    else:
        outdir = tempfile.mkdtemp(prefix="stride_test_")
        cleanup = not args.keep
    print("Output directory: %s%s" % (outdir, "" if cleanup else "  (kept)"))

    chk = Checker()
    try:
        if args.case in ("embryo", "all"):
            case_deconvolve(
                chk, stride, outdir, "GEO-seq_E7.5",
                sc_count=os.path.join(TEST_DIR, "scRNA_E7.5_gene_count.txt"),
                sc_celltype=os.path.join(TEST_DIR, "scRNA_E7.5_lineage.txt"),
                st_count=os.path.join(TEST_DIR, "GEO-seq_E7.5_gene_count.txt"),
                ntopics=args.ntopics)

        if args.case in ("heart", "all"):
            case_deconvolve(
                chk, stride, outdir, "Human_heart",
                sc_count=os.path.join(TEST_DIR, "Human_heart_scRNA_gene_count.h5"),
                sc_celltype=os.path.join(TEST_DIR, "Human_heart_scRNA_celltype_curated.txt"),
                st_count=os.path.join(TEST_DIR, "Human_heart_ST_gene_count.h5"),
                gene_use=os.path.join(TEST_DIR, "Human_heart_scRNA_markers.txt"),
                ntopics=args.ntopics)
            case_map(chk, stride, outdir, "Human_heart",
                     sc_celltype=os.path.join(TEST_DIR, "Human_heart_scRNA_celltype_curated.txt"))
    finally:
        if cleanup:
            shutil.rmtree(outdir, ignore_errors=True)

    print(BOLD + "\n== summary ==" + RESET)
    passed = chk.total - chk.failures
    colour = GREEN if chk.failures == 0 else RED
    print(colour + "  %d/%d checks passed" % (passed, chk.total) + RESET)
    return 0 if chk.failures == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
