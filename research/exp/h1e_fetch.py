#!/usr/bin/env python3
"""H1e data: Otaniemi (Espoo/Helsinki) 2018 EGS stimulation of OTN-3.

Sources (both CC BY 4.0, cite them):
  * Hydraulic data (1-min pressure and rate): Kukkonen et al. (2023) Geothermics, Mendeley Data
    https://doi.org/10.17632/ycdcw2rdfj.1   (file OTN3.xlsx, local Helsinki time = UTC+3 in summer)
  * Seismic catalogue: Leonhardt, Kwiatek, Martinez-Garzon, Heikkinen (2021), GFZ Data Services
    https://doi.org/10.5880/GFZ.4.2.2021.001  (UTC; 61,163 events, ML_HEL and Mw)

Writes compact files to research/data/otaniemi/:
  hydraulic_10min.csv : t_h (hours since 2018-06-04 00:00 UTC), rate_m3min (mean), whp_bar (mean), n (samples)
  catalogue.csv       : t_h, ML, Mw, cls
Raw downloads go to research/data/otaniemi/raw/ and can be deleted afterwards (re-run to fetch again).
"""
import csv, datetime as dt, io, os, sys, urllib.request, zipfile

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data", "otaniemi")
RAW = os.path.join(DATA, "raw")
T0 = dt.datetime(2018, 6, 4, 0, 0, 0)            # UTC origin
URL_XLSX = ("https://data.mendeley.com/public-files/datasets/ycdcw2rdfj/files/"
            "ce823c87-39c8-4116-9608-34bf39790d43/file_downloaded")
URL_ZIP = "https://datapub.gfz.de/download/10.5880.GFZ.4.2.2021.001mnod/2021-001_Leonhardt-et-al_data.zip"


def fetch(url, path):
    if not os.path.exists(path):
        print("downloading", url)
        with urllib.request.urlopen(url, timeout=300) as r, open(path, "wb") as f:
            f.write(r.read())


def hydraulic():
    import openpyxl
    p = os.path.join(RAW, "OTN3.xlsx")
    fetch(URL_XLSX, p)
    wb = openpyxl.load_workbook(p, read_only=True)
    ws = wb["Pumping_data"]
    bins = {}
    for i, r in enumerate(ws.iter_rows(values_only=True)):
        if i == 0 or r[0] is None or r[1] is None:
            continue
        d, tm = r[0], r[1]
        if not isinstance(d, dt.datetime) or not isinstance(tm, dt.time):
            continue
        local = dt.datetime(d.year, d.month, d.day, tm.hour, tm.minute, tm.second)
        utc = local - dt.timedelta(hours=3)          # EEST
        th = (utc - T0).total_seconds() / 3600.0
        try:
            pr, q = float(r[4]), float(r[5])
        except (TypeError, ValueError):
            continue
        k = int(th * 6)                              # 10-min bins
        b = bins.setdefault(k, [0.0, 0.0, 0])
        b[0] += q; b[1] += pr; b[2] += 1
    out = os.path.join(DATA, "hydraulic_10min.csv")
    with open(out, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["t_h", "rate_m3min", "whp_bar", "n"])
        for k in sorted(bins):
            q, pr, n = bins[k]
            w.writerow([f"{(k + 0.5) / 6:.4f}", f"{q / n:.5f}", f"{pr / n:.3f}", n])
    print("wrote", out, len(bins), "bins")


def catalogue():
    zp = os.path.join(RAW, "leonhardt2021.zip")
    fetch(URL_ZIP, zp)
    z = zipfile.ZipFile(zp)
    name = [n for n in z.namelist() if n.endswith(".csv")][0]
    txt = z.read(name).decode("latin-1").splitlines()
    hdr = next(i for i, l in enumerate(txt) if l.startswith("id;class"))
    out = os.path.join(DATA, "catalogue.csv")
    n = 0
    with open(out, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["t_h", "ML", "Mw", "cls"])
        for l in txt[hdr + 1:]:
            c = l.split(";")
            if len(c) < 11 or not c[3].strip():
                continue
            fl = lambda s: float(s.replace(",", "."))
            y, mo, d, h, mi, s = int(c[3]), int(c[4]), int(c[5]), int(c[6]), int(c[7]), fl(c[8])
            t = dt.datetime(y, mo, d, h, mi) + dt.timedelta(seconds=s)
            w.writerow([f"{(t - T0).total_seconds() / 3600.0:.5f}", f"{fl(c[9]):.3f}", f"{fl(c[10]):.3f}", c[1]])
            n += 1
    print("wrote", out, n, "events")


if __name__ == "__main__":
    os.makedirs(RAW, exist_ok=True)
    hydraulic()
    catalogue()
