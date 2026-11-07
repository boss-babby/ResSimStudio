#!/usr/bin/env python3
"""H8' data (Cycle 6): Otaniemi 2020 stimulation of OTN-2 (St1 Deep Heat, Helsinki).

Sources (both CC BY 4.0, cite them):
  * Hydraulic data (1-min WHP and rate, sheet 'Pumping Data'): Kukkonen et al. (2023) Geothermics, Mendeley Data
    https://doi.org/10.17632/ycdcw2rdfj.1  (file OTN2.xlsx, local Helsinki time = UTC+3 in May)
  * Seismic catalogue: Kwiatek, Martinez-Garzon & Karjalainen (2022), GFZ Data Services
    https://doi.org/10.5880/GFZ.4.2.2022.001  (Catalog_St1-OTN2-2020.csv, UTC, 6,318 events, ML)
Writes research/data/otaniemi2020/{hydraulic_10min.csv, catalogue.csv} (same format as exp/h1e_fetch.py);
events before the time origin (2020-05-05 00:00 UTC) are dropped.  Raw files go to raw/ (deletable).
"""
import csv, datetime as dt, os, urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data", "otaniemi2020")
RAW = os.path.join(DATA, "raw")
T0 = dt.datetime(2020, 5, 5, 0, 0, 0)
URL_XLSX = ("https://data.mendeley.com/public-files/datasets/ycdcw2rdfj/files/"
            "dd7fdb3c-10e8-4504-9999-63d29318dd7c/file_downloaded")
URL_CAT = "https://datapub.gfz.de/download/10.5880.GFZ.4.2.2022.001KheZV/Catalog_St1-OTN2-2020.csv"


def fetch(url, path):
    if not os.path.exists(path):
        print("downloading", url)
        req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
        with urllib.request.urlopen(req, timeout=300) as r, open(path, "wb") as f:
            f.write(r.read())


def hydraulic():
    import openpyxl
    p = os.path.join(RAW, "OTN2.xlsx")
    fetch(URL_XLSX, p)
    ws = openpyxl.load_workbook(p, read_only=True)["Pumping Data"]
    bins = {}
    for i, r in enumerate(ws.iter_rows(values_only=True)):
        if i == 0 or not isinstance(r[0], dt.datetime) or r[1] is None:
            continue
        tm = r[1]
        if isinstance(tm, str):
            hh, mm, ss = (int(x) for x in tm.split(":"))
        else:
            hh, mm, ss = tm.hour, tm.minute, tm.second
        local = dt.datetime(r[0].year, r[0].month, r[0].day, hh, mm, ss)
        th = (local - dt.timedelta(hours=3) - T0).total_seconds() / 3600.0
        try:
            pr, q = float(r[4]), float(r[5])
        except (TypeError, ValueError):
            continue
        b = bins.setdefault(int(th * 6), [0.0, 0.0, 0])
        b[0] += max(q, 0.0); b[1] += pr; b[2] += 1
    out = os.path.join(DATA, "hydraulic_10min.csv")
    with open(out, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["t_h", "rate_m3min", "whp_bar", "n"])
        for k in sorted(bins):
            q, pr, n = bins[k]
            w.writerow([f"{(k + 0.5) / 6:.4f}", f"{q / n:.5f}", f"{pr / n:.3f}", n])
    print("wrote", out, len(bins), "bins")


def catalogue():
    p = os.path.join(RAW, "Catalog_St1-OTN2-2020.csv")
    fetch(URL_CAT, p)
    out = os.path.join(DATA, "catalogue.csv")
    n = 0
    with open(p, encoding="utf-8-sig") as fi, open(out, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["t_h", "ML", "Mw", "cls"])
        for l in fi.read().splitlines()[1:]:
            c = l.split(";")
            if len(c) < 3:
                continue
            t = dt.datetime.strptime(c[1].rstrip("Z"), "%Y-%m-%dT%H:%M:%S.%f")
            th = (t - T0).total_seconds() / 3600.0
            if th < 0:
                continue
            ml = float(c[2].replace(",", "."))
            w.writerow([f"{th:.5f}", f"{ml:.3f}", "", "res"])
            n += 1
    print("wrote", out, n, "events")


if __name__ == "__main__":
    os.makedirs(RAW, exist_ok=True)
    hydraulic()
    catalogue()
