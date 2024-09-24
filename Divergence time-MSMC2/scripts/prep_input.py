import pandas as pd

def get_url(sheet_id, sheet_name):
    return f"https://docs.google.com/spreadsheets/d/{sheet_id}/gviz/tq?tqx=out:csv&sheet={sheet_name}"
sheet_id = "1zZgU2F9BYUbnaVicXPBpm-pAMtKYz7Y-cH5tzib_fYU"
sheet_name = "Blad1"

def download_sheet(url):
    return pd.read_csv(url)

url = get_url(sheet_id, sheet_name)

a = download_sheet(url)
outbase = "prep"
b=a[~a["Sample ID"].isna()]
filtercols = [x for x in b.columns if x.startswith("Filter") and "King" not in x]
keep = b[filtercols].sum(1)<1
for name, data  in b[keep].groupby("Clean_groups_K=10"):
    with open(f"{outbase}/{name}.group", 'w') as fh:
        print('\n'.join(data["Sample ID"]), file=fh)
