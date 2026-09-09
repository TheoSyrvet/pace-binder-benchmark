#!/usr/bin/env python3
import json, os, glob, csv, argparse
import numpy as np

BINDER_INDEX = {
    "spike":0,"glp1r":0,"il7ra":0,"cd28":0,
    "mdm2":1,"pdl1":1,"egfr":1,"insulinr":1,"stat3":1,"kit":1,
    "pdgfrb":1,"fgfr2":1,"trka":1,"ha_head":1,"clathrin":1,"ccr5":1,"cxcr4":1,
    "tnfa":2,"vegf":2,"ha_stem":2
}
BINDER_LEN = {
    "mdm2":76,"pdl1":131,"spike":130,"tnfa":118,"egfr":129,"insulinr":64,
    "vegf":109,"stat3":148,"kit":123,"glp1r":124,"pdgfrb":55,"fgfr2":80,
    "trka":90,"ha_head":88,"ha_stem":121,"il7ra":130,"clathrin":137,
    "cd28":69,"ccr5":149,"cxcr4":130
}
TOTAL_LEN = {
    "mdm2":161,"pdl1":240,"spike":324,"tnfa":382,"egfr":319,"insulinr":177,
    "vegf":297,"stat3":249,"kit":591,"glp1r":235,"pdgfrb":208,"fgfr2":214,
    "trka":191,"ha_head":299,"ha_stem":377,"il7ra":316,"clathrin":497,
    "cd28":187,"ccr5":441,"cxcr4":432
}

IPTM_THRESH  = 0.5
PTM_THRESH   = 0.55
PLDDT_THRESH = 0.8
IPAE_THRESH  = 10.85

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--target", required=True)
    parser.add_argument("--colabfold_dir", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    target = args.target
    bidx = BINDER_INDEX[target]
    blen = BINDER_LEN[target]
    tlen = TOTAL_LEN[target]

    json_files = glob.glob(f"{args.colabfold_dir}/*/*scores_rank_001*.json")
    print(f"[{target}] {len(json_files)} JSON trouvés", flush=True)

    rows = []
    for jf in sorted(json_files):
        try:
            d = json.load(open(jf))
            iptm       = d.get("iptm", 0)
            ptm        = d.get("ptm", 0)
            plddt      = d.get("plddt", [])
            pae        = d.get("pae", [])
            mean_plddt = sum(plddt) / len(plddt) if plddt else 0

            if pae:
                pae_arr = np.array(pae)
                n = pae_arr.shape[0]
                if bidx == 0:
                    b_start, b_end = 0, blen
                    t_start, t_end = blen, n
                elif bidx == 1:
                    t_actual = tlen - blen
                    t_start, t_end = 0, t_actual
                    b_start, b_end = t_actual, n
                else:
                    b_start, b_end = n - blen, n
                    t_start, t_end = 0, n - blen
                block1 = pae_arr[b_start:b_end, t_start:t_end]
                block2 = pae_arr[t_start:t_end, b_start:b_end]
                ipae   = float(np.mean([block1.mean(), block2.mean()]))
                ipsae  = ipae / 31.0
            else:
                ipae = 999; ipsae = 999

            if not (iptm > IPTM_THRESH and ptm > PTM_THRESH and
                    mean_plddt > PLDDT_THRESH and ipae < IPAE_THRESH):
                continue

            design_name = os.path.basename(os.path.dirname(jf))
            rows.append({
                "name":       design_name,
                "tool":       "rfdiff_mpnn",
                "iptm":       round(iptm, 4),
                "ptm":        round(ptm, 4),
                "mean_plddt": round(mean_plddt, 4),
                "i_pae":      round(ipae, 4),
                "ipsae":      round(ipsae, 4),
            })
        except Exception as e:
            print(f"  ERREUR {jf}: {e}", flush=True)

    print(f"[{target}] {len(rows)} designs passent les filtres AF2", flush=True)
    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    with open(args.output, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["name","tool","iptm","ptm","mean_plddt","i_pae","ipsae"])
        writer.writeheader()
        writer.writerows(rows)
    print(f"[{target}] CSV écrit -> {args.output}", flush=True)

if __name__ == "__main__":
    main()
