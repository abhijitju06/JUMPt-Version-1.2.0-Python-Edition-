#!/usr/bin/env python3
import argparse, os, re, sys

def _read_unit(params_path:str)->str:
    unit = "days"
    try:
        with open(params_path, "r", encoding="utf-8") as f:
            for raw in f:
                s = raw.strip()
                if not s or "=" not in s or s.startswith(("#","%",";")):
                    continue
                k, v = s.split("=", 1)
                k = k.strip().lower()
                v = re.split(r"[#%]", v, 1)[0].strip().strip(";").strip("'").strip('"').lower()
                if k in ("pulse_time_unit", "time_unit", "timeunit"):
                    unit = v
                    break
    except Exception:
        pass
    return unit

def main(params_path:str):
    unit = _read_unit(params_path)
    if unit.startswith("h"):   # hours
        import jumpt_hours as runner
    else:                      # default: days 
        import jumpt_days as runner
    runner.main(params_path)

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--params", default="JUMPt_main.params")
    args = ap.parse_args()
    main(args.params)
