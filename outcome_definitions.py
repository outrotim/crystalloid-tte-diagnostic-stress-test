"""Pure outcome-definition functions used in the Study 14 stress test.

This module contains no database access and accepts caller-supplied, tabular
measurements only. Missing composite components remain missing unless another
component establishes the event.
"""

from __future__ import annotations

import pandas as pd


def followup_end(time_zero, discharge_time):
    start = pd.Timestamp(time_zero)
    day30 = start + pd.Timedelta(days=30)
    return day30 if pd.isna(discharge_time) else min(pd.Timestamp(discharge_time), day30)


def death_in_followup(death_time, start, end) -> int:
    if pd.isna(death_time):
        return 0
    return int(pd.Timestamp(start) <= pd.Timestamp(death_time) <= pd.Timestamp(end))


def _clean_creatinine_labs(labs: pd.DataFrame) -> pd.DataFrame:
    required = {"charttime", "valuenum"}
    missing = required - set(labs.columns)
    if missing:
        raise ValueError(f"Creatinine data missing columns: {sorted(missing)}")
    data = labs.copy()
    data["charttime"] = pd.to_datetime(data["charttime"], errors="coerce")
    data["valuenum"] = pd.to_numeric(data["valuenum"], errors="coerce")
    return data.loc[
        data["charttime"].notna() & data["valuenum"].notna() & data["valuenum"].gt(0)
    ].sort_values("charttime")


def select_baseline_creatinine(labs: pd.DataFrame, time_zero, index_admittime):
    """Select measured baseline creatinine without a post-time-zero fallback."""
    data = _clean_creatinine_labs(labs)
    if "stable_prior" not in data:
        data["stable_prior"] = False
    t0 = pd.Timestamp(time_zero)
    admit = pd.Timestamp(index_admittime)
    prior = data.loc[
        data["stable_prior"].fillna(False).astype(bool)
        & data["charttime"].between(
            t0 - pd.Timedelta(days=365), t0 - pd.Timedelta(days=7), inclusive="both"
        )
    ]
    if not prior.empty:
        measured_at = prior["charttime"].max()
        tied = prior.loc[prior["charttime"].eq(measured_at), "valuenum"]
        return float(tied.median()), "stable_prior_7_365d", pd.Timestamp(measured_at)
    current = data.loc[data["charttime"].between(admit, t0, inclusive="both")]
    if not current.empty:
        minimum = current["valuenum"].min()
        measured_at = current.loc[current["valuenum"].eq(minimum), "charttime"].max()
        return float(minimum), "current_admission_pretreatment_min", pd.Timestamp(measured_at)
    return None, "missing", None


def persistent_renal_dysfunction(labs, baseline, start, end) -> dict:
    data = _clean_creatinine_labs(labs)
    if baseline is None or pd.isna(baseline) or float(baseline) <= 0:
        return {"persistent_renal_dysfunction": None, "final_creatinine": None,
                "final_creatinine_time": None, "peak_creatinine": None}
    eligible = data.loc[data["charttime"].between(
        pd.Timestamp(start), pd.Timestamp(end), inclusive="both"
    )]
    if eligible.empty:
        return {"persistent_renal_dysfunction": None, "final_creatinine": None,
                "final_creatinine_time": None, "peak_creatinine": None}
    final_time = eligible["charttime"].max()
    final_value = float(eligible.loc[eligible["charttime"].eq(final_time), "valuenum"].median())
    return {
        "persistent_renal_dysfunction": int(final_value >= 2.0 * float(baseline)),
        "final_creatinine": final_value,
        "final_creatinine_time": pd.Timestamp(final_time),
        "peak_creatinine": float(eligible["valuenum"].max()),
    }


def classify_new_rrt(rrt_times, start, end) -> dict:
    times = sorted(pd.Timestamp(value) for value in rrt_times if not pd.isna(value))
    t0, stop = pd.Timestamp(start), pd.Timestamp(end)
    if any(value < t0 for value in times):
        return {"prior_rrt": 1, "new_rrt": None, "first_new_rrt_time": None}
    eligible = [value for value in times if t0 <= value <= stop]
    return {"prior_rrt": 0, "new_rrt": int(bool(eligible)),
            "first_new_rrt_time": min(eligible) if eligible else None}


def classify_kdigo_creatinine(labs, baseline, start) -> dict:
    data = _clean_creatinine_labs(labs)
    if baseline is None or pd.isna(baseline) or float(baseline) <= 0:
        return {"aki": None, "aki_stage": None, "aki_first_time": None}
    t0 = pd.Timestamp(start)
    followup = data.loc[data["charttime"].between(t0, t0 + pd.Timedelta(days=7), inclusive="both")]
    if followup.empty:
        return {"aki": None, "aki_stage": None, "aki_first_time": None}
    series = pd.concat([
        pd.DataFrame({"charttime": [t0], "valuenum": [float(baseline)]}),
        followup[["charttime", "valuenum"]],
    ], ignore_index=True)
    series = series.groupby("charttime", as_index=False)["valuenum"].median().sort_values("charttime")
    max_stage, first_time = 0, None
    for row in series.itertuples(index=False):
        current_time, current = pd.Timestamp(row.charttime), float(row.valuenum)
        ratio = current / float(baseline)
        prior48 = series.loc[
            series["charttime"].between(current_time - pd.Timedelta(hours=48), current_time,
                                         inclusive="both")
            & series["charttime"].lt(current_time), "valuenum"
        ]
        absolute_rise = current - float(prior48.min()) if not prior48.empty else 0.0
        if ratio >= 3.0 or (current >= 4.0 and current - float(baseline) >= 0.3):
            stage = 3
        elif ratio >= 2.0:
            stage = 2
        elif ratio >= 1.5 or absolute_rise >= 0.3 - 1e-12:
            stage = 1
        else:
            stage = 0
        if stage > 0 and first_time is None:
            first_time = current_time
        max_stage = max(max_stage, stage)
    return {"aki": int(max_stage > 0), "aki_stage": int(max_stage),
            "aki_first_time": first_time}


def combine_make30(death, new_rrt, persistent):
    components = [death, new_rrt, persistent]
    if any(value == 1 for value in components if value is not None and not pd.isna(value)):
        return 1
    if all(value is not None and not pd.isna(value) and value == 0 for value in components):
        return 0
    return None
