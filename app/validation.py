"""请求体校验：所有越界 / 非法输入都汇总为明确的字段错误。

校验通过后返回规范化数据，供求解器直接使用；校验失败时返回完整的
字段错误列表（一次性报告所有问题，而不是只报第一个）。
"""

from __future__ import annotations

MIN_SEGMENTS = 6
MAX_SEGMENTS = 12
MIN_WINDOWS = 8
MAX_WINDOWS = 20
MAX_SEGMENT_LENGTH = 1_000_000          # 单段长度上限（与应变单位配套）
MAX_STRAIN_ABS = 1_000_000            # 整数微应变的绝对值上限
# 单窗长度加权应变和的理论上限：12 段 × 10^6（长度）× 10^6（应变）
MAX_ELONGATION_ABS = MAX_SEGMENTS * MAX_SEGMENT_LENGTH * MAX_STRAIN_ABS


def _is_int(value) -> bool:
    # JSON 的 true/false 在 Python 里是 bool，属于 int 的子类，必须排除。
    return isinstance(value, int) and not isinstance(value, bool)


def _err(errors: list[dict], field: str, message: str) -> None:
    errors.append({"field": field, "message": message})


def validate_payload(payload) -> tuple[dict | None, list[dict]]:
    """校验请求体。

    返回 (规范化数据, 错误列表)。存在任何错误时规范化数据为 None。
    规范化数据形如：
        {
            "segment_lengths": [int, ...],
            "strain_min": int, "strain_max": int,
            "windows": [{"start": int, "end": int, "min": int, "max": int}, ...],
        }
    """
    errors: list[dict] = []
    if not isinstance(payload, dict):
        return None, [{"field": "(body)", "message": "request body must be a JSON object"}]

    # ---- segment_lengths ----
    lengths = payload.get("segment_lengths")
    n = None  # 段数合法时才有确定的下标范围，用于窗口 start/end 检查
    if "segment_lengths" not in payload:
        _err(errors, "segment_lengths", "field is required")
    elif not isinstance(lengths, list):
        _err(errors, "segment_lengths", "must be an array of integers")
    elif not (MIN_SEGMENTS <= len(lengths) <= MAX_SEGMENTS):
        _err(
            errors,
            "segment_lengths",
            f"must contain {MIN_SEGMENTS}..{MAX_SEGMENTS} segments, got {len(lengths)}",
        )
    else:
        n = len(lengths)
        for i, v in enumerate(lengths):
            if not _is_int(v):
                _err(errors, f"segment_lengths[{i}]", "must be an integer")
            elif not (1 <= v <= MAX_SEGMENT_LENGTH):
                _err(errors, f"segment_lengths[{i}]", f"must be within 1..{MAX_SEGMENT_LENGTH}")

    # ---- strain_bounds ----
    bounds = payload.get("strain_bounds")
    strain_min = strain_max = None
    if "strain_bounds" not in payload:
        _err(errors, "strain_bounds", "field is required")
    elif not isinstance(bounds, dict):
        _err(errors, "strain_bounds", "must be an object with integer fields min and max")
    else:
        smin = bounds.get("min")
        smax = bounds.get("max")
        smin_ok = smax_ok = False
        if not _is_int(smin):
            _err(errors, "strain_bounds.min", "must be an integer")
        elif not (-MAX_STRAIN_ABS <= smin <= MAX_STRAIN_ABS):
            _err(errors, "strain_bounds.min", f"must be within -{MAX_STRAIN_ABS}..{MAX_STRAIN_ABS}")
        else:
            smin_ok = True
        if not _is_int(smax):
            _err(errors, "strain_bounds.max", "must be an integer")
        elif not (-MAX_STRAIN_ABS <= smax <= MAX_STRAIN_ABS):
            _err(errors, "strain_bounds.max", f"must be within -{MAX_STRAIN_ABS}..{MAX_STRAIN_ABS}")
        else:
            smax_ok = True
        if smin_ok and smax_ok:
            if smin > smax:
                _err(errors, "strain_bounds", "min must be <= max")
            else:
                strain_min, strain_max = smin, smax

    # ---- windows ----
    windows = payload.get("windows")
    if "windows" not in payload:
        _err(errors, "windows", "field is required")
    elif not isinstance(windows, list):
        _err(errors, "windows", "must be an array of observation windows")
    elif not (MIN_WINDOWS <= len(windows) <= MAX_WINDOWS):
        _err(
            errors,
            "windows",
            f"must contain {MIN_WINDOWS}..{MAX_WINDOWS} windows, got {len(windows)}",
        )
    else:
        for j, w in enumerate(windows):
            prefix = f"windows[{j}]"
            if not isinstance(w, dict):
                _err(errors, prefix, "must be an object with start/end/elongation")
                continue
            start = w.get("start")
            end = w.get("end")
            start_ok = end_ok = False
            if not _is_int(start):
                _err(errors, f"{prefix}.start", "must be an integer")
            elif n is not None and not (0 <= start < n):
                _err(errors, f"{prefix}.start", f"must be within 0..{n - 1}")
            else:
                start_ok = True
            if not _is_int(end):
                _err(errors, f"{prefix}.end", "must be an integer")
            elif n is not None and not (0 <= end < n):
                _err(errors, f"{prefix}.end", f"must be within 0..{n - 1}")
            else:
                end_ok = True
            if start_ok and end_ok and start > end:
                _err(errors, prefix, "start must be <= end")

            elong = w.get("elongation")
            if not isinstance(elong, dict):
                _err(errors, f"{prefix}.elongation", "must be an object with integer fields min and max")
                continue
            emin = elong.get("min")
            emax = elong.get("max")
            emin_ok = emax_ok = False
            if not _is_int(emin):
                _err(errors, f"{prefix}.elongation.min", "must be an integer")
            elif not (-MAX_ELONGATION_ABS <= emin <= MAX_ELONGATION_ABS):
                _err(errors, f"{prefix}.elongation.min", f"must be within ±{MAX_ELONGATION_ABS}")
            else:
                emin_ok = True
            if not _is_int(emax):
                _err(errors, f"{prefix}.elongation.max", "must be an integer")
            elif not (-MAX_ELONGATION_ABS <= emax <= MAX_ELONGATION_ABS):
                _err(errors, f"{prefix}.elongation.max", f"must be within ±{MAX_ELONGATION_ABS}")
            else:
                emax_ok = True
            if emin_ok and emax_ok and emin > emax:
                _err(errors, f"{prefix}.elongation", "min must be <= max")

    if errors:
        return None, errors

    data = {
        "segment_lengths": list(lengths),
        "strain_min": strain_min,
        "strain_max": strain_max,
        "windows": [
            {
                "start": w["start"],
                "end": w["end"],
                "min": w["elongation"]["min"],
                "max": w["elongation"]["max"],
            }
            for w in windows
        ],
    }
    return data, []
