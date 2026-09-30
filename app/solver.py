"""整数微应变联合反演器。

海缆检修组只能取得多个重叠标距（观测窗）的累计伸长量，本模块把
所有观测窗联合起来，反演出连续缆段的整数微应变序列。

可行解按以下顺序依次优化（字典序分层目标）：
1. 相邻段应变差绝对值的最大值（一级平滑指标）最小；
2. 全部相邻差绝对值之和（二级平滑指标）最小；
3. 应变序列本身的字典序最小。

每个阶段都是一个规模很小的整数线性规划（ILP），由 HiGHS 求到
可证明最优（mip_rel_gap=0）。求解器内部使用浮点，但最终结果在
返回前一律取整并用 Python 精确整数重新校验全部约束与两级指标，
保证对外暴露的每个数字都能由提交数据直接复核。
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.optimize import Bounds, LinearConstraint, milp

# 单个 ILP 的求解时间上限（秒）。问题规模极小（<=24 个变量、<=65 行），
# 正常求解以毫秒计，该上限只是兜底。
SOLVER_TIME_LIMIT = 10.0

_INT_TOL = 1e-6


class InfeasibleError(Exception):
    """观测彼此冲突：不存在同时满足所有观测窗闭区间的整数应变序列。"""


class SolverFailureError(Exception):
    """底层 MIP 求解器未能给出可通过精确整数复核的最优解。"""


@dataclass(frozen=True)
class InversionResult:
    """反演结果。所有字段均为精确整数。"""

    strains: list[int]            # 逐段整数微应变
    max_adjacent_diff: int        # 一级平滑指标：相邻段应变差绝对值的最大值
    sum_adjacent_diff: int        # 二级平滑指标：相邻段应变差绝对值之和
    window_sums: list[int]        # 逐窗回算的长度加权应变和


def _run_mip(c, rows, row_lo, row_hi, var_lo, var_hi):
    """调用 HiGHS 求解一个纯整数线性规划，返回原始结果对象。"""
    return milp(
        c=np.asarray(c, dtype=np.float64),
        integrality=np.ones(len(c), dtype=np.int64),
        bounds=Bounds(
            np.asarray(var_lo, dtype=np.float64),
            np.asarray(var_hi, dtype=np.float64),
        ),
        constraints=LinearConstraint(
            np.asarray(rows, dtype=np.float64),
            np.asarray(row_lo, dtype=np.float64),
            np.asarray(row_hi, dtype=np.float64),
        ),
        options={"time_limit": SOLVER_TIME_LIMIT, "mip_rel_gap": 0.0, "disp": False},
    )


def _require_optimal(res, stage: str):
    if res.status != 0 or res.x is None:
        raise SolverFailureError(
            f"{stage}: MIP solver did not reach a proven optimum "
            f"(status={res.status}, message={res.message!r})"
        )


def _as_exact_int(value: float, what: str) -> int:
    rounded = int(round(float(value)))
    if abs(float(value) - rounded) > _INT_TOL:
        raise SolverFailureError(f"{what} is not an integer: {value!r}")
    return rounded


def invert_strains(
    lengths: list[int],
    strain_min: int,
    strain_max: int,
    windows: list[dict],
) -> InversionResult:
    """联合反演整数微应变。

    参数
    ----
    lengths:    按顺序排列的各缆段长度（正整数，6..12 段）。
    strain_min: 统一应变闭区间下界（整数微应变）。
    strain_max: 统一应变闭区间上界（整数微应变）。
    windows:    观测窗列表，每项为 {"start", "end", "min", "max"}，
                表示第 start..end 段（含两端）的长度加权应变和
                sum(lengths[i] * strain[i]) 必须落入闭区间 [min, max]。

    异常
    ----
    InfeasibleError:    观测窗彼此冲突，无可行解。
    SolverFailureError: 求解器异常或结果未通过精确整数复核。
    """
    n = len(lengths)
    m = n - 1  # 相邻段对数
    span = strain_max - strain_min

    def window_rows(nvar: int):
        rows, lo, hi = [], [], []
        for w in windows:
            row = [0.0] * nvar
            for i in range(w["start"], w["end"] + 1):
                row[i] = float(lengths[i])
            rows.append(row)
            lo.append(float(w["min"]))
            hi.append(float(w["max"]))
        return rows, lo, hi

    # ---------- 阶段 1：最小化 max |x_{i+1} - x_i| ----------
    # 变量：x_0..x_{n-1}, M
    nvar = n + 1
    rows, lo, hi = window_rows(nvar)
    for i in range(m):
        # x_{i+1} - x_i - M <= 0
        row = [0.0] * nvar
        row[i] = -1.0
        row[i + 1] = 1.0
        row[n] = -1.0
        rows.append(row)
        lo.append(-np.inf)
        hi.append(0.0)
        # x_i - x_{i+1} - M <= 0
        row = [0.0] * nvar
        row[i] = 1.0
        row[i + 1] = -1.0
        row[n] = -1.0
        rows.append(row)
        lo.append(-np.inf)
        hi.append(0.0)
    c = [0.0] * n + [1.0]
    var_lo = [float(strain_min)] * n + [0.0]
    var_hi = [float(strain_max)] * n + [float(span)]
    res = _run_mip(c, rows, lo, hi, var_lo, var_hi)
    if res.status == 2:
        raise InfeasibleError(
            "observation windows are mutually conflicting: no integer strain "
            "assignment satisfies every window's closed interval"
        )
    _require_optimal(res, "stage-1 (minimize max adjacent diff)")
    max_diff = _as_exact_int(res.x[n], "max adjacent diff")

    # ---------- 阶段 2：在 max_diff 约束下最小化 sum |x_{i+1} - x_i| ----------
    # 变量：x_0..x_{n-1}, M, d_0..d_{m-1}，其中 d_i >= |x_{i+1} - x_i|
    nvar = n + 1 + m
    rows, lo, hi = window_rows(nvar)
    for i in range(m):
        for sign in (1.0, -1.0):
            # sign * (x_{i+1} - x_i) - M <= 0
            row = [0.0] * nvar
            row[i] = -sign
            row[i + 1] = sign
            row[n] = -1.0
            rows.append(row)
            lo.append(-np.inf)
            hi.append(0.0)
            # sign * (x_{i+1} - x_i) - d_i <= 0
            row = [0.0] * nvar
            row[i] = -sign
            row[i + 1] = sign
            row[n + 1 + i] = -1.0
            rows.append(row)
            lo.append(-np.inf)
            hi.append(0.0)
    c = [0.0] * (n + 1) + [1.0] * m
    var_lo = [float(strain_min)] * n + [0.0] + [0.0] * m
    var_hi = [float(strain_max)] * n + [float(max_diff)] + [float(span)] * m
    res = _run_mip(c, rows, lo, hi, var_lo, var_hi)
    _require_optimal(res, "stage-2 (minimize sum of adjacent diffs)")
    sum_diff = _as_exact_int(res.fun, "sum of adjacent diffs")

    # ---------- 阶段 3：在前两级最优的前提下按字典序最小化应变序列 ----------
    # 追加约束：sum(d_i) <= sum_diff
    rows.append([0.0] * (n + 1) + [1.0] * m)
    lo.append(-np.inf)
    hi.append(float(sum_diff))
    fixed_lo = list(var_lo)
    fixed_hi = list(var_hi)
    strains: list[int] = []
    for k in range(n):
        c = [0.0] * nvar
        c[k] = 1.0
        res = _run_mip(c, rows, lo, hi, fixed_lo, fixed_hi)
        _require_optimal(res, f"stage-3 (lexicographic, segment {k})")
        value = _as_exact_int(res.x[k], f"strain[{k}]")
        strains.append(value)
        fixed_lo[k] = float(value)
        fixed_hi[k] = float(value)

    # ---------- 精确整数复核（全部比较均为整数运算） ----------
    if any(not (strain_min <= x <= strain_max) for x in strains):
        raise SolverFailureError("strain out of bounds after integer rounding")
    window_sums: list[int] = []
    for w in windows:
        s = 0
        for i in range(w["start"], w["end"] + 1):
            s += lengths[i] * strains[i]
        if not (w["min"] <= s <= w["max"]):
            raise SolverFailureError("window constraint violated after integer rounding")
        window_sums.append(s)
    diffs = [abs(strains[i + 1] - strains[i]) for i in range(m)]
    if max(diffs) != max_diff or sum(diffs) != sum_diff:
        raise SolverFailureError("smoothing metrics mismatch after integer rounding")

    return InversionResult(
        strains=strains,
        max_adjacent_diff=max_diff,
        sum_adjacent_diff=sum_diff,
        window_sums=window_sums,
    )
