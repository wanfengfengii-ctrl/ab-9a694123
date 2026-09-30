"""求解器单元测试：与暴力枚举（按定义显然正确的参考实现）交叉验证。"""

from __future__ import annotations

import itertools
import random

import pytest

from app.solver import InfeasibleError, invert_strains


def brute_force_best(lengths, strain_min, strain_max, windows):
    """枚举全部整数应变序列，返回 (最优序列, max_diff, sum_diff) 或 None。

    目标顺序与求解器一致：先最小化相邻差最大值，再最小化相邻差之和，
    最后取字典序最小。
    """
    n = len(lengths)
    best_key = None
    best_xs = None
    for xs in itertools.product(range(strain_min, strain_max + 1), repeat=n):
        feasible = True
        for w in windows:
            s = sum(lengths[i] * xs[i] for i in range(w["start"], w["end"] + 1))
            if not (w["min"] <= s <= w["max"]):
                feasible = False
                break
        if not feasible:
            continue
        diffs = [abs(xs[i + 1] - xs[i]) for i in range(n - 1)]
        key = (max(diffs), sum(diffs), xs)
        if best_key is None or key < best_key:
            best_key = key
            best_xs = xs
    if best_key is None:
        return None
    return best_xs, best_key[0], best_key[1]


def make_windows_around(lengths, target, spans, slack):
    windows = []
    for a, b in spans:
        s = sum(lengths[i] * target[i] for i in range(a, b + 1))
        windows.append({"start": a, "end": b, "min": s - slack, "max": s + slack})
    return windows


def assert_matches_reference(lengths, smin, smax, windows):
    expected = brute_force_best(lengths, smin, smax, windows)
    if expected is None:
        with pytest.raises(InfeasibleError):
            invert_strains(lengths, smin, smax, windows)
        return
    result = invert_strains(lengths, smin, smax, windows)
    assert tuple(result.strains) == expected[0]
    assert result.max_adjacent_diff == expected[1]
    assert result.sum_adjacent_diff == expected[2]
    # 逐窗回算和必须能由提交数据直接复核
    for w, s in zip(windows, result.window_sums):
        assert s == sum(lengths[i] * result.strains[i] for i in range(w["start"], w["end"] + 1))
        assert w["min"] <= s <= w["max"]


def test_deterministic_case():
    lengths = [3, 1, 4, 1, 5, 9]
    target = [2, -1, 0, 1, 1, -2]
    spans = [(0, 2), (1, 3), (2, 4), (3, 5), (0, 3), (2, 5), (0, 5), (1, 4)]
    windows = make_windows_around(lengths, target, spans, slack=2)
    assert_matches_reference(lengths, -4, 4, windows)


def test_constant_strain_is_optimal_when_feasible():
    # 所有窗口都允许常数应变时，两级平滑指标都应为 0
    lengths = [10, 20, 30, 40, 50, 60]
    target = [3] * 6
    spans = [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (0, 2), (2, 4), (0, 5)]
    windows = make_windows_around(lengths, target, spans, slack=10)
    result = invert_strains(lengths, -10, 10, windows)
    assert result.max_adjacent_diff == 0
    assert result.sum_adjacent_diff == 0
    assert len(set(result.strains)) == 1


def test_lexicographic_tie_break():
    # 两个块之间必须过渡时，字典序最小的解应把较小的值尽量放前面
    lengths = [1] * 6
    windows = [
        {"start": 0, "end": 2, "min": 3, "max": 3},   # x0+x1+x2 == 3
        {"start": 3, "end": 5, "min": 0, "max": 0},   # x3+x4+x5 == 0
        {"start": 0, "end": 5, "min": 3, "max": 3},
        {"start": 0, "end": 0, "min": -10, "max": 10},
        {"start": 1, "end": 1, "min": -10, "max": 10},
        {"start": 2, "end": 2, "min": -10, "max": 10},
        {"start": 3, "end": 4, "min": -10, "max": 10},
        {"start": 4, "end": 5, "min": -10, "max": 10},
    ]
    assert_matches_reference(lengths, -10, 10, windows)
    result = invert_strains(lengths, -10, 10, windows)
    # 参考解：相邻差最大值 1、相邻差之和 1（两个等值块之间恰好过渡一次）
    assert result.max_adjacent_diff == 1
    assert result.sum_adjacent_diff == 1
    assert result.strains == [1, 1, 1, 0, 0, 0]


def test_infeasible_conflicting_windows():
    lengths = [10] * 6
    windows = [
        {"start": 0, "end": 2, "min": 400, "max": 500},
        {"start": 0, "end": 2, "min": -500, "max": -400},
    ] + [{"start": i, "end": i, "min": -1000, "max": 1000} for i in range(6)]
    with pytest.raises(InfeasibleError):
        invert_strains(lengths, -5, 5, windows)


def test_random_cases_against_brute_force():
    rng = random.Random(20260930)
    for case in range(12):
        n = 6
        lengths = [rng.randint(1, 4) for _ in range(n)]
        smin, smax = -3, 3
        spans = set()
        while len(spans) < 8:
            a = rng.randrange(n)
            b = rng.randrange(n)
            if a > b:
                a, b = b, a
            spans.add((a, b))
        windows = []
        for a, b in sorted(spans):
            # 围绕随机中心取小区间：既可能可行也可能不可行，两种路径都要走到
            lo_full = sum(lengths[i] * smin for i in range(a, b + 1))
            hi_full = sum(lengths[i] * smax for i in range(a, b + 1))
            center = rng.randint(lo_full, hi_full)
            half = rng.randint(0, 4)
            windows.append({"start": a, "end": b, "min": center - half, "max": center + half})
        assert_matches_reference(lengths, smin, smax, windows)


def test_large_case_consistency():
    # 12 段、20 个窗口的最大规模：不做暴力枚举，只验证结果自洽且可复核
    rng = random.Random(7)
    n = 12
    lengths = [rng.randint(50, 150) for _ in range(n)]
    target = [rng.randint(-20, 20) for _ in range(n)]
    spans = set()
    while len(spans) < 20:
        a = rng.randrange(n)
        b = rng.randrange(n)
        if a > b:
            a, b = b, a
        spans.add((a, b))
    windows = make_windows_around(lengths, target, sorted(spans), slack=30)
    result = invert_strains(lengths, -30, 30, windows)
    assert len(result.strains) == n
    assert len(result.window_sums) == 20
    assert all(-30 <= x <= 30 for x in result.strains)
    diffs = [abs(result.strains[i + 1] - result.strains[i]) for i in range(n - 1)]
    assert max(diffs) == result.max_adjacent_diff
    assert sum(diffs) == result.sum_adjacent_diff
    for w, s in zip(windows, result.window_sums):
        assert s == sum(lengths[i] * result.strains[i] for i in range(w["start"], w["end"] + 1))
        assert w["min"] <= s <= w["max"]
