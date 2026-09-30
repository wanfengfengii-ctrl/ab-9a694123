"""一次性验收服务（docker compose 服务名：verify）。

在 API 依赖健康之后依次完成：
  1. 代码测试（pytest 单元测试）
  2. 构建检查（字节码编译 + ASGI 应用装配 + OpenAPI 生成）
  3. 一组反演 API 冒烟（健康、成功反演及复核、字段错误、不可行冲突）

全部通过则以退出码 0 自行退出，任何一步失败以退出码 1 退出。
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
import time
import urllib.error
import urllib.request

REPO_ROOT = os.path.dirname(os.path.abspath(__file__))
API_BASE_URL = os.environ.get("API_BASE_URL", "http://localhost:8000").rstrip("/")
HEALTH_TIMEOUT_SECONDS = float(os.environ.get("VERIFY_HEALTH_TIMEOUT", "90"))


# ---------------------------------------------------------------- HTTP 工具

def http(method: str, path: str, payload=None) -> tuple[int, dict]:
    data = None if payload is None else json.dumps(payload).encode()
    req = urllib.request.Request(
        API_BASE_URL + path,
        data=data,
        method=method,
        headers={"Content-Type": "application/json"},
    )
    try:
        with urllib.request.urlopen(req, timeout=30) as resp:
            return resp.status, json.loads(resp.read().decode())
    except urllib.error.HTTPError as exc:
        return exc.code, json.loads(exc.read().decode())


# ---------------------------------------------------------------- 测试载荷

def feasible_payload() -> dict:
    lengths = [120, 80, 100, 100, 90, 110, 100, 95]
    target = [3, 3, 2, 2, 1, 1, 0, 0]
    spans = [(0, 2), (1, 3), (2, 4), (3, 5), (4, 6), (5, 7), (0, 3), (4, 7), (0, 7), (2, 5)]
    windows = []
    for a, b in spans:
        s = sum(lengths[i] * target[i] for i in range(a, b + 1))
        windows.append({"start": a, "end": b, "elongation": {"min": s - 25, "max": s + 25}})
    return {
        "segment_lengths": lengths,
        "strain_bounds": {"min": -50, "max": 50},
        "windows": windows,
    }


def infeasible_payload() -> dict:
    windows = [
        {"start": 0, "end": 2, "elongation": {"min": 400, "max": 500}},
        {"start": 0, "end": 2, "elongation": {"min": -500, "max": -400}},
    ] + [{"start": i, "end": i, "elongation": {"min": -1000, "max": 1000}} for i in range(6)]
    return {
        "segment_lengths": [100] * 6,
        "strain_bounds": {"min": -5, "max": 5},
        "windows": windows,
    }


# ---------------------------------------------------------------- 各验收步骤

def wait_for_api() -> None:
    deadline = time.monotonic() + HEALTH_TIMEOUT_SECONDS
    last = None
    while time.monotonic() < deadline:
        try:
            status, body = http("GET", "/health")
            if status == 200 and body.get("status") == "ok":
                print(f"API healthy at {API_BASE_URL}")
                return
            last = f"status={status} body={body}"
        except Exception as exc:  # noqa: BLE001 - 等待阶段容忍一切连接错误
            last = str(exc)
        time.sleep(2)
    raise RuntimeError(f"API did not become healthy within {HEALTH_TIMEOUT_SECONDS}s: {last}")


def run_unit_tests() -> None:
    proc = subprocess.run(
        [sys.executable, "-m", "pytest", "-q", "tests"],
        cwd=REPO_ROOT,
        capture_output=True,
        text=True,
    )
    sys.stdout.write(proc.stdout)
    if proc.returncode != 0:
        sys.stdout.write(proc.stderr)
        raise RuntimeError(f"pytest exited with code {proc.returncode}")


def run_build_check() -> None:
    proc = subprocess.run(
        [sys.executable, "-m", "compileall", "-q", "app", "tests", "verify.py"],
        cwd=REPO_ROOT,
        capture_output=True,
        text=True,
    )
    if proc.returncode != 0:
        raise RuntimeError(f"byte-compile failed:\n{proc.stderr}")
    proc = subprocess.run(
        [
            sys.executable,
            "-c",
            "from app.main import app; "
            "schema = app.openapi(); "
            "print('openapi paths:', ', '.join(sorted(schema['paths'])))",
        ],
        cwd=REPO_ROOT,
        capture_output=True,
        text=True,
    )
    if proc.returncode != 0:
        raise RuntimeError(f"ASGI app assembly failed:\n{proc.stderr}")
    print(proc.stdout.strip())


def smoke_health() -> None:
    status, body = http("GET", "/health")
    assert status == 200, f"expected 200, got {status}: {body}"
    assert body.get("status") == "ok", body


def smoke_invert_optimal() -> None:
    payload = feasible_payload()
    status, body = http("POST", "/api/v1/strain/invert", payload)
    assert status == 200, f"expected 200, got {status}: {body}"
    assert body["status"] == "optimal", body

    lengths = payload["segment_lengths"]
    strains = body["strains"]
    assert len(strains) == len(lengths), "应变序列长度与缆段数不一致"
    lo = payload["strain_bounds"]["min"]
    hi = payload["strain_bounds"]["max"]
    assert all(lo <= x <= hi for x in strains), "应变越出统一闭区间"

    # 逐窗回算证据：每个 weighted_sum 必须能由提交数据直接复核
    assert len(body["windows"]) == len(payload["windows"])
    for w_in, w_out in zip(payload["windows"], body["windows"]):
        s = sum(lengths[i] * strains[i] for i in range(w_in["start"], w_in["end"] + 1))
        assert w_out["weighted_sum"] == s, f"窗口 {w_out['index']} 回算和与提交数据不符"
        assert w_in["elongation"]["min"] <= s <= w_in["elongation"]["max"], f"窗口 {w_out['index']} 越界"

    # 两级平滑指标复核
    diffs = [abs(strains[i + 1] - strains[i]) for i in range(len(strains) - 1)]
    assert body["objectives"]["max_adjacent_diff"] == max(diffs), "一级平滑指标复核失败"
    assert body["objectives"]["sum_adjacent_diff"] == sum(diffs), "二级平滑指标复核失败"
    print("strains:", strains, "objectives:", body["objectives"])


def smoke_invert_invalid() -> None:
    # 段长值越界 + 窗口下标越界：两个字段错误应一次性同时报告
    payload = feasible_payload()
    payload["segment_lengths"][0] = 0          # 段长越界
    payload["windows"][0]["end"] = 99          # 窗口下标越界
    status, body = http("POST", "/api/v1/strain/invert", payload)
    assert status == 422, f"expected 422, got {status}: {body}"
    assert body["status"] == "invalid", body
    fields = {e["field"] for e in body["errors"]}
    assert "segment_lengths[0]" in fields, fields
    assert "windows[0].end" in fields, fields
    assert all(e.get("message") for e in body["errors"]), "字段错误缺少 message"

    # 段数越界也应得到明确字段错误
    payload = feasible_payload()
    payload["segment_lengths"] = payload["segment_lengths"][:5]
    status, body = http("POST", "/api/v1/strain/invert", payload)
    assert status == 422, f"expected 422, got {status}: {body}"
    assert body["status"] == "invalid", body
    assert any(e["field"] == "segment_lengths" for e in body["errors"]), body["errors"]


def smoke_invert_infeasible() -> None:
    status, body = http("POST", "/api/v1/strain/invert", infeasible_payload())
    assert status == 409, f"expected 409, got {status}: {body}"
    assert body["status"] == "infeasible", body


# ---------------------------------------------------------------- 主流程

def main() -> int:
    steps = [
        ("等待 API 健康", wait_for_api),
        ("代码测试 (pytest)", run_unit_tests),
        ("构建检查", run_build_check),
        ("冒烟: 健康检查", smoke_health),
        ("冒烟: 成功反演与复核", smoke_invert_optimal),
        ("冒烟: 字段错误", smoke_invert_invalid),
        ("冒烟: 不可行冲突", smoke_invert_infeasible),
    ]
    failures = []
    for name, fn in steps:
        print(f"\n=== {name} ===", flush=True)
        try:
            fn()
        except Exception as exc:  # noqa: BLE001 - 任何失败都记录并继续后续步骤
            failures.append(name)
            print(f"[FAIL] {name}: {exc}", flush=True)
        else:
            print(f"[PASS] {name}", flush=True)

    print("\n=== 验收汇总 ===", flush=True)
    if failures:
        print("失败步骤: " + ", ".join(failures), flush=True)
        return 1
    print("全部通过", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
