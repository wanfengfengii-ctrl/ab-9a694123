"""API 层测试：健康检查、成功反演、字段错误、不可行冲突。"""

from __future__ import annotations

from fastapi.testclient import TestClient

from app.main import app

client = TestClient(app)


def feasible_payload():
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


def infeasible_payload():
    windows = [
        {"start": 0, "end": 2, "elongation": {"min": 400, "max": 500}},
        {"start": 0, "end": 2, "elongation": {"min": -500, "max": -400}},
    ] + [
        {"start": i, "end": i, "elongation": {"min": -1000, "max": 1000}} for i in range(6)
    ]
    return {
        "segment_lengths": [100] * 6,
        "strain_bounds": {"min": -5, "max": 5},
        "windows": windows,
    }


def test_health():
    resp = client.get("/health")
    assert resp.status_code == 200
    assert resp.json() == {"status": "ok"}


def test_invert_optimal_and_self_verifiable():
    payload = feasible_payload()
    resp = client.post("/api/v1/strain/invert", json=payload)
    assert resp.status_code == 200, resp.text
    body = resp.json()
    assert body["status"] == "optimal"

    lengths = payload["segment_lengths"]
    strains = body["strains"]
    assert len(strains) == len(lengths)
    lo = payload["strain_bounds"]["min"]
    hi = payload["strain_bounds"]["max"]
    assert all(lo <= x <= hi for x in strains)

    # 逐窗回算：用提交数据独立重算，必须与响应一致且落入闭区间
    assert len(body["windows"]) == len(payload["windows"])
    for w_in, w_out in zip(payload["windows"], body["windows"]):
        s = sum(lengths[i] * strains[i] for i in range(w_in["start"], w_in["end"] + 1))
        assert w_out["weighted_sum"] == s
        assert w_in["elongation"]["min"] <= s <= w_in["elongation"]["max"]
        assert w_out["within_bounds"] is True

    # 两级平滑指标回算
    diffs = [abs(strains[i + 1] - strains[i]) for i in range(len(strains) - 1)]
    assert body["objectives"]["max_adjacent_diff"] == max(diffs)
    assert body["objectives"]["sum_adjacent_diff"] == sum(diffs)
    assert body["checks"]["strains_within_bounds"] is True
    assert body["checks"]["windows_within_bounds"] is True


def test_invert_invalid_segment_count():
    payload = feasible_payload()
    payload["segment_lengths"] = payload["segment_lengths"][:5]  # 只有 5 段
    resp = client.post("/api/v1/strain/invert", json=payload)
    assert resp.status_code == 422
    body = resp.json()
    assert body["status"] == "invalid"
    assert any(e["field"] == "segment_lengths" for e in body["errors"])
    assert all("field" in e and "message" in e for e in body["errors"])


def test_invert_invalid_window_index_and_type():
    payload = feasible_payload()
    payload["windows"][0]["end"] = 99            # 越界
    payload["windows"][1]["start"] = 2.5         # 非整数
    payload["windows"][2]["elongation"] = {"min": 10, "max": -10}  # min > max
    resp = client.post("/api/v1/strain/invert", json=payload)
    assert resp.status_code == 422
    fields = {e["field"] for e in resp.json()["errors"]}
    assert "windows[0].end" in fields
    assert "windows[1].start" in fields
    assert "windows[2].elongation" in fields


def test_invert_reports_multiple_field_errors_together():
    # 段长值越界与窗口下标越界应一次性同时报告
    payload = feasible_payload()
    payload["segment_lengths"][0] = 0
    payload["windows"][0]["end"] = 99
    resp = client.post("/api/v1/strain/invert", json=payload)
    assert resp.status_code == 422
    fields = {e["field"] for e in resp.json()["errors"]}
    assert "segment_lengths[0]" in fields
    assert "windows[0].end" in fields


def test_invert_invalid_strain_bounds():
    payload = feasible_payload()
    payload["strain_bounds"] = {"min": 10, "max": -10}
    resp = client.post("/api/v1/strain/invert", json=payload)
    assert resp.status_code == 422
    assert any(e["field"] == "strain_bounds" for e in resp.json()["errors"])


def test_invert_window_count_out_of_range():
    payload = feasible_payload()
    payload["windows"] = payload["windows"][:7]  # 只有 7 个窗口
    resp = client.post("/api/v1/strain/invert", json=payload)
    assert resp.status_code == 422
    assert any(e["field"] == "windows" for e in resp.json()["errors"])


def test_invert_infeasible_returns_409():
    resp = client.post("/api/v1/strain/invert", json=infeasible_payload())
    assert resp.status_code == 409
    body = resp.json()
    assert body["status"] == "infeasible"
    assert "detail" in body


def test_invert_rejects_non_object_body():
    resp = client.post("/api/v1/strain/invert", content=b"[1, 2, 3]", headers={"Content-Type": "application/json"})
    assert resp.status_code == 422
    assert resp.json()["status"] == "invalid"
