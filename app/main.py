"""海缆检修联合反演 API。

POST /api/v1/strain/invert
    输入按顺序排列的缆段长度、统一应变闭区间、一组观测窗
    （连续起止段 + 累计伸长量闭区间），输出逐段整数微应变、
    逐窗回算证据以及两级平滑指标。

GET /health
    健康检查，供 Docker Compose healthcheck 使用。
"""

from __future__ import annotations

from fastapi import FastAPI, Request
from fastapi.responses import JSONResponse

from app.solver import InfeasibleError, SolverFailureError, invert_strains
from app.validation import (
    MAX_ELONGATION_ABS,
    MAX_SEGMENT_LENGTH,
    MAX_SEGMENTS,
    MAX_STRAIN_ABS,
    MAX_WINDOWS,
    MIN_SEGMENTS,
    MIN_WINDOWS,
    validate_payload,
)

app = FastAPI(
    title="cable-strain-inversion",
    version="1.0.0",
    description="联合反演连续缆段整数微应变：重叠观测窗累计伸长量 -> 逐段应变",
)


@app.get("/health")
async def health():
    return {"status": "ok"}


@app.get("/")
async def index():
    return {
        "service": "cable-strain-inversion",
        "endpoints": {"invert": "/api/v1/strain/invert", "health": "/health"},
        "limits": {
            "segments": [MIN_SEGMENTS, MAX_SEGMENTS],
            "windows": [MIN_WINDOWS, MAX_WINDOWS],
            "segment_length": [1, MAX_SEGMENT_LENGTH],
            "strain": [-MAX_STRAIN_ABS, MAX_STRAIN_ABS],
            "elongation": [-MAX_ELONGATION_ABS, MAX_ELONGATION_ABS],
        },
    }


@app.post("/api/v1/strain/invert")
async def invert(request: Request):
    try:
        payload = await request.json()
    except Exception:
        return JSONResponse(
            status_code=422,
            content={
                "status": "invalid",
                "errors": [{"field": "(body)", "message": "request body must be valid JSON"}],
            },
        )

    data, errors = validate_payload(payload)
    if errors:
        return JSONResponse(status_code=422, content={"status": "invalid", "errors": errors})

    try:
        result = invert_strains(
            data["segment_lengths"],
            data["strain_min"],
            data["strain_max"],
            data["windows"],
        )
    except InfeasibleError as exc:
        return JSONResponse(status_code=409, content={"status": "infeasible", "detail": str(exc)})
    except SolverFailureError as exc:
        return JSONResponse(status_code=500, content={"status": "error", "detail": str(exc)})

    windows_out = []
    for idx, (w, s) in enumerate(zip(data["windows"], result.window_sums)):
        windows_out.append(
            {
                "index": idx,
                "start": w["start"],
                "end": w["end"],
                "elongation_min": w["min"],
                "elongation_max": w["max"],
                "weighted_sum": s,
                "within_bounds": w["min"] <= s <= w["max"],
            }
        )

    return {
        "status": "optimal",
        # 回显提交数据，便于调用方直接复核每个回算和与两级平滑指标。
        "segment_lengths": data["segment_lengths"],
        "strain_bounds": {"min": data["strain_min"], "max": data["strain_max"]},
        "strains": result.strains,
        "objectives": {
            "max_adjacent_diff": result.max_adjacent_diff,
            "sum_adjacent_diff": result.sum_adjacent_diff,
        },
        "windows": windows_out,
        "checks": {
            "strains_within_bounds": all(
                data["strain_min"] <= x <= data["strain_max"] for x in result.strains
            ),
            "windows_within_bounds": all(w["within_bounds"] for w in windows_out),
        },
    }
