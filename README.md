# 海缆检修联合反演服务（cable-strain-inversion）

海缆检修组只能取得多个**重叠标距**（观测窗）的**累计伸长量**。本服务把所有观测窗
联合起来，反演出连续缆段的**整数微应变**序列，避免逐窗均摊后得到彼此矛盾的局部结论。

## 数学模型

- 决策变量：逐段整数微应变 `x[0..n-1]`，全部落入统一应变闭区间 `[strain_min, strain_max]`。
- 约束：每个观测窗 `w = (start, end, [emin, emax])` 要求长度加权应变和
  `sum(lengths[i] * x[i], i = start..end)` 落入闭区间 `[emin, emax]`。
- 可行解按以下顺序依次优化（分层字典序目标）：
  1. **一级平滑指标**：相邻段应变差绝对值的最大值 `max |x[i+1] - x[i]|` 最小；
  2. **二级平滑指标**：全部相邻差绝对值之和 `sum |x[i+1] - x[i]|` 最小；
  3. **应变序列字典序**最小。

每个阶段是一个纯整数线性规划（ILP），由 HiGHS 求到可证明最优（`mip_rel_gap=0`）。
求解结果在返回前用 Python 精确整数重新校验全部约束与两级指标，因此响应中的每个
回算和与平滑指标都能由提交数据直接复核。

## 快速开始（Docker）

```bash
# 构建镜像、启动 API，并在 API 健康后运行一次性 verify 验收服务
docker compose up --build --exit-code-from verify

# verify 的退出码即验收结论（0 = 全部通过）
echo $?
```

- `api`：反演 API 服务，宿主机端口通过环境变量 `API_HOST_PORT` 配置（默认 `8000`），
  内置 `/health` 健康检查（Compose healthcheck 每 5s 探测一次）。
- `verify`：一次性服务，在 `api` 健康（`service_healthy`）后依次完成
  **代码测试**（pytest）、**构建检查**（字节码编译 + ASGI 装配 + OpenAPI 生成）、
  **一组反演 API 冒烟**（健康、成功反演与复核、字段错误、不可行冲突），
  随后自行退出并以退出码报告结论。

自定义宿主机端口：

```bash
API_HOST_PORT=9000 docker compose up --build
```

## API

### `POST /api/v1/strain/invert`

请求体：

```json
{
  "segment_lengths": [120, 80, 100, 100, 90, 110, 100, 95],
  "strain_bounds": {"min": -50, "max": 50},
  "windows": [
    {"start": 0, "end": 2, "elongation": {"min": 580, "max": 630}}
  ]
}
```

| 字段 | 约束 |
| --- | --- |
| `segment_lengths` | 整数数组，6..12 段，每段长度 1..1,000,000 |
| `strain_bounds.min/max` | 整数，-1,000,000..1,000,000，且 `min <= max` |
| `windows` | 数组，8..20 个观测窗 |
| `windows[i].start/end` | 整数，`0 <= start <= end < 段数`（连续起止段） |
| `windows[i].elongation.min/max` | 整数，±12,000,000,000,000 内，且 `min <= max` |

成功（`200`）：返回逐段应变、两级平滑指标、逐窗回算证据，并回显提交数据便于复核。

```json
{
  "status": "optimal",
  "segment_lengths": [120, 80, 100, 100, 90, 110, 100, 95],
  "strain_bounds": {"min": -50, "max": 50},
  "strains": [3, 3, 2, 2, 1, 1, 0, 0],
  "objectives": {"max_adjacent_diff": 1, "sum_adjacent_diff": 3},
  "windows": [
    {"index": 0, "start": 0, "end": 2,
     "elongation_min": 580, "elongation_max": 630,
     "weighted_sum": 600, "within_bounds": true}
  ],
  "checks": {"strains_within_bounds": true, "windows_within_bounds": true}
}
```

复核方法（全部精确整数运算）：

- 逐窗回算：`weighted_sum == sum(segment_lengths[i] * strains[i], i = start..end)`，
  且落在 `[elongation_min, elongation_max]` 内；
- 一级指标：`max_adjacent_diff == max(|strains[i+1] - strains[i]|)`；
- 二级指标：`sum_adjacent_diff == sum(|strains[i+1] - strains[i]|)`。

输入越界或窗口非法（`422`）——一次性报告全部字段错误：

```json
{"status": "invalid",
 "errors": [{"field": "windows[0].end", "message": "must be within 0..7"}]}
```

观测彼此冲突（`409`）：

```json
{"status": "infeasible",
 "detail": "observation windows are mutually conflicting: ..."}
```

### `GET /health`

返回 `{"status": "ok"}`，供 Docker 健康检查使用。

## 本地开发

```bash
python3 -m venv .venv
.venv/bin/pip install -r requirements.txt
.venv/bin/python -m pytest -q                 # 单元测试（含暴力枚举交叉验证）
.venv/bin/uvicorn app.main:app --port 8000    # 启动服务
API_BASE_URL=http://127.0.0.1:8000 .venv/bin/python verify.py   # 本地验收
```

## 代码结构

```
app/
  main.py        # FastAPI 应用：/health 与 /api/v1/strain/invert
  solver.py      # 三阶段 ILP 联合反演（HiGHS），精确整数复核
  validation.py  # 请求体校验，汇总全部字段错误
tests/
  test_solver.py # 与暴力枚举参考实现交叉验证 + 边界用例
  test_api.py    # API 行为测试（成功 / 字段错误 / 不可行）
verify.py        # 一次性验收：代码测试 + 构建检查 + API 冒烟
Dockerfile       # 单一镜像，同时用于 api 与 verify 服务
docker-compose.yml
```
