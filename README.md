# spline-motion-planner

Node.js 镜头样条运动规划库。它把非均匀位置控制点构造成向心 Catmull–Rom 样条，再通过弧长重参数化、曲率限速和加速度扫描生成稳定的等间隔时间采样。控制点在原始参数域中的疏密不会直接决定相机速度。

## 安装与测试

```bash
node --test
```

项目只使用 Node.js 内置模块，无需安装第三方依赖，ESM 入口为 `index.js`。

## 接口

```js
import { planSplineMotion } from './index.js';

const result = planSplineMotion(points, {
  duration: 4,
  sampleDt: 0.1,
  maxSpeed: 2,
  maxAcceleration: 4,
  maxCurvature: 5,
  closed: false,
  orientation: {
    up: [0, 1, 0]
  }
});
```

`points` 是 `[x, y, z][]`。时间参数可二选一：

- `duration`：总时长，单位秒。
- `times`：与原始控制点一一对应的严格递增时间戳，单位秒；库使用其总跨度作为目标时长。

`orientation` 支持：

- `{ up }`：沿路径切线观察，并用参考上方向抑制横滚。
- `{ targets, up }`：每个控制点给出观察目标点，按弧长做四元数球面插值。
- `{ quaternions }`：直接提供每个去重后控制点的 `[w, x, y, z]` 单位四元数。

## 单位

所有输入和输出采用一致的 SI 风格单位：

- 位置、弧长：米。
- 时间：秒。
- 速度：米/秒。
- 线性加速度：米/秒²，包含切向加速度和向心加速度的合矢量。
- 曲率：1/米，转弯半径为曲率的倒数。
- 四元数：`[w, x, y, z]`，无单位。
- 角度只在内部轴角计算中使用弧度。

## 数值方法

1. 使用向心 Catmull–Rom 样条生成 C1 位置路径；开放路径通过镜像幽灵点避免端点零切线，闭合路径按索引环绕并保持首尾切线连续。
2. 弧长积分使用自适应 Gauss–Legendre 求积：同时比较三点和五点规则，误差不满足相对容差时二分区间，深度上限为 16。
3. 建路径时预计算 `(u,s)` 弧长表，查询时用二分定位区间并做单调线性插值反解 `u(s)`，再以 `s` 作为运动规划自变量，因此参数速度不再受控制点密度影响。
4. 曲率由 `κ = |r'(u) × r''(u)| / |r'(u)|³` 计算。
5. 速度上限同时包含 `v ≤ maxSpeed` 和向心约束 `v²κ ≤ maxAcceleration`。
6. 在弧长网格上执行前向和后向扫描，使用二分求解切向加速度与向心加速度的合量约束；闭合轨道环绕扫描，直到接缝速度收敛。
7. 朝向使用最短路径四元数 slerp，显式处理接近 180 度的相对旋转；当观察方向与上方向平行时自动选择退化最小的备用世界轴。

## 边界情况

- 重复或零长度相邻控制点会被移除，原索引记录在 `duplicateControlPointIndexes`。
- 闭合路径末尾若重复首点，会自动删除重复闭合点。
- 精确 180 度折返会产生零切线和曲率奇点，抛出 `InfeasibleCurveError`；应加入中间控制点形成圆角。
- 几何曲率超过 `maxCurvature` 时抛出 `InfeasibleCurveError`，并在 `details` 中给出弧长和曲率。
- 有界速度、加速度和曲率下若请求时长小于物理可行的最小时间，抛出 `InfeasibleTimeError`。
- 所有采样结果按固定输入顺序和确定性算法生成，不使用随机数。

## 输出

`result.samples[]` 包含：

- `time`、`index`
- `arcLength`
- `position`、`tangent`
- `velocity`、`speed`
- `acceleration`、`accelerationMagnitude`
- `curvature`
- `quaternion`

顶层还返回 `duration`、`arcLength`、`arcLengthError`、`minimumDuration`、`speedScale` 和积分诊断信息。开放路径包含终点；闭合路径不重复周期起点，最后一个采样点位于一个周期之前。

## 测试覆盖

`test/planner.test.js` 使用内置 `node:test`，覆盖：

- 直线弧长误差、端点和速度约束。
- 含重复点的圆角折返。
- 圆轨半径、曲率、周期切线和速度稳定性。
- 密集非均匀控制点下的速度连续性。
- 急弯曲率限速和加速度约束。
- 循环闭合的不可达曲率与不可达时间。
- 精确 180 度折返奇点。
- 180 度四元数、退化上方向和路径级 180 度朝向参考。
- 非法朝向参考和零长度观察向量。
- 确定性与等间隔时间采样。
- 控制点时间戳总时长与弧长采样解耦。
