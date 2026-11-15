import assert from 'node:assert/strict';
import { test } from 'node:test';
import {
  InfeasibleCurveError,
  InfeasibleTimeError,
  buildLookQuaternion,
  planSplineMotion
} from '../index.js';
import { slerp } from '../lib/vectors.js';
import { CatmullRomPath } from '../lib/spline.js';
import { OrientationTrack } from '../lib/orientation.js';

const LIMITS = {
  maxSpeed: 2,
  maxAcceleration: 4,
  maxCurvature: 10,
  sampleDt: 0.1
};

function vectorDistance(a, b) {
  return Math.hypot(a[0] - b[0], a[1] - b[1], a[2] - b[2]);
}

function assertKinematicLimits(result, maxSpeed, maxAcceleration, maxCurvature) {
  for (const sample of result.samples) {
    assert.ok(sample.speed <= maxSpeed + 1e-7);
    assert.ok(sample.accelerationMagnitude <= maxAcceleration + 2e-7);
    assert.ok(sample.curvature <= maxCurvature + 1e-7);
    assert.equal(sample.quaternion.length, 4);
    assert.ok(sample.quaternion.every(Number.isFinite));
  }
}

test('straight line has accurate arc length, stable timing and bounded motion', () => {
  const result = planSplineMotion([[0, 0, 0], [1, 0, 0], [2, 0, 0]], LIMITS);

  assert.ok(Math.abs(result.arcLength - 2) < 1e-7);
  assert.ok(result.arcLengthError < 1e-7);
  assert.deepEqual(result.samples[0].position, [0, 0, 0]);
  assert.deepEqual(result.samples.at(-1).position, [2, 0, 0]);
  assert.equal(result.samples[0].speed, 0);
  assert.equal(result.samples.at(-1).speed, 0);

  for (let i = 1; i < result.samples.length - 1; i += 1) {
    const current = result.samples[i];
    const previous = result.samples[i - 1];
    assert.ok(Math.abs(current.arcLength - previous.arcLength - current.speed * (current.time - previous.time)) < 0.02);
    assert.ok(vectorDistance(current.position, [current.arcLength, 0, 0]) < 1e-5);
  }
  assertKinematicLimits(result, LIMITS.maxSpeed, LIMITS.maxAcceleration, LIMITS.maxCurvature);
});

test('rounded reversal accepts duplicate points and decelerates through the turn', () => {
  const result = planSplineMotion([
    [0, 0, 0],
    [0, 0, 0],
    [1, 0, 0],
    [1.5, 0.5, 0],
    [1, 1, 0],
    [0, 1, 0]
  ], { ...LIMITS, maxCurvature: 8, duration: 3 });

  assert.deepEqual(result.duplicateControlPointIndexes, [1]);
  assert.ok(result.arcLength > 2.9);
  const turn = result.samples.reduce((candidate, sample) => (
    sample.curvature > candidate.curvature ? sample : candidate
  ), result.samples[0]);
  const straight = result.samples[Math.floor(result.samples.length / 4)];
  assert.ok(turn.speed < straight.speed);
  assertKinematicLimits(result, LIMITS.maxSpeed, LIMITS.maxAcceleration, 8);
});

test('circle orbit stays on radius and maintains periodic tangent continuity', () => {
  const points = [];
  for (let i = 0; i < 16; i += 1) {
    const angle = 2 * Math.PI * i / 16;
    points.push([Math.cos(angle), Math.sin(angle), 0]);
  }
  const result = planSplineMotion(points, {
    ...LIMITS,
    closed: true,
    maxCurvature: 1.25,
    sampleDt: 0.15
  });

  for (const sample of result.samples) {
    const radius = Math.hypot(sample.position[0], sample.position[1]);
    assert.ok(Math.abs(radius - 1) < 0.03);
    assert.ok(Math.abs(sample.curvature - 1) < 0.25);
  }

  const speedVariance = Math.max(...result.samples.map((sample) => sample.speed))
    - Math.min(...result.samples.map((sample) => sample.speed));
  assert.ok(speedVariance < 0.12);
  assert.ok(result.samples.at(-1).time < result.duration);
  const closedPath = new CatmullRomPath(points, { closed: true });
  assert.ok(vectorDistance(closedPath.tangentAt(0), closedPath.tangentAt(1)) < 1e-9);
  const orientation = new OrientationTrack(closedPath);
  assert.deepEqual(orientation.at(0), orientation.at(closedPath.arcLength));
  assertKinematicLimits(result, LIMITS.maxSpeed, LIMITS.maxAcceleration, 1.25);
});

test('dense nonuniform control points do not create velocity jumps', () => {
  const densePoints = [
    [0, 0, 0],
    [0.2, 0.01, 0],
    [0.4, 0.02, 0],
    [0.45, 0.03, 0],
    [1, 0.1, 0],
    [2, 0.2, 0]
  ];
  const sparsePoints = [[0, 0, 0], [1, 0.1, 0], [2, 0.2, 0]];
  const dense = planSplineMotion(densePoints, { ...LIMITS, duration: 3 });
  const sparse = planSplineMotion(sparsePoints, { ...LIMITS, duration: 3 });

  assert.ok(Math.abs(dense.arcLength - sparse.arcLength) < 0.1);
  for (let i = 1; i < dense.samples.length; i += 1) {
    const current = dense.samples[i];
    const previous = dense.samples[i - 1];
    const dt = current.time - previous.time;
    if (dt > Number.EPSILON) {
      assert.ok(Math.abs(current.speed - previous.speed) / dt < 5);
    }
  }
  assertKinematicLimits(dense, LIMITS.maxSpeed, LIMITS.maxAcceleration, LIMITS.maxCurvature);
});

test('sharp bend slows according to curvature instead of following parameter density', () => {
  const result = planSplineMotion([
    [0, 0, 0],
    [1, 0, 0],
    [1, 0.5, 0],
    [1, 1, 0],
    [0, 1, 0]
  ], { ...LIMITS, maxCurvature: 10, duration: 3 });

  const bend = result.samples.reduce((candidate, sample) => (
    sample.curvature > candidate.curvature ? sample : candidate
  ), result.samples[0]);
  assert.ok(bend.curvature > 2);
  assert.ok(bend.speed < 1.2);
  assertKinematicLimits(result, LIMITS.maxSpeed, LIMITS.maxAcceleration, 10);
});

test('closed spline rejects curvature and acceleration requests outside bounded reachability', () => {
  const points = [[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0]];
  assert.throws(() => planSplineMotion(points, {
    ...LIMITS,
    closed: true,
    maxCurvature: 0.5
  }), InfeasibleCurveError);

  const feasible = planSplineMotion(points, { ...LIMITS, closed: true, maxCurvature: 10 });
  assert.throws(() => planSplineMotion(points, {
    ...LIMITS,
    closed: true,
    maxCurvature: 10,
    duration: feasible.minimumDuration * 0.5
  }), InfeasibleTimeError);
});

test('exact 180 degree foldback is reported as a geometric cusp', () => {
  assert.throws(() => planSplineMotion([
    [0, 0, 0],
    [1, 0, 0],
    [0, 0, 0]
  ], LIMITS), InfeasibleCurveError);
});

test('invalid orientation references and zero-length target vectors are rejected', () => {
  const points = [[0, 0, 0], [1, 0, 0], [2, 0, 0]];
  assert.throws(() => planSplineMotion(points, {
    ...LIMITS,
    orientation: { targets: [[0, 0, 0], [2, 0, 0], [3, 0, 0]] }
  }), /target must not coincide/);
  assert.throws(() => planSplineMotion(points, {
    ...LIMITS,
    orientation: { up: [0, Number.NaN, 0] }
  }), /orientation\.up/);
});

test('180 degree quaternion orientation and degenerate up vector are explicit', () => {
  const halfway = slerp([1, 0, 0, 0], [0, 1, 0, 0], 0.5);
  assert.ok(Math.abs(halfway[0] - Math.SQRT1_2) < 1e-10);
  assert.ok(Math.abs(halfway[1] - Math.SQRT1_2) < 1e-10);

  const q = buildLookQuaternion([0, 1, 0], [0, 1, 0]);
  assert.equal(q.length, 4);
  assert.ok(q.every(Number.isFinite));
  assert.ok(Math.abs(Math.hypot(...q) - 1) < 1e-10);
});

test('path-level orientation survives a 180 degree reference rotation', () => {
  const result = planSplineMotion([[0, 0, 0], [1, 0, 0], [2, 0, 0]], {
    ...LIMITS,
    duration: 2,
    orientation: {
      quaternions: [
        [1, 0, 0, 0],
        [0, 1, 0, 0],
        [0, 0, 1, 0]
      ]
    }
  });

  assert.ok(result.samples.every((sample) => (
    sample.quaternion.length === 4
    && sample.quaternion.every(Number.isFinite)
    && Math.abs(Math.hypot(...sample.quaternion) - 1) < 1e-10
  )));
});

test('planning is deterministic and produces uniform time samples', () => {
  const points = [[0, 0, 0], [1, 0.1, 0], [2, 0, 0], [3, 0.1, 0]];
  const first = planSplineMotion(points, { ...LIMITS, duration: 2.5, seed: 7 });
  const second = planSplineMotion(points, { ...LIMITS, duration: 2.5, seed: 99 });

  assert.equal(first.samples.length, second.samples.length);
  for (let i = 0; i < first.samples.length; i += 1) {
    assert.deepEqual(first.samples[i].position, second.samples[i].position);
    assert.deepEqual(first.samples[i].quaternion, second.samples[i].quaternion);
    if (i > 0) {
      const dt = first.samples[i].time - first.samples[i - 1].time;
      assert.ok(Math.abs(dt - LIMITS.sampleDt) < 1e-10 || i === first.samples.length - 1);
    }
  }
});

test('control-point timestamps set total duration but do not drive parameter speed', () => {
  const result = planSplineMotion([[0, 0, 0], [1, 0, 0], [2, 0, 0]], {
    ...LIMITS,
    times: [1, 7, 13]
  });

  assert.equal(result.duration, 12);
  assert.ok(Math.abs(result.arcLength - 2) < 1e-7);
  assert.ok(result.samples.every((sample) => sample.speed <= LIMITS.maxSpeed + 1e-7));
});
