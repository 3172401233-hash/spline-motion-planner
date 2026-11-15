import { ValidationError } from './errors.js';

export const add = (a, b) => [a[0] + b[0], a[1] + b[1], a[2] + b[2]];
export const subtract = (a, b) => [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
export const scale = (a, scalar) => [a[0] * scalar, a[1] * scalar, a[2] * scalar];
const negate4 = (a) => [-a[0], -a[1], -a[2], -a[3]];
export const dot = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
export const cross = (a, b) => [
  a[1] * b[2] - a[2] * b[1],
  a[2] * b[0] - a[0] * b[2],
  a[0] * b[1] - a[1] * b[0]
];

export function norm(a) {
  return Math.hypot(a[0], a[1], a[2]);
}

export function normalize(a) {
  const length = norm(a);
  if (length <= Number.EPSILON) {
    throw new ValidationError('Cannot normalize a zero-length vector.');
  }
  return scale(a, 1 / length);
}

export function clamp(value, min, max) {
  return Math.min(max, Math.max(min, value));
}

export function isFiniteNumber(value) {
  return typeof value === 'number' && Number.isFinite(value);
}

export function isFiniteVector(value) {
  return Array.isArray(value) && value.length === 3 && value.every(isFiniteNumber);
}

export function distance(a, b) {
  return norm(subtract(a, b));
}

export function normalizeQuaternion(q) {
  if (!Array.isArray(q) || q.length !== 4 || !q.every(isFiniteNumber)) {
    throw new ValidationError('Quaternions must be finite arrays [w, x, y, z].');
  }
  const length = Math.sqrt(q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3]);
  if (length <= Number.EPSILON) {
    throw new ValidationError('Cannot normalize a zero-length quaternion.');
  }
  return q.map((component) => component / length);
}

export function quaternionDot(a, b) {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2] + a[3] * b[3];
}

export function slerp(start, end, t) {
  const q0 = normalizeQuaternion(start);
  let q1 = normalizeQuaternion(end);
  let cosine = quaternionDot(q0, q1);

  if (cosine < 0) {
    q1 = negate4(q1);
    cosine = -cosine;
  }

  if (cosine > 0.9995) {
    const result = q0.map((value, index) => value + (q1[index] - value) * t);
    return normalizeQuaternion(result);
  }

  const theta = Math.acos(clamp(cosine, -1, 1));
  const sine = Math.sin(theta);
  const weight0 = Math.sin((1 - t) * theta) / sine;
  const weight1 = Math.sin(t * theta) / sine;
  const interpolated = q0.map((value, index) => value * weight0 + q1[index] * weight1);
  if (theta > Math.PI / 2 - 1e-4 && interpolated.some((value) => !Number.isFinite(value))) {
    const relativeVector = [
      q1[1] * q0[0] - q1[0] * q0[1] - q1[3] * q0[2] + q1[2] * q0[3],
      q1[2] * q0[0] + q1[3] * q0[1] - q1[0] * q0[2] - q1[1] * q0[3],
      q1[3] * q0[0] - q1[2] * q0[1] + q1[1] * q0[2] - q1[0] * q0[3]
    ];
    const axisLength = Math.hypot(relativeVector[0], relativeVector[1], relativeVector[2]);
    if (axisLength <= Number.EPSILON) {
      return t < 0.5 ? q0.slice() : q1.slice();
    }
    const axis = scale(relativeVector, 1 / axisLength);
    const halfAngle = theta * t;
    const halfSine = Math.sin(halfAngle);
    const partial = [
      Math.cos(halfAngle),
      axis[0] * halfSine,
      axis[1] * halfSine,
      axis[2] * halfSine
    ];
    return normalizeQuaternion([
      q0[0] * partial[0] - q0[1] * partial[1] - q0[2] * partial[2] - q0[3] * partial[3],
      q0[0] * partial[1] + q0[1] * partial[0] + q0[2] * partial[3] - q0[3] * partial[2],
      q0[0] * partial[2] - q0[1] * partial[3] + q0[2] * partial[0] + q0[3] * partial[1],
      q0[0] * partial[3] + q0[1] * partial[2] - q0[2] * partial[1] + q0[3] * partial[0]
    ]);
  }
  return normalizeQuaternion(interpolated);
}

function quaternionFromMatrix(m00, m01, m02, m10, m11, m12, m20, m21, m22) {
  const trace = m00 + m11 + m22;
  let w;
  let x;
  let y;
  let z;

  if (trace > 0) {
    const s = Math.sqrt(trace + 1) * 2;
    w = 0.25 * s;
    x = (m12 - m21) / s;
    y = (m20 - m02) / s;
    z = (m01 - m10) / s;
  } else if (m00 > m11 && m00 > m22) {
    const s = Math.sqrt(1 + m00 - m11 - m22) * 2;
    w = (m12 - m21) / s;
    x = 0.25 * s;
    y = (m10 + m01) / s;
    z = (m20 + m02) / s;
  } else if (m11 > m22) {
    const s = Math.sqrt(1 + m11 - m00 - m22) * 2;
    w = (m20 - m02) / s;
    x = (m10 + m01) / s;
    y = 0.25 * s;
    z = (m21 + m12) / s;
  } else {
    const s = Math.sqrt(1 + m22 - m00 - m11) * 2;
    w = (m01 - m10) / s;
    x = (m20 + m02) / s;
    y = (m21 + m12) / s;
    z = 0.25 * s;
  }

  return normalizeQuaternion([w, x, y, z]);
}

export function quaternionFromAxes(right, up, forward) {
  return quaternionFromMatrix(
    right[0], up[0], forward[0],
    right[1], up[1], forward[1],
    right[2], up[2], forward[2]
  );
}

export function rotateVector(q, v) {
  const [w, x, y, z] = q;
  const qv = [x, y, z];
  const uv = cross(qv, v);
  const uuv = cross(qv, uv);
  const term = add(scale(uv, 2 * w), scale(uuv, 2));
  return add(v, term);
}
