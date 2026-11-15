import { ValidationError } from './errors.js';
import {
  cross,
  dot,
  isFiniteVector,
  normalize,
  normalizeQuaternion,
  quaternionFromAxes,
  slerp,
  subtract
} from './vectors.js';

const WORLD_AXES = [
  [1, 0, 0],
  [0, 1, 0],
  [0, 0, 1]
];

function fallbackUp(forward) {
  let best = WORLD_AXES[0];
  let bestAlignment = Math.abs(dot(forward, best));
  for (const axis of WORLD_AXES.slice(1)) {
    const alignment = Math.abs(dot(forward, axis));
    if (alignment < bestAlignment) {
      best = axis;
      bestAlignment = alignment;
    }
  }
  return best;
}

export function buildLookQuaternion(forward, referenceUp = [0, 1, 0]) {
  const f = normalize(forward);
  let up = normalize(referenceUp);

  if (Math.abs(dot(f, up)) > 1 - 1e-9) {
    up = fallbackUp(f);
  }

  const right = normalize(cross(up, f));
  const correctedUp = cross(f, right);
  return quaternionFromAxes(right, correctedUp, f);
}

export function buildReferenceQuaternions(path, reference = {}) {
  const { quaternions, targets, up = [0, 1, 0] } = reference;
  if (!isFiniteVector(up)) {
    throw new ValidationError('orientation.up must be a finite [x, y, z] vector.');
  }
  if (quaternions !== undefined) {
    if (!Array.isArray(quaternions) || quaternions.length !== path.count) {
      throw new ValidationError('Orientation quaternions must contain one reference per non-duplicate control point.');
    }
    return quaternions.map((quaternion) => normalizeQuaternion(quaternion));
  }

  if (targets !== undefined) {
    if (!Array.isArray(targets) || targets.length !== path.count || !targets.every(isFiniteVector)) {
      throw new ValidationError('Orientation targets must contain one finite target point per non-duplicate control point.');
    }
    return path.points.map((point, index) => {
      const forward = subtract(targets[index], point);
      if (dot(forward, forward) <= Number.EPSILON) {
        throw new ValidationError('Orientation target must not coincide with its control point.');
      }
      return buildLookQuaternion(forward, up);
    });
  }

  return path.points.map((point, index) => {
    const u = path.closed
      ? index / path.count
      : index / Math.max(1, path.count - 1);
    return buildLookQuaternion(path.tangentAt(u), up);
  });
}

export class OrientationTrack {
  constructor(path, reference = {}) {
    this.path = path;
    this.quaternions = buildReferenceQuaternions(path, reference);
    this.arcStations = [0];
    let cumulative = 0;
    for (let i = 0; i < path.segmentCount - 1; i += 1) {
      cumulative += path.segmentLengths[i];
      this.arcStations.push(cumulative);
    }
    if (!path.closed && this.arcStations.length < path.count) {
      this.arcStations.push(path.arcLength);
    }
  }

  at(s) {
    const length = this.path.arcLength;
    const target = this.path.closed ? ((s % length) + length) % length : Math.max(0, Math.min(length, s));

    if (target >= length - 1e-12) {
      if (this.path.closed) {
        return this.quaternions[0].slice();
      }
      return this.quaternions[this.quaternions.length - 1].slice();
    }

    let start = 0;
    for (let i = 0; i < this.arcStations.length; i += 1) {
      if (target >= this.arcStations[i]) {
        start = i;
      }
    }
    if (this.path.closed && start >= this.quaternions.length - 1) {
      start = this.quaternions.length - 1;
    }

    let end;
    let segmentLength;
    if (this.path.closed && start === this.arcStations.length - 1) {
      end = 0;
      segmentLength = length - this.arcStations[start];
    } else {
      end = start + 1;
      segmentLength = this.arcStations[end] - this.arcStations[start];
    }

    const t = segmentLength <= 1e-12 ? 0 : (target - this.arcStations[start]) / segmentLength;
    return slerp(this.quaternions[start], this.quaternions[end], t);
  }
}
