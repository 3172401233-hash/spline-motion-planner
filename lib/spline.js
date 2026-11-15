import { adaptiveIntegrate } from './integrate.js';
import { InfeasibleCurveError, ValidationError } from './errors.js';
import {
  add,
  cross,
  distance,
  isFiniteNumber,
  isFiniteVector,
  norm,
  scale,
  subtract
} from './vectors.js';

const EPSILON = 1e-12;

export function validatePoints(points, closed) {
  if (!Array.isArray(points) || points.length < 2) {
    throw new ValidationError('At least two position control points are required.');
  }
  if (!points.every(isFiniteVector)) {
    throw new ValidationError('Position control points must all be finite [x, y, z] vectors.');
  }
  if (typeof closed !== 'boolean') {
    throw new ValidationError('Closed must be a boolean.');
  }

  const sourceIndexes = [0];
  const unique = [points[0]];
  const skippedZeroLength = [];

  for (let i = 1; i < points.length; i += 1) {
    const previous = points[unique.length - 1];
    const separation = distance(previous, points[i]);
    if (separation <= EPSILON) {
      skippedZeroLength.push(i);
    } else {
      unique.push(points[i]);
      sourceIndexes.push(i);
    }
  }

  let closingDuplicate = false;
  if (closed && unique.length > 1 && distance(unique[0], unique[unique.length - 1]) <= EPSILON) {
    unique.pop();
    sourceIndexes.pop();
    closingDuplicate = true;
  }

  if (unique.length < 2 || (closed && unique.length < 3)) {
    throw new ValidationError('After removing duplicate and zero-length control points, the path has no non-degenerate segment.');
  }

  return { points: unique, sourceIndexes, skippedZeroLength, closingDuplicate };
}

function wrapIndex(index, count, closed) {
  if (closed) {
    return ((index % count) + count) % count;
  }
  return Math.max(0, Math.min(count - 1, index));
}

function centripetalKnots(p0, p1) {
  const d = distance(p0, p1);
  return Math.sqrt(Math.max(d, EPSILON));
}

export class CatmullRomPath {
  constructor(points, { closed = false } = {}) {
    const cleaned = validatePoints(points, closed);
    this.closed = closed;
    this.points = cleaned.points;
    this.sourceIndexes = cleaned.sourceIndexes;
    this.duplicateIndexes = cleaned.skippedZeroLength;
    this.closingDuplicate = cleaned.closingDuplicate;
    this.count = this.points.length;
    this.segmentCount = closed ? this.count : this.count - 1;
    this.knots = this.#buildKnots();
    this.arcLength = 0;
    this.segmentLengths = [];
    this.segmentLengthErrors = [];
    this.totalLengthError = 0;
    this.#buildArcLength();
    this.arcTable = this.#buildArcTable();
  }

  #buildKnots() {
    const knots = [0];
    for (let i = 1; i <= this.segmentCount; i += 1) {
      const previousIndex = wrapIndex(i - 1, this.count, this.closed);
      const currentIndex = wrapIndex(i, this.count, this.closed);
      knots.push(knots[i - 1] + centripetalKnots(this.points[previousIndex], this.points[currentIndex]));
    }
    return knots;
  }

  #segmentPoint(segmentIndex, localU) {
    const i = segmentIndex;
    let p0;
    if (!this.closed && i === 0) {
      p0 = scale(subtract(scale(this.points[0], 2), this.points[1]), 1);
    } else {
      p0 = this.points[wrapIndex(i - 1, this.count, this.closed)];
    }
    const p1 = this.points[wrapIndex(i, this.count, this.closed)];
    const p2 = this.points[wrapIndex(i + 1, this.count, this.closed)];
    let p3;
    if (!this.closed && i === this.segmentCount - 1) {
      p3 = scale(subtract(scale(this.points[this.count - 1], 2), this.points[this.count - 2]), 1);
    } else {
      p3 = this.points[wrapIndex(i + 2, this.count, this.closed)];
    }

    const t0 = 0;
    const t1 = centripetalKnots(p0, p1);
    const t2 = t1 + centripetalKnots(p1, p2);
    const t3 = t2 + centripetalKnots(p2, p3);
    const t = t1 + (t2 - t1) * localU;

    const a1 = this.#interpolate(p0, p1, t0, t1, t);
    const a2 = this.#interpolate(p1, p2, t1, t2, t);
    const a3 = this.#interpolate(p2, p3, t2, t3, t);
    const b1 = this.#interpolate(a1, a2, t0, t2, t);
    const b2 = this.#interpolate(a2, a3, t1, t3, t);
    return this.#interpolate(b1, b2, t1, t2, t);
  }

  #interpolate(a, b, ta, tb, t) {
    const denominator = tb - ta;
    const alpha = (t - ta) / denominator;
    return add(a, scale(subtract(b, a), alpha));
  }

  positionAt(u) {
    const clamped = this.closed ? this.periodicU(u) : Math.max(0, Math.min(1, u));
    const scaled = clamped * this.segmentCount;
    const segmentIndex = Math.min(this.segmentCount - 1, Math.floor(scaled));
    const localU = scaled - segmentIndex;
    return this.#segmentPoint(segmentIndex, localU);
  }

  derivativeAt(u) {
    const h = 1e-5;
    if (!this.closed && u <= h) {
      const rightU = Math.min(1, u + h);
      return scale(subtract(this.positionAt(rightU), this.positionAt(u)), 1 / (rightU - u));
    }
    if (!this.closed && u >= 1 - h) {
      const leftU = Math.max(0, u - h);
      return scale(subtract(this.positionAt(u), this.positionAt(leftU)), 1 / (u - leftU));
    }

    const leftU = this.closed ? u - h : Math.max(0, u - h);
    const rightU = this.closed ? u + h : Math.min(1, u + h);
    const denominator = rightU - leftU;
    if (denominator <= EPSILON) {
      const finiteH = 1e-6;
      if (u <= EPSILON) {
        return scale(subtract(this.positionAt(finiteH), this.positionAt(0)), 1 / finiteH);
      }
      return scale(subtract(this.positionAt(1), this.positionAt(1 - finiteH)), 1 / finiteH);
    }
    return scale(subtract(this.positionAt(rightU), this.positionAt(leftU)), 1 / denominator);
  }

  secondDerivativeAt(u) {
    const h = 1e-4;
    let leftU = u - h;
    let rightU = u + h;
    let leftH = h;
    let rightH = h;
    if (!this.closed) {
      leftU = Math.max(0, leftU);
      rightU = Math.min(1, rightU);
      leftH = u - leftU;
      rightH = rightU - u;
    }

    const left = this.derivativeAt(leftU);
    const right = this.derivativeAt(rightU);
    if (leftH <= EPSILON || rightH <= EPSILON) {
      const finiteH = 1e-4;
      if (u <= EPSILON) {
        const d0 = this.derivativeAt(0);
        const d1 = this.derivativeAt(finiteH);
        return scale(subtract(d1, d0), 1 / finiteH);
      }
      const d0 = this.derivativeAt(1 - finiteH);
      const d1 = this.derivativeAt(1);
      return scale(subtract(d1, d0), 1 / finiteH);
    }
    return scale(subtract(right, left), 1 / (leftH + rightH));
  }

  tangentAt(u) {
    const d = this.derivativeAt(u);
    const length = norm(d);
    if (length <= EPSILON) {
      throw new InfeasibleCurveError('The path contains a zero-length tangent or an exact cusp; use a rounded turn instead.');
    }
    return scale(d, 1 / length);
  }

  curvatureAt(u) {
    const first = this.derivativeAt(u);
    const second = this.secondDerivativeAt(u);
    const numerator = norm(cross(first, second));
    const denominator = norm(first) ** 3;
    if (denominator <= EPSILON) {
      throw new InfeasibleCurveError('Curvature is undefined because the path has a zero-length tangent.');
    }
    return numerator / denominator;
  }

  #buildArcLength() {
    let total = 0;
    let totalError = 0;
    for (let i = 0; i < this.segmentCount; i += 1) {
      const result = adaptiveIntegrate((u) => {
        return norm(this.derivativeAt(u));
      }, i / this.segmentCount, (i + 1) / this.segmentCount, 1e-11);
      this.segmentLengths.push(result.value);
      this.segmentLengthErrors.push(result.error);
      total += result.value;
      totalError += result.error;
    }
    this.arcLength = total;
    this.totalLengthError = totalError;
  }

  #buildArcTable() {
    const subdivisions = 96;
    const table = [{ u: 0, s: 0 }];
    let cumulative = 0;

    for (let segment = 0; segment < this.segmentCount; segment += 1) {
      const u0 = segment / this.segmentCount;
      const u1 = (segment + 1) / this.segmentCount;
      for (let step = 1; step <= subdivisions; step += 1) {
        const u = u0 + (u1 - u0) * step / subdivisions;
        const previousU = table[table.length - 1].u;
        cumulative += adaptiveIntegrate(
          (parameter) => norm(this.derivativeAt(parameter)),
          previousU,
          u,
          1e-12
        ).value;
        table.push({ u, s: cumulative });
      }
    }

    table[table.length - 1].s = this.arcLength;
    return table;
  }

  sToU(arcLength) {
    const target = Math.max(0, Math.min(this.arcLength, arcLength));
    if (target <= EPSILON) {
      return 0;
    }
    if (target >= this.arcLength - EPSILON) {
      return 1;
    }
    let lowIndex = 0;
    let highIndex = this.arcTable.length - 1;
    while (lowIndex < highIndex) {
      const middleIndex = Math.floor((lowIndex + highIndex) / 2);
      if (this.arcTable[middleIndex].s < target) {
        lowIndex = middleIndex + 1;
      } else {
        highIndex = middleIndex;
      }
    }

    const upper = this.arcTable[lowIndex];
    const lower = this.arcTable[Math.max(0, lowIndex - 1)];
    if (upper.s === lower.s) {
      return upper.u;
    }

    const alpha = (target - lower.s) / (upper.s - lower.s);
    return lower.u + (upper.u - lower.u) * alpha;
  }

  pointAtS(s) {
    return this.positionAt(this.sToU(s));
  }

  periodicU(u) {
    if (!this.closed) {
      return Math.max(0, Math.min(1, u));
    }
    return ((u % 1) + 1) % 1;
  }
}

export function validateLimits(options) {
  const numericOptions = [
    ['maxSpeed', options.maxSpeed],
    ['maxAcceleration', options.maxAcceleration],
    ['maxCurvature', options.maxCurvature],
    ['sampleDt', options.sampleDt]
  ];
  for (const [name, value] of numericOptions) {
    if (!isFiniteNumber(value) || value <= 0) {
      throw new ValidationError(`${name} must be a finite positive number.`);
    }
  }
}
