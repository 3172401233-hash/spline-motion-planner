import { InfeasibleCurveError, InfeasibleTimeError, ValidationError } from './errors.js';
import { add, dot, scale, subtract } from './vectors.js';

function curvatureLimit(maxSpeed, maxAcceleration, curvature) {
  if (!Number.isFinite(curvature) || curvature < 0) {
    throw new InfeasibleCurveError('The spline contains undefined curvature.', { curvature });
  }
  if (curvature <= Number.EPSILON) {
    return maxSpeed;
  }
  const normalSpeed = Math.sqrt(maxAcceleration / curvature);
  return Math.min(maxSpeed, normalSpeed);
}

function segmentSpeedFromStart(startSpeed, segmentCurvature, distanceStep, maxAcceleration, upperBound) {
  const evaluate = (candidate) => {
    const tangential = Math.abs(candidate * candidate - startSpeed * startSpeed) / (2 * distanceStep);
    const normal = candidate * candidate * segmentCurvature;
    return Math.hypot(tangential, normal);
  };

  if (evaluate(upperBound) <= maxAcceleration + 1e-10) {
    return upperBound;
  }

  let low = startSpeed;
  let high = upperBound;
  for (let i = 0; i < 28; i += 1) {
    const middle = (low + high) / 2;
    if (evaluate(middle) <= maxAcceleration) {
      low = middle;
    } else {
      high = middle;
    }
  }
  return low;
}

function segmentSpeedFromEnd(endSpeed, segmentCurvature, distanceStep, maxAcceleration, upperBound) {
  const evaluate = (candidate) => {
    const tangential = Math.abs(endSpeed * endSpeed - candidate * candidate) / (2 * distanceStep);
    const normal = candidate * candidate * segmentCurvature;
    return Math.hypot(tangential, normal);
  };

  if (evaluate(upperBound) <= maxAcceleration + 1e-10) {
    return upperBound;
  }

  let low = endSpeed;
  let high = upperBound;
  for (let i = 0; i < 28; i += 1) {
    const middle = (low + high) / 2;
    if (evaluate(middle) <= maxAcceleration) {
      low = middle;
    } else {
      high = middle;
    }
  }
  return low;
}

export class SpeedProfile {
  constructor(path, { maxSpeed, maxAcceleration, maxCurvature, sampleDt }) {
    this.path = path;
    this.maxSpeed = maxSpeed;
    this.maxAcceleration = maxAcceleration;
    this.maxCurvature = maxCurvature;

    const desiredStep = Math.max(sampleDt * maxSpeed, path.arcLength / 1000);
    const intervalCount = Math.max(
      200,
      path.segmentCount * 8,
      Math.ceil(path.arcLength / desiredStep)
    );
    this.intervalCount = Math.min(intervalCount, 8000);
    this.distanceStep = path.arcLength / this.intervalCount;
    this.stations = [];
    this.curvatures = [];
    this.limits = [];
    this.speeds = [];

    for (let i = 0; i <= this.intervalCount; i += 1) {
      const s = path.arcLength * i / this.intervalCount;
      const u = path.sToU(s);
      const curvature = path.curvatureAt(u);
      if (curvature > maxCurvature + 1e-9) {
        throw new InfeasibleCurveError(
          'The spline curvature exceeds maxCurvature; add control points to round the turn.',
          { arcLength: s, curvature, maxCurvature }
        );
      }
      this.stations.push(s);
      this.curvatures.push(curvature);
      this.limits.push(curvatureLimit(maxSpeed, maxAcceleration, curvature));
    }

    this.segmentCurvatures = [];
    for (let i = 0; i < this.intervalCount; i += 1) {
      const midpointS = (this.stations[i] + this.stations[i + 1]) / 2;
      const midpointCurvature = path.curvatureAt(path.sToU(midpointS));
      const peakCurvature = this.#curvaturePeak(
        this.stations[i],
        this.stations[i + 1],
        this.curvatures[i],
        midpointCurvature,
        this.curvatures[i + 1],
        0
      );
      this.segmentCurvatures.push(peakCurvature);
    }

    this.#solveOpenProfile();
    this.timeStations = this.#buildTimeStations(1);
    this.minimumDuration = this.timeStations[this.timeStations.length - 1];
    this.speedScale = 1;
  }

  #solveOpenProfile() {
    const speeds = this.limits.slice();
    if (!this.path.closed) {
      speeds[0] = 0;
      speeds[speeds.length - 1] = 0;

      for (let i = 0; i < this.intervalCount; i += 1) {
        speeds[i + 1] = Math.min(
          speeds[i + 1],
          segmentSpeedFromStart(
            speeds[i],
            this.segmentCurvatures[i],
            this.distanceStep,
            this.maxAcceleration,
            this.limits[i + 1]
          )
        );
      }

      for (let i = this.intervalCount - 1; i >= 0; i -= 1) {
        speeds[i] = Math.min(
          speeds[i],
          segmentSpeedFromEnd(
            speeds[i + 1],
            this.segmentCurvatures[i],
            this.distanceStep,
            this.maxAcceleration,
            this.limits[i]
          )
        );
      }
    } else {
      for (let sweep = 0; sweep < 24; sweep += 1) {
        const before = speeds.slice();
        for (let step = 0; step < this.intervalCount; step += 1) {
          const i = step % this.intervalCount;
          const j = (i + 1) % this.intervalCount;
          speeds[j] = Math.min(
            speeds[j],
            segmentSpeedFromStart(
              speeds[i],
              this.segmentCurvatures[i],
              this.distanceStep,
              this.maxAcceleration,
              this.limits[j]
            )
          );
        }
        for (let step = this.intervalCount - 1; step >= 0; step -= 1) {
          const i = step % this.intervalCount;
          const j = (i + 1) % this.intervalCount;
          speeds[i] = Math.min(
            speeds[i],
            segmentSpeedFromEnd(
              speeds[j],
              this.segmentCurvatures[i],
              this.distanceStep,
              this.maxAcceleration,
              this.limits[i]
            )
          );
        }
        const change = Math.max(...speeds.map((value, index) => Math.abs(value - before[index])));
        if (change <= 1e-10) {
          break;
        }
      }
      speeds[speeds.length - 1] = speeds[0];
    }
    this.speeds = speeds;
  }

  #curvaturePeak(leftS, rightS, leftCurvature, middleCurvature, rightCurvature, depth) {
    const peak = Math.max(leftCurvature, middleCurvature, rightCurvature);
    if (peak > this.maxCurvature + 1e-9) {
      throw new InfeasibleCurveError(
        'The spline curvature exceeds maxCurvature; add control points to round the turn.',
        { arcLength: (leftS + rightS) / 2, curvature: peak, maxCurvature: this.maxCurvature }
      );
    }

    if (depth >= 8) {
      return peak;
    }

    const variation = Math.max(
      Math.abs(middleCurvature - leftCurvature),
      Math.abs(rightCurvature - middleCurvature)
    );
    if (variation <= 0.02 * Math.max(1, peak)) {
      return peak;
    }

    const leftMiddleS = (leftS + (leftS + rightS) / 2) / 2;
    const rightMiddleS = ((leftS + rightS) / 2 + rightS) / 2;
    const leftMiddleCurvature = this.path.curvatureAt(this.path.sToU(leftMiddleS));
    const rightMiddleCurvature = this.path.curvatureAt(this.path.sToU(rightMiddleS));
    const middleS = (leftS + rightS) / 2;

    return Math.max(
      this.#curvaturePeak(leftS, middleS, leftCurvature, leftMiddleCurvature, middleCurvature, depth + 1),
      this.#curvaturePeak(middleS, rightS, middleCurvature, rightMiddleCurvature, rightCurvature, depth + 1)
    );
  }

  #intervalTime(i, scale) {
    const v0 = this.speeds[i] * scale;
    const v1 = this.speeds[i + 1] * scale;
    if (v0 <= Number.EPSILON && v1 <= Number.EPSILON) {
      throw new InfeasibleCurveError('A non-zero path segment has zero reachable speed at both endpoints.');
    }
    if (v0 <= Number.EPSILON || v1 <= Number.EPSILON) {
      return 2 * this.distanceStep / Math.max(v0, v1);
    }
    return 2 * this.distanceStep / (v0 + v1);
  }

  #buildTimeStations(scale) {
    const times = [0];
    for (let i = 0; i < this.intervalCount; i += 1) {
      times.push(times[i] + this.#intervalTime(i, scale));
    }
    return times;
  }

  applyDuration(duration) {
    if (duration === undefined) {
      this.speedScale = 1;
      this.timeStations = this.#buildTimeStations(1);
      return this.timeStations[this.timeStations.length - 1];
    }
    if (!Number.isFinite(duration) || duration <= 0) {
      throw new ValidationError('duration must be a finite positive number of seconds.');
    }
    if (duration < this.minimumDuration - 1e-9) {
      throw new InfeasibleTimeError(
        'Requested duration is shorter than the minimum duration allowed by speed and acceleration limits.',
        { requestedDuration: duration, minimumDuration: this.minimumDuration }
      );
    }
    this.speedScale = this.minimumDuration / duration;
    this.timeStations = this.#buildTimeStations(this.speedScale);
    this.timeStations[this.timeStations.length - 1] = duration;
    return duration;
  }

  #intervalIndexAtS(s) {
    const target = this.path.closed
      ? ((s % this.path.arcLength) + this.path.arcLength) % this.path.arcLength
      : Math.max(0, Math.min(this.path.arcLength, s));
    const index = Math.min(this.intervalCount - 1, Math.floor(target / this.distanceStep));
    return { index, target };
  }

  speedAtS(s) {
    const { index, target } = this.#intervalIndexAtS(s);
    const s0 = this.stations[index];
    const s1 = this.stations[index + 1];
    const alpha = (target - s0) / (s1 - s0);
    const v0 = this.speeds[index] * this.speedScale;
    const v1 = this.speeds[index + 1] * this.speedScale;
    return v0 + (v1 - v0) * alpha;
  }

  sAtTime(time) {
    const duration = this.timeStations[this.timeStations.length - 1];
    const targetTime = this.path.closed
      ? ((time % duration) + duration) % duration
      : Math.max(0, Math.min(duration, time));

    let low = 0;
    let high = this.timeStations.length - 1;
    while (low < high) {
      const middle = Math.floor((low + high) / 2);
      if (this.timeStations[middle] <= targetTime) {
        low = middle + 1;
      } else {
        high = middle;
      }
    }
    const interval = Math.max(0, low - 1);
    const t0 = this.timeStations[interval];
    const t1 = this.timeStations[interval + 1];
    const alpha = t1 === t0 ? 0 : (targetTime - t0) / (t1 - t0);
    return this.stations[interval] + (this.stations[interval + 1] - this.stations[interval]) * alpha;
  }

  accelerationAtS(s) {
    const { index } = this.#intervalIndexAtS(s);
    const v0 = this.speeds[index] * this.speedScale;
    const v1 = this.speeds[index + 1] * this.speedScale;
    const tangential = (v1 * v1 - v0 * v0) / (2 * this.distanceStep);
    const u = this.path.sToU(s);
    const first = this.path.derivativeAt(u);
    const second = this.path.secondDerivativeAt(u);
    const firstLengthSquared = dot(first, first);
    const projection = scale(first, dot(first, second) / firstLengthSquared);
    const curvatureVector = scale(subtract(second, projection), 1 / firstLengthSquared);
    const speed = this.speedAtS(s);
    const tangent = scale(first, 1 / Math.sqrt(firstLengthSquared));
    return add(scale(tangent, tangential), scale(curvatureVector, speed * speed));
  }
}
