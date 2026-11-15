import { ValidationError } from './lib/errors.js';
import { OrientationTrack } from './lib/orientation.js';
import { CatmullRomPath, validateLimits } from './lib/spline.js';
import { SpeedProfile } from './lib/speed.js';
import { isFiniteNumber, isFiniteVector, norm, scale } from './lib/vectors.js';

function validateTimeInput(points, times, duration) {
  if (times === undefined && duration === undefined) {
    return undefined;
  }
  if (times !== undefined && duration !== undefined) {
    throw new ValidationError('Provide either times or duration, not both.');
  }
  if (duration !== undefined) {
    if (!isFiniteNumber(duration) || duration <= 0) {
      throw new ValidationError('duration must be a finite positive number of seconds.');
    }
    return duration;
  }

  if (!Array.isArray(times) || times.length !== points.length || !times.every(isFiniteNumber)) {
    throw new ValidationError('times must contain one finite timestamp for each original control point.');
  }
  for (let i = 1; i < times.length; i += 1) {
    if (times[i] <= times[i - 1]) {
      throw new ValidationError('times must be strictly increasing.');
    }
  }
  const span = times[times.length - 1] - times[0];
  if (span <= 0) {
    throw new ValidationError('The times span must be positive.');
  }
  return span;
}

export function planSplineMotion(points, options = {}) {
  if (!Array.isArray(points) || points.length < 2 || !points.every(isFiniteVector)) {
    throw new ValidationError('points must contain at least two finite [x, y, z] vectors.');
  }

  const closed = options.closed === true;
  const requestedDuration = validateTimeInput(points, options.times, options.duration);
  validateLimits(options);

  const path = new CatmullRomPath(points, { closed });
  const speedProfile = new SpeedProfile(path, options);
  const duration = speedProfile.applyDuration(requestedDuration);
  const orientationTrack = new OrientationTrack(path, options.orientation ?? {});

  const sampleCount = closed
    ? Math.max(1, Math.ceil(duration / options.sampleDt))
    : Math.max(2, Math.ceil(duration / options.sampleDt) + 1);
  const samples = [];
  let maxObservedAcceleration = 0;

  for (let index = 0; index < sampleCount; index += 1) {
    const time = closed
      ? index * duration / sampleCount
      : index * duration / (sampleCount - 1);
    const arcLength = speedProfile.sAtTime(time);
    const u = path.sToU(arcLength);
    const position = path.positionAt(u);
    const tangent = path.tangentAt(u);
    const speed = speedProfile.speedAtS(arcLength);
    const acceleration = speedProfile.accelerationAtS(arcLength);
    const accelerationMagnitude = norm(acceleration);
    maxObservedAcceleration = Math.max(maxObservedAcceleration, accelerationMagnitude);

    samples.push({
      index,
      time,
      arcLength,
      position,
      tangent,
      velocity: scale(tangent, speed),
      speed,
      acceleration,
      accelerationMagnitude,
      curvature: path.curvatureAt(u),
      quaternion: orientationTrack.at(arcLength)
    });
  }

  return {
    closed,
    duration,
    sampleCount,
    sampleDt: options.sampleDt,
    arcLength: path.arcLength,
    arcLengthError: path.totalLengthError,
    minimumDuration: speedProfile.minimumDuration,
    speedScale: speedProfile.speedScale,
    duplicateControlPointIndexes: path.duplicateIndexes,
    samples,
    diagnostics: {
      arcLengthIntegration: {
        method: 'adaptive Gauss-Legendre quadrature, 3-point and 5-point error estimate',
        segmentErrors: path.segmentLengthErrors
      },
      intervalCount: speedProfile.intervalCount,
      maxObservedAcceleration,
      maxObservedCurvature: Math.max(...speedProfile.curvatures),
      closingDuplicateRemoved: path.closingDuplicate
    }
  };
}

export { buildLookQuaternion } from './lib/orientation.js';
export {
  InfeasibleCurveError,
  InfeasibleTimeError,
  PlannerError,
  ValidationError
} from './lib/errors.js';
