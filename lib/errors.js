export class PlannerError extends Error {
  constructor(message) {
    super(message);
    this.name = 'PlannerError';
  }
}

export class ValidationError extends PlannerError {
  constructor(message) {
    super(message);
    this.name = 'ValidationError';
  }
}

export class InfeasibleCurveError extends PlannerError {
  constructor(message, details = {}) {
    super(message);
    this.name = 'InfeasibleCurveError';
    this.details = details;
  }
}

export class InfeasibleTimeError extends PlannerError {
  constructor(message, details = {}) {
    super(message);
    this.name = 'InfeasibleTimeError';
    this.details = details;
  }
}
