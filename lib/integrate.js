const GL3 = {
  nodes: [-Math.sqrt(3 / 5), 0, Math.sqrt(3 / 5)],
  weights: [5 / 9, 8 / 9, 5 / 9]
};

const GL5 = {
  nodes: [
    -0.906179845938664,
    -0.5384693101056831,
    0,
    0.5384693101056831,
    0.906179845938664
  ],
  weights: [
    0.2369268850561891,
    0.4786286704993665,
    0.5688888888888889,
    0.4786286704993665,
    0.2369268850561891
  ]
};

function gaussLegendre(f, a, b, rule) {
  const midpoint = (a + b) / 2;
  const halfLength = (b - a) / 2;
  let sum = 0;
  for (let i = 0; i < rule.nodes.length; i += 1) {
    const x = midpoint + halfLength * rule.nodes[i];
    sum += rule.weights[i] * f(x);
  }
  return halfLength * sum;
}

export function adaptiveIntegrate(f, a, b, tolerance = 1e-10, depth = 0) {
  const coarse = gaussLegendre(f, a, b, GL3);
  const fine = gaussLegendre(f, a, b, GL5);
  const error = Math.abs(fine - coarse);
  const allowedError = tolerance * Math.max(1, Math.abs(fine));

  if (depth >= 16 || error <= allowedError) {
    return { value: fine, error };
  }

  const middle = (a + b) / 2;
  const left = adaptiveIntegrate(f, a, middle, tolerance, depth + 1);
  const right = adaptiveIntegrate(f, middle, b, tolerance, depth + 1);
  return { value: left.value + right.value, error: left.error + right.error };
}
