import { KEYPOINT, MIN_CONFIDENCE } from '../shared/constants.js';

// ── Geometry helpers ──

function kp(keypoints, index) {
  const p = keypoints[index];
  return p && p.score >= MIN_CONFIDENCE ? p : null;
}

function midpoint(a, b) {
  return { x: (a.x + b.x) / 2, y: (a.y + b.y) / 2 };
}

function distance(a, b) {
  return Math.sqrt((a.x - b.x) ** 2 + (a.y - b.y) ** 2);
}

function angleDeg(a, b, c) {
  const ba = { x: a.x - b.x, y: a.y - b.y };
  const bc = { x: c.x - b.x, y: c.y - b.y };
  const dot = ba.x * bc.x + ba.y * bc.y;
  const cross = ba.x * bc.y - ba.y * bc.x;
  return Math.atan2(cross, dot) * (180 / Math.PI);
}

// ── Front camera metrics ──

export function calculateFrontMetrics(keypoints, imageWidth, imageHeight) {
  const nose = kp(keypoints, KEYPOINT.NOSE);
  const lEye = kp(keypoints, KEYPOINT.L_EYE);
  const rEye = kp(keypoints, KEYPOINT.R_EYE);
  const lEar = kp(keypoints, KEYPOINT.L_EAR);
  const rEar = kp(keypoints, KEYPOINT.R_EAR);
  const lShoulder = kp(keypoints, KEYPOINT.L_SHOULDER);
  const rShoulder = kp(keypoints, KEYPOINT.R_SHOULDER);
  const lHip = kp(keypoints, KEYPOINT.L_HIP);
  const rHip = kp(keypoints, KEYPOINT.R_HIP);

  if (!lShoulder || !rShoulder) return null;

  const shoulderMid = midpoint(lShoulder, rShoulder);
  const shoulderDist = distance(lShoulder, rShoulder);
  if (shoulderDist < 1) return null; // Too close together to be meaningful

  const metrics = {};

  // 1. Shoulder slope: angle of shoulder line from horizontal
  metrics.shoulderSlopeDeg = Math.atan2(
    rShoulder.y - lShoulder.y,
    rShoulder.x - lShoulder.x
  ) * (180 / Math.PI);

  // 2. Nose-to-shoulder vertical ratio (primary slouch indicator)
  if (nose) {
    metrics.noseShoulderVerticalRatio =
      (shoulderMid.y - nose.y) / shoulderDist;
  }

  // 3. Nose horizontal offset from center
  if (nose) {
    metrics.noseHorizontalOffset =
      (nose.x - shoulderMid.x) / shoulderDist;
  }

  // 4. Head tilt via eye line
  if (lEye && rEye) {
    metrics.headTiltDeg = Math.atan2(
      rEye.y - lEye.y,
      rEye.x - lEye.x
    ) * (180 / Math.PI);
  }

  // 5. Ear-shoulder vertical alignment
  if (lEar) {
    metrics.leftEarShoulderVertical =
      (lShoulder.y - lEar.y) / shoulderDist;
  }
  if (rEar) {
    metrics.rightEarShoulderVertical =
      (rShoulder.y - rEar.y) / shoulderDist;
  }

  // 6. Torso inclination (if hips visible)
  if (lHip && rHip) {
    const hipMid = midpoint(lHip, rHip);
    metrics.torsoInclinationDeg = Math.atan2(
      shoulderMid.x - hipMid.x,
      hipMid.y - shoulderMid.y
    ) * (180 / Math.PI);
  }

  metrics.shoulderWidth = shoulderDist;
  return metrics;
}

// ── Side camera metrics ──

export function calculateSideMetrics(keypoints, imageWidth, imageHeight) {
  const ear = kp(keypoints, KEYPOINT.L_EAR) || kp(keypoints, KEYPOINT.R_EAR);
  const shoulder = kp(keypoints, KEYPOINT.L_SHOULDER) || kp(keypoints, KEYPOINT.R_SHOULDER);
  const hip = kp(keypoints, KEYPOINT.L_HIP) || kp(keypoints, KEYPOINT.R_HIP);
  const nose = kp(keypoints, KEYPOINT.NOSE);

  if (!shoulder) return null;

  const metrics = {};

  // 1. Neck inclination (forward head posture gold standard)
  if (ear) {
    metrics.neckInclinationDeg = Math.atan2(
      Math.abs(ear.x - shoulder.x),
      shoulder.y - ear.y
    ) * (180 / Math.PI);
  }

  // 2. Torso inclination from side
  if (hip) {
    metrics.torsoInclinationDeg = Math.atan2(
      Math.abs(shoulder.x - hip.x),
      hip.y - shoulder.y
    ) * (180 / Math.PI);
  }

  // 3. Ear-shoulder-hip angle
  if (ear && hip) {
    metrics.earShoulderHipAngle = angleDeg(ear, shoulder, hip);
  }

  // 4. Head forward offset
  if (nose && hip) {
    const torsoLength = distance(shoulder, hip);
    if (torsoLength > 1) {
      metrics.headForwardOffset = (nose.x - shoulder.x) / torsoLength;
    }
  }

  return metrics;
}

// ── Baseline creation ──

export function createBaseline(metricsArray) {
  if (!metricsArray || metricsArray.length === 0) return null;

  const filtered = metricsArray.filter(m => m !== null);
  if (filtered.length === 0) return null;

  const baseline = {};
  const keys = Object.keys(filtered[0]).filter(k => !k.startsWith('_'));

  for (const key of keys) {
    const values = filtered.map(m => m[key]).filter(v => v !== undefined && v !== null);
    if (values.length > 0) {
      baseline[key] = values.reduce((a, b) => a + b, 0) / values.length;
    }
  }

  // Store standard deviation for natural variance
  baseline._stddev = {};
  for (const key of keys) {
    const values = filtered.map(m => m[key]).filter(v => v != null);
    if (values.length > 1) {
      const mean = baseline[key];
      const variance = values.reduce((sum, v) => sum + (v - mean) ** 2, 0) / values.length;
      baseline._stddev[key] = Math.sqrt(variance);
    }
  }

  return baseline;
}

// ── Scoring ──

function scoreFrontMetrics(current, baseline) {
  const scores = [];

  // Shoulder slope deviation
  if (current.shoulderSlopeDeg != null && baseline.shoulderSlopeDeg != null) {
    const deviation = Math.abs(current.shoulderSlopeDeg - baseline.shoulderSlopeDeg);
    scores.push({
      metric: 'shoulderSlope',
      score: Math.max(0, 100 - (deviation / 15) * 100),
      weight: 1.0,
      deviation
    });
  }

  // Vertical slouch (most important front metric)
  if (current.noseShoulderVerticalRatio != null && baseline.noseShoulderVerticalRatio != null) {
    const ratio = current.noseShoulderVerticalRatio / baseline.noseShoulderVerticalRatio;
    scores.push({
      metric: 'verticalSlouch',
      score: Math.max(0, Math.min(100, ((ratio - 0.7) / 0.3) * 100)),
      weight: 2.5,
      deviation: 1 - ratio
    });
  }

  // Nose horizontal offset
  if (current.noseHorizontalOffset != null && baseline.noseHorizontalOffset != null) {
    const deviation = Math.abs(current.noseHorizontalOffset - baseline.noseHorizontalOffset);
    scores.push({
      metric: 'headCentering',
      score: Math.max(0, 100 - (deviation / 0.3) * 100),
      weight: 0.8,
      deviation
    });
  }

  // Head tilt
  if (current.headTiltDeg != null && baseline.headTiltDeg != null) {
    const deviation = Math.abs(current.headTiltDeg - baseline.headTiltDeg);
    scores.push({
      metric: 'headTilt',
      score: Math.max(0, 100 - (deviation / 20) * 100),
      weight: 0.8,
      deviation
    });
  }

  // Ear-shoulder vertical alignment
  const earDeviations = [];
  if (current.leftEarShoulderVertical != null && baseline.leftEarShoulderVertical != null) {
    earDeviations.push(current.leftEarShoulderVertical / baseline.leftEarShoulderVertical);
  }
  if (current.rightEarShoulderVertical != null && baseline.rightEarShoulderVertical != null) {
    earDeviations.push(current.rightEarShoulderVertical / baseline.rightEarShoulderVertical);
  }
  if (earDeviations.length > 0) {
    const avgRatio = earDeviations.reduce((a, b) => a + b, 0) / earDeviations.length;
    scores.push({
      metric: 'earShoulderAlignment',
      score: Math.max(0, Math.min(100, ((avgRatio - 0.7) / 0.3) * 100)),
      weight: 1.5,
      deviation: 1 - avgRatio
    });
  }

  return scores;
}

function scoreSideMetrics(current, baseline) {
  const scores = [];

  // Neck inclination (highest weight - clinical gold standard)
  if (current.neckInclinationDeg != null && baseline.neckInclinationDeg != null) {
    const deviation = current.neckInclinationDeg - baseline.neckInclinationDeg;
    scores.push({
      metric: 'neckInclination',
      score: Math.max(0, 100 - (Math.max(0, deviation) / 20) * 100),
      weight: 3.0,
      deviation
    });
  }

  // Torso inclination
  if (current.torsoInclinationDeg != null && baseline.torsoInclinationDeg != null) {
    const deviation = current.torsoInclinationDeg - baseline.torsoInclinationDeg;
    scores.push({
      metric: 'torsoInclination',
      score: Math.max(0, 100 - (Math.abs(deviation) / 25) * 100),
      weight: 2.0,
      deviation
    });
  }

  // Head forward offset
  if (current.headForwardOffset != null && baseline.headForwardOffset != null) {
    const deviation = current.headForwardOffset - baseline.headForwardOffset;
    scores.push({
      metric: 'headForwardOffset',
      score: Math.max(0, 100 - (Math.max(0, deviation) / 0.25) * 100),
      weight: 2.5,
      deviation
    });
  }

  return scores;
}

export function evaluatePosture(poseData, calibration) {
  if (!calibration) {
    return { score: null, status: 'not_calibrated' };
  }

  const scores = [];
  let frontMetrics = null;
  let sideMetrics = null;

  if (poseData.front) {
    frontMetrics = calculateFrontMetrics(
      poseData.front.keypoints,
      poseData.front.width,
      poseData.front.height
    );
    if (frontMetrics && calibration.front) {
      scores.push(...scoreFrontMetrics(frontMetrics, calibration.front));
    }
  }

  if (poseData.side) {
    sideMetrics = calculateSideMetrics(
      poseData.side.keypoints,
      poseData.side.width,
      poseData.side.height
    );
    if (sideMetrics && calibration.side) {
      scores.push(...scoreSideMetrics(sideMetrics, calibration.side));
    }
  }

  if (scores.length === 0) {
    return { score: null, status: 'no_data', metrics: { frontMetrics, sideMetrics } };
  }

  const totalWeight = scores.reduce((sum, s) => sum + s.weight, 0);
  const weightedScore = scores.reduce((sum, s) => sum + s.score * s.weight, 0) / totalWeight;
  const finalScore = Math.round(Math.max(0, Math.min(100, weightedScore)));

  return {
    score: finalScore,
    status: finalScore >= 75 ? 'good' : finalScore >= 40 ? 'fair' : 'poor',
    metrics: { frontMetrics, sideMetrics },
    breakdown: scores,
    timestamp: Date.now()
  };
}
