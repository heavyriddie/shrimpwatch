// Companion page URL — hosted on GitHub Pages
export const COMPANION_URL = 'https://heavyriddie.github.io/shrimpwatch';

// MoveNet keypoint indices
export const KEYPOINT = {
  NOSE: 0,
  L_EYE: 1,
  R_EYE: 2,
  L_EAR: 3,
  R_EAR: 4,
  L_SHOULDER: 5,
  R_SHOULDER: 6,
  L_ELBOW: 7,
  R_ELBOW: 8,
  L_WRIST: 9,
  R_WRIST: 10,
  L_HIP: 11,
  R_HIP: 12,
  L_KNEE: 13,
  R_KNEE: 14,
  L_ANKLE: 15,
  R_ANKLE: 16
};

// Minimum keypoint confidence to use
export const MIN_CONFIDENCE = 0.3;

// Message types between service worker, offscreen doc, popup, content scripts
export const MSG = {
  // Service worker <-> Offscreen
  CAPTURE_AND_ANALYZE: 'CAPTURE_AND_ANALYZE',
  POSTURE_RESULT: 'POSTURE_RESULT',
  POSTURE_ERROR: 'POSTURE_ERROR',
  STOP_STREAMS: 'STOP_STREAMS',
  STREAMS_STOPPED: 'STREAMS_STOPPED',
  ENUMERATE_CAMERAS: 'ENUMERATE_CAMERAS',
  CAMERA_LIST: 'CAMERA_LIST',
  CAMERA_DISCONNECTED: 'CAMERA_DISCONNECTED',

  // Service worker <-> Content script
  BLUR: 'BLUR',
  UNBLUR: 'UNBLUR',
  SHRIMP_BANANAS: 'SHRIMP_BANANAS',
  SHRIMP_BANANAS_STOP: 'SHRIMP_BANANAS_STOP',
  BLUR_DISMISSED: 'BLUR_DISMISSED',

  // Popup <-> Service worker
  GET_STATUS: 'GET_STATUS',
  STATUS_RESPONSE: 'STATUS_RESPONSE',
  TOGGLE_MONITORING: 'TOGGLE_MONITORING',
  CHECK_NOW: 'CHECK_NOW',
  SAVE_SETTINGS: 'SAVE_SETTINGS',
  START_CALIBRATION: 'START_CALIBRATION',
  CALIBRATION_COMPLETE: 'CALIBRATION_COMPLETE'
};

// Alert modes
export const ALERT_MODE = {
  NOTIFICATION: 'notification',
  BLUR: 'blur',
  SHRIMP_BANANAS: 'shrimp_bananas'
};

// Default settings
export const DEFAULT_SETTINGS = {
  monitoringEnabled: false,
  checkIntervalSeconds: 60,
  keepStreamOpen: true,

  // Alert settings
  alertMode: ALERT_MODE.NOTIFICATION,
  notifyOnPoor: true,
  notifyOnRecovery: true,
  poorThreshold: 40,
  goodThreshold: 75,
  consecutivePoorBeforeAlert: 2,

  // Blur settings
  blurIntensity: 5,
  blurAutoRemove: true,

  // Camera settings
  frontCameraId: null,
  sideCameraId: null,
  captureWidth: 640,
  captureHeight: 480
};

// Default session state
export const DEFAULT_SESSION = {
  startedAt: null,
  consecutivePoorCount: 0,
  currentStreakType: null,
  currentStreakStartedAt: null,
  lastScore: null,
  lastCheckAt: null,
  isBlurActive: false,
  blurredTabId: null
};

// Posture messages - poor
export const POOR_POSTURE_MESSAGES = [
  "You're slouching! Sit up tall, shoulders back.",
  "Head's drifting forward. Tuck your chin, lift your chest.",
  "Time for a posture reset! Imagine a string pulling you up from the top of your head.",
  "Your back called. It wants its posture back.",
  "Shrimp detected! Unshrimp yourself.",
  "You're turning into a human question mark. Straighten up!",
  "Your spine is not a banana. Straighten it out!",
  "Slouch alert! Your future self will thank you for sitting up.",
];

// Posture messages - good
export const GOOD_POSTURE_MESSAGES = [
  "Great posture! Keep it up!",
  "Looking tall and proud! Nice work.",
  "Your spine is happy right now.",
  "Posture game: strong.",
  "The shrimp approves. Stay straight!",
];

// Notification sounds (Web Audio API frequencies for jingle)
export const SHRIMP_JINGLE_NOTES = [
  // "SHRIMP SHRIMP SHRIMP" - bouncy ascending pattern
  { freq: 523.25, duration: 0.15, gap: 0.05 },  // C5
  { freq: 659.25, duration: 0.15, gap: 0.05 },  // E5
  { freq: 783.99, duration: 0.15, gap: 0.1 },   // G5
  { freq: 523.25, duration: 0.15, gap: 0.05 },  // C5
  { freq: 659.25, duration: 0.15, gap: 0.05 },  // E5
  { freq: 783.99, duration: 0.15, gap: 0.1 },   // G5
  { freq: 523.25, duration: 0.15, gap: 0.05 },  // C5
  { freq: 659.25, duration: 0.15, gap: 0.05 },  // E5
  { freq: 1046.50, duration: 0.3, gap: 0 },     // C6 (finale!)
];
