import { MSG, ALERT_MODE } from '../shared/constants.js';

const $ = (sel) => document.querySelector(sel);

let nextCheckTimer = null;

// ── Init ──

document.addEventListener('DOMContentLoaded', async () => {
  await loadStatus();
  setupListeners();
  startNextCheckCountdown();
});

// ── Load status from service worker ──

async function loadStatus() {
  try {
    const status = await chrome.runtime.sendMessage({ target: 'background', type: MSG.GET_STATUS });
    if (!status) return;

    const { settings, session, hasCalibration, nextCheckAt } = status;

    // Monitoring toggle
    $('#monitorToggle').checked = settings.monitoringEnabled;

    // Score display
    if (session.lastScore != null) {
      updateScore(session.lastScore);
      $('#statsSection').style.display = 'flex';
    }

    // Next check countdown
    if (nextCheckAt && settings.monitoringEnabled) {
      $('#nextCheckRow').style.display = 'flex';
      updateNextCheckLabel(nextCheckAt);
    }

    // Load daily stats
    const { dailySummaries } = await chrome.storage.local.get('dailySummaries');
    const today = new Date().toISOString().slice(0, 10);
    const todayStats = dailySummaries?.[today];
    if (todayStats) {
      const goodPct = todayStats.totalChecks > 0
        ? Math.round((todayStats.goodChecks / todayStats.totalChecks) * 100) : 0;
      $('#todayPercent').textContent = goodPct + '%';
      $('#checksToday').textContent = todayStats.totalChecks;
    }

    // Streak
    if (session.currentStreakStartedAt && session.currentStreakType) {
      const mins = Math.round((Date.now() - session.currentStreakStartedAt) / 60000);
      $('#currentStreak').textContent = mins < 60 ? `${mins}m` : `${Math.floor(mins / 60)}h${mins % 60}m`;
    }

    // Calibration status
    if (!hasCalibration) {
      $('#checkNowBtn').disabled = true;
      $('#settingsBtn').textContent = 'Calibrate (Required)';
    }

    // Camera permission check
    try {
      const perm = await navigator.permissions.query({ name: 'camera' });
      if (perm.state === 'prompt') {
        $('#permissionBanner').style.display = 'block';
      }
    } catch (e) {
      // permissions API might not be available
    }

    updateStatusText(settings.monitoringEnabled, hasCalibration);
  } catch (e) {
    console.error('Failed to load status:', e);
    $('#statusText').textContent = 'Error loading status';
  }
}

function updateScore(score) {
  const circumference = 2 * Math.PI * 52; // r=52
  const offset = circumference - (score / 100) * circumference;

  let color = '#4CAF50'; // green
  let label = 'Good';
  if (score < 40) { color = '#F44336'; label = 'Poor'; }
  else if (score < 75) { color = '#FF9800'; label = 'Fair'; }

  $('#scoreValue').textContent = score;
  $('#scoreLabel').textContent = label;
  $('#scoreArc').setAttribute('stroke', color);
  $('#scoreArc').setAttribute('stroke-dashoffset', offset);
}

function updateNextCheckLabel(nextCheckAt) {
  const remaining = Math.max(0, Math.round((nextCheckAt - Date.now()) / 1000));
  $('#nextCheckLabel').textContent = `Next check in ${remaining}s`;
}

function startNextCheckCountdown() {
  if (nextCheckTimer) clearInterval(nextCheckTimer);
  nextCheckTimer = setInterval(async () => {
    const alarm = await chrome.alarms.get('posture-check');
    if (alarm) {
      updateNextCheckLabel(alarm.scheduledTime);
    }
  }, 1000);
}

function updateStatusText(monitoring, hasCalibration) {
  if (!hasCalibration) {
    $('#statusText').textContent = 'Calibration needed to start monitoring';
  } else if (monitoring) {
    $('#statusText').textContent = 'Monitoring active';
  } else {
    $('#statusText').textContent = 'Monitoring paused';
  }
}

// ── Event listeners ──

function setupListeners() {
  // Monitor toggle
  $('#monitorToggle').addEventListener('change', async (e) => {
    const result = await chrome.runtime.sendMessage({ target: 'background', type: MSG.TOGGLE_MONITORING });
    if (result) {
      const show = result.enabled;
      $('#nextCheckRow').style.display = show ? 'flex' : 'none';
      if (show) startNextCheckCountdown();
    }
  });

  // Check now
  $('#checkNowBtn').addEventListener('click', async () => {
    $('#checkNowBtn').textContent = 'Checking...';
    $('#checkNowBtn').disabled = true;
    await chrome.runtime.sendMessage({ target: 'background', type: MSG.CHECK_NOW });
    setTimeout(async () => {
      await loadStatus();
      $('#checkNowBtn').textContent = 'Check Now';
      $('#checkNowBtn').disabled = false;
    }, 2000);
  });

  // Settings / Calibrate — opens options page
  $('#settingsBtn').addEventListener('click', () => {
    chrome.runtime.openOptionsPage();
  });

  // Grant camera permission
  $('#grantCameraBtn').addEventListener('click', async () => {
    try {
      const stream = await navigator.mediaDevices.getUserMedia({ video: true });
      stream.getTracks().forEach(t => t.stop());
      $('#permissionBanner').style.display = 'none';
      await loadCameras();
    } catch (e) {
      console.error('Camera permission denied:', e);
      $('#permissionBanner').querySelector('p').textContent =
        'Permission denied. Check browser settings.';
    }
  });
}
