import { MSG, ALERT_MODE, DEFAULT_SETTINGS, COMPANION_URL } from '../shared/constants.js';
import { getSettings, saveSettings, getCalibration, saveCalibration, getDailySummaries, getRecentChecks } from '../shared/storage.js';

const $ = (sel) => document.querySelector(sel);

let calibrationStream = null;

// ── Init ──

document.addEventListener('DOMContentLoaded', async () => {
  await loadSettings();
  await loadCalibrationStatus();
  await loadAnalytics();
  setupListeners();
});

// ── Load settings ──

async function loadSettings() {
  const s = await getSettings();

  $('#checkInterval').value = s.checkIntervalSeconds;
  $('#keepStreamOpen').checked = s.keepStreamOpen;
  $('#alertMode').value = s.alertMode;
  $('#notifyOnPoor').checked = s.notifyOnPoor;
  $('#notifyOnRecovery').checked = s.notifyOnRecovery;
  $('#poorThreshold').value = s.poorThreshold;
  $('#poorThresholdVal').textContent = s.poorThreshold;
  $('#goodThreshold').value = s.goodThreshold;
  $('#goodThresholdVal').textContent = s.goodThreshold;
  $('#consecutivePoor').value = s.consecutivePoorBeforeAlert;
  $('#blurIntensity').value = s.blurIntensity;
  $('#blurIntensityVal').textContent = s.blurIntensity;
  $('#blurAutoRemove').checked = s.blurAutoRemove;

  updateAlertModeUI(s.alertMode);
}

function updateAlertModeUI(mode) {
  $('#blurSettings').style.display = (mode === 'blur') ? 'block' : 'none';
  $('#bananasPreview').style.display = (mode === 'shrimp_bananas') ? 'block' : 'none';
}

// ── Calibration ──

async function loadCalibrationStatus() {
  const cal = await getCalibration();
  if (cal) {
    const date = new Date(cal.calibratedAt).toLocaleString();
    const cams = [];
    if (cal.front) cams.push('Front');
    if (cal.side) cams.push('Side');
    $('#calibrationStatus').textContent =
      `Last calibrated: ${date} (${cams.join(' + ')} camera${cams.length > 1 ? 's' : ''})`;
    $('#startCalibrationBtn').textContent = 'Recalibrate';
  }
}

async function startCalibration() {
  $('#calibrationUI').style.display = 'block';
  $('#startCalibrationBtn').style.display = 'none';
  $('#stopCalibrationBtn').style.display = 'inline-block';
  $('#calibrationResult').style.display = 'none';

  // Enumerate cameras
  try {
    // Ensure we have permission first
    const testStream = await navigator.mediaDevices.getUserMedia({ video: true });
    testStream.getTracks().forEach(t => t.stop());

    const devices = await navigator.mediaDevices.enumerateDevices();
    const cameras = devices.filter(d => d.kind === 'videoinput');

    const select = $('#calibrationCameraSelect');
    select.innerHTML = '';
    for (const cam of cameras) {
      const opt = document.createElement('option');
      opt.value = cam.deviceId;
      opt.textContent = cam.label || `Camera ${cam.deviceId.slice(0, 8)}`;
      select.appendChild(opt);
    }

    // Start preview
    await startPreview(cameras[0]?.deviceId);
  } catch (e) {
    $('#calibrationGuide').querySelector('p').textContent =
      'Could not access camera. Please grant permission.';
  }
}

async function startPreview(deviceId) {
  stopPreview();
  if (!deviceId) return;

  try {
    calibrationStream = await navigator.mediaDevices.getUserMedia({
      video: {
        deviceId: { exact: deviceId },
        width: { ideal: 640 },
        height: { ideal: 480 }
      }
    });
    const video = $('#calibrationVideo');
    video.srcObject = calibrationStream;

    // Unflip for side camera
    const role = document.querySelector('input[name="calRole"]:checked').value;
    video.style.transform = role === 'front' ? 'scaleX(-1)' : 'none';
    $('#calibrationCanvas').style.transform = role === 'front' ? 'scaleX(-1)' : 'none';
  } catch (e) {
    console.error('Preview failed:', e);
  }
}

function stopPreview() {
  if (calibrationStream) {
    calibrationStream.getTracks().forEach(t => t.stop());
    calibrationStream = null;
  }
  const video = $('#calibrationVideo');
  video.srcObject = null;
}

async function captureBaseline() {
  const deviceId = $('#calibrationCameraSelect').value;
  const role = document.querySelector('input[name="calRole"]:checked').value;

  if (!deviceId) return;

  // Ask service worker to ensure offscreen document exists
  // (chrome.offscreen is only available in the service worker context)
  await chrome.runtime.sendMessage({ target: 'background', type: 'ENSURE_OFFSCREEN' });

  $('#captureProgress').style.display = 'block';
  $('#captureBaselineBtn').disabled = true;
  $('#captureBaselineBtn').textContent = 'Capturing... sit still!';

  // Animate progress
  let progress = 0;
  const progressTimer = setInterval(() => {
    progress = Math.min(progress + 3, 90);
    $('#captureProgressFill').style.width = progress + '%';
  }, 100);

  try {
    // Send calibration request to offscreen document
    const result = await chrome.runtime.sendMessage({
      target: 'offscreen',
      type: 'CALIBRATE',
      cameraRole: role,
      deviceId: deviceId,
      numFrames: 5,
      intervalMs: 500
    });

    clearInterval(progressTimer);
    $('#captureProgressFill').style.width = '100%';

    if (result?.baseline) {
      // Save calibration
      const existing = await getCalibration() || {};
      existing[role] = result.baseline;
      await saveCalibration(existing);

      // Also save camera ID
      const settingsKey = role === 'front' ? 'frontCameraId' : 'sideCameraId';
      await saveSettings({ [settingsKey]: deviceId });

      $('#calibrationResult').style.display = 'block';
      $('#calibrationResult').textContent =
        `${role === 'front' ? 'Front' : 'Side'} camera calibrated successfully!`;
      $('#calibrationResult').style.color = '#4CAF50';

      // Notify service worker
      chrome.runtime.sendMessage({ target: 'background', type: MSG.CALIBRATION_COMPLETE });

      await loadCalibrationStatus();
    } else {
      throw new Error(result?.error || 'No baseline returned');
    }
  } catch (e) {
    clearInterval(progressTimer);
    $('#calibrationResult').style.display = 'block';
    $('#calibrationResult').textContent = 'Calibration failed: ' + e.message;
    $('#calibrationResult').style.color = '#F44336';
  }

  $('#captureBaselineBtn').disabled = false;
  $('#captureBaselineBtn').textContent = 'Capture Baseline';
}

function stopCalibration() {
  stopPreview();
  $('#calibrationUI').style.display = 'none';
  $('#startCalibrationBtn').style.display = 'inline-block';
  $('#stopCalibrationBtn').style.display = 'none';
}

// ── Analytics ──

async function loadAnalytics() {
  const summaries = await getDailySummaries();
  const today = new Date().toISOString().slice(0, 10);
  const todayData = summaries[today];

  if (!todayData || todayData.totalChecks === 0) return;

  $('#analyticsSummary').querySelector('.no-data')?.remove();
  $('#analyticsDetails').style.display = 'block';
  $('#analyticsChart').style.display = 'block';

  // Fill stats
  $('#avgScore').textContent = todayData.averageScore;
  const goodPct = Math.round((todayData.goodChecks / todayData.totalChecks) * 100);
  $('#goodPercent').textContent = goodPct + '%';
  $('#totalChecksAnalytics').textContent = todayData.totalChecks;
  $('#selfCorrections').textContent = todayData.selfCorrections || 0;

  // Draw chart
  drawChart(summaries);
}

function drawChart(summaries) {
  const canvas = $('#analyticsChart');
  const ctx = canvas.getContext('2d');
  const dpr = window.devicePixelRatio || 1;

  canvas.width = canvas.offsetWidth * dpr;
  canvas.height = 200 * dpr;
  ctx.scale(dpr, dpr);

  const width = canvas.offsetWidth;
  const height = 200;
  const padding = { top: 20, right: 20, bottom: 30, left: 40 };

  // Get last 7 days
  const dates = [];
  for (let i = 6; i >= 0; i--) {
    const d = new Date();
    d.setDate(d.getDate() - i);
    dates.push(d.toISOString().slice(0, 10));
  }

  const scores = dates.map(d => summaries[d]?.averageScore ?? null);

  // Background
  ctx.fillStyle = '#fafafa';
  ctx.fillRect(0, 0, width, height);

  // Grid lines
  ctx.strokeStyle = '#eee';
  ctx.lineWidth = 1;
  for (let y = 0; y <= 100; y += 25) {
    const py = padding.top + ((100 - y) / 100) * (height - padding.top - padding.bottom);
    ctx.beginPath();
    ctx.moveTo(padding.left, py);
    ctx.lineTo(width - padding.right, py);
    ctx.stroke();

    ctx.fillStyle = '#bbb';
    ctx.font = '10px sans-serif';
    ctx.textAlign = 'right';
    ctx.fillText(y, padding.left - 6, py + 3);
  }

  // Color zones
  const chartH = height - padding.top - padding.bottom;
  const goodY = padding.top;
  const fairY = padding.top + (25 / 100) * chartH;
  const poorY = padding.top + (60 / 100) * chartH;

  ctx.globalAlpha = 0.08;
  ctx.fillStyle = '#4CAF50';
  ctx.fillRect(padding.left, goodY, width - padding.left - padding.right, fairY - goodY);
  ctx.fillStyle = '#FF9800';
  ctx.fillRect(padding.left, fairY, width - padding.left - padding.right, poorY - fairY);
  ctx.fillStyle = '#F44336';
  ctx.fillRect(padding.left, poorY, width - padding.left - padding.right, height - padding.bottom - poorY);
  ctx.globalAlpha = 1;

  // Plot line
  const pointsX = dates.map((_, i) =>
    padding.left + (i / (dates.length - 1)) * (width - padding.left - padding.right)
  );

  ctx.strokeStyle = '#FF6B35';
  ctx.lineWidth = 2.5;
  ctx.lineJoin = 'round';
  ctx.beginPath();

  let started = false;
  for (let i = 0; i < scores.length; i++) {
    if (scores[i] == null) continue;
    const x = pointsX[i];
    const y = padding.top + ((100 - scores[i]) / 100) * chartH;
    if (!started) { ctx.moveTo(x, y); started = true; }
    else ctx.lineTo(x, y);
  }
  ctx.stroke();

  // Points
  for (let i = 0; i < scores.length; i++) {
    if (scores[i] == null) continue;
    const x = pointsX[i];
    const y = padding.top + ((100 - scores[i]) / 100) * chartH;

    ctx.fillStyle = scores[i] >= 75 ? '#4CAF50' : scores[i] >= 40 ? '#FF9800' : '#F44336';
    ctx.beginPath();
    ctx.arc(x, y, 4, 0, Math.PI * 2);
    ctx.fill();
  }

  // X labels
  ctx.fillStyle = '#999';
  ctx.font = '10px sans-serif';
  ctx.textAlign = 'center';
  for (let i = 0; i < dates.length; i++) {
    const label = dates[i].slice(5); // MM-DD
    ctx.fillText(label, pointsX[i], height - 8);
  }
}

// ── Event listeners ──

function setupListeners() {
  // Calibration
  $('#startCalibrationBtn').addEventListener('click', startCalibration);
  $('#stopCalibrationBtn').addEventListener('click', stopCalibration);
  $('#captureBaselineBtn').addEventListener('click', captureBaseline);

  $('#calibrationCameraSelect').addEventListener('change', (e) => {
    startPreview(e.target.value);
  });

  document.querySelectorAll('input[name="calRole"]').forEach(radio => {
    radio.addEventListener('change', () => {
      const role = document.querySelector('input[name="calRole"]:checked').value;
      const video = $('#calibrationVideo');
      video.style.transform = role === 'front' ? 'scaleX(-1)' : 'none';
      $('#calibrationCanvas').style.transform = role === 'front' ? 'scaleX(-1)' : 'none';
      // Show phone camera option when side is selected
      $('#phoneCameraSection').style.display = role === 'side' ? 'block' : 'none';
    });
  });

  // Settings auto-save
  const settingsMap = {
    checkInterval: { key: 'checkIntervalSeconds', type: 'int' },
    keepStreamOpen: { key: 'keepStreamOpen', type: 'bool' },
    alertMode: { key: 'alertMode', type: 'str' },
    notifyOnPoor: { key: 'notifyOnPoor', type: 'bool' },
    notifyOnRecovery: { key: 'notifyOnRecovery', type: 'bool' },
    poorThreshold: { key: 'poorThreshold', type: 'int', display: 'poorThresholdVal' },
    goodThreshold: { key: 'goodThreshold', type: 'int', display: 'goodThresholdVal' },
    consecutivePoor: { key: 'consecutivePoorBeforeAlert', type: 'int' },
    blurIntensity: { key: 'blurIntensity', type: 'int', display: 'blurIntensityVal' },
    blurAutoRemove: { key: 'blurAutoRemove', type: 'bool' },
  };

  for (const [id, config] of Object.entries(settingsMap)) {
    const el = $(`#${id}`);
    const event = config.type === 'bool' ? 'change' : 'input';

    el.addEventListener(event, async () => {
      let value;
      if (config.type === 'bool') value = el.checked;
      else if (config.type === 'int') value = parseInt(el.value);
      else value = el.value;

      if (config.display) {
        $(`#${config.display}`).textContent = value;
      }

      await chrome.runtime.sendMessage({
        target: 'background',
        type: MSG.SAVE_SETTINGS,
        settings: { [config.key]: value }
      });

      // Update UI for alert mode changes
      if (id === 'alertMode') {
        updateAlertModeUI(value);
      }
    });
  }

  // Test SHRIMP GOES BANANAS
  $('#testBananasBtn').addEventListener('click', async () => {
    const [tab] = await chrome.tabs.query({ active: true, currentWindow: true });
    if (tab?.id) {
      try {
        await chrome.tabs.sendMessage(tab.id, { type: 'SHRIMP_BANANAS', score: 15 });
      } catch (e) {
        // If content script not loaded, try scripting API
        try {
          await chrome.scripting.executeScript({
            target: { tabId: tab.id },
            files: ['content/blur-overlay.js']
          });
          await chrome.tabs.sendMessage(tab.id, { type: 'SHRIMP_BANANAS', score: 15 });
        } catch (e2) {
          alert('Could not test on this tab. Try a regular web page.');
        }
      }
    }
  });

  // Clear data
  $('#clearDataBtn').addEventListener('click', async () => {
    if (confirm('Clear all ShrimpWatch data? This removes calibration, history, and settings.')) {
      await chrome.storage.local.clear();
      await chrome.storage.local.set({ settings: DEFAULT_SETTINGS });
      location.reload();
    }
  });

  // ── Phone camera WebRTC ──

  let phoneCheckInterval = null;

  $('#connectPhoneBtn').addEventListener('click', async () => {
    const btn = $('#connectPhoneBtn');
    const errEl = $('#phoneError');
    btn.disabled = true;
    btn.textContent = 'Setting up...';
    errEl.style.display = 'none';

    try {
      // Step 1: ensure offscreen document is created and ready
      await chrome.runtime.sendMessage({ target: 'background', type: 'ENSURE_OFFSCREEN' });

      // Step 2: send INIT_PEER directly (offscreen doc picks it up)
      const result = await chrome.runtime.sendMessage({ target: 'offscreen', type: 'INIT_PEER' });

      if (!result?.peerId) {
        throw new Error(result?.error || 'No peer ID returned');
      }

      const companionUrl = `${COMPANION_URL}/#peer=${result.peerId}`;
      const qrUrl = `https://api.qrserver.com/v1/create-qr-code/?size=200x200&data=${encodeURIComponent(companionUrl)}`;

      $('#qrCodeImg').src = qrUrl;
      $('#companionUrl').textContent = companionUrl;

      // Show QR code UI
      $('#phoneDisconnected').style.display = 'none';
      $('#phoneConnecting').style.display = 'block';

      // Poll for phone connection
      phoneCheckInterval = setInterval(async () => {
        try {
          const status = await chrome.runtime.sendMessage({ target: 'offscreen', type: 'PHONE_STATUS' });
          if (status?.connected) {
            clearInterval(phoneCheckInterval);
            phoneCheckInterval = null;
            $('#phoneConnecting').style.display = 'none';
            $('#phoneConnected').style.display = 'block';
          }
        } catch (e) { /* ignore */ }
      }, 2000);

      return; // success — don't reset button
    } catch (e) {
      console.error('Phone connection error:', e);
      errEl.textContent = 'Error: ' + e.message;
      errEl.style.display = 'block';
    }

    btn.disabled = false;
    btn.textContent = 'Connect Phone via QR Code';
  });

  $('#cancelPhoneBtn').addEventListener('click', () => {
    if (phoneCheckInterval) {
      clearInterval(phoneCheckInterval);
      phoneCheckInterval = null;
    }
    chrome.runtime.sendMessage({ target: 'background', type: 'DESTROY_PEER' });
    $('#phoneConnecting').style.display = 'none';
    $('#phoneDisconnected').style.display = 'block';
  });

  $('#disconnectPhoneBtn').addEventListener('click', () => {
    chrome.runtime.sendMessage({ target: 'background', type: 'DESTROY_PEER' });
    $('#phoneConnected').style.display = 'none';
    $('#phoneDisconnected').style.display = 'block';
  });

  // Check if phone is already connected
  (async () => {
    try {
      const status = await chrome.runtime.sendMessage({
        target: 'background', type: 'PHONE_STATUS'
      });
      if (status?.connected) {
        $('#phoneDisconnected').style.display = 'none';
        $('#phoneConnected').style.display = 'block';
      }
    } catch (e) { /* offscreen may not exist yet */ }
  })();
}
