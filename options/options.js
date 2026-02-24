import { MSG, ALERT_MODE, DEFAULT_SETTINGS, COMPANION_URL } from '../shared/constants.js';
import { getSettings, saveSettings, getCalibration, saveCalibration, getDailySummaries, getRecentChecks } from '../shared/storage.js';

const $ = (sel) => document.querySelector(sel);

let calibrationStream = null;

// ── Init ──

document.addEventListener('DOMContentLoaded', async () => {
  await loadSettings();
  await loadCalibrationStatus();
  await loadMonitoringStatus();
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

// ── Monitoring toggle ──

async function loadMonitoringStatus() {
  try {
    const status = await chrome.runtime.sendMessage({ target: 'background', type: MSG.GET_STATUS });
    if (!status) return;

    const enabled = status.settings.monitoringEnabled;
    $('#monitoringToggle').checked = enabled;
    updateMonitoringUI(enabled, status.hasCalibration);
  } catch (e) {
    $('#monitoringStatus').textContent = 'Could not load status';
  }
}

function updateMonitoringUI(enabled, hasCalibration) {
  const card = document.querySelector('.monitoring-card');
  if (enabled) {
    card.classList.add('active');
    $('#monitoringStatus').textContent = 'Active — checking your posture';
  } else if (!hasCalibration) {
    card.classList.remove('active');
    $('#monitoringStatus').textContent = 'Calibrate below to start monitoring';
  } else {
    card.classList.remove('active');
    $('#monitoringStatus').textContent = 'Paused';
  }
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
    // Show posture progression guide
    $('#postureProgression').style.display = 'block';
  }
}

async function startCalibration() {
  $('#calibrationUI').style.display = 'block';
  $('#startCalibrationBtn').style.display = 'none';
  $('#stopCalibrationBtn').style.display = 'inline-block';
  $('#calibrationResult').style.display = 'none';

  // Show the right section based on current role
  updateCalibrationView();

  // Enumerate system cameras for front view
  try {
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

    // Start preview if front camera is selected
    if (document.querySelector('input[name="calRole"]:checked').value === 'front') {
      await startPreview(cameras[0]?.deviceId);
    }
  } catch (e) {
    $('#calibrationGuide').querySelector('p').textContent =
      'Could not access camera. Please grant permission.';
  }
}

function updateCalibrationView() {
  const role = document.querySelector('input[name="calRole"]:checked').value;
  $('#frontCameraSection').style.display = role === 'front' ? 'block' : 'none';
  $('#sideCameraSection').style.display = role === 'side' ? 'block' : 'none';
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
    video.style.transform = 'scaleX(-1)';
    $('#calibrationCanvas').style.transform = 'scaleX(-1)';
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
  if (video) video.srcObject = null;
}

async function captureBaseline(role, deviceId, progressEl, progressFillEl, resultEl, btnEl) {
  if (!deviceId) return;

  await chrome.runtime.sendMessage({ target: 'background', type: 'ENSURE_OFFSCREEN' });

  progressEl.style.display = 'block';
  btnEl.disabled = true;
  const origText = btnEl.textContent;
  btnEl.textContent = 'Capturing... sit still!';

  let progress = 0;
  const progressTimer = setInterval(() => {
    progress = Math.min(progress + 3, 90);
    progressFillEl.style.width = progress + '%';
  }, 100);

  try {
    const result = await chrome.runtime.sendMessage({
      target: 'offscreen',
      type: 'CALIBRATE',
      cameraRole: role,
      deviceId: deviceId,
      numFrames: 5,
      intervalMs: 500
    });

    clearInterval(progressTimer);
    progressFillEl.style.width = '100%';

    if (result?.baseline) {
      const existing = await getCalibration() || {};
      existing[role] = result.baseline;
      existing.calibratedAt = Date.now();
      await saveCalibration(existing);

      const settingsKey = role === 'front' ? 'frontCameraId' : 'sideCameraId';
      await saveSettings({ [settingsKey]: deviceId });

      chrome.runtime.sendMessage({ target: 'background', type: MSG.CALIBRATION_COMPLETE });
      await loadCalibrationStatus();
      // Update monitoring toggle — calibration auto-starts monitoring
      $('#monitoringToggle').checked = true;
      updateMonitoringUI(true, true);

      // Show success with guidance
      const updated = await getCalibration();
      const hasFront = !!updated?.front;
      const hasSide = !!updated?.side;

      resultEl.style.display = 'block';
      resultEl.style.color = '#4CAF50';

      // Show posture progression guide
      $('#postureProgression').style.display = 'block';

      const otherDone = role === 'front' ? hasSide : hasFront;

      if (role === 'front' && otherDone) {
        resultEl.innerHTML =
          `<strong>Front camera calibrated!</strong> Both cameras are now active. ` +
          `ShrimpWatch will check your posture automatically and alert you if it deteriorates.`;
      } else if (role === 'front') {
        resultEl.innerHTML =
          `<strong>Front camera calibrated!</strong> Monitoring is active. ` +
          `For better accuracy, also add a side camera — switch to "Side Camera" above.`;
      } else if (role === 'side' && otherDone) {
        resultEl.innerHTML =
          `<strong>Side camera calibrated!</strong> Both cameras are now active. ` +
          `ShrimpWatch will check your posture automatically and alert you if it deteriorates.`;
      } else {
        resultEl.innerHTML =
          `<strong>Side camera calibrated!</strong> Monitoring is active. ` +
          `Don't forget to calibrate your front camera too — switch to "Front Camera" above.`;
      }
    } else {
      throw new Error(result?.error || 'No baseline returned');
    }
  } catch (e) {
    clearInterval(progressTimer);
    resultEl.style.display = 'block';
    resultEl.textContent = 'Calibration failed: ' + e.message;
    resultEl.style.color = '#F44336';
  }

  btnEl.disabled = false;
  btnEl.textContent = origText;
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
  // Monitoring toggle
  $('#monitoringToggle').addEventListener('change', async () => {
    const result = await chrome.runtime.sendMessage({ target: 'background', type: MSG.TOGGLE_MONITORING });
    if (result) {
      const cal = await getCalibration();
      updateMonitoringUI(result.enabled, !!cal);
    }
  });

  // Calibration
  $('#startCalibrationBtn').addEventListener('click', startCalibration);
  $('#stopCalibrationBtn').addEventListener('click', stopCalibration);

  // Front camera: capture baseline
  $('#captureBaselineBtn').addEventListener('click', () => {
    captureBaseline('front', $('#calibrationCameraSelect').value,
      $('#captureProgress'), $('#captureProgressFill'),
      $('#calibrationResult'), $('#captureBaselineBtn'));
  });

  // Side camera: capture baseline from phone
  $('#captureSideBaselineBtn').addEventListener('click', () => {
    captureBaseline('side', '__webrtc__',
      $('#sideProgress'), $('#sideProgressFill'),
      $('#sideCalibrationResult'), $('#captureSideBaselineBtn'));
  });

  $('#calibrationCameraSelect').addEventListener('change', (e) => {
    startPreview(e.target.value);
  });

  document.querySelectorAll('input[name="calRole"]').forEach(radio => {
    radio.addEventListener('change', () => {
      updateCalibrationView();
      const role = document.querySelector('input[name="calRole"]:checked').value;
      if (role === 'front') {
        const deviceId = $('#calibrationCameraSelect').value;
        if (deviceId) startPreview(deviceId);
      } else {
        stopPreview();
      }
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
    const btn = $('#testBananasBtn');
    btn.disabled = true;
    btn.textContent = 'Launching...';

    try {
      // Find any regular web page tab
      const allTabs = await chrome.tabs.query({ currentWindow: true });
      let tab = allTabs.find(t => /^https?:\/\//.test(t.url || ''));

      if (!tab) {
        // Open a new tab and wait for it to fully load
        tab = await chrome.tabs.create({ url: 'https://example.com', active: true });
        await new Promise((resolve) => {
          const listener = (tabId, info) => {
            if (tabId === tab.id && info.status === 'complete') {
              chrome.tabs.onUpdated.removeListener(listener);
              resolve();
            }
          };
          chrome.tabs.onUpdated.addListener(listener);
          setTimeout(resolve, 5000); // fallback timeout
        });
      } else {
        await chrome.tabs.update(tab.id, { active: true });
      }

      // Inject content script and trigger bananas with retry
      const sendBananas = async (retries = 3) => {
        for (let i = 0; i < retries; i++) {
          try {
            await chrome.scripting.executeScript({
              target: { tabId: tab.id },
              files: ['content/blur-overlay.js']
            });
          } catch (e) { /* might already be loaded */ }
          // Give the content script time to register its listener
          await new Promise(r => setTimeout(r, 500));
          try {
            await chrome.tabs.sendMessage(tab.id, { type: 'SHRIMP_BANANAS', score: 15 });
            return; // success
          } catch (e) {
            if (i === retries - 1) throw e;
          }
        }
      };
      await sendBananas();
    } catch (e) {
      alert('Could not test: ' + e.message);
    }

    btn.disabled = false;
    btn.textContent = 'Test SHRIMP GOES BANANAS';
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
