import { MSG, ALERT_MODE, DEFAULT_SETTINGS, DEFAULT_SESSION, POOR_POSTURE_MESSAGES, GOOD_POSTURE_MESSAGES } from '../shared/constants.js';
import { getSettings, saveSettings, getCalibration, getSession, saveSession, addCheck, updateDailySummary } from '../shared/storage.js';
import { sendToActiveTab } from '../shared/messages.js';

const ALARM_NAME = 'posture-check';

// ── Offscreen document management ──

async function ensureOffscreenDocument() {
  const contexts = await chrome.runtime.getContexts({
    contextTypes: ['OFFSCREEN_DOCUMENT']
  });
  if (contexts.length === 0) {
    await chrome.offscreen.createDocument({
      url: 'offscreen/offscreen.html',
      reasons: ['USER_MEDIA'],
      justification: 'Camera access for posture analysis'
    });
    console.log('[ShrimpWatch] Offscreen document created');
    // Wait for scripts to load and message listener to register
    await waitForOffscreen();
  }
}

async function waitForOffscreen(retries = 10) {
  for (let i = 0; i < retries; i++) {
    try {
      const resp = await chrome.runtime.sendMessage({ target: 'offscreen', type: 'PING' });
      if (resp?.type === 'PONG') return;
    } catch (e) { /* not ready yet */ }
    await new Promise(r => setTimeout(r, 300));
  }
  console.warn('[ShrimpWatch] Offscreen document may not be ready');
}

// ── Alarm management ──

async function ensureAlarm(intervalSeconds) {
  const existing = await chrome.alarms.get(ALARM_NAME);
  const intervalMinutes = Math.max(0.5, intervalSeconds / 60); // Chrome minimum is 30s
  if (!existing || Math.abs(existing.periodInMinutes - intervalMinutes) > 0.01) {
    await chrome.alarms.clear(ALARM_NAME);
    chrome.alarms.create(ALARM_NAME, {
      delayInMinutes: intervalMinutes,
      periodInMinutes: intervalMinutes
    });
    console.log(`[ShrimpWatch] Alarm set: every ${intervalSeconds}s`);
  }
}

async function startMonitoring() {
  const settings = await getSettings();
  await ensureOffscreenDocument();
  await ensureAlarm(settings.checkIntervalSeconds);
  await saveSession({ startedAt: Date.now() });
  chrome.action.setBadgeText({ text: 'ON' });
  chrome.action.setBadgeBackgroundColor({ color: '#4CAF50' });
  console.log('[ShrimpWatch] Monitoring started');
}

async function stopMonitoring() {
  await chrome.alarms.clear(ALARM_NAME);
  chrome.action.setBadgeText({ text: 'OFF' });
  chrome.action.setBadgeBackgroundColor({ color: '#9E9E9E' });

  const session = await getSession();
  if (session.isBlurActive && session.blurredTabId) {
    try {
      await chrome.tabs.sendMessage(session.blurredTabId, { type: MSG.UNBLUR });
      await chrome.tabs.sendMessage(session.blurredTabId, { type: MSG.SHRIMP_BANANAS_STOP });
    } catch (e) { /* tab may not exist */ }
  }
  await saveSession({ ...DEFAULT_SESSION });
  console.log('[ShrimpWatch] Monitoring stopped');
}

// ── Posture check ──

async function performPostureCheck() {
  const settings = await getSettings();
  const calibration = await getCalibration();

  if (!calibration) {
    console.warn('[ShrimpWatch] No calibration data, skipping check');
    return;
  }

  await ensureOffscreenDocument();

  const cameras = {};
  if (settings.frontCameraId) cameras.front = settings.frontCameraId;
  if (settings.sideCameraId) cameras.side = settings.sideCameraId;

  if (Object.keys(cameras).length === 0) {
    console.warn('[ShrimpWatch] No cameras configured');
    return;
  }

  try {
    const result = await chrome.runtime.sendMessage({
      target: 'offscreen',
      type: MSG.CAPTURE_AND_ANALYZE,
      cameras,
      calibration
    });

    if (result?.type === MSG.POSTURE_ERROR) {
      console.error('[ShrimpWatch] Posture check error:', result.error);
      return;
    }

    if (result?.score != null) {
      await handlePostureResult(result, settings);
    }
  } catch (err) {
    console.error('[ShrimpWatch] Check failed:', err);
  }
}

async function handlePostureResult(result, settings) {
  const session = await getSession();
  const { score, status } = result;

  // Store result
  await addCheck({ timestamp: Date.now(), score, status });
  await updateDailySummary(score, status);

  // Update badge
  chrome.action.setBadgeText({ text: String(score) });
  chrome.action.setBadgeBackgroundColor({
    color: status === 'good' ? '#4CAF50' : status === 'fair' ? '#FF9800' : '#F44336'
  });

  const previousConsecutivePoor = session.consecutivePoorCount;

  if (status === 'poor') {
    const newCount = session.consecutivePoorCount + 1;
    await saveSession({
      consecutivePoorCount: newCount,
      lastScore: score,
      lastCheckAt: Date.now(),
      currentStreakType: 'poor',
      currentStreakStartedAt: session.currentStreakType === 'poor'
        ? session.currentStreakStartedAt : Date.now()
    });

    if (newCount >= settings.consecutivePoorBeforeAlert && settings.notifyOnPoor) {
      await triggerAlert(settings, score);
    }
  } else {
    const wasRecovery = previousConsecutivePoor >= settings.consecutivePoorBeforeAlert;

    // Remove blur/bananas if posture improved
    if (session.isBlurActive && settings.blurAutoRemove) {
      await removeAlert(session);
    }

    await saveSession({
      consecutivePoorCount: 0,
      lastScore: score,
      lastCheckAt: Date.now(),
      currentStreakType: status === 'good' ? 'good' : 'fair',
      currentStreakStartedAt: session.currentStreakType === 'good' || session.currentStreakType === 'fair'
        ? session.currentStreakStartedAt : Date.now()
    });

    // Positive reinforcement on recovery
    if (wasRecovery && settings.notifyOnRecovery && status === 'good') {
      const msg = GOOD_POSTURE_MESSAGES[Math.floor(Math.random() * GOOD_POSTURE_MESSAGES.length)];
      chrome.notifications.create('posture-praise-' + Date.now(), {
        type: 'basic',
        iconUrl: 'icons/icon128.png',
        title: 'ShrimpWatch: Nice recovery!',
        message: msg,
        priority: 0,
        silent: true
      });
    }
  }
}

async function triggerAlert(settings, score) {
  const mode = settings.alertMode;

  if (mode === ALERT_MODE.NOTIFICATION || mode === ALERT_MODE.BLUR) {
    // Always show notification
    const msg = POOR_POSTURE_MESSAGES[Math.floor(Math.random() * POOR_POSTURE_MESSAGES.length)];
    chrome.notifications.create('posture-alert', {
      type: 'basic',
      iconUrl: 'icons/icon128.png',
      title: `ShrimpWatch: Score ${score}/100`,
      message: msg,
      buttons: [
        { title: 'I fixed it!' },
        { title: 'Snooze 5 min' }
      ],
      priority: 2,
      requireInteraction: true
    });
  }

  if (mode === ALERT_MODE.BLUR) {
    const result = await sendToActiveTab(MSG.BLUR, { level: settings.blurIntensity });
    if (result?.success) {
      const [tab] = await chrome.tabs.query({ active: true, currentWindow: true });
      await saveSession({ isBlurActive: true, blurredTabId: tab?.id });
    }
  }

  if (mode === ALERT_MODE.SHRIMP_BANANAS) {
    const result = await sendToActiveTab(MSG.SHRIMP_BANANAS, { score });
    if (result?.success) {
      const [tab] = await chrome.tabs.query({ active: true, currentWindow: true });
      await saveSession({ isBlurActive: true, blurredTabId: tab?.id });
    }
    // Also show notification
    chrome.notifications.create('posture-alert', {
      type: 'basic',
      iconUrl: 'icons/icon128.png',
      title: 'SHRIMP GOES BANANAS!',
      message: `Score: ${score}/100. THE SHRIMP IS UPSET. SIT UP STRAIGHT!`,
      buttons: [
        { title: 'I fixed it!' },
        { title: 'Snooze 5 min' }
      ],
      priority: 2,
      requireInteraction: true
    });
  }
}

async function removeAlert(session) {
  if (session.blurredTabId) {
    try {
      await chrome.tabs.sendMessage(session.blurredTabId, { type: MSG.UNBLUR });
      await chrome.tabs.sendMessage(session.blurredTabId, { type: MSG.SHRIMP_BANANAS_STOP });
    } catch (e) { /* tab may not exist */ }
  }
  await saveSession({ isBlurActive: false, blurredTabId: null });
}

// ── Event listeners ──

// Alarm handler
chrome.alarms.onAlarm.addListener((alarm) => {
  if (alarm.name === ALARM_NAME) {
    performPostureCheck();
  }
});

// Install / startup
chrome.runtime.onInstalled.addListener(async (details) => {
  if (details.reason === 'install') {
    await chrome.storage.local.set({
      settings: DEFAULT_SETTINGS,
      calibration: null,
      recentChecks: [],
      dailySummaries: {},
      session: DEFAULT_SESSION
    });
    chrome.action.setBadgeText({ text: 'NEW' });
    chrome.action.setBadgeBackgroundColor({ color: '#2196F3' });
    // Open options page for onboarding
    chrome.runtime.openOptionsPage();
  }
});

chrome.runtime.onStartup.addListener(async () => {
  const settings = await getSettings();
  const calibration = await getCalibration();
  if (calibration && settings.monitoringEnabled) {
    await startMonitoring();
  }
});

// Notification button handler
chrome.notifications.onButtonClicked.addListener(async (id, buttonIndex) => {
  if (id === 'posture-alert') {
    const settings = await getSettings();
    if (buttonIndex === 0) {
      // "I fixed it!" - do a quick recheck in 10 seconds
      chrome.alarms.create('recheck', { delayInMinutes: 10 / 60 });
      // Remove any blur/bananas
      const session = await getSession();
      if (session.isBlurActive) await removeAlert(session);
    } else if (buttonIndex === 1) {
      // "Snooze 5 min" - restart alarm with delay
      await chrome.alarms.clear(ALARM_NAME);
      chrome.alarms.create(ALARM_NAME, {
        delayInMinutes: 5,
        periodInMinutes: Math.max(0.5, settings.checkIntervalSeconds / 60)
      });
      const session = await getSession();
      if (session.isBlurActive) await removeAlert(session);
    }
    chrome.notifications.clear(id);
  }
});

// Recheck alarm
chrome.alarms.onAlarm.addListener((alarm) => {
  if (alarm.name === 'recheck') {
    performPostureCheck();
  }
});

// Message handler from popup/options/content scripts
chrome.runtime.onMessage.addListener((message, sender, sendResponse) => {
  if (message.target && message.target !== 'background') return false;

  switch (message.type) {
    case MSG.GET_STATUS:
      (async () => {
        const settings = await getSettings();
        const session = await getSession();
        const calibration = await getCalibration();
        const alarm = await chrome.alarms.get(ALARM_NAME);
        sendResponse({
          type: MSG.STATUS_RESPONSE,
          settings,
          session,
          hasCalibration: !!calibration,
          calibratedAt: calibration?.calibratedAt,
          alarmScheduled: !!alarm,
          nextCheckAt: alarm?.scheduledTime
        });
      })();
      return true;

    case MSG.TOGGLE_MONITORING:
      (async () => {
        const settings = await getSettings();
        const newEnabled = !settings.monitoringEnabled;
        await saveSettings({ monitoringEnabled: newEnabled });
        if (newEnabled) {
          await startMonitoring();
        } else {
          await stopMonitoring();
        }
        sendResponse({ enabled: newEnabled });
      })();
      return true;

    case MSG.CHECK_NOW:
      performPostureCheck().then(() => sendResponse({ done: true }));
      return true;

    case MSG.SAVE_SETTINGS:
      (async () => {
        await saveSettings(message.settings);
        const settings = await getSettings();
        // Update alarm interval if monitoring is active
        if (settings.monitoringEnabled) {
          await ensureAlarm(settings.checkIntervalSeconds);
        }
        sendResponse({ saved: true });
      })();
      return true;

    case MSG.BLUR_DISMISSED:
      saveSession({ isBlurActive: false, blurredTabId: null });
      return false;

    case MSG.CALIBRATION_COMPLETE:
      (async () => {
        // If this is the first calibration and monitoring wasn't on, start it
        const settings = await getSettings();
        if (!settings.monitoringEnabled) {
          await saveSettings({ monitoringEnabled: true });
          await startMonitoring();
        }
        sendResponse({ started: true });
      })();
      return true;

    case 'ENSURE_OFFSCREEN':
      ensureOffscreenDocument().then(() => sendResponse({ ready: true }));
      return true;

    case 'INIT_PEER':
      (async () => {
        await ensureOffscreenDocument();
        const result = await chrome.runtime.sendMessage({
          target: 'offscreen', type: 'INIT_PEER'
        });
        sendResponse(result);
      })();
      return true;

    case 'DESTROY_PEER':
      chrome.runtime.sendMessage({ target: 'offscreen', type: 'DESTROY_PEER' });
      sendResponse({ done: true });
      return false;

    case 'PHONE_STATUS':
      (async () => {
        await ensureOffscreenDocument();
        const result = await chrome.runtime.sendMessage({
          target: 'offscreen', type: 'PHONE_STATUS'
        });
        sendResponse(result);
      })();
      return true;

    case 'PHONE_CAMERA_CONNECTED':
      // Phone connected via WebRTC, update badge
      console.log('[ShrimpWatch] Phone camera connected via WebRTC');
      return false;

    default:
      return false;
  }
});
