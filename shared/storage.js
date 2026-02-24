import { DEFAULT_SETTINGS, DEFAULT_SESSION } from './constants.js';

// Typed wrappers around chrome.storage.local

export async function getSettings() {
  const { settings } = await chrome.storage.local.get('settings');
  return { ...DEFAULT_SETTINGS, ...settings };
}

export async function saveSettings(settings) {
  const current = await getSettings();
  await chrome.storage.local.set({ settings: { ...current, ...settings } });
}

export async function getCalibration() {
  const { calibration } = await chrome.storage.local.get('calibration');
  return calibration || null;
}

export async function saveCalibration(calibration) {
  await chrome.storage.local.set({
    calibration: { ...calibration, calibratedAt: Date.now(), version: 1 }
  });
}

export async function getSession() {
  const { session } = await chrome.storage.local.get('session');
  return { ...DEFAULT_SESSION, ...session };
}

export async function saveSession(session) {
  const current = await getSession();
  await chrome.storage.local.set({ session: { ...current, ...session } });
}

export async function getRecentChecks() {
  const { recentChecks } = await chrome.storage.local.get('recentChecks');
  return recentChecks || [];
}

export async function addCheck(check) {
  const checks = await getRecentChecks();
  checks.push(check);
  // Keep last 200 entries
  if (checks.length > 200) checks.splice(0, checks.length - 200);
  await chrome.storage.local.set({ recentChecks: checks });
}

export async function getDailySummaries() {
  const { dailySummaries } = await chrome.storage.local.get('dailySummaries');
  return dailySummaries || {};
}

export async function updateDailySummary(score, status) {
  const summaries = await getDailySummaries();
  const today = new Date().toISOString().slice(0, 10);
  const hour = new Date().getHours();

  if (!summaries[today]) {
    summaries[today] = {
      totalChecks: 0,
      goodChecks: 0,
      fairChecks: 0,
      poorChecks: 0,
      totalScore: 0,
      averageScore: 0,
      longestGoodStreakMinutes: 0,
      alertsSent: 0,
      selfCorrections: 0,
      hourlyScores: {},
      hourlyCounts: {}
    };
  }

  const s = summaries[today];
  s.totalChecks++;
  s.totalScore += score;
  s.averageScore = Math.round(s.totalScore / s.totalChecks);

  if (status === 'good') s.goodChecks++;
  else if (status === 'fair') s.fairChecks++;
  else if (status === 'poor') s.poorChecks++;

  // Hourly tracking
  if (!s.hourlyScores[hour]) {
    s.hourlyScores[hour] = 0;
    s.hourlyCounts[hour] = 0;
  }
  s.hourlyScores[hour] += score;
  s.hourlyCounts[hour]++;

  // Prune summaries older than 90 days
  const cutoff = new Date();
  cutoff.setDate(cutoff.getDate() - 90);
  const cutoffStr = cutoff.toISOString().slice(0, 10);
  for (const date of Object.keys(summaries)) {
    if (date < cutoffStr) delete summaries[date];
  }

  await chrome.storage.local.set({ dailySummaries: summaries });
}

export async function getCameras() {
  const { cameras } = await chrome.storage.local.get('cameras');
  return cameras || { lastEnumeration: [], lastEnumeratedAt: null };
}

export async function saveCameras(cameraList) {
  await chrome.storage.local.set({
    cameras: { lastEnumeration: cameraList, lastEnumeratedAt: Date.now() }
  });
}
