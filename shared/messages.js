// Message helper functions

export function sendToOffscreen(type, data = {}) {
  return chrome.runtime.sendMessage({ target: 'offscreen', type, ...data });
}

export function sendToBackground(type, data = {}) {
  return chrome.runtime.sendMessage({ target: 'background', type, ...data });
}

export async function sendToActiveTab(type, data = {}) {
  const [tab] = await chrome.tabs.query({ active: true, currentWindow: true });
  if (tab?.id) {
    try {
      return await chrome.tabs.sendMessage(tab.id, { type, ...data });
    } catch (e) {
      // Tab might not have content script loaded (chrome:// pages, etc.)
      console.warn('Could not send to tab:', e.message);
      return null;
    }
  }
  return null;
}
