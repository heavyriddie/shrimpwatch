// Message helper functions

export function sendToOffscreen(type, data = {}) {
  return chrome.runtime.sendMessage({ target: 'offscreen', type, ...data });
}

export function sendToBackground(type, data = {}) {
  return chrome.runtime.sendMessage({ target: 'background', type, ...data });
}

export async function sendToActiveTab(type, data = {}) {
  // Find a suitable tab — prefer active tab, fall back to any regular web page
  const [activeTab] = await chrome.tabs.query({ active: true, currentWindow: true });
  let tab = activeTab;

  // Skip chrome://, chrome-extension://, edge://, about: pages
  if (tab && /^(chrome|chrome-extension|edge|about|devtools):/.test(tab.url || '')) {
    const [webTab] = await chrome.tabs.query({ currentWindow: true, url: ['http://*/*', 'https://*/*'] });
    tab = webTab || null;
  }

  if (!tab?.id) return null;

  try {
    return await chrome.tabs.sendMessage(tab.id, { type, ...data });
  } catch (e) {
    // Content script not loaded — inject it and retry
    try {
      await chrome.scripting.executeScript({
        target: { tabId: tab.id },
        files: ['content/blur-overlay.js']
      });
      return await chrome.tabs.sendMessage(tab.id, { type, ...data });
    } catch (e2) {
      console.warn('Could not send to tab:', e2.message);
      return null;
    }
  }
}
