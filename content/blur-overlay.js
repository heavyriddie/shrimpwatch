// Content script: blur overlay and SHRIMP GOES BANANAS mode
// Injected into all web pages

(function () {
  let blurOverlay = null;
  let bananasOverlay = null;
  let bananasInterval = null;
  let audioCtx = null;

  // ── Shrimp SVG (inline, used in bananas mode) ──
  const SHRIMP_SVG = `<svg viewBox="0 0 200 200" xmlns="http://www.w3.org/2000/svg">
    <g transform="translate(100,100)">
      <!-- Body -->
      <ellipse cx="0" cy="10" rx="45" ry="30" fill="#FF6B35" stroke="#D4522A" stroke-width="2"/>
      <!-- Tail segments -->
      <ellipse cx="-35" cy="30" rx="20" ry="12" fill="#FF8C5A" stroke="#D4522A" stroke-width="1.5" transform="rotate(-30 -35 30)"/>
      <ellipse cx="-50" cy="45" rx="15" ry="10" fill="#FFA07A" stroke="#D4522A" stroke-width="1.5" transform="rotate(-50 -50 45)"/>
      <path d="M-58 52 Q-75 65 -60 72 Q-45 65 -55 55" fill="#FFB899" stroke="#D4522A" stroke-width="1.5"/>
      <!-- Head -->
      <circle cx="35" cy="-5" r="22" fill="#FF8C5A" stroke="#D4522A" stroke-width="2"/>
      <!-- Eye -->
      <circle cx="42" cy="-12" r="6" fill="white" stroke="#333" stroke-width="1.5"/>
      <circle cx="44" cy="-13" r="3" fill="#333"/>
      <circle cx="45" cy="-14" r="1" fill="white"/>
      <!-- Antennae -->
      <path d="M45 -25 Q55 -50 65 -45" fill="none" stroke="#D4522A" stroke-width="2" stroke-linecap="round"/>
      <path d="M40 -25 Q45 -55 55 -55" fill="none" stroke="#D4522A" stroke-width="2" stroke-linecap="round"/>
      <!-- Legs -->
      <line x1="10" y1="35" x2="15" y2="55" stroke="#D4522A" stroke-width="2" stroke-linecap="round"/>
      <line x1="-5" y1="38" x2="-8" y2="58" stroke="#D4522A" stroke-width="2" stroke-linecap="round"/>
      <line x1="-20" y1="37" x2="-25" y2="55" stroke="#D4522A" stroke-width="2" stroke-linecap="round"/>
      <line x1="25" y1="30" x2="32" y2="48" stroke="#D4522A" stroke-width="2" stroke-linecap="round"/>
      <!-- Claws -->
      <path d="M50 5 Q65 -5 60 5 Q65 15 50 10" fill="#FF6B35" stroke="#D4522A" stroke-width="1.5"/>
    </g>
  </svg>`;

  // ── Blur overlay ──

  function showBlur(level) {
    if (blurOverlay) return;
    blurOverlay = document.createElement('div');
    blurOverlay.id = 'shrimpwatch-blur';
    blurOverlay.innerHTML = `
      <div style="
        position:fixed; inset:0; z-index:2147483646;
        backdrop-filter:blur(${level || 5}px);
        -webkit-backdrop-filter:blur(${level || 5}px);
        display:flex; align-items:center; justify-content:center;
        background:rgba(0,0,0,0.08);
        font-family:-apple-system,BlinkMacSystemFont,sans-serif;
      ">
        <div style="
          background:white; border-radius:20px; padding:40px 56px;
          text-align:center; box-shadow:0 12px 48px rgba(0,0,0,0.2);
          max-width:400px;
        ">
          <div style="font-size:56px; margin-bottom:12px;">🦐</div>
          <div style="font-size:20px; font-weight:700; color:#333; margin-bottom:8px;">
            Sit up straight!
          </div>
          <div style="font-size:14px; color:#666; margin-bottom:20px;">
            Fix your posture to unblur this page
          </div>
          <button id="shrimpwatch-dismiss-blur" style="
            padding:10px 28px; border:none; border-radius:10px;
            background:#f0f0f0; cursor:pointer; font-size:13px;
            font-weight:500; color:#555;
            transition:background 0.2s;
          ">Dismiss</button>
        </div>
      </div>
    `;
    document.documentElement.appendChild(blurOverlay);
    blurOverlay.querySelector('#shrimpwatch-dismiss-blur').addEventListener('click', () => {
      removeBlur();
      chrome.runtime.sendMessage({ type: 'BLUR_DISMISSED' });
    });
  }

  function removeBlur() {
    if (blurOverlay) {
      blurOverlay.remove();
      blurOverlay = null;
    }
  }

  // ── SHRIMP GOES BANANAS mode ──

  function playJingle() {
    try {
      audioCtx = new (window.AudioContext || window.webkitAudioContext)();
      const notes = [
        { freq: 523.25, dur: 0.12, gap: 0.04 },  // C5 - "SHRIMP"
        { freq: 659.25, dur: 0.12, gap: 0.04 },  // E5 - "SHRIMP"
        { freq: 783.99, dur: 0.12, gap: 0.08 },  // G5 - "SHRIMP"
        { freq: 523.25, dur: 0.12, gap: 0.04 },  // C5
        { freq: 659.25, dur: 0.12, gap: 0.04 },  // E5
        { freq: 783.99, dur: 0.12, gap: 0.08 },  // G5
        { freq: 523.25, dur: 0.12, gap: 0.04 },  // C5
        { freq: 659.25, dur: 0.12, gap: 0.04 },  // E5
        { freq: 1046.50, dur: 0.35, gap: 0 },    // C6 - finale!
      ];

      let time = audioCtx.currentTime + 0.05;
      for (const note of notes) {
        const osc = audioCtx.createOscillator();
        const gain = audioCtx.createGain();
        osc.type = 'square';
        osc.frequency.value = note.freq;
        gain.gain.setValueAtTime(0.15, time);
        gain.gain.exponentialRampToValueAtTime(0.01, time + note.dur);
        osc.connect(gain);
        gain.connect(audioCtx.destination);
        osc.start(time);
        osc.stop(time + note.dur);
        time += note.dur + note.gap;
      }
    } catch (e) {
      console.warn('[ShrimpWatch] Could not play jingle:', e);
    }
  }

  function startBananas(score) {
    if (bananasOverlay) return;

    bananasOverlay = document.createElement('div');
    bananasOverlay.id = 'shrimpwatch-bananas';

    // Create the chaos
    bananasOverlay.innerHTML = `
      <div id="sw-bananas-bg" style="
        position:fixed; inset:0; z-index:2147483647;
        display:flex; flex-direction:column; align-items:center; justify-content:center;
        font-family:'Comic Sans MS','Chalkboard SE',cursive,-apple-system,sans-serif;
        overflow:hidden;
        background:rgba(255,50,0,0.85);
        animation:sw-bg-flash 0.3s ease-in-out infinite alternate;
      ">
        <style>
          @keyframes sw-bg-flash {
            0% { background:rgba(255,50,0,0.85); }
            50% { background:rgba(255,200,0,0.9); }
            100% { background:rgba(255,100,0,0.85); }
          }
          @keyframes sw-shrimp-pulse {
            0% { transform:scale(0.6) rotate(-15deg); }
            50% { transform:scale(1.4) rotate(15deg); }
            100% { transform:scale(0.6) rotate(-15deg); }
          }
          @keyframes sw-text-bounce {
            0% { transform:scale(1) translateY(0); }
            25% { transform:scale(1.3) translateY(-20px); }
            50% { transform:scale(0.9) translateY(0); }
            75% { transform:scale(1.2) translateY(-10px); }
            100% { transform:scale(1) translateY(0); }
          }
          @keyframes sw-shake {
            0%, 100% { transform:translateX(0); }
            10% { transform:translateX(-10px) rotate(-2deg); }
            20% { transform:translateX(10px) rotate(2deg); }
            30% { transform:translateX(-10px) rotate(-1deg); }
            40% { transform:translateX(10px) rotate(1deg); }
            50% { transform:translateX(-5px); }
            60% { transform:translateX(5px); }
            70% { transform:translateX(-5px); }
            80% { transform:translateX(5px); }
            90% { transform:translateX(-2px); }
          }
          @keyframes sw-word-pop {
            0%, 100% { transform:scale(1); opacity:1; }
            50% { transform:scale(1.5); opacity:0.8; }
          }
          #sw-bananas-shrimp {
            animation:sw-shrimp-pulse 0.4s ease-in-out infinite;
            filter:drop-shadow(0 0 30px rgba(255,200,0,0.8));
          }
          #sw-bananas-title {
            animation:sw-text-bounce 0.5s ease-in-out infinite;
          }
          #sw-bananas-container {
            animation:sw-shake 0.5s ease-in-out infinite;
          }
          .sw-shrimp-word {
            display:inline-block;
            animation:sw-word-pop 0.3s ease-in-out infinite;
            text-shadow:4px 4px 0 #8B0000, -2px -2px 0 #FF4500;
          }
          .sw-shrimp-word:nth-child(2) { animation-delay:0.1s; }
          .sw-shrimp-word:nth-child(3) { animation-delay:0.2s; }
        </style>

        <div id="sw-bananas-container" style="text-align:center;">
          <!-- The shrimp -->
          <div id="sw-bananas-shrimp" style="width:200px; height:200px; margin:0 auto 20px;">
            ${SHRIMP_SVG}
          </div>

          <!-- SHRIMP SHRIMP SHRIMP -->
          <div style="margin-bottom:16px;">
            <span class="sw-shrimp-word" style="font-size:48px; font-weight:900; color:white; letter-spacing:4px;">SHRIMP </span>
            <span class="sw-shrimp-word" style="font-size:48px; font-weight:900; color:white; letter-spacing:4px;">SHRIMP </span>
            <span class="sw-shrimp-word" style="font-size:48px; font-weight:900; color:white; letter-spacing:4px;">SHRIMP</span>
          </div>

          <!-- SHRIMP GOES BANANAS -->
          <div id="sw-bananas-title" style="
            font-size:64px; font-weight:900; color:#FFD700;
            text-shadow:4px 4px 0 #8B0000, -2px -2px 0 #FF4500, 0 0 40px rgba(255,215,0,0.5);
            letter-spacing:6px; margin-bottom:24px;
          ">
            SHRIMP GOES BANANAS
          </div>

          <!-- Score -->
          <div style="font-size:24px; color:white; font-weight:700; margin-bottom:24px; text-shadow:2px 2px 0 rgba(0,0,0,0.3);">
            Posture Score: ${score}/100 — SIT UP NOW!
          </div>

          <!-- Dismiss -->
          <button id="sw-bananas-dismiss" style="
            padding:16px 48px; border:none; border-radius:16px;
            background:white; cursor:pointer; font-size:18px;
            font-weight:800; color:#FF4500;
            box-shadow:0 6px 24px rgba(0,0,0,0.3);
            text-transform:uppercase; letter-spacing:2px;
            transition:transform 0.1s;
          ">OK OK I'LL SIT UP</button>
        </div>
      </div>
    `;

    document.documentElement.appendChild(bananasOverlay);

    // Play the jingle
    playJingle();

    // Replay jingle every 2 seconds
    bananasInterval = setInterval(playJingle, 2000);

    bananasOverlay.querySelector('#sw-bananas-dismiss').addEventListener('click', () => {
      stopBananas();
      chrome.runtime.sendMessage({ type: 'BLUR_DISMISSED' });
    });
  }

  function stopBananas() {
    if (bananasInterval) {
      clearInterval(bananasInterval);
      bananasInterval = null;
    }
    if (audioCtx) {
      audioCtx.close().catch(() => {});
      audioCtx = null;
    }
    if (bananasOverlay) {
      bananasOverlay.remove();
      bananasOverlay = null;
    }
  }

  // ── Message listener ──

  chrome.runtime.onMessage.addListener((message, sender, sendResponse) => {
    switch (message.type) {
      case 'BLUR':
        showBlur(message.level);
        sendResponse({ success: true });
        break;

      case 'UNBLUR':
        removeBlur();
        sendResponse({ success: true });
        break;

      case 'SHRIMP_BANANAS':
        startBananas(message.score || '??');
        sendResponse({ success: true });
        break;

      case 'SHRIMP_BANANAS_STOP':
        stopBananas();
        sendResponse({ success: true });
        break;

      default:
        break;
    }
  });
})();
