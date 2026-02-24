// Offscreen document: camera capture + MoveNet inference
// This file is bundled by webpack into offscreen/offscreen-bundle.js

import * as poseDetection from '@tensorflow-models/pose-detection';
import '@tensorflow/tfjs-core';
import '@tensorflow/tfjs-converter';
import '@tensorflow/tfjs-backend-webgl';

import { Peer } from 'peerjs';

import { evaluatePosture, calculateFrontMetrics, calculateSideMetrics, createBaseline } from '../offscreen/posture-engine.js';
import { MSG } from '../shared/constants.js';

let detector = null;
const streams = {};  // { front: MediaStream, side: MediaStream }
const videos = {};   // { front: HTMLVideoElement, side: HTMLVideoElement }
const canvases = {}; // { front: HTMLCanvasElement, side: HTMLCanvasElement }

// ── WebRTC (phone companion) ──

let peerInstance = null;
let peerDataConn = null;

async function initDetector() {
  if (detector) return detector;

  // Use WebGL backend
  const tf = await import('@tensorflow/tfjs-core');
  await tf.setBackend('webgl');
  await tf.ready();

  detector = await poseDetection.createDetector(
    poseDetection.SupportedModels.MoveNet,
    {
      modelType: poseDetection.movenet.modelType.SINGLEPOSE_LIGHTNING,
      enableSmoothing: false
    }
  );

  console.log('[ShrimpWatch] MoveNet Lightning loaded');
  return detector;
}

async function getOrCreateStream(role, deviceId) {
  if (streams[role] && streams[role].active) {
    return streams[role];
  }

  const stream = await navigator.mediaDevices.getUserMedia({
    video: {
      deviceId: { exact: deviceId },
      width: { ideal: 640 },
      height: { ideal: 480 },
      frameRate: { ideal: 5 }
    }
  });

  streams[role] = stream;

  const video = document.createElement('video');
  video.srcObject = stream;
  video.autoplay = true;
  video.playsInline = true;
  await new Promise((resolve, reject) => {
    video.onloadeddata = resolve;
    video.onerror = reject;
    setTimeout(() => reject(new Error('Video load timeout')), 10000);
  });
  videos[role] = video;

  const canvas = document.createElement('canvas');
  canvas.width = video.videoWidth || 640;
  canvas.height = video.videoHeight || 480;
  canvases[role] = canvas;

  // Listen for track ending (camera disconnected)
  stream.getVideoTracks()[0].addEventListener('ended', () => {
    console.warn(`[ShrimpWatch] Camera ${role} disconnected`);
    chrome.runtime.sendMessage({ type: MSG.CAMERA_DISCONNECTED, role });
    delete streams[role];
    delete videos[role];
    delete canvases[role];
  });

  return stream;
}

async function captureAndAnalyze(cameras, calibration) {
  const det = await initDetector();
  const poseData = {};

  // Check if phone side camera is connected via WebRTC
  // If so, add it as the 'side' role even if no local deviceId
  const roles = { ...cameras };
  if (!roles.side && isPhoneConnected()) {
    roles.side = '__webrtc__';
  }

  for (const [role, deviceId] of Object.entries(roles)) {
    if (!deviceId) continue;

    try {
      // For WebRTC streams, skip getUserMedia — stream already exists
      if (deviceId !== '__webrtc__') {
        await getOrCreateStream(role, deviceId);
      }

      const video = videos[role];
      const canvas = canvases[role];
      if (!video || !canvas) continue;

      const ctx = canvas.getContext('2d');
      ctx.drawImage(video, 0, 0, canvas.width, canvas.height);
      const poses = await det.estimatePoses(canvas);

      if (poses.length > 0 && poses[0].keypoints) {
        poseData[role] = {
          keypoints: poses[0].keypoints,
          width: canvas.width,
          height: canvas.height
        };
      }
    } catch (err) {
      console.error(`[ShrimpWatch] Error capturing from ${role} camera:`, err);
    }
  }

  return evaluatePosture(poseData, calibration);
}

async function captureCalibrationFrames(cameraRole, deviceId, numFrames, intervalMs) {
  const det = await initDetector();
  // If deviceId is __webrtc__, the phone stream is already in streams/videos/canvases
  if (deviceId !== '__webrtc__') {
    await getOrCreateStream(cameraRole, deviceId);
  }
  if (!videos[cameraRole] || !canvases[cameraRole]) {
    throw new Error(`No ${cameraRole} camera stream available`);
  }

  const metricsList = [];
  for (let i = 0; i < numFrames; i++) {
    const video = videos[cameraRole];
    const canvas = canvases[cameraRole];
    const ctx = canvas.getContext('2d');

    ctx.drawImage(video, 0, 0, canvas.width, canvas.height);
    const poses = await det.estimatePoses(canvas);

    if (poses.length > 0 && poses[0].keypoints) {
      const metrics = cameraRole === 'front'
        ? calculateFrontMetrics(poses[0].keypoints, canvas.width, canvas.height)
        : calculateSideMetrics(poses[0].keypoints, canvas.width, canvas.height);
      metricsList.push(metrics);
    }

    if (i < numFrames - 1) {
      await new Promise(r => setTimeout(r, intervalMs));
    }
  }

  return createBaseline(metricsList);
}

function stopStreams(roles) {
  const toStop = roles || Object.keys(streams);
  for (const role of toStop) {
    if (streams[role]) {
      streams[role].getTracks().forEach(t => t.stop());
      delete streams[role];
      delete videos[role];
      delete canvases[role];
    }
  }
}

async function enumerateCameras() {
  const devices = await navigator.mediaDevices.enumerateDevices();
  return devices
    .filter(d => d.kind === 'videoinput')
    .map(d => ({ deviceId: d.deviceId, label: d.label || `Camera ${d.deviceId.slice(0, 8)}` }));
}

// ── WebRTC: receive phone camera stream ──

function initPeer() {
  return new Promise((resolve, reject) => {
    if (peerInstance && !peerInstance.destroyed) {
      resolve(peerInstance.id);
      return;
    }

    // Generate a stable-ish peer ID based on extension ID
    const peerId = 'shrimpwatch-' + Math.random().toString(36).slice(2, 10);
    peerInstance = new Peer(peerId);

    peerInstance.on('open', (id) => {
      console.log('[ShrimpWatch] WebRTC peer ready:', id);

      // Listen for incoming calls from phone companion
      peerInstance.on('call', (call) => {
        console.log('[ShrimpWatch] Incoming call from phone companion');
        // Answer with no stream (we only receive)
        call.answer();

        call.on('stream', (remoteStream) => {
          console.log('[ShrimpWatch] Phone camera stream received');

          // Store as the 'side' stream
          streams['side'] = remoteStream;

          const video = document.createElement('video');
          video.srcObject = remoteStream;
          video.autoplay = true;
          video.playsInline = true;
          video.onloadeddata = () => {
            videos['side'] = video;
            const canvas = document.createElement('canvas');
            canvas.width = video.videoWidth || 640;
            canvas.height = video.videoHeight || 480;
            canvases['side'] = canvas;

            // Notify service worker
            chrome.runtime.sendMessage({
              type: 'PHONE_CAMERA_CONNECTED',
              target: 'background'
            });

            // Send confirmation to phone via data channel
            if (peerDataConn) {
              peerDataConn.send({ type: 'connected' });
            }
          };

          remoteStream.getVideoTracks()[0].addEventListener('ended', () => {
            console.warn('[ShrimpWatch] Phone camera stream ended');
            delete streams['side'];
            delete videos['side'];
            delete canvases['side'];
            chrome.runtime.sendMessage({
              type: MSG.CAMERA_DISCONNECTED,
              role: 'side'
            });
          });
        });

        call.on('close', () => {
          console.log('[ShrimpWatch] Phone call closed');
          delete streams['side'];
          delete videos['side'];
          delete canvases['side'];
        });
      });

      // Listen for data connections
      peerInstance.on('connection', (conn) => {
        peerDataConn = conn;
        conn.on('open', () => {
          console.log('[ShrimpWatch] Data channel open with phone');
        });
      });

      resolve(id);
    });

    peerInstance.on('error', (err) => {
      console.error('[ShrimpWatch] Peer error:', err);
      reject(err);
    });

    // Timeout after 15 seconds
    setTimeout(() => reject(new Error('Peer connection timeout')), 15000);
  });
}

function destroyPeer() {
  if (peerDataConn) {
    peerDataConn.close();
    peerDataConn = null;
  }
  if (peerInstance) {
    peerInstance.destroy();
    peerInstance = null;
  }
  // Clean up side stream if it was from WebRTC
  if (streams['side'] && !streams['side']._isLocal) {
    delete streams['side'];
    delete videos['side'];
    delete canvases['side'];
  }
}

function isPhoneConnected() {
  return !!(streams['side'] && videos['side'] && streams['side'].active);
}

// ── Message handler ──

chrome.runtime.onMessage.addListener((message, sender, sendResponse) => {
  if (message.target && message.target !== 'offscreen') return false;

  switch (message.type) {
    case MSG.CAPTURE_AND_ANALYZE:
      captureAndAnalyze(message.cameras, message.calibration)
        .then(result => sendResponse({ type: MSG.POSTURE_RESULT, ...result }))
        .catch(err => sendResponse({ type: MSG.POSTURE_ERROR, error: err.message }));
      return true; // async

    case MSG.STOP_STREAMS:
      stopStreams(message.roles);
      sendResponse({ type: MSG.STREAMS_STOPPED });
      return false;

    case MSG.ENUMERATE_CAMERAS:
      enumerateCameras()
        .then(cameras => sendResponse({ type: MSG.CAMERA_LIST, cameras }))
        .catch(err => sendResponse({ type: MSG.CAMERA_LIST, cameras: [], error: err.message }));
      return true; // async

    case 'CALIBRATE':
      captureCalibrationFrames(message.cameraRole, message.deviceId, message.numFrames || 5, message.intervalMs || 500)
        .then(baseline => sendResponse({ type: 'CALIBRATION_RESULT', baseline, cameraRole: message.cameraRole }))
        .catch(err => sendResponse({ type: MSG.POSTURE_ERROR, error: err.message }));
      return true; // async

    case 'INIT_MODEL':
      initDetector()
        .then(() => sendResponse({ type: 'MODEL_READY' }))
        .catch(err => sendResponse({ type: MSG.POSTURE_ERROR, error: err.message }));
      return true; // async

    case 'INIT_PEER':
      initPeer()
        .then(peerId => sendResponse({ type: 'PEER_READY', peerId }))
        .catch(err => sendResponse({ type: 'PEER_ERROR', error: err.message }));
      return true; // async

    case 'DESTROY_PEER':
      destroyPeer();
      sendResponse({ type: 'PEER_DESTROYED' });
      return false;

    case 'PHONE_STATUS':
      sendResponse({ type: 'PHONE_STATUS', connected: isPhoneConnected() });
      return false;

    case 'PING':
      sendResponse({ type: 'PONG' });
      return false;

    default:
      return false;
  }
});

// Signal that offscreen document is ready
console.log('[ShrimpWatch] Offscreen document loaded');
