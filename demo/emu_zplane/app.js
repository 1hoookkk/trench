/**
 * app.js — Entry point. Wires DOM to all modules.
 */

import { Cascade }   from './cascade.js';
import * as interp   from './interpolate.js';
import { loadPills } from './pills.js';
import { runParity } from './parity.js';
import {
    makeImpulse, makeSineBurst, makeWhiteNoise, makeClick,
    renderOffline, toAudioBuffer,
} from './render_audio.js';
import { realFft, magnitudeDb } from './fft.js';
import { drawResponse, drawImpulse } from './plot.js';

// ---------------------------------------------------------------------------
// State
// ---------------------------------------------------------------------------

let audioCtx    = null;
let sourceNode  = null;
let cascade     = new Cascade();
let pills       = null;
let currentPill = null;
let stimType    = 'impulse';
let isPlaying   = false;
let debounceTimer = null;

const SR_OFFLINE = 44100; // for rendering (no AudioContext needed initially)

// DOM refs
const pillSelect      = document.getElementById('pill-select');
const morphSlider     = document.getElementById('morph-slider');
const qSlider         = document.getElementById('q-slider');
const btnPlay         = document.getElementById('btn-play');
const btnStop         = document.getElementById('btn-stop');
const btnParity       = document.getElementById('btn-parity');
const stimImpulse     = document.getElementById('stim-impulse');
const stimSine        = document.getElementById('stim-sine');
const stimNoise       = document.getElementById('stim-noise');
const stimClick       = document.getElementById('stim-click');
const canvasResponse  = document.getElementById('canvas-response');
const statusDiv       = document.getElementById('status');

// ---------------------------------------------------------------------------
// AudioContext
// ---------------------------------------------------------------------------

function getAudioCtx() {
    if (!audioCtx) audioCtx = new AudioContext({ sampleRate: SR_OFFLINE });
    if (audioCtx.state === 'suspended') audioCtx.resume();
    return audioCtx;
}

// ---------------------------------------------------------------------------
// Stimulus
// ---------------------------------------------------------------------------

const STIM_LEN = 4096;

function makeStimulus(sr) {
    switch (stimType) {
        case 'sine':   return makeSineBurst(440, STIM_LEN, sr);
        case 'noise':  return makeWhiteNoise(STIM_LEN, sr);
        case 'click':  return makeClick(STIM_LEN, sr);
        default:       return makeImpulse(1024, sr);  // impulse; shorter for FFT
    }
}

// ---------------------------------------------------------------------------
// Compile + render pipeline
// ---------------------------------------------------------------------------

function compileCascade(sr) {
    if (!currentPill) return;
    const morph   = parseFloat(morphSlider.value);
    const q       = parseFloat(qSlider.value) * 2 - 1; // 0..1 → -1..1
    const secs    = interp.sectionsAtMorph(currentPill, morph);
    const qed     = interp.applyLiveQ(secs, q);
    const target  = interp.compileSections(qed, sr);

    cascade = new Cascade();
    cascade.setTargets(target, 1);
    cascade.setBoost(4.0, 1);
}

function renderAndPlot() {
    if (!currentPill) return;
    const sr       = audioCtx ? audioCtx.sampleRate : SR_OFFLINE;
    compileCascade(sr);

    // Impulse for FFT
    const impulse  = makeImpulse(1024, sr);
    const fftCasc  = new Cascade();
    {
        const morph  = parseFloat(morphSlider.value);
        const q      = parseFloat(qSlider.value) * 2 - 1;
        const secs   = interp.sectionsAtMorph(currentPill, morph);
        const qed    = interp.applyLiveQ(secs, q);
        const target = interp.compileSections(qed, sr);
        fftCasc.setTargets(target, 1);
        fftCasc.setBoost(4.0, 1);
    }
    const rendered = renderOffline(fftCasc, impulse);

    // Pad to 1024
    const padded = new Float32Array(1024);
    padded.set(rendered.slice(0, 1024));

    const { re, im } = realFft(padded);
    const db         = magnitudeDb(re, im);
    drawResponse(canvasResponse, db, sr);
}

// ---------------------------------------------------------------------------
// Playback
// ---------------------------------------------------------------------------

function stopPlayback() {
    if (sourceNode) {
        try { sourceNode.stop(); } catch (_) {}
        sourceNode.disconnect();
        sourceNode = null;
    }
    isPlaying = false;
    btnPlay.textContent = 'Play';
}

async function startPlayback() {
    stopPlayback();
    const ctx = getAudioCtx();
    if (!currentPill) return;

    compileCascade(ctx.sampleRate);
    const stim    = makeStimulus(ctx.sampleRate);
    const out     = renderOffline(cascade, stim);
    const abuf    = toAudioBuffer(ctx, out);

    sourceNode    = ctx.createBufferSource();
    sourceNode.buffer = abuf;
    sourceNode.loop   = true;
    sourceNode.connect(ctx.destination);
    sourceNode.start();
    isPlaying     = true;
    btnPlay.textContent = 'Playing';
}

// ---------------------------------------------------------------------------
// Debounced update
// ---------------------------------------------------------------------------

function scheduleUpdate() {
    clearTimeout(debounceTimer);
    debounceTimer = setTimeout(() => {
        renderAndPlot();
        if (isPlaying) startPlayback();
    }, 30);
}

// ---------------------------------------------------------------------------
// Parity
// ---------------------------------------------------------------------------

async function doParity() {
    statusDiv.textContent = 'Running parity…';
    try {
        const results = await runParity(interp, SR_OFFLINE);
        const lines = results.map(r => {
            const icon = r.ok ? 'OK' : 'FAIL';
            return `${icon}  ${r.label}  maxAbsErr=${r.maxAbsErr.toExponential(3)}`;
        });
        const allOk = results.every(r => r.ok);
        statusDiv.textContent = lines.join('\n') + (allOk ? '\nAll corners green.' : '\nSome corners failed!');
        statusDiv.style.color = allOk ? '#2f4' : '#f44';
    } catch (e) {
        statusDiv.textContent = `Parity error: ${e.message}`;
        statusDiv.style.color = '#f44';
    }
}

// ---------------------------------------------------------------------------
// Init
// ---------------------------------------------------------------------------

async function init() {
    statusDiv.textContent = 'Loading pills…';
    try {
        pills = await loadPills();
    } catch (e) {
        statusDiv.textContent = `Failed to load pills: ${e.message}`;
        return;
    }

    // Populate select
    pillSelect.innerHTML = '';
    for (const name of pills.names) {
        const opt   = document.createElement('option');
        opt.value   = name;
        opt.textContent = name;
        if (name === 'hedz') opt.selected = true;
        pillSelect.appendChild(opt);
    }

    currentPill = pills.getPill(pillSelect.value || pills.names[0]);
    statusDiv.textContent = `Loaded ${pills.names.length} pills.`;
    statusDiv.style.color = '#8b8';
    renderAndPlot();

    // Pill change
    pillSelect.addEventListener('change', () => {
        currentPill = pills.getPill(pillSelect.value);
        scheduleUpdate();
    });

    // Sliders
    morphSlider.addEventListener('input', scheduleUpdate);
    qSlider.addEventListener('input', scheduleUpdate);

    // Stimulus buttons
    function setStim(type, btn) {
        stimType = type;
        [stimImpulse, stimSine, stimNoise, stimClick].forEach(b => b.classList.remove('active'));
        btn.classList.add('active');
        scheduleUpdate();
    }
    stimImpulse.addEventListener('click', () => setStim('impulse', stimImpulse));
    stimSine.addEventListener('click',    () => setStim('sine', stimSine));
    stimNoise.addEventListener('click',   () => setStim('noise', stimNoise));
    stimClick.addEventListener('click',   () => setStim('click', stimClick));

    // Play / Stop
    btnPlay.addEventListener('click', () => {
        startPlayback();
    });
    btnStop.addEventListener('click', () => {
        stopPlayback();
    });

    // Parity
    btnParity.addEventListener('click', doParity);

    // Activate impulse by default
    stimImpulse.classList.add('active');
}

init().catch(e => {
    statusDiv.textContent = `Init error: ${e.message}`;
    console.error(e);
});
