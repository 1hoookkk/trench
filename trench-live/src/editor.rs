use baseview::{WindowHandle, WindowOpenOptions, WindowScalePolicy};
use crossbeam::atomic::AtomicCell;
use nih_plug::params::persist::PersistentField;
use nih_plug::prelude::*;
use raw_window_handle::{HasRawDisplayHandle, HasRawWindowHandle, RawWindowHandle};
use serde::{Deserialize, Serialize};
use std::{
    num::{NonZeroIsize, NonZeroU32},
    ptr::NonNull,
    sync::{
        atomic::{AtomicBool, Ordering},
        Arc,
    },
};

use crate::TrenchLiveParams;
use trench_core::Cartridge;

pub const UI_WIDTH: u32 = 420;
pub const UI_HEIGHT: u32 = 260;

// ── Editor state ──

#[derive(Debug, Serialize, Deserialize)]
pub struct LiveEditorState {
    #[serde(with = "nih_plug::params::persist::serialize_atomic_cell")]
    size: AtomicCell<(u32, u32)>,
    #[serde(skip)]
    open: AtomicBool,
}

impl LiveEditorState {
    pub fn new() -> Arc<Self> {
        Arc::new(Self {
            size: AtomicCell::new((UI_WIDTH, UI_HEIGHT)),
            open: AtomicBool::new(false),
        })
    }
}

impl<'a> PersistentField<'a, LiveEditorState> for Arc<LiveEditorState> {
    fn set(&self, new_value: LiveEditorState) {
        self.size.store(new_value.size.load());
    }
    fn map<F, R>(&self, f: F) -> R
    where
        F: Fn(&LiveEditorState) -> R,
    {
        f(self)
    }
}

// ── Editor ──

pub struct LiveEditor {
    params: Arc<TrenchLiveParams>,
    state: Arc<LiveEditorState>,
    cartridge: Arc<std::sync::Mutex<Option<Cartridge>>>,
    body_name: Arc<AtomicCell<[u8; 64]>>,
}

impl LiveEditor {
    pub fn new(
        params: Arc<TrenchLiveParams>,
        state: Arc<LiveEditorState>,
        cartridge: Arc<std::sync::Mutex<Option<Cartridge>>>,
        body_name: Arc<AtomicCell<[u8; 64]>>,
    ) -> Self {
        Self { params, state, cartridge, body_name }
    }
}

impl Editor for LiveEditor {
    fn spawn(
        &self,
        parent: ParentWindowHandle,
        context: Arc<dyn GuiContext>,
    ) -> Box<dyn std::any::Any + Send> {
        let params = self.params.clone();
        let state = self.state.clone();
        let cartridge = self.cartridge.clone();
        let body_name = self.body_name.clone();
        let gui_context = context.clone();

        let window = baseview::Window::open_parented(
            &ParentAdapter(parent),
            WindowOpenOptions {
                title: "TRENCH LIVE".to_string(),
                size: baseview::Size::new(UI_WIDTH as f64, UI_HEIGHT as f64),
                scale: WindowScalePolicy::ScaleFactor(1.0),
            },
            move |window: &mut baseview::Window<'_>| -> LiveWindow {
                let target = baseview_to_softbuffer(window);
                let ctx = softbuffer::Context::new(target.clone()).expect("softbuffer ctx");
                let mut surface = softbuffer::Surface::new(&ctx, target).expect("softbuffer surface");
                surface.resize(
                    NonZeroU32::new(UI_WIDTH).unwrap(),
                    NonZeroU32::new(UI_HEIGHT).unwrap(),
                ).unwrap();

                LiveWindow {
                    gui_context: gui_context.clone(),
                    params: params.clone(),
                    cartridge: cartridge.clone(),
                    body_name: body_name.clone(),
                    _ctx: Some(ctx),
                    surface: Some(surface),
                    fb: vec![0xFF1A1A1A; (UI_WIDTH * UI_HEIGHT) as usize],
                    drag: None,
                    last_x: 0.0,
                    last_y: 0.0,
                }
            },
        );

        self.state.open.store(true, Ordering::Release);
        Box::new(EditorHandle { state, window })
    }

    fn size(&self) -> (u32, u32) {
        (UI_WIDTH, UI_HEIGHT)
    }
    fn set_scale_factor(&self, _: f32) -> bool { false }
    fn param_value_changed(&self, _: &str, _: f32) {}
    fn param_modulation_changed(&self, _: &str, _: f32) {}
    fn param_values_changed(&self) {}
}

struct EditorHandle {
    state: Arc<LiveEditorState>,
    window: WindowHandle,
}
unsafe impl Send for EditorHandle {}
impl Drop for EditorHandle {
    fn drop(&mut self) {
        self.state.open.store(false, Ordering::Release);
        self.window.close();
    }
}

// ── Window ──

const BG: u32      = 0xFF1A1A1A;
const SCOPE_BG: u32 = 0xFF0A140A;
const GRID: u32     = 0xFF1A2A1A;
const ZERO: u32     = 0xFF335533;
const CURVE: u32    = 0xFFAAE0DE;
const SL_BG: u32    = 0xFF333333;
const SL_FG: u32    = 0xFF00CC88;
const DIM: u32      = 0xFF666666;
const BRIGHT: u32   = 0xFFCCCCCC;
const NAME_C: u32   = 0xFF00CC88;

const SX: u32 = 10; const SY: u32 = 10;
const SW: u32 = 400; const SH: u32 = 140;
const MY: u32 = 162; const QY: u32 = 196;
const SLX: u32 = 60; const SLW: u32 = 280; const SLH: u32 = 18;
const NY: u32 = 235;

#[derive(Clone, Copy)]
enum Drag { Morph, Q }

struct LiveWindow {
    gui_context: Arc<dyn GuiContext>,
    params: Arc<TrenchLiveParams>,
    cartridge: Arc<std::sync::Mutex<Option<Cartridge>>>,
    body_name: Arc<AtomicCell<[u8; 64]>>,
    _ctx: Option<softbuffer::Context<SbAdapter>>,
    surface: Option<softbuffer::Surface<SbAdapter, SbAdapter>>,
    fb: Vec<u32>,
    drag: Option<Drag>,
    last_x: f64,
    last_y: f64,
}

impl LiveWindow {
    fn px(&mut self, x: u32, y: u32, c: u32) {
        if x < UI_WIDTH && y < UI_HEIGHT {
            self.fb[(y * UI_WIDTH + x) as usize] = c;
        }
    }

    fn rect(&mut self, x: u32, y: u32, w: u32, h: u32, c: u32) {
        for dy in 0..h { for dx in 0..w { self.px(x+dx, y+dy, c); } }
    }

    fn text3x5(&mut self, x: u32, y: u32, s: &str, c: u32) {
        static G: &[(char, [u8;5])] = &[
            ('0',[7,5,5,5,7]),('1',[2,6,2,2,7]),('2',[7,1,7,4,7]),('3',[7,1,7,1,7]),
            ('4',[5,5,7,1,1]),('5',[7,4,7,1,7]),('6',[7,4,7,5,7]),('7',[7,1,2,2,2]),
            ('8',[7,5,7,5,7]),('9',[7,5,7,1,7]),('%',[5,1,2,4,5]),('.',[0,0,0,0,2]),
            (' ',[0,0,0,0,0]),('-',[0,0,7,0,0]),('M',[5,7,7,5,5]),('O',[7,5,5,5,7]),
            ('R',[6,5,6,5,5]),('P',[6,5,6,4,4]),('H',[5,5,7,5,5]),('Q',[7,5,5,7,1]),
            (':',[0,2,0,2,0]),('A',[2,5,7,5,5]),('B',[6,5,6,5,6]),('C',[7,4,4,4,7]),
            ('D',[6,5,5,5,6]),('E',[7,4,7,4,7]),('F',[7,4,6,4,4]),('G',[7,4,5,5,7]),
            ('I',[7,2,2,2,7]),('K',[5,5,6,5,5]),('L',[4,4,4,4,7]),('N',[5,7,7,5,5]),
            ('S',[7,4,7,1,7]),('T',[7,2,2,2,2]),('U',[5,5,5,5,7]),('V',[5,5,5,5,2]),
            ('W',[5,5,7,7,5]),('X',[5,5,2,5,5]),('Y',[5,5,2,2,2]),('_',[0,0,0,0,7]),
        ];
        let mut cx = x;
        for ch in s.chars() {
            let u = ch.to_ascii_uppercase();
            if let Some((_, g)) = G.iter().find(|(k,_)| *k == u || *k == ch) {
                for (r, bits) in g.iter().enumerate() {
                    for col in 0..3 {
                        if bits & (1 << (2-col)) != 0 { self.px(cx+col, y+r as u32, c); }
                    }
                }
            }
            cx += 4;
        }
    }

    fn render(&mut self) {
        self.fb.fill(BG);

        // Scope
        self.rect(SX, SY, SW, SH, SCOPE_BG);
        for &f in &[100.0_f64, 1000.0, 10000.0] {
            let t = (f / 20.0).ln() / (20000.0_f64 / 20.0).ln();
            let x = SX + (t * SW as f64) as u32;
            for y in SY..SY+SH { self.px(x, y, GRID); }
        }
        for &db in &[-40.0, -20.0, 0.0, 20.0] {
            let t = 1.0 - (db - (-60.0)) / (40.0 - (-60.0));
            let y = SY + (t * SH as f64) as u32;
            let col = if db == 0.0 { ZERO } else { GRID };
            for x in SX..SX+SW { self.px(x, y, col); }
        }

        // Curve — compute y-values into local vec first to avoid borrow conflict
        let morph = self.params.morph.value() as f64;
        let q = self.params.q.value() as f64;
        let curve_ys: Vec<u32> = {
            let guard = self.cartridge.lock();
            if let Ok(ref g) = guard {
                if let Some(ref cart) = **g {
                    let coeffs = cart.interpolate(morph, q);
                    (0..SW).map(|i| {
                        let t = i as f64 / (SW - 1) as f64;
                        let freq = 20.0 * (1000.0_f64).powf(t);
                        let w = 2.0 * std::f64::consts::PI * freq / 44100.0;
                        let (zr, zi) = (w.cos(), -w.sin());
                        let (z2r, z2i) = ((2.0*w).cos(), -(2.0*w).sin());
                        let mut mag = 1.0_f64;
                        for stage in &coeffs {
                            let (c0,c1,c2,c3,c4) = (stage[0],stage[1],stage[2],stage[3],stage[4]);
                            if c0.abs() < 1e-12 && c1.abs() < 1e-12 { continue; }
                            let nr = c0 + c1*zr + c2*z2r;
                            let ni = c1*zi + c2*z2i;
                            let dr = 1.0 + c3*zr + c4*z2r;
                            let di = c3*zi + c4*z2i;
                            mag *= (nr*nr+ni*ni).sqrt() / (dr*dr+di*di).sqrt().max(1e-30);
                        }
                        let db = (20.0 * mag.max(1e-30).log10()).clamp(-60.0, 40.0);
                        let yt = 1.0 - (db - (-60.0)) / (40.0 - (-60.0));
                        SY + (yt * SH as f64).clamp(0.0, (SH-1) as f64) as u32
                    }).collect()
                } else { vec![] }
            } else { vec![] }
        };
        let mut prev_y: Option<u32> = None;
        for (i, &py) in curve_ys.iter().enumerate() {
            let px = SX + i as u32;
            self.px(px, py, CURVE);
            if let Some(prev) = prev_y {
                let (a, b) = if prev < py { (prev, py) } else { (py, prev) };
                for yy in a..=b { self.px(px, yy, CURVE); }
            }
            prev_y = Some(py);
        }

        // Morph slider
        let mv = self.params.morph.value();
        self.text3x5(10, MY+5, "MORPH", DIM);
        self.rect(SLX, MY, SLW, SLH, SL_BG);
        self.rect(SLX, MY, (mv * SLW as f32) as u32, SLH, SL_FG);
        self.text3x5(SLX+SLW+8, MY+5, &format!("{}%", (mv*100.0) as i32), BRIGHT);

        // Q slider
        let qv = self.params.q.value();
        self.text3x5(10, QY+5, "Q", DIM);
        self.rect(SLX, QY, SLW, SLH, SL_BG);
        self.rect(SLX, QY, (qv * SLW as f32) as u32, SLH, SL_FG);
        self.text3x5(SLX+SLW+8, QY+5, &format!("{}%", (qv*100.0) as i32), BRIGHT);

        // Body name
        let buf = self.body_name.load();
        let len = buf.iter().position(|&b| b == 0).unwrap_or(64);
        let name = std::str::from_utf8(&buf[..len]).unwrap_or("---");
        self.text3x5(10, NY, name, NAME_C);

        // Blit
        if let Some(ref mut s) = self.surface {
            if let Ok(mut b) = s.buffer_mut() {
                b.copy_from_slice(&self.fb);
                let _ = b.present();
            }
        }
    }

    fn slider_hit(&self, x: f64, y: f64) -> Option<Drag> {
        if x >= SLX as f64 && x <= (SLX+SLW) as f64 {
            if y >= MY as f64 && y <= (MY+SLH) as f64 { return Some(Drag::Morph); }
            if y >= QY as f64 && y <= (QY+SLH) as f64 { return Some(Drag::Q); }
        }
        None
    }

    fn apply_drag(&self, x: f64) {
        let t = ((x - SLX as f64) / SLW as f64).clamp(0.0, 1.0) as f32;
        match self.drag {
            Some(Drag::Morph) => {
                let ptr = self.params.morph.as_ptr();
                unsafe {
                    self.gui_context.raw_begin_set_parameter(ptr);
                    self.gui_context.raw_set_parameter_normalized(ptr, t);
                    self.gui_context.raw_end_set_parameter(ptr);
                }
            }
            Some(Drag::Q) => {
                let ptr = self.params.q.as_ptr();
                unsafe {
                    self.gui_context.raw_begin_set_parameter(ptr);
                    self.gui_context.raw_set_parameter_normalized(ptr, t);
                    self.gui_context.raw_end_set_parameter(ptr);
                }
            }
            None => {}
        }
    }
}

impl baseview::WindowHandler for LiveWindow {
    fn on_frame(&mut self, _: &mut baseview::Window) {
        self.render();
    }
    fn on_event(&mut self, _: &mut baseview::Window, event: baseview::Event) -> baseview::EventStatus {
        match event {
            baseview::Event::Mouse(baseview::MouseEvent::CursorMoved { position, .. }) => {
                self.last_x = position.x;
                self.last_y = position.y;
                if self.drag.is_some() { self.apply_drag(position.x); }
                baseview::EventStatus::Captured
            }
            baseview::Event::Mouse(baseview::MouseEvent::ButtonPressed { button: baseview::MouseButton::Left, .. }) => {
                self.drag = self.slider_hit(self.last_x, self.last_y);
                if self.drag.is_some() { self.apply_drag(self.last_x); }
                baseview::EventStatus::Captured
            }
            baseview::Event::Mouse(baseview::MouseEvent::ButtonReleased { button: baseview::MouseButton::Left, .. }) => {
                self.drag = None;
                baseview::EventStatus::Captured
            }
            _ => baseview::EventStatus::Ignored,
        }
    }
}

// ── Window handle adapters (copied from trench-plugin) ──

struct ParentAdapter(nih_plug::editor::ParentWindowHandle);
unsafe impl HasRawWindowHandle for ParentAdapter {
    fn raw_window_handle(&self) -> RawWindowHandle {
        match self.0 {
            nih_plug::editor::ParentWindowHandle::X11Window(w) => {
                let mut h = raw_window_handle::XcbWindowHandle::empty(); h.window = w;
                RawWindowHandle::Xcb(h)
            }
            nih_plug::editor::ParentWindowHandle::AppKitNsView(v) => {
                let mut h = raw_window_handle::AppKitWindowHandle::empty(); h.ns_view = v;
                RawWindowHandle::AppKit(h)
            }
            nih_plug::editor::ParentWindowHandle::Win32Hwnd(hwnd) => {
                let mut h = raw_window_handle::Win32WindowHandle::empty(); h.hwnd = hwnd;
                RawWindowHandle::Win32(h)
            }
        }
    }
}

type SbAdapter = SoftbufferAdapter;
#[derive(Clone)]
struct SoftbufferAdapter {
    dh: raw_window_handle_06::RawDisplayHandle,
    wh: raw_window_handle_06::RawWindowHandle,
}
impl raw_window_handle_06::HasDisplayHandle for SoftbufferAdapter {
    fn display_handle(&self) -> Result<raw_window_handle_06::DisplayHandle<'_>, raw_window_handle_06::HandleError> {
        unsafe { Ok(raw_window_handle_06::DisplayHandle::borrow_raw(self.dh)) }
    }
}
impl raw_window_handle_06::HasWindowHandle for SoftbufferAdapter {
    fn window_handle(&self) -> Result<raw_window_handle_06::WindowHandle<'_>, raw_window_handle_06::HandleError> {
        unsafe { Ok(raw_window_handle_06::WindowHandle::borrow_raw(self.wh)) }
    }
}

fn baseview_to_softbuffer(window: &baseview::Window<'_>) -> SoftbufferAdapter {
    let rdh = window.raw_display_handle();
    let rwh = window.raw_window_handle();
    SoftbufferAdapter {
        dh: match rdh {
            raw_window_handle::RawDisplayHandle::Windows(_) =>
                raw_window_handle_06::RawDisplayHandle::Windows(raw_window_handle_06::WindowsDisplayHandle::new()),
            raw_window_handle::RawDisplayHandle::Xlib(h) =>
                raw_window_handle_06::RawDisplayHandle::Xlib(raw_window_handle_06::XlibDisplayHandle::new(NonNull::new(h.display), h.screen)),
            raw_window_handle::RawDisplayHandle::AppKit(_) =>
                raw_window_handle_06::RawDisplayHandle::AppKit(raw_window_handle_06::AppKitDisplayHandle::new()),
            _ => todo!(),
        },
        wh: match rwh {
            raw_window_handle::RawWindowHandle::Win32(h) => {
                let mut r = raw_window_handle_06::Win32WindowHandle::new(NonZeroIsize::new(h.hwnd as isize).unwrap());
                r.hinstance = NonZeroIsize::new(h.hinstance as isize);
                raw_window_handle_06::RawWindowHandle::Win32(r)
            }
            raw_window_handle::RawWindowHandle::Xlib(h) =>
                raw_window_handle_06::RawWindowHandle::Xlib(raw_window_handle_06::XlibWindowHandle::new(h.window)),
            raw_window_handle::RawWindowHandle::AppKit(h) =>
                raw_window_handle_06::RawWindowHandle::AppKit(raw_window_handle_06::AppKitWindowHandle::new(NonNull::new(h.ns_view).unwrap())),
            _ => todo!(),
        },
    }
}
