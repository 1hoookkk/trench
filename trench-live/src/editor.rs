//! TRENCH LIVE editor — vizia GUI thread only.
//!
//! This is a scaffold for the R1 spike. The layout is a structural
//! placeholder for the v1 faceplate: 4-body selector row, LCD
//! placeholder, MORPH / Q / SPACE / TRIG rollers, REV / ENV buttons.
//! MORPH and Q are wired to the live plugin params. SPACE / TRIG /
//! REV / ENV are inert placeholders until Task R2 adds the
//! corresponding DSP stages and params.
//!
//! AUDIO-THREAD NOTE: this module runs on the GUI thread. Nothing in
//! here may be called from `process()`. The DSP doctrine in
//! `lib.rs::process` stands untouched.
use nih_plug::prelude::Editor;
use nih_plug_vizia::vizia::prelude::*;
use nih_plug_vizia::widgets::*;
use nih_plug_vizia::{assets, create_vizia_editor, ViziaState, ViziaTheming};
use std::sync::Arc;

use crate::TrenchLiveParams;

#[derive(Lens)]
struct Data {
    params: Arc<TrenchLiveParams>,
}

impl Model for Data {}

pub(crate) fn default_state() -> Arc<ViziaState> {
    ViziaState::new(|| (340, 510))
}

pub(crate) fn create(
    params: Arc<TrenchLiveParams>,
    editor_state: Arc<ViziaState>,
) -> Option<Box<dyn Editor>> {
    create_vizia_editor(editor_state, ViziaTheming::Custom, move |cx, _| {
        assets::register_noto_sans_light(cx);
        assets::register_noto_sans_thin(cx);

        Data {
            params: params.clone(),
        }
        .build(cx);

        VStack::new(cx, |cx| {
            // --- Body selector row (stub — Task R2+) -----------------
            HStack::new(cx, |cx| {
                for i in 0..4 {
                    Label::new(cx, &format!("B{}", i + 1))
                        .width(Pixels(60.0))
                        .height(Pixels(40.0))
                        .background_color(Color::rgb(40, 40, 40))
                        .border_color(Color::rgb(80, 80, 80))
                        .border_width(Pixels(1.0))
                        .child_space(Stretch(1.0));
                }
            })
            .height(Pixels(40.0))
            .col_between(Pixels(8.0));

            // --- LCD placeholder rectangle ----------------------------
            // Pale mint body per identity doctrine. The real spectrum
            // renderer lands later — this rectangle stakes the space.
            Element::new(cx)
                .width(Stretch(1.0))
                .height(Pixels(100.0))
                .background_color(Color::rgb(180, 230, 200))
                .border_color(Color::black())
                .border_width(Pixels(2.0));

            // --- Roller row 1: MORPH + Q (wired) ---------------------
            HStack::new(cx, |cx| {
                VStack::new(cx, |cx| {
                    Label::new(cx, "MORPH").height(Pixels(20.0));
                    ParamSlider::new(cx, Data::params, |p| &p.morph);
                })
                .row_between(Pixels(4.0));

                VStack::new(cx, |cx| {
                    Label::new(cx, "Q").height(Pixels(20.0));
                    ParamSlider::new(cx, Data::params, |p| &p.q);
                })
                .row_between(Pixels(4.0));
            })
            .col_between(Pixels(12.0));

            // --- Roller row 2: SPACE + TRIG (stubs) ------------------
            HStack::new(cx, |cx| {
                VStack::new(cx, |cx| {
                    Label::new(cx, "SPACE").height(Pixels(20.0));
                    // TODO Task R2: ParamSlider::new(cx, Data::params, |p| &p.space);
                    Element::new(cx)
                        .height(Pixels(24.0))
                        .background_color(Color::rgb(60, 60, 60))
                        .border_color(Color::rgb(30, 30, 30))
                        .border_width(Pixels(1.0));
                })
                .row_between(Pixels(4.0));

                VStack::new(cx, |cx| {
                    Label::new(cx, "TRIG").height(Pixels(20.0));
                    // TODO Task R2: ParamSlider::new(cx, Data::params, |p| &p.trig);
                    Element::new(cx)
                        .height(Pixels(24.0))
                        .background_color(Color::rgb(60, 60, 60))
                        .border_color(Color::rgb(30, 30, 30))
                        .border_width(Pixels(1.0));
                })
                .row_between(Pixels(4.0));
            })
            .col_between(Pixels(12.0));

            // --- Button row: REV + ENV (stubs) -----------------------
            HStack::new(cx, |cx| {
                // TODO Task R2: wire to reverb on/off param.
                Label::new(cx, "REV")
                    .width(Pixels(80.0))
                    .height(Pixels(32.0))
                    .background_color(Color::rgb(180, 40, 40))
                    .border_color(Color::black())
                    .border_width(Pixels(1.0))
                    .child_space(Stretch(1.0));
                // TODO Task R2: wire to envelope on/off param.
                Label::new(cx, "ENV")
                    .width(Pixels(80.0))
                    .height(Pixels(32.0))
                    .background_color(Color::rgb(60, 60, 60))
                    .border_color(Color::black())
                    .border_width(Pixels(1.0))
                    .child_space(Stretch(1.0));
            })
            .height(Pixels(32.0))
            .col_between(Pixels(12.0));

            ResizeHandle::new(cx);
        })
        // transparent-case dark base — identity refinement is later work.
        .background_color(Color::rgb(25, 25, 25))
        .child_space(Pixels(12.0))
        .row_between(Pixels(12.0));
    })
}
