use std::fs;
fn main() {
    let p = std::env::args().nth(1).expect("path");
    let text = fs::read_to_string(&p).unwrap();
    let cart = trench_core::Cartridge::from_json(&text).expect("load");
    println!("loaded: name={}, corners={}", cart.name, cart.corners().len());
}
