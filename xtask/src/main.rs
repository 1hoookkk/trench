use std::process::Command;

fn main() -> nih_plug_xtask::Result<()> {
    let mut args = std::env::args().skip(1);
    if let Some(cmd) = args.next() {
        if cmd == "release-gate" {
            let python = std::env::var("PYTHON").unwrap_or_else(|_| "python".to_string());
            let mut child = Command::new(python);
            child.arg("tools/run_shipping_release_gate.py");
            for arg in args {
                child.arg(arg);
            }
            let status = child.status()?;
            if !status.success() {
                std::process::exit(status.code().unwrap_or(1));
            }
            return Ok(());
        }
    }
    nih_plug_xtask::main()
}
