#!/usr/bin/env bash
# Linux -> Windows x64 GNU cross-build. Requires matching Rust target std + MinGW.
set -euo pipefail
export CARGO_TARGET_X86_64_PC_WINDOWS_GNU_LINKER="${CARGO_TARGET_X86_64_PC_WINDOWS_GNU_LINKER:-x86_64-w64-mingw32-gcc}"
export CC_x86_64_pc_windows_gnu="${CC_x86_64_pc_windows_gnu:-x86_64-w64-mingw32-gcc}"
export AR_x86_64_pc_windows_gnu="${AR_x86_64_pc_windows_gnu:-x86_64-w64-mingw32-ar}"
cargo build --locked --release -p hx-app --target x86_64-pc-windows-gnu
python3 scripts/package-windows.py --target x86_64-pc-windows-gnu
